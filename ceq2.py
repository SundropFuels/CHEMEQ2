import numpy as np
import pandas as pd
import cantera as ct
from collections import OrderedDict
import sys
import copy
import time
import scipy.optimize as spo

class ceq2Exception(Exception):
    def __init__(self,value):
        self.value = value

    def __str__(self):
        return repr(self.value)

class NaNError(ceq2Exception):
    pass

class BadInputError(ceq2Exception):
    pass

class ConvectionNotSetup(ceq2Exception):
    pass


def increasing(x):
    """Determines if a numpy array is strictly increasing"""
    x_shift = x[1:]
    a = x[:-1]
    return (x_shift-a>0).all()


class ChemEQ2Solver:

    def __init__(self, ct_phase, mode = 'isothermal'):

        #should check on the phase here
        self.ct_phase = ct_phase
        self.conv_eps = 1E-1
        self.c_factor = 10.0
        self.dt_eps = 1E-2
        self.T = self.ct_phase.T
        self.P = self.ct_phase.P
        self.stability_adjustment = False
        self.mode = mode
        self.recursion_depth = 0
        self.iterations = 0
        self.cutoff_moles = 1E-20
        print "Solver set up."



    def initialize(self, t):
        self._init_t(t)
        self._init_y()
        self.dt_max = np.min(self.t[1:] - self.t[:-1])
        self.dt_min = 1E-10			#put a lower bound on this -- may cause some problems, which means we really need a temperature ramp-up
        print "Solver initialized"

    def _init_t(self, t):

        if isinstance(t, np.ndarray):
            #should check that t makes sense  #!#
            self.t = t
        else:
            raise BadInputError, "t must be a numpy array right now"
        

        if self.t[0] != 0.0:
            raise BadInputError, "time must start at 0"

        if not increasing(self.t):
            raise BadInputError, "time must be strictly increasing"
        

        


    def _init_y(self):
        #data = {key:np.zeros(len(self.t)) for key in self.ct_phase.species_names}
        data = OrderedDict()
        for key in self.ct_phase.species_names:
            data[key] = np.zeros(len(self.t))
            data[key][0] = self.ct_phase.X[self.ct_phase.species_index(key)]		#building off of mole fractions, with a basis scalable by volumetric flowrate

        self.species_len = len(self.ct_phase.species_names)

        

        #Temperature and pressure
        data['T'] = np.zeros(len(self.t))
        data['P'] = np.zeros(len(self.t))

        data['T'][0] = self.T
        data['P'][0] = self.P

        self.T_index = self.species_len
              
        self.y = pd.DataFrame(data = data, index = self.t)

        self.y0 = self.y.iloc[0,:].values    #This is a numpy array now

        #print self.y0        


    def solve(self, Nc):
        self.Nc = Nc
        self.t_now = self.t[0]
        self.yc_history = pd.DataFrame(index = range(0,self.Nc), columns = ['yc'])
        self.estimate_dt()
        for t in self.t[1:]:   #skip 0
            #self.c_step = t
            #self.estimate_dt()
            
            while self.t_now < t:
                #print "Solving major step for time %s" % t
                y_last = self.y0
                t_last = self.t_now
                #self.recursion_depth = 0
                self.solve_ts()

                if self.converged() and self.stable():
                    self.t_now += self.dt
                    self.adjust_dt()
                    y = self.yc  #these may be slow, but necessary
                    self.y0 = y  #ditto
         
                else: 
                    self.adjust_dt()		#cannot generalize this because of the stability adjustment
                    self.stability_adjustment = False

                
                    
            
            y_int = (y - y_last)/(self.t_now - t_last)*(t-t_last) + y_last
            self.y.loc[t,:] = y_int
            
        print "Solved in %s iterations" % self.iterations

    def solve_ts(self):

        self.calc_yp()
        self.yc = copy.deepcopy(self.yp)
        sys.stdout.write("current time:\t%0.5e  current dt:\t%0.5e  current temp:\t%0.5e\r" % (self.t_now, self.ct_phase.mean_molecular_weight*np.sum(self.y0[:self.species_len]), self.y0[self.T_index]))
        self.H_in = self.ct_phase.enthalpy_mass
        N = 0
        #print "Current timestep:\t%s" % self.t_now
        #print "Current dt:\t%s" % self.dt
        
        while N < self.Nc:
            self.calc_yc(self.yc)
            self.yc_history['yc'][N] = copy.deepcopy(self.yc)
            N += 1
        #print self.recursion_depth
        self.iterations += 1
        
        

    def y_pc(self, y0, q, p):
        """Calculates either the predictor or the corrector, depending on the terms fed to it"""
        
        return y0 + (self.dt * (q-p*y0))/(1.0 + self.pade(self.dt*p)*p*self.dt)

    def calc_yp(self):
        #self.y0 holds the number of moles at this timestep for each specie in kmol
        self.ct_phase.TPX = self.y0[self.T_index], self.y0[self.T_index+1], self.ct_str(self.y0) 
        #Calculate q0 and p0
        #need to strip off the molar flows first
        #V = np.sum(self.y0[:self.species_len]) * ct.gas_constant*self.y0[self.T_index]/self.y0[self.T_index+1]
        V = self.ct_phase.volume_mole*np.sum(self.y0[:self.species_len])
        self.q0 = self.ct_phase.creation_rates*V
        self.p0 = self.ct_phase.destruction_rates*V/self.y0[:self.species_len]
        self.p0[np.logical_not(np.isfinite(self.p0))] = 0.0
        self.yp = self.y_pc(self.y0[:self.species_len], self.q0, self.p0)
        #now need to solve the energy balance
        self.yp_dT = getattr(self, 'energy_balance_%s' % self.mode)(self.y0)
        
        self.yp = np.append(self.yp,self.y0[self.T_index] + self.dt * self.yp_dT)
        self.yp = np.append(self.yp,self.P)
        

    def calc_yc(self, yp):
        self.ct_phase.TPX = self.y0[self.T_index], self.y0[self.T_index+1], self.ct_str(yp)
        #V = np.sum(yp[:self.species_len])*ct.gas_constant*yp[self.T_index]/yp[self.T_index+1]
        V = self.ct_phase.volume_mole*np.sum(yp[:self.species_len])
        qp = self.ct_phase.creation_rates*V
        pp = self.ct_phase.destruction_rates*V/yp[:self.species_len]
        pp[np.logical_not(np.isfinite(pp))] = 0.0
        p_bar = 0.5*(pp+self.p0)
        alpha_bar = self.pade(p_bar*self.dt)
        q_tilde = alpha_bar * qp + (1.0-alpha_bar)*self.q0
        self.yc = self.y_pc(self.y0[:self.species_len], q_tilde, p_bar)
        
        self.yc = np.append(self.yc, self.y0[self.T_index] + self.dt * 0.5 * (self.yp_dT+ getattr(self, 'energy_balance_%s' % self.mode)(yp)))
        
        self.yc = np.append(self.yc, self.P)
        

    def pade(self, r_inv):
        r = 1.0/r_inv
        n = 180.0*r**3 + 60.0*r**2 + 11*r + 1
        d = 360.0*r**3 + 60.0*r**2 + 12*r + 1
        a = n/d
        a[np.logical_not(np.isfinite(a))] = 0.5
        return a

    def converged(self):
        """Determine if the current step is converged.  Store the convergence check."""
        #first check that we have finite values of everything -- doing this here is a convenient place to catch this error
        if not np.isfinite(self.yc).all():		#should also check >0 -- do we install something to put a floor on the levels to avoid negatives?
            raise NaNError, "The solver has encountered a non-finite solution.  Check the equations and try again."
        cc = np.abs(self.yc[self.yc>self.cutoff_moles]-self.yp[self.yc>self.cutoff_moles])/(self.conv_eps * self.yc[self.yc>self.cutoff_moles])

        self.sigma = np.max(cc[np.isfinite(cc)])
        return self.sigma <= 1.0

    def stable(self):
        """This checks the stability criterion"""
        #Do nothing right now -- will require some work to do the lagging and to adjust the time step correctly
        return True
        if self.Nc >= 3:
            cc = np.abs(self.yc_history['yc'][self.Nc-1] - self.yc_history['yc'][self.Nc-2]) - np.abs(self.yc_history['yc'][self.Nc-2] - self.yc_history['yc'][self.Nc-3])
        else:
            cc = np.zeros(3)  # just to make the final function workable
        retval = np.max(cc) < 0
        self.stabilty_adjustment = not retval
        #cc = np.abs(self.yc - self.yc_lag1) - np.abs(self.yc_lag1-self.yc_lag2)
        return retval

    def adjust_dt(self):


        #sqrt_sigma = self.newton_sqrt(self.sigma, 3)
        #print np.power(self.sigma, 0.5)

        #if self.recursion_depth > 500:
        #    self.dt *= 0.001	#Go small, and go small fast
        if self.stability_adjustment:
            print "here"
            d = np.abs(self.yc_history['yc'][self.Nc-1] - self.yc_history['yc'][self.Nc-2])+0.001
            n = np.abs(self.yc_history['yc'][self.Nc-2] - self.yc_history['yc'][self.Nc-3])
            r = n/d
            self.dt = self.dt * np.max(r[np.isfinite(r)])
        else:
            if self.sigma > 0:

                #self.dt = self.dt * (1/np.power(self.sigma,0.5))   #This could be too slow -- may want to replace the sqrt function with a 3-time newton iteration
                self.dt = self.dt * (1/self.newton_sqrt(self.sigma*self.c_factor, 3))  #This will fail the current tests (probably)


        if self.dt > self.dt_max:
            self.dt = self.dt_max

        #if self.dt > self.c_step - self.t_now:
        #    self.dt = self.c_step - self.t_now

    def estimate_dt(self):
        V = np.sum(self.y0[:self.species_len])*ct.gas_constant*self.y0[self.T_index]/self.y0[self.T_index+1]
        r = self.y0[:self.species_len]/(self.ct_phase.net_production_rates*V)
        r[np.logical_not(np.isfinite(r))] = 10000000	#where net production rate is 0, ignore
        r[r==0.0] = 10000000
        self.dt = self.dt_eps * np.min(np.abs(r))
        if self.dt > self.dt_max:
            self.dt = self.dt_max
        
            
    def newton_sqrt(self, val, iterations):
        N=0
        ans = val
        while N<iterations:
            ans = 0.5*(ans+val/ans)
            N += 1
        return ans

    def ebal(self, T):
        self.ct_phase.TPX = T, self.P, self.ct_str(self.yc)
        return self.H_in - self.ct_phase.enthalpy_mass

    def energy_balance_isothermal(self, y):
        return 0.0


    def energy_balance_adiabatic(self, y):
        #Assumes the ct_phase is already properly set -- I could do this explicitly, but that would take more processor time
        n = np.sum(y[:self.species_len])
        return -1*(np.sum(self.ct_phase.net_production_rates*self.ct_phase.volume_mole*n*self.ct_phase.partial_molar_enthalpies))/(n*self.ct_phase.cp_mole)
        



    def energy_balance_convective(self, y):
        try:
            Q = self.h * self.SA * (self.Tw - y[self.T_index])
        except AttributeError:
            raise ConvectionNotSetup, "The convection problem has not been properly set up -- you need to specify an h and a surface area"


    def ct_str(self, y):
        #y needs to be a data row (pandas Series object)
        s = ""

        for name in self.ct_phase.species_names:
            s += "%s:%s," % (name, y[self.ct_phase.species_index(name)])
        s = s[:-1]
        return s


