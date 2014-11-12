import numpy as np
import pandas as pd
import cantera as ct
from collections import OrderedDict
import sys


class ceq2Exception(Exception):
    def __init__(self,value):
        self.value = value

    def __str__(self):
        return repr(self.value)

class NaNError(ceq2Exception):
    pass

class BadInputError(ceq2Exception):
    pass


def increasing(x):
    """Determines if a numpy array is strictly increasing"""
    x_shift = x[1:]
    a = x[:-1]
    return (x_shift-a>0).all()


class ChemEQ2Solver:

    def __init__(self, ct_phase):

        #should check on the phase here
        self.ct_phase = ct_phase
        self.conv_eps = 1E-2
        self.dt_eps = 1E-2
        self.T = self.ct_phase.T
        self.P = self.ct_phase.P
        self.stability_adjustment = False
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
        
        self.y = pd.DataFrame(data = data, index = self.t)
        self.y0 = self.y.iloc[0,:].values
        #print self.y0        


    def solve(self, Nc):
        self.Nc = Nc
        self.t_now = self.t[0]
        
        for t in self.t[1:]:   #skip 0
            self.estimate_dt()
            while self.t_now < t:
                #print "Solving major step for time %s" % t
                y_last = self.y0
                t_last = self.t_now
                self.recursion_depth = 0
                y = self.solve_ts()
                self.y0 = y		#set the current value to the most recent timestep
            
            #print "current y:\t%s" % y
            #actually, probably need to do some interpolation here if t_now > t
            y_int = (y - y_last)/(self.t_now - t_last)*(t-t_last) + y_last
            self.y.loc[t,:] = y_int


    def solve_ts(self):
        self.yp = self.calc_yp()
        self.yc = self.yp
        sys.stdout.write("current time:\t%s\r" % (self.t_now))
        N = 0
        #print "Current timestep:\t%s" % self.t_now
        #print "Current dt:\t%s" % self.dt
        self.yc_history = {}
        while N < self.Nc:
            self.yc = self.calc_yc(self.yc)
            self.yc_history[N] = self.yc
            N += 1
        #print self.recursion_depth
        if self.converged() and self.stable():
            self.t_now += self.dt
            self.adjust_dt()
            return self.yc
        else: 
            self.recursion_depth += 1
            self.adjust_dt()
            self.stability_adjustment = False
            return self.solve_ts()

    def y_pc(self, y0, q, p):
        """Calculates either the predictor or the corrector, depending on the terms fed to it"""
        n = self.dt * (q-p*y0)
        d = 1.0 + self.pade(self.dt*p)*p*self.dt
        return y0 + n/d

    def calc_yp(self):
        #self.y0 holds the number of moles at this timestep for each specie in kmol
        self.ct_phase.TPX = self.T, self.P, self.ct_str(self.y0) 
        #Calculate q0 and p0
        V = np.sum(self.y0) * ct.gas_constant*self.T/self.P
        self.q0 = self.ct_phase.creation_rates*V
        self.p0 = self.ct_phase.destruction_rates*V/self.y0
        self.p0[np.logical_not(np.isfinite(self.p0))] = 0.0
        
        return self.y_pc(self.y0, self.q0, self.p0)

    def calc_yc(self, yp):
        self.ct_phase.TPX = self.T, self.P, self.ct_str(yp)
        V = np.sum(yp)*ct.gas_constant*self.T/self.P
        qp = self.ct_phase.creation_rates*V
        pp = self.ct_phase.destruction_rates*V/yp
        pp[np.logical_not(np.isfinite(pp))] = 0.0
        p_bar = 0.5*(pp+self.p0)
        alpha_bar = self.pade(p_bar*self.dt)
        q_tilde = alpha_bar * qp + (1.0-alpha_bar)*self.q0
        return self.y_pc(self.y0, q_tilde, p_bar)

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
        cc = np.abs(self.yc-self.yp)/(self.conv_eps * self.yc)

        self.sigma = np.max(cc[np.isfinite(cc)])
        return self.sigma <= 1.0

    def stable(self):
        """This checks the stability criterion"""
        #Do nothing right now -- will require some work to do the lagging and to adjust the time step correctly
        return True
        if self.Nc >= 3:
            cc = np.abs(self.yc_history[self.Nc-1] - self.yc_history[self.Nc-2]) - np.abs(self.yc_history[self.Nc-2] - self.yc_history[self.Nc-3])
        else:
            cc = np.zeros(3)  # just to make the final function workable
        retval = np.max(cc) < 0
        self.stabilty_adjustment = not retval
        #cc = np.abs(self.yc - self.yc_lag1) - np.abs(self.yc_lag1-self.yc_lag2)
        return retval

    def adjust_dt(self):


        #sqrt_sigma = self.newton_sqrt(self.sigma, 3)
        #print np.power(self.sigma, 0.5)

        if self.recursion_depth > 500:
            self.dt *= 0.001	#Go small, and go small fast
        elif self.stability_adjustment:
            d = np.abs(self.yc_history[self.Nc-1] - self.yc_history[self.Nc-2])+0.001
            n = np.abs(self.yc_history[self.Nc-2] - self.yc_history[self.Nc-3])
            r = n/d
            self.dt = self.dt * np.max(r[np.isfinite(r)])
        else:
            if self.sigma > 0:
                self.dt = self.dt * (1/self.newton_sqrt(self.sigma,3))
                #self.dt = self.dt * (1/np.power(self.sigma,0.5))   #This could be too slow -- may want to replace the sqrt function with a 3-time newton iteration
            

        if self.dt > self.dt_max:
            self.dt = self.dt_max


    def estimate_dt(self):
        V = np.sum(self.y0)*ct.gas_constant*self.T/self.P
        r = self.y0/(self.ct_phase.net_production_rates*V)
        r[np.logical_not(np.isfinite(r))] = 10000000	#where net production rate is 0, ignore
        r[r==0.0] = 10000000
        self.dt = self.dt_eps * np.min(np.abs(r))
        if self.dt > self.dt_max:
            self.dt = self.dt_max
        #if self.dt < self.dt_min:
        #    self.dt = self.dt_min
        #print self.dt
    
    def newton_sqrt(self, val, iterations):
        N=0
        ans = val
        while N<iterations:
            ans = 0.5*(ans+val/ans)
            N += 1
        return ans

    def ct_str(self, y):
        #y needs to be a data row (pandas Series object)
        s = ""
        for name in self.ct_phase.species_names:
            s += "%s:%s," % (name, y[self.ct_phase.species_index(name)])
        s = s[:-1]
        return s


