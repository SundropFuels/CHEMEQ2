
import ceq2_1D as kin
import numpy as np
import pandas as pd
import cantera as ct
from collections import OrderedDict
import unitConversion as uc

if __name__ == "__main__":

    #Define the size of the system

    conv = uc.UnitConverter()
    D = conv.convert_units(1.5, 'in', 'm')
    h = 30.0		#W/m^2-K
    Tw = conv.convert_units(1400, 'C', 'K')
    L = conv.convert_units(24, 'in', 'm')
    eps = [0.0, 0.3]
    N2_ratio = [0.0, 15.0]
    H2O_ratio = [0.0, 1.0, 3.0]
    CO2_ratio = [0.0, 1.0]
     
    factor = 1.0/24.0/60.0/1000.0
    #
    flow = [1.0*factor, 3.0*factor, 5.0*factor] 			#SLPM to kmol/s
    

    z = np.arange(0, L, 0.001)
    ph = ct.Solution('POLIMI_PRF_PAH_HT_1407.cti')
    
    Tw = [1300, 1400, 1500]
    P = [conv.convert_units(75.0, 'psig', 'Pa')]
    
    for pressure in P:
        
        for temp in Tw:
            for ep in eps:
                for N2r in N2_ratio:
                    for CO2r in CO2_ratio:
                        for H2Or in H2O_ratio:   
                            for CH4_flow in flow:
                                Y = 'CH4:%s,N2:%s,H2O:%s,CO2:%s' % (CH4_flow,CH4_flow*N2r,CH4_flow*H2Or,CH4_flow*CO2r)
                                ph.TPX =273.15+25, pressure, Y
                                rx = kin.Tube(D = D, h = h, Tw = temp, eps = ep)
                                solver = kin.ChemEQ2Solver(ct_phase = ph, tube = rx, mode = 'convective')
                                solver.initialize(z, CH4_flow*(1+N2r+H2Or+CO2r))
                                solver.conv_eps = 1E-1
                                solver.solve(Nc=3)

                                outframe = solver.y.copy()
                                outframe['Enthalpy'] = np.zeros(len(solver.y['H2']))
                                for i in solver.y.index:
                                    ph.TPX = 273.15+temp, pressure, solver.ct_str(solver.y.loc[i,:])
                                    outframe.loc[i,'Enthalpy'] = ph.enthalpy_mass

                                filename = "./hw_lab/methane_crack_hotwall_T%s_P%s_CH4%0.2f_N2%s_H2O%s_CO2%s_eps%s.csv" % (temp, pressure/101325.0,CH4_flow,N2r,H2Or,CO2r,ep)
                                outframe.to_csv(filename)
    
    
