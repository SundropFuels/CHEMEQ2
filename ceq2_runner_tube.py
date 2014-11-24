
import ceq2_1D as kin
import numpy as np
import pandas as pd
import cantera as ct
from collections import OrderedDict
import unitConversion as uc

if __name__ == "__main__":

    #Define the size of the system

    conv = uc.UnitConverter()
    D = conv.convert_units(1.0, 'in', 'm')
    h = 30.0		#W/m^2-K
    Tw = conv.convert_units(1400, 'C', 'K')
    eps = 0.0
    

    #
    flow = 10/24.0/60.0			#SLPM to mol/s

    z = np.arange(0, 1.0, 0.001)
    ph = ct.Solution('POLIMI_PRF_PAH_HT_1407.cti')
    
    Tw = [1200, 1300, 1400, 1500, 1600]
    P = [101325*5.0]
    Y = 'CH4:%s' % (flow)
    for pressure in P:
        
        for temp in Tw:
            ph.TPX =273.15+25, pressure, Y
            rx = kin.tube = Tube(D = D, h = h, Tw = temp, eps = eps)
            solver = kin.ChemEQ2Solver(ct_phase = ph, tube = rx, mode = 'convective')
            solver.initialize(z, flow)
            solver.conv_eps = 1E-2
            solver.solve(Nc=4)

            outframe = solver.y.copy()
            outframe['Enthalpy'] = np.zeros(len(solver.y['H2']))
            for i in solver.y.index:
                ph.TPX = 273.15+temp, pressure, solver.ct_str(solver.y.loc[i,:])
                outframe.loc[i,'Enthalpy'] = ph.enthalpy_mass

            filename = "methane_crack_hotwall_%s_%s.csv" % (temp, pressure/101325.0)
            outframe.to_csv(filename)
    
    
