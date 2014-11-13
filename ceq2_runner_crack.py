
import ceq2 as kin
import numpy as np
import pandas as pd
import cantera as ct
from collections import OrderedDict

if __name__ == "__main__":
    t = np.arange(0, 0.5, 0.001)
    
    ph = ct.Solution('POLIMI_PRF_PAH_HT_1407.cti')
    #T = range(1300,1500, 100)
    #P = [101325.0*5, 101325.0*10]
    T = [1500,1600]
    P = [101325*2.]
    Y = 'CH4:0.5, H2:0.5'
    for pressure in P:
        
        for temp in T:
            ph.TPX =273.15+temp, pressure, Y

            solver = kin.ChemEQ2Solver(ct_phase = ph)
            solver.initialize(t)
            solver.solve(Nc=3)
            #want to add the enthalpies -- this is cludgy -- I should do it as I solve, but that will take longer -- branch and build that later
            outframe = solver.y.copy()
            outframe['Enthalpy'] = np.zeros(len(solver.y['H2']))
            for i in solver.y.index:
                ph.TPX = 273.15+temp, pressure, solver.ct_str(solver.y.loc[i,:])
                outframe.loc[i,'Enthalpy'] = ph.enthalpy_mass

            filename = "methane_crack_%s_%s.csv" % (temp, pressure/101325.0)
            outframe.to_csv(filename)
    
    
