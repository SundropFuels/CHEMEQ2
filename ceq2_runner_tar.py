
import ceq2 as kin
import numpy as np
import pandas as pd
import cantera as ct
from collections import OrderedDict

if __name__ == "__main__":
    t = np.arange(0, 4.0, 0.01)
    
    ph = ct.Solution('POLIMI_TOT_1407.cti')
    #T = range(1300,1500, 100)
    #P = [101325.0*5, 101325.0*10]
    T = [1300,1200,1400]
    P = [101325*5.0]
    Y = 'H2:0.132,CO2:0.0708,CO:0.233,H2O:0.30,CH4:0.0850,C2H2:0.0023,C2H4:0.0241,C2H6:0.0016,C6H6:0.0007,C10H8:0.0013'
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

            filename = "tarcrack_highT_%s_%s_silva.csv" % (temp, pressure/101325.0)
            outframe.to_csv(filename)
    
    
