
import ceq2 as kin
import numpy as np
import pandas as pd
import cantera as ct
from collections import OrderedDict

if __name__ == "__main__":
    t = np.arange(0, 600, 0.1)
    
    ph = ct.Solution('completeMechanism.xml')
    T = range(550,1050, 50)
    P = range(101325.0*5, 101325.0*10)
    Y = 'LIGC:0.1,LIGO:0.1,LIGH:0.1,CELL:0.39,ASH:0.01,HCE:0.3,H2O:1.0'
    for pressure in P:
        
        for temp in T:
            ph.TPY =273.15+temp, pressure, Y

            solver = kin.ChemEQ2Solver(ct_phase = ph)
            solver.initialize(t)
            solver.solve(Nc=4)
            #want to add the enthalpies -- this is cludgy -- I should do it as I solve, but that will take longer -- branch and build that later
            outframe = solver.y
            outframe['Enthalpy'] = np.zeros(len(solver.y['H2']))
            for i in solver.y.index:
                ph.TPX = 273.15+temp, pressure, solver.ct_str(solver.y.loc[i,:])
                outframe.loc(i,'Enthalpy') = ph.enthalpy_mole

            filename = "cracking_%s_%1.0.csv" % temp
            solver.y.to_csv(filename)
    
    
