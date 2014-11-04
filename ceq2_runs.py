
import ceq2 as kin
import numpy as np
import pandas as pd
import cantera as ct
from collections import OrderedDict

if __name__ == "__main__":
    t = np.arange(0, 0.1, 0.001)

    ph = ct.Solution('POLIMI_TOT_1407_reduced.cti')
    ph.TPX = 1200, 101325, 'CH4:1.0'

    solver = kin.ChemEQ2Solver(ct_phase = ph)
    solver.initialize(t)
    solver.solve(Nc=4)
        
    t_p = np.arange(0,0.1,0.01)
    print solver.y.loc[t_p, ['CH4','H2','C2H2','C2H4']]
    
