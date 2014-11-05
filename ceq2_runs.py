
import ceq2 as kin
import numpy as np
import pandas as pd
import cantera as ct
from collections import OrderedDict

if __name__ == "__main__":
    t = np.arange(0, 0.1, 0.001)

    ph = ct.Solution('POLIMI_TOT_1407.cti')
    temps = range(1000,1450,50)
    pressures = [101325, 101325*2, 101325*3]
    for T in temps:
        for P in pressures:
            ph.TPX = T+273.15, P, 'CH4:1.0'

            solver = kin.ChemEQ2Solver(ct_phase = ph)
            solver.initialize(t)
            solver.solve(Nc=4)
            filename = 'CH4_crack_no_oxidant_%s_%s' % (T, P/101325.0)
            solver.y.to_csv(filename)
    

    
