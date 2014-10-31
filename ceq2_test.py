"""Unit tests for the CHEMEQ2 python solver"""

import unittest
import ceq2 as kin
import numpy as np
import pandas as pd
import cantera as ct



#test solver creation
#	-new solver with test phase
#	-initialize solver t - evenly spaced
#	-initialize solver t - unevenly spaced
#	-bad t input - type
#	-bad t input -- not beginning with 0
#	-bad t input -- not monotonically increasing

class SolverCreationTests(unittest.TestCase):
    def testCorrectCreation(self):
        ph = ct.Solution('test.cti')
        ph.TPX = 1000, 101325, 'HE:0.5,AR:0.5'
        solver = kin.ChemEQ2Solver(ct_phase = ph)
        self.assertEqual(solver.ct_phase, ph)

    def testCorrectInitializationEvenSpacing(self):
        t = np.arange(0, 100, 1)
        ph = ct.Solution('test.cti')
        ph.TPX = 1000, 101325, 'HE:0.5,AR:0.5'
        solver = kin.ChemEQ2Solver(ct_phase = ph)
        solver.initialize(t)
        self.assertTrue((solver.t==t).all())
        d = {'AR':np.zeros(len(t)), 'HE':np.zeros(len(t)), 'ARHE':np.zeros(len(t))}
        d['AR'][0] = 0.5
        d['HE'][0] = 0.5
        y = pd.DataFrame(data = d, index = t)
                        
        self.assertTrue((solver.y-y<1E-12).all().all())
       

    def testCorrectInitializationUnevenSpacing(self):
        t = np.array([0,1,5,23])
        ph = ct.Solution('test.cti')
        ph.TPX = 1000, 101325, 'HE:0.5,AR:0.5'
        solver = kin.ChemEQ2Solver(ct_phase = ph)
        solver.initialize(t)
        self.assertTrue((solver.t==t).all())
        d = {'AR':np.zeros(len(t)), 'HE':np.zeros(len(t)), 'ARHE':np.zeros(len(t))}
        d['AR'][0] = 0.5
        d['HE'][0] = 0.5
        y = pd.DataFrame(data = d, index = t)
                        
        self.assertTrue((solver.y-y<1E-12).all().all())
       

    def testBadt_type(self):
        t = 23.0
        ph = ct.Solution('test.cti')
        ph.TPX = 1000, 101325, 'HE:0.5,AR:0.5'
        solver = kin.ChemEQ2Solver(ph)
        
        self.assertRaises(kin.BadInputError, solver.initialize, t)


    def testBadt_nonzeroinitial(self):
        t = np.array([1,34,56])
        ph = ct.Solution('test.cti')
        ph.TPX = 1000, 101325, 'HE:0.5,AR:0.5'
        solver = kin.ChemEQ2Solver(ph)
        
        self.assertRaises(kin.BadInputError, solver.initialize, t)

    def testBadt_nonmonotonic(self):
        t = np.array([1,2,1.5])
        ph = ct.Solution('test.cti')
        ph.TPX = 1000, 101325, 'HE:0.5,AR:0.5'
        solver = kin.ChemEQ2Solver(ph)
        
        self.assertRaises(kin.BadInputError, solver.initialize, t)

#test ys calculation
#	-correct calculation

class ysCalculationTests(unittest.TestCase):
    def testCorrectys(self):
        t = np.arange(0, 100, 1)
        ph = ct.Solution('test.cti')
        ph.TPX = 1000, 101325, 'HE:0.5,AR:0.5'
        solver = kin.ChemEQ2Solver(ct_phase = ph)
        solver.initialize(t)
        #for the purposes of testing, we need to jury-rig the system a little
        solver.dt = 0.1
        y = np.array([0, 0.2, 0.6])
        q = np.array([0.1, 0.0, 0.1])
        p = np.array([0.2, 0.1, 0.0])
        ans = np.array([0,0.20199,0.6])
        self.assertTrue((np.abs(ans-solver.y_pc(y,q,p))<1E-4).all())


#test yp calculation
#	-correct yp

#class ypCalculationTests(unittest.TestCase):
#    def testCorrectyp(self):
#        self.assertTrue(False)

#test yc calculation
#	-correct yc

#class ycCalculationTests(unittest.TestCase):
#    def testCorrectyc(self):
#        self.assertTrue(False)

#THESE ABOVE SHOULD PASS IF THE MAIN SOLVER PASSES -- SKIPPING FOR NOW, AS THEY ARE SIMPLE COMPOSITE FUNCTIONS

#test convergence calculation
#	-correct convergence True
#	-correct convergence False
#	-bad inputs -- kick out of program if this is the case

class ConvergenceCheckTests(unittest.TestCase):
    def testConvergenceTrue(self):
        t = np.arange(0, 100, 1)
        ph = ct.Solution('test.cti')
        ph.TPX = 1000, 101325, 'HE:0.5,AR:0.5'
        solver = kin.ChemEQ2Solver(ct_phase = ph)
        solver.initialize(t)
        solver.yp = np.array([0.10001, 0.200002, 0.300001])
        solver.yc = np.array([0.1, 0.2, 0.3])
        self.assertTrue(solver.converged())
    
    def testConvergenceFalse(self):
        t = np.arange(0, 100, 1)
        ph = ct.Solution('test.cti')
        ph.TPX = 1000, 101325, 'HE:0.5,AR:0.5'
        solver = kin.ChemEQ2Solver(ct_phase = ph)
        solver.initialize(t)
        solver.yp = np.array([0.4, 0.200002, 0.300001])
        solver.yc = np.array([0.1, 0.2, 0.3])
        self.assertTrue(not solver.converged())

    def testBadInputConvergence(self):
        t = np.arange(0, 100, 1)
        ph = ct.Solution('test.cti')
        ph.TPX = 1000, 101325, 'HE:0.5,AR:0.5'
        solver = kin.ChemEQ2Solver(ct_phase = ph)
        solver.initialize(t)
        solver.yp = np.array([0.0, 0.0, 1.0])
        solver.yc = np.array([np.nan, 0.2, 0.3])
        self.assertRaises(kin.NaNError, solver.converged)


#test stability calculation
#	-correct stability

class StabilityCheckTests(unittest.TestCase):
    def testStabilityCorrect(self):
        self.assertTrue(False)

#adjust dt
#	-sigma <=1
#	-sigma >1

class AdjustdtTests(unittest.TestCase):
    def testAdjustdtCorrect_lowsigma(self):
        t = np.arange(0, 100, 1)
        ph = ct.Solution('test.cti')
        ph.TPX = 1000, 101325, 'HE:0.5,AR:0.5'
        solver = kin.ChemEQ2Solver(ct_phase = ph)
        solver.initialize(t)
        solver.dt = 0.1
        solver.sigma = 0.2
        solver.adjust_dt()
        self.assertAlmostEqual(0.224107, solver.dt,4)

    def testAdjustdtCorrect_highsigma(self):
        t = np.arange(0, 100, 1)
        ph = ct.Solution('test.cti')
        ph.TPX = 1000, 101325, 'HE:0.5,AR:0.5'
        solver = kin.ChemEQ2Solver(ct_phase = ph)
        solver.initialize(t)
        solver.dt = 0.1
        solver.sigma = 1.3
        solver.adjust_dt()
        self.assertAlmostEqual(0.088206, solver.dt,4)

#test pade estimation
#	-pdt < 1
#	-pdt = 0
#	-pdt >> 1

class PadeEstimationTests(unittest.TestCase):
    def testPadeEstimation_small(self):
        
        t = np.arange(0, 100, 1)
        ph = ct.Solution('test.cti')
        ph.TPX = 1000, 101325, 'HE:0.5,AR:0.5'
        solver = kin.ChemEQ2Solver(ct_phase = ph)
        solver.initialize(t)
        
        p = np.array([0.911392, 0.508332, 0])
        arg = np.array([10, 0.1, 0])
        self.assertTrue((np.abs(p-solver.pade(arg))<1E-3).all())
        

#correct solution
#	-match analytical solution of "Ar-He" fake problem

class CorrectSolutionTests(unittest.TestCase):
    def testAnalyticalSolution(self):
        t = np.arange(0, 10, 1)
        ph = ct.Solution('test.cti')
        ph.TPX = 300, 101325, 'HE:0.5,AR:0.5'
        solver = kin.ChemEQ2Solver(ct_phase = ph)
        solver.initialize(t)
        solver.solve(Nc=4)
        Ar = np.array([0.020312,0.02001, 0.019433, 0.018627, 0.017651, 0.016566, 0.015428, 0.014283, 0.013166, 0.012102, 0.011104])
        He = np.array([0.020312,0.02001, 0.019433, 0.018627, 0.017651, 0.016566, 0.015428, 0.014283, 0.013166, 0.012102, 0.011104])
        ArHe = np.array([0.0, 0.000302, 0.000879, 0.001685, 0.002661, 0.003746, 0.004884, 0.006029, 0.007146, 0.00821, 0.009208])
        d = {'AR':Ar, 'HE':He, 'ARHE':ArHe}
        y_check = pd.DataFrame(data = d, index = t)
        self.assertTrue((y==y_check).all().all())


if __name__ == "__main__":
    unittest.main()
        
        

