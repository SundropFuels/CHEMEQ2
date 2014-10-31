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
        self.assertRaises(kin.BadInputError, kin.ChemEQ2Solver, ph)


    def testBadt_nonzeroinitial(self):
        t = np.array([1,34,56])
        ph = ct.Solution('test.cti')
        ph.TPX = 1000, 101325, 'HE:0.5,AR:0.5'
        self.assertRaises(kin.BadInputError, kin.ChemEQ2Solver, ph)

    def testBadt_nonmonotonic(self):
        t = np.array([1,2,1.5])
        ph = ct.Solution('test.cti')
        ph.TPX = 1000, 101325, 'HE:0.5,AR:0.5'
        self.assertRaises(kin.BadInputError, kin.ChemEQ2Solver, ph)

#test ys calculation
#	-correct calculation
#	-bad inputs

class ysCalculationTests(unittest.TestCase):
    def testCorrectys(self):
        self.assertTrue(False)

    def testBadInputsys(self):
        self.assertTrue(False)

#test yp calculation
#	-correct yp
#	-bad inputs

class ypCalculationTests(unittest.TestCase):
    def testCorrectyp(self):
        self.assertTrue(False)

    def testBadInputsyp(self):
        self.assertTrue(False)

#test yc calculation
#	-correct yc
#	-bad inputs

class ycCalculationTests(unittest.TestCase):
    def testCorrectyc(self):
        self.assertTrue(False)

    def testBadInputsyc(self):
        self.assertTrue(False)

#test convergence calculation
#	-correct convergence True
#	-correct convergence False
#	-bad inputs

class ConvergenceCheckTests(unittest.TestCase):
    def testConvergenceTrue(self):
        self.assertTrue(False)
    
    def testConvergenceFalse(self):
        self.assertTrue(False)

    def testBadInputConvergence(self):
        self.assertTrue(False)

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
        self.assertTrue(False)

    def testAdjustdtCorrect_highsigma(self):
        self.assertTrue(False)

#test pade estimation
#	-pdt < 1
#	-pdt = 0
#	-pdt >> 1

class PadeEstimationTests(unittest.TestCase):
    def testPadeEstimation_small(self):
        self.assertTrue(False)
        #ph = ct.Solution('test.cti')
        #ph.TPX = 1000, 101325, 'He:0.5,Ar:0.5'
        #solver = kin.ChemEQ2Solver(ct_phase = ph)

    def testPadeEstimation_zero(self):
        self.assertTrue(False)

    def testPadeEstimation_big(self):
        self.assertTrue(False)

#correct solution
#	-match analytical solution of "Ar-He" fake problem

class CorrectSolutionTests(unittest.TestCase):
    def testAnalyticalSolution(self):
        self.assertTrue(False)


if __name__ == "__main__":
    unittest.main()
        
        

