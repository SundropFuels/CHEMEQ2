"""Unit tests for the CHEMEQ2 python solver"""

import unittest
import ceq2 as kin
import cantera as ct

#test solver creation
#	-new solver with test phase
#	-initialize solver t
#	-initialize solver y

class SolverCreationTests(unittest.TestCase):
    def testCorrectCreation(self):
        self.assertTrue(False)

    def testCorrectInitialization(self):
        self.assertTrue(False)

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



        
        

