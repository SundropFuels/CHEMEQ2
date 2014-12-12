import unitConversion as uc
import numpy as np


class ThermoError(Exception):
    def __init__(self,value):
        self.value = value

    def __str__(self):
        return repr(self.value)

class ThermoGenerator:
    """Generates thermodynamic data for various species"""

    NASA_poly = {}

    NASA_poly['H2_high'] = [3.33727920E+00,-4.94024731E-05,4.99456778E-07, -1.79566394E-10, 2.00255376E-14, -9.50158922E+02,-3.20502331E+00]   
    NASA_poly['H2_low'] = [2.34433112E+00, 7.98052075E-03,-1.94781510E-05, 2.01572094E-08, -7.37611761E-12, -9.17935173E+02, 6.83010238E-01]

    NASA_poly['H2O_high'] = [3.03399249E+00, 2.17691804E-03, -1.64072518E-07, -9.70419870E-11, 1.68200992E-14, -3.00042971E+04, 4.96677010E+00]
    NASA_poly['H2O_low'] = [4.19864056E+00, -2.03643410E-03, 6.52040211E-06,-5.48797062E-09, 1.77197817E-12, -3.02937267E+04, -8.49032208E-01]

    NASA_poly['CO_high'] = [2.71518561E+00, 2.06252743E-03,-9.98825771E-07, 2.30053008E-10,-2.03647716E-14,-1.41518724E+04, 7.81868772E+00]
    NASA_poly['CO_low'] = [3.57953347E+00, -6.10353680E-04,  1.01681433E-06, 9.07005884E-10,-9.04424499E-13,-1.43440860E+04, 3.50840928E+00]
    
    NASA_poly['CO2_high'] = [3.85746092E+00, 4.1437026E-3, 2.21481404E-6, 5.23490188E-10, -4.72084164E-14, -4.87591660E+04, 2.27163806E+00]
    NASA_poly['CO2_low'] = [2.35677352E+00, 8.98459677E-03, -7.12356269E-6, 2.45919022E-9, -1.43699548E-13, -4.83719697E4, 9.90105222E00]

    NASA_poly['CH4_high'] = [7.48514950E-2, 1.33909467E-02, -5.73285809E-06, 1.22292535E-9, -1.01815230E-13, -9.46834459E3, 1.84373180E1]
    NASA_poly['CH4_low'] = [5.14987613E00, -1.36709788E-2, 4.91800599E-05, -4.84743026E-08, 1.66693956E-11, -1.02466476E4, -4.64130376E0]

    NASA_poly['Ar_high'] = [0.02500000E2, 0.00000000E0, 0.00000000E0, 0.00000000E0, 0.000000000E0, -0.0745375E4, 0.04366000E2]
    NASA_poly['Ar_low'] = [0.02500000E2, 0.00, 0.00, 0.00, 0.00, -0.0745375E4, 0.04366000E2]

    NASA_poly['N2_high'] = [0.0292664E2, 0.14879768E-2, -0.0568476E-5, 0.10097938E-9, -0.06753351E-13, -0.09227977E4, 0.05980528E2]
    NASA_poly['N2_low'] = [0.03298677E2, 0.14082404E-2, -0.03963222E-4, 0.05641515E-7, -0.02444854E-10, -0.10208999E4, 0.03950372E2]

    NASA_poly['C2H2_high'] = [4.14756964E0, 5.96166664E-3, -2.37294852E-6, 4.6712171E-10, 3.61235213E-14, 2.59359992E4, -1.23028121E0]
    NASA_poly['C2H2_low'] = [8.08681094E-1, 2.33615629E-2, -3.55171815E-5, 2.80152437E-8, -8.50072974E-12, 2.6428987E4, 1.39397051E1]

    NASA_poly['C2H4_high'] = [2.03611116E0, 1.46454151E-2, -6.71077915E-6, 1.47222923E-9, -1.25706061E-13, 4.93988614E3, 1.03053693E1]
    NASA_poly['C2H4_low'] = [3.95920138E0, -7.57052247E-3, 5.70990292E-5, -6.91588753E-8, 2.69884373E-11, 5.08977593E3, 4.09733096E0]

    NASA_poly['C2H6_high'] = [1.07188150E0, 2.16852677E-2, -1.00256067E-5,2.21412001E-9, -1.90002890E-13, -1.14263932E4, 1.51156107E1]
    NASA_poly['C2H6_low'] = [4.29142492E0, -5.50154270E-3, 5.99438288E-5, -7.08466285E-8, 2.68685771E-11, -1.15222055E4, 2.66682316E0]

    NASA_poly['C3H8_high'] = [0.75341368E+01, 0.18872239E-01,-0.62718491E-05, 0.91475649E-09,-0.47838069E-13, -0.16467516E+05, -0.17892349E+02]
    NASA_poly['C3H8_low'] = [.93355381E0, 0.26424579E-01, 0.61059727E-05,-0.21977499E-07, 0.95149253E-11,-0.13958520E+05, 0.19201691E+02]
    
    NASA_poly['C3H6_high'] = [6.038704990E+00, 1.629638950E-02, -5.821306240E-06, 9.359364830E-10, -5.586029030E-14, -7.765950920E+02, -8.438243220E+00]
    NASA_poly['C3H6_low'] = [3.834645240E+00, 3.290784050E-03, 5.052281840E-05, -6.662514180E-08, 2.637075850E-11, 7.538382950E+02, 7.534109950E+00]

    NASA_poly['C6H6_high'] = [1.107717080E+01, 2.070678950E-02, -7.516251000E-06, 1.222094160E-09, -7.353125130E-14, 4.309883950E+03, -4.001169500E+01]
    NASA_poly['C6H6_low'] = [5.034696640E-01, 1.851423630E-02, 7.378644090E-05, -1.181061270E-07,5.071825270E-11, 8.552662930E+03, 2.164817960E+01]

    NASA_poly['C7H8_high'] = [1.293947500E+01, 2.669215580E-02, -9.684201080E-06,   1.573921400E-09,  -9.466704820E-14, -6.770357690E+02,  -4.672553020E+01]
    NASA_poly['C7H8_low'] = [1.611914000E+00, 2.111889020E-02, 8.532214530E-05, -1.325668760E-07, 5.594061090E-11, 4.096519760E+03,   2.029736140E+01]

    NASA_poly['C10H8_high'] = [1.861298990E+01, 3.044941410E-02, -1.112247990E-05, 1.816154060E-09, -1.096012240E-13, 8.915529440E+03, -8.002304790E+01]
    NASA_poly['C10H8_low'] = [-1.049193260E+00,  4.629706110E-02, 7.075922030E-05, -1.384081860E-07, 6.204757480E-11, 1.598463880E+04, 3.021215710E+01]
    
    NASA_poly['cellulose_high'] = []
    NASA_poly['cellulose_low'] = []

    NASA_poly['hemicellulose_high'] = []
    NASA_poly['hemicellulose_low'] = []

    NASA_poly['lig_C_high'] = []
    NASA_poly['lig_C_low'] = []

    NASA_poly['lig_O_high'] = []
    NASA_poly['lig_O_low'] = []

    NASA_poly['lig_H_high'] = []
    NASA_poly['lig_H_low'] = []

    R = (8.314, "J/mol/K")

    def __init__(self):
        """Creates a thermodynamic property generator"""
        self.unit_converter = uc.UnitConverter()

    def calc_heat_capacity(self, species, T, units):
        """Returns the heat capacity at temperature T (value, units) in the given units"""
        #This will always work in J/mol-K as the baseline, and then make the conversion later
        if '%s_high' % species not in ThermoGenerator.NASA_poly.keys() or '%s_low' % species not in ThermoGenerator.NASA_poly.keys():
            raise ThermoError, "The given species is not in the NASA polynomial database"

        temp = self.unit_converter.convert_units(T[0], T[1], 'K')
        if temp < 1000:
            qual = 'low'
        else:
            qual = 'high'
        polys = ThermoGenerator.NASA_poly['%s_%s' % (species, qual)]
        Cp_SI = ThermoGenerator.R[0] * (polys[0] + polys[1] * temp + polys[2] * np.power(temp,2) + polys[3] * np.power(temp,3) + polys[4] * np.power(temp,4))
        return self.unit_converter.convert_units(Cp_SI, ThermoGenerator.R[1], units)

    def calc_enthalpy(self, species, T, units):
        """Returns the specific enthalpy at temperature T in the given units"""
        if '%s_high' % species not in ThermoGenerator.NASA_poly.keys() or '%s_low' % species not in ThermoGenerator.NASA_poly.keys():
            raise ThermoError, "The given species is not in the NASA polynomial database"

        temp = self.unit_converter.convert_units(T[0], T[1], 'K')
        if temp < 1000:
            qual = 'low'
        else:
            qual = 'high'
        polys = ThermoGenerator.NASA_poly['%s_%s' % (species, qual)]
        H_SI = ThermoGenerator.R[0] * temp*(polys[0] + polys[1] * temp/2.0 + polys[2] * np.power(temp,2)/3.0 + polys[3] * np.power(temp,3)/4.0 + polys[4] * np.power(temp,4)/5.0+polys[5]/temp)
        return self.unit_converter.convert_units(H_SI, '%s*K' % ThermoGenerator.R[1], units)

    def calc_entropy(self, species, T, units):
        """Returns the specific entropy at temperature T in the given units"""
        if '%s_high' % species not in ThermoGenerator.NASA_poly.keys() or '%s_low' % species not in ThermoGenerator.NASA_poly.keys():
            raise ThermoError, "The given species is not in the NASA polynomial database"

        temp = self.unit_converter.convert_units(T[0], T[1], 'K')
        if temp < 1000:
            qual = 'low'
        else:
            qual = 'high'
        polys = ThermoGenerator.NASA_poly['%s_%s' % (species, qual)]
        S_SI = ThermoGenerator.R[0] * (polys[0]*np.log(temp) + polys[1] * temp + polys[2] * np.power(temp,2)/2.0 + polys[3] * np.power(temp,3)/3.0 + polys[4] * np.power(temp,4)/4.0+polys[6])
        return self.unit_converter.convert_units(S_SI, ThermoGenerator.R[1], units)

if __name__ == "__main__":
    Therm = ThermoGenerator()
    print Therm.calc_enthalpy('CO', (1500, 'K'), 'kJ/mol')

