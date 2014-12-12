
import ceq2_1D as kin
import numpy as np
import pandas as pd
import cantera as ct
from collections import OrderedDict
import unitConversion as uc
import LabFunctionLib as lfl

if __name__ == "__main__":

    #Define the size of the system

    conv = uc.UnitConverter()
    D = conv.convert_units(1.5, 'in', 'm')
    h = 30.0		#W/m^2-K
    
    eps = 0.0
    cantera_filename = 'POLIMI_PRF_PAH_HT_1407.cti'
    MFC_stdT = [70, "F"]
    MFC_stdP = [101325, "Pa"]

    #Build the streams to be analyzed
    ###1. Methane  ###
    methane_feed = lfl.Stream('CH4_feed', flowrate = (np.ones(1)*5.4, "L/min"), composition = {'CH4':1.0}, basis = "std_gas_volume", cantera_filename = cantera_filename)
    methane_feed.temperature = [25, 'C']
    methane_feed.pressure = [50, 'psig']
    methane_feed.std_temperature = MFC_stdT
    methane_feed.std_pressure = MFC_stdP


    ###2. CO2      ###

    CO2_feed = lfl.Stream('CO2_feed', flowrate = (np.ones(1)*4.5, "L/min"), composition = {'CO2':1.0}, basis = "std_gas_volume", cantera_filename = cantera_filename)
    CO2_feed.temperature = [25, 'C']
    CO2_feed.pressure = [50, 'psig']
    CO2_feed.std_temperature = MFC_stdT
    CO2_feed.std_pressure = MFC_stdP


    ###3. N2       ###

    N2_feed = lfl.Stream('N2_feed', flowrate = (np.ones(1)*47.4, "L/min"), composition = {'N2':1.0}, basis = "std_gas_volume", cantera_filename = cantera_filename)
    N2_feed.temperature = [25, 'C']
    N2_feed.pressure = [50, 'psig']
    N2_feed.std_temperature = MFC_stdT
    N2_feed.std_pressure = MFC_stdP

    ###4. Ar       ###

    Ar_feed = lfl.Stream('Ar_feed', flowrate = (np.ones(1)*2.0, "L/min"), composition = {'Ar':1.0}, basis = "std_gas_volume", cantera_filename = cantera_filename)
    Ar_feed.temperature = [25, 'C']
    Ar_feed.pressure = [50, 'psig']
    Ar_feed.std_temperature = MFC_stdT
    Ar_feed.std_pressure = MFC_stdP


    ###5. H2O      ###

    steam_feed = lfl.Stream('H2O_feed', flowrate = (np.ones(1)*32.0, "g/min"), composition = {'H2O':1.0}, basis = "mass", cantera_filename = cantera_filename)
    steam_feed.temperature = [500, 'C']
    steam_feed.pressure = [50, 'psig']
    steam_feed.std_temperature = MFC_stdT
    steam_feed.std_pressure = MFC_stdP

    #

    inlet_streams = [methane_feed, CO2_feed, N2_feed, Ar_feed, steam_feed]
    mixed_flow = lfl.Mixer('inlet_mix', inlets = inlet_streams)

    #Set up the solver stream
    inlet_temp = conv.convert_units(mixed_flow.outlets[0].temperature[0], mixed_flow.outlets[0].temperature[1], 'K')
    inlet_press = conv.convert_units(mixed_flow.outlets[0].pressure[0], mixed_flow.outlets[0].pressure[1], 'Pa')
    Y = ""
    for specie in mixed_flow.outlets[0].composition:
        if specie == 'Ar':
            s = 'AR'
        else:
            s = specie
        Y += '%s:%s,' % (s, mixed_flow.outlets[0].calcSpeciesMolarFlowrate(specie)[0])
    print Y
    Y = Y[:-1]
    inlet_flow = 0.0
    for s in ['CO2', 'CH4', 'H2O', 'N2', 'Ar']:
         inlet_flow += mixed_flow.outlets[0].calcSpeciesMolarFlowrate(s)[0]


    inlet_flow = inlet_flow/1000.0              #convert to kmol/s from mol/s

    #flow = 10/24.0/60.0/1000			#SLPM to kmol/s

    z = np.arange(0, 24.0/39.37, 0.001)
    ph = ct.Solution('POLIMI_PRF_PAH_HT_1407.cti')
    
    Tw = [1250]
    P = inlet_press
    
    #Y = 'CH4:%s' % (flow)
    for temp in Tw:
        ph.TPX =inlet_temp, inlet_press, Y
        rx = kin.Tube(D = D, h = h, Tw = temp, eps = eps)
        solver = kin.ChemEQ2Solver(ct_phase = ph, tube = rx, mode = 'convective')
        solver.initialize(z, inlet_flow)
        solver.conv_eps = 1E-2
        solver.solve(Nc=4)
        outframe = solver.y.copy()
        outframe['Enthalpy'] = np.zeros(len(solver.y['H2']))
        for i in solver.y.index:
            ph.TPX = inlet_temp, inlet_press, solver.ct_str(solver.y.loc[i,:])
            outframe.loc[i,'Enthalpy'] = ph.enthalpy_mass
        filename = "Longmont_experimental_matrix_point1.csv"
        outframe.to_csv(filename)
    
    
