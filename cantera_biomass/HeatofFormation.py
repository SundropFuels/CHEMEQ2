'''
Created on Mar 11, 2011

@author: JalvingJH
'''

#This script calculates Heats of Formation of the species defined in the .cti file in cal/mol.  It writes the species to a .csv
#It is currently not set to take in arguments

#def heatofFormation(species,inputFile):

from Cantera import *
from Cantera.constants import GasConst_cal_mol_K as R
import numpy as np    
    
Temp = 298.15
f = open('Heat_of_Formation.csv','w')

gas = importPhase('GasifierSpecies.cti','gas')
g_names = gas.speciesNames()
fractions = np.ones(len(g_names))
gas.setTemperature(Temp)
gas.setMoleFractions(fractions)


writeCSV(f, ['Component','Heat of Formation [cal/mol]','Temperature [K]'] )
gas.setTemperature(Temp)
for i in range(len(g_names)):
    species = g_names[i]
    index = gas.speciesIndex(species)
    HoF = (gas.enthalpies_RT()[index]*R*Temp)#/gas.molarMasses()[index]  #cal/mol


    writeCSV(f, [species,HoF,Temp])
    print species, HoF
        
        #return HoF
    
    
