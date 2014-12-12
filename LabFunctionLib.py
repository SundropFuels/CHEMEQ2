"""LabFunctionLib.py
   General data processing and handling functions for experimental analysis
   Chris Perkins
   2012-01-23

   Version history:
   Version 0.1 - First porting over from R libraries

"""

import numpy as np
import db_toolbox as SQL
import dataFrame_pd as df
import datetime
import unitConversion as uc
import Thermo
import matplotlib.pyplot as plt
import collections
import element_parser as ep
import scipy.optimize as spo
try:
    import Cantera as ct
except ImportError:
    import cantera as ct
import distutils.version as vch
#we need to set a global version parameter for the cantera version
if vch.LooseVersion(ct.__version__) < vch.LooseVersion("2.1.1"):
    ct_api = "old"
else:
    ct_api = "new"

import time
import pandas as pd
import run_eqiv as eqv
import scipy.stats as st


conv = uc.UnitConverter()

class lflException(df.dfException):
    pass

class SQLInterfaceError(lflException):
    pass

class TimeError(lflException):
    pass

class InterpArgumentError(lflException):
    pass

class SQLInterfaceError(lflException):
    pass

class UnitConversionError(lflException):
    pass

class NoUnitsError(UnitConversionError):
    pass

class SpeciesNotDefinedError(lflException):
    pass

class BadStreamError(lflException):
    pass

class ConversionError(lflException):
    pass

class NoInletOutletFlowrateError(lflException):
    pass

class BadCTTransInputError(lflException):
    pass

class ts_data(df.Dataframe):
    """General timeseries data class"""
    def __init__(self, start, end, data = None, units_dict = None):
        df.Dataframe.__init__(self, data = data, units_dict = units_dict)
        if type(start) !=  datetime.datetime or type(end) !=  datetime.datetime:
            raise TimeError, "The start and end times for a timeseries must be datetime.datetime objects"
        if start > end:
            raise TimeError, "The start time should happen before the end time"

        self.start = start
        self.end = end
        self.uncertainty_dict = {}  #The uncertainty dictionary -- need to add a function to fill this from an SQL table

        self.avgs = {}
        self.stdevs = {}
   

    def reinitialize_data(self, data):
        self.__init__(self.start, self.end, data = data, units_dict = self.units)

    def SQL_load(self, db_interface, table = "", glossary = 'tag_glossary_tbl', timeseries_col = 'ts'):

        if not isinstance(db_interface, SQL.db_interface):
            raise SQLInterfaceError, "The passed interface to the database is not valid"

        start_condition = "%s >=  '%s'" % (timeseries_col, self.start.strftime("%Y-%m-%d %H:%M:%S"))
        end_condition = "%s <=  '%s'" % (timeseries_col, self.end.strftime("%Y-%m-%d %H:%M:%S"))
        try:
            self.SQL_load_data(db_interface, table, conditions = [start_condition, end_condition])
        except SQL.NoDBInUseError:
            raise SQLInterfaceError, "No database selected for the passed interface"
        except SQL.InterfaceNotConnectedError:
            raise SQLInterfaceError, "Passed interface is not connected to a server"
        counter = np.arange(len(self))
        self['counter'] = counter

        for i in self.columns:
            try:
                q=SQL.select_Query(objects=['units'], table=glossary, condition_list=["simple_name='%s'" % i])
                
                self.units[i]=db_interface.query(q)[0]['units']
            except IndexError:
                self.units[i]=None

        #Need to add the units load to this and to the unittests

    def SQL_db_upload(self, db_interface, table = ""):
        if not isinstance(db_interface, SQL.db_interface):
            raise SQLInterfaceError, "The passed interface to the database is not valid"

        try:
            self.SQL_upload_data(db_interface, table)
        except SQL.NoDBInUseError:
            raise SQLInterfaceError, "No database selected for the passed interface"
        except SQL.InterfaceNotConnectedError:
            raise SQLInterfaceError, "Passed interface is not connected to a server"
        except df.dfException, e:
            raise SQLInterfaceError, "%s" % e

    def SQL_uc_load(self, db_interface, table = ""):
        pass

    def interpolate_col(self, x, y):
        """Interpolates between None or np.nan values on a given column ('y') relative to 'x' and adds the interpolated column as 'Name'_interp"""
        if x not in self.columns or y not in self.columns:
            raise InterpArgumentError, "The given column name is not in the data frame -- cannot interpolate"

        interp_col = self[y].copy()
        #Scroll forward to find the first non-nan/non-None value
        i = 0
        try:
            
            while i < len(self[y]) and (self[y][i] ==  None or np.isnan(self[y][i])):
                i+= 1
            if i ==  len(self[y]):
                #The column was empty -- nothing to do
                return interp_col
            begin = i
            while i < len(self[y]):
                i +=  1
                if i ==  len(self[y]):
                    #Trailing nones/nans
                    return interp_col
                elif self[y][i] ==  None or np.isnan(self[y][i]):
                    continue
                elif i ==  len(self[y]):
                    #Trailing nones/nans
                    return interp_col
                elif i - begin ==  1:
                    #There were no nans or Nones in the mix, move forward
                    begin = i
                else:
                    #Now interpolate between the values
                    end = i
                    for j in range(begin+1,end):
                        interp_col[j] = (self[y][end] - self[y][begin])/(self[x][end] - self[x][begin])*(self[x][j] - self[x][begin]) + self[y][begin]
                    begin = i
        except IndexError:
            return interp_col
        return interp_col

    def get_ms_from_csv(self, csvfile):
        ppmlist = ['Benzene_MS', 'Hydrogen Sulfide_MS', 'Napthalene_MS', 'Toluene_MS']
        for i in self.data:
            if i.endswith('_MS'):
                for j in range(len(self[i])): #Clear list
                    self[i][j] = None
                with open(csvfile) as f:
                    colnames = f.readline().split(',')
                    timecol = colnames.index('Time&Date')
                    icol = colnames.index(i.replace('_MS', ''))
                    for line in f:
                        row = line.split(',')
                        time = datetime.datetime.strptime(row[timecol],'%Y/%m/%d %H:%M:%S')
                        if i in ppmlist:
                            conc = float(row[icol])/10000
                        else:
                            conc = float(row[icol])
                        self[i][np.where(self['ts'] == time)[0]] = conc
                    
    def calc_rel_time(self):
        rel_time = []
        for i in self['ts']:
            rel_time.append((i-self['ts'][0]).seconds)
        rel_time = np.array(rel_time)
        return rel_time

    def rid_nones(self, colname):
        for i in range(len(colname)-1):
            if colname[i-1] == None and colname[i]!= None:
                firstval = i
            if colname[i+1] == None and colname[i]!= None:
                lastval = i
        colname[0:firstval] = colname[firstval]
        colname[lastval+1:] = colname[lastval]
        return(colname)

    def generate_averages_and_stdevs(self, cols=None):
        """Generates numbers for the averages and standard deviations of the indicated columns, and stores them in member dicts"""
        if cols is None:
            cols = self.data.keys()

        for key in cols:
            try:

                #by default, I ignore nan and inf values
                
                self.avgs[key] = self[key][np.isfinite(self[key].astype(float))].mean(dtype='float64')
                self.stdevs[key] = self[key][np.isfinite(self[key].astype(float))].std(dtype='float64')

            except KeyError:
                raise lflException, "%s is not a key in the dataframe" % key
            except ZeroDivisionError:
                self.avgs[key] = np.nan
                self.stdevs[key] = np.nan
	
    def interp_ga(self):
        self['rel_time'] = self.calc_rel_time()
        self.units['rel_time']='s'
        for i in self.data:
            if i.endswith('_MS') or i.endswith('_GC'):

                colname = i
                try:
                    interpcol = self.interpolate_col('rel_time', colname)
                    interpcol = self.rid_nones(interpcol)
                    self[colname] = interpcol
                except IndexError: pass
                if self.units[i]=='ppm':
                    self[i]=self[i]/10000
                    self.units[i]='%'

    def calc_product_flowrates(self):
        """Calculate product gas flow rates in moles/s"""
        gaslist = []
        for i in self.data:
            if i.endswith('_MS'):
                gaslist.append(i)
        for i in gaslist:
            product_flowrate = self.data[i]/100*self.data['outlet_flowrate']
            self.data[i+'_flowrate'] = product_flowrate

    def calc_inst_conversion(self):
        """Calculates instantaneous Carbon conversion"""
        self['carbon_conversion'] = self['carbon_out_total']/self['carbon_in']

    def calc_average_feed_rate(self):
        return np.average(self.data['mass_flow_brush_feeder'])

    def calc_std_feed_rate(self):
        return np.std(self.data['mass_flow_brush_feeder'])
        
    def calc_total_c_in(self):
        return np.sum(self.data['carbon_in'])

    def calc_total_c_out(self):
        return np.sum(self.data['carbon_out_total'])

    def calc_carbon_x(self):
        return np.sum(self.data['carbon_out_total'])/np.sum(self.data['carbon_in'])

    def convert_col_units(self, column, units):
        """Will convert the given column into the new units from the existing units"""
        
        if column not in self.columns:
            raise df.NoColumnError, "The requested column is not in the data frame."

        if self.units ==  {}:
            raise NoUnitsError, "The units have not been defined for this column"

        conv = uc.UnitConverter()

        try:
            self[column] = conv.convert_units(self[column], self.units[column], units)
        except uc.UnitNotFoundError:
            raise UnitConversionError, "The unit was not in the database"
        except uc.InconsistentUnitError:
            raise UnitConversionError, "The units are not dimensionally consistent with the current units"

        self.units[column] = units    

    def get_val(self, col_name):
        if col_name in self.units.keys():
            return (self[col_name], self.units[col_name])
        else:
            return (self[col_name], "")

    def list_keys(self):
        """Generates list of column names in their respective categories"""
        keylist = collections.OrderedDict({})
        keylist['MFCs'] = []
        keylist['Temperatures'] = []
        keylist['Pressures'] = []
        keylist['Setpoints'] = []
        keylist['Mass Spec'] = []
        keylist['GC'] = []
        keylist['NDIR'] = []
        keylist['Product Flow Rates'] = []
        keylist['Carbon Flows'] = []
        
        for i in self.data:
            if i.endswith('_MS'):
                keylist['Mass Spec'].append(i)
            if i.endswith('_GC'):
                keylist['GC'].append(i)
            if i.startswith('temp_'):
                keylist['Temperatures'].append(i)
            if i.startswith('pressure_'):
                keylist['Pressures'].append(i)
            if i.startswith('mass_'):
                keylist['MFCs'].append(i)
            if i.endswith('_flowrate'):
                keylist['Product Flow Rates'].append(i)
            if i.startswith('carbon_') or i.endswith('_carbon_out'):
                keylist['Carbon Flows'].append(i)
            if i.startswith('setpoint_'):
                keylist['Setpoints'].append(i)
            if i.startswith('NDIR_'):
                keylist['NDIR'].append(i)
        for k in keylist:
            keylist[k].sort()
        return keylist

class Stream:
    species = {}
    species['CO'] = {'C':1,'O':1, 'name':'Carbon Monoxide', 'phase':'g'}
    species['CO2'] = {'C':1,'O':2, 'name':'Carbon Dioxide', 'phase':'g'}
    species['CH4'] = {'C':1, 'H':4, 'name':'Methane','phase':'g'}
    species['H2'] = {'H':2, 'name':'Hydrogen','phase':'g'}
    species['H2S'] = {'H':2, 'S':1, 'name':'Hydrogen Sulfide','phase':'g'}
    species['C2H2'] = {'C':2,'H':2, 'name':'Acetylene','phase':'g'}
    species['C2H4'] = {'C':2,'H':4, 'name':'Ethylene','phase':'g'}
    species['C2H6'] = {'C':2,'H':6, 'name':'Ethane','phase':'g'}
    species['C3H8'] = {'C':3, 'H':8, 'name':'Propane','phase':'g'}
    species['C3H6'] = {'C':3, 'H':6, 'name':'Propylene','phase':'g'}
    species['C6H6'] = {'C':6, 'H':6, 'name':'Benzene','phase':'g'}
    species['C7H8'] = {'C':7, 'H':8, 'name':'Toluene','phase':'g'}
    species['C10H8'] = {'C':10, 'H':8, 'name':'Napthalene', 'phase':'g'}
    species['C4H8'] = {'C':4, 'H':8, 'name':'1-Butene', 'phase':'g'}
    species['C4H10'] = {'C':4, 'H':10, 'name':'n-Butane', 'phase':'g'}
    species['CH3CHCH3CH3'] = {'C':4, 'H':10, 'name':'i-Butane','phase':'g'}
    species['C6H4CH3CH3'] = {'C':8, 'H':10, 'name':'o-xylene', 'phase':'g'}
    species['C6H5CH2CH3'] = {'C':8, 'H':10, 'name':'ethyl benzene', 'phase':'g'}
    species['H2O'] = {'H':2,'O':1, 'name':'Water', 'phase':'g'}
    species['Ar'] = {'Ar':1, 'name':'Argon', 'phase':'g'}
    species['N2'] = {'N':2, 'name':'Nitrogen', 'phase':'g'}
    species['CELL'] = {'C':6, 'H':10, 'O':5, 'name':'Cellulose', 'phase':'s'}
    species['HCE'] = {'C':5, 'H':8, 'O':4, 'name':'Hemicellulose', 'phase':'s'}
    species['LIGH'] = {'C':22, 'H':28, 'O':9, 'name':'Lig-H', 'phase':'s'}
    species['LIGO'] = {'C':20, 'H':22, 'O':10, 'name':'Lig-O', 'phase':'s'}
    species['LIGC'] = {'C':15, 'H':14, 'O':4, 'name':'Lig-C', 'phase':'s'}
    species['O2'] = {'O':2, 'name':'Oxygen','phase':'g'}
    
    ct_trans = {}
    ct_trans['Ar'] = {'AR':1.0}
    ct_trans['biomass'] = {'CELL':0.3932, 'LIGC':0.10, 'LIGH':0.10, 'LIGO':0.10, 'HCE':0.30}   #Still missing ASH - will correct in a bit

    names = {}
    for i in species.keys():       
        names[species[i]['name']] = i

    def __init__(self, name, flowrate = None, composition = None,
                 basis = "molar", temperature = None, pressure = None,
                 density = None, compressible = None, std_temperature = (25.0, 'C'), std_pressure = (101325.0, 'Pa'), mode = None, cantera_filename = 'cantera_biomass/GasifierSpecies.cti'):

        if composition is None:
            composition = {}
        self.name = name
        self.mode = mode
        self.length = None
        self.temperature = None
        self.pressure = None
        self.composition = None
        #need to average temperature, pressure, and composition if they are arrays
        self.set_temperature(temperature)
        self.set_pressure(pressure)
        self.set_composition(composition)
        self.basis = basis
               
        self.density = density #(value, units)
        self.flowrate = flowrate #(value, units)
        self.compressible = compressible
        self.std_temperature = std_temperature
        self.std_pressure = std_pressure

        self.special_species = {}
        self.enthalpy_reserve = [0.0, 'W']  #This is a cludge fix -- it will be used to handle liquid water until I make streams n-phase capable

        self.cantera_helper = CanteraHelper()
        self.ctphase = self.cantera_helper.importPhase(cantera_filename, 'gas')
        #self.ctphase = ct.importPhase('cantera_biomass/GasifierSpecies.cti','gas')

    def convert_scalar_to_vector(self, length):
        if self.mode == "scalar":
            self.mode = "vector"
            if self.temperature is not None and not isinstance(self.temperature[0], pd.Series) and not isinstance(self.temperature[0],np.ndarray):
                self.temperature[0] = pd.Series(np.ones(length)*self.temperature[0])
            if self.pressure is not None and not isinstance(self.pressure[0], pd.Series) and not isinstance(self.pressure[0],np.ndarray):
                self.pressure[0] = pd.Series(np.ones(length)*self.pressure[0])
            if self.composition is not None:

                comp = {}
                for k in self.composition:
                    if not isinstance(self.composition[k], pd.Series) and not isinstance(self.composition[k],np.ndarray):
                        self.composition[k] = pd.Series(np.ones(length)*self.composition[k])
        
    def set_temperature(self, temperature):
        
        if temperature is not None:
            if isinstance(temperature[0], pd.Series) or isinstance(temperature[0], np.ndarray):
                if self.mode is None or self.mode == "vector":
                    if self.length is None or len(temperature[0]) == self.length:
                        self.mode = "vector"
                        self.length = len(temperature[0])
                    else:
                        raise Exception, "Length of new temperature is not consistent with length of vector stream"
                
                elif self.mode == "scalar":
                    self.length = len(temperature[0])
                    self.convert_scalar_to_vector(self.length)
            
            else:
                if self.mode is None or self.mode == "scalar":
                    self.mode = "scalar"
                elif self.mode == "vector":
                    temperature[0] = pd.Series(temperature[0]*np.ones(self.length))

            self.temperature = [temperature[0], temperature[1]]
        else:
            self.temperature = None

    def set_pressure(self, pressure):
        if pressure is not None:
            if isinstance(pressure[0], pd.Series) or isinstance(pressure[0], np.ndarray):
                if self.mode is None or self.mode == "vector":
                    if self.length is None or len(pressure[0]) == self.length:
                        self.mode = "vector"
                        self.length = len(pressure[0])
                    else:
                        raise Exception, "Length of new pressure vector is not consistent with length of vector stream"
                elif self.mode == "scalar":
                    self.length = len(pressure[0])
                    self.convert_scalar_to_vector(self.length)
            
            else:
                if self.mode is None or self.mode == "scalar":
                    self.mode = "scalar"
                elif self.mode == "vector":
                    pressure[0] = pd.Series(pressure[0]*np.ones(self.length))
                
            
            self.pressure = [pressure[0], pressure[1]]
        else:
            self.pressure = None

    def set_composition(self, composition):
        """Sets the composition"""
        if composition == {}:
            self.composition = {}
        else:

            #dict checking
            if not isinstance(composition, dict):
                raise Exception, "The composition must be provided as a {species:fraction} dictionary"

            #mode checking
            if isinstance(composition[composition.keys()[0]], pd.Series) or isinstance(composition[composition.keys()[0]], np.ndarray):
                if self.mode is None or self.mode == "vector":
                    self.mode = "vector"
                elif self.mode == "scalar":
                    self.length = len(composition[composition.keys()[0]])
                    self.convert_scalar_to_vector(self.length)

            else:
                if self.mode == "vector":
                    #convert the composition to a vector to allow for simple scalar assignment
                    for k in composition:
                        composition[k] = pd.Series(composition[k]*np.ones(self.length))
                elif self.mode is None or self.mode == "scalar":
                    self.mode = "scalar"
                else:
                    raise Exception, "Trying to set a scalar composition to a vector valued stream"
            
            getattr(self, "_set_composition_%s" % self.mode)(composition)

    def _set_composition_scalar(self, composition):
        """Internal scalar function for setting the composition"""
        self.composition = composition

    def _set_composition_vector(self, composition):
        """Internal vector function for setting the composition"""
        #len checking first
        l = len(composition[composition.keys()[0]])
        s = np.zeros(l)
        
        for k in composition:
            if len(composition[k])!=l:
                raise Exception, "Not all of the fractional composition vectors are of the same length"
            s += composition[k]
        
        if not (s==np.ones(l)).all():
            pass
            #raise Exception, "The compositions do not all add up to 1.0 in the composition setup"
            #!!!#FIX
        self.composition = composition

    def gas_volumetric_flowrate(self, units):
        """Returns the gas volumetric flowrate, in the desired units"""
        conv = uc.UnitConverter()
        if self.basis == "gas_volume":
            return conv.convert_units(self.flowrate[0], self.flowrate[1], units)
        elif self.basis == "std_gas_volume":
            T = conv.convert_units(self.temperature[0], self.temperature[1], 'K')
            p = conv.convert_units(self.pressure[0], self.pressure[1], 'Pa')
            std_T = conv.convert_units(self.std_temperature[0], self.std_temperature[1], 'K')
            std_P = conv.convert_units(self.std_pressure[0], self.std_pressure[1], 'Pa')
            return conv.convert_units(self.flowrate[0], self.flowrate[1], units)*T/std_T*std_P/p
        elif self.basis == "molar":
            f = conv.convert_units(self.flowrate[0], self.flowrate[1], 'mol/s')
            T = conv.convert_units(self.temperature[0], self.temperature[1], 'K')
            p = conv.convert_units(self.pressure[0], self.pressure[1], 'Pa')
            return conv.convert_units(8.314*f*T/p, 'm^3/s', units)
        elif self.basis == "mass":
            #convert to a molar flow first ###!!!###
            val = conv.convert_units(self.flowrate[0], self.flowrate[1], 'g/s')
            T = conv.convert_units(self.temperature[0], self.temperature[1], 'K')
            p = conv.convert_units(self.pressure[0], self.pressure[1], 'Pa')
            avgMWinv = 0.0
            for species in self.composition.keys():
                if Stream.species[species]['phase'] == 'g':
                    if species in SpecialMolecule.MW.keys():
                        MW = SpecialMolecule.MW[species]
                    else:
                        MW = 0.0
                        breakdown = ep.parse_species(species)
                        try:
                            for ele, v in breakdown.items():
                                MW +=  v*Element.MW[ele]
                        except KeyError:
                            raise BadStreamError, "%s does not have an entry in the Element molecular weight dictionary" % ele
                    avgMWinv += self.composition[species]/MW

            avgMW = 1.0/avgMWinv 
            return val/avgMW*8.314*T/p


    def elementalFactor(self, element):
        """Returns the units of element/basis unit of a given element in the feed"""
        factor = 0.0
        for (specie, fraction) in self.composition.items():
            try:
                spe_dict = Stream.species[specie]
                
            except KeyError:
                try:
                    spe_dict = self.special_species[specie]
                    
                except KeyError:
                    raise SpeciesNotDefinedError, "%s does not have an elemental breakdown definition...check special_species{}" % specie
            try:
                
                factor +=  fraction*spe_dict[element]
            except KeyError:
                pass #It's okay to pass here, because not all species will have the specified element
        
        return factor

    def speciesFactor(self, specie):
        """Returns the units of species/basis unit of a given element in the feed"""
        try:
            return self.composition[specie]
        except KeyError:
            pass

    def calcElementalMolarFlowrate(self, element):
        """Calculates the molar feed rates (mol/s) of an element for a given stream flowrate"""
        if self.basis ==  "molar":
            #First need to convert to mol/s
            conv = uc.UnitConverter()
            val = conv.convert_units(self.flowrate[0], self.flowrate[1], 'mol/s')
            return val*self.elementalFactor(element)
        elif self.basis ==  "mass":
            #First need to convert to g/s
            conv = uc.UnitConverter()
            val = conv.convert_units(self.flowrate[0], self.flowrate[1], 'g/s')
            
            return val*self.elementalFactor(element)/Element.MW[element]
        elif self.basis ==  "gas_volume":
            #First need to convert to m^3/s
            conv = uc.UnitConverter()
            val = conv.convert_units(self.flowrate[0], self.flowrate[1], 'm^3/s')
            p = conv.convert_units(self.pressure[0], self.pressure[1], 'kg/s^2/m')
            T = conv.convert_units(self.temperature[0], self.temperature[1], 'K')
            #Need to convert the temperature to K
            return val*p/(8.314*T)*self.elementalFactor(element)
        elif self.basis == "std_gas_volume":
            #First need to convert to Sm^3/s
            conv = uc.UnitConverter()
            val = conv.convert_units(self.flowrate[0], self.flowrate[1], 'm^3/s')
            p = conv.convert_units(self.std_pressure[0], self.std_pressure[1], 'kg/s^2/m')
            T = conv.convert_units(self.std_temperature[0], self.std_temperature[1], 'K')
            return val*p/(8.314*T)*self.elementalFactor(element)

        elif self.basis ==  "liquid_volume":
            #First need to convert the density to the same units as the liquid
            pass
        else:
            pass

    def calcSpeciesMolarFlowrate(self, species):
        """Calculates the molar feed rate (mol/s) of a species for a given stream flowrate"""
        if species not in self.composition.keys():
            return np.zeros(len(self.flowrate[0]))

        if self.basis ==  "molar":
            conv = uc.UnitConverter()
            val = conv.convert_units(self.flowrate[0], self.flowrate[1], 'mol/s')
            return val*self.composition[species]
        elif self.basis ==  "mass":
            #Need to convert to molar amounts based on the molecular weight
            conv = uc.UnitConverter()
            val = conv.convert_units(self.flowrate[0], self.flowrate[1], 'g/s')
            
            if species in SpecialMolecule.MW.keys():
                MW = SpecialMolecule.MW[species]
            else:
                MW = 0
                breakdown = ep.parse_species(species)
                try:
                    for ele, v in breakdown.items():
                        MW +=  v*Element.MW[ele]
                except KeyError:
                    raise BadStreamError, "%s does not have an entry in the Element molecular weight dictionary" % ele
            return val*self.composition[species]/MW
        elif self.basis ==  "gas_volume":
            
            conv = uc.UnitConverter()
            val = conv.convert_units(self.flowrate[0], self.flowrate[1], 'm^3/s')
            p = conv.convert_units(self.pressure[0], self.pressure[1], 'kg/s^2/m')
            T = conv.convert_units(self.temperature[0], self.temperature[1], 'K')
            return val*self.composition[species]*p/(8.314*T)
        elif self.basis == "std_gas_volume":
            conv = uc.UnitConverter()
            val = conv.convert_units(self.flowrate[0], self.flowrate[1], 'm^3/s')
            
            p = conv.convert_units(self.std_pressure[0], self.std_pressure[1], 'Pa')
            T = conv.convert_units(self.std_temperature[0], self.std_temperature[1], 'K')
            return val*self.composition[species]*p/(8.314*T)

        else:
            raise lflException, "Basis not recognized!"

    def calcSpeciesVolumetricFlowrate(self, species):    
        """Calculates the volumetric feed rate (m^3/s) of a species for a given stream flowrate"""
        #This needs to be generalized to admit liquid states -- also, unit tests must be written!
        if species not in self.composition.keys():
            return np.zeros(len(self.flowrate[0]))

        if self.basis ==  "molar":
            conv = uc.UnitConverter()
            val = conv.convert_units(self.flowrate[0], self.flowrate[1], 'mol/s')
            T = conv.convert_units(self.temperature[0], self.temperature[1], 'K')
            P = conv.convert_units(self.pressure[0], self.pressure[1], 'kg/s^2/m')
            return val*self.composition[species]*8.314*T/P
        elif self.basis ==  "mass":
            #Need to convert to molar amounts based on the molecular weight
            conv = uc.UnitConverter()
            val = conv.convert_units(self.flowrate[0], self.flowrate[1], 'g/s')
            
            if species in SpecialMolecule.MW.keys():
                MW = SpecialMolecule.MW[species]
            else:
                MW = 0
                breakdown = ep.parse_species(species)
                try:
                    for ele, v in breakdown.items():
                        MW +=  v*Element.MW[ele]
                except KeyError:
                    raise BadStreamError, "%s does not have an entry in the Element molecular weight dictionary" % ele
            molar = val*self.composition[species]/MW
            T = conv.convert_units(self.temperature[0], self.temperature[1], 'K')
            P = conv.convert_units(self.pressure[0], self.pressure[1], 'kg/s^2/m')
            return molar * T*8.314/P
        elif self.basis ==  "gas_volume":
            conv = uc.UnitConverter()
            val = conv.convert_units(self.flowrate[0], self.flowrate[1], 'm^3/s')          
            return val*self.composition[species]
        elif self.basis == "std_gas_volume":
            conv = uc.UnitConverter()
            val = conv.convert_units(self.flowrate[0], self.flowrate[1], 'm^3/s')
            T = conv.convert_units(self.temperature[0], self.temperature[1], 'K')
            p = conv.convert_units(self.pressure[0], self.pressure[1], 'Pa')
            T_std = conv.convert_units(self.std_temperature[0], self.std_temperature[1], 'K')
            p_std = conv.convert_units(self.std_pressure[0], self.std_pressure[1], 'Pa')
            val *= T/T_std * p_std/p
            return val*self.composition[species]

        else:
            raise lflException, "Basis not recognized!"

    def _parse_species(self, species_str):
        """Parses a molecule into into a dictionary of {element: # of atoms}"""

        parsed_elements = {}
        current_element = ""
        current_number_str = ""

        i = 0
        if not species_str[0].isalpha():
            raise BadCharacterError, "A molecule must start with an alphabetical character"

        while i < len(species_str):
            current_char = species_str[i]
        
            if current_char.isalpha():
                if i+1 ==  len(species_str) or species_str[i+1].isupper():
                    #Allows for single character names, like CH4
                    current_element = "".join([current_element, current_char])
                    if current_element in parsed_elements.keys():
                        parsed_elements[current_element] +=  1
                    else:
                        parsed_elements[current_element] = 1
                    current_element = ""

                else:
                
                    current_element = "".join([current_element, current_char])
                i +=  1
                continue

            elif current_char.isdigit():
                #we have gotten to the end of an element name
                while i < len(species_str) and species_str[i].isdigit():
                    current_char = species_str[i]
                    current_number_str = "".join([current_number_str, current_char])
                    i +=  1
                if current_number_str ==  '':
                    raise BadCharacterError, "Each element must have a number associated with it"
                if current_element in parsed_elements.keys():
                    parsed_elements[current_element] +=  int(current_number_str)
                else:
                    parsed_elements[current_element] = int(current_number_str)
            
                current_element = ""
                current_number_str = ""
            else:
                raise BadCharacterError, "A molecule can only contain alphabetical and numerical characters"

        return parsed_elements

    def get_enthalpy(self, units):
        if not isinstance(units, str):
            raise UnitConversionError, "The provided units must be a string for the enthalpy function"
        
        self._calc_enthalpy()
        conv = uc.UnitConverter()
        return conv.convert_units(self.enthalpy[0], self.enthalpy[1], units) + conv.convert_units(self.enthalpy_reserve[0], self.enthalpy_reserve[1], units)

    def get_entropy(self, units):
        if not isinstance(units, str):
            raise UnitConversionError, "The provided units must be a string for the enthalpy function"
        
        self._calc_entropy()
        conv = uc.UnitConverter()
        return conv.convert_units(self.entropy[0], self.entropy[1], units)

    def change_biomass_composition(self, bio_comp):
        #need to do some error checking here
        if not isinstance(bio_comp, dict):
            raise BadCTTransInputError, "The new biomass composition must be a dictionary"

        for key in bio_comp:
            if key not in ["CELL","LIGH","LIGC","LIGO","HCE","ASH"]:
                raise BadCTTransInputError, "%s is not a valid biomass component" % key

        for value in bio_comp.values():
            try:
                int(value)
            except ValueError:
                raise BadCTTransInputError, "All compositional values must be numeric"

            if not value >= 0.0:
                raise BadCTTransInputError, "All compositional values must be greater than zero"

        Stream.ct_trans['biomass'] = bio_comp
            


    def ct_setcomp(self, composition = None):
        """Returns a string that can be used to input the stream's composition into a Cantera object"""
        # Missing transport data for H2S necessary to work with Cantera.  Prevent Stream object from putting H2S into Cantera composition, for now.  Shouldn't have too much of an impact in calculations, as H2S is a trace species.
        ignore_species = ['H2S']
        
        if composition == None:
            composition = self.composition
            
        for specie in ignore_species:
            if specie in composition:
                del(composition[specie])
                
        #First, need to format everything into the Cantera Species
        if self.basis == 'molar' or self.basis == 'gas_volume' or self.basis == 'std_gas_volume':
            #This is the easiest case -- just pull all the compositional values and set the phase appropriately
            set_string = ""
            for specie in composition.keys():
                if specie in Stream.ct_trans:
                    for sub in Stream.ct_trans[specie]:
                        set_string += '%s:%s, ' % (sub, Stream.ct_trans[specie][sub]*composition[specie])  #This, of course, assumes that the basis in ct_trans is the same as in composition
                else:
                    set_string += '%s:%s, ' % (specie, composition[specie])
            #set the phase composition
            set_string = set_string[:-2] #Remove last ', ' from set_string build
            self.cantera_helper.setMoleFractions(self.ctphase, set_string)
            #self.ctphase.setMoleFractions(set_string)
            

        elif self.basis == 'mass':
            #This is also pretty easy -- just pull all the compositional values and set the phase appropriately
            set_string = ""
            for specie in composition.keys():
                if specie in Stream.ct_trans:
                    for sub in Stream.ct_trans[specie]:
                        set_string += '%s:%s, ' % (sub, Stream.ct_trans[specie][sub]*composition[specie])  #This, of course, assumes that the basis in ct_trans is the same as in composition
                else:
                    set_string += '%s:%s, ' % (specie, composition[specie])
            #set the phase composition
            set_string = set_string[:-2] #Remove last ', ' from set_string build
            self.cantera_helper.setMassFractions(self.ctphase, set_string)
            #self.ctphase.setMassFractions(set_string)

        else:
            raise lflException, '%s is not a valid stream basis' % self.basis

    def _calc_enthalpy(self):
        """Calculates the stream enthalpy and stores it in self.enthalpy"""
        
        if self.temperature==None:
            raise BadStreamError, '%s Stream temperature needs to be defined to calculate enthalpy.' % self.name
        if self.pressure==None:
            raise BadStreamError, '%s Stream pressure needs to be defined to calculate enthalpy.' % self.name
        if self.composition ==None:
            raise BadStreamError, '%s Stream composition needs to be defined to calculate enthalpy.' % self.name
        conv = uc.UnitConverter()
        
        #get the specific enthalpy -- for the vector case will need to loop through and build the enthalpy vectors
        getattr(self, "_calc_spec_enthalpy_%s" % self.mode)()

        #self.ct_setcomp(self.composition)
        #self.ctphase.set(T = conv.convert_units(self.temperature[0], self.temperature[1], 'K'), P = conv.convert_units(self.pressure[0], self.pressure[1], 'Pa'))
         
        #Cantera output is J/kmol or J/kg, so conversions must follow this for molar and mass flow rates.
        if self.basis == 'molar':
            #convert to kmol/s:
            
            flow = conv.convert_units(self.flowrate[0], self.flowrate[1], 'kmol/s')

        elif self.basis == 'mass':
            #convert to kg/s
            flow = conv.convert_units(self.flowrate[0], self.flowrate[1], 'kg/s')
                       
            
        elif self.basis ==  "gas_volume":
            val = conv.convert_units(self.flowrate[0], self.flowrate[1], 'm^3/s')
            p = conv.convert_units(self.pressure[0], self.pressure[1], 'Pa')
            T = conv.convert_units(self.temperature[0], self.temperature[1], 'K')
            flow =  val*p/(8.314*T)/1000
            
        elif self.basis == "std_gas_volume":
            val = conv.convert_units(self.flowrate[0], self.flowrate[1], 'm^3/s')
            p = conv.convert_units(self.std_pressure[0], self.std_pressure[1], 'Pa')
            T = conv.convert_units(self.std_temperature[0], self.std_temperature[1], 'K')
            flow =  val*p/(8.314*T)/1000

        enthalpy = flow*self.spec_enthalpy
        self.enthalpy = (enthalpy, 'J/s')

    def _calc_spec_enthalpy_scalar(self):
        """Internal function for calculating the specific enthalpy for scalar valued streams"""
        self.ct_setcomp(self.composition)
        self.cantera_helper.set(self.ctphase, T = conv.convert_units(self.temperature[0], self.temperature[1], 'K'), P = conv.convert_units(self.pressure[0], self.pressure[1], 'Pa')) 
        #self.ctphase.set(T = conv.convert_units(self.temperature[0], self.temperature[1], 'K'), P = conv.convert_units(self.pressure[0], self.pressure[1], 'Pa'))
        if self.basis == "mass":
            self.spec_enthalpy = self.cantera_helper.enthalpy_mass(self.ctphase)
            #self.spec_enthalpy = self.ctphase.enthalpy_mass()	#Units will be J/kg
        else:
            self.spec_enthalpy = self.cantera_helper.enthalpy_mole(self.ctphase)		#Everything else is a molar case
	    #self.spec_enthalpy = self.ctphase.enthalpy_mole()

    def _calc_spec_enthalpy_vector(self):
        """Internal function for calculating the specific enthalpy for vector valued streams"""
        #Vector streams need to loop through the temperature, pressure, and composition lines to set the Cantera phase one by one
        self.spec_enthalpy = np.zeros(len(self.temperature[0]))
        self.spec_enthalpy[:] = np.nan				#appropriate until values are filled in
        if self.basis == "mass":
            kw = "mass"
        else:
            kw = "mole"
        for i in range(0, len(self.temperature[0])):
            #need to build a composition dictionary to set the Cantera phase
            comp = {}
            for k in self.composition:
                comp[k] = self.composition[k][i]
            self.ct_setcomp(comp)
            self.cantera_helper.set(self.ctphase, T = conv.convert_units(self.temperature[0][i], self.temperature[1], 'K'), P = conv.convert_units(self.pressure[0][i], self.pressure[1], 'Pa'))
            #self.ctphase.set(T = conv.convert_units(self.temperature[0][i], self.temperature[1], 'K'), P = conv.convert_units(self.pressure[0][i], self.pressure[1], 'Pa'))
            self.spec_enthalpy[i] = getattr(self.cantera_helper, "enthalpy_%s" % kw)(self.ctphase)
            #self.spec_enthalpy[i] = getattr(self.ctphase, "enthalpy_%s" % kw)()


    def _calc_entropy(self):
        """Calculates the stream entropy and stores it in self.entropy"""
        conv = uc.UnitConverter()
        if self.temperature==None:
            raise BadStreamError, 'Stream temperature is not defined.'
        if self.pressure==None:
            raise BadStreamError, 'Stream pressure is not defined.'

        #set the Cantera phase
        self.ct_setcomp(self.composition)
        self.cantera_helper.set(self.ctphase, T = conv.convert_units(self.temperature[0], self.temperature[1], 'K'), P = conv.convert_units(self.pressure[0], self.pressure[1], 'Pa'))
        #self.ctphase.set(T = conv.convert_units(self.temperature[0], self.temperature[1], 'K'), P = conv.convert_units(self.pressure[0], self.pressure[1], 'Pa'))

        if self.basis == 'molar':
            #convert to kmol/s:
            flow = conv.convert_units(self.flowrate[0], self.flowrate[1], 'kmol/s')
            self.entropy = [flow*self.cantera_helper.entropy_mole(self.ctphase), 'J/K/s']
            #self.entropy = [flow*self.ctphase.entropy_mole(), 'J/K/s']

        elif self.basis == 'mass':
            #convert to kg/s
            flow = conv.convert_units(self.flowrate[0], self.flowrate[1], 'kg/s')
            self.entropy = [flow*self.cantera_helper.entropy_mass(self.ctphase), 'J/K/s']
            #self.entropy = [flow*self.ctphase.entropy_mass(), 'J/K/s']

        elif self.basis ==  "gas_volume":
                        
            val = conv.convert_units(self.flowrate[0], self.flowrate[1], 'm^3/s')
            p = conv.convert_units(self.pressure[0], self.pressure[1], 'kg/s^2/m')
            T = conv.convert_units(self.temperature[0], self.temperature[1], 'K')
            flow =  val*p/(8.314*T)
            self.entropy = [flow*self.cantera_helper.entropy_mole(self.ctphase), 'J/K/s']
            #self.entropy = [flow*self.ctphase.entropy_mole(), 'J/K/s']

        elif self.basis == "std_gas_volume":
            
            val = conv.convert_units(self.flowrate[0], self.flowrate[1], 'm^3/s')
            
            p = conv.convert_units(self.std_pressure[0], self.std_pressure[1], 'Pa')
            T = conv.convert_units(self.std_temperature[0], self.std_temperature[1], 'K')
            flow =  val*p/(8.314*T)
            self.entropy = [flow*self.cantera_helper.entropy_mole(self.ctphase), 'J/K/s']
            #self.entropy = [flow*self.ctphase.entropy_mole(), 'J/K/s']
        
    def convert_to_mass_basis(self):
        if self.basis == 'mass':
            #do nothing
            pass
        
        elif self.basis == 'gas_volume':
            #I'm going to be sneaky here and convert to molar first, then let the execution drop downward
            conv = uc.UnitConverter()
            p = conv.convert_units(self.pressure[0], self.pressure[1], 'Pa')
            T = conv.convert_units(self.temperature[0], self.temperature[1], 'K')
            self.flowrate = [conv.convert_units(self.flowrate[0], self.flowrate[1], 'm^3/s')*p/(8.314*T), 'mol/s']
            self.basis = 'molar'

        elif self.basis == 'std_gas_volume':
            conv = uc.UnitConverter()
            p = conv.convert_units(self.std_pressure[0], self.std_pressure[1], 'Pa')
            T = conv.convert_units(self.std_temperature[0], self.std_temperature[1], 'K')
            self.flowrate = [conv.convert_units(self.flowrate[0], self.flowrate[1], 'm^3/s')*p/(8.314*T), 'mol/s']
            self.basis = 'molar'

        if self.basis == 'molar':
            #calculate the molar compositions from the mass compositions
            
            conv = uc.UnitConverter()
            #calculate the average molecular weight first
            species_dict = {}
            for species in self.composition:
            
                if species in SpecialMolecule.MW.keys():
                    MW = SpecialMolecule.MW[species]
                else:
                    MW = 0
                    breakdown = ep.parse_species(species)
                    try:
                        for ele, v in breakdown.items():
                            MW +=  v*Element.MW[ele]
                    except KeyError:
                        raise BadStreamError, "%s does not have an entry in the Element molecular weight dictionary" % ele
                species_dict[species] = MW*self.composition[species]
            avg_MW = sum(species_dict.values())
            for species in species_dict:
                self.composition[species] = species_dict[species]/avg_MW
            self.flowrate = [avg_MW*conv.convert_units(self.flowrate[0],self.flowrate[1], 'mol/s'), 'g/s']
            self.basis = 'mass'

class CanteraHelper:
    """This class is a wrapper to allow the cantera functions to work in the same way as before"""
    def __init__(self):
        if vch.LooseVersion(ct.__version__) < vch.LooseVersion("2.1.1"):
            self.api = "old"
        else:
            self.api = "new"

    def importPhase(self, filename, phasename=None):
        return getattr(self, "_importPhase_%s" % self.api)(filename, phasename)

    def setMoleFractions(self, phase, mole_fraction_string):
        return getattr(self, "_sety_%s" % self.api)(phase, mole_fraction_string)

    def setMassFractions(self, phase, mass_fraction_string):
        return getattr(self, "_setw_%s" % self.api)(phase, mass_fraction_string) 

    def set(self, phase, T, P, X=None, Y=None):
        return getattr(self, "_set_%s" % self.api)(phase,T,P,X,Y)

    def enthalpy_mole(self, phase):
        return getattr(self, "_enthalpy_mole_%s" % self.api)(phase)

    def enthalpy_mass(self, phase):
        return getattr(self, "_enthalpy_mass_%s" % self.api)(phase)

    def entropy_mole(self, phase):
        return getattr(self, "_entropy_mole_%s" % self.api)(phase)

    def entropy_mass(self, phase):
        return getattr(self, "_entropy_mass_%s" % self.api)(phase)
    


    def _importPhase_old(self, filename, phasename):
        if phasename is not None:
            return ct.importPhase(filename, phasename) 
        else:
            return ct.importPhase(filename)

    def _importPhase_new(self, filename, phasename):
        if phasename is not None:
            return ct.Solution(filename, phasename)
        else:
           return ct.Solution(filename)

    def _sety_old(self, phase, y_string):
        phase.setMoleFractions(y_string)

    def _sety_new(self, phase, y_string):
        phase.X = y_string

    def _setw_old(self, phase, w_string):
        phase.setMassFractions(w_string)
    
    def _setw_new(self, phase, w_string):
        #print w_string
        phase.Y = w_string

    def _set_old(self, phase, T, P, X, Y):
        if X is not None:
            phase.set(T=T, P=P, X=X)
        elif Y is not None:
            phase.set(T=T, P=P, Y=Y)
        else:
            phase.set(T=T, P=P)

    def _set_new(self, phase, T,P,X,Y):
        if X is not None:
            phase.TPX = T,P,X
        elif Y is not None:
            phase.TPY = T,P,Y
        else:
            phase.TP = T,P

    def _enthalpy_mole_old(self, phase):
        return phase.enthalpy_mole()

    def _enthalpy_mole_new(self, phase):
        return phase.enthalpy_mole

    def _enthalpy_mass_old(self, phase):
        return phase.enthalpy_mass()

    def _enthalpy_mass_new(self, phase):
        return phase.enthalpy_mass

    def _entropy_mole_old(self, phase):
        return phase.entropy_mole()

    def _entropy_mole_new(self, phase):
        return phase.entropy_mole

    def _entropy_mass_old(self, phase):
        return phase.entropy_mass()

    def _entropy_mass_new(self, phase):
        return phase.entropy_mass


class ProcessObject:
    
    def __init__(self, inlets, outlets = None):
        
        #Inlet Checks.  Maybe in future don't need to have inlets fully defined if outlets are, but not sure how to implement now.
        if type(inlets) != list:
            raise BadStreamError, 'Inlets must be input as a list of Stream objects.'
        for inlet in inlets:
            if not isinstance(inlet, Stream):
                raise BadStreamError, 'Inlet %s is not a Stream object.' %inlet.name
        #!!!FIX!!! Should check that all inlets are of the same length            
 
        #Outlet Checks
        if outlets is not None:
            if type(outlets) != list:
                raise BadStreamError, 'Outlets must be input as a list of Stream objects.'
            for outlet in outlets:
                if not isinstance(outlet, Stream):
                    raise BadStreamError, 'Outlet %s is not a Stream object.' %outlet.name

    #There has GOT to be a way to generalize this with decorators or similar
    
        self.inlets = inlets
        self.outlets = outlets
        
    def totalInletEnthalpy(self, units):
        H = 0
        for stream in self.inlets:
            
            H += stream.get_enthalpy(units)
        return H

    def totalOutletEnthalpy(self, units):
        H = 0
        for stream in self.outlets:
            
            H += stream.get_enthalpy(units)
        return H

    def totalInletEntropy(self,units):
        S = 0
        for stream in self.inlets:
            S += stream.get_entropy(units)
        return S

    def totalOutletEntropy(self,units):
        S = 0
        for stream in self.outlets:
            S += stream.get_entropy(units)
        return S

    def deltaH(self, units):
        #print "Total Inlet Enthalpy:\t%s" % self.totalInletEnthalpy(units)[0]
        #print "Total Outlet Enthalpy:\t%s" % self.totalOutletEnthalpy(units)[0]
        return self.totalOutletEnthalpy(units) - self.totalInletEnthalpy(units)
    
    def deltaS(self, units):
        return self.totalOutletEntropy(units) - self.totalInletEnthalpy(units)
            
class Mixer(ProcessObject):
    
    """A mixer blends multiple inlet streams into one outlet stream"""

    def __init__(self, name, outlet_pressure = None, temp_method = 'default', **kwargs):
        ProcessObject.__init__(self, **kwargs)
        self.temp_method = temp_method
        
        self.outlets = [Stream(name = '%s_outlet' % name)]
        
        #Need to solve for the outlet stream pressure, temperature, and composition
        #pressure is easy -- it is assumed that the pressure drops to the LOWEST stream pressure entering the mixer if an outlet pressure is not specified
        self._calc_outlet_pressure(outlet_pressure)
	
        self._calc_outlet_flowrate()
        self._calc_outlet_temperature()

        #need to drag along the enthalpy reserve from the inlet streams -- will be fixed when streams get n-phase capability
        conv = uc.UnitConverter()
        for inlet in self.inlets:
            self.outlets[0].enthalpy_reserve[0] += conv.convert_units(inlet.enthalpy_reserve[0], inlet.enthalpy_reserve[1], 'W')


    def recalc(self, outlet_pressure = None):
        self._calc_outlet_pressure(outlet_pressure)
        self._calc_outlet_flowrate()
        self._calc_outlet_temperature()

    def _calc_outlet_pressure(self, outlet_pressure):
        conv = uc.UnitConverter()
        if outlet_pressure is None:

            minP = None
            mode = "scalar"
            for stream in self.inlets:
                if stream.mode == "vector":
                    mode = "vector"
        
            getattr(self, "_calc_outlet_pressure_%s" % mode)()

        else:
            self.outlets[0].set_pressure(outlet_pressure)

    def _calc_outlet_pressure_scalar(self):
        conv = uc.UnitConverter()
        minP = None
        
        for stream in self.inlets:
            if minP is None or conv.convert_units(stream.pressure[0], stream.pressure[1], 'Pa') < minP:
                minP = conv.convert_units(stream.pressure[0], stream.pressure[1], 'Pa')
        
        self.outlets[0].set_pressure((minP, 'Pa'))

    def _calc_outlet_pressure_vector(self):
        
        conv = uc.UnitConverter()
        minP = np.zeros(self.inlets[0].length)
        minP[:] = np.nan
        
        for i in range(0,self.inlets[0].length):
            for stream in self.inlets:
                
                if np.isnan(minP[i]) or conv.convert_units(stream.pressure[0][i], stream.pressure[1], 'Pa') < minP[i]:
                    minP[i] = conv.convert_units(stream.pressure[0][i], stream.pressure[1], 'Pa')
        self.outlets[0].set_pressure((minP, 'Pa'))
        

    def _calc_outlet_flowrate(self):
        conv = uc.UnitConverter()
        #Need to put everything on a consistent basis - if everything the same, just add them all together; otherwise use mass
        #Also, need to set the composition of the outlet
        c = True
        basis_fl_dict = {'molar':'mol/s', 'mass':'kg/s', 'gas_volume':'m^3/s', 'std_gas_volume':'m^3/s'}
        basis_choice = self.inlets[0].basis
        
        #if we have vectors for outlet flowrates, convert all the streams to vector streams
        mode = None
        for inlet in self.inlets:
            
            if isinstance(inlet.flowrate[0], pd.Series) or isinstance(inlet.flowrate[0], np.ndarray):
                mode = "vector"
                length = len(inlet.flowrate[0])
        if mode == "vector":
            for inlet in self.inlets:
                inlet.convert_scalar_to_vector(length)

        for inlet in self.inlets:
            if inlet.basis != basis_choice or inlet.basis == 'gas_volume' or inlet.basis == 'std_gas_volume':
                c = False

        if not c:
            #convert all streams to a mass basis
            for inlet in self.inlets:
                if inlet.basis != 'mass':
                    inlet.convert_to_mass_basis()
            basis_choice = 'mass'
        
        self.outlets[0].basis = basis_choice
        fl_sum = 0
        for inlet in self.inlets:
            fl_sum += conv.convert_units(inlet.flowrate[0], inlet.flowrate[1], basis_fl_dict[basis_choice])
        self.outlets[0].flowrate = (fl_sum, basis_fl_dict[basis_choice])
        #need to generate a total species list for compositional matching - then composition is simply the sum of the streams for that species divided by the total flowrate
        species_list = []
        for inlet in self.inlets:
            for species in inlet.composition: 
                if species not in species_list:
                    species_list.append(species)
        composition = {}
        for species in species_list:
            spec_sum = 0            
            for inlet in self.inlets:
                if species in inlet.composition: 
                    spec_sum += conv.convert_units(inlet.flowrate[0], inlet.flowrate[1], basis_fl_dict[basis_choice])*inlet.composition[species]
            composition[species] = spec_sum/fl_sum
        
        self.outlets[0].set_composition(composition)

        #Just driving all volume directly to mass now -- should work
        #!!!#Need a function that converts a scalar stream to a vector stream

    def enth_func(self, T):
        self.outlets[0].set_temperature((T, 'K'))      
        d_H = self.deltaH('J/s')
        return d_H
    
    def _calc_outlet_temperature(self):
        #Solving the equation dH = 0 (not a heat exchanger, so Q and W are both 0)
        #guess a temperature for the outlet as a mean of the inlet temperatures
        
        conv = uc.UnitConverter()
        temp_sum = 0.0
        for inlet in self.inlets:
            temp_sum += conv.convert_units(inlet.temperature[0], inlet.temperature[1], 'K')
        
        temp_avg = temp_sum/len(self.inlets)
                
        if self.outlets[0].mode == "vector": 
            if self.temp_method == 'default':

                outlet_temp = spo.newton_krylov(F = self.enth_func, xin = temp_avg, f_tol = 1E-3)
            elif self.temp_method == 'fast_mean':
                #Create streams for each of the inlet streams at the mean value of flowrate, temperature, and pressure
                temp_streams = []
                for inlet in self.inlets:
                    T = [st.nanmean(inlet.temperature[0]), inlet.temperature[1]]
                    P = [st.nanmean(inlet.pressure[0]), inlet.pressure[1]]
                    F = [st.nanmean(inlet.flowrate[0]), inlet.flowrate[1]]
                    #need an average inlet composition
                    comp = {}
                    for specie in inlet.composition:
                        if F[0] != 0:                        
                            comp[specie] = st.nanmean(inlet.composition[specie]*inlet.flowrate[0])/F[0]
                        else:
                            comp[specie] = st.nanmean(inlet.composition[specie])

                    temp_streams.append(Stream(name = inlet.name, flowrate = F, temperature = T, pressure = P, basis = inlet.basis, composition = comp))  #switch composition to inlet.composition to get old behavior
                

                temp_m = Mixer(name = 'mean_mix', inlets = temp_streams)
                
                outlet_temp = temp_m.outlets[0].temperature[0]*np.ones(len(self.inlets[0].flowrate[0]))
                

        elif self.outlets[0].mode == "scalar":
            
            minT = np.inf
            maxT = -np.inf
            for inlet in self.inlets:
                T1 = conv.convert_units(inlet.temperature[0], inlet.temperature[1], 'K')
                if T1 < minT:
                    minT = T1
                if T1 > maxT:
                    maxT = T1
               
                   
           
            
            outlet_temp = spo.bisect(f = self.enth_func,a = 1.0, b=1500.0)
            
        self.outlets[0].set_temperature((outlet_temp, 'K'))
        

        
class Reactor(ProcessObject):
    """Reactor Class..."""
    def __init__(self, temperature = None, pressure = None, **kwargs):

        ProcessObject.__init__(self,**kwargs)
        self.temperature = temperature
        self.pressure = pressure

    def calc_species_generation(self):
        pass
    def calc_species_consumption(self):
        pass
    def calc_enthalpy_change(self, units):
        return self.deltaH(units)
    def calc_entropy_change(self, units):
        return self.deltaS(units)

class SDFIdealGasifier(Reactor):
    """Idealized SDF gasifier -- 100% conversion"""
    def __init__(self, y_CO2_out, y_CH4_out, **kwargs):
        Reactor.__init__(self, **kwargs)
        self.y_CO2 = y_CO2_out
        self.y_CH4 = y_CH4_out

        #need to convert the biomass in the inlet streams to the components
        for inlet in self.inlets:
            if 'biomass' in inlet.composition:
                for key in Stream.ct_trans['biomass']:
                    inlet.composition[key] = Stream.ct_trans['biomass'][key]*inlet.composition['biomass']
                del inlet.composition['biomass']
                  


    def generate_outlet_stream(self):
        #calculate the molar flowrates of each of the species
	species = []
        for inlet in self.inlets:
            for key in inlet.composition:
                if key not in species:
                    species.append(key)
            

        flowrates = {}
        for key in species:
            for inlet in self.inlets:
                 try:
                     flowrates[key] += inlet.calcSpeciesMolarFlowrate(key)
                 except KeyError:
                     flowrates[ key] = inlet.calcSpeciesMolarFlowrate(key)

        inc_spec = ["CELL", "HCE", "LIGC", "LIGO", "LIGH", "H2", "CO", "CO2", "CH4", "H2O", "N2", "Ar"]
        for spec in inc_spec:
            if spec not in flowrates:
                 flowrates[spec] = np.zeros(len(self.inlets[0].flowrate[0]))

        #for spec in flowrates:
        #    print "%s:\t%s" % (spec, flowrates[spec][0])


        
        #calculate grouped parameters
        xi1 = flowrates["CELL"]
        xi2 = flowrates["HCE"]
        xi3 = flowrates["LIGC"]
        xi4 = flowrates["LIGO"]
        xi5 = flowrates["LIGH"]
        #phi0 = flowrates["H2"] + flowrates["CO"] + flowrates["CO2"] + flowrates["CH4"]
        #phi1 = 12*xi1 + 10*xi2 + 33*xi3 + 41*xi4 + 49*xi5
        xi6 = 0.0 #setting this to zero - essentially turning off water-gas shift(self.y_CO2*(phi0+phi1)+3*self.y_CO2*flowrates["CH4"] + 3*self.y_CH4*flowrates["CO2"] - flowrates["CO2"])/(1+3*self.y_CH4-self.y_CO2)
        xi7 = 0.0 #setting this to zero - essentially turning off steam-reforming(self.y_CH4/self.y_CO2)*(flowrates["CO2"]+xi6) - flowrates["CH4"]     
   
        #now do the mass balances
        outs = {}
        outs['CELL'] = flowrates['CELL'] - xi1
        outs['HCE'] = flowrates['HCE'] - xi2
        outs['LIGC'] = flowrates['LIGC'] - xi3
        outs['LIGO'] = flowrates['LIGO'] - xi4
	outs['LIGH'] = flowrates['LIGH'] - xi5
        outs['H2'] = flowrates['H2'] +6*xi1 + 5*xi2 + 18*xi3 + 21*xi4 + 27*xi5 + xi6 - 3*xi7
        outs['CO'] = flowrates['CO'] +6*xi1 + 5*xi2 + 15*xi3 + 20*xi4 + 22*xi5 - xi6 - xi7
        outs['CO2'] = flowrates['CO2'] + xi6
        outs['CH4'] = flowrates['CH4'] + xi7
        outs['H2O'] = flowrates['H2O'] -xi1 -xi2 - 11*xi3 - 10*xi4 - 13*xi5 - xi6 + xi7
        outs['N2'] = flowrates['N2']
        outs['Ar'] = flowrates['Ar']

        #for out in outs:
        #     print "%s:\t%s" % (out, outs[out][0])


        #set up the outlet stream
        fr = 0.0
        for specie in outs:
            fr += outs[specie]

        composition = {}
        for specie in outs:
            composition[specie] = outs[specie]/fr

        outlet = Stream('max_dH_outlet', flowrate = [fr, 'mol/s'], temperature = self.temperature, pressure = self.pressure, composition = composition, basis = "molar")

        self.outlets = [outlet,]
     

        

        
class Condenser(ProcessObject):
    def __init__(self, **kwargs):
        ProcessObject.__init__(**kwargs)
        
class PhysicalConstants:
    """Class to hold global physical constant variables"""
    std_temp = 298   #In K
    std_pressure = 100000 #In Pa

class Element:
    """A physical element.  Will have properties of the elements, plus a handy reference library"""
    Elements = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au','Hg', 'Tl', 'Pb','Bi', 'Po', 'At', 'Rn', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu']
    MW = {}
    MW['C'] = 12.0107
    MW['H'] = 1.00794
    MW['O'] = 15.9994
    MW['N'] = 14.0067
    MW['S'] = 32.065
    MW['Cl'] = 35.453
    MW['Ar'] = 39.948

class Molecule:
    def __init__(self, formula):
        self.formula = formula
        self.form_dict = ep.parse_species(self.formula)

    def MW(self):
        MW = 0.0
        
        for element in self.form_dict:
            MW += Element.MW[element] * self.form_dict[element]
        return MW

class SpecialMolecule(Molecule):
    """This class holds specially designated molcules, so that I can calculate molecular weights without too much of a problem"""
    MW = {}
    MW['LIGH'] = 436
    MW['LIGO'] = 422
    MW['LIGC'] = 258
    MW['CELL'] = 162
    MW['HCE'] = 132

class ProcTS(ts_data):

    """Timeseries data for chemical processes"""
    def __init__(self, start, end, data = None, units_dict = None):
        ts_data.__init__(self, start, end, data, units_dict)
        
        self.inlet_streams = []
        self.outlet_streams = []
        self.proc_elements = []
        self.proc_species = []
        self.inert_species = []

    

    def generate_inlet_outlet_elemental_flows(self, name_qualifier = None):
        """Generates the standard inlet and outlet elemental flows for elements in self.proc_elements"""
        if self.inlet_streams ==  [] or self.outlet_streams ==  []:
            raise BadStreamError, "Inlet/outlet streams not yet defined"
        
        if name_qualifier ==  None:
            outlet_name = 'outlet'
            inlet_name = 'inlet'

        else:
            if type(name_qualifier) !=  str:
                raise lflException, "The new column name qualifier must be a string"
            else:
                outlet_name = 'outlet_%s' % name_qualifier
       
        for element in self.proc_elements:
            if element not in Element.Elements:
                raise BadStreamError, "The given element: %s is not a physical element" % element

            self['%s_%s' % (element,inlet_name)] = self.calc_elemental_molar_feedrates(self.inlet_streams, element)
            self['%s_%s' % (element,outlet_name)] = self.calc_elemental_molar_feedrates(self.outlet_streams, element)
            self.units['%s_%s' % (element,inlet_name)] = 'mol/s'
            self.units['%s_%s' % (element,outlet_name)] = 'mol/s'


    def calc_elemental_molar_feedrates(self, stream_list, element):
        """Returns the elemental molar feedrate in a given set of streams"""
        tot = 0        
        

        for stream in stream_list:
            
            if not isinstance(stream, Stream):
                
                raise BadStreamError, "%s is not a stream, which is required for elemental calculations" % stream
            
            tot +=  stream.calcElementalMolarFlowrate(element)

        return tot

    def outlet_stream_from_tracer(self, inlet_streams, basis, tracer_species, outlet_tracer_concentration, name='outlet'):
        """Generates an outlet stream, on a determined basis, from a mass balance over a list of inlet streams on the tracer species"""
        #inlet_streams is a list of streams, basis is either "Mass" or "Molar", tracer_species is a string, and tracer_concentration is a numpy array
        #generate the total amount of tracer species coming into the system

        total_tracer_in = np.zeros(len(inlet_streams[0].flowrate[0]))
        #a len function in stream would be nice!
        
        for stream in inlet_streams:
                    
            if basis == "Molar" or basis == "Mass":         #Mass not working yet

                total_tracer_in += getattr(stream, "calcSpecies%sFlowrate" % basis)(tracer_species)


        #try to break the error around outlet_tracer_concentration == 0
        outlet_tracer_concentration[outlet_tracer_concentration<=0] = np.nan
  
        total_outlet_flowrate = total_tracer_in/outlet_tracer_concentration
        
        #create the return Stream
        return Stream(name, flowrate = (total_outlet_flowrate,"mol/s"), basis = "molar")

    def generate_inlet_outlet_species_flows(self, name_qualifier = None):
        """Creates inlet and outlet flowrates for the species of interest in the self.proc_species list"""
        if self.inlet_streams ==  [] or self.outlet_streams ==  []:
            raise BadStreamError, "Inlet/outlet streams not yet defined"

        if name_qualifier ==  None:
            outlet_name = 'outlet'
            inlet_name = 'inlet'

        else:
            if type(name_qualifier) !=  str:
                raise lflException, "The new column name qualifier must be a string"
            else:
                outlet_name = 'outlet_%s' % name_qualifier
                inlet_name = 'inlet_%s' % name_qualifier
        
        for specie in self.proc_species:
        
            self['%s_%s' % (specie, inlet_name)] = self.calc_species_molar_feedrates(self.inlet_streams, specie)
            self['%s_%s' % (specie,outlet_name)] = self.calc_species_molar_feedrates(self.outlet_streams, specie)
            self.units['%s_%s' % (specie, inlet_name)] = 'mol/s'
            self.units['%s_%s' % (specie,outlet_name)] = 'mol/s'

    def calc_species_molar_feedrates(self, stream_list, specie):
        """Returns the species molar feedrates in a given set of streams"""
        tot = 0
        for stream in stream_list:
            if not isinstance(stream, Stream):

                raise BadStreamError, "%s is not a stream, which is required for elemental calculations" % stream
            
            tot +=  stream.calcSpeciesMolarFlowrate(specie)
        return tot

    def _calc_conversion(self, inlet_carbon = None, outlet_carbon = None):
        """Calculates a conversion from the given lists of names"""
        if inlet_carbon is None or outlet_carbon is None:
            pass #raise an error here on a problematic list


        #OK, I know, this seems rather unsophisticated and stupid.  However, it is REQUIRED to have simple functions that return an array for uncertainty determination purposes.
        try:
            conv_val = outlet_carbon/inlet_carbon
        except KeyError:
            raise df.NoColumnError, "The desired column is not in the data frame"
        return conv_val

    def _calc_normalized_comp(self,species_to_normalize, excluded_list = None):
        """Calculates a normalized composition given a list of inerts to exclude"""
        try:
            if excluded_list is None:
                return self[species_to_normalize]
            arr_list = []
            for spec in excluded_list:
                arr_list.append(self[spec])
            return self[species_to_normalize]/(1.0-sum(arr_list))
        except KeyError:
            raise df.NoColumnError, "One of the columns in the excluded list or the species to normalize was not defined before attempting normalization"
        except ValueError:
            raise lflException, "You tried to divide by zero!"

    def inlet_enthalpy(self, units):
        """Returns the inlet enthalpy in the process"""
        H = 0.0
        for stream in self.inlet_streams:
            H +=  stream.get_enthalpy(units)
        return (H, units)

    def outlet_enthalpy(self, units):
        """Returns the outlet enthalpy in the process"""
        H = 0.0
        for stream in self.outlet_streams:
            H +=  stream.get_enthalpy(units)
        return (H, units)

    def enthalpy_change(self, units):
        """Returns the enthalpy change in the process"""
        #Because we cannot measure the C coming out, we need to make up for it in an outlet stream
        #In a similar manner, we need to account for the water that is not measured and set the outlet temperature

        #Another way to accomplish this would be to multiply the conversion against the total dH_max -- this would give an approximate answer that is probably just as good
        #And that is how we are going to calculate this approximation

        #Create a reactor object with inlet and outlet streams
        r = Reactor(inlets = self.inlet_streams, outlets = self.outlet_streams)

        return r.calc_enthalpy_change(units)
        #return r.deltaH(units) -- may need to use this if the other does not work
    
    def generate_enthalpy_change(self, units):
        """Generates a new column with enthalpy change in the data series"""
        self['delta_H'] = self.enthalpy_change(units)
              
        self.units['delta_H'] = units
        
    def collected_enthalpy(self, streams, units):
        """Returns the enthalpy for a collection of streams"""
        H = 0.0
        for stream in streams:
            H +=  stream.get_enthalpy(units)
        return (H, units)

    def enthalpy_diff(self, stream_list_1, stream_list_2, units):
        return self.collected_enthalpy(stream_list_2, units) - self.collected_enthalpy(stream_list_1, units)

class GasifierProcTS(ProcTS):
    """Timeseries data for the gasification chemical process"""
    def __init__(self, start, end, data = None, units_dict = None):
        ProcTS.__init__(self, start, end, data, units_dict)

    

    def generate_carbon_conversions(self, inlet_qualifier = "", qualifier = ""):
        """This calculates the four "standard" conversions/yields that we look for in the biomass gasification on an instantaneous basis"""
        try:
            #CO Yield
            self['CO_yield%s' % qualifier] = self._calc_conversion(self['C_inlet%s' % inlet_qualifier]-self['CO2_inlet%s'% inlet_qualifier],self['CO_outlet%s' % qualifier])
            self.units['CO_yield%s' % qualifier] = None
            #filter these values to throw out the high ones
            self.filter_vals('CO_yield%s' % qualifier, 1.00, 'high')
            self.filter_vals('CO_yield%s' % qualifier, 0.00, 'low')
       
            #Good conversion
            self['X_good%s' % qualifier] = self._calc_conversion(self['C_inlet%s' % inlet_qualifier]-self['CO2_inlet%s' % inlet_qualifier],sum([self['CO_outlet%s' % qualifier],self['CO2_outlet%s' % qualifier]])-self['CO2_inlet%s' % inlet_qualifier])
            self.units['X_good%s' % qualifier] = None
            self.filter_vals('X_good%s' % qualifier, 1.00, 'high')
            self.filter_vals('X_good%s' % qualifier, 0.00, 'low')

            #Standard conversion
            self['X_std%s' % qualifier] = self._calc_conversion(self['C_inlet%s' % inlet_qualifier]-self['CO2_inlet%s' % inlet_qualifier],sum([self['CO_outlet%s' % qualifier],self['CO2_outlet%s' % qualifier],self['CH4_outlet%s' % qualifier]])-self['CO2_inlet%s' % inlet_qualifier])
            self.units['X_std%s' % qualifier] = None
            self.filter_vals('X_std%s' % qualifier, 1.00, 'high')
            self.filter_vals('X_std%s' % qualifier, 0.00, 'low')

            #Total conversion to gas
            self['X_tot%s' % qualifier] = self._calc_conversion(self['C_inlet%s' % inlet_qualifier]-self['CO2_inlet%s' % inlet_qualifier],self['C_outlet%s' % qualifier]-self['CO2_inlet%s' % inlet_qualifier])
            self.units['X_tot%s' % qualifier] = None
            self.filter_vals('X_tot%s' % qualifier, 1.00, 'high')
            self.filter_vals('X_tot%s' % qualifier, 0.00, 'low')

        except KeyError:
            raise ConversionError, "The necessary columns to calculate CO yield have not been generated yet."

    def generate_normalized_compositions(self, name_qualifier = ""):
        """Generate normalized species from the list of process species"""
        norm_flow = sum([self["%s_outlet%s" % (spec, name_qualifier)] for spec in self.proc_species]) - sum([self["%s_outlet%s" % (spec, name_qualifier)] for spec in self.inert_species])

        for species in self.proc_species:
            self['%s_normalized' % species] = self["%s_outlet%s" % (species, "%s" % name_qualifier)]/norm_flow

    def calc_tar_rate(self, exit_stream, name_qualifier = "", tar_list = ['C6H6', 'C7H8', 'C10H8'], inclusive_tar_list = ['C2H2', 'C2H4', 'C2H6', 'C3H8', 'C3H6', 'C6H6', 'C7H8', 'C10H8'], norm_excl_names = ['Ar_MS', 'N2_MS']):
        """Calculate the level of tar in mg/Nm^3 exiting the gasifier"""
        #I'm just going to start by assuming that the outlet flowrate is molar -- I can make it more general later
	
        #Need to normalized tar to remove N2 and Ar
        factor = 1.0
        for norm_name in norm_excl_names:
            factor -= self[norm_name]/100.0
        outlet_vol_rate = exit_stream.flowrate[0] * 0.0224 * factor	#Nm^3/s, assuming mol/s for basis of original flowrate -- make it more general later

        total_tar = np.zeros(len(self.index))
        total_tar_incl = np.zeros(len(self.index))
        #total tar mass rate
        
        for molecule in tar_list:
            
            total_tar += Molecule(molecule).MW() * self['%s_outlet%s' % (molecule, name_qualifier)]

        for molecule in inclusive_tar_list:
            
            total_tar_incl += Molecule(molecule).MW() * self['%s_outlet%s' % (molecule, name_qualifier)]

        self['tar_loading'] = total_tar/outlet_vol_rate*1000.0
        self['tar_loading_incl'] = total_tar_incl/outlet_vol_rate*1000.0
        self.units['tar_loading'] = 'mg/m^3'
        self.units['tar_loading_incl'] = 'mg/m^3'

    def calc_space_time(self, reactor_vol, excluded_species, temp_method = 'fast_mean'):
        """Calculates the inlet space time of the reactor based on the inlet streams"""
        #Right now excluded_species must be in their own stream, find a way to fix this if I can... Maybe create temporary streams without the excluded species.
        conv = uc.UnitConverter()
        vol = conv.convert_units(reactor_vol[0], reactor_vol[1], 'm^3')
        
        excl_inlets = []
        for stream in self.inlet_streams:
            stream_flow = 0
            non_excl_species = {}
            for species in stream.composition:
                if species in excluded_species:
                    excl_inlets.append(stream)
        
        
        #need a mixer that works on both the FULL streams as well as the gas-only
        
        temp_inlets = [i for i in self.inlet_streams if i not in excl_inlets] 

        mix = Mixer('inlet_mix', inlets = temp_inlets, temp_method = temp_method)


        #mix.recalc()

        
        V_dot = mix.outlets[0].gas_volumetric_flowrate('m^3/s')
        self['volumetric_inlet_gas_only'] = V_dot
        self.units['volumetric_inlet_gas_only'] = 'm^3/s'

        self['T_cupmix_gas_only'] = mix.outlets[0].temperature[0]
        self.units['T_cupmix_gas_only'] = mix.outlets[0].temperature[1]

        tau = vol/V_dot
        #return tau

        self['space_time'] = tau
        self.units['space_time'] = 's'

    def calc_max_dH(self, temperature, pressure, units = 'W'):
        """Calculates the potential enthalpy change given the wall temperature and the feedrate of biomass"""
        #Create a Reactor from the inlet streams
        r = SDFIdealGasifier(y_CO2_out = 0.065, y_CH4_out = 0.0029, temperature = temperature, pressure = pressure, inlets = self.inlet_streams)
        r.generate_outlet_stream()
        self['dH_max'] = r.calc_enthalpy_change(units)
        self.units['dH_max'] = units
        #print self['dH_max'][0]

        #Calculate the Reactor's outlet stream using fixed mass balances (this could use an extensible mass balance package, but c'est la vie

        #Determine the potential change in enthalpy using built-in functions

    def calc_normalized_max_dH(self, tube_length, tube_diameter):
        """Calculates the dH_max/A measurement"""
        conv = uc.UnitConverter()
        A_lat = conv.convert_units(tube_length[0], tube_length[1], 'm') * conv.convert_units(tube_diameter[0], tube_length[1], 'm') * 3.14159
        self['dH_max/A'] = self['dH_max']/A_lat
        self.units['dH_max/A'] = "%s/m^2" % self['dH_max'].units

    def calc_min_RT(self, tube_length, tube_diameter, exit_temperature_tag, exit_pressure_tag, mode = "implicit"):
        """Calculates a minimum residence time -- needs a unit test!!!"""
        conv = uc.UnitConverter()
        V = conv.convert_units(tube_length[0], tube_length[1], 'm') * (conv.convert_units(tube_diameter[0],tube_diameter[1],'m'))**2.0*np.pi/4.0
        self['t_min'] = getattr(self, "_%s_minRT" % mode)(V, exit_temperature_tag, exit_pressure_tag)
        self.units['t_min'] = 's'

    def _implicit_minRT(self, V, T, P):
        #sum all of the of molar species
	conv = uc.UnitConverter()
        ndot = self.outlet_streams[0].flowrate[0]	#This is naturally in mol/s
        ndot -= self['H2O_outlet']			#subtract off the outlet water flow
	ndot += self['H2O_inlet']			#add on the inlet water flow
        Vdot = ndot * conv.convert_units(self[T], self.units[T], 'K')*8.314/conv.convert_units(self[P], self.units[P], 'Pa')
        return V/Vdot

    def _EERC_minRT(self, V, T, P):
        conv = uc.UnitConverter()
        ndot = self.outlet_streams[0].flowrate[0]
        ndot -= self['H2O_outlet']
        ndot += (1-self['H2O_MS']/100.0)/(1-self['ai_outlet_moisture']/100.0)*self['ai_outlet_moisture']/100.0*self.outlet_streams[0].flowrate[0]
        Vdot = ndot * conv.convert_units(self[T], self.units[T], 'K')*8.314/conv.convert_units(self[P],self.units[P], 'Pa')
        return V/Vdot

    def calc_dimensionless_numbers(self):
        """Calculates Re, Gr, Ri, ... at the inlet and outlet of the tube"""
	pass


    def generate_enthalpy_change(self, units):
        """Calculates the enthalpy change by multiplying the max change over the conversion"""
        conv = uc.UnitConverter()
        self['delta_H'] = conv.convert_units(self['dH_max'], self.units['dH_max'], units)*self['X_tot']


    def calc_min_residence_time(self):
        """Calculates the minimum bound on the residence time, assuming complete conversion and heat up at the instant materials enter the reactor"""
        pass


    def calc_min_residence_time(self, tube_length, tube_diameter):
        """Calculates the minimum bound on the residence time, assuming the conversion observed and heat up at the instant materials enter the reactor"""
        Vdot = self.outlet_streams[0].gas_volumetric_flowrate('m^3/s')
	conv = uc.UnitConverter()
        Vrx = conv.convert_units(tube_length[0], tube_length[1], 'm')*(conv.convert_units(tube_diameter[0],tube_diameter[1], 'm'))**2/4.0*3.14159
        self['t_min'] = Vrx/Vdot
        self.units['t_min'] = 's'

    def calc_optical_thickness(self, tubeD, density, particle_size):
        """Calculates the optical thickness of the inlet mixture.  Particle size should be a dictionary with d## keys"""
        conv = uc.UnitConverter()

        #Need to first get the solids flowrate in
        mdot = 0.0
	for stream in self.inlet_streams:
            for species in stream.composition:
                if species not in Stream.species or Stream.species[species]['phase'] == 's':  #OK, this is a huge stretch -- I'm assuming that if what we are looking for is not in the species dictionary, it is a solid
                    if stream.basis == 'mass':
                        mdot += conv.convert_units(stream.flowrate[0]*stream.composition[species],stream.flowrate[1],'kg/s')
                    #Need to add something for 'molar' or 'volume' later, but this will get it working

	rho = conv.convert_units(density[0], density[1], 'kg/m^3') #passed in for now, but should be able to add this as a material property later and get it through stream introspection
        Vdot = conv.convert_units(self['volumetric_inlet_gas_only'], self.units['volumetric_inlet_gas_only'], 'm^3/s')
        if tubeD is None:
            tubeD = [np.nan, 'm']
        if tubeD[0] is None:
            tubeD[0] = np.nan
        D = conv.convert_units(tubeD[0], tubeD[1], 'm')
        for dp in particle_size:
            if particle_size[dp][0] is None:
                particle_size[dp][0] = np.nan
            self['optical_thickness_%s' % dp] = 1.5*mdot/rho/conv.convert_units(particle_size[dp][0], particle_size[dp][1], 'm')/Vdot*D

    def generate_C_mass_balance(self):
        if 'C_inlet' not in self.columns or 'C_outlet' not in self.columns:
            raise NoInletOutletFlowrateError, "The inlet and outlet flowrates of carbon must be solved for first"
        self['C_gas_mass_balance'] = (self['C_inlet'] - self['C_outlet'])/self['C_inlet']
        self['C_gas_mass_balance'][np.logical_not(np.isfinite(self['C_gas_mass_balance']))] = np.nan

    def generate_CH4_yield(self):
        if 'C_inlet' not in self.columns or 'CH4_outlet' not in self.columns:
            raise NoInletOutletFlowrateError, "The inlet and outlet flowrates of carbon/CH4 must be solved for first"
        self['CH4_yield'] = (self['CH4_outlet']-self['CH4_inlet'])/(self['C_inlet']-self['CH4_inlet'] - self['CO2_inlet'])




class RunInformation:
    """Container class for information about experimental runs -- wrapper class around the dictionary object"""
    def __init__(self):
        self.info = {}
        self.start = None
        self.end = None

    def SQL_load(self, interface, table, run_id):
        """Load the data in from the given table into member objects"""
        if isinstance(run_id, int):
            query = SQL.select_Query(objects = ['*'], table = table, condition_list = ['run_id=%s' % run_id])
        else:
            query = SQL.select_Query(objects = ['*'], table = table, condition_list = ["run_id='%s'" % run_id])

        results = interface.query(query)

        #results will be a list of dicts
        for key, value in results[0].items():
            self.info[key] = value               #we'll use a glossary entry to convert the field names to variable names in the program, so I can choose whatever I want

    def __getitem__(self, index):
        try:
            return self.info[index]
        except KeyError:
            raise lflException, "The desired key is not in the run information container"

    def __setitem__(self, index, value):
        try:
            self.info[index] = value
        except KeyError:
            raise lflException, "The desired key is not in the run information container"
        
if __name__ == '__main__':

    pass
   
