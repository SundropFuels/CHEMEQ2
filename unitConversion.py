"""Unit conversion module for Python"""

import numpy as np

class UnitError(Exception):
    pass

class BadCharacterError(UnitError):
    pass

class EmptyUnitError(UnitError):
    pass

class UnitNotFoundError(UnitError):
    pass

class InconsistentUnitError(UnitError):
    pass

class UnitConverter:
    """This class converts arbitrary fundamental units from one to the other"""

    mass_dict = {'kg':1.0,'g':1000.0, 'lb':2.20462,'mg':1000000.0}
    time_dict = {'hr':1.0/3600.0, 'min':1.0/60.0, 's':1.0, 'day':1.0/(24.0*3600.0)}
    length_dict = {'m':1.0, 'in':39.37, 'mm':1000.0, 'cm':100.0, 'ft' : 39.37/12.0}
    mole_dict = {'mol':1.0, 'kmol':1.0/1000.0, 'SL':24.159}
    temperature_dict = {'K': 1.0, 'C': 1.0, 'R': 1.8, 'F': 1.8}
    money_dict = {'$':1.0}
    percent_dict = {'%':1.0, 'ppm':10000.0, 'fraction':1.0/100.0}

    energy_dict = {'J':1.0, 'kJ':0.001, 'MJ':1E-6, 'cal':0.239005736138, 'kcal':2.39005736138E-4, 'Btu':9.47817120313E-4, 'MMBtu':9.47817120313E-10}
    energy_units = 'kg*m^2/s^2'    

    pressure_dict = {'Pa':1, 'psia':1.4503773773E-4, 'psig':1.4503773773E-4, 'bar':1E-5, 'atm':9.86923266716E-6, 'mmHg':0.007500616827, 'mmH2O':0.101971621298, 'inH2O':0.00401463076}
    pressure_units = 'kg/s^2/m'

    volume_dict = {'L':1000.0, 'gal':264.172052358, 'mL':1000000.0}
    volume_units = 'm^3'

    power_dict = {'W':1, 'kW':0.001, 'MW':1E-6, 'mW':1000}
    power_units = 'kg*m^2/s^3'

    units = [mass_dict, time_dict, length_dict, mole_dict, temperature_dict, money_dict, percent_dict]
    derived_units = [(energy_units, energy_dict),(pressure_units,pressure_dict),(volume_units,volume_dict),(power_units,power_dict)]

    non_abs_temp_factors = {'K':0, 'C':-273.15, 'R':0, 'F':-459.67}
    temp_units = ['R', 'C', 'K', 'F']

    non_abs_pressure_factors = {'psig':0, 'psia':14.6959488}
    pressure_units = ['psig', 'psia']
    
    def __init__(self):
        self.unit_list = []

    
    def _derived_unit_replace(self, unit):
        """Takes derived units (N, J, etc.) and replaces them with new text representing the fundamental units"""
        #Search through the dictionary and see if the given unit is a derived unit
        #Return text for the fundamental unit  #THIS IS NOT WORKING YET!!!!!!!!!!!!
        if unit in UnitConverter.derived_units.keys():
            return UnitConverter.derived_units[unit]
        else:
            return unit

    def convert_units(self, val, from_str, to_str):
        from_parsed = self._parse_inputstr(from_str)
        to_parsed = self._parse_inputstr(to_str)
        
        if from_str == 'psig':
            val = val + (UnitConverter.non_abs_pressure_factors['psia'] - UnitConverter.non_abs_pressure_factors['psig'])
            
        if len(from_parsed) == 1:
            if from_parsed[0][0] in UnitConverter.temp_units and from_parsed[0][1] == 1:
                if to_parsed[0][0] in UnitConverter.temp_units and len(to_parsed) == 1 and to_parsed[0][1] == 1:
                    ret_val = val - UnitConverter.non_abs_temp_factors[from_parsed[0][0]]
                    ret_val = ret_val / UnitConverter.temperature_dict[from_parsed[0][0]]
                    ret_val = ret_val * UnitConverter.temperature_dict[to_parsed[0][0]]
                    ret_val = ret_val + UnitConverter.non_abs_temp_factors[to_parsed[0][0]]
                    return ret_val
                else:
                    raise InconsistentUnitError, "The units were not consistent for this temperature conversion"
            #Do nothing if it is not a direct temperature conversion; we could put gauge pressure conversions here, as well

        found = False

        consistency_list = [0] * len(UnitConverter.units)

        factor = 1.0
        for (unit, exponent) in from_parsed:
            
            
            #Check if the unit is in the derived_unit dictionaries first

            for (base_units, unit_dict) in UnitConverter.derived_units:
                for (key,value) in unit_dict.items():
                    if unit == key:
                       
                        factor = factor / np.power(value, float(exponent))
                        found = True
                        
                        parsed_sub = self._parse_inputstr(base_units)
                        for st, exponent2 in parsed_sub:
                            from_parsed.append((st, float(exponent2)*float(exponent)))

                if found == True:
                    break

            
            for unit_dict in UnitConverter.units:
                for (key, value) in unit_dict.items():
                    if unit == key:
                        factor = factor / np.power(value, float(exponent))
                        found = True
                        consistency_list[UnitConverter.units.index(unit_dict)] += exponent


                        break
                if found == True:
                    break

            if found == False:
                raise UnitNotFoundError, "The unit %s is not in the database" % unit
            found = False

        

        for (unit, exponent) in to_parsed:
            
            for (base_units, unit_dict) in UnitConverter.derived_units:
                for (key,value) in unit_dict.items():
                    if unit == key:
                        
                        factor = factor * np.power(value, float(exponent))
                        found = True
                        
                        parsed_sub = self._parse_inputstr(base_units)
                        for st, exponent2 in parsed_sub:
                            to_parsed.append((st, float(exponent2)*float(exponent)))

                if found == True:
                    break

            

            for unit_dict in UnitConverter.units:
                for (key, value) in unit_dict.items():
                    if unit == key:
                        factor = factor * np.power(value, float(exponent))
                        found = True
                        consistency_list[UnitConverter.units.index(unit_dict)] += -exponent
                        break
                if found == True:
                    break


            
            if found == False:
                raise UnitNotFoundError, "The unit %s is not in the database" % unit
            found = False

       

        for checksum in consistency_list:
            if checksum != 0:
                raise InconsistentUnitError, "The to and from units do not match: %s vs %s" % (from_str, to_str)

        #print consistency_list

        
        
               
        if to_str == 'psig':
            return val * float(factor) - (UnitConverter.non_abs_pressure_factors['psia'] - UnitConverter.non_abs_pressure_factors['psig'])
        else:
            
            return val * float(factor)

        

    def _parse_inputstr(self, i_string):
        """Parser for arbitrary symbolic mathematical expressions of units"""

        operators = ['*','/','^','(',')']
        number_sym = ['.', '-']
        parsed_units = []      
        sub_parsed = None
        last_char = ""
        current_unit = ""
        current_exp_str = ""
        current_exp = 1
        i = 0
        while i < len(i_string):
            current_char = i_string[i]
            if current_char not in operators:
                
                if not current_char.isalpha() and current_char != '$' and current_char != '%':
                    raise BadCharacterError, "Unit names can only contain alphabetical characters"
                else:
                    current_unit = "".join([current_unit, current_char])
                    
                    last_char = current_char
                    i += 1
                    
            else: 
                if current_char == '*' or current_char == '/':
                    
                    
                    if sub_parsed is not None:
                        for (u,e) in sub_parsed:
                            e = e * current_exp
                            parsed_units.append((u,e))
                        sub_parsed = None
                    else:
                        if current_unit == "":
                            raise EmptyUnitError, "Cannot add an empty unit"
                        parsed_units.append((current_unit, current_exp))
                    current_unit = ""
                    current_exp = 1
                    if current_char == '/':
                        current_exp = -1
                    i += 1
                
                elif current_char == '^':
                                       
                    while (i+1) < len(i_string) and (i_string[i+1] in number_sym or i_string[i+1].isdigit()):
                        i += 1
                        current_exp_str = "".join([current_exp_str, i_string[i]])
                    try:
                        current_exp = current_exp * float(current_exp_str)
                        current_exp_str = ""
                    except ValueError:
                        raise BadExponentError, "The exponent for %s was formatted incorrectly" % current_unit
                    if (i+1) < len(i_string) and i_string[i+1] not in operators:
                        #There are still more characters to come, so I want to make sure the next character is an operator
                        raise BadCharacterError, "An improper character followed the exponent string"
                    i += 1

                elif current_char == '(':
                    #get the substring
                    paren_level = 1
                    j = i
                    while paren_level > 0:
                        j += 1
                        if i_string[j] == '(':
                            paren_level += 1
                        elif i_string[j] == ')':
                            paren_level -= 1
                    
                    sub_parsed = self._parse_inputstr(i_string[i+1:j])
                    i = j + 1      #move to the next i
                    

        if sub_parsed is not None:
            for (u,e) in sub_parsed:
                e = e * current_exp
                parsed_units.append((u,e))
            sub_parsed = None
        else:
            if current_unit == "":
                raise EmptyUnitError, "Cannot add an empty unit"
            parsed_units.append((current_unit, current_exp))
        return parsed_units

if __name__ == '__main__':
   conv = UnitConverter()
   print "3 m^3 to L: %s" % conv.convert_units(3, 'm^3', 'L')
   

                
                    
