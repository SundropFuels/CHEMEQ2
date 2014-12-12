import numpy as np
import pandas as pd
import db_toolbox as SQL
import datetime
import csv
import unitConversion as uc

class dfException(Exception):
    def __init__(self,value):
        self.value = value

    def __str__(self):
        return repr(self.value)

class GlossaryError(Exception):
    pass

class BadGlossaryTagError(GlossaryError):
    pass

class BadGlossaryTypeError(GlossaryError):
    pass

class NoColumnError(dfException):
    pass

class BadArgumentError(dfException):
    pass

class UncertaintyError(dfException):
    pass

class BadUnitError(dfException):
    pass

class dfSQLError(dfException):
    pass

class Dataframe(pd.DataFrame):

    def __init__(self, data = None, units_dict = None, index = None, **kwargs):
        """Constructor puts together the basic data structure, which includes a dictionary of arrays -- this is essentially a subclassed pandas DataFrame now"""
        #should check whether **kwargs contains either data or index and drop those keys, so they are not called twice; also, I should write the tests first
        
        pd.DataFrame.__init__(self, data, index = index, **kwargs)
       
        #will only need to explicitly handle errors that pandas does NOT handle

        if units_dict is None:
            units_dict = {}

        for key in units_dict.keys():
            if key not in self.columns:
                raise NoColumnError, "Unit given for non-existent column"

        for value in units_dict.values():
            if type(value) != str and value is not None:
                raise BadUnitError, "Units must be given as strings"

        self.units = units_dict.copy()
   
    def _set_item(self, key, value):
        """Overriding the inherent set item function to allow addition of a column when the index is empty"""
        if len(self) == 0:
            #call your reinitialize function on the dictionary
            self.reinitialize_data(data = {key:value})
        else:
            pd.DataFrame._set_item(self, key, value)


    def val_units(self, col_name, units = None):
        if col_name not in self.columns:
            raise NoColumnError, "The specified column is not in the dataframe"

        if units is None:
            units = self.units[col_name]
        conv = uc.UnitConverter()
        return (conv.convert_units(self[col_name],self.units[col_name],units), units)


    def SQL_load_data(self, db_interface, table = "", conditions = None):
        """Allows the user to load the dataframe directly from a MySQL database; must pass interface"""
        if conditions is None:
            conditions = []
        if not isinstance(db_interface, SQL.db_interface):
            raise dfSQLError, "The passed interface is not a db_toolbox interface"

        query = SQL.select_Query(objects = ['*'], table = table, condition_list = conditions)
        

        #try:
        #This should now return a pandas Dataframe object, from which we can initialize our own object
        s = db_interface.query(query, return_type = 'pandas_dataframe')
           
        #except SQL.DBToolboxError:
        #    raise dfSQLError, "There was an error in using the SQL interface"

                
        #try this:
        self.reinitialize_data(s)
      

        
       
    def reinitialize_data(self, data):
        
        self.__init__(data = data, units_dict = self.units)


    def SQL_upload_data(self, db_interface, table = "", conditions = None, index_col = "ts", max_query_length = 2000):
        """Allows the user to upload data to a MySQL database; tries an INSERT command first, then an UPDATE if an error is received"""
        if conditions is None:
            conditions = []

        if not isinstance(db_interface, SQL.db_interface):
            raise dfSQLError, "The passed interface is not a db_toolbox interface"


        #Because I INSERT the entire table anyway when there is a conflict, instead of UPDATE we will just DELETE all of the rows to INSERT
	#As this is conditional in SQL, it will fail silently when the row is not in the database (or should)

        #build the DELETE query
        delete_condition = ""
        for index in range(0, len(self)):
            delete_condition += "%s = '%s' OR " % (index_col, self[index_col][index])
        delete_condition = delete_condition[:-4]

        q = SQL.delete_Query(table = table, conditions = [delete_condition,])
        #print q.getQuery()[:1000]
        db_interface.query(q)
       
        #Now build the INSERT query, broken up into max query lengths
        num_ranges = int(len(self)/max_query_length)
        
        rmin = [i for i in range(0,num_ranges*max_query_length,max_query_length)]
        rmax = [i for i in range(max_query_length, (num_ranges+1)*max_query_length,max_query_length)]

        if num_ranges == 0:
            rmin = [0]
            rmax = [len(self)]
        
       
        if len(self) > max_query_length:
            rmin.append(num_ranges*max_query_length)
            rmax.append(len(self))
        #This is a bit cludgy, but it is what I am doing
        
        q_vals_set = {}
        for i in rmin:  #just need the right number of qvals in the set
            key_map = {}
            q_vals = {}
            for key in self.columns:
                if ' ' in key:
                    key_map[key] = '`' + key + '`'
                else:
                    key_map[key] = key
                q_vals[key_map[key]] = []
            q_vals_set[i] = q_vals             #use a dict for easy access later
            

        i = 0
        for (range_min, range_max) in zip(rmin, rmax):
            print "uploading pass %s." % (i+1)
            for index in range(range_min,range_max):
            
            
                for key in self.columns:
                
                    #print "key: %s\tvalue: %s" % (key, self[key][index])
                    q_vals_set[range_min][key_map[key]].append("'" + str(self[key][index]) + "'")
                    if q_vals_set[range_min][key_map[key]][index-max_query_length*i] == "'nan'":
                        q_vals_set[range_min][key_map[key]][index-max_query_length*i] = 'NULL'
                
            try:
                query = SQL.multiple_insert_Query(objects = q_vals_set[range_min], table = table)
                
                #f1 = open('sql_output.txt', 'w')
                #f1.write(query.getQuery())
                #f1.close()
                db_interface.query(query)
                query = SQL.commit_Query()
                db_interface.query(query)

            except SQL.DBToolboxError:
                raise dfSQLError, "Whoa...can't load that into the table, brah.  Don't know why. %s"
            i+=1

    def write_csv(self, filename, mode = "new", col_list = None):
        """Writes the dataframe to a csv file"""
        mode_dict = {"new":"w", "append":"a"}
        if not isinstance(filename, str):
            raise dfException, "The filename must be a string"
        outputfile = open(filename, mode_dict[mode])
        writer = csv.writer(outputfile, delimiter = ',', quoting=csv.QUOTE_MINIMAL)
        if col_list is not None:
            for key in col_list:
                if key not in self.columns:
                    raise dfException, "%s is not a column in the dataframe -- cannot write" % key
            head = col_list
        else:
            head = self.columns
        writer.writerow(head)
        for index in range(0,len(self)):
            row = []
            for key in head:
                row.append(self[key][index])
            writer.writerow(row)
        outputfile.close()

    def glossary_replace(self, glossary = None):
        """Replace the names of the keys in the data frame with new keys in the glossary"""
        if glossary == None:
            glossary = {}

        if type(glossary) != dict:
            raise BadGlossaryTypeError, "The glossary must be a dictionary"

        


        for key in glossary.keys():
            if key not in self.columns:
                #print key
                raise BadGlossaryTagError, "There is a tag in the glossary that is not in the dataframe -- %s" % key
            


            if key != glossary[key]:
                self[glossary[key]] = self[key]
                del self[key]

    def set_units(self, units_dict = None):
        if units_dict == None:
            units_dict = {}

        for key in units_dict.keys():
            if key not in self.columns:
                raise dfException, "Unit given for non-existent column: %s" % key

        for value in units_dict.values():
            if type(value) != str and value is not None:
                raise dfException, "Units must be given as strings"

        self.units = units_dict.copy()

    def filter_vals(self, column, v, action):
        
        if type(column) != str or not isinstance(v, (int, long, complex, float)) or type(action) != str:
            raise BadArgumentError, "The filter function received an argument not consistent with its parameter list"

        if action not in ('high', 'low'):
            raise BadArgumentError, "The action for a filter function must be either 'high' or 'low'"

        if column not in self.columns:
            raise NoColumnError, "The specified column is not in the dataframe"

        ret_val = self[column].copy()












        if action == 'high':
            ret_val[ret_val>v] = np.nan
        else:
            ret_val[ret_val<v] = np.nan

        return ret_val

    ###HERE -- there are specific functions to do this within pandas, I should rely on those
    #def replace_None_w_nan(self, column):
    #    if column not in self.columns:
    #        raise NoColumnError, "The specified column is not in the dataframe"
    #    
    #    if self[column].dtype != 'int64':
    #        
    #        self[column][np.equal(self[column],None)] = np.nan
    def replace_None_w_nan(self, column):
        if column not in self.columns:
            raise NoColumnError, "The specified column is not in the dataframe"
        
        if self[column].dtype != 'int8' and self[column].dtype != 'int16' and self[column].dtype != 'int32' and self[column].dtype != 'int64':
            
            self[column][self[column].isnull()] = np.nan


    def replace_None_w_nan_all(self):
        for column in self.columns:
            
            self.replace_None_w_nan(column)


    def finite_vals(self, key):
        """Returns only the finite values -- will build on this later to add the ability to return indices as well, to allow for pairs indexed on finite only"""
        if key not in self.columns:
            raise NoColumnError, "The specified column is not in the dataframe"
        return (self[key][self[key].notnull()]).values

    def finite_set(self, key, cols = None):
        """Returns a dataframe built up from only the finite values for key"""
        array_dict = {}
        array_dict[key] = (self[key][self[key].notnull()]).values
        if cols is not None:
            for col in cols:
                array_dict[col] = (self[col][self[key].notnull()]).values
        return Dataframe(data = array_dict)

    


"""
if __name__ == "__main__":
    print "testing this object -- really need to add unit tests, you know"
    frame = Dataframe()
    interface = SQL.db_interface(host = "localhost", user = "root", passwd = "udflyer87")
    interface.connect()
    q = SQL.use_Query('trial')
    interface.query(q)
    frame.SQL_load_data(interface, table = "lv_data")
    print frame['TE_1']
"""

if __name__ == "__main__":
    print "testing the distribution object"
    normal = NormalDistribution(0,1)
    print "Value of cumulative probability halfway: %s" % normal.cum[5000,1]
    (left1, right1) = binary_search_1D(normal.cum[:,0],-1.0)
    (left2, right2) = binary_search_1D(normal.cum[:,0],1.0)
  
    print normal.cum[left2,1] - normal.cum[left1,1]
    print normal.sample_N(20)


