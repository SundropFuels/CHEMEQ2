"""Database toolbox for working with mySQL
Classes and methods for using the mySQLdb API

Chris Perkins
10-31-2010

Version 0.1 - Initial implementations

"""
#! /usr/bin/env python
import MySQLdb
import pandas.io.sql as pd_sql

class DBToolboxError(Exception):
    pass

class NoDBInUseError(DBToolboxError):
    pass

class InterfaceNotConnectedError(DBToolboxError):
    pass

class db_interface:
    
    def __init__(self, host = "", user = "", passwd = ""):
        self.parameters = {}
        self['host'] = host
        self['user'] = user
        self['passwd'] = passwd
        self.connected = False
        self.active_db = None        


    def connect(self):
        self.connection = MySQLdb.connect(host = self['host'], user = self['user'], passwd = self['passwd'])
        self.cursor = self.connection.cursor(MySQLdb.cursors.DictCursor)
        
        self.connected = True
        #self.connection.autocommit(True)

    def __setitem__(self,key,item):
        self.parameters[key] = item

    def __getitem__(self,key):
        return self.parameters[key]

    def query(self, myQuery, return_type = 'default'):
        
        if self.connected == False:
            raise InterfaceNotConnectedError, "Interface not connected to a database server"        

        if self.active_db == None and not (isinstance(myQuery, show_Query) or isinstance(myQuery, use_Query)):
            raise NoDBInUseError, "You cannot query a database that is not in use"
        if return_type == 'default' or return_type == 'dict_of_dicts':
            try:
                self.cursor.execute(myQuery.getQuery())
                result = self.cursor.fetchall()
            except MySQLdb.Error, e:
                #print "This is the error I got: %s" % e
                raise DBToolboxError
        elif return_type == 'pandas_dataframe':
            try:
                #call the pandas library to get a dataframe from the request
                result = pd_sql.read_frame(myQuery.getQuery(), self.connection)
                
            except MySQLdb.Error, e:
                #print "This is the error I got: %s" % e
                raise DBToolboxError
        
        if isinstance(myQuery, use_Query):
            self.active_db = myQuery.getQuery().lstrip('USE ').rstrip(';')

                    
        return result

    def close(self):
        self.connection.close()
        self.connected = False


class Query:
    def __init__(self):
        self.query = []
        self.items = {}
        
    def getQuery(self):
                
        q =  ' '.join(e for e in self.query)
        return q + ';'

    def setObjects(self, obj_string):
        self['objects'] = obj_string

    def __setitem__(self, key, item):
        self.items[key] = item

    def __getitem__(self, key):
        return self.items[key]

class oneWord_Query(Query):
    def __init__(self):
        Query.__init__(self)

    def getQuery(self):
        try:
            self.query = []
            self.query.extend([self['command']])
        except KeyError:
            print "Error - Missing Query Key"
        return Query.getQuery(self)

class commit_Query(oneWord_Query):
    def __init__(self):
        oneWord_Query.__init__(self)
        self['command'] = 'COMMIT'

class short_Query(Query):
    def __init__(self, objects = []):
        Query.__init__(self)
        self['objects'] = objects

    def getQuery(self):
        try:
            self.query = []
            self.query.extend([self['command'],self['objects']])
        except KeyError:
            print "Error - Missing Query Key"
        return Query.getQuery(self)

class use_Query(short_Query):
    def __init__(self, objects = []):
        short_Query.__init__(self, objects = objects)
        self['command'] = "USE"
        
class show_Query(short_Query):
    def __init__(self, objects = []):
        short_Query.__init__(self, objects = objects)
        self['command'] = "SHOW"

class describe_Query(short_Query):
    def __init__(self, objects = []):
        short_Query.__init__(self, objects = objects)
        self['command'] = "DESCRIBE"

class set_Query(Query):
    """
        Set a database variable
    """
    def __init__(self,variable,objects = []):
        Query.__init__(self)
        self['command'] = ["SELECT",":="]
        self.setObjects(objects)
        self.setVariable(variable)
    
    def setVariable(self,variable):
        self['variable'] = variable
        
    def getQuery(self):
        try:
            self.query.extend([self['command'][0],self['variable'],self['command'][1],self['objects']])
            #print self.query
        except KeyError:
            print "Error - Missing Query Key"
        return Query.getQuery(self)   
     
class call_Query(Query):
    """
        Call a stored procedure
    """
    def __init__(self,procedure,objects = []):
        Query.__init__(self)
        self['command'] = "CALL"
        self.setObjects(objects)
        self.setProcedure(procedure)
        
    def setObjects(self, objects):
        self['objects'] =  "(" + ",".join(obj for obj in objects) + ")"  
    
    def setProcedure(self,procedure):
        self['procedure'] = procedure
    
    def getQuery(self):
        try:
            self.query.extend([self['command'],self['procedure'],self['objects']])
        except KeyError:
            print "Error - Missing Query Key"
        return Query.getQuery(self)       
     
class select_Query(Query):
    def __init__(self, objects = [], table = "", condition_list = []):
        Query.__init__(self)
        
        self['command'] = "SELECT"
        self.setObjects(objects)
        self.setTable(table)
        self.setConditions(condition_list)


    def setObjects(self, obj_list):
        self['objects'] =  ",".join(obj for obj in obj_list) 
           
    def setTable(self, table):
        self['table'] = "FROM " + table

    def setConditions(self, condition_list):
        if condition_list != []:
            self['conditions'] = "WHERE " + " AND ".join(cond for cond in condition_list)

    def getQuery(self):
        """Stiches together query from conditions"""
        self.query = []
        try:
            self.query.extend([self['command'],self['objects'],self['table']])
        except KeyError:
            print "Error - Missing Query Key"
        try:
            self.query.append(self['conditions'])
        except KeyError:
            pass
        return Query.getQuery(self)

class insert_Query(Query):
    def __init__(self, objects = {}, table = ""):
        Query.__init__(self)
        self['command'] = "INSERT INTO"
        self.setObjects(objects)
        self.setTable(table)
        


    def setObjects(self, objects):
        self['objects'] =  "(" + ",".join(obj for obj in objects.keys()) + ")"
        #self['values'] = "VALUES ('" + "','".join(val for val in objects.values()) + "')"
        
        vals = []
        for v in objects.values():
            if v == 'NULL':
                vals.append(v)
            else:
                vals.append("'" + v + "'")

        self['values'] = "VALUES (" + ",".join(val for val in vals) + ")"
           
    def setTable(self, table):
        self['table'] = table

           

    def getQuery(self):
        """Stiches together query from conditions"""
        self.query = []
        try:
            self.query.extend([self['command'],self['table'],self['objects'],self['values']])
        except KeyError:
            print "Error - Missing Query Key"
            
        
        return Query.getQuery(self)

class multiple_insert_Query(insert_Query):
    def __init__(self, **kwargs):
        insert_Query.__init__(self, **kwargs)

    def setObjects(self, objects):
        #make sure the objects are the right length and type
        if not isinstance(objects, dict):
            raise DBToolboxError, "objects must be a dictionary of lists"
        for key in objects:
            if not isinstance(objects[key], list):
                raise DBToolboxError, "objects must be a dictionary of lists"

        l = len(objects[objects.keys()[0]])
        for key in objects:
            if len(objects[key]) != l:
                raise DBToolBoxError, "The length of the lists in objects for key %s is inconsistent with the first item" % key

        

        self['objects'] = "(" + ",".join(obj for obj in objects.keys()) + ")"
        #now need to create the lists of specific rows to join
        vals = []
        for i in range(0, l):
            row = []
            for key in objects:
                row.append(str(objects[key][i]))
            row = "(" + ",".join(row) + ")"
            vals.append(row)
        self['values'] = "VALUES " + ",".join(vals) 


class insert_udk_Query(insert_Query):
    "On duplicate key update type insert"
    def __init__(self, **kwargs):
        insert_Query.__init__(self, **kwargs)

    def getQuery(self):
        """Stiches together query from conditions"""
        a = insert_Query.getQuery(self)   #just stowing the response, not needed
        self.query.append('ON DUPLICATE KEY UPDATE')
        return Query.getQuery(self)

class delete_Query(Query):
    def __init__(self, table = "", conditions = []):
        Query.__init__(self)
        self['command'] = "DELETE FROM"
        self.setTable(table)
        self.setConditions(conditions)

    def setTable(self, table):
        self['table'] = table

    def setConditions(self, conditions):
        self['conditions'] = "WHERE " + " AND ".join(cond for cond in conditions)
            
    def getQuery(self):
        """Stiches together query from conditions"""
        try:
            self.query.extend([self['command'],self['table'],self['conditions']])
        except KeyError:
            print "Error - Missing Query Key"
        return Query.getQuery(self)

class replace_Query(Query):
    #If no primary key added, this works exactly like INSERT
    def __init__(self, objects = {}, table = ""):
        Query.__init__(self)
        self['command'] = "REPLACE INTO"
        self.setObjects(objects)
        self.setTable(table)
        


    def setObjects(self, objects):
        self['objects'] =  "(" + ",".join(obj for obj in objects.keys()) + ")"
        #self['values'] = "VALUES ('" + "','".join(val for val in objects.values()) + "')" 
        vals = []
        for v in objects.values():
            if v == 'NULL':
                vals.append(v)
            else:
                vals.append("'" + v + "'")

        self['values'] = "VALUES (" + ",".join(val for val in vals) + ")"   

    def setTable(self, table):
        self['table'] = table

           

    def getQuery(self):
        """Stiches together query from conditions"""
        self.query = []
        try:
            self.query.extend([self['command'],self['table'],self['objects'],self['values']])
        except KeyError:
            print "Error - Missing Query Key"          
        
        return Query.getQuery(self)

class update_Query(Query):
    def __init__(self,objects = {},table = "",conditions = []):
        Query.__init__(self)
        self['command'] = ["UPDATE","SET"]
        self.setObjects(objects)
        self.setTable(table)
        self.setConditions(conditions)
        
    def setObjects(self, objects):
        vals = {}
        for k,v in objects.items():
            if v == 'NULL':
                vals[k] = v
            else:
                vals[k] = ("'" + v + "'")

        
        self['objects'] =  ",".join("%s = %s"%(k,v) for k,v in vals.items())
        
    def setTable(self, table):
        self['table'] = table

    def setConditions(self, conditions):
        if conditions != []:
            self['conditions'] = "WHERE " + " AND ".join(cond for cond in conditions)
            
    def getQuery(self):
        """Stiches together query from conditions"""
        self.query = []
        try:
            self.query.extend([self['command'][0],self['table'],self['command'][1],self['objects'],self['conditions']])
        except KeyError:
            print "Error - Missing Query Key"
            
        
        return Query.getQuery(self)        
        
        
if __name__ == '__main__':
    interface = db_interface()
    interface['host'] = 'localhost'
    interface['user'] = 'root'
    interface['passwd'] = 'Sparticus'        

    use_query = use_Query('crashcourse')
    use_querytext = use_query.getQuery()
    print use_querytext
    
    interface.connect()
    interface.query(use_query)
    
    #Select Query
    selectQuery = select_Query(objects = ['*'], table = 'orders',condition_list = ['order_num = %s'%'20005'])
    results = interface.query(selectQuery)
    #print results
#    for row in results:
#        for (k,v) in row.items():
#            print k,v
            
    #Update Query
    updateQuery = update_Query(objects = {'customer_state':'CO','customer_city':'Longmont'},table = 'customers',
                               conditions = ['customer_name = Jack'])
    
    
    #Call Query
    callQuery = call_Query(procedure = 'add_duplicate_samples',objects = ['2','1','1'])
    setQuery = set_Query(variable = '@last_run_id',objects = 'last_insert_id()')
    #print selectQuery.getQuery()
    #print updateQuery.getQuery()
    #print callQuery.getQuery()
