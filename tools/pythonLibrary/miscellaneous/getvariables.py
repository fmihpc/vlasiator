from vlsvreader import *

def get_variables(variable, fileNames, cellids):
   #Create a dictionary for holding cell ids
   cellidvariables = {}
   for i in cellids:
      cellidvariables[i] = []
   #Read variables
   for f in fileNames:
      vlsvReader = VlsvFile(f)
      variables = vlsvReader.read_variables_for_cellids(name=variable, cellids=cellids)
      for i in xrange(len(cellids)):
         cellidvariables[cellids[i]].append(variables[i])
   #Return cell id variables
   return cellidvariables
