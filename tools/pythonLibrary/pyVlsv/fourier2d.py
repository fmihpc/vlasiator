import numpy as np
from getcellid import *

def get_cell_indices( bounds, cellid ):
   cellid = (int)(cellid - 1)
   # Get xmin, xmax, etc
   xmin = bounds[0]
   ymin = bounds[0 + 3]
   zmin = bounds[0 + 3 + 3]
   
   xmax = bounds[1]
   ymax = bounds[1 + 3]
   zmax = bounds[1 + 3 + 3]

   xcells = bounds[2]
   ycells = bounds[2 + 3]
   zcells = bounds[2 + 3 + 3]

   # Get cell lengths:
   cell_lengths = np.array([(xmax - xmin)/(float)(xcells), (ymax - ymin)/(float)(ycells), (zmax - zmin)/(float)(zcells)])

   cellindices = np.zeros(3)

   cellindices[0] = cellid%(int)(xcells)
   cellindices[1] = (cellid%ycells)/(int)(xcells)
   cellindices[2] = cellid/(int)(xcells*ycells)

   return cellindices

#def get_cell_id_from_indices( bounds, cellindices ):
#   # Get xmin, xmax, etc
#   xmin = bounds[0]
#   ymin = bounds[0 + 3]
#   zmin = bounds[0 + 3 + 3]
#   
#   xmax = bounds[1]
#   ymax = bounds[1 + 3]
#   zmax = bounds[1 + 3 + 3]
#
#   xcells = bounds[2]
#   ycells = bounds[2 + 3]
#   zcells = bounds[2 + 3 + 3]
#
#   return cellindices[0] + cellindices[1]
   

#def get_3d_cellids( bounds, BBOX ):
#   # Get xmin, xmax, etc
#   xmin = bounds[0]
#   ymin = bounds[0 + 3]
#   zmin = bounds[0 + 3 + 3]
#   
#   xmax = bounds[1]
#   ymax = bounds[1 + 3]
#   zmax = bounds[1 + 3 + 3]
#
#   xcells = bounds[2]
#   ycells = bounds[2 + 3]
#   zcells = bounds[2 + 3 + 3]
#
#   # Get cell lengths:
#   cell_lengths = np.array([(xmax - xmin)/(float)(xcells), (ymax - ymin)/(float)(ycells), (zmax - zmin)/(float)(zcells)])
#
#   minCellId = get_cell_id(bounds, [BBOX[0], BBOX[2], BBOX[4]])
#   maxCellId = 

def get_2d_array( fileNames, variables, BBOX ):
   '''Fetches the cell's data from filenames
      :param fileNames         List of the files
      :param variables         list of variables
      :param BBOX              boundary box, (=[xmin, xmax, ymin, ymax, zmin, zmax])
      :returns array of variables
   '''
   cellids = []

   fileNames = np.atleast_1d(fileNames)

   # Get the cell ids inside the boundary box:
   vlsvReader = VlsvFile(fileNames[0])
   # Get xmax, xmin and xcells_ini
   xmax = vlsvReader.read_parameter(name="xmax")
   xmin = vlsvReader.read_parameter(name="xmin")
   xcells = vlsvReader.read_parameter(name="xcells_ini")
   # Do the same for y
   ymax = vlsvReader.read_parameter(name="ymax")
   ymin = vlsvReader.read_parameter(name="ymin")
   ycells = vlsvReader.read_parameter(name="ycells_ini")
   # And for z
   zmax = vlsvReader.read_parameter(name="zmax")
   zmin = vlsvReader.read_parameter(name="zmin")
   zcells = vlsvReader.read_parameter(name="zcells_ini")

   # Get bounds for reading cell ids from coordinates
   bounds = [xmin, xmax, xcells, ymin, ymax, ycells, zmin, zmax, zcells]

   # Get cell lengths
   cell_lengths = np.array([(xmax - xmin)/(float)(xcells), (ymax - ymin)/(float)(ycells), (zmax - zmin)/(float)(zcells)])

   array = {}
   #Get cell ids:
   #Iterate through the BBOX coordinates:
   minCoordinates = [BBOX[0], BBOX[2], BBOX[4]]
   maxCoordinates = [BBOX[1], BBOX[3], BBOX[5]]
   coordinates = minCoordinates

   xindex = 0
   yindex = 0
   zindex = 0

   while coordinates[2] <= maxCoordinates[2]:
      while coordinates[1] <= maxCoordinates[1]:
         while coordinates[0] <= maxCoordinates[0]:
            cellid = get_cell_id( bounds, coordinates )
            cellids.append(cellid)
            # Move to the next cell in x direction
            coordinates[0] = coordinates[0] + cell_lengths[0]
         # Move to the next cell in y direction
         coordinates[1] = coordinates[1] + cell_lengths[1]
         coordinates[0] = minCoordinates[0]
      # Move to the next cell in z direction
      coordinates[2] = coordinates[2] + cell_lengths[2]
      coordinates[1] = minCoordinates[1]
      coordinates[0] = minCoordinates[0]
   # Create an array for holding variables:
   maxIndices = get_cell_indices( bounds, max(cellids) )
   minIndices = get_cell_indices( bounds, min(cellids) )
   if (maxIndices[0] - minIndices[0])*(maxIndices[1] - minIndices[1]) != len(cellids):
      print "BAD CELLIDS"

   # Create an array for holding variables
   two_d_variables = []
   for i in xrange(maxIndices[1] - minIndices[1]):
      two_d_variables.append(np.zeros((maxIndices[0] - minIndices[0])))
   two_d_variables = np.array(two_d_variables)

   # Input positions of the cell ids
   positions = {}
   for i in cellids:
      positions[i] = get_cell_indices( bounds, i ) - minIndices

   variable_dict = {}
   for i in variables:
      variable_dict[i] = []

   # Read variables:
   for f in fileNames:
      # Open file
      vlsvReader = VlsvFile(f)
      # Read variables:
      for i in variables:
         # Read in the variable
         variableArray = vlsvReader.read_variables_for_cellids(i, cellids)
         # Input into correct positions:
         for j in xrange(len(cellids)):
            position = positions[cellids[j]]
            two_d_variables[position[1]][position[0]] = variableArray[j]
         # Input variables:
         variable_dict[i].append(np.array(two_d_variables))
   # Return variables:
   return np.array(variable_dict)

def fourier_2d_array():
   print "test"


















