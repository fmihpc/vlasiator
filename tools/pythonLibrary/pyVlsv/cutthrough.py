from vlsvreader import *
import numpy as np

def cut_through(fileName, point1, point2):
   '''Creates an array of cut-through of some variable from a vlsv file
   :param fileName         Name of the vlsv file
   :param point1           The starting coordinates of the line
   :param point2           The ending coordinates of the line
   :returns                Array [cellids, coordinates, distances] of the cut through
   '''
   # Transform point1 and point2 into numpy array:
   point1 = np.array(point1)
   point2 = np.array(point2)
   # Open the file
   vlsvReader = VlsvFile(fileName)
   # Get parameters from the file to determine a good length between points (step length):
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
   #Calculate cell lengths:
   cell_lengths = np.array([(xmax - xmin)/(float)(xcells), (ymax - ymin)/(float)(ycells), (zmax - zmin)/(float)(zcells)])

   # Get coordinates between the points:
   coordinateList = []
   unitVec = (point2 - point1)/np.linalg.norm(point2 - point1)
   # Get a vector of cell length's length in the direction of (point2-point1)
   oneVec = unitVec * np.min(cell_lengths)
   # Get the number of points:
   numberOfPoints = (int)(np.linalg.norm(point2 - point1) / np.linalg.norm(oneVec)) + 1
   # Some routine checks:
   if numberOfPoints <= 2:
      print "ERROR, TOO SMALL A DISTANCE BETWEEN POINT1 AND POINT2"
      return
   # Input coordinates:
   for i in xrange((int)(np.linalg.norm(point2 - point1) / np.linalg.norm(oneVec)) + 1):
      coordinateList.append( point1 + oneVec*i )
   # Get the cell ids at the coordinates:
   cellids = []
   for coordinate in coordinateList:
      # Get cell indices:
      cellindices = np.array([(int)((coordinate[0] - xmin)/(float)(cell_lengths[0])), (int)((coordinate[1] - ymin)/(float)(cell_lengths[1])), (int)((coordinate[2] - zmin)/(float)(cell_lengths[2]))])
      # Get the cell id:
      cellid = cellindices[0] + cellindices[1] * xcells + cellindices[2] * xcells * ycells + 1
      # Append the cell id:
      cellids.append(cellid)
      # Get the distances to coordinates from the starting point:
      distances = []
      for i in coordinateList:
         # Get the distance
         distances.append(np.linalg.norm(i - point1))
   # Get the variables:
   
   # Return the cut through:
   return [cellids, coordinateList, distances]

