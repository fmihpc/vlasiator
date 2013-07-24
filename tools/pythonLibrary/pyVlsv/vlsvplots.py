from vlsvreader import *
import pylab as pl
import numpy as np


def draw_plots_by_cellid(vlsvReader, variables1, variables2, cellids, coordinates=[], distances=[]):
   '''Draws a plot of given variables for the given cell ids
      :param vlsvReader    some VlsvFile with a file open
      :param variables1    some dictionary of variables to be plotted as the x-axis
      :param variables2    some dictionary of variables to be plotted as the y-axis
      :param cellids       some list of cellids
      :param coordinates   some list of coordinates for the cell ids (Only used if variables include the variable \"coordinates\")
      :param distances     some list of distances from the starting cell id (Only used if variables include the variable \"distances\")
   '''
   # Number of variables to be plotted in y-axis
   xplotnum = 0
   for var in variables1.iteritems():
      xplotnum = xplotnum + len(np.atleast_1d(var[1]))
   # Number of variables to be plotted in x-axis
   yplotnum = 0
   for var in variables2.iteritems():
      yplotnum = yplotnum + len(np.atleast_1d(var[1]))
   # Plot the variables:
   index = 1
   # Insert variables into the plot:
   # Note: variables1 is a python dictionary dict()
   for variable1 in variables1.iteritems():
      for variable2 in variables2.iteritems():
         # Read the variables:
         def read_plot_var(variable, cellids, coordinates, distances, vlsvReader):
            if variable != "CellID" and variable != "coordinates" and variable != "distance":
               # Get the variable from a file:
               variablelist = vlsvReader.read_variables_for_cellids(variable, cellids)
            elif variable == "CellID":
               variablelist = np.array(cellids)
            elif variable == "distance":
               variablelist = np.array(distances)
            elif variable == "coordinates":
               variablelist = np.array(coordinates)
            else:
               print "Bad variable1 name " + str(variable)
               variablelist = []
            return variablelist
         variable1list = read_plot_var(variable=variable1[0], cellids=cellids, coordinates=coordinates, distances=distances, vlsvReader=vlsvReader)
         variable2list = read_plot_var(variable=variable2[0], cellids=cellids, coordinates=coordinates, distances=distances, vlsvReader=vlsvReader)

         # Plot the variable:
         for i in np.atleast_1d(variable1[1]):
            # Get x plot coordinates:
            x = []
            for variableVector in variable1list:
               x.append(np.atleast_1d(variableVector)[i])
            for j in np.atleast_1d(variable2[1]):
               # Get y plot coordinates:
               y = []
               for variableVector in variable2list:
                  y.append(np.atleast_1d(variableVector)[j])
               # Set subplot
               pl.subplot(yplotnum,xplotnum,index)
               # Set labels:
               pl.xlabel(variable1[0] + "[" + str(i) + "]")
               pl.ylabel(variable2[0] + "[" + str(j) + "]")
               # Set x and y limits:
               xlength = max(x) - min(x)
               xmin = min(x) - 0.05*xlength
               xmax = max(x) + 0.05*xlength
               ylength = max(y) - min(y)
               ymin = min(y) - 0.05*ylength
               ymax = max(y) + 0.05*ylength
               pl.xlim((xmin, xmax))
               pl.ylim((ymin, ymax))
               # Plot the variable:
               pl.plot(x,y)
               #Increment the index:
               index = index + 1
   # Show the plot:
   pl.ion()
   pl.show()


def take_cut_through( fileName, variables1, variables2, point1, point2 ):
   '''Creates a plot of cut-through of some variable from a vlsv file
   :param fileName         Name of the vlsv file
   :param variables1       Variables to be plotted in the x-axis
   :param variables2       Variables to be plotted in the y-axis
   :param point1           The starting coordinates of the line
   :param point2           The ending coordinates of the line
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
   ###########################################################################################
   # Plot the variables:
   draw_plots_by_cellid(vlsvReader=vlsvReader, variables1=variables1, variables2=variables2, cellids=cellids, coordinates=coordinateList, distances=distances)
