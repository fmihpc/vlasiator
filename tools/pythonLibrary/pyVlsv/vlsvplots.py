from vlsvreader import *
import pylab as pl
import numpy as np
import scipy as sc
from scipy import optimize
from plotfunctions import *
from fourierplot import *

class nullfit:
   nulldata=[]

#def draw_by_cellid_array(vlsvReader, variables1, variables2, cellids, coordinates=[], distances=[], fitFunction=nullfit):
#   '''Returns variables for the given cell ids as well as fitted function if specified
#      :param vlsvReader    some VlsvFile with a file open
#      :param variables1    some dictionary of variables to be plotted as the x-axis
#      :param variables2    some dictionary of variables to be plotted as the y-axis
#      :param cellids       some list of cellids
#      :param coordinates   some list of coordinates for the cell ids (Only used if variables include the variable \"coordinates\")
#      :param distances     some list of distances from the starting cell id (Only used if variables include the variable \"distances\")
#      :param fitFunction   some function for making a fit for the x, y plot
#      :returns             list of variables fetched from the vlsvReader as well as the fitted variables in dictionary form: np.array([variables1list, variables2list, variables1listfitted, variables2listfitted])
#   '''
#   if fitFunction != nullfit:
#      # Get function:
#      fetchedFunction = False
#      global functionList
#      # Try to get the function by name first
#      for i in functionList.iteritems():
#         if fitFunction == i[0]:
#            functionClass = functionList[fitFunction]
#            fetchedFunction = True
#      # If the user didn't fetch function by name, then it's simply the class itself:
#      if fetchedFunction == False:
#         functionClass = fitFunction
#      function = functionClass.get_function()
#      parameters = functionClass.get_parameters()
#
#   # Number of variables in y-axis
#   xplotnum = 0
#   for var in variables1.iteritems():
#      xplotnum = xplotnum + len(np.atleast_1d(var[1]))
#   # Number of variables in x-axis
#   yplotnum = 0
#   for var in variables2.iteritems():
#      yplotnum = yplotnum + len(np.atleast_1d(var[1]))
#   index = 1
#   # Insert variables:
#   variable1list = dict()
#   variable2list = dict()
#   variable1listfitted = dict()
#   variable2listfitted = dict()
#   for i in variables1.iteritems():
#      variable1list[i[0]] = []
#      variable1listfitted[i[0]] = []
#      for j in variables2.iteritems():
#         variable2list[j[0]] = []
#         variable2listfitted[j[0]] = []
#   def read_plot_var(variable, cellids, coordinates, distances, vlsvReader):
#      if variable != "CellID" and variable != "coordinates" and variable != "distance":
#         # Get the variable from a file:
#         variablelist = vlsvReader.read_variables_for_cellids(variable, cellids)
#      elif variable == "CellID":
#         variablelist = np.array(cellids)
#      elif variable == "distance":
#         variablelist = np.array(distances)
#      elif variable == "coordinates":
#         variablelist = np.array(coordinates)
#      else:
#         print "Bad variable1 name " + str(variable)
#         variablelist = []
#      return variablelist
#   # Note: variables1 is a python dictionary dict()
#   for variable1 in variables1.iteritems():
#      variable1list[variable1[0]].append(read_plot_var(variable=variable1[0], cellids=cellids, coordinates=coordinates, distances=distances, vlsvReader=vlsvReader))
#      for variable2 in variables2.iteritems():
#         # Read the variables:
#         variable2list[variable2[0]].append(read_plot_var(variable=variable2[0], cellids=cellids, coordinates=coordinates, distances=distances, vlsvReader=vlsvReader))
#         # Get the fitted functions:
#         # Fit a curve into the plot if the user wants to:
#         if fitFunction != nullfit:
#            # Fit the function into raw data:
#            fit = optimize.leastsq(function, np.ones(parameters), args=(x,y))
#            # Get arguments
#            fitargs = []
#            for k in xrange(parameters):
#               fitargs.append(fit[0][k])
#            # Get the X coordinates for the fit plot:
#            X = min(x) + np.arange(100*len(x)).astype(float) / (float)(100*len(x)) * max(x)
#            # Get the Y coordinates for the fit plot:
#            Y = (-1)*function(fitargs, X, 0)
#            
#         else:
#            variables1listfitted = []
#            variables2listfitted = []
#   # Return the fitted and normal variables:
#   return [variables1list, variables2list, variables1listfitted, variables2listfitted]

def draw_plots_by_cellid(vlsvReader, variables1, variables2, cellids, coordinates=[], distances=[], fitFunction=nullfit):
   '''Draws a plot of given variables for the given cell ids
      :param vlsvReader    some VlsvFile with a file open
      :param variables1    some dictionary of variables to be plotted as the x-axis
      :param variables2    some dictionary of variables to be plotted as the y-axis
      :param cellids       some list of cellids
      :param coordinates   some list of coordinates for the cell ids (Only used if variables include the variable \"coordinates\")
      :param distances     some list of distances from the starting cell id (Only used if variables include the variable \"distances\")
      :param fitFunction   some function for making a fit for the x, y plot
   '''
   if fitFunction != nullfit:
      # Get function:
      fetchedFunction = False
      global functionList
      # Try to get the function by name first
      for i in functionList.iteritems():
         if fitFunction == i[0]:
            functionClass = functionList[fitFunction]
            fetchedFunction = True
      # If the user didn't fetch function by name, then it's simply the class itself:
      if fetchedFunction == False:
         functionClass = fitFunction
      function = functionClass.get_function()
      parameters = functionClass.get_parameters()

   # Number of variables to be plotted in y-axis
   xplotnum = 0
   for var in variables1.iteritems():
      xplotnum = xplotnum + len(np.atleast_1d(var[1]))
   # Number of variables to be plotted in x-axis
   yplotnum = 0
   for var in variables2.iteritems():
      yplotnum = yplotnum + len(np.atleast_1d(var[1]))
   # Plot the variables:
   pl.figure()
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
               # Check the vector's index ( 3 stands for the magnitude )
               if i != 3:
                  # Get the vector by index
                  x.append(np.atleast_1d(variableVector)[i])
               else:
                  # Get the vector's magnitude:
                  x.append(np.linalg.norm(np.atleast_1d(variableVector)))
            for j in np.atleast_1d(variable2[1]):
               # Get y plot coordinates:
               y = []
               for variableVector in variable2list:
                  if j != 3:
                     # Get the vector by index
                     y.append(np.atleast_1d(variableVector)[j])
                  else:
                     # Get the vector's magnitude:
                     y.append(np.linalg.norm(np.atleast_1d(variableVector)))
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
               pl.plot(x,y,linestyle='-',marker='.',markersize=2, color='b')
               # Fit a curve into the plot if the user wants to:
               if fitFunction != nullfit:
                  # Fit the function into raw data:
                  fit = optimize.leastsq(function, np.ones(parameters), args=(x,y))
                  # Get arguments
                  fitargs = []
                  for k in xrange(parameters):
                     fitargs.append(fit[0][k])
                  # Get the X coordinates for the fit plot:
                  X = min(x) + np.arange(100*len(x)).astype(float) / (float)(100*len(x)) * max(x)
                  # Get the Y coordinates for the fit plot:
                  Y = (-1)*function(fitargs, X, 0)
                  pl.plot(X,Y,color='r')
               #Increment the index:
               index = index + 1
   # Show the plot:
   pl.tight_layout()
   pl.ion()
   pl.show()

def take_cut_through_array( fileName, point1, point2 ):
   '''Creates an array of cut-through of some variable from a vlsv file
   :param fileName         Name of the vlsv file
   :param point1           The starting coordinates of the line
   :param point2           The ending coordinates of the line
   :returns                Numpy array np.array([cellids, coordinates, distances]) of the cut through
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
   # Return the cut through:
   return [cellids, coordinateList, distances]

def take_cut_through( fileName, variables1, variables2, point1, point2, fitFunction=nullfit):
   '''Creates a plot of cut-through of some variable from a vlsv file
   :param fileName         Name of the vlsv file
   :param variables1       Variables to be plotted in the x-axis
   :param variables2       Variables to be plotted in the y-axis
   :param point1           The starting coordinates of the line
   :param point2           The ending coordinates of the line
   '''
   # Get the cell ids, coordinates and distances from the cut through:
   cutThrough = take_cut_through_array( fileName=fileName, point1=point1, point2=point2 )
   cellids = cutThrough[0]
   coordinateList = cutThrough[1]
   distances = cutThrough[2]
   # Plot the variables:
   draw_plots_by_cellid(vlsvReader=VlsvFile(fileName), variables1=variables1, variables2=variables2, cellids=cellids, coordinates=coordinateList, distances=distances, fitFunction=fitFunction)

def time_evolution_array( fileNames, dt, cellid, variables, forceConstantAmplitude=False, fitFunction=nullfit, kaiserwindowparameter=0 ):
   ''' Plots the time evolution of some cell and fits a fourier series in the plot
       :param fileNames      Name of the files
       :param dt             The time step between frames in the files
       :param cellid         The cell's id
       :param variables      Variables along the y-axis
       :param variables      Some dictionary with variables as key and indexes as values to be plotted along the y-axis
       :param forceConstantAmplitude If true, forces constant amplitudes for any possible waves in the plots (OPTIONAL)
       :param fitFunction    fit function for forcing constant amplitudes, first degree polynomial by default (OPTIONAL)
       :returns              array with fourier series of the variables (and frequencies) [[t2, y2], freq] as a function of t as well as the variables np.array([fourier_variables, y])
       Return graph:
        Array
->       [0]                               |            [1]
   fourier_variables                       |             y
->   [0]     |        [1]     |        [2]          |     ..
fitFunc0     |    fitFunc1    |     fitFunc2        |     ..
->   [0]                  |                    [1]
t2 & y2                   |                    frequencies
-> [0] | [1]
   t2  |  y2
   '''
   # Read and record files
   t = []
   y = []
   for i in variables.iteritems():
      for j in np.atleast_1d(i[1]):
         y.append([])
   index=0
   for f in fileNames:
      # Index of y variable:
      yindex = 0
      # Get time
      t.append(dt*index)
      # Open file for reading
      vlsvReader = VlsvFile(f)
      # Get the variable(s)
      for i in variables.iteritems():
         readVariable = vlsvReader.read_variable(cellid=cellid, name=i[0])
         for j in np.atleast_1d(i[1]):
            y[yindex].append( np.atleast_1d(readVariable)[j] )
            yindex = yindex + 1
      # Increment index
      index = index + 1
   # Check if the user wants to force constant amplitudes
   if forceConstantAmplitude==False:
      fourier_variables = []
      for i in xrange(len(np.atleast_1d(y))):
         fourier_variables.append(fourier_array(t, y[i], kaiserwindowparameter=kaiserwindowparameter))
   else:
      fitFunction = np.atleast_1d(fitFunction)
      fourier_variables=[]
      for j in fitFunction:
         # Fit polynomial of second degree into the data
         # Get fit function data:
         if j==nullfit:
            functionClass = functionList["polynomialfirst"]
         else:
            # Get function:
            fetchedFunction = False
            global functionList
            # Try to get the function by name first
            for i in functionList.iteritems():
               if j == i[0]:
                  functionClass = functionList[j]
                  fetchedFunction = True
            # If the user didn't fetch function by name, then it's simply the class itself:
            if fetchedFunction == False:
               functionClass = j
         # Get the function and the num of parameters
         function = functionClass.get_function()
         parameters = functionClass.get_parameters()
         for i in xrange(len(np.atleast_1d(y))):
            # Fit a polynomial into the data
            fit = optimize.leastsq(function, np.ones(parameters), args=(t,y[i]))
            y_fitted = (-1)*function(fit[0], t, 0)
            # Create a new array y2 which has a forced constant amplitude for the (possible) waves:
            y2 = y[i] - y_fitted
            # Plot the data  with fourier fit
            fourier_variables.append(fourier_array(t, y2, kaiserwindowparameter=kaiserwindowparameter))
   # Return the fourier array:
   return np.array([fourier_variables, y])

def plot_time_evolution( fileNames, dt, cellid, variables, forceConstantAmplitude=False, fitFunction=nullfit, saveplot="none", showplots=True, savedata=False, kaiserwindowparameter=0 ):
   ''' Plots the time evolution of some cell and fits a fourier series in the plot
       :param fileNames      Name of the files
       :param dt             The time step between frames in the files
       :param cellid         The cell's id
       :param variables      Some dictionary with variables as key and indexes as values to be plotted along the y-axis
       :param forceConstantAmplitude If true, forces constant amplitudes for any possible waves in the plots
       :param fitFunction    fit function for forcing constant amplitudes, first degree polynomial by default
   '''
   # Read and record files
   t = []
   #y = [[] for i in xrange(len(np.atleast_1d(variables)))]
   y = []
   names = []
   for i in variables.iteritems():
      for j in np.atleast_1d(i[1]):
         y.append([])
         names.append(i[0] + "[" + str(j) + "]")
   index=0
   for f in fileNames:
      yindex = 0
      # Get time
      t.append(dt*index)
      # Open file for reading
      vlsvReader = VlsvFile(f)
      # Get the variable(s)
      for i in variables.iteritems():
         readVariable = vlsvReader.read_variable(cellid=cellid, name=i[0])
         for j in np.atleast_1d(i[1]):
            y[yindex].append( np.atleast_1d(readVariable)[j] )
            yindex = yindex + 1
      # Increment index
      index = index + 1
   # Plot:
   # Create a new plot window:
   if forceConstantAmplitude==False:
      pl.figure()
      #subplotnums=[[len(y)*2,1,1],[len(y)*2,1,2]]
      plotsPerVariable = (len(np.atleast_1d(kaiserwindowparameter)) + 1)
      subplotnums=[[len(y)*plotsPerVariable, 1, k + 1] for k in range(plotsPerVariable)]
      for i in xrange(len(np.atleast_1d(y))):
         if savedata == True:
            save=saveplot
         else:
            save="none"
         plot_fourier(t, y[i], subplotnums=subplotnums, savedata=save, kaiserwindowparameter=kaiserwindowparameter)
         for k in xrange(len(subplotnums)):
            subplotnums[k][2] = subplotnums[k][2] + plotsPerVariable
         # Save the plot if user wants to:
         pl.tight_layout()
         if saveplot != "none":
            pl.savefig(saveplot + ".ps")
   else:
      fitFunction = np.atleast_1d(fitFunction)
      plotnumber = 0
      for j in fitFunction:
         pl.figure()
         # Fit polynomial of second degree into the data
         # Get fit function data:
         if j==nullfit:
            functionClass = functionList["polynomialfirst"]
         else:
            # Get function:
            fetchedFunction = False
            global functionList
            # Try to get the function by name first
            for i in functionList.iteritems():
               if j == i[0]:
                  functionClass = functionList[j]
                  fetchedFunction = True
            # If the user didn't fetch function by name, then it's simply the class itself:
            if fetchedFunction == False:
               functionClass = j
         function = functionClass.get_function()
         parameters = functionClass.get_parameters()
   
         # Declare subplotnum (For plotting multiple plots in one window)
         plotsPerVariable = (len(np.atleast_1d(kaiserwindowparameter)) + 2)
         numberofVariables = len(y)
         subplotnums=[[numberofVariables*plotsPerVariable, 1, k + 2] for k in range(plotsPerVariable - 1)]
         #subplotnums=[[len(y)*3,1,2],[len(y)*3,1,3]]
         for i in xrange(len(np.atleast_1d(y))):
            # Fit a polynomial into the data
            if functionClass != functionList["emptyfit"]:
               fit = optimize.leastsq(function, np.ones(parameters), args=(t,y[i]))
               y_fitted = (-1)*function(fit[0], t, 0)
            else:
               y_fitted = np.zeros(len(y[i]))
            # Plot t, y:
            pl.subplot(numberofVariables*plotsPerVariable,1,i*plotsPerVariable+1)
            pl.plot(t,y[i], '.')
            # Set title:
            if i == 0:
               pl.title(np.atleast_1d(fileNames)[0] + " - " + np.atleast_1d(fileNames)[len(np.atleast_1d(fileNames)) - 1])
            pl.plot(t,y_fitted)
            from scipy.interpolate import InterpolatedUnivariateSpline
            s = InterpolatedUnivariateSpline(t,y[i])
            pl.plot(t, s(t))
            pl.legend(["data", "fitted_curve", "spline"])
            # Create a new array y2 which has a forced constant amplitude for the (possible) waves:
            y2 = y[i] - y_fitted
            # Plot the data  with fourier fit
            pl.ylabel(names[i])
            if savedata == True:
               save=saveplot
            else:
               save="none"
            plot_fourier(t, y2, subplotnums=subplotnums, savedata=save, kaiserwindowparameter=kaiserwindowparameter)
            for k in xrange(len(subplotnums)):
               subplotnums[k][2] = subplotnums[k][2] + plotsPerVariable
            #subplotnums[0][2] = subplotnums[0][2] + 3
            #subplotnums[1][2] = subplotnums[1][2] + 3
         # Save the plot if user wants to:
         pl.tight_layout()
         if saveplot != "none":
            pl.savefig(saveplot + str(plotnumber) + ".png")
            plotnumber = plotnumber + 1
         if showplots != True:
            pl.close()
   # Show the plots
   if showplots == True:
      pl.ion()
      pl.show()
   else:
      pl.close()

