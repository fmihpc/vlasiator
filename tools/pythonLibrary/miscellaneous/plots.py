import pylab as pl
import numpy as np
from plotfunctions import *
from vlsvreader import *
from vlsvplots import *
from scipy.interpolate import interp1d

def get_function( fitfunction ):
   # Get function:
   fetchedFunction = False
   global functionList
   # Try to get the function by name first
   for i in functionList.iteritems():
      if fitfunction == i[0]:
         functionClass = functionList[fitfunction]
         fetchedFunction = True
   # If the user didn't fetch function by name, then it's simply the class itself:
   if fetchedFunction == False:
      functionClass = fitfunction
   function = functionClass.get_function()
   return function

def get_parameters( fitfunction ):
   # Get function:
   fetchedFunction = False
   global functionList
   # Try to get the function by name first
   for i in functionList.iteritems():
      if fitfunction == i[0]:
         functionClass = functionList[fitfunction]
         fetchedFunction = True
   # If the user didn't fetch function by name, then it's simply the class itself:
   if fetchedFunction == False:
      functionClass = fitfunction
   parameters = functionClass.get_parameters()
   return parameters

def fit_function_into_data( x, y, fitfunction ):
   function = get_function(fitfunction)
   parameters = get_parameters(fitfunction)
   fit = optimize.leastsq(function, np.ones(parameters), args=(x,y))
   print "Parameters: " + str(fit[0])
   y_fitted = (-1)*function(fit[0], x, 0)
   return y_fitted

def plot_variables( x, y, showplots=False, fitfunction="nullfit", cubicspline=True ):
   pl.plot(x,y)
   if cubicspline == True:
      f = interp1d(x, y, kind='cubic')
      xlength = max(x) - min(x)
      steps = np.arange(len(x)*3) / (float)(len(x)*3)
      xnew = min(x) + steps*xlength
      pl.plot(xnew, f(xnew))
   # Fit a function:
   if fitfunction != "nullfit":
      y_fitted = fit_function_into_data(x, y, fitfunction)
      pl.plot(x,y_fitted)
   if showplots == True:
      pl.ion()
      pl.show()

#def plot_cut_through_var_vs_distance( fileName, variables, point1, point2 ):
#   variables = dict()
#   cutThrough = take_cut_through_array( fileName=fileName, point1=point1, point2=point2)
#   cellids = cutThrough[0]
#   coordinates = cutThrough[1]
#   distances = cutThrough[2]
#   vlsvReader = VlsvFile("bulk.0000519.vlsv")
#   for i in variables:
#      variables[i] = vlsvReader.read_variables_for_cellids("rho", cellids)
#      plot_variables( x = distances, y = variables["rho"], showplots = False, fitfunction="nullfit" )

def take_cut_through_variables( fileName, variableNames, point1, point2 ):
   variables = dict()
   cutThrough = take_cut_through_array( fileName=fileName, point1=point1, point2=point2)
   cellids = cutThrough[0]
   coordinates = cutThrough[1]
   distances = cutThrough[2]
   vlsvReader = VlsvFile(fileName)
   for i in variableNames:
      variables[i] = vlsvReader.read_variables_for_cellids(i, cellids)
   variables["cellids"] = cellids
   variables["coordinates"] = coordinates
   variables["distances"] = distances
   return variables

def plot_from_list( list_of_coordinates, fitfunction="nullfit", cubicspline=True, xticks=-1, yticks=-1 ):
   '''Plots the given list. The list should be in the following format:
      [[x1,y1], [x2,y2], [x3,y3], [x4,y4], .. ]
      :param list_of_coordinates    List of coordinates to be plotted
   '''
   index = 1
   fig = pl.figure()
   for i in list_of_coordinates:
      #title="(" + i[2] + "," + i[3] + ")"
      pl.subplot(len(list_of_coordinates), 1, index)
      #ax = fig.add_subplot(len(list_of_coordinates), 1, index)
      x = i[0]
      y = i[1]
      pl.xlabel(i[2])
      pl.ylabel(i[3])
      xlength = max(x) - min(x)
      ylength = max(y) - min(y)
      pl.xlim([min(x) - 0.01*xlength, max(x) + 0.01*xlength])
      pl.ylim([min(y) - 0.05*ylength, max(y) + 0.05*ylength])
      plot_variables(x,y,False,fitfunction,cubicspline=cubicspline)
      pl.ticklabel_format(style='sci', axis='y', scilimits=(0,3))
      # Set y ticks
      if yticks > 0 and len(i) == 4:
         ticks = min(y) + np.arange(yticks + 1) / (float)(yticks)*ylength
         print len(i)
         from decimal import *
         getcontext().prec = 2
         # Get two decimals
         ticks = [float(Decimal(j)/Decimal(1)) for j in ticks]
         pl.yticks(ticks)
         print ticks
      # Set x ticks
      if xticks > 0 and len(i) == 4:
         ticks = min(x) + np.arange(xticks + 1) / (float)(xticks)*xlength
         from decimal import *
         getcontext().prec = 2
         # Get two decimals
         ticks = [float(Decimal(j)/Decimal(1)) for j in ticks]
         pl.xticks(ticks)
      index = index + 1
   #pl.tight_layout()
   pl.ion()
   pl.show()

















