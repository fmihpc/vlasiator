import sys
import numpy as np
import pylab as pl
from useroptions import *
from vlsvreader import *

def make_moving_spacecraft_plot_scalars( filename, variables_to_plot, x_variable, starting_coordinates, time_limit ):
   '''
   Creates a moving spacecraft that plots the given variables as a function of time
   :param filename              Name of the vlsv or silo file
   :param variables_to_plot     Array of variables to be plotted
   :param x_variable            Name of the variable to be plotted against the x-axis (Note: other than the variables to be plotted x_variable can be \"time\", \"xcoordinates\", \"ycoordinates\", \"zcoordinates\" and \"cellid\"
   :param starting_coordinates  Starting coordinates for the spacecraft
   :param time_limit            The amount of time in seconds that the spacecraft will move (or until hits a boundary or gets stuck)
   '''
   if time_limit <= 0:
      print "TIME_LIMIT MUST BE GREATER THAN ZERO!"
      return
   mySpaceCraft = SpaceCraft(time=0, dt=0.5, dt_global=0.8, filenames=np.array([filename]), variables=np.array(variables_to_plot), coordinates=np.array(starting_coordinates), cellid=-1)
   while True:
      if mySpaceCraft.record_next_value() == False:
         print "HIT BOUNDARY OR GOT STUCK"
         break
      if mySpaceCraft.get_time() >= time_limit:
         print "TIME LIMIT HIT"
         break
   i = 0
   for varName in variables_to_plot:
      mySpaceCraft.plot_values(variable1=varName, variable2=x_variable, variable1index=-1, variable2index=-1)
   return

def make_moving_spacecraft_plot( filename, variables_to_plot, x_variable, starting_coordinates, time_limit ):
   '''
   Creates a moving spacecraft that plots the given variables as a function of time
   :param filename              Name of the vlsv or silo file
   :param variables_to_plot     Dictionary of variables to be plotted, the key must be name of the variable and the value should be the index
   :param x_variable            Name of the variable to be plotted against the x-axis (Note: other than the variables to be plotted x_variable can be \"time\", \"xcoordinates\", \"ycoordinates\", \"zcoordinates\" and \"cellid\"
   :param starting_coordinates  Starting coordinates for the spacecraft
   :param time_limit            The amount of time in seconds that the spacecraft will move (or until hits a boundary or gets stuck)
   '''
   if len(x_variable) != 1:
      print "Bad x_variable size!"
      return
   if time_limit <= 0:
      print "TIME_LIMIT MUST BE GREATER THAN ZERO!"
      return
   # Make sure x_variable is a dictionary
   x_var = dict(x_variable)
   # Get variable names:
   variableNames = []
   for i in x_var.iteritems():
      x_variableName = i[0]
      # Check length
      if len([i[1]]) != 1:
         print "Bad x_variable index length!"
         return
   for i in variables_to_plot.iteritems():
      variableNames.append(i[0])
   # Create spacecraft
   mySpaceCraft = SpaceCraft(time=0, dt=0.5, dt_global=0.8, filenames=np.array([filename]), variables=np.array(variableNames), coordinates=np.array(starting_coordinates), cellid=-1)
   # Move spacecraft
   while True:
      if mySpaceCraft.record_next_value() == False:
         print "HIT BOUNDARY OR GOT STUCK"
         break
      if mySpaceCraft.get_time() >= time_limit:
         print "TIME LIMIT HIT"
         break
   # Get number of plots:
   num_of_plots = 0
   for i in variables_to_plot.iteritems():
      indexes = np.array(i[1])
      for j in xrange(len(indexes)):
         num_of_plots = num_of_plots + 1
   # Plot results
   plotindice = 1
   for i in variables_to_plot.iteritems():
      varName = i[0]
      indexes = np.array(i[1])
      for j in xrange(len(indexes)):
         # Plot all indexes
         mySpaceCraft.plot_values(variable1=varName, variable2=x_variableName, variable1index=j, variable2index=-1, showplot=False, subplotnum=[num_of_plots,1,plotindice])
         plotindice = plotindice + 1
   # Show the plots
   pl.ion()
   pl.show()
   return


class SpaceCraft(object):
   '''
   Class for spacecraft in plasma
   '''
   __cellid=-1
   __coordinates=np.zeros(3)
   __time=0
   __dt=0
   __dt_global=0
   __time_index=0
   __filenames=np.array([""])
   __vlsvfiles=np.array([])
   __values=dict()
   __xmax=1
   __ymax=1
   __zmax=1
   __xmin=1
   __ymin=1
   __zmin=1
   __xcells=1
   __ycells=1
   __zcells=1
   __velocity=np.zeros(3)
   __cell_lengths=np.zeros(3)
   __variables=[]
   def __init__(self, time, dt, dt_global, filenames, variables, coordinates=np.zeros(3), cellid=-1):
      if (coordinates == np.zeros(3)).all() and cellid == -1:
         print "YOU MUST SPECIFY EITHER COORDINATES OR CELLID TO SPACECRAFT INIT"
      if (coordinates != np.zeros(3)).all() and cellid != -1:
         print "YOU MUST SPECIFY EITHER COORDINATES OR CELLID TO SPACECRAFT INIT"
      self.__time = time
      self.__dt = dt
      self.__dt_global = dt_global
      self.__time_index = 0
      # Note: Expecting a sorted numpy array of strings
      self.__filenames = filenames
      # Create a list of vlsv file objects for the files
      self.__vlsvfiles = np.array([VlsvFile(f) for f in filenames])
      if len(self.__vlsvfiles) == 0:
         print "BAD LENGTH IN SPACECRAFT INIT!"
      # Get xmax, xmin and xcells_ini
      self.__xmax = self.__vlsvfiles[0].read_parameter(name="xmax")
      self.__xmin = self.__vlsvfiles[0].read_parameter(name="xmin")
      self.__xcells = self.__vlsvfiles[0].read_parameter(name="xcells_ini")
      # Do the same for y
      self.__ymax = self.__vlsvfiles[0].read_parameter(name="ymax")
      self.__ymin = self.__vlsvfiles[0].read_parameter(name="ymin")
      self.__ycells = self.__vlsvfiles[0].read_parameter(name="ycells_ini")
      # And for z
      self.__zmax = self.__vlsvfiles[0].read_parameter(name="zmax")
      self.__zmin = self.__vlsvfiles[0].read_parameter(name="zmin")
      self.__zcells = self.__vlsvfiles[0].read_parameter(name="zcells_ini")
      
      # Get cell lengths
      self.__cell_lengths = np.array([(self.__xmax - self.__xmin)/(float)(self.__xcells), (self.__ymax - self.__ymin)/(float)(self.__ycells), (self.__zmax - self.__zmin)/(float)(self.__zcells)])
      # Get the cell id and coordinates:
      if cellid != -1:
         self.__cellid = cellid
         self.__coordinates = self.__get_coordinates(cellid)
      else:
         self.__cellid = self.__get_cellid(coordinates)
         self.__coordinates = coordinates
      # Get the bulk velocity within the cell:
      rho_v = self.__vlsvfiles[0].read_variable("rho_v", self.__cellid)
      rho = self.__vlsvfiles[0].read_variable("rho", self.__cellid)
      if len(rho_v) != 3:
         print "BAD RHO_V LENGTH"
      # Check if we're in 2d:
      if self.__zcells == 1:
         rho_v[2] = 0
      self.__velocity = np.array([rho_v[0]/(float)(rho), rho_v[1]/(float)(rho), rho_v[2]/(float)(rho)])
      # Get variable:
      self.__variables = variables
      # Get value arrays for the variables
      self.__values = dict()
      for varName in self.__variables:
         self.__values[varName] = []
      # Set other value arrays:
      self.__values["time"] = []
      self.__values["cellid"] = []
      self.__values["xcoordinates"] = []
      self.__values["ycoordinates"] = []
      self.__values["zcoordinates"] = []
      #self.__coordinateslist = []
      # Update variables: (Appends variables at the space craft's current location to list of variables)
      self.record_values()
      # Print for debugging:
      print "xmin:" + str(self.__xmin)
      print "ymin:" + str(self.__ymin)
      print "zmin:" + str(self.__zmin)
      print "xmax:" + str(self.__xmax)
      print "ymax:" + str(self.__ymax)
      print "zmax:" + str(self.__zmax)
      print "xcells:" + str(self.__xcells)
      print "ycells:" + str(self.__ycells)
      print "zcells:" + str(self.__zcells)
      print "cell_lengths:" + str(self.__cell_lengths)
      print "time:" + str(self.__time)
      print "timeindex:" + str(self.__time_index)
      print "dt:" + str(self.__dt)
      print "filenames:" + str(self.__filenames)
      print "values:" + str(self.__values)
      print "variables:" + str(self.__variables)
      print "velocity:" + str(self.__velocity)
      print "coordinates:" + str(self.__coordinates)
      print "cellid: " + str(self.__cellid)
      print "values: "
      print self.__values

   def __get_coordinates(self, cellid):
      '''Gets the coordinates of a given cell id
      '''
      # Get the global id first: (Equals cell id minus 1)
      globalid = cellid - 1
      # Get the cell indices:
      cellindices = np.array([globalid%self.__xcells, ((int)(globalid)/(int)(self.__xcells))%self.__ycells, (int)(globalid)/(int)(self.__xcells*self.__ycells)])
      # Get the cell coordinates:
      cellcoordinates = np.array([self.__xmin + cellindices[0]*self.__cell_lengths[0], self.__ymin + cellindices[1]*self.__cell_lengths[1], self.__zmin + cellindices[2]*self.__cell_lengths[2]])
      # Return the cell coordinates:
      return cellcoordinates

   def __get_cellid(self, coordinates):
      '''Gets the cell id of given coordinates
      '''
      # Get cell indices:
      cellindices = np.array([(int)((coordinates[0] - self.__xmin)/(float)(self.__cell_lengths[0])), (int)((coordinates[1] - self.__ymin)/(float)(self.__cell_lengths[1])), (int)((coordinates[2] - self.__zmin)/(float)(self.__cell_lengths[2]))])
      # Get the cell id:
      cellid = cellindices[0] + cellindices[1] * self.__xcells + cellindices[2] * self.__xcells * self.__ycells + 1
      return cellid

   def get_time(self):
      return self.__time

   def reset_variables(self, time=-1, dt=-1, dt_global=-1, coordinates=np.zeros(3), cellid=-1):
      '''Resets variables (and cellid/coordinates if specified)
      '''
      if (coordinates != np.zeros(3)).all() and cellid != -1:
         print "YOU MUST SPECIFY EITHER COORDINATES OR CELLID TO SPACECRAFT INIT, NOT BOTH"
      # Get the cell id and coordinates:
      if cellid != -1:
         self.__cellid = cellid
         self.__coordinates = self.__get_coordinates(cellid)
      elif (coordinates != np.zeros(3)).any():
         self.__cellid = self.__get_cellid(coordinates)
         self.__coordinates = coordinates
      if time != -1:
         self.__time=time
      if dt != -1:
         self.__dt=dt
      if dt_global != -1:
         self.__dt_global=dt_global
      self.__time_index=0
      for varName in self.__variables:
         self.__values[varName] = []
      # Set other value arrays:
      self.__values["time"] = []
      self.__values["cellid"] = []
      self.__values["xcoordinates"] = []
      self.__values["ycoordinates"] = []
      self.__values["zcoordinates"] = []
      # Get the bulk velocity within the cell:
      rho_v = self.__vlsvfiles[0].read_variable("rho_v", self.__cellid)
      rho = self.__vlsvfiles[0].read_variable("rho", self.__cellid)
      if len(rho_v) != 3:
         print "BAD RHO_V LENGTH"
      # Check if we're in 2d:
      if self.__zcells == 1:
         rho_v[2] = 0
      if rho == 0:
         print "BAD RHO: " + str(rho)
         print "CELLID: " + str(self.__cellid)
         rho = 0.00000001
      self.__velocity = np.array([rho_v[0]/(float)(rho), rho_v[1]/(float)(rho), rho_v[2]/(float)(rho)])
      self.__cellidlist = []
      self.__coordinateslist = []
      # Update variables: (Appends variables at the space craft's current location to list of variables)
      self.record_values()


   def increase_time(self, updatevariables=False, t_step=-1):
      '''Increases the time variable for space craft and moves the space craft
         :param updatevariables=False      Updates variables for the spacecraft if true
         :param t_step=-1                      The amount of time to step (if -1 then equals the member variable dt)
      '''
      # Get the current time
      currentTime = self.__time
      # Get dt
      if t_step == -1:
         dt = self.__dt
      else:
         # user set dt:
         dt = t_step
      # Move the spacecraft:
      self.__coordinates = self.__coordinates + self.__velocity * dt
      # Check if the cell has changed:
      newcellid = self.__get_cellid(self.__coordinates)
      # The chances are that by moving dt time the cell won't change, so time variable is the only thing changing
      onlyTimeVariableChanged = True
      if newcellid != self.__cellid:
         # Set the new cell id
         self.__cellid = newcellid
         # Check if the cell id is within bounds:
         if newcellid >= self.__xcells*self.__ycells*self.__zcells or newcellid <= 0:
            print "CELL ID OUT OF BOUNDS, EXITING increase_time FUNCTION"
            return False
         # Update velocity
         rho_v = self.__vlsvfiles[self.__time_index].read_variable("rho_v", self.__cellid)
         # Check if we're in 2d:
         if self.__zcells == 1:
            rho_v[2] = 0
         rho = self.__vlsvfiles[self.__time_index].read_variable("rho", self.__cellid)
         if( rho == 0 ):
            print "BAD RHO IN INCREASE_TIME (SPACE CRAFT)"
            print "rho: " + str(rho)
            print "rho_v: " + str(rho_v)
            print "newcellid: " + str(newcellid)
            return False
         # self.__velocity = np.array([(float)(rho_v[self.__time_index])/(float)(rho), (float)(rho_v[self.__time_index])/(float)(rho), (float)(rho_v[self.__time_index])/(float)(rho)])a
         # Update velocity:
         self.__velocity = rho_v / (float)(rho)
         # Changed velocity cell so other variables except time has changed, too
         onlyTimeVariableChanged = False
      # Check if the user wants to update variables:
      if updatevariables == True:
         # Update the values:
         self.record_values(updateTimeOnly=onlyTimeVariableChanged)
      # Increase time:
      self.__time = self.__time + dt
      return True
      # TODO: INCREASE TIMEINDEX HERE
      # if self.__time >= self.__dt_global*self.__time_index:
      #    self.__time_index = self.__time_index + 1
      #    if self.__time_index >= len(self.__vlsvfiles):
      #       print "Timelimit reached!"

   def record_next_value(self):
      '''Calculates the time needed to get to the next cell and calls increase_time with the calculated time step
      '''
      # Get the coordinates:
      coordinates = self.__coordinates
      # Get the current cell's bounds (The coordinates in x, y, z direction where the cell the space craft is in ends and new cell begins)
      min_bounds = self.__get_coordinates(self.__cellid)
      max_bounds = np.array([min_bounds[i] + self.__cell_lengths[i] for i in range(0,3)])
      # Get the smallest time step required to get to the next cell
      smallest_dt = sys.float_info.max
      for i in xrange(len(min_bounds)):
         # Get the time step required to get to the x, y, z min and max boundary and get the smallest NON-NEGATIVE time step:
         if self.__velocity[i] != 0:
            calculated_dt = ((min_bounds[i] - coordinates[i]) - self.__cell_lengths[i]*0.0001) / self.__velocity[i]
            # Edit smallest dt if calculated_dt is smaller:
            if smallest_dt > calculated_dt and calculated_dt > 0:
               smallest_dt = calculated_dt
      for i in xrange(len(max_bounds)):
         # Get the time step required to get to the x, y, z min and max boundary and get the smallest NON-NEGATIVE time step:
         if self.__velocity[i] != 0:
            calculated_dt = ((max_bounds[i] - coordinates[i]) + self.__cell_lengths[i]*0.0001) / self.__velocity[i]
            # Edit smallest dt if calculated_dt is smaller:
            if smallest_dt > calculated_dt and calculated_dt > 0:
               smallest_dt = calculated_dt
      # Check if we found any good values:
      if smallest_dt == sys.float_info.max:
         print "COULDNT FIND A GOOD DT IN record_next_value"
         print "Coordinates: " + str(coordinates)
         print "Velocity: " + str(self.__velocity)
         return False
      # Do the update:
      return self.increase_time(updatevariables=True, t_step=smallest_dt)

   def record_values(self, updateTimeOnly=False):
      '''Records the current variable values of the space craft into a private list (to be plotted later)
      '''
      # Saves variable values!
      # Append the time variable:
      self.__values["time"].append([self.__time])
      self.__values["cellid"].append([self.__cellid])
      self.__values["xcoordinates"].append([self.__coordinates[0]])
      self.__values["ycoordinates"].append([self.__coordinates[1]])
      self.__values["zcoordinates"].append([self.__coordinates[2]])
      #self.__timevariable.append(self.__time)
      #self.__cellidlist.append(self.__cellid)
      #self.__coordinateslist.append(self.__coordinates)
      # Set other variables if updateTimeOnly is false:
      if updateTimeOnly == False:
         for varName in self.__variables:
            # Append the value at spacecraft's position:
            self.__values[varName].append(np.array(self.__vlsvfiles[self.__time_index].read_variable(varName, self.__cellid)))
      else:
         for i in self.__variables:
            # Append the value at spacecraft's position:
            # Note: Appends the latest value again so if list is [[5,2], [0, 3]] then after appending it looks like this: [[5,2,2], [0,3,3]]
            self.__values[varName].append(np.array(self.__values[varName][len(self.__values[varName])-1]))

   def plot_values(self, variable1, variable2="time", variable1index=-1, variable2index=-1, showplot=True, subplotnum=[1,1,1], saveplot=False, fromIndex=-1, toIndex=-1):
      '''Makes a plot of the given variables
      '''
      # Get the index of the variables:
      if (variable1 in self.__values) == False:
         print "Variable " + variable1 + " not found!"
         return
      if (variable2 in self.__values) == False:
         print "Variable " + variable2 + " not found!"
         return
      # Input variables:
      if fromIndex == -1 and toIndex == -1:
         var1 = np.array(self.__values[variable1])
         var2 = np.array(self.__values[variable2])
      else:
         # User input
         if toIndex >= len(self.__values[variable1]):
            print "Bad index (index out of bounds)!"
            return
         var1 = np.array([self.__values[variable1][i] for i in xrange(fromIndex, toIndex+1)])
         var2 = np.array([self.__values[variable2][i] for i in xrange(fromIndex, toIndex+1)])
      if len(var1) == 0 or len(var2) == 0:
         print "Bad variable length"
         return
      
      # Input variable 1
      if len(np.atleast_1d(var1[0])) != 1:
         # The variable is a vector
         if variable1index == -1:
            print "Variable is of size " + str(len(var1[0])) + " so can't input variable1index -1"
            return
         if variable1index >= len(var1[0]):
            print "Bad variable1index! Index out of bounds!"
            return
         plot_variable1 = []
         # Input vector's values for index variable1index
         for i in xrange(len(var1)):
            plot_variable1.append(var1[i][variable1index])
         # Transform to numpy array
         plot_variable1 = np.array(plot_variable1)
      else:
         print "This called"
         print var1
         # The variable is a scalar
         plot_variable1 = np.array(var1)
      
      # Input variable 2
      if len(np.atleast_1d(var2[0])) != 1:
         # The variable is a vector
         if variable2index == 0:
            print "Variable is of size " + str(len(var1[0])) + " so can't input variable1index -1"
            return
         if variable2index >= len(var2[0]):
            print "Bad variable1index! Index out of bounds!"
            return
         plot_variable2 = []
         # Input vector's values for index variable1index
         for i in xrange(len(var2)):
            plot_variable2.append(var2[i][variable2index])
         # Transform to numpy array
         plot_variable2 = np.array(plot_variable2)
      else:
         # The variable is a scalar
         plot_variable2 = np.array(var2)
      # Plot it:
      if showplot == False:
         if len(subplotnum) != 3:
            print "Bad subplotnum length"
            return
         pl.subplot(subplotnum[0], subplotnum[1], subplotnum[2])
      pl.ion()
      pl.plot(plot_variable2, plot_variable1)
      pl.xlabel(variable2 + "[" + str(variable2index) + "]")
      pl.ylabel(variable1 + "[" + str(variable1index) + "]")
      # Set x and y limits:
      xlength = max(plot_variable2) - min(plot_variable2)
      xmin = min(plot_variable2) - 0.01*xlength
      xmax = max(plot_variable2) + 0.01*xlength
      ylength = max(plot_variable1) - min(plot_variable1)
      ymin = min(plot_variable1) - 0.01*ylength
      ymax = max(plot_variable1) + 0.01*ylength
      
      pl.xlim((xmin, xmax))
      pl.ylim((ymin, ymax))
      if showplot == True:
         pl.show()
      return

   def plot_route(self):
      pl.ion()
      pl.plot(self.__values["xcoordinates"], self.__values["ycoordinates"], '.')
      xmin = self.__xmin
      xmax = self.__xmax
      ymin = self.__ymin
      ymax = self.__ymax
      pl.xlim((xmin, ymax))
      pl.ylim((ymin, ymax))
      pl.show()

   def list_variables(self):
      print "xmin:" + str(self.__xmin)
      print "ymin:" + str(self.__ymin)
      print "zmin:" + str(self.__zmin)
      print "xmax:" + str(self.__xmax)
      print "ymax:" + str(self.__ymax)
      print "zmax:" + str(self.__zmax)
      print "xcells:" + str(self.__xcells)
      print "ycells:" + str(self.__ycells)
      print "zcells:" + str(self.__zcells)
      print "time:" + str(self.__time)
      print "timeindex:" + str(self.__time_index)
      print "dt:" + str(self.__dt)
      print "filenames:" + str(self.__filenames)
      print "values:" + str(self.__values)

