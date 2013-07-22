import sys
import numpy as np
import pylab as pl
from useroptions import *
from oldvlsvreader import *

class SpaceCraft(object):
   '''
   Class for spacecraft in plasma
   '''
   _cellid=-1
   _coordinates=np.zeros(3)
   _time=0
   _dt=0
   _dt_global=0
   _time_index=0
   _filenames=np.array([""])
   _vlsvfiles=np.array([])
   _values=dict()
   _xmax=1
   _ymax=1
   _zmax=1
   _xmin=1
   _ymin=1
   _zmin=1
   _xcells=1
   _ycells=1
   _zcells=1
   _velocity=np.zeros(3)
   _cell_lengths=np.zeros(3)
   _variables=[]
   def __init__(self, time, dt, dt_global, filenames, variables, coordinates=np.zeros(3), cellid=-1):
      if (coordinates == np.zeros(3)).all() and cellid == -1:
         print "YOU MUST SPECIFY EITHER COORDINATES OR CELLID TO SPACECRAFT INIT"
      if (coordinates != np.zeros(3)).all() and cellid != -1:
         print "YOU MUST SPECIFY EITHER COORDINATES OR CELLID TO SPACECRAFT INIT"
      self._time = time
      self._dt = dt
      self._dt_global = dt_global
      self._time_index = 0
      # Note: Expecting a sorted numpy array of strings
      self._filenames = filenames
      # Create a list of vlsv file objects for the files
      self._vlsvfiles = np.array([VlsvFile(f) for f in filenames])
      # Create a list of values for the variables:
      self._values = dict()
      if len(self._vlsvfiles) == 0:
         print "BAD LENGTH IN SPACECRAFT INIT!"
      # Get xmax, xmin and xcells_ini
      self._xmax = self._vlsvfiles[0].read(name="xmax",tag="PARAMETERS")
      self._xmin = self._vlsvfiles[0].read(name="xmin",tag="PARAMETERS")
      self._xcells = self._vlsvfiles[0].read(name="xcells_ini",tag="PARAMETERS")
      # Do the same for y
      self._ymax = self._vlsvfiles[0].read(name="ymax",tag="PARAMETERS")
      self._ymin = self._vlsvfiles[0].read(name="ymin",tag="PARAMETERS")
      self._ycells = self._vlsvfiles[0].read(name="ycells_ini",tag="PARAMETERS")
      # And for z
      self._zmax = self._vlsvfiles[0].read(name="zmax",tag="PARAMETERS")
      self._zmin = self._vlsvfiles[0].read(name="zmin",tag="PARAMETERS")
      self._zcells = self._vlsvfiles[0].read(name="zcells_ini",tag="PARAMETERS")
      # Get cell lengths
      self._cell_lengths = np.array([(self._xmax - self._xmin)/(float)(self._xcells), (self._ymax - self._ymin)/(float)(self._ycells), (self._zmax - self._zmin)/(float)(self._zcells)])
      # Get the cell id and coordinates:
      if cellid != -1:
         self._cellid = cellid
         self._coordinates = self.get_coordinates(cellid)
      else:
         self._cellid = self.get_cellid(coordinates)
         self._coordinates = coordinates
      # Get the bulk velocity within the cell:
      rho_v = self._vlsvfiles[0].read_variable("rho_v", self._cellid)
      rho = self._vlsvfiles[0].read_variable("rho", self._cellid)
      if len(rho_v) != 3:
         print "BAD RHO_V LENGTH"
      # Check if we're in 2d:
      if self._zcells == 1:
         rho_v[2] = 0
      self._velocity = np.array([rho_v[0]/(float)(rho), rho_v[1]/(float)(rho), rho_v[2]/(float)(rho)])
      # Get variable:
      self._variables = variables
      # Get value arrays for the variables
      self._values = dict()
      for i in xrange(len(self._variables)):
         self._values[self._variables[i]] = []
      
      #self._coordinateslist = []
      # Update variables: (Appends variables at the space craft's current location to list of variables)
      self.record_values()
      # Print for debugging:
      print "xmin:" + str(self._xmin)
      print "ymin:" + str(self._ymin)
      print "zmin:" + str(self._zmin)
      print "xmax:" + str(self._xmax)
      print "ymax:" + str(self._ymax)
      print "zmax:" + str(self._zmax)
      print "xcells:" + str(self._xcells)
      print "ycells:" + str(self._ycells)
      print "zcells:" + str(self._zcells)
      print "cell_lengths:" + str(self._cell_lengths)
      print "time:" + str(self._time)
      print "timeindex:" + str(self._time_index)
      print "dt:" + str(self._dt)
      print "filenames:" + str(self._filenames)
      print "values:" + str(self._values)
      print "variables:" + str(self._variables)
      print "timevariable:" + str(self._timevariable)
      print "velocity:" + str(self._velocity)
      print "coordinates:" + str(self._coordinates)
      print "cellid: " + str(self._cellid)

   def resert_variables(self, time, dt, dt_global, coordinates=np.zeros(3), cellid=-1):
      '''Resets variables (and cellid/coordinates if specified)
      '''
      if (coordinates != np.zeros(3)).all() and cellid != -1:
         print "YOU MUST SPECIFY EITHER COORDINATES OR CELLID TO SPACECRAFT INIT, NOT BOTH"
      # Get the cell id and coordinates:
      if cellid != -1:
         self._cellid = cellid
         self._coordinates = self.get_coordinates(cellid)
      elif (coordinates == np.zeros(3)).all():
         self._cellid = self.get_cellid(coordinates)
         self._coordinates = coordinates
      self._time=time
      self._dt=dt
      self._dt_global=dt_global
      self._time_index=0
      for i in xrange(len(self._values)):
         self._values[i] = []
      self._variables=[]
      self._timevariable=[]
      # Get the bulk velocity within the cell:
      rho_v = self._vlsvfiles[0].read_variable("rho_v", self._cellid)
      rho = self._vlsvfiles[0].read_variable("rho", self._cellid)
      if len(rho_v) != 3:
         print "BAD RHO_V LENGTH"
      # Check if we're in 2d:
      if self._zcells == 1:
         rho_v[2] = 0
      self._velocity = np.array([rho_v[0]/(float)(rho), rho_v[1]/(float)(rho), rho_v[2]/(float)(rho)])
      self._cellidlist = []
      self._coordinateslist = []


   def get_coordinates(self, cellid):
      '''Gets the coordinates of a given cell id
      '''
      # Get the global id first: (Equals cell id minus 1)
      globalid = cellid - 1
      # Get the cell indices:
      cellindices = np.array([globalid%self._xcells, ((int)(globalid)/(int)(self._xcells))%self._ycells, (int)(globalid)/(int)(self._xcells*self._ycells)])
      # Get the cell coordinates:
      cellcoordinates = np.array([self._xmin + cellindices[0]*self._cell_lengths[0], self._ymin + cellindices[1]*self._cell_lengths[1], self._zmin + cellindices[2]*self._cell_lengths[2]])
      # Return the cell coordinates:
      return cellcoordinates

   def get_cellid(self, coordinates):
      '''Gets the cell id of given coordinates
      '''
      # Get cell indices:
      cellindices = np.array([(int)((coordinates[0] - self._xmin)/(float)(self._cell_lengths[0])), (int)((coordinates[1] - self._ymin)/(float)(self._cell_lengths[1])), (int)((coordinates[2] - self._zmin)/(float)(self._cell_lengths[2]))])
      # Get the cell id:
      cellid = cellindices[0] + cellindices[1] * self._xcells + cellindices[2] * self._xcells * self._ycells + 1
      return cellid

   def increase_time(self, updatevariables=False, t_step=-1):
      '''Increases the time variable for space craft and moves the space craft
         :param updatevariables=False      Updates variables for the spacecraft if true
         :param t_step=-1                      The amount of time to step (if -1 then equals the member variable dt)
      '''
      # Get the current time
      currentTime = self._time
      # Get dt
      if t_step == -1:
         dt = self._dt
      else:
         # user set dt:
         dt = t_step
      # Move the spacecraft:
      self._coordinates = self._coordinates + self._velocity * dt
      # Check if the cell has changed:
      newcellid = self.get_cellid(self._coordinates)
      # DEBUGGING:
      print "new cellid:" + str(newcellid)
      print "new coordinates:" + str(self._coordinates)
      # The chances are that by moving dt time the cell won't change, so time variable is the only thing changing
      onlyTimeVariableChanged = True
      if newcellid != self._cellid:
         # Set the new cell id
         self._cellid = newcellid
         # Check if the cell id is within bounds:
         if newcellid >= self._xcells*self._ycells*self._zcells or newcellid <= 0:
            print "CELL ID OUT OF BOUNDS, EXITING increase_time FUNCTION"
            return False
         # Update velocity
         rho_v = self._vlsvfiles[self._time_index].read_variable("rho_v", self._cellid)
         # Check if we're in 2d:
         if self._zcells == 1:
            rho_v[2] = 0
         rho = self._vlsvfiles[self._time_index].read_variable("rho", self._cellid)
         if( rho == 0 ):
            print "BAD RHO IN INCREASE_TIME (SPACE CRAFT)"
            print "rho: " + str(rho)
            print "rho_v: " + str(rho_v)
            print "newcellid: " + str(newcellid)
            return False
         # self._velocity = np.array([(float)(rho_v[self._time_index])/(float)(rho), (float)(rho_v[self._time_index])/(float)(rho), (float)(rho_v[self._time_index])/(float)(rho)])a
         # Update velocity:
         self._velocity = rho_v / (float)(rho)
         # Changed velocity cell so other variables except time has changed, too
         onlyTimeVariableChanged = False
      # Check if the user wants to update variables:
      if updatevariables == True:
         # Update the values:
         self.record_values(updateTimeOnly=onlyTimeVariableChanged)
      # Increase time:
      self._time = self._time + dt
      return True
      # TODO: INCREASE TIMEINDEX HERE
      # if self._time >= self._dt_global*self._time_index:
      #    self._time_index = self._time_index + 1
      #    if self._time_index >= len(self._vlsvfiles):
      #       print "Timelimit reached!"

   def record_next_value(self):
      '''Calculates the time needed to get to the next cell and calls increase_time with the calculated time step
      '''
      # Get the coordinates:
      coordinates = self._coordinates
      # Get the current cell's bounds (The coordinates in x, y, z direction where the cell the space craft is in ends and new cell begins)
      min_bounds = self.get_coordinates(self._cellid)
      max_bounds = np.array([min_bounds[i] + self._cell_lengths[i] for i in range(0,3)])
      all_bounds = np.append(min_bounds, max_bounds)
      # Get the smallest time step required to get to the next cell
      smallest_dt = sys.float_info.max
      for i in xrange(len(min_bounds)):
         # Get the time step required to get to the x, y, z min and max boundary and get the smallest NON-NEGATIVE time step:
         if self._velocity[i] != 0:
            calculated_dt = ((min_bounds[i] - coordinates[i]) - self._cell_lengths[i]*0.0001) / self._velocity[i]
            # Edit smallest dt if calculated_dt is smaller:
            if smallest_dt > calculated_dt and calculated_dt > 0:
               smallest_dt = calculated_dt
      for i in xrange(len(max_bounds)):
         # Get the time step required to get to the x, y, z min and max boundary and get the smallest NON-NEGATIVE time step:
         if self._velocity[i] != 0:
            calculated_dt = ((max_bounds[i] - coordinates[i]) + self._cell_lengths[i]*0.0001) / self._velocity[i]
            # Edit smallest dt if calculated_dt is smaller:
            if smallest_dt > calculated_dt and calculated_dt > 0:
               smallest_dt = calculated_dt
      # Check if we found any good values:
      if smallest_dt == sys.float_info.max:
         print "COULDNT FIND A GOOD DT IN record_next_value"
         print "Coordinates: " + str(coordinates)
         print "Velocity: " + str(self._velocity)
         return False
      # Do the update:
      return self.increase_time(updatevariables=True, t_step=smallest_dt)

   def record_values(self, updateTimeOnly=False):
      '''Records the current variable values of the space craft into a private list (to be plotted later)
      '''
      # Saves variable values!
      # Append the time variable:
      self._timevariable.append(self._time)
      self._cellidlist.append(self._cellid)
      self._coordinateslist.append(self._coordinates)
      # Set other variables if updateTimeOnly is false:
      if updateTimeOnly == False:
         for i in xrange(len(self._variables)):
            # Append the value at spacecraft's position:
            self._values[i].append(self._vlsvfiles[self._time_index].read_variable(self._variables[i], self._cellid))
      else:
         for i in xrange(len(self._variables)):
            # Append the value at spacecraft's position:
            # Note: Appends the latest value again so if list is [[5,2], [0, 3]] then after appending it looks like this: [[5,2,2], [0,3,3]]
            self._values[i].append(self._values[i][len(self._values[i])-1])

   def plot_values(self, variable1, variable2="time"):
      '''Makes a plot of the given variables
      '''
      # Get the index of the variables:
      index1 = -1
      index2 = -1
      for i in xrange(len(self._variables)):
         if variable1 == self._variables[i]:
            index1 = i
         if variable2 == self._variables[i]:
            index2 = i
      # Check if indexes were found:
      #plot_variable1 = np.array([])
      if index1 == -1:
         print "ERROR, COULDNT FIND VARIABLE " + str(variable1)
         return;
      else:
         plot_variable1 = np.array(self._values[index1])
      #plot_variable2 = np.array([])
      if index2 == -1:
         if variable2 != "time":
            print "ERROR, COULDNT FIND VARIABLE " + str(variable2)
            return;
         else:
            plot_variable2 = np.array(self._timevariable)
      else:
         plot_variable2 = np.array(self._values[index2])
      # Plot it:
      print "Lengths:" + str(len(plot_variable1)) + " " + str(len(plot_variable2))
      print "Variables:"
      print plot_variable1
      print plot_variable2
      pl.plot(plot_variable2, plot_variable1)
      pl.xlabel(variable2)
      pl.ylabel(variable1)
      pl.show()


   def list_variables(self):
      print "xmin:" + str(self._xmin)
      print "ymin:" + str(self._ymin)
      print "zmin:" + str(self._zmin)
      print "xmax:" + str(self._xmax)
      print "ymax:" + str(self._ymax)
      print "zmax:" + str(self._zmax)
      print "xcells:" + str(self._xcells)
      print "ycells:" + str(self._ycells)
      print "zcells:" + str(self._zcells)
      print "time:" + str(self._time)
      print "timeindex:" + str(self._time_index)
      print "dt:" + str(self._dt)
      print "filenames:" + str(self._filenames)
      print "values:" + str(self._values)

