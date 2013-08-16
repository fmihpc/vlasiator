from vlsvreader import *
import numpy as np
import pylab as pl

def get_pitch_angles( fileName, cellid ):
   # Open the file
   vlsvReader = VlsvFile( fileName )
   # Read the velocity cells:
   velocity_cell_data = vlsvReader.read_velocity_cells(cellid)
   # Calculate the pitch angles for the data:
   B = vlsvReader.read_variable("B", cellid)
   B_unit = B / np.linalg.norm(B)
   pitch_angles = []
   avgs = []
   for i in velocity_cell_data.iteritems():
      # Get the velocity:
      vcellid = i[0]
      avg = i[1]
      v = vlsvReader.get_velocity_cell_coordinates(vcellid)
      # Get the angle
      angle = np.arccos(v.dot(B_unit) / np.linalg.norm(v))
      # Input pitch angle and avgs:
      pitch_angles.append(angle)
      avgs.append(avg)
   # Put the angles into plottable array:
   pitch_angles_array = []
   avgs_array = []
   for i in pitch_angles.iteritems():
      

#      # Get the parallel component:
#      v_par = v * B_unit
#      # Get the perpendicular component:
#      v_per = v - v_par
#      angle = np.arctan(np.linalg.norm(v_per) / np.linalg.norm(v_par))

