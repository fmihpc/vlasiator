from vlsvreader import *
import numpy as np

def get_cell_id( bounds, coordinates ):
   '''Returns the cell id of of given coordinates
      :param bounds         bounds of the simulation (=[xmin,xmax,xcells,ymin,ymax,ycells,zmin,zmax,zcells])
      :param coordinates    coordinates of the cell
      :returns cell id
   '''
   # Get xmin, xmax, etc
   xmin = bounds[0]
   ymin = bounds[0 + 3]
   zmin = bounds[0 + 3 + 3]

   xmax = bounds[1]
   ymax = bounds[1 + 3]
   zmax = bounds[1 + 3 + 3]

   xcells = (int)(bounds[2])
   ycells = (int)(bounds[2 + 3])
   zcells = (int)(bounds[2 + 3 + 3])

   # Get cell lengths:
   cell_lengths = np.array([(xmax - xmin)/(float)(xcells), (ymax - ymin)/(float)(ycells), (zmax - zmin)/(float)(zcells)])
   # Get cell indices:
   cellindices = np.array([(int)((coordinates[0] - xmin)/(float)(cell_lengths[0])), (int)((coordinates[1] - ymin)/(float)(cell_lengths[1])), (int)((coordinates[2] - zmin)/(float)(cell_lengths[2]))])
   # Get the cell id:
   cellid = cellindices[0] + cellindices[1] * xcells + cellindices[2] * xcells * ycells + 1
   return cellid

def get_cell_id_from_file( fileName, coordinates ):
   '''Returns the cell id of of given coordinates
      :param fileName       Name of the file where VlsvFile reads the bounds of the simulation
      :param coordinates    coordinates of the cell
      :returns cell id
      Example:
      get_cell_id("bulk.0000042.vlsv", [42, 4e6, 0])
   '''
   # Create a list of vlsv file objects for the files
   vlsvReader = VlsvFile(fileName)
   # Get xmax, xmin and xcells_ini
   xmax = vlsvReader.read_parameter(name="xmax")
   xmin = vlsvReader.read_parameter(name="xmin")
   xcells = (int)(vlsvReader.read_parameter(name="xcells_ini"))
   # Do the same for y
   ymax = vlsvReader.read_parameter(name="ymax")
   ymin = vlsvReader.read_parameter(name="ymin")
   ycells = (int)(vlsvReader.read_parameter(name="ycells_ini"))
   # And for z
   zmax = vlsvReader.read_parameter(name="zmax")
   zmin = vlsvReader.read_parameter(name="zmin")
   zcells = (int)(vlsvReader.read_parameter(name="zcells_ini"))

   # Get cell lengths:
   cell_lengths = np.array([(xmax - xmin)/(float)(xcells), (ymax - ymin)/(float)(ycells), (zmax - zmin)/(float)(zcells)])
   # Get cell indices:
   cellindices = np.array([(int)((coordinates[0] - xmin)/(float)(cell_lengths[0])), (int)((coordinates[1] - ymin)/(float)(cell_lengths[1])), (int)((coordinates[2] - zmin)/(float)(cell_lengths[2]))])
   # Get the cell id:
   cellid = cellindices[0] + cellindices[1] * xcells + cellindices[2] * xcells * ycells + 1
   return cellid


def get_cell_coordinates_from_file( fileName, cellid ):
   '''Returns the cell id of of given coordinates
      :param fileName       Name of the file where VlsvFile reads the bounds of the simulation
      :param cellid         ID of the cell
      :returns coordinates
      Example:
      coordinates = get_cell_id("bulk.0000042.vlsv", 14923)
   '''
   # Create a list of vlsv file objects for the files
   vlsvReader = VlsvFile(fileName)
   # Get xmax, xmin and xcells_ini
   xmax = vlsvReader.read_parameter(name="xmax")
   xmin = vlsvReader.read_parameter(name="xmin")
   xcells = (int)(vlsvReader.read_parameter(name="xcells_ini"))
   # Do the same for y
   ymax = vlsvReader.read_parameter(name="ymax")
   ymin = vlsvReader.read_parameter(name="ymin")
   ycells = (int)(vlsvReader.read_parameter(name="ycells_ini"))
   # And for z
   zmax = vlsvReader.read_parameter(name="zmax")
   zmin = vlsvReader.read_parameter(name="zmin")
   zcells = (int)(vlsvReader.read_parameter(name="zcells_ini"))

   # Get cell lengths:
   cell_lengths = np.array([(xmax - xmin)/(float)(xcells), (ymax - ymin)/(float)(ycells), (zmax - zmin)/(float)(zcells)])

   cellindices = np.zeros(3)

   # Get cell indices:
   cellid = (int)(cellid - 1)
   cellindices[0] = (int)(cellid)%(int)(xcells)
   #cellindices[1] = ((int)(cellid - cellindices[0])%(int)(ycells))/xcells
   #cellindices[2] = (cellid - cellindices[1]*xcells - cellindices[0])/(xcells*ycells)
   #cellindices[1] = ((int)(cellid)%(int)(ycells))/(int)(xcells)
   #cellindices[2] = (int)(cellid)/(int)(xcells*ycells)
   cellindices[1] = ((int)(cellid)/(int)(xcells))%(int)(ycells)
   cellindices[2] = (int)(cellid)/(int)(xcells*ycells)

   # Get cell coordinates:
   cellcoordinates = np.zeros(3)
   cellcoordinates[0] = xmin + (cellindices[0] + 0.5) * cell_lengths[0]
   cellcoordinates[1] = ymin + (cellindices[1] + 0.5) * cell_lengths[1]
   cellcoordinates[2] = zmin + (cellindices[2] + 0.5) * cell_lengths[2]
   # Return the coordinates:
   return np.array(cellcoordinates)




