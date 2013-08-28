import vlsvreader
from numpy import mgrid, empty, sin, pi, ravel
from tvtk.api import tvtk
from mayavi import mlab
#from mayavi.api import engine
import numpy as np


class MayaviPlots:
   '''Class for constructing plots with MayaVi
   '''
   def __init__(self, vlsvReader):
      print "Constructing mayavi plot"
      self.__vlsvReader = vlsvReader

   def __picker_callback( self, picker ):
      """ Picker callback: this get called when on pick events.
      """
      # Find which data point corresponds to the point picked:
      # we have to account for the fact that each data point is
      # represented by a glyph with several points
      #point_id = picker.point_id/glyph_points.shape[0]
      point_id = picker.cell_id
      # If the no points have been selected, we have '-1'
      print point_id
#      if point_id != -1:
#         # Retrieve the coordinnates coorresponding to that data
#         # point
#         x, y, z = x1[point_id], y1[point_id], z1[point_id]
#         print x
#         print y
#         print z
#         # Move the outline to the data point.
#         outline.bounds = (x-0.1, x+0.1,
#                       y-0.1, y+0.1,
#                       z-0.1, z+0.1)
   
   
   def __generate_grid( self, mins, lengths, cells, datas, names ):
      ''' Generates a grid from given data
          :param mins           An array of minimum coordinates for the grid for ex. [-100, 0, 0]
          :param lengths        An array of cell lengths (the cell's lengths in x, y, z direction)
          :param cells          An array of number of cells in x, y, z direction
          :param datas          Scalar data for the grid e.g. array([ cell1Rho, cell2Rho, cell3Rho, cell4Rho, .., cellNRho ])
      '''
      figure = mlab.gcf()
      mlab.clf()
      figure.scene.disable_render = True
      # Create nodes
      #x, y, z = mgrid[0:0.1:5j, 0:0.1:5j, 0:0.1*2.0/5.0:2j]
      x, y, z = mgrid[mins[0]:lengths[0]*(cells[0]+1):(cells[0]+1)*complex(0,1), mins[1]:lengths[1]*(cells[1]+1):(cells[1]+1)*complex(0,1), mins[2]:lengths[2]*(cells[2]+1):(cells[2]+1)*complex(0,1)]
      # Cell coordinates:
      x2 = 0.1*0.5 + np.arange(4)/4.0*0.1
      y2 = 0.1*0.5 + np.arange(4)/4.0*0.1
      z2 = 0.1*2.0/5.0*0.5 + np.arange(1)
      
      # Create points for the nodes:
      pts = empty(z.shape + (3,), dtype=float)
      pts[...,0] = x
      pts[...,1] = y
      pts[...,2] = z
      
      # Input scalars
      scalars = np.array(datas)
      # Input vectors
      #vectors = empty(z.shape + (3,), dtype=float)
      #vectors[...,0] = (4 - y*2)
      #vectors[...,1] = (x*3 - 12)
      #vectors[...,2] = sin(z*pi)
      
      # We reorder the points, scalars and vectors so this is as per VTK's
      # requirement of x first, y next and z last.
      pts = pts.transpose(2, 1, 0, 3).copy()
      pts.shape = pts.size/3, 3
      scalars = scalars.T.copy()
      #vectors = vectors.transpose(2, 1, 0, 3).copy()
      #vectors.shape = vectors.size/3, 3
      
      # Create the dataset.
      sg = tvtk.StructuredGrid(dimensions=x.shape, points=pts)
      sg.cell_data.scalars = ravel(scalars.copy())
      sg.cell_data.scalars.name = names
      #sg.point_data.vectors = vectors
      #sg.point_data.vectors.name = 'velocity'
      
      
      # Visualize the data
      d = mlab.pipeline.add_dataset(sg)
      iso = mlab.pipeline.surface(d)
      picker = figure.on_mouse_pick( self.__picker_callback, type='cell' )
      picker.tolerance = 0
      #iso.contour.maximum_contour = 75.0
      #vec = mlab.pipeline.vectors(d)
      #vec.glyph.mask_input_points = True
      #vec.glyph.glyph.scale_factor = 1.5
      figure.scene.disable_render = False
      mlab.show()

   def __generate_velocity_grid( self, cellid ):
      '''Generates a velocity grid from a given spatial cell id
         :param cellid           The spatial cell's ID
      '''
      figure = mlab.gcf()
      mlab.clf()
      figure.scene.disable_render = True
      # Create nodes




   def load_grid( self, variable ):
      ''' Creates a grid and inputs scalar variables from a vlsv file
          :param variable        Name of the variable to plot
      '''
      # Get the cell params:
      mins = np.array([self.__vlsvReader.read_parameter("xmin"), self.__vlsvReader.read_parameter("ymin"), self.__vlsvReader.read_parameter("zmin")])
      cells = np.array([self.__vlsvReader.read_parameter("xcells_ini"), self.__vlsvReader.read_parameter("ycells_ini"), self.__vlsvReader.read_parameter("zcells_ini")])
      maxs = np.array([self.__vlsvReader.read_parameter("xmax"), self.__vlsvReader.read_parameter("ymax"), self.__vlsvReader.read_parameter("zmax")])
      lengths = (maxs - mins) / cells.astype(float)
      # Get the variables:
      index_for_cellid_dict = self.__vlsvReader.get_cellid_locations()
      variable_array = self.__vlsvReader.read_variables( name=variable )
      # Sort the dictionary by cell id
      import operator
      sorted_index_for_cellid_dict = sorted(index_for_cellid_dict.iteritems(), key=operator.itemgetter(0))
      # Add the variable values:
      variable_array_sorted = []
      for i in sorted_index_for_cellid_dict:
         variable_array_sorted.append(variable_array[i[1]])
      # Draw the grid:
      self.__generate_grid( mins=mins, lengths=lengths, cells=cells, datas=variable_array_sorted, variable )
   












