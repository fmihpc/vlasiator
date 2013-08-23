import struct
import xml.etree.ElementTree as ET
import ast
import numpy as np


class VlsvFile(object):
   ''' Class for reading VLSV files
   ''' 
   def __init__(self, file_name):
      self.__file_name = file_name
      self.__xml_root = ET.fromstring("<VLSV></VLSV>")
      self.__fileindex_for_cellid=dict()
      self.__read_xml_footer()
      # Check if the file is using new or old vlsv format
      self.__uses_new_vlsv_format = False
      for child in self.__xml_root:
         if child.tag == "PARAMETER":
           if child.attrib["name"] == "version":
             self.__uses_new_vlsv_format = True
      self.__read_fileindex_for_cellid()
      # Read parameters (Note: Reading the spatial cell locations and storing them will anyway take the most time and memory):
      self.__vxblocks = (int)(self.read_parameter("vxblocks_ini"))
      self.__vyblocks = (int)(self.read_parameter("vyblocks_ini"))
      self.__vzblocks = (int)(self.read_parameter("vzblocks_ini"))

      self.__xcells = (int)(self.read_parameter("xcells_ini"))
      self.__ycells = (int)(self.read_parameter("ycells_ini"))
      self.__zcells = (int)(self.read_parameter("zcells_ini"))

      self.__xmin = self.read_parameter("xmin")
      self.__ymin = self.read_parameter("ymin")
      self.__zmin = self.read_parameter("zmin")
      self.__xmax = self.read_parameter("xmax")
      self.__ymax = self.read_parameter("ymax")
      self.__zmax = self.read_parameter("zmax")

      self.__vxmin = self.read_parameter("vxmin")
      self.__vymin = self.read_parameter("vymin")
      self.__vzmin = self.read_parameter("vzmin")
      self.__vxmax = self.read_parameter("vxmax")
      self.__vymax = self.read_parameter("vymax")
      self.__vzmax = self.read_parameter("vzmax")

      self.__dx = (self.__xmax - self.__xmin) / (float)(self.__xcells)
      self.__dy = (self.__ymax - self.__ymin) / (float)(self.__xcells)
      self.__dz = (self.__zmax - self.__zmin) / (float)(self.__xcells)

      # Velocity cell lengths
      velocity_cells_per_direction = 4
      self.__dvx = ((self.__vxmax - self.__vxmin) / (float)(self.__vxblocks)) / (float)(velocity_cells_per_direction)
      self.__dvy = ((self.__vymax - self.__vymin) / (float)(self.__vyblocks)) / (float)(velocity_cells_per_direction)
      self.__dvz = ((self.__vzmax - self.__vzmin) / (float)(self.__vzblocks)) / (float)(velocity_cells_per_direction)


   def __read_xml_footer(self):
      ''' Reads in the XML footer of the VLSV file and store all the content
      ''' 
      max_xml_size = 1000000
      fptr = open(self.__file_name,"rb")
      #(endianness,) = struct.unpack("c", fptr.read(1))
      fptr.seek(8)
      (offset,) = struct.unpack("Q", fptr.read(8))
      fptr.seek(offset)
      xml_data = fptr.read(max_xml_size)
      fptr.close() 
      (xml_string,) = struct.unpack("%ds" % len(xml_data), xml_data)
      self.__xml_root = ET.fromstring(xml_string)

   def __read_fileindex_for_cellid_old(self):
      """ Read in the cell ids and create an internal dictionary to give the index of an arbitrary cellID
      """
      cellids=self.read(name="SpatialGrid",tag="MESH")
      for index,cellid in enumerate(np.atleast_1d(cellids)):
         self.__fileindex_for_cellid[cellid]=index

   def __read_fileindex_for_cellid(self):
      """ Read in the cell ids and create an internal dictionary to give the index of an arbitrary cellID
      """
      if self.__uses_new_vlsv_format == True:
         cellids=self.read(mesh="SpatialGrid",name="CellID", tag="VARIABLE")
         for index,cellid in enumerate(cellids):
            self.__fileindex_for_cellid[cellid]=index
      else:
         # Uses old format
         self.__read_fileindex_for_cellid_old()

   def __list_old(self):
      ''' Print out a description of the content of the file. Useful
         for interactive usage
      '''
      print "tag = PARAMETERS"
      for child in self.__xml_root:
         if child.tag == "PARAMETERS":
            print "   ", child.attrib["name"], " = ", child.attrib["value"]
      print "tag = VARIABLE"
      for child in self.__xml_root:
         if child.tag == "VARIABLE":
            print "   ", child.attrib["name"]
      print "Other:"
      for child in self.__xml_root:
         if child.tag != "PARAMETERS" and child.tag != "VARIABLE":
            print "    tag = ", child.tag, " name = ", child.attrib["name"]
      return

   def __read_parameter_old(self, name):
      return self.read(name=name, tag="PARAMETERS")

   def __read_blocks_old(self, cell_id):
      ''' Read raw block data from the open file. Not so useful yet as 
         there is no information on the location of each block.
      
      Arguments:
      :param cell_id Cell ID of the cell whose velocity blocks are read
      :returns numpy array with blocks in the cell. Empty if cell has no stored blocks.
      '''
         
      #these two arrays are in the same order: 
      #list of cells for which dist function is saved
      cells_with_blocks = self.read("SpatialGrid","CELLSWITHBLOCKS")
      #number of blocks in each cell for which data is stored
      blocks_per_cell = self.read("SpatialGrid","BLOCKSPERCELL")
      (cells_with_blocks_index,) = np.where(cells_with_blocks == cell_id)

      if len(cells_with_blocks_index) == 0:
         #block data did not exist
         return []

      num_of_blocks = blocks_per_cell[cells_with_blocks_index[0]]

      for child in self.__xml_root:
         if child.attrib["name"] == "avgs":
            vector_size = ast.literal_eval(child.attrib["vectorsize"])
            #array_size = ast.literal_eval(child.attrib["arraysize"])
            element_size = ast.literal_eval(child.attrib["datasize"])
            datatype = child.attrib["datatype"]
            
            offset = ast.literal_eval(child.text)
            
            for i in range(0, cells_with_blocks_index[0]):
               offset += blocks_per_cell[i]*vector_size*element_size

            fptr = open(self.__file_name,"rb")
            fptr.seek(offset)
            if datatype == "float" and element_size == 4:
               data = np.fromfile(fptr, dtype = np.float32, count = vector_size*num_of_blocks)
            if datatype == "float" and element_size == 8:
               data = np.fromfile(fptr, dtype = np.float64, count = vector_size*num_of_blocks)
            fptr.close()
            return data.reshape(num_of_blocks, vector_size)   

   
      return []

   def list(self):
      ''' Print out a description of the content of the file. Useful
         for interactive usage
      '''
      if self.__uses_new_vlsv_format == False:
         return self.__list_old()
      print "tag = PARAMETER"
      for child in self.__xml_root:
         if child.tag == "PARAMETER":
            print "   ", child.attrib["name"]
      print "tag = VARIABLE"
      for child in self.__xml_root:
         if child.tag == "VARIABLE":
            print "   ", child.attrib["name"]
      print "tag = MESH"
      for child in self.__xml_root:
         if child.tag == "MESH":
            print "   ", child.attrib["name"]
      print "Other:"
      for child in self.__xml_root:
         if child.tag != "PARAMETER" and child.tag != "VARIABLE" and child.tag != "MESH":
            print "    tag = ", child.tag, " mesh = ", child.attrib["mesh"]

   def get_cellid_locations(self):
      return self.__fileindex_for_cellid

   def read(self, name="", tag="", mesh="", read_single_cellid=-1):
      ''' Read data from the open vlsv file. 
      
      Arguments:
      :param name Name of the data array
      :param tag  Tag of the data array.
      :param read_single_cellid  If -1 then all data is read. If nonzero then only the vector for the specified cell id is read
      :returns numpy array with the data

      '''
      if tag == "" and name == "" and tag == "":
         print "Bad arguments at read"
      #TODO, read_single_cellid should perhaps be an list/numpy array with cellids that are read in. This could be more efficient to 
      #     study multiple cells, e.g., along a line
      for child in self.__xml_root:
         if tag != "":
            if child.tag != tag:
               continue
         if name != "":
            if child.attrib["name"] != name:
               continue
         if mesh != "":
            if child.attrib["mesh"] != mesh:
               continue
         if child.tag == tag and child.attrib["name"] == name:
            vector_size = ast.literal_eval(child.attrib["vectorsize"])
            array_size = ast.literal_eval(child.attrib["arraysize"])
            element_size = ast.literal_eval(child.attrib["datasize"])
            datatype = child.attrib["datatype"]            
            offset = ast.literal_eval(child.text)
            if read_single_cellid >= 0:
               offset=offset+self.__fileindex_for_cellid[read_single_cellid]*element_size*vector_size
               array_size=1

            fptr = open(self.__file_name, "rb")
            fptr.seek(offset)

            if datatype == "float" and element_size == 4:
               data = np.fromfile(fptr, dtype = np.float32, count=vector_size*array_size)
            if datatype == "float" and element_size == 8:
               data = np.fromfile(fptr, dtype=np.float64, count=vector_size*array_size)
            if datatype == "int" and element_size == 4:
               data = np.fromfile(fptr, dtype=np.int32, count=vector_size*array_size)
            if datatype == "int" and element_size == 8:
               data = np.fromfile(fptr, dtype=np.int64, count=vector_size*array_size)
            if datatype == "uint" and element_size == 4:
               data = np.fromfile(fptr, dtype=np.uint32, count=vector_size*array_size)
            if datatype == "uint" and element_size == 8:
               data = np.fromfile(fptr, dtype=np.uint64, count=vector_size*array_size)
            fptr.close() 

            if vector_size > 1:
               data=data.reshape(array_size, vector_size)
            
            if array_size == 1:
               return data[0]
            else:
               return data

   def read_variables(self, name):
      ''' Read variables from the open vlsv file. 
      
      Arguments:
      :param name Name of the variable
      :returns numpy array with the data

      '''
      return self.read(mesh="SpatialGrid", name=name, tag="VARIABLE", read_single_cellid=-1)

   def read_variables_for_cellids(self, name, cellids, index=3):
      ''' Read variables from the open vlsv file. 
      
      Arguments:
      :param name Name of the variable
      :param cellids List of cellids
      :returns numpy array with the data
      '''
      # Read the variable:
      variablelist = self.read_variables(name)
      #Pick the variables with the cell ids in the list:
      returnvariablelist = []
      for cellid in cellids:
         if index == 3:
            returnvariablelist.append(np.array(variablelist[self.__fileindex_for_cellid[cellid]]))
         else:
            returnvariablelist.append(np.array(variablelist[self.__fileindex_for_cellid[cellid]][index]))
      # Return the variables:
      return np.array(returnvariablelist)

   def get_cellid(self, coordinates):
      ''' Returns the cell id at given coordinates

      Arguments:
      :param coordinates        The cell's coordinates
      :returns the cell id
      '''
      # Get xmax, xmin and xcells_ini
      xmax = self.read_parameter(name="xmax")
      xmin = self.read_parameter(name="xmin")
      xcells = (int)(self.read_parameter(name="xcells_ini"))
      # Do the same for y
      ymax = self.read_parameter(name="ymax")
      ymin = self.read_parameter(name="ymin")
      ycells = (int)(self.read_parameter(name="ycells_ini"))
      # And for z
      zmax = self.read_parameter(name="zmax")
      zmin = self.read_parameter(name="zmin")
      zcells = (int)(self.read_parameter(name="zcells_ini"))
   
      # Get cell lengths:
      cell_lengths = np.array([(xmax - xmin)/(float)(xcells), (ymax - ymin)/(float)(ycells), (zmax - zmin)/(float)(zcells)])
   
      # Get cell indices:
      cellindices = np.array([(int)((coordinates[0] - xmin)/(float)(cell_lengths[0])), (int)((coordinates[1] - ymin)/(float)(cell_lengths[1])), (int)((coordinates[2] - zmin)/(float)(cell_lengths[2]))])
      # Get the cell id:
      cellid = cellindices[0] + cellindices[1] * xcells + cellindices[2] * xcells * ycells + 1
      return cellid

   def get_cell_coordinates(self, cellid):
      ''' Returns a given cell's coordinates as a numpy array

      Arguments:
      :param cellid            The cell's ID
      :returns a numpy array with the coordinates
      '''
      # Get xmax, xmin and xcells_ini
      xmax = self.read_parameter(name="xmax")
      xmin = self.read_parameter(name="xmin")
      xcells = (int)(self.read_parameter(name="xcells_ini"))
      # Do the same for y
      ymax = self.read_parameter(name="ymax")
      ymin = self.read_parameter(name="ymin")
      ycells = (int)(self.read_parameter(name="ycells_ini"))
      # And for z
      zmax = self.read_parameter(name="zmax")
      zmin = self.read_parameter(name="zmin")
      zcells = (int)(self.read_parameter(name="zcells_ini"))
      # Get cell lengths:
      cell_lengths = np.array([(xmax - xmin)/(float)(xcells), (ymax - ymin)/(float)(ycells), (zmax - zmin)/(float)(zcells)])
      # Get cell indices:
      cellid = (int)(cellid - 1)
      cellindices = np.zeros(3)
      cellindices[0] = (int)(cellid)%(int)(xcells)
      cellindices[1] = ((int)(cellid)/(int)(xcells))%(int)(ycells)
      cellindices[2] = (int)(cellid)/(int)(xcells*ycells)
   
      # Get cell coordinates:
      cellcoordinates = np.zeros(3)
      cellcoordinates[0] = xmin + (cellindices[0] + 0.5) * cell_lengths[0]
      cellcoordinates[1] = ymin + (cellindices[1] + 0.5) * cell_lengths[1]
      cellcoordinates[2] = zmin + (cellindices[2] + 0.5) * cell_lengths[2]
      # Return the coordinates:
      return np.array(cellcoordinates)

   def get_velocity_cell_coordinates(self, vcellids):
      ''' Returns a given velocity cell's coordinates as a numpy array

      Arguments:
      :param vcellid       The velocity cell's ID
      :return a numpy array with the coordinates
      '''
      vcellids = np.atleast_1d(vcellids)
      # Get block ids:
      blocks = vcellids.astype(int) / 64
      # Get block coordinates:
      blockIndicesX = np.remainder(blocks.astype(int), (int)(self.__vxblocks))
      blockIndicesY = np.remainder(blocks.astype(int)/(int)(self.__vxblocks), (int)(self.__vyblocks))
      blockIndicesZ = blocks.astype(int)/(int)(self.__vxblocks*self.__vyblocks)
      blockCoordinatesX = blockIndicesX.astype(float) * self.__dvx * 4 + self.__vxmin
      blockCoordinatesY = blockIndicesY.astype(float) * self.__dvy * 4 + self.__vymin
      blockCoordinatesZ = blockIndicesZ.astype(float) * self.__dvz * 4 + self.__vzmin
      # Get cell indices:
      cellids = np.remainder(vcellids.astype(int), (int)(64))
      cellIndicesX = np.remainder(cellids.astype(int), (int)(4))
      cellIndicesY = np.remainder((cellids.astype(int)/(int)(4)).astype(int), (int)(4))
      cellIndicesZ = cellids.astype(int)/(int)(16)
      # Get cell coordinates:
      cellCoordinates = np.array([blockCoordinatesX.astype(float) + (cellIndicesX.astype(float) + 0.5) * self.__dvx,
                                  blockCoordinatesY.astype(float) + (cellIndicesY.astype(float) + 0.5) * self.__dvy,
                                  blockCoordinatesZ.astype(float) + (cellIndicesZ.astype(float) + 0.5) * self.__dvz])

#      # Get block id:
#      block = (int)(vcellid) / 64
#      # Get block coordinates:
#      blockIndices = np.array( [(int)(block)%(int)(self.__vxblocks), (int)((int)(block)/(int)(self.__vxblocks))%(int)(self.__vyblocks), (int)(block)/(int)(self.__vxblocks*self.__vyblocks)] )
#      blockCoordinates = np.array([self.__vxmin + 4 * self.__dvx * blockIndices[0],
#                                   self.__vymin + 4 * self.__dvy * blockIndices[1],
#                                   self.__vzmin + 4 * self.__dvz * blockIndices[2]])
#      #vcellid = 64 * velocity_block_id + kv*4*4 + jv*4 + iv
#      # Get cell coordinates:
#      cellIndices = np.array([(int)((int)(vcellid)%64)%4, (int)(((int)((int)(vcellid)%64)/4))%4, (int)((int)(vcellid)%64)/(int)(4*4)])
#      cellCoordinates = np.array([(cellIndices[0] + 0.5) * self.__dvx, (cellIndices[1] + 0.5) * self.__dvy, (cellIndices[2] + 0.5) * self.__dvz])
#      # Get the coordinates:
#      vcellCoordinates = blockCoordinates + cellCoordinates
      # Return cell coordinates:
      return cellCoordinates.transpose()


   def read_variable(self, name, cellid):
      ''' Read a variable of a given cell from the open vlsv file. 
      
      Arguments:
      :param name Name of the variable
      :param cellid Cell's cell id
      :returns numpy array with the data

      '''
      return self.read(mesh="SpatialGrid", name=name, tag="VARIABLE", read_single_cellid=cellid)

   def read_parameter(self, name):
      ''' Read a parameter from the vlsv file

      :param name   Name of the parameter
      '''
      if self.__uses_new_vlsv_format == True:
         return self.read(name=name, tag="PARAMETER")
      else:
         # Uses old format
         return self.__read_parameter_old(name=name)

   def read_blocks(self, cell_id):
      ''' Read raw block data from the open file. Not so useful yet as 
         there is no information on the location of each block.
      
      Arguments:
      :param cell_id Cell ID of the cell whose velocity blocks are read
      :returns numpy array with blocks in the cell. Empty if cell has no stored blocks.
      '''
      if self.__uses_new_vlsv_format == False:
         # Uses old format
         return self.__read_blocks_old(cell_id)
      #these two arrays are in the same order: 
      #list of cells for which dist function is saved
      cells_with_blocks = self.read("","CELLSWITHBLOCKS")
      #number of blocks in each cell for which data is stored
      blocks_per_cell = self.read("","BLOCKSPERCELL")
      (cells_with_blocks_index,) = np.where(cells_with_blocks == cell_id)

      if len(cells_with_blocks_index) == 0:
         #block data did not exist
         return []

      num_of_blocks = blocks_per_cell[cells_with_blocks_index[0]]

      for child in self.__xml_root:
         if child.attrib["name"] == "avgs":
            vector_size = ast.literal_eval(child.attrib["vectorsize"])
            #array_size = ast.literal_eval(child.attrib["arraysize"])
            element_size = ast.literal_eval(child.attrib["datasize"])
            datatype = child.attrib["datatype"]
            
            offset = ast.literal_eval(child.text)
            
            for i in range(0, cells_with_blocks_index[0]):
               offset += blocks_per_cell[i]*vector_size*element_size

            fptr = open(self.__file_name,"rb")
            fptr.seek(offset)
            if datatype == "float" and element_size == 4:
               data = np.fromfile(fptr, dtype = np.float32, count = vector_size*num_of_blocks)
            if datatype == "float" and element_size == 8:
               data = np.fromfile(fptr, dtype = np.float64, count = vector_size*num_of_blocks)
            fptr.close()
            return data.reshape(num_of_blocks, vector_size)   

   
      return []


#getVelocityBlockCoordinates( cellStruct, blockId, blockCoordinates )
#void getVelocityBlockCoordinates(const CellStructure & cellStruct, const uint64_t block, array<Real, 3> & coordinates ) {
#   //First get indices:
#   array<uint64_t, 3> blockIndices;
#   blockIndices[0] = block % cellStruct.vcell_bounds[0];
#   blockIndices[1] = (block / cellStruct.vcell_bounds[0]) % cellStruct.vcell_bounds[1];
#   blockIndices[2] = block / (cellStruct.vcell_bounds[0] * cellStruct.vcell_bounds[1]);
#   //Store the coordinates:
#   for( int i = 0; i < 3; ++i ) {
#      coordinates[i] = cellStruct.min_vcoordinates[i] + cellStruct.vblock_length[i] * blockIndices[i];
#   }
#   return;
#}
#if( vlsvReader.readParameter( "vxblocks_ini", vcell_bounds[0] ) == false ) cerr << "FAILED TO READ PARAMETER AT " << __FILE__ << " " << __LINE__ << endl;
#//y-direction
#if( vlsvReader.readParameter( "vyblocks_ini", vcell_bounds[1] ) == false ) cerr << "FAILED TO READ PARAMETER AT " << __FILE__ << " " << __LINE__ << endl;
#//z-direction
#if( vlsvReader.readParameter( "vzblocks_ini", vcell_bounds[2] ) == false ) cerr << "FAILED TO READ PARAMETER AT " << __FILE__ << " " << __LINE__ << endl;
#for child in self.__xml_root:
#   if tag != "":
#      if child.tag != tag:
#         continue
#   if name != "":
#      if child.attrib["name"] != name:
#         continue
#   if mesh != "":
#      if child.attrib["mesh"] != mesh:
#         continue
#   if child.tag == tag and child.attrib["name"] == name:
#      vector_size = ast.literal_eval(child.attrib["vectorsize"])
#      array_size = ast.literal_eval(child.attrib["arraysize"])
#      element_size = ast.literal_eval(child.attrib["datasize"])
#      datatype = child.attrib["datatype"]            
#      offset = ast.literal_eval(child.text)
#      if read_single_cellid >= 0:
#         offset=offset+self.__fileindex_for_cellid[read_single_cellid]*element_size*vector_size
#         array_size=1
#
#      fptr = open(self.__file_name, "rb")
#      fptr.seek(offset)
#
#      if datatype == "float" and element_size == 4:
#         data = np.fromfile(fptr, dtype = np.float32, count=vector_size*array_size)
#      if datatype == "float" and element_size == 8:
#         data = np.fromfile(fptr, dtype=np.float64, count=vector_size*array_size)
#      if datatype == "int" and element_size == 4:
#         data = np.fromfile(fptr, dtype=np.int32, count=vector_size*array_size)
#      if datatype == "int" and element_size == 8:
#         data = np.fromfile(fptr, dtype=np.int64, count=vector_size*array_size)
#      if datatype == "uint" and element_size == 4:
#         data = np.fromfile(fptr, dtype=np.uint32, count=vector_size*array_size)
#      if datatype == "uint" and element_size == 8:
#         data = np.fromfile(fptr, dtype=np.uint64, count=vector_size*array_size)
#      fptr.close() 
#
#      if vector_size > 1:
#         data=data.reshape(array_size, vector_size)
#      
#      if array_size == 1:
#         return data[0]
#      else:
#         return data

   def __read_velocity_cells_old( self, cellid, cells_with_blocks, blocks_per_cell, cells_with_blocks_index  ):
      # Read in the coordinates:
      #block_coordinates = self.read(name="",tag="BLOCKCOORDINATES")
      # Navigate to the correct position:
      offset = 0
      for i in xrange(0, cells_with_blocks_index[0]):
         offset += blocks_per_cell[i]

      num_of_blocks = np.atleast_1d(blocks_per_cell)[cells_with_blocks_index[0]]

      # Read in avgs and velocity cell ids:
      for child in self.__xml_root:
         # Read in avgs
         if child.attrib["name"] == "avgs":
            vector_size = ast.literal_eval(child.attrib["vectorsize"])
            #array_size = ast.literal_eval(child.attrib["arraysize"])
            element_size = ast.literal_eval(child.attrib["datasize"])
            datatype = child.attrib["datatype"]

            # Navigate to the correct position
            offset_avgs = offset * vector_size * element_size + ast.literal_eval(child.text)
#            for i in range(0, cells_with_blocks_index[0]):
#               offset_avgs += blocks_per_cell[i]*vector_size*element_size

            fptr = open(self.__file_name,"rb")
            fptr.seek(offset_avgs)
            if datatype == "float" and element_size == 4:
               data_avgs = np.fromfile(fptr, dtype = np.float32, count = vector_size*num_of_blocks)
            if datatype == "float" and element_size == 8:
               data_avgs = np.fromfile(fptr, dtype = np.float64, count = vector_size*num_of_blocks)
            fptr.close()
            data_avgs = data_avgs.reshape(num_of_blocks, vector_size)
         # Read in block coordinates:
         if child.attrib["name"] == "SpatialGrid" and child.tag == "BLOCKCOORDINATES":
            vector_size = ast.literal_eval(child.attrib["vectorsize"])
            #array_size = ast.literal_eval(child.attrib["arraysize"])
            element_size = ast.literal_eval(child.attrib["datasize"])
            datatype = child.attrib["datatype"]

            offset_block_coordinates = offset * vector_size * element_size + ast.literal_eval(child.text)

            fptr = open(self.__file_name,"rb")
            fptr.seek(offset_block_coordinates)
            if datatype == "float" and element_size == 4:
               data_block_coordinates = np.fromfile(fptr, dtype = np.float32, count = vector_size*num_of_blocks)
            if datatype == "float" and element_size == 8:
               data_block_coordinates = np.fromfile(fptr, dtype = np.float64, count = vector_size*num_of_blocks)
            fptr.close()

            data_block_coordinates = data_block_coordinates.reshape(num_of_blocks, vector_size)

      # Check to make sure the sizes match (just some extra debugging)
      if len(data_avgs) != len(data_block_coordinates):
         print "BAD DATA SIZES"
      # Make a dictionary (hash map) out of velocity cell ids and avgs:
      velocity_cells = {}
      array_size = len(data_avgs)

      # Construct velocity cells:
      velocity_cell_ids = []
      for kv in xrange(4):
         for jv in xrange(4):
            for iv in xrange(4):
               velocity_cell_ids.append(kv*16 + jv*4 + iv)
#      velocity_cell_ids_temp = np.array(velocity_cell_ids)
#      # Construct velocity blocks:
#      dv = np.array([self.__dvx, self.__dvy, self.__dvz])
#      v_mins = np.array([self.__vx_min, self.__vy_min, self.__vz_min])
#      velocity_block_indices = np.array(np.floor(data_block_coordinates - v_mins) / (4*dv))
#      velocity_block_ids = velocity_block_indices * np.array([1, self.__vx_blocks, self.__vx_blocks*self.__vy_blocks])
#      velocity_cell_ids = np.array([


      for i in xrange(array_size):
         block_coordinate = data_block_coordinates[i]
         # The minimum corner coordinates of the blocks
         vx = block_coordinate[0]
         vy = block_coordinate[1]
         vz = block_coordinate[2]
         # The diff in blocks
#         dvx = block_coordinate[3]
#         dvy = block_coordinate[4]
#         dvz = block_coordinate[5]
         avgs = data_avgs[i]
         # Get the velocity cell id (First transform coordinates to block indices, then block indices to block id and then block id to velocity cell ids):
         velocity_block_indices = np.array([np.floor((vx - self.__vxmin) / (4*self.__dvx)), np.floor((vy - self.__vymin) / (4*self.__dvy)), np.floor((vz - self.__vzmin) / (4*self.__dvz))])
         velocity_block_id = velocity_block_indices[0] + velocity_block_indices[1] * self.__vxblocks + velocity_block_indices[2] * self.__vxblocks * self.__vyblocks
         avgIndex = 0

         for j in velocity_cell_ids + 64*velocity_block_id:
            velocity_cells[(int)(j)] = avgs[avgIndex]
            avgIndex = avgIndex + 1
#         for kv in xrange(4):
#            for jv in xrange(4):
#               for iv in xrange(4):
#                  #appending = np.array([iv, jv, kv])
#                  #velocity_cell_indices = velocity_block_indices * 4 + appending
#                  #Note: There  are 64 velocity cells in every block and 4 in every direction (4 times 4 times 4 = 64)
#                  #vcellid = 64 * velocity_block_id + kv*4*4 + jv*4 + iv
#                  velocity_cells[(int)(64 * velocity_block_id + kv*16 + jv*4 + iv)] = avgs[avgIndex]
#                  avgIndex = avgIndex + 1
      return velocity_cells
      
#const velocity_block_indices_t indices = {{
#   (unsigned int) np.floor((vx - SpatialCell::vx_min) / SpatialCell::block_dvx),
#   (unsigned int) np.floor((vy - SpatialCell::vy_min) / SpatialCell::block_dvy),
#   (unsigned int) np.floor((vz - SpatialCell::vz_min) / SpatialCell::block_dvz)
#}};
#indices[0] = cell % block_vx_length;
#indices[1] = (cell / block_vx_length) % block_vy_length;
#indices[2] = cell / (block_vx_length * block_vy_length);
#static unsigned int get_velocity_block(
#   const Real vx,
#   const Real vy,
#   const Real vz
#) {
#   if (vx < SpatialCell::vx_min || vx >= SpatialCell::vx_max
#       || vy < SpatialCell::vy_min || vy >= SpatialCell::vy_max
#       || vz < SpatialCell::vz_min || vz >= SpatialCell::vz_max) {
#      return error_velocity_block;
#   }
#
#   const velocity_block_indices_t indices = {{
#      (unsigned int) np.floor((vx - SpatialCell::vx_min) / SpatialCell::block_dvx),
#      (unsigned int) np.floor((vy - SpatialCell::vy_min) / SpatialCell::block_dvy),
#      (unsigned int) np.floor((vz - SpatialCell::vz_min) / SpatialCell::block_dvz)
#   }};
#
#   return SpatialCell::get_velocity_block(indices);
#}

   def read_velocity_cells(self, cellid):
      ''' Read velocity cells from a spatial cell
      
      Arguments:
      :param cellid Cell ID of the cell whose velocity cells are read
      :returns numpy array with blocks in the cell. Empty if cell has no stored blocks.
      '''
      #these two arrays are in the same order: 
      #list of cells for which dist function is saved
      cells_with_blocks = self.read("SpatialGrid","CELLSWITHBLOCKS")
      #number of blocks in each cell for which data is stored
      blocks_per_cell = self.read("SpatialGrid","BLOCKSPERCELL")
      (cells_with_blocks_index,) = np.where(cells_with_blocks == cellid)

      if len(cells_with_blocks_index) == 0:
         #block data did not exist
         print "Cell does not have velocity distribution"
         return []

      num_of_blocks = np.atleast_1d(blocks_per_cell)[cells_with_blocks_index[0]]

      # Check for the old library
      if self.__uses_new_vlsv_format == False:
         # Uses old format
         return self.__read_velocity_cells_old(cellid=cellid, cells_with_blocks=cells_with_blocks, blocks_per_cell=blocks_per_cell, cells_with_blocks_index=cells_with_blocks_index)
      else:
         # Uses new format:
         print "NOT IMPLEMENTED YET"
         return self.__read_velocity_cells_new(cellid=cellid, cells_with_blocks=cells_with_blocks, blocks_per_cell=blocks_per_cell, cells_with_blocks_index=cells_with_blocks_index)










