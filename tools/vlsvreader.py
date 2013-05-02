import struct
import sys
import os
import xml.etree.ElementTree as ET
import ast
import numpy as np
    

class File(object):

    def __init__(self,file_name):
        self.__file_name=file_name
        self.__parse_xml_footer()


    def __parse_xml_footer(self):
        max_xml_size=1000000;
        f=open(self.__file_name,"rb")
        (endianness,)=struct.unpack("c",f.read(1))
        f.seek(8)
        (offset,)=struct.unpack("Q",f.read(8))
        f.seek(offset)
        xml_data=f.read(max_xml_size)
        f.close() 
        (xml_string,)=struct.unpack("%ds" % len(xml_data),xml_data)
        self.__xml_root=ET.fromstring(xml_string)


    def ls(self):
        print "tag = PARAMETERS"
        for child in self.__xml_root:
            if child.tag=="PARAMETERS":
                print "    ",child.attrib["name"]," = ",child.attrib["value"]
        print "tag = VARIABLE"
        for child in self.__xml_root:
            if child.tag=="VARIABLE":
                print "    ",child.attrib["name"]
        print "Other:"
        for child in self.__xml_root:
            if child.tag!="PARAMETERS" and child.tag!="VARIABLE":
                print "     tag = ",child.tag, " name = ", child.attrib["name"]

    
    def read(self,name,tag="VARIABLE"):
        for child in self.__xml_root:
            if child.tag==tag and child.attrib["name"]==name:
                vector_size=ast.literal_eval(child.attrib["vectorsize"])
                array_size=ast.literal_eval(child.attrib["arraysize"])
                element_size=ast.literal_eval(child.attrib["datasize"])
                datatype=child.attrib["datatype"]                
                offset=ast.literal_eval(child.text)

                fd=open(self.__file_name,"rb")
                fd.seek(offset)

                if datatype=="float" and element_size==4:
                    data=np.fromfile(fd,dtype=np.float32,count=vector_size*array_size);
                if datatype=="float" and element_size==8:
                    data=np.fromfile(fd,dtype=np.float64,count=vector_size*array_size);
                if datatype=="int" and element_size==4:
                    data=np.fromfile(fd,dtype=np.int32,count=vector_size*array_size);
                if datatype=="int" and element_size==8:
                    data=np.fromfile(fd,dtype=np.int64,count=vector_size*array_size);
                if datatype=="uint" and element_size==4:
                    data=np.fromfile(fd,dtype=np.uint32,count=vector_size*array_size);
                if datatype=="uint" and element_size==8:
                    data=np.fromfile(fd,dtype=np.uint64,count=vector_size*array_size);
                fd.close() 

                if vector_size > 1:
                    return data.reshape(array_size,vector_size)
                else:
                    if len(data) == 1:
                        return data[0]
                    else:
                        return data


    def read_blocks(self,cellID):
        #these two arrays are in the same order: 
        #list of cells for which dist function is saved
        cells_with_blocks=self.read("SpatialGrid","CELLSWITHBLOCKS")
        #number of blocks in each cell for which data is stored
        blocks_per_cell=self.read("SpatialGrid","BLOCKSPERCELL")
        (cells_with_blocks_index,)=np.where(cells_with_blocks==cellID)

        if len(cells_with_blocks_index)==0:
            #block data did not exist
            return []

        num_of_blocks=blocks_per_cell[cells_with_blocks_index[0]]

        for child in self.__xml_root:
            if child.attrib["name"]=="avgs":
                vector_size=ast.literal_eval(child.attrib["vectorsize"])
                array_size=ast.literal_eval(child.attrib["arraysize"])
                element_size=ast.literal_eval(child.attrib["datasize"])
                datatype=child.attrib["datatype"]
                
                offset=ast.literal_eval(child.text)
                
                for i in range(0,cells_with_blocks_index[0]):
                    offset+=blocks_per_cell[i]*vector_size*element_size

                fd=open(self.__file_name,"rb")
                fd.seek(offset)
                if datatype=="float" and element_size==4:
                    data=np.fromfile(fd,dtype=np.float32,count=vector_size*num_of_blocks)
                if datatype=="float" and element_size==8:
                    data=np.fromfile(fd,dtype=np.float64,count=vector_size*num_of_blocks)
                fd.close()
                return data.reshape(num_of_blocks,vector_size)    

    
        return []
            
