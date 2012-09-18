/*
This file is part of Vlasiator.

Copyright 2010, 2011, 2012 Finnish Meteorological Institute

Vlasiator is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3
as published by the Free Software Foundation.

Vlasiator is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef VLSVWRITER2_H
#define VLSVWRITER2_H

#include <stdint.h>
#include <mpi.h>

#include "mpiconversion.h"
#include "muxml.h"

struct WriteUnit {
    char *array;     //  pointer to data to be written
    MPI_Datatype mpiType;     //  type of data that is written
    uint64_t begin;  // where to begin writing in file compared to current offset (bytes)
    uint64_t amount; // how much to write (elements)
};

class VLSVWriter {
 public:
   VLSVWriter();
   ~VLSVWriter();
   
   bool close();
   bool endMultiwrite(const std::string& tagName,const std::string& arrayName,const std::map<std::string,std::string>& attribs);
   bool open(const std::string& fname,MPI_Comm comm,const int& MASTERRANK);
   
   bool startMultiwrite(const std::string& dataType,const uint64_t& arraySize,const uint64_t& vectorSize,const uint64_t& dataSize);
   
   template<typename T> bool multiwriteArray(const uint64_t& amount,T* array);
   template<typename T> 
   bool writeArray(const std::string& tagName,const std::string& arrayName,const std::map<std::string,std::string>& attribs,
                   const uint64_t& arraySize,const uint64_t& vectorSize,T* array);
   bool writeArray(const std::string& tagName,const std::string& arrayName,
                   const std::map<std::string,std::string>& attribs,
                   const uint64_t& arraySize,const uint64_t& vectorSize,
                   const std::string& dataType,const uint64_t& dataSize,void* array);
   template<typename T> 
   bool writeArrayMaster(const std::string& tagName,const std::string& arrayName,const std::map<std::string,std::string>& attribs,
                         const uint64_t& arraySize,const uint64_t& vectorSize,T* array);
   
 private:
   MPI_Comm comm;            /**< MPI communicator which is writing to the file.*/
   int myRank;               /**< Rank of this process in communicator comm.*/
   int N_processes;          /**< Number of processes in communicator comm.*/
   int masterRank;           /**< Rank of master process in communicator comm.*/
   
   bool fileOpen;            /**< If true, a file has been successfully opened for writing.*/
   std::string fileName;     /**< Name of the output file.*/
   MPI_File fileptr;         /**< MPI file pointer to the output file.*/
   uint64_t myBytes;         /**< Number of bytes this process is writing to the current array.*/
   MPI_Offset offset;        /**< MPI offset into output file for this process.*/
    
   MuXML* xmlWriter;         /**< Pointer to XML writer, used for writing a footer to the VLSV file.*/
   
   // Variables needed for multi-write mode. Multi-write mode 
   // is used for arrays which are too large to be buffered.
   
    std::string dataType;
    uint64_t arraySize;
    uint64_t dataSize;
    uint64_t vectorSize;

    //store the multiwrite requests in this vector
    std::vector<WriteUnit> multiWriteUnits;
    
    
    uint64_t* bytesPerProcess; //array with N_processes elements. Used to gather myBytes
    MPI_Offset* offsets; //array with N_processes elements. Used to scatter offsets 
   
    template<typename T> std::string arrayDataType();
};

// Templates which give the string datatype corresponding to the given datatype:
template<> inline std::string VLSVWriter::arrayDataType<int8_t>() {return "int";}
template<> inline std::string VLSVWriter::arrayDataType<int16_t>() {return "int";}
template<> inline std::string VLSVWriter::arrayDataType<int32_t>() {return "int";}
template<> inline std::string VLSVWriter::arrayDataType<int64_t>() {return "int";}
template<> inline std::string VLSVWriter::arrayDataType<uint8_t>() {return "uint";}
template<> inline std::string VLSVWriter::arrayDataType<uint16_t>() {return "uint";}
template<> inline std::string VLSVWriter::arrayDataType<uint32_t>() {return "uint";}
template<> inline std::string VLSVWriter::arrayDataType<uint64_t>() {return "uint";}
template<> inline std::string VLSVWriter::arrayDataType<float>() {return "float";}
template<> inline std::string VLSVWriter::arrayDataType<double>() {return "float";}
template<> inline std::string VLSVWriter::arrayDataType<long double>() {return "float";}

template<typename T> inline bool VLSVWriter::multiwriteArray(const uint64_t& amount,T* array) {
   bool success = true;
   //record the write, but only do it in the end routine
   WriteUnit writeUnit;   
   writeUnit.array=(char*)array;
   writeUnit.mpiType=MPI_Type<T>();
   writeUnit.amount=amount*vectorSize;
   multiWriteUnits.push_back(writeUnit);
   return success; 
}



//FIXME, the actual implementation should be in this, or in the other non-template writeArray. Not in both....
template<typename T>
inline bool VLSVWriter::writeArray(const std::string& tagName,const std::string& arrayName,
                                   const std::map<std::string,std::string>& attribs,
				   const uint64_t& arraySize,const uint64_t& vectorSize,T* array) {
    bool success = true;

    //master gathers amount of data each process writes, computes offsets, and scatters them to others

    myBytes = arraySize * vectorSize *sizeof(T);
    MPI_Gather(&myBytes,sizeof(uint64_t),MPI_BYTE,
               bytesPerProcess,sizeof(uint64_t),MPI_BYTE,
               masterRank,comm);
    
    if (myRank == masterRank) {
        offsets[0]=offset; //rank 0 handles the starting point of this block of data
        for(int i=1;i<N_processes;i++)
            offsets[i]=offsets[i-1]+bytesPerProcess[i-1];
    }


    //scatter offsets so that everybody has the correct offset
    MPI_Scatter(offsets,sizeof(MPI_Offset),MPI_BYTE,&offset,sizeof(MPI_Offset),MPI_BYTE,
                masterRank,comm);    
     
    // Write this process's data:
#ifdef NO_WRITE_AT_ALL
    //avoid crayXt memory leak, do not use collective routine        
    if (MPI_File_write_at(fileptr,offset,array,arraySize*vectorSize,
                              MPI_Type<T>(),MPI_STATUS_IGNORE)!= MPI_SUCCESS)
        success = false;
#else
    if (MPI_File_write_at_all(fileptr,offset,array,arraySize*vectorSize,
                              MPI_Type<T>(),MPI_STATUS_IGNORE)!= MPI_SUCCESS)
        success = false;
#endif
    
   // Master writes footer tag:
   if (myRank == masterRank) {
      uint64_t totalBytes = 0;
      for(int i=0;i<N_processes;i++)
          totalBytes+=bytesPerProcess[i];
      
      
      XMLNode* root = xmlWriter->getRoot();
      XMLNode* xmlnode = xmlWriter->find("VLSV",root);
      XMLNode* node = xmlWriter->addNode(xmlnode,tagName,offset);
      for (std::map<std::string,std::string>::const_iterator it=attribs.begin(); it!=attribs.end(); ++it) 
	xmlWriter->addAttribute(node,it->first,it->second);
      xmlWriter->addAttribute(node,"name",arrayName);
      xmlWriter->addAttribute(node,"vectorsize",vectorSize);
      xmlWriter->addAttribute(node,"arraysize",totalBytes/sizeof(T)/vectorSize);
      xmlWriter->addAttribute(node,"datatype",arrayDataType<T>());
      xmlWriter->addAttribute(node,"datasize",sizeof(T));
      
   
       //move forward by the amount of written data
      offset += totalBytes ;
   }
   return success;
}



template<typename T>
inline bool VLSVWriter::writeArrayMaster(const std::string& tagName,const std::string& arrayName,
                                         const std::map<std::string,std::string>& attribs,
                                         const uint64_t& arraySize,const uint64_t& vectorSize,T* array) {
   bool success = true;
   if (myRank != masterRank) {
      success=false;
      return success;
   }
   myBytes = arraySize * vectorSize *sizeof(T);
   
   // Write this process's data:
   if (MPI_File_write_at(fileptr,offset,array,arraySize*vectorSize,
                         MPI_Type<T>(),MPI_STATUS_IGNORE)!= MPI_SUCCESS)
      success = false;

   uint64_t totalBytes = myBytes;

   XMLNode* root = xmlWriter->getRoot();
   XMLNode* xmlnode = xmlWriter->find("VLSV",root);
   XMLNode* node = xmlWriter->addNode(xmlnode,tagName,offset);
   for (std::map<std::string,std::string>::const_iterator it=attribs.begin(); it!=attribs.end(); ++it) 
      xmlWriter->addAttribute(node,it->first,it->second);
   xmlWriter->addAttribute(node,"name",arrayName);
   xmlWriter->addAttribute(node,"vectorsize",vectorSize);
   xmlWriter->addAttribute(node,"arraysize",totalBytes/sizeof(T)/vectorSize);
   xmlWriter->addAttribute(node,"datatype",arrayDataType<T>());
   xmlWriter->addAttribute(node,"datasize",sizeof(T));
   
   //move forward by the amount of written data
   offset += totalBytes ;
   return success;
}

#endif
