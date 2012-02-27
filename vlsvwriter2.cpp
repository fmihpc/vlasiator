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

#include <cstdlib>
#include <iostream>
#include <fstream>

#include "mpiconversion.h"
#include "vlsvwriter2.h"
#include "vlscommon.h"

using namespace std;

VLSVWriter::VLSVWriter() {
   fileOpen = false;
   offset = 0;
   xmlWriter = NULL;
}

VLSVWriter::~VLSVWriter() {
   delete xmlWriter;
   xmlWriter = NULL;
}

bool VLSVWriter::close() {
    if (fileOpen == false) return false; 
    MPI_File_close(&fileptr);
    
    // Master process appends footer to the end of binary file:
    if (myrank == masterRank) {

        fstream footer;
        footer.open(fileName.c_str(),fstream::out | fstream::app);
        footer.seekp(0,std::ios_base::end);
        
        // Get put position
        uint64_t footerOffset = footer.tellp();
        
        // Write footer:
        xmlWriter->print(footer);
        delete xmlWriter;
        xmlWriter = NULL;
        footer.close();
        
        // Write header position to the beginning of binary file:
        footer.open(fileName.c_str(),fstream::in | fstream::out | fstream::binary | fstream::ate);
        char* ptr = reinterpret_cast<char*>(&footerOffset);
        footer.seekp(sizeof(uint64_t));
        footer.write(ptr,sizeof(uint64_t));
        footer.close();
    }
    
    if (myrank == masterRank) {
        //delete some arrays
        delete[] offsets ;
        delete[] bytesPerProcess;
    }
    MPI_Barrier(comm);
    fileOpen = false;
    return true;
}

bool VLSVWriter::open(const std::string& fname,MPI_Comm comm,const int& MASTER_RANK) {
    //FIXME: should be MPI_COMM_DUP
    this->comm = comm;
    masterRank = MASTER_RANK;
    
    MPI_Comm_rank(comm,&myrank);
    MPI_Comm_size(comm,&N_processes);
    // All processes in communicator comm open the same file:
    int accessMode = (MPI_MODE_WRONLY | MPI_MODE_CREATE);
    MPI_Info MPIinfo = MPI_INFO_NULL;
    
    MPI_File_delete(const_cast<char*>(fname.c_str()),MPI_INFO_NULL);
    if (MPI_File_open(comm,const_cast<char*>(fname.c_str()),accessMode,MPIinfo,&fileptr)    
        != MPI_SUCCESS) 
        return false;
    MPI_File_set_view(fileptr,0,MPI_BYTE,MPI_BYTE,const_cast<char*>("native"),MPI_INFO_NULL);
    fileOpen = true;
    fileName = fname;
    offset = 0; //offset set to 0 when opening a new file

    
    //only master rank needs these arrays
    if (myrank == masterRank) {
        offsets=new MPI_Offset[N_processes];
        bytesPerProcess=new uint64_t[N_processes];    
    }
    
    // Master opens a separate file for writing the footer:
    if (myrank == masterRank) {
        xmlWriter = new MuXML();
        XMLNode* root = xmlWriter->getRoot();
        xmlWriter->addNode(root,"VLSV","");
    }
    
    // Master writes 2 64bit integers to the start of file:
    bool success = true;
    if (myrank == masterRank) {
        // Write file endianness to the first byte:
        uint64_t endianness = 0;
        unsigned char* ptr = reinterpret_cast<unsigned char*>(&endianness);
        ptr[0] = detectEndianness();
        if (MPI_File_write_at(fileptr,0,&endianness,1,MPI_Type<uint64_t>(),MPI_STATUS_IGNORE) 
			!= MPI_SUCCESS) 
			success = false;
        if (MPI_File_write_at(fileptr,8,&endianness,1,MPI_Type<uint64_t>(),MPI_STATUS_IGNORE)
		    != MPI_SUCCESS) 
		    success = false;	
        offset += 2*sizeof(uint64_t); //only master rank keeps a running count
    }
    
return success;
}

bool VLSVWriter::startMultiwrite(const std::string& dataType,const uint64_t& arraySize,
                                 const uint64_t& vectorSize,const uint64_t& dataSize) {
    bool success = true;
    
   // Calculate the amount of data written by this process in bytes, 
   // and send offset to process myrank+1:
   myBytes = arraySize * vectorSize * dataSize;
   MPI_Gather(&myBytes,sizeof(uint64_t),MPI_BYTE,
	      bytesPerProcess,sizeof(uint64_t),MPI_BYTE,
	      0,comm);
    
   if (myrank == 0) {
      offsets[0]=offset; //rank 0 handles the starting point of this block of data
      for(int i=1;i<N_processes;i++)
	offsets[i]=offsets[i-1]+bytesPerProcess[i-1];
   }
   
   //scatter offsets so that everybody has the correct offset
   MPI_Scatter(offsets,sizeof(MPI_Offset),MPI_BYTE,&offset,sizeof(MPI_Offset),MPI_BYTE,
	       0,comm);
      
   this->dataType = dataType;
   this->arraySize = arraySize;
   this->vectorSize = vectorSize;
   this->dataSize = dataSize;
   
   multiWriteUnits.clear();
   
   return success;
}
bool VLSVWriter::endMultiwrite(const std::string& tagName,const std::string& arrayName,
    const std::map<std::string,std::string>& attribs) {
    bool success = true;
    MPI_Barrier(comm);
    
    if (multiWriteUnits.size() >0 ){
       //create data type for data that we use to do all writes at once
       MPI_Aint* displacements = new MPI_Aint[multiWriteUnits.size()];
       MPI_Datatype* types = new MPI_Datatype[multiWriteUnits.size()];
       int* blockLengths = new int[multiWriteUnits.size()];
       for (uint i=0; i<multiWriteUnits.size(); ++i){
	  displacements[i] = multiWriteUnits[i].array - multiWriteUnits[0].array;
	  types[i]         = multiWriteUnits[0].mpiType;
	  blockLengths[i]  = multiWriteUnits[0].amount;
       }
       MPI_Datatype outputType;
       MPI_Type_create_struct(multiWriteUnits.size(),blockLengths,displacements,types,&outputType);
       MPI_Type_commit(&outputType);
       
       //write out actual data with one collective call       
       if(MPI_File_write_at_all(fileptr,offset,multiWriteUnits[0].array,
				1,outputType,MPI_STATUS_IGNORE) != MPI_SUCCESS)
	 success = false;
       //free type
       MPI_Type_free(&outputType);
       delete[] displacements;
       delete[] types;
       delete[] blockLengths;
    }
    else{
        //we have no data to write, write out zero data to participate in collective call
        if(MPI_File_write_at_all(fileptr,offset,NULL,
                                 0,MPI_BYTE,MPI_STATUS_IGNORE) != MPI_SUCCESS)
	 success = false;
    }
    
   // Master writes footer tag:
   if (myrank == masterRank) {
       uint64_t totalBytes = 0;
       for(int i=0;i<N_processes;i++)
           totalBytes+=bytesPerProcess[i];
       
       XMLNode* root = xmlWriter->getRoot();
       XMLNode* xmlnode = xmlWriter->find("VLSV",root);
       XMLNode* node = xmlWriter->addNode(xmlnode,tagName,offset);
       for (std::map<std::string,std::string>::const_iterator it=attribs.begin(); 
            it!=attribs.end(); ++it)
           xmlWriter->addAttribute(node,it->first,it->second);
       xmlWriter->addAttribute(node,"name",arrayName);
       xmlWriter->addAttribute(node,"vectorsize",vectorSize);
       xmlWriter->addAttribute(node,"arraysize",totalBytes/dataSize/vectorSize);
       xmlWriter->addAttribute(node,"datatype",dataType);
       xmlWriter->addAttribute(node,"datasize",dataSize);

       offset +=totalBytes; //update offset
   }
   return success;
}
  

bool VLSVWriter::writeArray(const std::string& tagName,const std::string& arrayName,
                            const std::map<std::string,std::string>& attribs,
                            const uint64_t& arraySize,const uint64_t& vectorSize,
                            const std::string& dataType,const uint64_t& dataSize,void* array) {
    bool success = true;
    // All processes except the master receive the offset from process with rank = myrank-1:
    // Calculate the amo unt of data written by this process in bytes, and 
   // send the next process its offset:
    myBytes = arraySize * vectorSize * dataSize;
    MPI_Gather(&myBytes,sizeof(uint64_t),MPI_BYTE,
               bytesPerProcess,sizeof(uint64_t),MPI_BYTE,
               0,comm);
    
    if (myrank == 0) {
        offsets[0]=offset; //rank 0 handles the starting point of this block of data
        for(int i=1;i<N_processes;i++)
            offsets[i]=offsets[i-1]+bytesPerProcess[i-1];
    }

   //scatter offsets so that everybody has the correct offset
   MPI_Scatter(offsets,sizeof(MPI_Offset),MPI_BYTE,&offset,sizeof(MPI_Offset),MPI_BYTE,0,comm);

    // Write this process's data:
   if (MPI_File_write_at_all(fileptr,offset,array,myBytes,MPI_BYTE,MPI_STATUS_IGNORE) != MPI_SUCCESS) {
      return false;
   }
   
   // Master writes footer tag:
   if (myrank == masterRank) {
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
       xmlWriter->addAttribute(node,"arraysize",totalBytes/dataSize/vectorSize);
       xmlWriter->addAttribute(node,"datatype",dataType);
       xmlWriter->addAttribute(node,"datasize",dataSize);
       //move forward by the amount of written data
       offset += totalBytes ;
   }
   return success;
}












