#include <cstdlib>

#include "mpiconversion.h"
#include "vlsvwriter2.h"
#include "vlscommon.h"

#include <fstream>
#include <iostream>

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
   MPI_Wait(&mpiRequest,MPI_STATUS_IGNORE);
   MPI_Cancel(&recvRequest);
   
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
   
   fileOpen = false;
   return true;
}

bool VLSVWriter::endMultiwrite(const std::string& tagName,const std::string& arrayName,const std::map<std::string,std::string>& attribs) {
   bool success = true;
   
   // Wait for nonblocking communications to complete:
   uint64_t totalBytes;
   if (myrank == masterRank) {
      if (MPI_Wait(&recvRequest,MPI_STATUS_IGNORE) != MPI_SUCCESS) success = false;
      totalBytes = offsetIn-offset;
      int source = N_processes-1;
      MPI_Irecv(&offsetIn,1,MPI_Type<MPI_Offset>(),source,0,comm,&recvRequest);
   }
   
   // Master writes footer tag:
   if (myrank == masterRank) {
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
      offset = offsetIn;
   }
   
   return success;
}

bool VLSVWriter::open(const std::string& fname,MPI_Comm comm,const int& MASTER_RANK) {
   this->comm = comm;
   masterRank = MASTER_RANK;
   MPI_Comm_rank(comm,&myrank);
   MPI_Comm_size(comm,&N_processes);
   
   // All processes in communicator comm open the same file:
   int accessMode = (MPI_MODE_WRONLY | MPI_MODE_CREATE);
   MPI_Info MPIinfo = MPI_INFO_NULL;
   
   MPI_File_delete(const_cast<char*>(fname.c_str()),MPI_INFO_NULL);
   if (MPI_File_open(comm,const_cast<char*>(fname.c_str()),accessMode,MPIinfo,&fileptr) != MPI_SUCCESS) return false;
   MPI_File_set_view(fileptr,0,MPI_BYTE,MPI_BYTE,const_cast<char*>("native"),MPI_INFO_NULL);
   fileOpen = true;
   fileName = fname;

   // Post a receive for data offset:
   mpiRequest = MPI_REQUEST_NULL;
   int source = myrank - 1;
   if (source == -1 ) source = N_processes-1;
   MPI_Irecv(&offsetIn,1,MPI_Type<MPI_Offset>(),source,0,comm,&recvRequest);
   
   // Master opens a separate file for writing the footer:
   if (myrank == masterRank) {
      xmlWriter = new MuXML();
      XMLNode* root = xmlWriter->getRoot();
      XMLNode* xmlnode = xmlWriter->addNode(root,"VLSV","");
   }
   
   // Master writes 2 64bit integers to the start of file:
   bool success = true;
   if (myrank == masterRank) {
      // Write file endianness to the first byte:
      uint64_t endianness = 0;
      unsigned char* ptr = reinterpret_cast<unsigned char*>(&endianness);
      ptr[0] = detectEndianness();
      if (MPI_File_write_at(fileptr,0,&endianness,1,MPI_Type<uint64_t>(),MPI_STATUS_IGNORE) != MPI_SUCCESS) success = false;
      if (MPI_File_write_at(fileptr,8,&endianness,1,MPI_Type<uint64_t>(),MPI_STATUS_IGNORE) != MPI_SUCCESS) success = false;
      offset += 2*sizeof(uint64_t);
   }
   
   return success;
}

bool VLSVWriter::startMultiwrite(const std::string& dataType,const uint64_t& arraySize,const uint64_t& vectorSize,const uint64_t& dataSize) {
   bool success = true;
   // All processes except the master receive the offset from process with rank = myrank-1
   // (master receives offset in endMultiWrite):
   if (myrank != masterRank) {
      if (MPI_Wait(&recvRequest,MPI_STATUS_IGNORE) != MPI_SUCCESS) success = false;
      offset = offsetIn;
      int source = myrank - 1;
      MPI_Irecv(&offsetIn,1,MPI_Type<MPI_Offset>(),source,0,comm,&recvRequest);
   }
   
   // Calculate the amount of data written by this process in bytes, 
   // and send offset to process myrank+1:
   myBytes = arraySize * vectorSize * dataSize;
   offsetOut = offset + myBytes;
   mpiRequest = MPI_REQUEST_NULL;
   int target = myrank+1;
   if (target == N_processes) target = 0;
   if (MPI_Wait(&mpiRequest,MPI_STATUS_IGNORE) != MPI_SUCCESS) success = false;
   if (MPI_Isend(&offsetOut,1,MPI_Type<MPI_Offset>(),target,0,comm,&mpiRequest) != MPI_SUCCESS) success = false;
   
   this->dataType = dataType;
   this->arraySize = arraySize;
   this->vectorSize = vectorSize;
   this->dataSize = dataSize;
   return success;
}

bool VLSVWriter::writeArray(const std::string& tagName,const std::string& arrayName,const std::map<std::string,std::string>& attribs,const uint64_t& arraySize,const uint64_t& vectorSize,const std::string& dataType,const uint64_t& dataSize,char* array) {
   bool success = true;
   // All processes except the master receive the offset from process with rank = myrank-1:
   if (myrank != masterRank) {
      if (MPI_Wait(&recvRequest,MPI_STATUS_IGNORE) != MPI_SUCCESS) success = false;
      offset = offsetIn;
      int source = myrank - 1;
      MPI_Irecv(&offsetIn,1,MPI_Type<MPI_Offset>(),source,0,comm,&recvRequest);
   }
   
   // Calculate the amount of data written by this process in bytes, and 
   // send the next process its offset:
   myBytes = arraySize * vectorSize * dataSize;
   offsetOut = offset + myBytes;

   int target = myrank+1;
   if (target == N_processes) target = 0;
   if (MPI_Wait(&mpiRequest,MPI_STATUS_IGNORE) != MPI_SUCCESS) success = false;
   if (MPI_Isend(&offsetOut,1,MPI_Type<MPI_Offset>(),target,0,comm,&mpiRequest) != MPI_SUCCESS) success = false;
   
   // Write this process's data:
   if (MPI_File_write_at_all(fileptr,offset,array,myBytes,MPI_BYTE,MPI_STATUS_IGNORE) != MPI_SUCCESS) success = false;

   // Master writes footer tag:
   if (myrank == masterRank) {
      if (MPI_Wait(&recvRequest,MPI_STATUS_IGNORE) != MPI_SUCCESS) success = false;
      uint64_t totalBytes = offsetIn-offset;
      int source = N_processes-1;
      MPI_Irecv(&offsetIn,1,MPI_Type<MPI_Offset>(),source,0,comm,&recvRequest);
      
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
      
      offset = offsetIn;
   }
   return success;
}












