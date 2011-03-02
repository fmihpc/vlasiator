#ifndef VLSVWRITER2_H
#define VLSVWRITER2_H

#include <stdint.h>
#include <mpi.h>

#include "mpiconversion.h"
#include "muxml.h"

class VLSVWriter {
 public:
   VLSVWriter();
   ~VLSVWriter();
   
   bool close();
   bool endMultiwrite(const std::string& tagName,const std::string& arrayName,const std::map<std::string,std::string>& attribs);
   bool open(const std::string& fname,MPI_Comm comm,const int& MASTER_RANK);
   
   bool startMultiwrite(const std::string& dataType,const uint64_t& arraySize,const uint64_t& vectorSize,const uint64_t& dataSize);
   
   template<typename T> bool multiwriteArray(const uint64_t& begin,const uint64_t& amount,T* array);
   template<typename T> 
     bool writeArray(const std::string& tagName,const std::string& arrayName,const std::map<std::string,std::string>& attribs,
		     const uint64_t& arraySize,const uint64_t& vectorSize,T* array);
   bool writeArray(const std::string& tagName,const std::string& arrayName,const std::map<std::string,std::string>& attribs,
		   const uint64_t& arraySize,const uint64_t& vectorSize,const std::string& dataType,const uint64_t& dataSize,char* array);
 private:
   MPI_Comm comm;            /**< MPI communicator which is writing to the file.*/
   int myrank;               /**< Rank of this process in communicator comm.*/
   int N_processes;          /**< Number of processes in communicator comm.*/
   int masterRank;           /**< Rank of master process in communicator comm.*/
   
   bool fileOpen;            /**< If true, a file has been successfully opened for writing.*/
   std::string fileName;     /**< Name of the output file.*/
   MPI_File fileptr;         /**< MPI file pointer to the output file.*/
   uint64_t myBytes;         /**< Number of bytes this process is writing to the current array.*/
   MPI_Offset offset;        /**< MPI offset into output file for this process.*/
   MPI_Offset offsetIn;
   MPI_Offset offsetOut;
   MPI_Request mpiRequest;
   MPI_Request recvRequest;
   MuXML* xmlWriter;         /**< Pointer to XML writer, used for writing a footer to the VLSV file.*/
   
   // Variables needed for multi-write mode. Multi-write mode 
   // is used for arrays which are too large to be buffered.
   
   std::string dataType;
   uint64_t arraySize;
   uint64_t dataSize;
   uint64_t vectorSize;
   
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

template<typename T> inline bool VLSVWriter::multiwriteArray(const uint64_t& begin,const uint64_t& amount,T* array) {
   bool success = true;
   if (MPI_File_write_at(fileptr,offset+begin*vectorSize*sizeof(T),array,amount*vectorSize,MPI_Type<T>(),MPI_STATUS_IGNORE) != MPI_SUCCESS) success = false;
   return success;
}

template<typename T>
inline bool VLSVWriter::writeArray(const std::string& tagName,const std::string& arrayName,const std::map<std::string,std::string>& attribs,
				   const uint64_t& arraySize,const uint64_t& vectorSize,T* array) {
   bool success = true;
   // All processes except the master receive the offset from process with rank = myrank-1, 
   if (myrank != masterRank) {
      if (MPI_Wait(&recvRequest,MPI_STATUS_IGNORE) != MPI_SUCCESS) success = false;
      offset = offsetIn;
      int source = myrank - 1;
      MPI_Irecv(&offsetIn,1,MPI_Type<MPI_Offset>(),source,0,comm,&recvRequest);
   }
   
   // Calculate the amount of data written by this process in bytes, and 
   // send the next process its offset:
   uint64_t myBytes = arraySize * vectorSize * sizeof(T);
   offsetOut = offset+myBytes;

   int target = myrank+1;
   if (target == N_processes) target = 0;
   if (MPI_Wait(&mpiRequest,MPI_STATUS_IGNORE) != MPI_SUCCESS) success = false;
   if (MPI_Isend(&offsetOut,1,MPI_Type<MPI_Offset>(),target,0,comm,&mpiRequest) != MPI_SUCCESS) success = false;
   
   // Write this process's data:
   if (MPI_File_write_at_all(fileptr,offset,array,arraySize*vectorSize,MPI_Type<T>(),MPI_STATUS_IGNORE) != MPI_SUCCESS) success = false;

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
      xmlWriter->addAttribute(node,"arraysize",totalBytes/sizeof(T)/vectorSize);
      xmlWriter->addAttribute(node,"datatype",arrayDataType<T>());
      xmlWriter->addAttribute(node,"datasize",sizeof(T));
      
      offset = offsetIn;
   }
   return success;
}

#endif
