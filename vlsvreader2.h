#ifndef VLSVREADER2_H
#define VLSVREADER2_H

#include <mpi.h>
#include <stdint.h>
#include <list>
#include <fstream>
#include "muxml.h"
#include "vlscommon.h"

class VLSVReader {
 public:
   VLSVReader();
   virtual ~VLSVReader();
   
   virtual bool close();
   virtual bool getArrayInfo(const std::string& tagName,const std::string& arrayName,uint64_t& arraySize,
			     uint64_t& vectorSize,VLSV::datatype& dataType,uint64_t& byteSize) const;
   virtual bool getArrayInfo(const std::string& tagName,const std::string& arrayName,const std::string& meshName,
			     uint64_t& arraySize,uint64_t& vectorSize,VLSV::datatype& dataType,uint64_t& byteSize) const;
   virtual bool getBlockVariableNames(const std::string& meshName,std::list<std::string>& varNames) const;
   virtual bool getMeshNames(std::list<std::string>& meshNames) const;
   virtual bool getArrayInfo(const std::string& tagName,const std::list<std::pair<std::string,std::string> >& attribs,
			     uint64_t& arraySize,uint64_t& vectorSize,VLSV::datatype& dataType,uint64_t& byteSize) const;
   virtual bool getVariableNames(const std::string& meshName,std::list<std::string>& varNames) const;
   virtual bool loadArray(const std::string& tagName,const std::list<std::pair<std::string,std::string> >& attribs);
   virtual bool open(const std::string& fname);
   virtual bool readArray(const std::string& tagName,const std::string& arrayName,const uint64_t& begin,const uint64_t& amount,char* buffer);
   virtual bool readArray(const std::string& tagName,const std::list<std::pair<std::string,std::string> >& attribs,
			  const uint64_t& begin,const uint64_t& end,char* buffer);
   virtual bool readArray(const std::string& tagName,const std::string& arrayName,const std::list<std::pair<std::string,std::string> >& attribs,
			  const uint64_t& begin,const uint64_t& end,char* buffer);
   
 protected:
   unsigned char endiannessFile;   /**< Endianness in VLSV file.*/
   unsigned char endiannessReader; /**< Endianness of computer which reads the data.*/
   std::fstream filein;            /**< Input file stream.*/
   std::string fileName;           /**< Name of the input file.*/
   bool fileOpen;                  /**< If true, a file is currently open.*/
   bool swapIntEndianness;         /**< If true, endianness should be swapped on read data (not implemented yet).*/
   MuXML xmlReader;                /**< XML reader used to parse VLSV footer.*/
   
   /** Struct used to store information on the currently open array.*/
   struct ArrayOpen {
      std::streamoff offset;
      std::string tagName;
      std::string arrayName;
      VLSV::datatype dataType;
      uint64_t arraySize;
      uint64_t vectorSize;
      uint64_t dataSize;
   } arrayOpen;
};

class VLSVParReader: public VLSVReader {
 public:
   VLSVParReader();
   ~VLSVParReader();
   
   bool close();
   bool getArrayInfo(const std::string& tagName,const std::list<std::pair<std::string,std::string> >& attribs,
		     uint64_t& arraySize,uint64_t& vectorSize,VLSV::datatype& dataType,uint64_t& byteSize);
   bool getArrayInfoMaster(const std::string& tagName,const std::list<std::pair<std::string,std::string> >& attribs,
			   uint64_t& arraySize,uint64_t& vectorSize,VLSV::datatype& dataType,uint64_t& dataSize);
   bool multiReadStart(const std::string& tagName,const std::list<std::pair<std::string,std::string> >& attribs);
   bool multiReadAddUnit(const uint64_t& amount,char* buffer);
   bool multiReadEnd(const uint64_t& offset);
   bool open(const std::string& fname,MPI_Comm comm,const int& masterRank,MPI_Info mpiInfo);
   bool readArrayMaster(const std::string& tagName,const std::list<std::pair<std::string,std::string> >& attribs,
			const uint64_t& begin,const uint64_t& amount,char* buffer);
   bool readArray(const std::string& tagName,const std::list<std::pair<std::string,std::string> >& attribs,
		  const uint64_t& begin,const uint64_t& amount,char* buffer);
   
 private:
   MPI_Comm comm;                  /**< MPI communicator used to read the file.*/
   MPI_File filePtr;               /**< MPI file pointer to input file.*/
   int masterRank;                 /**< MPI rank of master process.*/
   int myRank;                     /**< MPI rank of this process in communicator comm.*/
   int processes;                  /**< Number of MPI processes in communicator comm.*/
   
   bool getArrayInfo(const std::string& tagName,const std::list<std::pair<std::string,std::string> >& attribs);
   
   MPI_Datatype multiReadVectorType;
   MPI_Datatype multiReadArrayType;
   std::map<char*,uint64_t> multiReadUnits; // buffer,begin,amount
};

#endif
