#ifndef VLSVREADER2_H
#define VLSVREADER2_H

#include <stdint.h>
#include <list>
#include <fstream>
#include "muxml.h"
#include "vlscommon.h"

class VLSVReader {
 public:
   VLSVReader();
   ~VLSVReader();
   
   bool close();
   bool getArrayInfo(const std::string& tagName,const std::string& arrayName,uint64_t& arraySize,uint64_t& vectorSize,VLSV::datatype& dataType,uint64_t& byteSize) const;
   bool getArrayInfo(const std::string& tagName,const std::string& arrayName,const std::string& meshName,uint64_t& arraySize,uint64_t& vectorSize,VLSV::datatype& dataType,uint64_t& byteSize) const;
   bool getBlockVariableNames(const std::string& meshName,std::list<std::string>& varNames) const;
   bool getMeshNames(std::list<std::string>& meshNames) const;
   bool getArrayInfo(const std::string& tagName,const std::list<std::pair<std::string,std::string> >& attribs,uint64_t& arraySize,uint64_t& vectorSize,VLSV::datatype& dataType,uint64_t& byteSize) const;
   bool getVariableNames(const std::string& meshName,std::list<std::string>& varNames) const;
   bool open(const std::string& fname);
   bool readArray(const std::string& tagName,const std::string& arrayName,const uint64_t& begin,const uint64_t& amount,char* buffer);
   bool readArray(const std::string& tagName,const std::string& arrayName,const std::list<std::pair<std::string,std::string> >& attribs,const uint64_t& begin,const uint64_t& end,char* buffer);
   
 private:
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

#endif
