#ifndef VLSWRITER_H
#define VLSWRITER_H

#include <map>
#include <mpi.h>
#include <vector>
#include <stdint.h>
#include <cmath>

#include "definitions.h"
#include "mpifile.h"
#include "cell_spatial.h"
#include "datareducer.h"

class VlsWriter {
 public:
   
   VlsWriter();
   ~VlsWriter();
   
   bool close();
   bool flushBuffer();
   bool open(MPI_Comm comm,const std::string& fname);
   bool reserveSpatCellCoordBuffer(const unsigned int& N_cells,const DataReducer* const dr);
   bool sync();
   bool writeHeader(MPI_Comm comm,const int& masterRank);   
   bool writeSpatCellCoordEntry(cuint& cellID,const SpatialCell& cell,DataReducer* const dr);
   bool writeSpatCellCoordEntryBuffered(cuint& cellID,const SpatialCell& cell,DataReducer* const dr);
   bool writeSpatCellCoordEntryEndMarker(MPI_Comm comm,const int& masterRank);
   bool writeStaticVariableDesc(MPI_Comm comm,const int& masterRank,DataReducer* const dr);

 private:
   size_t bufferPointer;                /**< Next free position in VlsWriter::byteArray where data can 
					 * be written.*/
   size_t bufferSize;                   /**< The size of VlsWriter::byteArray.*/
   unsigned char* byteArray;            /**< Byte array which is used for buffered writing to a vlsv file.*/
   MPIFile mpiFile;                     /**< MPIFile which is used for I/O.*/
   unsigned int sizeByteArray;
   
   template<typename T> void appendHeaderElement(std::vector<unsigned char>& byteArray,const unsigned char& typeID,const T& element);
   void appendHeaderElement(std::vector<unsigned char>& byteArray,const unsigned char& typeID,const unsigned char* const array,const unsigned int& arraySize);
};
   
// **********************************
// ****** TEMPLATE DEFINITIONS ******
// **********************************

/** A specialized template of appendHeaderElement for strings.
 * @param byteArray Byte array in which the given header element is inserted.
 * @param typeID The ID of the header element.
 * @param element The value of the header element.
 * @see VlsWriter::appendHeaderElement.
 */
template<>
inline void VlsWriter::appendHeaderElement(std::vector<unsigned char>& byteArray,const unsigned char& typeID,const std::string& element) {
   const unsigned char byteID = static_cast<unsigned char>(typeID);
   byteArray.push_back(element.size()+1);
   byteArray.push_back(byteID);
   unsigned char* ptr = reinterpret_cast<unsigned char*>(const_cast<char*>(&(element[0])));
   for (size_t i=0; i<element.size(); ++i) byteArray.push_back(ptr[i]);
   byteArray.push_back(static_cast<unsigned char>('\0'));
}

/** Append the given header element into the given byte array. What is actually inserted into 
 * byte array is the <byte length>-<header ID>-<value> triplet.
 * @param byteArray Byte array in which the given header element is inserted.
 * @param typeID The ID of the header element.
 * @param element The value of the header element.
 */
template<typename T> 
inline void VlsWriter::appendHeaderElement(std::vector<unsigned char>& byteArray,const unsigned char& typeID,const T& element) {
   const unsigned char byteID = static_cast<unsigned char>(typeID);
   byteArray.push_back(sizeof(element));
   byteArray.push_back(byteID);
   unsigned char* ptr = reinterpret_cast<unsigned char*>(const_cast<T*>(&element));
   for (size_t i=0; i<sizeof(element); ++i) byteArray.push_back(ptr[i]);
}

#endif
