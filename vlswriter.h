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

#ifdef PARGRID
namespace ID {
   typedef uint type;
}
#endif

class VlsWriter {
 public:
   
   VlsWriter();
   ~VlsWriter();
   
   bool close();
   bool flushBuffer();
   bool open(MPI_Comm comm,const std::string& fname);
   bool reserveSpatCellCoordBuffer(const unsigned int& N_cells,const DataReducer* const dr);
   bool setBytesPerCellGID(const unsigned char& bytes);
   bool setWriteSpatNbrsLists(const bool& writeLists);
   bool sync();
   bool writeHeader(MPI_Comm comm,const int& masterRank);   
   bool writeSpatCellCoordEntryEndMarker(MPI_Comm comm,const int& masterRank);
   bool writeStaticVariableDesc(MPI_Comm comm,const int& masterRank,DataReducer* const dr);

   #ifdef PARGRID
   bool writeSpatCellCoordEntry(const ID::type& cellID,const SpatialCell& cell,DataReducer* const dr);
   bool writeSpatCellCoordEntryBuffered(const ID::type& cellID,const SpatialCell& cell,DataReducer* const dr,const std::vector<ID::type>& nbrs,cuchar& refLevel);
   bool writeVelocityBlockEntryBuffered(const ID::type& cellID,const SpatialCell& cell);
   bool writeVelocityBlockEntryHeaderBuffered(const ID::type& cellID,cuint& N_blocks);
   #else
   bool writeSpatCellCoordEntry(const uint64_t& cellID,const SpatialCell& cell,DataReducer* const dr);
   bool writeSpatCellCoordEntryBuffered(const uint64_t& cellID,const SpatialCell& cell,DataReducer* const dr,const std::vector<uint64_t>& nbrs,cuchar& refLevel);
   bool writeVelocityBlockEntryBuffered(const uint64_t& cellID,const SpatialCell& cell);
   bool writeVelocityBlockEntryHeaderBuffered(const uint64_t& cellID,cuint& N_blocks);
   #endif
   
 private:
   unsigned int bufferedCells;           /**< How many cells are currently in buffer.*/
   unsigned int bufferSize;              /**< Size of buffer, measured in cells.*/
   unsigned char bytesPerCellGID;        /**< The size of spatial cell global ID field, in bytes. This value is written to header.*/
   unsigned char bytesPerSpatNbrListSize; /**< The byte size of a field giving the byte size spatial cell neighbour list entry.*/
   std::vector<unsigned char> byteArray;
   MPIFile mpiFile;                     /**< MPIFile which is used for I/O.*/
   unsigned char N_dimensions;          /**< Dimensionality (spatial) of the data. Allowed values are 1, 2, or 3.*/
   bool writeSpatNbrList;               /**< If true, spatial neighbour lists are written. This value is written to header.*/
   
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
