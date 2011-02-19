#include <cstdlib>
#include <iostream>
#include <vector>
#include <limits>

#include "vlscommon.h"
#include "vlswriter.h"
#include "mpifile.h"

using namespace std;

VlsWriter::VlsWriter() {
   bufferedCells = 0;
   bufferSize = 0;
   bytesPerCellGID = numeric_limits<unsigned char>::max();
   bytesPerSpatNbrListSize = 1;
   N_dimensions = 3;    // This should be given by user
   writeSpatNbrList = false;
}

VlsWriter::~VlsWriter() { }

void VlsWriter::appendHeaderElement(std::vector<unsigned char>& byteArray,const unsigned char& typeID,
				    const unsigned char* const array,const unsigned int& arraySize) {
   byteArray.push_back(arraySize+1);
   byteArray.push_back(typeID);
   for (unsigned int i=0; i<arraySize; ++i) byteArray.push_back(array[i]);
}

/** Close a Vlasov file which has been previously opened by calling 
 * VlsWriter::openRead or VlsWriter::openWrite.
 * @return If true, the file was closed without errors.
 * @see VlsWriter::openRead
 * @see VlsWriter::openWrite
 */
bool VlsWriter::close() {
   byteArray.resize(0);
   bufferedCells = 0;
   bufferSize = 0;
   return mpiFile.close();
}

bool VlsWriter::flushBuffer() {
   // Write buffer to file:
   bool rvalue = true;
   if (byteArray.size() == 0 || bufferedCells == 0) return rvalue;
   
   int myrank;
   MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
   
   if (mpiFile.write(byteArray.size(),&(byteArray[0])) == false) rvalue = false;
   if (mpiFile.getCount<size_t>() != byteArray.size()) rvalue = false;
   // Clear buffer:
   bufferedCells = 0;
   byteArray.clear();
   return rvalue;
}

/** Open a new Vlasov file writing data.
 * @param comm MPI communicator in which the file is opened.
 * @param fname The name of the file.
 * @return If true, the file was opened successfully.
 */
bool VlsWriter::open(MPI_Comm comm,const std::string& fname) {
   int accessMode = (MPI_MODE_WRONLY | MPI_MODE_SEQUENTIAL | MPI_MODE_CREATE);
   bool rvalue = true;
   if (mpiFile.open(comm,fname,MPI_INFO_NULL,accessMode,true) == false) rvalue = false;
   if (mpiFile.resetPosition() == false) rvalue = false;
   return rvalue;
}

bool VlsWriter::reserveSpatCellCoordBuffer(const unsigned int& N_cells,const DataReducer* const dr) {
   bufferedCells = 0;
   bufferSize = N_cells;
   byteArray.clear();
   return true;
}

bool VlsWriter::setBytesPerCellGID(const unsigned char& bytes) {
   if (bytes > sizeof(unsigned long long int)) return false;
   bytesPerCellGID = bytes;
   return true;
}

bool VlsWriter::setWriteSpatNbrsLists(const bool& writeLists) {
   writeSpatNbrList = writeLists;
   return true;
}

bool VlsWriter::sync() {
   return mpiFile.resetPosition();
}

/** Write header into Vlasov output file. After this function 
 * returns, each MPI process in the same communicator as the 
 * process which wrote the header has to set its file view.
 * 
 * The header consists of <size> <type> <element> triplets, where 
 * <size> is a byte containing the sum of the sizes of <type> and 
 * <element> fields in bytes, <type> is a byte that identifies the 
 * element, and <element> is the actual value of the header element.
 * For example, <element> might be the name of a coordinate, and the 
 * coordinate is identified by the value of <type>.
 * 
 * @return If true, the header was written successfully.
 */
bool VlsWriter::writeHeader(MPI_Comm comm,const int& masterRank) {
   bool rvalue = true;
   unsigned char floatEndianness = VlsHeader::LITTLE_END; // Default to little-endian
   unsigned char intEndianness   = VlsHeader::LITTLE_END;

   // Sanity checks:
   if (bytesPerCellGID > sizeof(unsigned long long int)) {
      cerr << "VlsWriter::writeHeader ERROR: bytesPerCellGID value not set!" << endl;
      return false;
   }
   
   // Detect endianness (only important on master process):
   const int tmpValue = 1;
   const unsigned char* const tmpPtr = reinterpret_cast<const unsigned char*>(&tmpValue);
   if (tmpPtr[0] != 1) {
      floatEndianness = VlsHeader::BIG_END;
      intEndianness   = VlsHeader::BIG_END;
   }
   
   // Only master process within the given communicator writes header:
   int myrank;
   MPI_Comm_rank(comm,&myrank);
   if (myrank == masterRank) {
      // Construct a byte array containing the header:
      vector<unsigned char> byteArray;
      unsigned char typeID;

      // Floating point and integer endianness:
      appendHeaderElement(byteArray,VlsHeader::ENDIANNESS_FLOAT,floatEndianness);
      appendHeaderElement(byteArray,VlsHeader::ENDIANNESS_INT  ,intEndianness);
      
      // file version number:
      const string versionNumber = "1.0.0";
      appendHeaderElement(byteArray,VlsHeader::VERSION,versionNumber);
   
      // bytes per global ID (this should be given elsewhere):
      #ifdef PARGRID
         appendHeaderElement(byteArray,VlsHeader::BYTES_PER_CELL_GID,(unsigned char)sizeof(ID::type));
      #else
         appendHeaderElement(byteArray,VlsHeader::BYTES_PER_CELL_GID,(unsigned char)sizeof(uint64_t));
      #endif
   
      // bytes per coordinate value:
      appendHeaderElement(byteArray,VlsHeader::BYTES_PER_CELL_CRD,(unsigned char)sizeof(Real));
   
      // Number of dimensions:
      appendHeaderElement(byteArray,VlsHeader::DIMENSIONS,(unsigned char)3);
      
      // Byte size of variable name size field:
      appendHeaderElement(byteArray,VlsHeader::BYTES_PER_VARNAME_SIZE,DataReducer::getNameSizeEntryByteSize());

      // Does file contain spatial neighbour lists:
      appendHeaderElement(byteArray,VlsHeader::CONTAINS_SPAT_NBR_LIST,(unsigned char)writeSpatNbrList);

      // Byte size of field giving the size of spatial neighbour list entry:
      appendHeaderElement(byteArray,VlsHeader::BYTES_PER_SPAT_NBR_SIZE,bytesPerSpatNbrListSize);
      
      // Append header end marker and write the byte array to disk:
      byteArray.push_back(0);
      if (mpiFile.write(byteArray.size(),&(byteArray[0])) == false) rvalue = false;
   }
   
   // Master process has to let other processes in communicator comm know 
   // if the header was written successfully:
   unsigned char headerRead;
   if (myrank == masterRank) {
      if (rvalue == true) headerRead = 0;
      else headerRead = 1;
   }
   if (MPI_Bcast(&headerRead,1,MPI_UNSIGNED_CHAR,masterRank,comm) != MPI_SUCCESS) rvalue = false;
   if (headerRead > 0) rvalue = false;
   
   // Sync file position pointer on all processes on communicator comm:
   MPI_Barrier(comm);
   if (mpiFile.resetPosition() == false) rvalue = false;
   return rvalue;
}

#ifdef PARGRID
   bool VlsWriter::writeSpatCellCoordEntry(const ID::type& cellID,const SpatialCell& cell,DataReducer* dr) {
#else
   bool VlsWriter::writeSpatCellCoordEntry(const uint64_t& cellID,const SpatialCell& cell,DataReducer* dr) {
#endif
   // Reduce cell data and add to byte array varByteData:
   unsigned char* varByteData = NULL;
   unsigned int varDataByteSize = 0;
   if (dr != NULL) {
      if (dr->reduceData(cell) == false) {
	 cerr << "VlsWriter::writeSpatCellCoordEntry ERROR: DataReducer failed to reduce data!" << endl;
	 return false;
      }
      varDataByteSize = dr->getByteSize();
      varByteData = new unsigned char[varDataByteSize];
      dr->appendReducedData(cell,varByteData);
   }
   
   // Write Version_1_0_0 data:
   const unsigned int SIZE = 6*sizeof(Real) + sizeof(cellID) + varDataByteSize;
   byteArray.reserve(SIZE);
      
   unsigned int index = 0;
   unsigned char* ptr;
   
   // Append cell global ID to byteArray:
   #ifndef PARGRID
      uint tmpID = static_cast<uint>(cellID);
      ptr = reinterpret_cast<unsigned char*>(const_cast<uint*>(&tmpID));
      for (int i=0; i<sizeof(cellID); ++i) byteArray.push_back(ptr[i]);
   #else
      ptr = reinterpret_cast<unsigned char*>(const_cast<uint*>(&cellID));
      for (int i=0; i<sizeof(tmpID); ++i) byteArray.push_back(ptr[i]);
   #endif
   
   // Append cell x,y,z,dx,dy,dz to byteArray:
   ptr = reinterpret_cast<unsigned char*>(cell.cpu_cellParams);
   for (int i=0; i<6*sizeof(Real); ++i) byteArray.push_back(ptr[i]);
      
   // Append static-size variable data from varByteData to byteArray:
   for (unsigned int i=0; i<varDataByteSize; ++i) byteArray.push_back(varByteData[i]);
   delete varByteData;
   varByteData = NULL;
   
   // Write byteArray to file:
   bool rvalue = true; //mpiFile.write(SIZE,byteArray);
   if (mpiFile.write(byteArray.size(),&(byteArray[0])) == false) rvalue = false;
   if (mpiFile.getCount<unsigned int>() != SIZE) rvalue = false;
   return rvalue;
}

#ifdef PARGRID
bool VlsWriter::writeSpatCellCoordEntryBuffered(const ID::type& cellID,const SpatialCell& cell,DataReducer* dr,const vector<ID::type>& nbrs,cuchar& refLevel) {
   typedef ID::type GIDtype;
#else
bool VlsWriter::writeSpatCellCoordEntryBuffered(const uint64_t& cellID,const SpatialCell& cell,DataReducer* dr,const vector<uint64_t>& nbrs,cuchar& refLevel) {
   typedef uint64_t GIDtype;
#endif
   bool rvalue = true;
   // Check that cell global ID is sane:
   if (cellID == numeric_limits<GIDtype>::max()) {
      cerr << "VlsWriter::writeSpatCellCoordEntryBuffered ERROR: invalid cellID!" << endl;
      return false;
   }
   
   // Append cell global ID to buffer:
   const unsigned char* ptr;
   ptr = reinterpret_cast<const unsigned char*>(&cellID);
   for (int i=0; i<sizeof(cellID); ++i) byteArray.push_back(ptr[i]);
   // Append cell x,y,z,dx,dy,dz to buffer:
   ptr = reinterpret_cast<const unsigned char*>(cell.cpu_cellParams);
   for (int i=0; i<6*sizeof(Real); ++i) byteArray.push_back(ptr[i]);
   
   // Reduce cell data and add to buffer:
   if (dr != NULL) {
      if (dr->reduceData(cell) == false) {
	 cerr << "VlsWriter::writeSpatCellCoordEntry ERROR: DataReducer failed to reduce data!" << endl;
	 rvalue = false;
      } else {
	if (dr->appendReducedData(cell,byteArray) == false) {
	   cerr << "VlsWriter::writeSpatCellCoordEntry ERROR: DataReducer failed to append data!" << endl;
	}
      }
   }
   
   if (writeSpatNbrList == true) {
      // Calculate the byte size of neighbour list entry:
      if (bytesPerCellGID*nbrs.size() > numeric_limits<unsigned char>::max()) {
	 cerr << "VlsWriter::writeSpatCellCoordEntry ERROR: spat nbr list size field too small!" << endl;
	 rvalue = false;
      }
      const unsigned char size = bytesPerCellGID*nbrs.size();
      // Append byte size of neighbour list entry:
      ptr = reinterpret_cast<const unsigned char*>(&size);
      for (int i=0; i<sizeof(size); ++i) byteArray.push_back(ptr[i]);
      // Append refinement level:
      byteArray.push_back(refLevel);
      // Append neighbour global IDs:
      ptr = reinterpret_cast<const unsigned char*>(&(nbrs[0]));
      for (unsigned int i=0; i<bytesPerCellGID*nbrs.size(); ++i) byteArray.push_back(ptr[i]);
   }
   
   // If buffer is full, flush it to file:
   ++bufferedCells;
   if (bufferedCells == bufferSize)
     if (flushBuffer() == false) rvalue = false;
   return rvalue;
}

bool VlsWriter::writeSpatCellCoordEntryEndMarker(MPI_Comm comm,const int& masterRank) {
   // Only master process writes end marker:
   int myrank;
   MPI_Comm_rank(comm,&myrank);
   if (myrank != masterRank) return true;
   #ifdef PARGRID
     //cuint END_VALUE = numeric_limits<uint>::max();
     const ID::type END_VALUE = numeric_limits<ID::type>::max();
   #else
     const uint64_t END_VALUE = numeric_limits<uint64_t>::max();
   #endif
   mpiFile << END_VALUE;
   return true;
}

bool VlsWriter::writeStaticVariableDesc(MPI_Comm comm,const int& masterRank,DataReducer* const dr) {
   // Only master process writes variable descriptions:
   int myrank; 
   MPI_Comm_rank(comm,&myrank);

   bool success = true;
   if (myrank == masterRank) {
      // If NULL ptr was given, write empty variable description, i.e. write name size = 0.
      // However, the width of the name field must correspond to the value written to the header.
      if (dr == NULL) {
	 unsigned char zeroArray[] = {0,0,0,0,0,0,0,0};
	 if (DataReducer::getNameSizeEntryByteSize() > 8) success = false;
	 if (success == true) if (mpiFile.write(DataReducer::getNameSizeEntryByteSize(),zeroArray) == false) success = false;
      } else {
	 // Request variable description from DataReducer and write it:
	 unsigned char* tmpArr = NULL;
	 unsigned int tmpArrSize = 0;
	 if (dr->getDescription(tmpArr,tmpArrSize) == false) {
	    cerr << "VlsWriter::writeStaticVariableDesc ERROR: Failed to get variable description from DataReducer!" << endl;
	    unsigned char zeroArray[] = {0,0,0,0,0,0,0,0};
	    if (DataReducer::getNameSizeEntryByteSize() > 8) success = false;
	    if (success == true) if (mpiFile.write(DataReducer::getNameSizeEntryByteSize(),zeroArray) == false) success = false;
	    delete tmpArr;
	    goto exitWrite;
	 }
   
	 if (mpiFile.write(tmpArrSize,tmpArr) == false) success = false;
	 delete tmpArr;
	 tmpArr = NULL;
      }
   }
   
exitWrite:
   MPI_Barrier(comm);
   
   // Let all processes know if master wrote descriptions successfully:
   unsigned char descrWritten;
   if (myrank == masterRank) {
      if (success == true) descrWritten = 0;
      else descrWritten = 1;
   }
   if (MPI_Bcast(&descrWritten,1,MPI_UNSIGNED_CHAR,masterRank,comm) != MPI_SUCCESS) success = false;
   if (descrWritten > 0) success = false;
   
   // Wait that variable descriptions have been written and sync file position pointers:
   if (mpiFile.resetPosition() == false) success = false;
   return success;
}

#ifdef PARGRID
bool VlsWriter::writeVelocityBlockEntryBuffered(const ID::type& cellID,const SpatialCell& cell) {
#else
bool VlsWriter::writeVelocityBlockEntryBuffered(const uint64_t& cellID,const SpatialCell& cell) {
#endif
   bool rvalue = true;
   const char* ptr;
   for (uint b=0; b<cell.N_blocks; ++b) {
      // Write velocity block ID:
      ptr = reinterpret_cast<const char*>(&b);
      for (uint i=0; i<sizeof(b); ++i) byteArray.push_back(ptr[i]);
      // Write velocity block coordinates and size:
      ptr = reinterpret_cast<const char*>(&(cell.cpu_blockParams[b*SIZE_BLOCKPARAMS+BlockParams::VXCRD]));
      for (uint i=0; i<6*sizeof(Real); ++i) byteArray.push_back(ptr[i]);
      // Write velocity block neighbours:
      ptr = reinterpret_cast<const char*>(&(cell.cpu_nbrsVel[b*SIZE_NBRS_VEL+NbrsVel::VXNEG]));
      for (uint i=0; i<6*sizeof(uint); ++i) byteArray.push_back(ptr[i]);
      // Write distribution function data:
      ptr = reinterpret_cast<const char*>(&(cell.cpu_avgs[b*SIZE_VELBLOCK]));
      for (uint i=0; i<WID3*sizeof(Real); ++i) byteArray.push_back(ptr[i]);
   }
   
   ++bufferedCells;
   if (bufferedCells == bufferSize)
     if (flushBuffer() == false) rvalue = false;
   return rvalue;
}
   
#ifdef PARGRID
bool VlsWriter::writeVelocityBlockEntryHeaderBuffered(const ID::type& cellID,cuint& N_blocks) {
#else
bool VlsWriter::writeVelocityBlockEntryHeaderBuffered(const uint64_t& cellID,cuint& N_blocks) {
#endif
   const char* ptr;
   // Write spatial cell global ID:
   ptr = reinterpret_cast<const char*>(&cellID);
   for (uint i=0; i<sizeof(cellID); ++i) byteArray.push_back(ptr[i]);
   // Write number of velocity blocks:
   ptr = reinterpret_cast<const char*>(&N_blocks);
   for (uint i=0; i<sizeof(N_blocks); ++i) byteArray.push_back(ptr[i]);
}


