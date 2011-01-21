#include <cstdlib>
#include <iostream>
#include <vector>
#include <limits>

#include "vlscommon.h"
#include "vlswriter.h"
#include "mpifile.h"

using namespace std;

VlsWriter::VlsWriter() {
   bufferPointer = 0;
   bufferSize = 0;
   byteArray = NULL;
   sizeByteArray = 0;
}

VlsWriter::~VlsWriter() {
   delete byteArray;
   byteArray = NULL;
}

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
bool VlsWriter::close() {return mpiFile.close();}

bool VlsWriter::flushBuffer() {
   // Write buffer to file:
   bool rvalue = true;
   if (bufferPointer == 0 || bufferSize == 0) return rvalue;
   if (mpiFile.write(bufferPointer,byteArray) == false) rvalue = false;
   if (mpiFile.getCount<unsigned int>() != bufferPointer) rvalue = false;
   bufferPointer = 0;
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
   delete byteArray;
   
   // Calculate byte size of one coordinate entry, and reserve large enough array:
   size_t entrySize = sizeof(uint) + 6*sizeof(Real);
   if (dr != NULL) entrySize += dr->getByteSize();   
   bufferSize = N_cells * entrySize;
   bufferPointer = 0;
   
   byteArray = new unsigned char[bufferSize];
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
   
      // bytes per global ID:
      appendHeaderElement(byteArray,VlsHeader::BYTES_PER_CELL_GID,(unsigned char)sizeof(uint));
   
      // bytes per coordinate value:
      appendHeaderElement(byteArray,VlsHeader::BYTES_PER_CELL_CRD,(unsigned char)sizeof(Real));
   
      // Number of dimensions:
      appendHeaderElement(byteArray,VlsHeader::DIMENSIONS,(unsigned char)3);
      
      // Byte size of variable name size field:
      appendHeaderElement(byteArray,VlsHeader::BYTES_PER_VARNAME_SIZE,DataReducer::getNameSizeEntryByteSize());

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

bool VlsWriter::writeSpatCellCoordEntry(cuint& cellID,const SpatialCell& cell,DataReducer* dr) {
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
   const unsigned int SIZE = 6*sizeof(Real) + sizeof(uint) + varDataByteSize;
   char byteArray[SIZE];

   unsigned int index = 0;
   unsigned char* ptr;
   
   // Append cell global ID to byteArray:
   ptr = reinterpret_cast<unsigned char*>(const_cast<uint*>(&cellID));
   for (int i=0; i<sizeof(cellID); ++i) {byteArray[index] = ptr[i]; ++index;}
   
   // Append cell x,y,z,dx,dy,dz to byteArray:
   ptr = reinterpret_cast<unsigned char*>(cell.cpu_cellParams);
   for (int i=0; i<6*sizeof(Real); ++i) {byteArray[index] = ptr[i]; ++index;}   
   
   // Append static-size variable data from varByteData to byteArray:
   for (unsigned int i=0; i<varDataByteSize; ++i) {byteArray[index] = varByteData[i]; ++index;}
   delete varByteData;
   varByteData = NULL;
   
   // Write byteArray to file:
   bool rvalue = true; //mpiFile.write(SIZE,byteArray);
   if (mpiFile.write(SIZE,byteArray) == false) rvalue = false;
   if (mpiFile.getCount<unsigned int>() != SIZE) rvalue = false;
   return rvalue;
}

bool VlsWriter::writeSpatCellCoordEntryBuffered(cuint& cellID,const SpatialCell& cell,DataReducer* dr) {
   bool rvalue = true;
   // Check that buffer has space for the entry. If it doesn't, 
   // flush buffer to file before appending the given data:
   size_t entrySize = sizeof(uint) + 6*sizeof(Real);
   if (dr != NULL) entrySize += dr->getByteSize();
   if (bufferPointer + entrySize > bufferSize) {
      if (flushBuffer() == false) rvalue = false;
   }
   
   // Append cell global ID to buffer:
   const unsigned char* ptr;
   ptr = reinterpret_cast<const unsigned char*>(&cellID);
   for (int i=0; i<sizeof(cellID); ++i) {byteArray[bufferPointer] = ptr[i]; ++bufferPointer;}
   // Append cell x,y,z,dx,dy,dz to buffer:
   ptr = reinterpret_cast<const unsigned char*>(cell.cpu_cellParams);
   for (int i=0; i<6*sizeof(Real); ++i) {byteArray[bufferPointer] = ptr[i]; ++bufferPointer;}
   
   // Reduce cell data and add to buffer:
   if (dr != NULL) {
      if (dr->reduceData(cell) == false) {
	 cerr << "VlsWriter::writeSpatCellCoordEntry ERROR: DataReducer failed to reduce data!" << endl;
	 rvalue = false;
      } else {
	 dr->appendReducedData(cell,byteArray + bufferPointer);
	 bufferPointer += dr->getByteSize();
      }
   }
   return rvalue;
}

bool VlsWriter::writeSpatCellCoordEntryEndMarker(MPI_Comm comm,const int& masterRank) {
   // Only master process writes end marker:
   int myrank;
   MPI_Comm_rank(comm,&myrank);
   if (myrank != masterRank) return true;
   cerr << "Writing end marker" << endl;
   cuint END_VALUE = numeric_limits<uint>::max();
   cerr << "\t End marker value = '" << (cuint)END_VALUE << "'" << endl;
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




