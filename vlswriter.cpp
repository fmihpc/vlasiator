#include <cstdlib>
#include <iostream>
#include <vector>

#include "vlswriter.h"
#include "mpifile.h"

using namespace std;

MPI_Datatype MPIUINT64 = MPI_UNSIGNED_LONG_LONG;

VlsHeader::Real (*ptrCellCRD)(const unsigned char* const);
VlsHeader::Int (*ptrCellGID)(const unsigned char* const);
VlsHeader::Int (*ptrVarNameSize)(const unsigned char* const);

VlsHeader::Int convUInt1(const unsigned char* const ptr) {return *(reinterpret_cast<uint8_t*>(const_cast<unsigned char* const>(ptr)));}
VlsHeader::Int convUInt2(const unsigned char* const ptr) {return *(reinterpret_cast<uint16_t*>(const_cast<unsigned char* const>(ptr)));}
VlsHeader::Int convUInt4(const unsigned char* const ptr) {return *(reinterpret_cast<uint32_t*>(const_cast<unsigned char* const>(ptr)));}
VlsHeader::Int convUInt8(const unsigned char* const ptr) {return *(reinterpret_cast<uint64_t*>(const_cast<unsigned char* const>(ptr)));}

VlsHeader::Real convReal4(const unsigned char* const ptr) {return *(reinterpret_cast<float*>(const_cast<unsigned char* const>(ptr)));}
VlsHeader::Real convReal8(const unsigned char* const ptr) {return *(reinterpret_cast<double*>(const_cast<unsigned char* const>(ptr)));}
VlsHeader::Real convReal16(const unsigned char* const ptr) {return *(reinterpret_cast<long double*>(const_cast<unsigned char* const>(ptr)));}

enum VersionNumber {
   Unknown,
   Version_1_0_0
};

VersionNumber versionNumber = Unknown;

VersionNumber getVersionNumber(const std::string& stringVersion) {
   if (stringVersion == "1.0.0") return Version_1_0_0;
   return Unknown;
}

VlsWriter::VarDesc::VarDesc(const std::string& name,const unsigned char& typeID,const unsigned char& elementBytes): 
name(name),typeID(typeID),elementBytes(elementBytes) {
   byteSize = elementBytes;
   if (typeID == VlsHeader::NULLVARIABLE)  {byteSize *= 0; elements = 0;}
   else if (typeID == VlsHeader::SCALAR)   {byteSize *= 1; elements = 1;}
   else if (typeID == VlsHeader::VECTOR2)  {byteSize *= 2; elements = 2;}
   else if (typeID == VlsHeader::VECTOR3)  {byteSize *= 3; elements = 3;}
   else if (typeID == VlsHeader::TENSOR22) {byteSize *= 4; elements = 4;}
   else if (typeID == VlsHeader::TENSOR23) {byteSize *= 6; elements = 6;}
   else if (typeID == VlsHeader::TENSOR32) {byteSize *= 6; elements = 6;}
   else if (typeID == VlsHeader::TENSOR33) {byteSize *= 9; elements = 9;}
}

VlsWriter::VlsWriter() {
   byteArray = NULL;
   N_dimensions = 0;
   sizeByteArray = 0;
   sizeCellCoordEntry = 0;
   sizeCellCRD = 0;
   sizeCellGID = 0;
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

bool VlsWriter::broadcastVariableValues(MPI_Comm comm,const int& masterRank) {
   // Might be more efficient to pack the values of all these variables into a single 
   // MPI_Datatype, but for now this is easier to use:
   bool success = true;
   if (MPI_Bcast(&sizeCellCRD         ,1,MPI_UNSIGNED_CHAR,masterRank,comm) != MPI_SUCCESS) success = false;
   if (MPI_Bcast(&sizeCellGID         ,1,MPI_UNSIGNED_CHAR,masterRank,comm) != MPI_SUCCESS) success = false;
   if (MPI_Bcast(&N_dimensions        ,1,MPI_UNSIGNED_CHAR,masterRank,comm) != MPI_SUCCESS) success = false;
   if (MPI_Bcast(&varNameFieldByteSize,1,MPI_UNSIGNED_CHAR,masterRank,comm) != MPI_SUCCESS) success = false;

   // Broadcast version number:
   int vnumber = (int)versionNumber;
   if (MPI_Bcast(&vnumber,1,MPI_INT,masterRank,comm) != MPI_SUCCESS) success = false;
   versionNumber = static_cast<VersionNumber>(vnumber);
   return success;
}

void VlsWriter::calculateDerivedVariables() {
   // Calculate the size of spatial cell coordinate entry:
   sizeCellCoordEntry = sizeCellGID + N_dimensions*2*sizeCellCRD;
   
   for (vector<VarDesc>::const_iterator it=staticSizeVars.begin(); it!=staticSizeVars.end(); ++it) 
     sizeCellCoordEntry += it->byteSize;
}

bool VlsWriter::calculateVariableSizes() {
   const unsigned char* ptr;
   map<unsigned char,vector<unsigned char> >::const_iterator it;
   // Version number:
   it = headerElements.find(VlsHeader::VERSION);
   if (it == headerElements.end()) return false;
   ptr = &(it->second[0]);
   versionNumber = getVersionNumber(reinterpret_cast<const char*>(ptr));
   
   // Size of spatial cell coordinate entry:
   it = headerElements.find(VlsHeader::BYTES_PER_CELL_CRD);
   if (it == headerElements.end()) return false;
   sizeCellCRD = it->second[0];
   
   // Size of spatial cell global ID entry:
   it = headerElements.find(VlsHeader::BYTES_PER_CELL_GID);
   if (it == headerElements.end()) return false;
   sizeCellGID = it->second[0];
   
   // Number of spatial dimensions:
   it = headerElements.find(VlsHeader::DIMENSIONS);
   if (it == headerElements.end()) return false;
   N_dimensions = it->second[0];
   
   // Size of entry containing the size of variable name:
   it = headerElements.find(VlsHeader::BYTES_PER_VARNAME_SIZE);
   if (it == headerElements.end()) return false;
   varNameFieldByteSize = it->second[0];
   
   headerElements.clear();
}

/** Close a Vlasov file which has been previously opened by calling 
 * VlsWriter::openRead or VlsWriter::openWrite.
 * @return If true, the file was closed without errors.
 * @see VlsWriter::openRead
 * @see VlsWriter::openWrite
 */
bool VlsWriter::close() {return mpiFile.close();}

size_t VlsWriter::getNumberOfStaticVars() const {return staticSizeVars.size();}

unsigned int VlsWriter::getStaticVarElements(const size_t& varID) const {
   if (varID >= staticSizeVars.size()) return 0;
   return staticSizeVars[varID].elements;
}

std::string VlsWriter::getStaticVarName(const size_t& varID) const {
   if (varID >= staticSizeVars.size()) return string("");
   return staticSizeVars[varID].name;
}

unsigned char VlsWriter::getStaticVarType(const size_t& varID) const {
   if (varID >= staticSizeVars.size()) return VlsHeader::NULLVARIABLE;
   return staticSizeVars[varID].typeID;
}

/** Open a Vlasov file for reading data.
 * @param comm MPI communicator in which the file is opened.
 * @param fname The name of the file.
 * @return If true, the file was opened successfully.
 */
bool VlsWriter::openRead(MPI_Comm comm,const std::string& fname) {
   int accessMode = (MPI_MODE_RDONLY | MPI_MODE_SEQUENTIAL);
   bool rvalue = true;
   if (mpiFile.open(comm,fname,MPI_INFO_NULL,accessMode,false) == false) rvalue = false;
   if (mpiFile.resetPosition() == false) rvalue = false;
   return rvalue;
}

/** Open a new Vlasov file writing data.
 * @param comm MPI communicator in which the file is opened.
 * @param fname The name of the file.
 * @return If true, the file was opened successfully.
 */
bool VlsWriter::openWrite(MPI_Comm comm,const std::string& fname) {
   int accessMode = (MPI_MODE_WRONLY | MPI_MODE_SEQUENTIAL | MPI_MODE_CREATE);
   bool rvalue = true;
   if (mpiFile.open(comm,fname,MPI_INFO_NULL,accessMode,true) == false) rvalue = false;
   if (mpiFile.resetPosition() == false) rvalue = false;
   return rvalue;
}

/** Read the header of Vlasov file. After calling this function the values 
 * of header elements can be queried with VlsWriter::getHeaderElement member 
 * function.
 * @param comm The MPI communicator in which the file header is read. Note 
 * that only the master process within this communicator reads the header, 
 * all other processes return immediately with value false.
 * @param masterRank The rank of the master process in communicator comm.
 * @return If true, this process read the header successfully.
 * @see VlsWriter::getHeaderElement
 */
bool VlsWriter::readHeader(MPI_Comm comm,const int& masterRank) {
   int N_processes;
   MPI_Comm_size(comm,&N_processes);
   
   bool rvalue = true;
   headerElements.clear();
   // Only master process within the given communicator reads header:
   int myrank;
   MPI_Comm_rank(comm,&myrank);

   if (myrank == masterRank) {
      // Read <byte length>-<ID byte>-<value> triplets until length=0 is encountered:
      unsigned char length;
      unsigned char byteArray[255];
      unsigned char byteID;
      mpiFile >> length;
      while (length > 0 && rvalue == true) {
	 mpiFile >> byteID;
	 if (mpiFile.read(length,byteArray) == false) rvalue = false;
	 
	 // Insert the read element into headerElements:
	 vector<unsigned char> tmp(length);
	 for (unsigned char i=0; i<length; ++i) tmp[i] = byteArray[i];
	 headerElements[byteID] = tmp;
	 
	 mpiFile >> length;
      }
      
      // Determine values for internal variables:
      if (rvalue == true) calculateVariableSizes();
   }
   MPI_Barrier(comm);
   
   // Let all processes know whether or not the master process read header correctly:
   unsigned char headerRead;
   if (myrank == masterRank) {
      if (rvalue == true) headerRead = 0;
      else headerRead = 1;
   }
   if (MPI_Bcast(&headerRead,1,MPI_UNSIGNED_CHAR,masterRank,comm) != MPI_SUCCESS) rvalue = false;
   if (headerRead > 0) rvalue = false;
   
   // Since only the master process parses the header, the values of 
   // internal variables need to be broadcasted to other processes:
   if (rvalue == true) if (broadcastVariableValues(comm,masterRank) == false) rvalue = false;
   
   // Set int and float converters:
   if (rvalue == true) if (setInternalPointers() == false) rvalue = false;
   
   // Read static-size variable descriptions:
   if (rvalue == true) if (readStaticVariableDesc(comm,masterRank) == false) rvalue = false;

   // Calculate the values of internal variables that are derived 
   // from broadcasted values:
   if (rvalue == true) calculateDerivedVariables();
   
   // Sync file position pointer on all processes on communicator comm:
   mpiFile.resetPosition();
   return rvalue;
}

bool VlsWriter::readSpatCellCoordEntry() {
   // Make sure byteArray is large enough:
   if (sizeByteArray < sizeCellCoordEntry) {
      delete byteArray;
      byteArray = new unsigned char[sizeCellCoordEntry];
      sizeByteArray = sizeCellCoordEntry;
   }
   
   // Attempt to read spatial cell coordinate entry:   
   if (mpiFile.read(sizeCellCoordEntry,byteArray) == false) return false;
   if (mpiFile.getCount<unsigned int>() < sizeCellCoordEntry) return false; // End of file
   
   // Store values to internal variables. How this is done depends on 
   // the file version:
   if (versionNumber > Unknown) {
      spatCellGID = ptrCellGID(byteArray);
      xcrd        = ptrCellCRD(byteArray + sizeCellGID);
      ycrd        = ptrCellCRD(byteArray + sizeCellGID + 1*sizeCellCRD);
      zcrd        = ptrCellCRD(byteArray + sizeCellGID + 2*sizeCellCRD);
      dx          = ptrCellCRD(byteArray + sizeCellGID + 3*sizeCellCRD);
      dy          = ptrCellCRD(byteArray + sizeCellGID + 4*sizeCellCRD);
      dz          = ptrCellCRD(byteArray + sizeCellGID + 5*sizeCellCRD);
   }
   return true;
}

bool VlsWriter::readStaticVariableDesc(MPI_Comm comm,const int& masterRank) {
   staticSizeVars.clear();
   if (varNameFieldByteSize == 0) return false;
   
   // Only master process reads variable descriptions:
   bool success = true;
   int myrank;
   MPI_Comm_rank(comm,&myrank);
   
   // Reserve large enough byte array:
   delete byteArray;
   byteArray = new unsigned char[65536];
   sizeByteArray = 65536;
   
   // Read variable descriptions while the size of variable name > 0, 
   // and broadcast the descriptions to all processes:
   do {
      if (myrank == masterRank) {
	 // Read the byte size of variable name:
	 if (mpiFile.read(varNameFieldByteSize,byteArray) == false) {success = false; varNameSize = 0;}
	 if (mpiFile.getCount<unsigned int>() != varNameFieldByteSize) {success = false; varNameSize = 0;}
	 varNameSize = ptrVarNameSize(byteArray);
	 
	 // If variable exists, read its description:
	 if (varNameSize > 0) {
	    unsigned long int bytes = varNameSize + 2; // size(name) + typeID + elementBytes
	    if (mpiFile.read(bytes,byteArray) == false) {success = false; varNameSize = 0;}
	    if (mpiFile.getCount<unsigned long int>() != bytes) {success = false; varNameSize = 0;}
	 }
      }
      
      // Broadcast variable name size and variable byte description, if size > 0:
      MPI_Bcast(&varNameSize,1,MPIUINT64,masterRank,comm);
      if (varNameSize > 0) {
	 if (MPI_Bcast(byteArray,varNameSize+2,MPI_UNSIGNED_CHAR,masterRank,comm) != MPI_SUCCESS) {
	    success = false;
	    cerr << "VlsWriter::readStaticVariableDesc ERROR: Failed to Bcast variable description!" << endl;
	 }
	 // Parse description:
	 string varName;
	 unsigned char varTypeID;
	 unsigned char varElementSize;
	 varName = reinterpret_cast<char*>(byteArray);
	 varTypeID = byteArray[varName.size() + 1];      // +1 from end-of-character
	 varElementSize = byteArray[varName.size() + 2];
	 
	 staticSizeVars.push_back(VarDesc(varName,varTypeID,varElementSize));
      }
   } while (varNameSize > 0);
   
   // Calculate static size variable byte offsets:
   if (staticSizeVars.size() > 0) {
      staticSizeVars[0].offset = sizeCellGID + 2*N_dimensions*sizeCellCRD;
      for (size_t i=1; i<staticSizeVars.size(); ++i) {
	staticSizeVars[i].offset = staticSizeVars[i-1].offset + staticSizeVars[i-1].byteSize;
      }
   }
   
   // Let all processes know if the master process read variable 
   // descriptions successfully:
   unsigned char variablesRead;
   if (myrank == masterRank) {
      if (success == true) variablesRead = 0;
      else variablesRead = 1;
   }
   if (MPI_Bcast(&variablesRead,1,MPI_UNSIGNED_CHAR,masterRank,comm) != MPI_SUCCESS) success = false;
   if (variablesRead > 0) success = false;
   
   delete byteArray;
   byteArray = NULL;
   sizeByteArray = 0;
   
   return success;
}

bool VlsWriter::setInternalPointers() {
   bool rvalue = true;
   // Set converter for spatial cell coordinate value:
   if (sizeCellCRD == 4) ptrCellCRD = &convReal4;
   else if (sizeCellCRD == 8) ptrCellCRD = &convReal8;
   else if (sizeCellCRD == 16) ptrCellCRD = &convReal16;
   else {ptrCellCRD = NULL; rvalue = false;}
   
   // Set converter for spatial cell global ID value:
   if (sizeCellGID == 1) ptrCellGID = &convUInt1;
   else if (sizeCellGID == 2) ptrCellGID = &convUInt2;
   else if (sizeCellGID == 4) ptrCellGID = &convUInt4;
   else if (sizeCellGID == 8) ptrCellGID = &convUInt8;
   else {ptrCellGID = NULL; rvalue = false;}

   // Set converter for variable name size field:
   if (varNameFieldByteSize == 1) ptrVarNameSize = &convUInt1;
   else if (varNameFieldByteSize == 2) ptrVarNameSize = &convUInt2;
   else if (varNameFieldByteSize == 4) ptrVarNameSize = &convUInt4;
   else if (varNameFieldByteSize == 8) ptrVarNameSize = &convUInt8;
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
   
   // Only master process within the given communicator writes header:
   int myrank;
   MPI_Comm_rank(comm,&myrank);
   if (myrank == masterRank) {
      // Construct a byte array containing the header:
      vector<unsigned char> byteArray;
      unsigned char typeID;
   
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
      appendHeaderElement(byteArray,VlsHeader::BYTES_PER_VARNAME_SIZE,DataReducer::getNameByteSize());

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
	 if (DataReducer::getNameByteSize() > 8) success = false;
	 if (success == true) if (mpiFile.write(DataReducer::getNameByteSize(),zeroArray) == false) success = false;
      } else {
	 // Request variable description from DataReducer and write it:
	 unsigned char* tmpArr = NULL;
	 unsigned int tmpArrSize = 0;
	 if (dr->getDescription(tmpArr,tmpArrSize) == false) {
	    cerr << "VlsWriter::writeStaticVariableDesc ERROR: Failed to get variable description from DataReducer!" << endl;
	    unsigned char zeroArray[] = {0,0,0,0,0,0,0,0};
	    if (DataReducer::getNameByteSize() > 8) success = false;
	    if (success == true) if (mpiFile.write(DataReducer::getNameByteSize(),zeroArray) == false) success = false;
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




