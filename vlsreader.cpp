#include <cstdlib>
#include <iostream>
#include <limits>
#include "vlsreader.h"

using namespace std;

typedef float Real4;
typedef double Real8;
typedef long double Real16;

static bool fileFloatsLittleEndian = true;
static bool fileIntsLittleEndian = true;
static bool readerFloatsLittleEndian;
static bool readerIntsLittleEndian;

static unsigned int sizeBaNbrList = 0;
static unsigned int sizeBaDynData = 0;
static unsigned char* baNbrList = NULL; // Byte array for reading neighbour lists
static unsigned char* baDynData = NULL; // Byte array for dynamic-size data

static unsigned int xNegNbrs = 0;
static unsigned int xPosNbrs = 0;
static unsigned int yNegNbrs = 0;
static unsigned int yPosNbrs = 0;
static unsigned int zNegNbrs = 0;
static unsigned int zPosNbrs = 0;

static VlsHeader::UInt nbrIDArray[24];

unsigned int countNbrListSize(const char& nbrStatus) {
   char byte = 1;

   xNegNbrs = 1;
   xPosNbrs = 1;
   yNegNbrs = 1;
   yPosNbrs = 1;
   zNegNbrs = 1;
   zPosNbrs = 1;
   if ((nbrStatus & byte) == 1) xNegNbrs = 4;
   byte = byte << 1;
   if ((nbrStatus & byte) == 1) xPosNbrs = 4;
   byte = byte << 1;
   if ((nbrStatus & byte) == 1) yNegNbrs = 4;
   byte = byte << 1;
   if ((nbrStatus & byte) == 1) yPosNbrs = 4;
   byte = byte << 1;
   if ((nbrStatus & byte) == 1) zNegNbrs = 4;
   byte = byte << 1;
   if ((nbrStatus & byte) == 1) zPosNbrs = 4;
   byte = byte << 1;
   
   return xNegNbrs + xPosNbrs + yNegNbrs + yPosNbrs + zNegNbrs + zPosNbrs;
}

void VlsReader::parseNeighbourIndices() {
   unsigned int offset = 1;
   for (int i=0; i<24; ++i) nbrIDArray[i] = spatCellGID;
   // Parse -x neighbour IDs:
   nbrIDArray[ 0] = getUInt(baNbrList+offset,sizeCellGID,swapIntEndian);
   if (xNegNbrs == 4) {
      nbrIDArray[ 1] = getUInt(baNbrList+offset+1*sizeCellGID,sizeCellGID,swapIntEndian);
      nbrIDArray[ 2] = getUInt(baNbrList+offset+2*sizeCellGID,sizeCellGID,swapIntEndian);
      nbrIDArray[ 3] = getUInt(baNbrList+offset+3*sizeCellGID,sizeCellGID,swapIntEndian);
      offset += 4*sizeCellGID;
   } else
     offset += sizeCellGID;
   // Parse +x neighbour IDs:
   nbrIDArray[ 4] = getUInt(baNbrList+offset,sizeCellGID,swapIntEndian);
   if (xPosNbrs == 4) {
      nbrIDArray[ 5] = getUInt(baNbrList+offset+1*sizeCellGID,sizeCellGID,swapIntEndian);
      nbrIDArray[ 6] = getUInt(baNbrList+offset+2*sizeCellGID,sizeCellGID,swapIntEndian);
      nbrIDArray[ 7] = getUInt(baNbrList+offset+3*sizeCellGID,sizeCellGID,swapIntEndian);
      offset += 4*sizeCellGID;
   } else
     offset += sizeCellGID;
   // Parse -y neighbour IDs:
   nbrIDArray[ 8] = getUInt(baNbrList+offset,sizeCellGID,swapIntEndian);
   if (yNegNbrs == 4) {
      nbrIDArray[ 9] = getUInt(baNbrList+offset+1*sizeCellGID,sizeCellGID,swapIntEndian);
      nbrIDArray[10] = getUInt(baNbrList+offset+2*sizeCellGID,sizeCellGID,swapIntEndian);
      nbrIDArray[11] = getUInt(baNbrList+offset+3*sizeCellGID,sizeCellGID,swapIntEndian);
      offset += 4*sizeCellGID;
   } else
     offset += sizeCellGID;
   // Parse +y neighbour IDs:
   nbrIDArray[12] = getUInt(baNbrList+offset,sizeCellGID,swapIntEndian);
   if (yPosNbrs == 4) {
      nbrIDArray[13] = getUInt(baNbrList+offset+1*sizeCellGID,sizeCellGID,swapIntEndian);
      nbrIDArray[14] = getUInt(baNbrList+offset+2*sizeCellGID,sizeCellGID,swapIntEndian);
      nbrIDArray[15] = getUInt(baNbrList+offset+3*sizeCellGID,sizeCellGID,swapIntEndian);
      offset += 4*sizeCellGID;
   } else
     offset += sizeCellGID;
   // Parse -z neighbour IDs:
   nbrIDArray[16] = getUInt(baNbrList+offset,sizeCellGID,swapIntEndian);
   if (yPosNbrs == 4) {
      nbrIDArray[17] = getUInt(baNbrList+offset+1*sizeCellGID,sizeCellGID,swapIntEndian);
      nbrIDArray[18] = getUInt(baNbrList+offset+2*sizeCellGID,sizeCellGID,swapIntEndian);
      nbrIDArray[19] = getUInt(baNbrList+offset+3*sizeCellGID,sizeCellGID,swapIntEndian);
      offset += 4*sizeCellGID;
   } else
     offset += sizeCellGID;
   // Parse +z neighbour IDs:
   nbrIDArray[20] = getUInt(baNbrList+offset,sizeCellGID,swapIntEndian);
   if (yPosNbrs == 4) {
      nbrIDArray[21] = getUInt(baNbrList+offset+1*sizeCellGID,sizeCellGID,swapIntEndian);
      nbrIDArray[22] = getUInt(baNbrList+offset+2*sizeCellGID,sizeCellGID,swapIntEndian);
      nbrIDArray[23] = getUInt(baNbrList+offset+3*sizeCellGID,sizeCellGID,swapIntEndian);
      offset += 4*sizeCellGID;
   } else
     offset += sizeCellGID;
}

VlsHeader::Real (*ptrCellCRD)(const unsigned char* const,const bool&);
VlsHeader::Int (*ptrCellGID)(const unsigned char* const,const bool&);
VlsHeader::Int (*ptrVarNameSize)(const unsigned char* const,const bool&);

VlsHeader::Int VlsReader::getUInt(const unsigned char* const ptr,const int& size,const bool& swapIntEndian) const {
   switch (size) {
    case 1:
      return convUInt1(ptr,swapIntEndian);
      break;
    case 2:
      return convUInt2(ptr,swapIntEndian);
      break;
    case 4:
      return convUInt4(ptr,swapIntEndian);
      break;
    case 8:
      return convUInt8(ptr,swapIntEndian);
      break;
    default:
      return numeric_limits<VlsHeader::Int>::max();
      break;
   }
}

VlsHeader::Real VlsReader::getReal(const unsigned char* const ptr,const int& size,const bool& swapFloatEndian) const {
   switch (size) {
    case 4:
      return convReal4(ptr,swapFloatEndian);
      break;
    case 8:
      return convReal8(ptr,swapFloatEndian);
      break;
    case 16:
      return convReal16(ptr,swapFloatEndian);
      break;
    default:
      return NAN;
      break;
   }
}

enum VersionNumber {
   Unknown,
   Version_1_0_0
};

VersionNumber versionNumber = Unknown;

VersionNumber getVersionNumber(const std::string& stringVersion) {
   if (stringVersion == "1.0.0") return Version_1_0_0;
   return Unknown;
}

// ********************************************************
// ***** DEFINITIONS FOR CONSTRUCTORS AND DESTRUCTORS *****
// ********************************************************

VlsReader::VlsReader(): byteArray(NULL),fileptr(NULL) { 
   bytesPerSpatNbrListSize = 1;
   N_dimensions = 0;
   readDynamicData = false;
   readSpatNbrList = false;
   sizeCellCoordEntry = 0;
   sizeCellCRD = 0;
   sizeCellGID = 0;
   swapFloatEndian = false;
   swapIntEndian = false;
   
   // By default VLSV files are written in little-endian for both 
   // ints and floats:
   fileFloatsLittleEndian = true;
   fileIntsLittleEndian = true;
   
   // Detect the endianness on the machine in which vlsreader is run. 
   // For now I don't know how to test floating point endianness, so I use 
   // the same endianness for both ints and floats:
   const int number = 1;
   const char* const ptr = reinterpret_cast<const char*>(&number);
   if (ptr[0] == 1) {
      readerFloatsLittleEndian = true;
      readerIntsLittleEndian = true;
   } else {
      readerFloatsLittleEndian = false;
      readerIntsLittleEndian = false;
   }
}

VlsReader::~VlsReader() { 
   delete byteArray;
   delete fileptr;
   delete baDynData;
   delete baNbrList;
   byteArray = NULL;
   fileptr = NULL;
   baDynData = NULL;
   baNbrList = NULL;
}

VlsReader::VarDesc::VarDesc(const std::string& name,const unsigned char& typeID,const unsigned char& elementBytes):
name(name),typeID(typeID),elementBytes(elementBytes) {
   byteSize = elementBytes;
   if (typeID == VlsVariable::NULLVARIABLE)  {byteSize *= 0; elements = 0;}
   else if (typeID == VlsVariable::SCALAR)   {byteSize *= 1; elements = 1;}
   else if (typeID == VlsVariable::VECTOR2)  {byteSize *= 2; elements = 2;}
   else if (typeID == VlsVariable::VECTOR3)  {byteSize *= 3; elements = 3;}
   else if (typeID == VlsVariable::TENSOR22) {byteSize *= 4; elements = 4;}
   else if (typeID == VlsVariable::TENSOR23) {byteSize *= 6; elements = 6;}
   else if (typeID == VlsVariable::TENSOR32) {byteSize *= 6; elements = 6;}
   else if (typeID == VlsVariable::TENSOR33) {byteSize *= 9; elements = 9;}
}

// ********************************************
// ***** DEFINITIONS FOR MEMBER FUNCTIONS *****
// ********************************************

void VlsReader::calculateDerivedVariables() {
   // Calculate the size of spatial cell coordinate entry:
   sizeCellCoordEntry = sizeCellGID + N_dimensions*2*sizeCellCRD;
   
   for (vector<VarDesc>::const_iterator it=staticSizeVars.begin(); it!=staticSizeVars.end(); ++it)
     sizeCellCoordEntry += it->byteSize;
}

bool VlsReader::calculateVariableSizes() {
   bool success = true;
   const unsigned char* ptr;
   map<unsigned char,vector<unsigned char> >::const_iterator it;
   
   // Determine file endianness (if it was written):
   it = headerElements.find(VlsHeader::ENDIANNESS_FLOAT);
   if (it != headerElements.end()) {
      if (it->second[0] == VlsHeader::LITTLE_END) fileFloatsLittleEndian = true;
      else if (it->second[0] == VlsHeader::BIG_END) fileFloatsLittleEndian = false;
      else {
	 cerr << "VlsReader: erroneous value in ENDIANNESS_FLOAT!" << endl;
	 success = false;
      }
   }
   if (fileFloatsLittleEndian != readerFloatsLittleEndian) swapFloatEndian = true;
   
   it = headerElements.find(VlsHeader::ENDIANNESS_INT);
   if (it != headerElements.end()) {
      if (it->second[0] == VlsHeader::LITTLE_END) fileIntsLittleEndian = true;
      else if (it->second[0] == VlsHeader::BIG_END) fileIntsLittleEndian = false;
      else {
	 cerr << "VlsReader: erroneous value in ENDIANNESS_INT!" << endl;
	 success = false;
      }
   }
   if (fileIntsLittleEndian != readerIntsLittleEndian) swapIntEndian = true;

   // Make sure that endianness is taken into account when parsing 
   // the values of header elements!
   
   // Version number:
   it = headerElements.find(VlsHeader::VERSION);
   if (it == headerElements.end()) return false;
   ptr = &(it->second[0]);
   versionNumber = getVersionNumber(reinterpret_cast<const char*>(ptr));
   
   // Size of spatial cell coordinate entry:
   it = headerElements.find(VlsHeader::BYTES_PER_CELL_CRD);
   if (it == headerElements.end()) return false;
   sizeCellCRD = getUInt(&(it->second[0]),it->second.size(),swapIntEndian);
   
   // Size of spatial cell global ID entry:
   it = headerElements.find(VlsHeader::BYTES_PER_CELL_GID);
   if (it == headerElements.end()) return false;
   sizeCellGID = getUInt(&(it->second[0]),it->second.size(),swapIntEndian);
   
   // Number of spatial dimensions:
   it = headerElements.find(VlsHeader::DIMENSIONS);
   if (it == headerElements.end()) return false;
   N_dimensions = getUInt(&(it->second[0]),it->second.size(),swapIntEndian);
   
   // Size of entry containing the size of variable name:
   it = headerElements.find(VlsHeader::BYTES_PER_VARNAME_SIZE);
   if (it == headerElements.end()) return false;
   varNameFieldByteSize = getUInt(&(it->second[0]),it->second.size(),swapIntEndian);

   // Determine if file contains spatial neighbour lists:
   it = headerElements.find(VlsHeader::CONTAINS_SPAT_NBR_LIST);
   if (it != headerElements.end()) {
      VlsHeader::UInt tmp = getUInt(&(it->second[0]),it->second.size(),swapIntEndian);
      if (tmp != 0) readSpatNbrList = true;
   }
   
   // Determine if file contains dynamic-size data:
   it = headerElements.find(VlsHeader::CONTAINS_DYNAMIC_DATA);
   if (it != headerElements.end()) {
      VlsHeader::UInt tmp = getUInt(&(it->second[0]),it->second.size(),swapIntEndian);
      if (tmp != 0) readDynamicData = true;
   }

   // Size of field giving the byte size of spatial cell neighbour list entry:
   it = headerElements.find(VlsHeader::BYTES_PER_SPAT_NBR_SIZE);
   if (it != headerElements.end()) {
      bytesPerSpatNbrListSize = getUInt(&(it->second[0]),it->second.size(),swapIntEndian);
   }
   
   //headerElements.clear();
   return success;
}

/** Close a Vlasov file which has been previously opened by calling
 * VlsReader::open.
 * @return If true, the file was closed without errors.
 * @see VlsReader::open
 */
bool VlsReader::close() {
   if (fileptr != NULL) fileptr->close();
   delete fileptr;
   fileptr = NULL;
   
   delete byteArray;
   delete baDynData;
   delete baNbrList;
   byteArray = NULL;
   baDynData = NULL;
   baNbrList = NULL;
   sizeBaDynData = 0;
   sizeBaNbrList = 0;
   sizeByteArray = 0;
   return true;
}

VlsHeader::UInt VlsReader::getNeighbourID(const unsigned int& nbr) const {
   if (nbr >= sizeSpatNbrList/sizeCellGID) return numeric_limits<VlsHeader::UInt>::max();
   return getUInt(baNbrList+nbr*sizeCellGID,sizeCellGID,swapIntEndian);
}

VlsHeader::UInt VlsReader::getNumberOfSpatialNbrs() const {
   return sizeSpatNbrList/sizeCellGID;
}

size_t VlsReader::getNumberOfStaticVars() const {return staticSizeVars.size();}

VlsHeader::UInt VlsReader::getRefinementLevel() const {
   return cellRefLevel;
}

unsigned int VlsReader::getStaticVarElements(const size_t& varID) const {
   if (varID >= staticSizeVars.size()) return 0;
   return staticSizeVars[varID].elements;
}

std::string VlsReader::getStaticVarName(const size_t& varID) const {
   if (varID >= staticSizeVars.size()) return string("");
   return staticSizeVars[varID].name;
}
   
unsigned char VlsReader::getStaticVarType(const size_t& varID) const {
   if (varID >= staticSizeVars.size()) return VlsVariable::NULLVARIABLE;
   return staticSizeVars[varID].typeID;
}

/** Open a Vlasov file for reading data.
 * @param fname The name of the file.
 * @return If true, the file was opened successfully.
 */
bool VlsReader::open(const std::string& fname) {
   fileptr = new fstream;
   fileptr->open(fname.c_str(), fstream::in);
   if (fileptr->good() == false) return false;
   return true;
}

/** Read the header of Vlasov file. After calling this function the values 
 * of header elements can be queried with VlsReader::getHeaderElement member
 * function.
 * @return If true, this process read the header successfully.
 * @see VlsReader::getHeaderElement
 */
bool VlsReader::readHeader() {
   bool rvalue = true;
   headerElements.clear();

   // Read <byte length>-<ID byte>-<value> triplets until length=0 is encountered:
   unsigned char length;
   unsigned char byteArray[255];
   unsigned char byteID;
   fileptr->read(reinterpret_cast<char*>(&length),1);
   while (length > 0 && rvalue == true) {
      fileptr->read(reinterpret_cast<char*>(&byteID),1);
      if (fileptr->good() == false) rvalue = false;
      fileptr->read(reinterpret_cast<char*>(byteArray),length);
      if (fileptr->good() == false) rvalue = false;

      // Insert the read element into headerElements:
      vector<unsigned char> tmp(length);
      for (unsigned char i=0; i<length; ++i) tmp[i] = byteArray[i];
      headerElements[byteID] = tmp;
      
      fileptr->read(reinterpret_cast<char*>(&length),1);
      if (fileptr->good() == false) rvalue = false;
   }
   // Determine values for internal variables:
   if (rvalue == true) calculateVariableSizes();
   // Set int and float converters:
   if (rvalue == true) if (setInternalPointers() == false) rvalue = false;
   // Read static-size variable descriptions:
   if (rvalue == true) if (readStaticVariableDesc() == false) rvalue = false;
   // Calculate the values of internal variables that are derived 
   // from broadcasted values:
   if (rvalue == true) calculateDerivedVariables();
   
   return rvalue;
}

bool VlsReader::readSpatCellCoordEntry() {
   // Make sure byteArray is large enough:
   if (sizeByteArray < sizeCellCoordEntry) {
      delete byteArray;
      byteArray = new unsigned char[sizeCellCoordEntry];
      sizeByteArray = sizeCellCoordEntry;
   }
       
   // Attempt to read spatial cell GID:
   fileptr->read(reinterpret_cast<char*>(byteArray),sizeCellGID);
   if (fileptr->gcount() < sizeCellGID) return false;
   if (fileptr->good() == false) return false;
   spatCellGID = ptrCellGID(byteArray,swapIntEndian);

   // If all bits in GID are one, the GID marks the end of cell coordinate part:
   bool endMarkerFound = true;
   for (int i=0; i<sizeCellGID; ++i) if (byteArray[i] != 255) {endMarkerFound = false; break;}
   if (endMarkerFound == true) return false;

   // Read the rest of the cell entry:
   fileptr->read(reinterpret_cast<char*>(byteArray),sizeCellCoordEntry-sizeCellGID);
   if (fileptr->gcount() != (sizeCellCoordEntry-sizeCellGID)) return false;
   if (fileptr->good() == false) return false;
   
   // Store values to internal variables, taking care of endianness:
   if (versionNumber > Unknown) {
      xcrd        = ptrCellCRD(byteArray                ,swapFloatEndian);
      ycrd        = ptrCellCRD(byteArray + 1*sizeCellCRD,swapFloatEndian);
      zcrd        = ptrCellCRD(byteArray + 2*sizeCellCRD,swapFloatEndian);
      dx          = ptrCellCRD(byteArray + 3*sizeCellCRD,swapFloatEndian);
      dy          = ptrCellCRD(byteArray + 4*sizeCellCRD,swapFloatEndian);
      dz          = ptrCellCRD(byteArray + 5*sizeCellCRD,swapFloatEndian);
   }

   // Read spatial neighbour list (if it is included):
   if (readSpatNbrList == true) {
      // Read the byte size of spatial neighbour list:
      unsigned char tmp;
      fileptr->read(reinterpret_cast<char*>(&tmp),1);
      if (fileptr->gcount() != 1) return false;
      sizeSpatNbrList = tmp;
      
      // Check that nbr list byte array is large enough:
      if (sizeBaNbrList < sizeSpatNbrList) {
	 delete baNbrList;
	 baNbrList = new unsigned char[sizeSpatNbrList];
	 sizeBaNbrList = sizeSpatNbrList;
      }
      // Read cell refinement level:
      fileptr->read(reinterpret_cast<char*>(&cellRefLevel),1);
      if (fileptr->gcount() != 1) return false;
      // Read neighbour list:
      fileptr->read(reinterpret_cast<char*>(baNbrList),sizeSpatNbrList);
      if (fileptr->gcount() != sizeSpatNbrList) return false;
   }
   
   // Read dynamic-size data if it is included:
   if (readDynamicData == true) {
      
   }
   
   return true;
}

bool VlsReader::readStaticVariableDesc() {
   staticSizeVars.clear();
   if (varNameFieldByteSize == 0) return false;
   
   // Reserve large enough byte array:
   bool success = true;
   delete byteArray;
   byteArray = new unsigned char[65536];
   sizeByteArray = 65536;

   // Read variable descriptions while the size of variable name > 0:
   do {
      // Read the byte size of variable name:
      fileptr->read(reinterpret_cast<char*>(byteArray),varNameFieldByteSize);
      if (fileptr->gcount() != varNameFieldByteSize) {success = false; varNameSize = 0;}
      if (fileptr->good() == false) {success = false; varNameSize = 0;}      
      varNameSize = ptrVarNameSize(byteArray,swapIntEndian);
      
      // If variable exists, read its description:
      if (varNameSize > 0) {
	 unsigned long int bytes = varNameSize + 2; // size(name) + typeID + elementBytes
	 fileptr->read(reinterpret_cast<char*>(byteArray),bytes);
	 if (fileptr->gcount() != bytes) {success = false; varNameSize = 0;}
	 if (fileptr->good() == false) {success = false; varNameSize = 0;}
      }
   
      // Parse description:
      string varName;
      unsigned char varTypeID;
      unsigned char varElementSize;
      varName = reinterpret_cast<char*>(byteArray);
      varTypeID      = getUInt(byteArray + varName.size() + 1,1,swapIntEndian); // +1 from end-of-character
      varElementSize = getUInt(byteArray + varName.size() + 2,1,swapIntEndian);
      staticSizeVars.push_back(VarDesc(varName,varTypeID,varElementSize));
   } while (varNameSize > 0);
   
   // Calculate static size variable byte offsets:
   if (staticSizeVars.size() > 0) {
      staticSizeVars[0].offset = 2*N_dimensions*sizeCellCRD;
      for (size_t i=1; i<staticSizeVars.size(); ++i) {
	 staticSizeVars[i].offset = staticSizeVars[i-1].offset + staticSizeVars[i-1].byteSize;
      }
   }
   
   delete byteArray;
   byteArray = NULL;
   sizeByteArray = 0;
   return success;
}
   
bool VlsReader::setInternalPointers() {
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




