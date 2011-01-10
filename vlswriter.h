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

/** A namespace which defines ID numbers for vlsv file header elements.
 * You should not change the values of existing constants as this will break 
 * backward compability with older file versions, you can only add new values.
 */
namespace VlsHeader {
   typedef uint64_t Int;
   typedef double Real;
   
   const unsigned char BYTES_PER_CELL_CRD     = 0; /**< Bytes per spatial cell coordinate value (>0).*/
   const unsigned char BYTES_PER_CELL_GID     = 1; /**< Bytes per spatial cell global ID (>0).*/
   const unsigned char DIMENSIONS             = 2; /**< Spatial dimensionality of data (1, 2, or 3).*/
   const unsigned char VERSION                = 3; /**< Version number of the file (string).*/
   const unsigned char BYTES_PER_VARNAME_SIZE = 4; /**< Byte size of field that gives the byte size of variable name entry.*/
   const unsigned char STATIC_SIZE_VARIABLES  = 5; /**< Description of static-size spatial cell variables, i.e. variables having the 
						    * same size for each spatial cell.*/
   
   const unsigned char NULLVARIABLE = 0; /**< Variable is null, i.e. an array with zero elements.*/
   const unsigned char SCALAR   = 1;     /**< Variable is a scalar, i.e. an array with one element.*/
   const unsigned char VECTOR2  = 2;     /**< Variable is a 2D vector, i.e. an array with two elements.*/
   const unsigned char VECTOR3  = 3;     /**< Variable is a 3D vector, i.e. an array with three elements.*/
   const unsigned char TENSOR22 = 4;     /**< Variable is a 2x2 tensor, i.e. an array with four elements.*/
   const unsigned char TENSOR23 = 5;     /**< Variable is a 2x3 tensor, i.e. an array with five elements.*/
   const unsigned char TENSOR32 = 6;     /**< Variable is a 3x2 tensor, i.e. an array with five elements.*/
   const unsigned char TENSOR33 = 7;     /**< Variable is a 3x3 tensor, i.e. an array with six elements.*/
   
   const unsigned char BEGIN_HEADER               = 0;
   const unsigned char BEGIN_STATIC_VARIABLE_DESC = 1;
   const unsigned char BEGIN_COORDs               = 2;
}

/**
 * Definition of vlsv file header:
 * Header consists of <byte size> <header ID> <value> triplets, where 
 * <byte size> and <header ID> are unsigned chars, i.e. bytes. <value> is a 
 * byte array of length <byte size>. <header ID> gives the ID number of the 
 * entry, and corresponds to one (and just one) numeric value of constants 
 * defined in namespace VlsHeader.
 * 
 * Note that MPI must have been initialized prior to using this class.
 */
class VlsWriter {
 public:
   
   VlsWriter();
   ~VlsWriter();
   
   bool close();
   template<typename T> T getCrdX() const;
   template<typename T> T getCrdY() const;
   template<typename T> T getCrdZ() const;
   template<typename T> T getDx() const;
   template<typename T> T getDy() const;
   template<typename T> T getDz() const;
   template<typename T> bool getHeaderElement(const unsigned char& ID,T& value) const;
   size_t getNumberOfStaticVars() const;
   template<typename T> T getSpatCellGID() const;
   template<typename T> T getStaticVar(const size_t& varID,const unsigned int& element) const;
   unsigned int getStaticVarElements(const size_t& varID) const;
   std::string getStaticVarName(const size_t& varID) const;
   unsigned char getStaticVarType(const size_t& varID) const;
   
   bool openRead(MPI_Comm comm,const std::string& fname);
   bool openWrite(MPI_Comm comm,const std::string& fname);
   bool readHeader(MPI_Comm comm,const int& masterRank);
   bool readSpatCellCoordEntry();
   bool readStaticVariableDesc(MPI_Comm comm,const int& masterRank);
   bool sync();
   bool writeHeader(MPI_Comm comm,const int& masterRank);   
   bool writeSpatCellCoordEntry(cuint& cellID,const SpatialCell& cell,DataReducer* const dr);
   bool writeStaticVariableDesc(MPI_Comm comm,const int& masterRank,DataReducer* const dr);

 private:
   struct VarDesc {
      std::string name;
      unsigned char typeID;
      unsigned char elements;
      unsigned char elementBytes;
      VlsHeader::Int byteSize;
      VlsHeader::Int offset;
      
      VarDesc(const std::string& name,const unsigned char& typeID,const unsigned char& elementBytes);
   };
   
   unsigned char* byteArray;
   VlsHeader::Real dx;
   VlsHeader::Real dy;
   VlsHeader::Real dz;
   std::map<unsigned char,std::vector<unsigned char> > headerElements; /**< A map which contains header element IDs and their values
									* after VlsWriter::readHeader has been called.*/
   MPIFile mpiFile;                  /** MPIFile which is used for I/O.*/
   unsigned char N_dimensions;
   unsigned int sizeByteArray;
   unsigned int sizeCellCoordEntry;  /**< The size of spatial cell coordinate entry, in bytes.*/
   unsigned char sizeCellCRD;        /**< Size of spatial cell coordinate value entry, in bytes.*/
   unsigned char sizeCellGID;        /**< Size of spatial cell global ID entry, in bytes.*/
   VlsHeader::Int  spatCellGID;
   std::vector<VarDesc> staticSizeVars;
   unsigned char varNameFieldByteSize;
   VlsHeader::Int varNameSize;
   VlsHeader::Real xcrd;
   VlsHeader::Real ycrd;
   VlsHeader::Real zcrd;
   
   template<typename T> void appendHeaderElement(std::vector<unsigned char>& byteArray,const unsigned char& typeID,const T& element);
   void appendHeaderElement(std::vector<unsigned char>& byteArray,const unsigned char& typeID,const unsigned char* const array,const unsigned int& arraySize);
   bool broadcastVariableValues(MPI_Comm comm,const int& masterRank);
   void calculateDerivedVariables();
   bool calculateVariableSizes();
   bool setInternalPointers();
};
   
// **********************************
// ****** TEMPLATE DEFINITIONS ******
// **********************************

template<typename T> T VlsWriter::getCrdX() const {return static_cast<T>(xcrd);}
template<typename T> T VlsWriter::getCrdY() const {return static_cast<T>(ycrd);}
template<typename T> T VlsWriter::getCrdZ() const {return static_cast<T>(zcrd);}
template<typename T> T VlsWriter::getDx() const {return static_cast<T>(dx);}
template<typename T> T VlsWriter::getDy() const {return static_cast<T>(dy);}
template<typename T> T VlsWriter::getDz() const {return static_cast<T>(dz);}

template<typename T> T VlsWriter::getSpatCellGID() const {return static_cast<T>(spatCellGID);}

template<typename T> T VlsWriter::getStaticVar(const size_t& varID,const unsigned int& element) const {
   if (varID >= staticSizeVars.size()) return NAN;
   if (element >= staticSizeVars[varID].elements) return NAN;
   /*
   if (staticSizeVars[varID].byteSize == 4) {
      return static_cast<T>(*(reinterpret_cast<float*>(&(byteArray[staticSizeVars[varID].offset + element*4]))));
   } else if (staticSizeVars[varID].byteSize == 8) { 
      return static_cast<T>(*(reinterpret_cast<double*>(&(byteArray[staticSizeVars[varID].offset + element*8]))));
   } else if (staticSizeVars[varID].byteSize == 16) {
      return static_cast<T>(*(reinterpret_cast<long double*>(&(byteArray[staticSizeVars[varID].offset + element*16]))));
   } else {
      return NAN;
   }
   */
   if (staticSizeVars[varID].elementBytes == 4) {
      return static_cast<T>(*(reinterpret_cast<float*>(&(byteArray[staticSizeVars[varID].offset + element*4]))));
   } else if (staticSizeVars[varID].elementBytes == 8) {
      return static_cast<T>(*(reinterpret_cast<double*>(&(byteArray[staticSizeVars[varID].offset + element*8]))));
   } else if (staticSizeVars[varID].elementBytes == 16) {
      return static_cast<T>(*(reinterpret_cast<long double*>(&(byteArray[staticSizeVars[varID].offset + element*16]))));
   } else {
      return NAN;
   }
}

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

/** Attempt to fetch the value of header element with the given ID. The byte size
 * parameter 'value' must match the byte size of the header element in file. For example,
 * if the given header element is to the file with four bytes, then sizeof(value) must also be 
 * equal to four. The only exception to this rule are when reading strings 
 * (see the specialized template for strings). The byte array read from file is type casted 
 * into primitive type represented by parameter 'value'.
 * @param ID The ID of the requested header element.
 * @param value A parameter in which the value of the header element is written.
 * @return If true, the value of the header element was found. Otherwise the file 
 * does not contain the given header element, or the byte size of value does not match the 
 * byte size of stored value.
 */
template<typename T>
inline bool VlsWriter::getHeaderElement(const unsigned char& ID,T& value) const {
   // Check that the given element has been read, and that value has the correct size:
   std::map<unsigned char,std::vector<unsigned char> >::const_iterator it = headerElements.find(ID);
   if (it == headerElements.end()) return false;
   if (it->second.size() != sizeof(value)) return false;
   // Cast the byte array into same type as value:
   T* ptr = reinterpret_cast<T*>(const_cast<unsigned char*>(&((it->second)[0])));
   value = *ptr;
   return true;
}

/** Attempt to fetch the value of header element with the given ID as a string (char array).
 * @param ID The ID of the requested header element.
 * @param value A string in which the value of the header element is written. The header element 
 * is written as a byte array, and each of those bytes is type casted into chars and stored to 
 * value.
 * @return If true, the value of the header element was found. Otherwise the file
 * does not contain the given header element.
 */
template<>
inline bool VlsWriter::getHeaderElement(const unsigned char& ID,std::string& value) const {
   // Check that the given element has been read:
   std::map<unsigned char,std::vector<unsigned char> >::const_iterator it = headerElements.find(ID);
   if (it == headerElements.end()) return false;
   // Cast each byte into char:
   value = reinterpret_cast<char*>(const_cast<unsigned char*>(&((it->second)[0])));
   return true;
}

#endif
