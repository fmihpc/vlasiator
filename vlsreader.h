/*
This file is part of Vlasiator.

Copyright 2010, 2011 Finnish Meteorological Institute

Vlasiator is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3
as published by the Free Software Foundation.

Vlasiator is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef VLSREADER_H
#define VLSREADER_H

#include <cmath>
#include <fstream>
#include <map>
#include <vector>

#include "vlscommon.h"

class VlsReader {
 public:
   VlsReader();
   ~VlsReader();
   
   bool close();
   template<typename T> T getCrdX() const;
   template<typename T> T getCrdY() const;
   template<typename T> T getCrdZ() const;
   template<typename T> T getDx() const;
   template<typename T> T getDy() const;
   template<typename T> T getDz() const;
   template<typename T> bool getHeaderElement(const unsigned char& ID,T& value) const;
   VlsHeader::UInt getNeighbourID(const unsigned int& nbr) const;
   VlsHeader::UInt getNumberOfSpatialNbrs() const;
   size_t getNumberOfStaticVars() const;
   VlsHeader::UInt getRefinementLevel() const;
   template<typename T> T getSpatCellGID() const;
   template<typename T> T getStaticVar(const size_t& varID,const unsigned int& element) const;
   unsigned int getStaticVarElements(const size_t& varID) const;
   std::string getStaticVarName(const size_t& varID) const;
   unsigned char getStaticVarType(const size_t& varID) const;
   
   bool open(const std::string& fname);
   bool readHeader();
   bool readSpatCellCoordEntry();
   
 private:
   /** Description of a static-size variable as read from a vlsv file.*/
   struct VarDesc {
      VlsHeader::Int byteSize;    /**< Size of the variable data array, in bytes.
				   * Equals to the number of elements times VarDesc::elementBytes.*/
      VlsHeader::Int offset;      /**< Static-size variable data for each cell is stored as a
				   * continuous byte array. This is an offset relative to the
				   * start of the array, which tells where the data of this
				   * variable begins.*/
      unsigned char elements;     /**< Number of elements in the variable data array.
				   Deduced from VarDesc::typeID.*/
      unsigned char elementBytes; /**< Size of an element in variable data array, in bytes.*/
      std::string name;           /**< The name of the variable.*/
      unsigned char typeID;       /**< Type ID of the variable, corresponds to a value defined
				   in namespace VlsHeader.*/
      
      VarDesc(const std::string& name,const unsigned char& typeID,const unsigned char& elementBytes);
   };

   unsigned char* byteArray;
   VlsHeader::UInt bytesPerSpatNbrListSize; /**< Byte size of field that gives the byte size of spatial neighbour 
					     * list entry. Defaults to unit value.*/
   unsigned char cellRefLevel;
   VlsHeader::Real dx;                  /**< Width of a cell in x-direction.*/
   VlsHeader::Real dy;                  /**< Width of a cell in y-direction.*/
   VlsHeader::Real dz;                  /**< Width of a cell in z-direction.*/
   std::fstream* fileptr;
   std::map<unsigned char,std::vector<unsigned char> > headerElements; /**< A map which contains header element IDs and their values
									* after VlsWriter::readHeader has been called.*/
   unsigned char N_dimensions;          /**< Dimensionality of the data stored in the vlsv file.*/
   bool readDynamicData;                /**< If true, cell entries in VLSV file contain dynamic-size data.
					 * Defaults to false.*/
   bool readSpatNbrList;                /**< If true, cell entries in VLSV file contain spatial neighbour lists.
					 * Defaults to false.*/
   unsigned int sizeByteArray;
   unsigned int sizeCellCoordEntry;     /**< The size of spatial cell coordinate entry, in bytes.*/
   unsigned char sizeCellCRD;           /**< Size of spatial cell coordinate value entry, in bytes.*/
   unsigned char sizeCellGID;           /**< Size of spatial cell global ID entry, in bytes.*/
   VlsHeader::UInt sizeSpatNbrList;     /**< Byte size of spatial neighbour list entry.*/
   VlsHeader::Int spatCellGID;          /**< Global ID of the cell which was read from the vlsv file.*/
   std::vector<VarDesc> staticSizeVars; /**< Container for the descriptions of static-size variables
					 * stored in the vlsv file.*/
   bool swapFloatEndian;                /**< If true, floating point byte order should be swapped when
					 * reading values from file.*/
   bool swapIntEndian;                  /**< If true, integer byte order should be swapped when
					 * reading values from file.*/
   unsigned char varNameFieldByteSize;  /**< The size (in bytes) of an entry which gives the size of
					 * a variable name.*/
   VlsHeader::Int varNameSize;          /**< The byte size of a variable name. Since variable names are
					 * just char arrays, this equals to the number of characters in
					 * the name (including the end-of-character).*/
   VlsHeader::Real xcrd;                /**< X-coordinate of the lower left corner of a cell which
					 * was read from a vlsv file.*/
   VlsHeader::Real ycrd;                /**< Y-coordinate of the lower left corner of a cell which
					 * was read from a vlsv file.*/
   VlsHeader::Real zcrd;                /**< Z-coordinate of the lower left corner of a cell which
7					 * was read from a vlsv file.*/
   
   void calculateDerivedVariables();
   bool calculateVariableSizes();
   VlsHeader::Int getUInt(const unsigned char* const ptr,const int& size,const bool& swapIntEndian) const;
   VlsHeader::Real getReal(const unsigned char* const ptr,const int& size,const bool& swapFloatEndian) const;
   void parseNeighbourIndices();
   bool readStaticVariableDesc();
   bool setInternalPointers();
};

// **********************************
// ****** TEMPLATE DEFINITIONS ******
// **********************************

template<typename T> T VlsReader::getCrdX() const {return static_cast<T>(xcrd);}
template<typename T> T VlsReader::getCrdY() const {return static_cast<T>(ycrd);}
template<typename T> T VlsReader::getCrdZ() const {return static_cast<T>(zcrd);}
template<typename T> T VlsReader::getDx() const {return static_cast<T>(dx);}
template<typename T> T VlsReader::getDy() const {return static_cast<T>(dy);}
template<typename T> T VlsReader::getDz() const {return static_cast<T>(dz);}

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
inline bool VlsReader::getHeaderElement(const unsigned char& ID,T& value) const {
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
inline bool VlsReader::getHeaderElement(const unsigned char& ID,std::string& value) const {
   // Check that the given element has been read:
   std::map<unsigned char,std::vector<unsigned char> >::const_iterator it = headerElements.find(ID);
   if (it == headerElements.end()) return false;
   // Cast each byte into char:
   value = reinterpret_cast<char*>(const_cast<unsigned char*>(&((it->second)[0])));
   return true;
}

template<typename T> 
T VlsReader::getSpatCellGID() const {return static_cast<T>(spatCellGID);}

template<typename T> 
T VlsReader::getStaticVar(const size_t& varID,const unsigned int& element) const {
   if (varID >= staticSizeVars.size()) return NAN;
   if (element >= staticSizeVars[varID].elements) return NAN;
   // Return the value of the requested element, taking care of endianness:
   return getReal(&(byteArray[staticSizeVars[varID].offset + element*staticSizeVars[varID].elementBytes]),
		  staticSizeVars[varID].elementBytes,swapFloatEndian);
}

#endif
