/*
This file is part of Vlasiator.

Copyright 2010, 2011, 2012 Finnish Meteorological Institute

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


/*! \file vlsvreader2.h
 * \brief Classes for reading in vlsv files
 */


#ifndef VLSVREADER2_H
#define VLSVREADER2_H

#ifndef TOOL_NOT_PARALLEL
   #include <mpi.h>
#endif
#include <stdint.h>
#include <list>
#include <fstream>
#include "muxml.h"
#include "vlscommon.h"


/*!
\brief A serial non-MPI class for reading in vlsv files.

*/
class VLSVReader {
 public:
   VLSVReader();
   virtual ~VLSVReader();
   /*!
     \brief Close file
     \return Returns true for success, and false for error
   */ 
   virtual bool close();
   /*!
     \brief Get information on array identified by its tag and name 
     \param[in] tagName Name of tag
     \param[in] arrayName Name of array (attribute name)
     \param[out] arraySize Number of vector elements in the array
     \param[out] vectorSize Number of scalar elements in each vector element 
     \param[out] dataType Datatype of each scalar element in the vectors ("int","uint", or "float")
     \param[out] byteSize The size of each scalar element in bytes.
     \return Returns true for success, and false for error
   */   
   virtual bool getArrayInfo(const std::string& tagName,const std::string& arrayName,uint64_t& arraySize,
			     uint64_t& vectorSize,VLSV::datatype& dataType,uint64_t& byteSize) const;
   /*!
     \brief Get information on array ientified by its tag, name meshName
     \param[in] tagName Name of tag
     \param[in] arrayName Name of array (attribute "name")
     \param[in] meshName Name of mesh on which the data is defind (attribute "mesh" )
     \param[out] arraySize Number of vector elements in the array
     \param[out] vectorSize Number of scalar elements in each vector element 
     \param[out] dataType Datatype of each scalar element in the vectors ("int","uint", or "float")
     \param[out] byteSize The size of each scalar element in bytes.
     \return Returns true for success, and false for error
   */
   virtual bool getArrayInfo(const std::string& tagName,const std::string& arrayName,const std::string& meshName,
			     uint64_t& arraySize,uint64_t& vectorSize,VLSV::datatype& dataType,uint64_t& byteSize) const;

   /*!
     \brief Get information on array identified by attributes
     \param[in] tagName Name of tag
     \param[in] attribs List of attribute pairs. Key can be, e.g., "name" and "mesh.
     \param[out] arraySize Number of vector elements in the array
     \param[out] vectorSize Number of scalar elements in each vector element 
     \param[out] dataType Datatype of each scalar element in the vectors ("int","uint", or "float")
     \param[out] byteSize The size of each scalar element in bytes.
     \return Returns true for success, and false for error
   */
   virtual bool getArrayInfo(const std::string& tagName,const std::list<std::pair<std::string,std::string> >& attribs,
			     uint64_t& arraySize,uint64_t& vectorSize,VLSV::datatype& dataType,uint64_t& byteSize) const;


   /*!
     \brief Get the name of all arrays with a "MESH" tag
     \param[out] meshNames List of values for "name" attributes for all "MESH" tags
     \return Returns true for success, and false for error
   */
   virtual bool getMeshNames(std::list<std::string>& meshNames) const;
   /*!
     \brief Get the name of all arrays with a "VARIABLE" tag and a specified "mesh" attribute
     \param[in] meshNames "mesh" atribute value
     \param[out] varNames List of values for "name" attributes for all matching arrays
     \return Returns true for success, and false for error
   */
   virtual bool getVariableNames(const std::string& meshName,std::list<std::string>& varNames) const;
   /*!
     \brief Get the name of all arrays with a "BLOCKVARIABLE" tag and a specified "mesh" attribute
     \param[in] meshNames "mesh" atribute value
     \param[out] varNames List of values for "name" attributes for all matching arrays
     \return Returns true for success, and false for error
   */

   virtual bool getBlockVariableNames(const std::string& meshName,std::list<std::string>& varNames) const;
/*!     
    \brief Open file      
    \param[in] fname Name of file that is to be opened
     \return Returns true for success, and false for error
  */
   virtual bool open(const std::string& fname);


   /*!
     \brief Read array
     \param[in] tagName Name of tag
     \param[in] arrayName Name of array (value of "name" attrib).
     \param[in] begin Location in array where read is to begin, in units of array vectors.
     \param[in] amount Number of vectors that is read
     \param[out] Pointer to array where the read data is placed
     \return Returns true for success, and false for error 
   */
   virtual bool readArray(const std::string& tagName,const std::string& arrayName,const uint64_t& begin,
                          const uint64_t& amount,char* buffer);
   /*!
     \brief Read array
     \param[in] tagName Name of tag
     \param[in] attribs List of attribute pairs. Key can be, e.g., "mesh".
     \param[in] begin Location in array where read is to begin, in units of array vectors.
     \param[in] amount Number of vectors that is read
     \param[out] Pointer to array where the read data is placed
     \return Returns true for success, and false for error 
   */
   virtual bool readArray(const std::string& tagName,const std::list<std::pair<std::string,std::string> >& attribs,
			  const uint64_t& begin,const uint64_t& amount,char* buffer);

   /*!
     \brief Read array
     \param[in] tagName Name of tag
     \param[in] arrayName Name of array (value of "name" attrib).
     \param[in] attribs List of attribute pairs. Key can be, e.g., "mesh".
     \param[in] begin Location in array where read is to begin, in units of array vectors.
     \param[in] amount Number of vectors that is read
     \param[out] Pointer to array where the read data is placed
     \return Returns true for success, and false for error 
   */
   virtual bool readArray(const std::string& tagName,const std::string& arrayName,
                          const std::list<std::pair<std::string,std::string> >& attribs,
			  const uint64_t& begin,const uint64_t& amount,char* buffer);
   
 protected:
   unsigned char endiannessFile;   /**< Endianness in VLSV file.*/
   unsigned char endiannessReader; /**< Endianness of computer which reads the data.*/
   std::fstream filein;            /**< Input file stream.*/
   std::string fileName;           /**< Name of the input file.*/
   bool fileOpen;                  /**< If true, a file is currently open.*/
   bool swapIntEndianness;         /**< If true, endianness should be swapped on read data (not implemented yet).*/
   MuXML xmlReader;                /**< XML reader used to parse VLSV footer.*/


   
   virtual bool loadArray(const std::string& tagName,const std::list<std::pair<std::string,std::string> >& attribs);
   
   /*! Struct used to store information on the currently open array.*/
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

#ifndef TOOL_NOT_PARALLEL
/*!
\brief A parallel MPI class for reading in vlsv files.

*/
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
#endif // #ifndef TOOL_NOT_PARALLEL
#endif // #ifndef VLSVREADER2_H
