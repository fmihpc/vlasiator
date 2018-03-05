/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://www.physics.helsinki.fi/vlasiator/
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * File:   mesh_data.h
 * Author: sandroos
 *
 * Created on June 2, 2015, 10:09 AM
 */

#ifndef MESH_DATA_H
#define MESH_DATA_H

#include <cstdlib>

namespace mesh {

   class MeshData {
   public:
      MeshData(): byteSize(0),vectorSize(0),dataPointer(NULL),meshSize(0) { }
      MeshData(const MeshData& mdc): byteSize(mdc.byteSize),vectorSize(mdc.vectorSize),dataType(mdc.dataType),meshSize(mdc.meshSize) {
         dataPointer = new char[byteSize*vectorSize*meshSize];
         for (size_t i=0; i<byteSize*vectorSize*meshSize; ++i) dataPointer[i] = mdc.dataPointer[i];
      }
      ~MeshData() {delete [] dataPointer; dataPointer = NULL;}

      template<typename T> T* getData();
      size_t getDataSize() const {return byteSize;}
      const std::string& getDataType() const {return dataType;}
      size_t getVectorSize() const {return vectorSize;}
      template<typename T> T index(const T& i,const T& j,const T& k) const;
      void reallocate() {
         delete [] dataPointer;
         dataPointer = new char[meshSize*vectorSize*byteSize];
      }
      template<typename T> bool setDataSize(const size_t& vectorSize,const std::string& datatype);
      
      /** Set the number of cells in the mesh.
       * @param meshSize Number of cells in the mesh.
       * @return If true, number of cells was set successfully.*/
      bool setMeshSize(const size_t& meshSize) {this->meshSize = meshSize; return true;}

   private:
      size_t byteSize;           /**< Byte size of one vector element.*/
      size_t vectorSize;         /**< Data vector size, number of elements per cell.*/
      std::string dataType;      /**< String representation of the stored datatype.*/
      char* dataPointer;         /**< Pointer to the data array.*/
      size_t meshSize;           /**< Number of cells in mesh.*/
   };

   /** Get pointer to the data. The template parameter must be of the same 
    * datatype that was passed to the setDataSize member function.
    * @return Pointer to the data.*/
   template<typename T> inline
   T* MeshData::getData() {
      return reinterpret_cast<T*>(dataPointer);
   }

   /** Set the data vector metadata and allocate memory. A vector is stored in 
    * each cell in the mesh. The vector has vectorSize elements, and the datatype 
    * of each element is given in the template parameter.
    * @param vectorSize Size of the data vector stored in each cell.
    * @return If true, memory for the data was allocated successfully.*/
   template<typename T> inline
   bool MeshData::setDataSize(const size_t& vectorSize,const std::string& datatype) {
      delete [] dataPointer;
      byteSize = sizeof(T);
      this->vectorSize = vectorSize;
      dataPointer = new char[meshSize*vectorSize*byteSize];

      this->dataType = datatype;
      return true;
   }

} // namespace mesh

#endif

