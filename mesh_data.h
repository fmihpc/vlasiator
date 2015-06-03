/*
 * This file is part of Vlasiator.
 * 
 * Copyright 2015 Finnish Meteorological Institute
 * File:   mesh_data.h
 * Author: sandroos
 *
 * Created on June 2, 2015, 10:09 AM
 */

#ifndef MESH_DATA_H
#define	MESH_DATA_H

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
      template<typename T> bool setDataSize(const size_t& vectorSize);
      
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
   bool MeshData::setDataSize(const size_t& vectorSize) {
      delete [] dataPointer;
      byteSize = sizeof(T);
      this->vectorSize = vectorSize;
      dataPointer = new char[meshSize*vectorSize*byteSize];
      
#warning FIXME
      dataType = "float";
      return true;
   }

} // namespace mesh

#endif

