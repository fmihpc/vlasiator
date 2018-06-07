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
 * File:   mesh_data_container.h
 * Author: sandroos
 *
 * Created on June 2, 2015, 11:18 AM
 */

#ifndef MESH_DATA_CONTAINER_H
#define	MESH_DATA_CONTAINER_H

#include <cstdlib>
#include <map>
#include <vector>

#include "mesh_data.h"

#ifndef AMR
   #include "velocity_mesh_old.h"
#else
   #include "velocity_mesh_amr.h"
#endif

namespace mesh {

   class MeshDataContainer {
      public:
         MeshDataContainer();
         
         template<typename T> size_t addData(const std::string& name,const size_t& vectorSize,const std::string& datatype);
         template<typename T> T* getData(const size_t& dataID);
         template<typename T> T* getData(const std::string& name);
         size_t getDataSize(const size_t& dataID) const {return meshData[dataID].getDataSize();}
         const std::string& getDataType(const size_t& dataID) const {return meshData[dataID].getDataType();}
         size_t getMeshSize() const {return mesh.size();}
         size_t getVectorSize(const size_t& dataID) const {return meshData[dataID].getVectorSize();}
         std::string getName(const size_t& dataID) const {
            for (std::map<std::string,size_t>::const_iterator it=meshDataNames.begin(); it!=meshDataNames.end(); ++it) {
               if (it->second == dataID) return it->first;
            }
            return "";
         }
         unsigned int getLocalID(const CellID& cellID) const {return mesh.getLocalID(cellID);}
         bool initialize(const std::string& name);
         void reallocate();
         size_t size() const {return meshData.size();}

      private:

         bool initialized;
         vmesh::VelocityMesh<CellID,unsigned int> mesh;
         std::vector<mesh::MeshData> meshData;
         std::map<std::string,size_t> meshDataNames;
         size_t meshID;
   };
   
   template<typename T> inline
   size_t MeshDataContainer::addData(const std::string& name,const size_t& vectorSize,const std::string& datatype) {
      // Return with an invalid data ID if class has not been initialized
      if (initialized == false) return std::numeric_limits<size_t>::max();
      
      // Check that data doesn't already exist
      std::map<std::string,size_t>::const_iterator it = meshDataNames.find(name);
      if (it != meshDataNames.end()) return it->second;

      size_t rvalue = meshData.size();
      meshDataNames[name] = rvalue;
      meshData.push_back(mesh::MeshData());
      meshData[rvalue].setMeshSize(mesh.size());
      meshData[rvalue].setDataSize<T>(vectorSize,datatype);
      return rvalue;
   }
   
   template<typename T> inline
   T* MeshDataContainer::getData(const size_t& dataID) {
      // Return with a NULL pointer if class has not been initialized
      if (initialized == false) return NULL;

      // Return pointer to the data
      return meshData[dataID].getData<T>();
   }
   
   template<typename T> inline
   T* MeshDataContainer::getData(const std::string& name) {
      // Return with a NULL pointer if class has not been initialized
      if (initialized == false) return NULL;
      
      // If data doesn't exist return a NULL pointer.
      std::map<std::string,size_t>::const_iterator it = meshDataNames.find(name);
      if (it == meshDataNames.end()) return NULL;
      
      // Return pointer to the data
      return meshData[it->second].getData<T>();      
   }
   
} // namespace mesh
#endif

