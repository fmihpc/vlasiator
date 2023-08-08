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
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 * 
 * File:   mesh_data_container.cpp
 * Author: sandroos
 *
 * Created on June 2, 2015, 11:18 AM
 */

#include "mesh_data_container.h"
#include "object_wrapper.h"

using namespace std;

extern ARCH_MANAGED GridParameters meshParams;

namespace mesh {

   MeshDataContainer::MeshDataContainer(): initialized(false) { }
   
   bool MeshDataContainer::initialize(const std::string& name) {
      vmesh::MeshParameters vmeshParams;
      vmeshParams.name = name;
      vmeshParams.max_velocity_blocks = numeric_limits<vmesh::LocalID>::max();
      vmeshParams.meshLimits[0] = meshParams.xmin;
      vmeshParams.meshLimits[1] = meshParams.xmax;
      vmeshParams.meshLimits[2] = meshParams.zmin;
      vmeshParams.meshLimits[3] = meshParams.zmax;
      vmeshParams.meshLimits[4] = meshParams.zmin;
      vmeshParams.meshLimits[5] = meshParams.zmax;
      vmeshParams.gridLength[0] = meshParams.xcells_ini;
      vmeshParams.gridLength[1] = meshParams.ycells_ini;
      vmeshParams.gridLength[2] = meshParams.zcells_ini;
      vmeshParams.blockLength[0] = 1;
      vmeshParams.blockLength[1] = 1;
      vmeshParams.blockLength[2] = 1;
      vmeshParams.refLevelMaxAllowed = 0;
      
      meshID = getObjectWrapper().velocityMeshes.size();
      getObjectWrapper().velocityMeshes.push_back(vmeshParams);
      mesh.initialize(meshID,getObjectWrapper().velocityMeshes);
      initialized = true;
      return initialized;
   }
   
   void MeshDataContainer::reallocate() {
      if (initialized == false) {
         cerr << "(MeshDataContainer) ERROR: Class has not been initialized, exiting at ";
         cerr << __FILE__ << ":" << __LINE__ << endl;
         exit(1);
      }
      
      // Recalculate local IDs
      mesh.setGrid(getLocalCells());

      // Reallocate memory 
      for (size_t i=0; i<meshData.size(); ++i) {
         meshData[i].setMeshSize(mesh.size());
         meshData[i].reallocate();
      }
   }

} // namespace mesh
