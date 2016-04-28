/* This file is part of Vlasiator.
 * 
 * File:   mesh_data_container.cpp
 * Author: sandroos
 *
 * Created on June 2, 2015, 11:18 AM
 * 
 * Copyright 2015 Finnish Meteorological Institute
 */

#include "mesh_data_container.h"
#include "object_wrapper.h"

using namespace std;

namespace mesh {

   MeshDataContainer::MeshDataContainer(): initialized(false) { }
   
   bool MeshDataContainer::initialize() {
      vmesh::MeshParameters meshParams;
      meshParams.name = "SpatialGrid";
      meshParams.max_velocity_blocks = numeric_limits<vmesh::LocalID>::max();
      meshParams.meshLimits[0] = Parameters::xmin;
      meshParams.meshLimits[1] = Parameters::xmax;
      meshParams.meshLimits[2] = Parameters::zmin;
      meshParams.meshLimits[3] = Parameters::zmax;
      meshParams.meshLimits[4] = Parameters::zmin;
      meshParams.meshLimits[5] = Parameters::zmax;
      meshParams.gridLength[0] = Parameters::xcells_ini;
      meshParams.gridLength[1] = Parameters::ycells_ini;
      meshParams.gridLength[2] = Parameters::zcells_ini;
      meshParams.blockLength[0] = 1;
      meshParams.blockLength[1] = 1;
      meshParams.blockLength[2] = 1;
      meshParams.refLevelMaxAllowed = 0;
      
      meshID = getObjectWrapper().velocityMeshes.size();
      getObjectWrapper().velocityMeshes.push_back(meshParams);
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
