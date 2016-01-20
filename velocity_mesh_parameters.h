/* This file is part of Vlasiator.
 * 
 * File:   velocity_mesh_parameters.h
 * Author: sandroos
 *
 * Created on April 10, 2015, 12:44 PM
 * 
 * Copyright 2015 Finnish Meteorological Institute
 */

#ifndef VELOCITY_MESH_PARAMETERS_H
#define	VELOCITY_MESH_PARAMETERS_H

#include <vector>
#include "definitions.h"

namespace vmesh {
   
   /** Wrapper for mesh parameters. The Project class (projects/project.cpp) reads 
    * one or more velocity meshes from the configuration file and stores them to 
    * vector ObjectWrapper::velocityMeshes. The velocity meshes stored in each spatial 
    * cell (one mesh per particle population) store a mesh ID, which is an index 
    * to ObjectWrapper::velocityMeshes. Many "get" functions in VelocityMesh are 
    * wrapper functions, which return the values stored in MeshParameters. This allows 
    * different particle populations to use the same velocity mesh (parameters) if needed.
    * 
    * While the variables have been named 
    * for velocity mesh, this struct can be used for spatial meshes as well.*/
   struct MeshParameters {
      std::string name;                         /**< Name of the mesh (unique).*/
      vmesh::LocalID max_velocity_blocks;       /**< Maximum valid block local ID.*/
      Real meshLimits[6];                       /**< Velocity mesh bounding box limits vx_min,vx_max,...,vz_max.*/
      vmesh::LocalID gridLength[3];             /**< Number of blocks in mesh per coordinate at base grid level.*/
      vmesh::LocalID blockLength[3];            /**< Number of phase-space cells per coordinate in block.*/
      uint8_t refLevelMaxAllowed;               /**< Maximum refinement level allowed, 0=no refinement.*/
      
      // ***** DERIVED PARAMETERS, CALCULATED BY VELOCITY MESH ***** //
      bool initialized;                         /**< If true, variables in this struct contain sensible values.*/
      Real meshMinLimits[3];                    /**< Minimum coordinate values of the grid bounding box.*/
      Real meshMaxLimits[3];                    /**< Maximum coordinate values of the grid bounding box.*/
      Real blockSize[3];                        /**< Size of a block at base grid level.*/
      Real cellSize[3];                         /**< Size of a cell in a block at base grid level.*/
      Real gridSize[3];                         /**< Physical size of the grid bounding box.*/

      // ***** DERIVED PARAMETERS SPECIFIC TO AMR ***** //
      std::vector<vmesh::GlobalID> offsets;     /**< Block global ID offsets for each refinement level.*/
      std::vector<Real> blockSizes;             /**< Velocity block sizes (dvx,dvy,dvz) for each refinement level.
                                                 * This vector is initialized to size 3*(refLevelMaxAllowed+1)
                                                 * in VelocityMesh::initialize (AMR mesh).*/
      std::vector<Real> cellSizes;              /**< Velocity block phase-space cell sizes (dvx,dvy,dvz) for each 
                                                 * refinement level. This vector is initialized to size 
                                                 * 3*(refLevelMaxAllowed+1) in VelocityMesh::initialize (AMR mesh).*/
      std::vector<vmesh::LocalID> gridLengths;  /**< Velocity grid lengths for each refinement level.
                                                 * This vector is initialized to size 3*(refLevelMaxAllowed+1) 
                                                 * in VelocityMesh::initialize (AMR mesh).*/

      MeshParameters() {
         initialized = false;
      }
   };

} // namespace vmesh

#endif	/* VELOCITY_MESH_PARAMETERS_H */

