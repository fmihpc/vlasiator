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
 * File:   velocity_mesh_parameters.h
 * Author: sandroos
 *
 * Created on April 10, 2015, 12:44 PM
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

      // ***** DERIVED PARAMETERS SPECIFIC TO VAMR ***** //
      std::vector<vmesh::GlobalID> offsets;     /**< Block global ID offsets for each refinement level.*/
      std::vector<Real> blockSizes;             /**< Velocity block sizes (dvx,dvy,dvz) for each refinement level.
                                                 * This vector is initialized to size 3*(refLevelMaxAllowed+1)
                                                 * in VelocityMesh::initialize (VAMR mesh).*/
      std::vector<Real> cellSizes;              /**< Velocity block phase-space cell sizes (dvx,dvy,dvz) for each 
                                                 * refinement level. This vector is initialized to size 
                                                 * 3*(refLevelMaxAllowed+1) in VelocityMesh::initialize (VAMR mesh).*/
      std::vector<vmesh::LocalID> gridLengths;  /**< Velocity grid lengths for each refinement level.
                                                 * This vector is initialized to size 3*(refLevelMaxAllowed+1) 
                                                 * in VelocityMesh::initialize (VAMR mesh).*/

      MeshParameters() {
         initialized = false;
      }
   };

} // namespace vmesh

#endif	/* VELOCITY_MESH_PARAMETERS_H */

