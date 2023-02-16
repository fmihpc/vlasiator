/*
 * This file is part of Vlasiator.
 * Copyright 2010-2023 Finnish Meteorological Institute & University of Helsinki
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
 * Author: sandroos, mbattarbee
 *
 * Created on April 10, 2015, 12:44 PM
 */

#ifndef VELOCITY_MESH_PARAMETERS_H
#define	VELOCITY_MESH_PARAMETERS_H

#include <vector>
#include "definitions.h"
#include <iostream>

#ifdef __CUDACC__
#define CUDA_HOSTDEV __host__ __device__
#define CUDA_DEV __device__
#else
#define CUDA_HOSTDEV
#define CUDA_DEV
#endif

#ifdef USE_CUDA
   #include "include/splitvector/splitvec.h"
#endif

namespace vmesh {

   /** Wrapper for mesh parameters. The object wrapper reads one or more velocity meshes
    * from the configuration file and stores them to the mesh vector
    * MeshWrapper::velocityMeshes. The particle species store a mesh ID, which is an index
    * to MeshWrapper::velocityMeshes. Many "get" functions in VelocityMesh are
    * wrapper functions, which return the values stored in MeshParameters.
    */
   struct MeshParameters {
      std::string name;                         /**< Name of the mesh (unique).*/
      vmesh::LocalID max_velocity_blocks;       /**< Maximum valid block local ID.*/
      Real meshLimits[6];                       /**< Velocity mesh bounding box limits vx_min,vx_max,...,vz_max.*/
      vmesh::LocalID gridLength[3];             /**< Number of blocks in mesh per coordinate at base grid level.*/
      vmesh::LocalID blockLength[3];            /**< Number of phase-space cells per coordinate in block.*/
      uint8_t refLevelMaxAllowed;               /**< Maximum refinement level allowed, 0=no refinement.*/

      // ***** DERIVED PARAMETERS, CALCULATED BY INITVELOCITYMESHES ***** //
      bool initialized;                         /**< If true, variables in this struct contain sensible values.*/
      Real meshMinLimits[3];                    /**< Minimum coordinate values of the grid bounding box.*/
      Real meshMaxLimits[3];                    /**< Maximum coordinate values of the grid bounding box.*/
      Real blockSize[3];                        /**< Size of a block at base grid level.*/
      Real cellSize[3];                         /**< Size of a cell in a block at base grid level.*/
      Real gridSize[3];                         /**< Physical size of the grid bounding box.*/

#ifdef VAMR
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
#endif
      MeshParameters() {
         initialized = false;
      }
   };

   struct MeshWrapper {
      MeshWrapper() {
#ifdef USE_CUDA
         velocityMeshes = new split::SplitVector<vmesh::MeshParameters>(1);
#else
         velocityMeshes = new std::vector<vmesh::MeshParameters>(1);
#endif
         velocityMeshes->clear();
      }
      //Don't bother with copy constructor or destructors
      //    MeshWrapper(const MeshWrapper& ow);
      //    MeshWrapper& operator=(const MeshWrapper& ow);

   public:
      /**< Parameters for velocity mesh(es) as a vector*/
#ifdef USE_CUDA
      split::SplitVector<vmesh::MeshParameters> *velocityMeshes;
#else
      std::vector<vmesh::MeshParameters> *velocityMeshes;
#endif

      void initVelocityMeshes(const uint nMeshes);  /**< Pre-calculate more helper parameters for velocity meshes. */
      void printVelocityMesh(const uint meshIndex); /**< debug purposes, print contents of mesh. */
      void uploadMeshWrapper();                     /**< Send a copy of the MeshWrapper into GPU memory */
   };

   void allocMeshWrapper();
   MeshWrapper* host_getMeshWrapper();
   CUDA_HOSTDEV MeshWrapper* dev_getMeshWrapper();

   // Caller, inlined into other compilation units, will call either host or device getter
   CUDA_HOSTDEV inline MeshWrapper* getMeshWrapper() {
      #ifndef __CUDA_ARCH__
      return host_getMeshWrapper();
      #else
      return dev_getMeshWrapper();
      #endif
   }

} // namespace vmesh

#endif	/* VELOCITY_MESH_PARAMETERS_H */
