/*
 * This file is part of Vlasiator.
 * Copyright 2010-2024 Finnish Meteorological Institute and University of Helsinki
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
#define VELOCITY_MESH_PARAMETERS_H

#include "definitions.h"
#include <array>
#include <iostream>
#include <vector>

#include "arch/arch_device_api.h"
#ifdef USE_GPU
   #include "arch/gpu_base.hpp"
#include "include/splitvector/splitvec.h"
#endif

// One per particle population
#define MAX_VMESH_PARAMETERS_COUNT 32

namespace vmesh {

   /** Wrapper for mesh parameters. The object wrapper reads one or more velocity meshes
    * from the configuration file and stores them to the mesh vector
    * MeshWrapper::velocityMeshes. The particle species store a mesh ID, which is an index
    * to MeshWrapper::velocityMeshes. Many "get" functions in VelocityMesh are
    * wrapper functions, which return the values stored in MeshParameters.
    */
   struct MeshParameters {
      std::string name;                   /**< Name of the mesh (unique).*/
      vmesh::LocalID max_velocity_blocks; /**< Maximum valid block local ID.*/
      Real meshLimits[6];                 /**< Velocity mesh bounding box limits vx_min,vx_max,...,vz_max.*/
      vmesh::LocalID gridLength[3];       /**< Number of blocks in mesh per coordinate at base grid level.*/
      vmesh::LocalID blockLength[3];      /**< Number of phase-space cells per coordinate in block.*/

      // ***** DERIVED PARAMETERS, CALCULATED BY INITVELOCITYMESHES ***** //
      bool initialized;      /**< If true, variables in this struct contain sensible values.*/
      Real meshMinLimits[3]; /**< Minimum coordinate values of the grid bounding box.*/
      Real meshMaxLimits[3]; /**< Maximum coordinate values of the grid bounding box.*/
      Real blockSize[3];     /**< Size of a block at base grid level.*/
      Real cellSize[3];      /**< Size of a cell in a block at base grid level.*/
      Real gridSize[3];      /**< Physical size of the grid bounding box.*/

      MeshParameters() { initialized = false; }
   };

   struct MeshWrapper {
      MeshWrapper() {
         velocityMeshesCreation = new std::vector<vmesh::MeshParameters>(1);
         velocityMeshesCreation->clear();
      }
      ~MeshWrapper() {
         delete velocityMeshes;
         delete velocityMeshesCreation;
      }
      MeshWrapper(const MeshWrapper& other) { velocityMeshesCreation = new std::vector<vmesh::MeshParameters>(*(other.velocityMeshesCreation)); }
      MeshWrapper& operator=(const MeshWrapper& other) {
         delete velocityMeshes;
         delete velocityMeshesCreation;
         velocityMeshesCreation = new std::vector<vmesh::MeshParameters>(*(other.velocityMeshesCreation));
         return *this;
      }
      std::vector<vmesh::MeshParameters>* velocityMeshesCreation;
      // We also need an array so we can copy this data into direct GPU-device memory.
      // On the CPU side we actually reserve enough room for
      // MAX_VMESH_PARAMETERS_COUNT MeshParameters.
      std::array<vmesh::MeshParameters, MAX_VMESH_PARAMETERS_COUNT>* velocityMeshes;
      void initVelocityMeshes(const uint nMeshes); /**< Pre-calculate more helper parameters for velocity meshes. */
      void uploadMeshWrapper();                    /**< Send a copy of the MeshWrapper into GPU memory */
   };

   void allocateMeshWrapper();
   MeshWrapper* host_getMeshWrapper();
   #ifdef USE_GPU
   // To avoid using relocatable code builds, we use static instances of the GPU constant memory for the mesh wrapper.
   // This means that each compilation unit will use its own. To make sure all instances are initialized with the same
   // address, we use the Ctor of a static object to register all intances so that the allocated memory pointer
   // could be copied to all of them.
   __device__ __constant__ MeshWrapper* meshWrapperDevInstance;
   ARCH_DEV static MeshWrapper* gpu_getMeshWrapper() { return meshWrapperDevInstance; };
   // Static object and corresponding instance whose Ctor is used register all the instances of device meshWrapperDev
   // symbols.
   struct meshWrapperDevRegistor {
      meshWrapperDevRegistor(MeshWrapper*&);
   };
   static meshWrapperDevRegistor meshWrapperDevRegistorInstance(meshWrapperDevInstance);

   void deallocateMeshWrapper(); /**< Deallocate GPU memory */
   #endif

   // Caller, inlined into other compilation units, will call either host or device getter
   ARCH_HOSTDEV inline MeshWrapper* getMeshWrapper() {
   #if defined(USE_GPU)
#if (defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__))
      return gpu_getMeshWrapper();
      #else
      return host_getMeshWrapper();
      #endif
#else
      return host_getMeshWrapper();
   #endif
   }

   ARCH_HOSTDEV inline void printVelocityMesh(const uint meshIndex) {
      vmesh::MeshParameters* vMesh = &((*(getMeshWrapper()->velocityMeshes))[meshIndex]);
      printf("\nPrintout of velocity mesh %d \n", meshIndex);
      // printf("Meshwrapper address 0x%lx\n",getMeshWrapper());
      // printf("array of meshes address 0x%lx\n",&(getMeshWrapper()->velocityMeshes));
      // printf("Mesh address 0x%lx\n",vMesh);
      printf("Mesh size\n");
      printf(" %d %d %d \n", vMesh->gridLength[0], vMesh->gridLength[1], vMesh->gridLength[2]);
      printf("Block size\n");
      printf(" %d %d %d \n", vMesh->blockLength[0], vMesh->blockLength[1], vMesh->blockLength[2]);
      printf("Mesh limits \n");
      printf(" %f %f %f %f \n", vMesh->meshMinLimits[0], vMesh->meshLimits[0], vMesh->meshMaxLimits[0], vMesh->meshLimits[1]);
      printf(" %f %f %f %f \n", vMesh->meshMinLimits[1], vMesh->meshLimits[2], vMesh->meshMaxLimits[1], vMesh->meshLimits[3]);
      printf(" %f %f %f %f \n", vMesh->meshMinLimits[2], vMesh->meshLimits[4], vMesh->meshMaxLimits[2], vMesh->meshLimits[5]);
      printf("Derived mesh parameters \n");
      printf(" gridSize %f %f %f \n", vMesh->gridSize[0], vMesh->gridSize[1], vMesh->gridSize[2]);
      printf(" blockSize %f %f %f \n", vMesh->blockSize[0], vMesh->blockSize[1], vMesh->blockSize[2]);
      printf(" cellSize %f %f %f \n", vMesh->cellSize[0], vMesh->cellSize[1], vMesh->cellSize[2]);
      printf(" max velocity blocks %d \n\n", vMesh->max_velocity_blocks);
   }

} // namespace vmesh

#endif  /* VELOCITY_MESH_PARAMETERS_H */
