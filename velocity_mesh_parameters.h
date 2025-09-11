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
#define	VELOCITY_MESH_PARAMETERS_H

#include <vector>
#include <array>
#include "definitions.h"
#include <iostream>

#include "arch/arch_device_api.h"
#ifdef USE_GPU
   #include "include/splitvector/splitvec.h"
   #include "arch/gpu_base.hpp"
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

      // TODO these should be const'd
      const std::string name;                         /**< Name of the mesh (unique).*/
      const std::array<Real, 6> meshLimits;                       /**< Velocity mesh bounding box limits vx_min,vx_max,...,vz_max.*/
      const std::array<uint32_t, 6> hiResRange;             // Min/max x,y,z indices (x_min, x_max, y_min etc.) with double resolution
      const std::array<uint32_t, 3> gridLength;             /**< Number of blocks in mesh per coordinate at base grid level.*/
      const std::array<uint32_t, 3> blockLength;            /**< Number of phase-space cells per coordinate in block.*/

      // ***** DERIVED PARAMETERS, CALCULATED BY INITVELOCITYMESHES ***** //
      const vmesh::LocalID max_velocity_blocks;       /**< Maximum valid block local ID.*/

      // TODO should these be functions instead?
      const std::array<Real, 3> meshMinLimits;                    /**< Minimum coordinate values of the grid bounding box.*/
      const std::array<Real, 3> meshMaxLimits;                    /**< Maximum coordinate values of the grid bounding box.*/
      const std::array<Real, 3> gridSize;                         /**< Physical size of the grid bounding box.*/

      const std::array<Real, 3> blockSize;                        /**< Size of a block at base grid level.*/

      // Based on blocksize
      //const std::array<Real, 3> cellSize;                         /**< Size of a cell in a block at base grid level.*/

      MeshParameters(std::string_view name, std::array<Real, 6> meshLimits, std::array<uint32_t, 6> hiResRange, std::array<uint32_t, 3> gridLength, std::array<uint32_t, 3> blockLength);

      ARCH_HOSTDEV std::array<uint32_t, 3> getIndices(const vmesh::GlobalID& globalID) const {
         if (globalID >= INVALID_GLOBALID) {
            return {INVALID_VEL_BLOCK_INDEX, INVALID_VEL_BLOCK_INDEX, INVALID_VEL_BLOCK_INDEX};
         }

         return {
            globalID % gridLength[0],
            (globalID / gridLength[0]) % gridLength[1],
            globalID / (gridLength[0] * gridLength[1])
         };
      }

      //[[deprecated]]
      ARCH_HOSTDEV Real getBlockDx(int idx) const {
         return blockSize[idx];
      }

      // Assumption: cell size in coordinate i only depends on the grid coordinate x_i
      ARCH_HOSTDEV Real getBlockDx(uint32_t cellIndex, int idx) const {
         return blockSize[idx] * (cellIndex >= hiResRange[2 * idx] && cellIndex < hiResRange[2 * idx + 1] ? 0.5 : 1.0);
      }

      ARCH_HOSTDEV Real getBlockDxFromID(const vmesh::GlobalID globalID, int idx) const {
         return getBlockDx(getIndices(globalID)[idx], idx);
      }

      //[[deprecated]]
      ARCH_HOSTDEV Real getCellDx(int idx) const {
         return blockSize[idx] / blockLength[idx];
      }

      ARCH_HOSTDEV Real getCellDx(uint32_t cellIndex, int idx) const {
         // I _think_ basic division works here
         return getBlockDx(cellIndex / blockLength[idx], idx) / blockLength[idx];
      }

      // This guy should probably check the cells block and then its dx
      ARCH_HOSTDEV Real getCellDxFromID(const vmesh::GlobalID globalID, int idx) const {
         return blockSize[idx] / blockLength[idx];
      }

      ARCH_HOSTDEV bool getBlockSize(const vmesh::GlobalID globalID, Real size[3]) const {
         for (int i = 0; i < 3; ++i) {
            size[i] = getBlockDxFromID(globalID, i);
         }
         return true;
      }

      ARCH_HOSTDEV bool getCellSize(const vmesh::GlobalID globalID, Real size[3]) const {
         for (int i = 0; i < 3; ++i) {
            size[i] = getCellDx(globalID, i);
         }
         return true;
      }
   };

   struct MeshWrapper {
      MeshWrapper() {
         velocityMeshesCreation = new std::vector<vmesh::MeshParameters>;
      }
      ~MeshWrapper() {
         delete velocityMeshes;
         delete velocityMeshesCreation;
      }
      MeshWrapper(const MeshWrapper& other) {
         velocityMeshesCreation = new std::vector<vmesh::MeshParameters>(*(other.velocityMeshesCreation));
      }
      MeshWrapper& operator=(const MeshWrapper& other) {
         delete velocityMeshes;
         delete velocityMeshesCreation;
         velocityMeshesCreation = new std::vector<vmesh::MeshParameters>(*(other.velocityMeshesCreation));
         return *this;
      }

      // TODO bounds checking?
      ARCH_HOSTDEV vmesh::MeshParameters& at(int index) {
         return velocityMeshes[index];
      }

      ARCH_HOSTDEV const vmesh::MeshParameters& at(int index) const {
         return velocityMeshes[index];
      }

      std::vector<vmesh::MeshParameters> *velocityMeshesCreation;
      // We also need an array so we can copy this data into direct GPU-device memory.
      // On the CPU side we actually reserve enough room for
      // MAX_VMESH_PARAMETERS_COUNT MeshParameters.
      vmesh::MeshParameters* velocityMeshes;
      //std::array<vmesh::MeshParameters,MAX_VMESH_PARAMETERS_COUNT> *velocityMeshes;
      void initVelocityMeshes(const uint nMeshes);  /**< Pre-calculate more helper parameters for velocity meshes. */
      void uploadMeshWrapper();   /**< Send a copy of the MeshWrapper into GPU memory */
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
      vmesh::MeshParameters *vMesh = &(getMeshWrapper()->at(meshIndex));
      printf("\nPrintout of velocity mesh %d \n",meshIndex);
      // printf("Meshwrapper address 0x%lx\n",getMeshWrapper());
      // printf("array of meshes address 0x%lx\n",&(getMeshWrapper()->velocityMeshes));
      // printf("Mesh address 0x%lx\n",vMesh);
      printf("Mesh size\n");
      printf(" %d %d %d \n",vMesh->gridLength[0],vMesh->gridLength[1],vMesh->gridLength[2]);
      printf("Block size\n");
      printf(" %d %d %d \n",vMesh->blockLength[0],vMesh->blockLength[1],vMesh->blockLength[2]);
      printf("Mesh limits \n");
      printf(" %f %f %f %f \n",vMesh->meshMinLimits[0],vMesh->meshLimits[0],vMesh->meshMaxLimits[0],vMesh->meshLimits[1]);
      printf(" %f %f %f %f \n",vMesh->meshMinLimits[1],vMesh->meshLimits[2],vMesh->meshMaxLimits[1],vMesh->meshLimits[3]);
      printf(" %f %f %f %f \n",vMesh->meshMinLimits[2],vMesh->meshLimits[4],vMesh->meshMaxLimits[2],vMesh->meshLimits[5]);
      printf("Derived mesh parameters \n");
      printf(" gridSize %f %f %f \n",vMesh->gridSize[0],vMesh->gridSize[1],vMesh->gridSize[2]);
      // This will be nonsense with variable size
      // printf(" blockSize %f %f %f \n",vMesh->blockSize[0],vMesh->blockSize[1],vMesh->blockSize[2]);
      // printf(" cellSize %f %f %f \n",vMesh->cellSize[0],vMesh->cellSize[1],vMesh->cellSize[2]);
      printf(" max velocity blocks %d \n\n",vMesh->max_velocity_blocks);
   }

} // namespace vmesh

#endif	/* VELOCITY_MESH_PARAMETERS_H */
