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
 */
/*!
Spatial cell class for Vlasiator that supports a variable number of velocity blocks.
*/

#ifndef VLASIATOR_SPATIAL_CELL_GPU_HPP
#define VLASIATOR_SPATIAL_CELL_GPU_HPP

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <mpi.h>
#include <limits>
#include <stdint.h>
#include <vector>
#include <array>
#include <unordered_map>
#include <set>
#include <map>
#include <phiprof.hpp>
#include <tuple>

#include "../memoryallocation.h"
#include "../common.h"
#include "../parameters.h"
#include "../definitions.h"

#include "velocity_mesh_gpu.h"
#include "velocity_block_container.h"

#ifdef DEBUG_VLASIATOR
   #ifndef DEBUG_SPATIAL_CELL
   #define DEBUG_SPATIAL_CELL
   #endif
#endif

namespace spatial_cell {

   /** GPU mini-kernel for resizing a vmesh on-device */
   __global__ static void resize_vmesh_ondevice_kernel (
      vmesh::VelocityMesh *vmesh,
      vmesh::LocalID nBlocks
      ) {
      vmesh->device_setNewSize(nBlocks);
   }

   /** GPU kernel for scaling a particle population */
   __global__ static void __launch_bounds__(WID3,4) population_scale_kernel (
      vmesh::LocalID nBlocks,
      vmesh::VelocityMesh *vmesh,
      vmesh::VelocityBlockContainer *blockContainer,
      const Real factor
      ) {
      const int blocki = blockIdx.x;
      const int i = threadIdx.x;
      const int j = threadIdx.y;
      const int k = threadIdx.z;
      const uint ti = k*WID2 + j*WID + i;
      // loop over whole velocity space and scale the values
      const uint blockLID = blocki;
      // Pointer to target block data
      Realf* data = blockContainer->getData(blockLID);
      // Scale value
      data[ti] = data[ti] * factor;
   }
   /** GPU kernel for adding a particle population to another with a scaling factor
    This kernel increments existing blocks, creates new ones if the block does
    not exist, and so is not parallel-safe.
   */
   __global__ static void __launch_bounds__(WID3,4) population_increment_kernel (
      vmesh::LocalID nBlocks,
      vmesh::VelocityMesh *vmesh,
      vmesh::VelocityBlockContainer *blockContainer,
      vmesh::VelocityMesh *otherVmesh,
      vmesh::VelocityBlockContainer *otherBlockContainer,
      const Real factor
      // GPUTODO: This could gather into a vector GIDs and (invalidGIDs) of only those GIDs which need to be added
      // and call another kernel to do just that?
      ) {
      //const int gpuBlocks = gridDim.x;
      //const int blocki = blockIdx.x;
      const int i = threadIdx.x;
      const int j = threadIdx.y;
      const int k = threadIdx.z;
      const uint ti = k*WID2 + j*WID + i;
      // if (gpuBlocks != 1) {
      //    if (ti==0 && blocki==0) {
      //       printf("Warning! Calling population_increment_new_kernel from parallel region!\n");
      //    }
      // }
      // for (vmesh::LocalID incLID=blocki; incLID<nBlocks; incLID += gpuBlocks) {
      for (vmesh::LocalID incLID=0; incLID<nBlocks; incLID++) {
         const Realf* fromData = otherBlockContainer->getData(incLID);
         // Global ID of the block containing incoming data
         const vmesh::GlobalID GID = otherVmesh->getGlobalID(incLID);
         // Get local ID of the target block. If the block doesn't exist, create it.
         __shared__ vmesh::LocalID writeLID;
         #ifdef USE_WARPACCESSORS
         vmesh::LocalID toLID = vmesh->warpGetLocalID(GID,ti);
         #else
         vmesh::LocalID toLID = vmesh->getLocalID(GID);
         #endif
         if (toLID == vmesh->invalidLocalID()) {
            #ifdef USE_WARPACCESSORS
            bool created = vmesh->warpPush_back(GID, ti);
            #else
            __shared__ bool created;
            if (ti==0) {
               created = vmesh->push_back(GID);
            }
            __syncthreads();
            #endif
            // Thread zero must create new block
            if (ti==0) {
               if (!created) {
                  assert(0 && "Error in incrementing blockContainer in population_increment_kernel!");
               }
               toLID = blockContainer->push_back();
               Real* parameters = blockContainer->getParameters(toLID);
               vmesh->getBlockInfo(GID, parameters+BlockParams::VXCRD);
               // make new LID available to all threads
               writeLID = toLID;
            }
            __syncthreads();
            Realf* toData = blockContainer->getData(writeLID);
            if (created) {
               // Write values from source cells
               toData[ti] = fromData[ti] * factor;
            }
         } else {
            // Increment with values from source cells
            Realf* toData = blockContainer->getData(toLID);
            toData[ti] += fromData[ti] * factor;
         }
      } // for-loop over velocity blocks
   }
   /** GPU kernel for replacing one particle population with another
   */
   __global__ static void __launch_bounds__(WID3,4) population_replace_kernel (
      vmesh::VelocityMesh *vmesh,
      vmesh::VelocityBlockContainer *blockContainer,
      vmesh::VelocityMesh *otherVmesh,
      vmesh::VelocityBlockContainer *otherBlockContainer
      ) {
      //const int gpuBlocks = gridDim.x; // incoming LID count - use same GID-LID-pairs as incoming
      const vmesh::LocalID LID = blockIdx.x;
      const int i = threadIdx.x;
      const int j = threadIdx.y;
      const int k = threadIdx.z;
      const uint ti = k*WID2 + j*WID + i;
      // Assumes vmesh and VBC vmesh have correct size before kernel is launched.
      // They are set before calling this kernel in setNewsizeClear()
      size_t newSize = otherVmesh->size();
      #ifdef DEBUG_SPATIAL_CELL
      if (ti==0) { // Check sizes
         if (blockContainer->size() != newSize) {
            printf("Incorrect VBC size in population replace kernel!\n");
            blockContainer->setNewSize(newSize);
         }
         if (vmesh->size() != newSize) {
            printf("Incorrect vmesh size in population replace kernel!\n");
            vmesh->device_setNewSize(newSize);
         }
      }
      __syncthreads();
      #endif
      // Global ID of the block containing incoming data
      const vmesh::GlobalID GID = otherVmesh->getGlobalID(LID);
      // Create block in vmesh
      #ifdef USE_WARPACCESSORS
      vmesh->warpPlaceBlock(GID,LID,ti);
      #else
      if (ti==0) {
         vmesh->placeBlock(GID,LID);
      }
      __syncthreads();
      #endif
      // Write values from source cells
      const Realf* fromData = otherBlockContainer->getData(LID);
      Realf* toData = blockContainer->getData(LID);
      toData[ti] = fromData[ti];
      // copy over also blockParameters
      const Real* fromParameters = otherBlockContainer->getParameters(LID);
      Real* toParameters = blockContainer->getParameters(LID);
      if (ti < BlockParams::N_VELOCITY_BLOCK_PARAMS) {
         toParameters[ti] = fromParameters[ti];
      }
   }
/** GPU kernel for clearing the vmesh globalToLocalMap and setting the
    localToGlobalMap size to newSize. If the maps do not have
    sufficient capacity, raises a reallocation flag and returns early.
    This allows operating on-device without page faults.
 */
__global__ static void resize_and_empty_kernel (
   vmesh::VelocityMesh *vmesh,
   vmesh::VelocityBlockContainer *blockContainer,
   vmesh::LocalID newSize
   ) {
   const int ti = threadIdx.x;
   const int blockSize = blockDim.x;
   if (ti==0) {
      // assert checks already happen in the actual setNewSize calls
      // assert(vmesh->capacity() >= newSize && "Insufficient vmesh capacity in resize_and_empty_kernel!");
      // assert(blockContainer->capacity() >= newSize && "Insufficient VBC capacity in resize_and_empty_kernel!");
      vmesh->device_setNewSize(newSize);
      blockContainer->setNewSize(newSize);
   }
   __syncthreads();
   // check map sizepower
   Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID> *globalToLocalMap = vmesh->gpu_expose_map();
   Hashinator::Info *info = globalToLocalMap->expose_mapinfo<false>();
   // Set map fill to zero. Target goal is to have empty map. Even if we need
   // to recapacitate the map, this will lead to a clean new allocation without
   // copying old key-value-pairs over.
   info->fill=0;
   const size_t len = globalToLocalMap->bucket_count();
   const vmesh::GlobalID emptybucket = globalToLocalMap->get_emptybucket();
   Hashinator::hash_pair<vmesh::GlobalID, vmesh::LocalID>* dst = globalToLocalMap->expose_bucketdata<false>();
   for (size_t i = ti; i < len; i+=blockSize) {
      if (dst[i].first != emptybucket) {
         dst[i].first = emptybucket;
      }
   }
}

   /** Wrapper for variables needed for each particle species.
    *  Change order if you know what you are doing.
    * All Real fields should be consecutive, as they are communicated as a block.
    *
    */
   struct Population {
      Real RHO;
      Real V[3];
      Real RHO_R;
      Real V_R[3];
      Real RHO_V;
      Real V_V[3];
      Real P[6];
      Real P_R[6];
      Real P_V[6];
      Real RHOLOSSADJUST = 0.0;      /*!< Counter for particle number loss from the destroying blocks in blockadjustment*/
      Real max_dt[2];                                                /**< Element[0] is max_r_dt, element[1] max_v_dt.*/
      Real velocityBlockMinValue;
      vmesh::LocalID reservation = 0; /* Guidance on vector size reservation */

      uint ACCSUBCYCLES;        /*!< number of subcyles for each cell*/
      vmesh::LocalID N_blocks;                                       /**< Number of velocity blocks, used when receiving velocity
                                                                      * mesh from remote neighbors using MPI.*/
      vmesh::VelocityMesh *vmesh;     /**< Velocity mesh. Contains all velocity blocks that exist
                                      * in this spatial cell. Cells are identified by their unique
                                      * global IDs.*/
      vmesh::VelocityBlockContainer *blockContainer;  /**< Velocity block data.*/
      /* pointers to device copies of vmesh and vbc */
      vmesh::VelocityMesh *dev_vmesh;
      vmesh::VelocityBlockContainer *dev_blockContainer;

      /**< Temporary storage of acceleration transform intersections and sybcycling dt.*/
      Real intersection_z,intersection_z_di,intersection_z_dj,intersection_z_dk;
      Real intersection_x,intersection_x_di,intersection_x_dj,intersection_x_dk;
      Real intersection_y,intersection_y_di,intersection_y_dj,intersection_y_dk;
      Real subcycleDt;

      // Constructor, destructor
      Population() {
         vmesh = new vmesh::VelocityMesh();
         blockContainer = new vmesh::VelocityBlockContainer();
         dev_vmesh = 0;
         dev_blockContainer = 0;
         // Host registers seem to break in multi-gpu per node runs
         // CHK_ERR(gpuHostRegister(&vmesh, sizeof(vmesh::VelocityMesh*),gpuHostRegisterPortable));
         // CHK_ERR(gpuHostRegister(&blockContainer, sizeof(vmesh::VelocityBlockContainer*),gpuHostRegisterPortable));
         // CHK_ERR(gpuHostRegister(vmesh, sizeof(vmesh::VelocityMesh),gpuHostRegisterPortable));
         // CHK_ERR(gpuHostRegister(blockContainer, sizeof(vmesh::VelocityBlockContainer),gpuHostRegisterPortable));
         gpuStream_t stream = gpu_getStream();
         CHK_ERR(gpuMallocAsync((void**)&dev_vmesh, sizeof(vmesh::VelocityMesh),stream));
         CHK_ERR(gpuMallocAsync((void**)&dev_blockContainer, sizeof(vmesh::VelocityBlockContainer),stream));
         CHK_ERR(gpuMemcpyAsync(dev_vmesh, vmesh, sizeof(vmesh::VelocityMesh), gpuMemcpyHostToDevice,stream));
         CHK_ERR(gpuMemcpyAsync(dev_blockContainer, blockContainer, sizeof(vmesh::VelocityBlockContainer), gpuMemcpyHostToDevice,stream));
         // Set values to zero in case of zero-block populations
         RHO = RHO_R = RHO_V = RHOLOSSADJUST = velocityBlockMinValue = ACCSUBCYCLES = N_blocks = 0;
         for (uint i=0; i<2; ++i) {
            max_dt[i] = 0;
         }
         for (uint i=0; i<3; ++i) {
            V[i] = V_R[i] = V_V[i] = 0;
         }
         for (uint i=0; i<6; i++) {
            P[i] = P_R[i] = P_V[i] = 0;
         }
      }
      ~Population() {
         if (dev_vmesh) {
            CHK_ERR(gpuFree(dev_vmesh));
            dev_vmesh=0;
         }
         if (dev_blockContainer) {
            CHK_ERR(gpuFree(dev_blockContainer));
            dev_blockContainer=0;
         }
         delete vmesh;
         delete blockContainer;
      }
      Population(const Population& other) {
         vmesh = new vmesh::VelocityMesh(*(other.vmesh));
         blockContainer = new vmesh::VelocityBlockContainer(*(other.blockContainer));
         dev_vmesh = 0;
         dev_blockContainer = 0;
         // Host registers seem to break in multi-gpu per node runs
         // CHK_ERR(gpuHostRegister(&vmesh, sizeof(vmesh::VelocityMesh*),gpuHostRegisterPortable));
         // CHK_ERR(gpuHostRegister(&blockContainer, sizeof(vmesh::VelocityBlockContainer*),gpuHostRegisterPortable));
         // CHK_ERR(gpuHostRegister(vmesh, sizeof(vmesh::VelocityMesh),gpuHostRegisterPortable));
         // CHK_ERR(gpuHostRegister(blockContainer, sizeof(vmesh::VelocityBlockContainer),gpuHostRegisterPortable));
         gpuStream_t stream = gpu_getStream();
         CHK_ERR(gpuMallocAsync((void**)&dev_vmesh, sizeof(vmesh::VelocityMesh),stream));
         CHK_ERR(gpuMallocAsync((void**)&dev_blockContainer, sizeof(vmesh::VelocityBlockContainer),stream));
         CHK_ERR(gpuMemcpyAsync(dev_vmesh, vmesh, sizeof(vmesh::VelocityMesh), gpuMemcpyHostToDevice,stream));
         CHK_ERR(gpuMemcpyAsync(dev_blockContainer, blockContainer, sizeof(vmesh::VelocityBlockContainer), gpuMemcpyHostToDevice,stream));

         RHO = other.RHO;
         RHO_R = other.RHO_R;
         RHO_V = other.RHO_V;
         RHOLOSSADJUST = other.RHOLOSSADJUST;
         velocityBlockMinValue = other.velocityBlockMinValue;
         ACCSUBCYCLES = other.ACCSUBCYCLES;
         N_blocks = other.N_blocks;
         reservation = other.reservation;
         for (uint i=0; i<2; ++i) {
            max_dt[i] = other.max_dt[i];
         }
         for (uint i=0; i<3; ++i) {
            V[i] = other.V[i];
            V_R[i] = other.V_R[i];
            V_V[i] = other.V_V[i];
         }
         for (uint i=0; i<6; i++) {
            P[i] = other.P[i];
            P_R[i] = other.P_R[i];
            P_V[i] = other.P_V[i];
         }
         vmesh->updateCachedSize();
         blockContainer->updateCachedSize();
         vmesh->updateCachedCapacity();
         blockContainer->updateCachedCapacity();
      }
      const Population& operator=(const Population& other) {
         gpuStream_t stream = gpu_getStream();
         const vmesh::LocalID newSize = other.vmesh->size();
         ResizeClear(newSize); // Updates cached values too

         if (newSize > 0) {
            dim3 block(WID,WID,WID);
            population_replace_kernel<<<newSize, block, 0, stream>>> (
               dev_vmesh,
               dev_blockContainer,
               other.dev_vmesh,
               other.dev_blockContainer
               );
            CHK_ERR( gpuPeekAtLastError() );
         }

         #ifdef DEBUG_SPATIAL_CELL
         vmesh->check();
         #endif

         RHO = other.RHO;
         RHO_R = other.RHO_R;
         RHO_V = other.RHO_V;
         RHOLOSSADJUST = other.RHOLOSSADJUST;
         velocityBlockMinValue = other.velocityBlockMinValue;
         ACCSUBCYCLES = other.ACCSUBCYCLES;
         N_blocks = newSize;
         reservation = other.reservation;
         for (uint i=0; i<2; ++i) {
            max_dt[i] = other.max_dt[i];
         }
         for (uint i=0; i<3; ++i) {
            V[i] = other.V[i];
            V_R[i] = other.V_R[i];
            V_V[i] = other.V_V[i];
         }
         for (uint i=0; i<6; i++) {
            P[i] = other.P[i];
            P_R[i] = other.P_R[i];
            P_V[i] = other.P_V[i];
         }
         return *this;
      }

      void Upload() {
         gpuStream_t stream = gpu_getStream();
         CHK_ERR( gpuMemcpyAsync(dev_vmesh, vmesh, sizeof(vmesh::VelocityMesh), gpuMemcpyHostToDevice, stream) );
         CHK_ERR( gpuMemcpyAsync(dev_blockContainer, blockContainer, sizeof(vmesh::VelocityBlockContainer), gpuMemcpyHostToDevice, stream) );
         //CHK_ERR( gpuStreamSynchronize(stream) );
      }

      void ResizeClear(const uint newSize) {
         // Clears the vmesh globalToLocalMap. Ensures the vmesh localToGlobalMap is of the requested size
         // and that the VBC has the correct size, but does not alter contents of these.
         gpuStream_t stream = gpu_getStream();

         const bool reallocated1 = blockContainer->setNewCapacity(newSize);
         const bool reallocated2 = vmesh->setNewCapacity(newSize);
         // vmesh->print_sizes();
         if (reallocated1 || reallocated2) { // Beware short-circuit evaluation, don't place the recapacitations inside this check!
            Upload();
         }

         // The following kernel tries to resize the vmesh localToGlobalMap,
         // clears the vmesh GlobalToLocalMap, and resizes the velocity block container.
         // Contents of the localToGlobalMap or the VBC are not edited.
         resize_and_empty_kernel<<<1, Hashinator::defaults::MAX_BLOCKSIZE, 0, stream>>> (
            dev_vmesh,
            dev_blockContainer,
            newSize
            );
         CHK_ERR( gpuPeekAtLastError() );
         vmesh->setNewCachedSize(newSize);
         blockContainer->setNewCachedSize(newSize);
         // CHK_ERR( gpuStreamSynchronize(stream) );
      }

      void Scale(creal factor) {
         RHO *= factor;
         RHO_R *= factor;
         RHO_V *= factor;
         for (uint i=0; i<3; ++i) {
            P[i] *= factor;
            P_R[i] *= factor;
            P_V[i] *= factor;
         }
         // Now loop over whole velocity space and scale the values
         vmesh::LocalID nBlocks = vmesh->size();
         gpuStream_t stream = gpu_getStream();
         if (nBlocks > 0) {
            dim3 block(WID,WID,WID);
            population_scale_kernel<<<nBlocks, block, 0, stream>>> (
               nBlocks,
               dev_vmesh,
               dev_blockContainer,
               factor
               );
            CHK_ERR( gpuPeekAtLastError() );
            CHK_ERR( gpuStreamSynchronize(stream) );
         }
      }
      void Increment(const Population& other, creal factor) {
         // Note: moments will be invalidated.
         // Ensure the vmesh and VBC are large enough
         if (factor==0) {
            // Nothing to add
            return;
         }
         gpuStream_t stream = gpu_getStream();
         vmesh::LocalID nBlocks = (other.vmesh)->size();
         vmesh::LocalID nExistingBlocks = vmesh->size();
         // GPUTODO: Find union size first and allocate based on that?
         vmesh->setNewCapacity(nExistingBlocks + nBlocks + 1);
         blockContainer->setNewCapacity(nExistingBlocks + nBlocks + 1);
         Upload(); // GPUTODO: Only upload if new capacity required re-allocation. VBC setNewCapacity
         // already returns bool, vmesh not yet.
         CHK_ERR( gpuStreamSynchronize(stream) );
         // Loop over the whole velocity space, and add scaled values with
         // a kernel. Addition of new blocks is not block-parallel-safe.
         if (nBlocks > 0) {
            dim3 block(WID,WID,WID);
            // Now serial
            population_increment_kernel<<<1, block, 0, stream>>> (
               nBlocks,
               dev_vmesh,
               dev_blockContainer,
               other.dev_vmesh,
               other.dev_blockContainer,
               factor
               );
            CHK_ERR( gpuPeekAtLastError() );
            CHK_ERR( gpuStreamSynchronize(stream) );
         }
         vmesh->updateCachedSize();
         blockContainer->updateCachedSize();
      }

   };

   /** GPU kernel for populating block data and parameters based on list of
       globalIDs and avgs.
   */
   template <typename fileReal> __global__ void add_blocks_from_buffer_kernel (
      const vmesh::VelocityMesh *vmesh,
      vmesh::VelocityBlockContainer *blockContainer,
      const vmesh::LocalID startLID,
      const vmesh::GlobalID* gpuInitBlocks,
      const fileReal* gpuInitBuffer,
      const uint nBlocks
      ) {
      const int blocki = blockIdx.x;
      //const int warpSize = blockDim.x*blockDim.y*blockDim.z;
      const uint ti = threadIdx.z*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;
      Real* parameters = blockContainer->getParameters(startLID);
      Realf *cellBlockData = blockContainer->getData(startLID);
      const uint index = blocki;
      {
         // Copy in cell data, perform conversion float<->double if necessary
         cellBlockData[index*WID3 + ti] = (Realf)gpuInitBuffer[index*WID3 + ti];
         // Set block parameters
         if (ti==0) {
            vmesh::GlobalID GID = gpuInitBlocks[index];
            vmesh->getBlockInfo(GID, parameters + index*BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VXCRD);
         }
         __syncthreads();
      }
   }

   typedef std::array<unsigned int, 3> velocity_cell_indices_t;   /**< Defines the indices of a velocity cell in a velocity block.
                                                                   * Indices start from 0 and the first value is the index in x direction.
                                                                   * Note: these are the (i,j,k) indices of the cell within the block.
                                                                   * Valid values are ([0,WID[,[0,WID[,[0,WID[).*/

   typedef std::array<vmesh::LocalID,3> velocity_block_indices_t; /**< Defines the indices of a velocity block in the velocity grid.
                                                                   * Indices start from 0 and the first value is the index in x direction.
                                                                   * Note: these are the (i,j,k) indices of the block.
                                                                   * Valid values are ([0,vx_length[,[0,vy_length[,[0,vz_length[).*/

   class SpatialCell {
   public:
      SpatialCell();
      ~SpatialCell();
      SpatialCell(const SpatialCell& other);
      const SpatialCell& operator=(const SpatialCell& other);

      void setReservation(const uint popID, const vmesh::LocalID reservationsize, bool force=false);
      vmesh::LocalID getReservation(const uint popID) const;
      void applyReservation(const uint popID);

      vmesh::GlobalID find_velocity_block(vmesh::GlobalID cellIndices[3],const uint popID);
      Realf* get_data(const uint popID);
      const Realf* get_data(const uint popID) const;
      Realf* get_data(const vmesh::LocalID& blockLID,const uint popID);
      const Realf* get_data(const vmesh::LocalID& blockLID,const uint popID) const;
      Real* get_block_parameters(const uint popID);
      const Real* get_block_parameters(const uint popID) const;
      Real* get_block_parameters(const vmesh::LocalID& blockLID,const uint popID);
      const Real* get_block_parameters(const vmesh::LocalID& blockLID,const uint popID) const;

      Realf* dev_get_data(const uint popID);
      const Realf* dev_get_data(const uint popID) const;
      // Realf* dev_get_data(const vmesh::LocalID& blockLID,const uint popID);
      // const Realf* dev_get_data(const vmesh::LocalID& blockLID,const uint popID) const;
      Real* dev_get_block_parameters(const uint popID);
      const Real* dev_get_block_parameters(const uint popID) const;
      // Real* dev_get_block_parameters(const vmesh::LocalID& blockLID,const uint popID);
      // const Real* dev_get_block_parameters(const vmesh::LocalID& blockLID,const uint popID) const;

      Real* get_cell_parameters();
      const Real* get_cell_parameters() const;

      vmesh::LocalID get_number_of_velocity_blocks(const uint popID) const;
      vmesh::LocalID get_number_of_all_velocity_blocks() const;
      int get_number_of_populations() const;
      void debug_population_check(const uint popID) const;
      void debug_population_check(const uint popID, const vmesh::LocalID blockLID) const;

      Population & get_population(const uint popID);
      const Population & get_population(const uint popID) const;
      void set_population(const Population& pop, cuint popID);
      void scale_population(creal factor, cuint popID);
      void increment_population(const Population& pop, creal factor, cuint popID);
      void increment_mass_loss(cuint popID, Real increment);

      std::vector<Population>& get_populations();
      const std::vector<Population>& get_populations() const;

      const Real& get_max_r_dt(const uint popID) const;
      const Real& get_max_v_dt(const uint popID) const;

      const vmesh::LocalID* get_velocity_grid_length(const uint popID);
      const vmesh::GlobalID* get_velocity_grid(const uint popID);
      const Real* get_velocity_grid_block_size(const uint popID);
      const Real* get_velocity_grid_cell_size(const uint popID);
      void get_velocity_block_coordinates(const uint popID,const vmesh::GlobalID& globalID,Real* coords);
      velocity_block_indices_t get_velocity_block_indices(const uint popID,const vmesh::GlobalID globalID);
      vmesh::GlobalID get_velocity_block(const uint popID,vmesh::GlobalID blockIndices[3]) const;
      vmesh::GlobalID get_velocity_block(const uint popID,const velocity_block_indices_t indices) const;
      vmesh::GlobalID get_velocity_block(const uint popID,const Real* coords) const;
      vmesh::GlobalID get_velocity_block(const uint popID,const Real vx,const Real vy,const Real vz) const;
      vmesh::GlobalID get_velocity_block_global_id(const vmesh::LocalID& blockLID,const uint popID) const;
      vmesh::LocalID get_velocity_block_local_id(const vmesh::GlobalID& blockGID,const uint popID) const;
      void get_velocity_block_size(const uint popID,const vmesh::GlobalID block,Real size[3]);
      Real get_velocity_block_vx_min(const uint popID,const vmesh::GlobalID block) const;
      Real get_velocity_block_vx_max(const uint popID,const vmesh::GlobalID block) const;
      Real get_velocity_block_vy_min(const uint popID,const vmesh::GlobalID block) const;
      Real get_velocity_block_vy_max(const uint popID,const vmesh::GlobalID block) const;
      Real get_velocity_block_vz_min(const uint popID,const vmesh::GlobalID block) const;
      Real get_velocity_block_vz_max(const uint popID,const vmesh::GlobalID block) const;

      static unsigned int invalid_block_index();
      static vmesh::GlobalID invalid_global_id();
      static vmesh::LocalID invalid_local_id();

      size_t count(const vmesh::GlobalID& block,const uint popID) const;

      void printMeshSizes();
      static bool setCommunicatedSpecies(const uint popID);

      // Following functions adjust velocity blocks stored on the cell //
      void adjustSingleCellVelocityBlocks(const uint popID, bool doDeleteEmpty=false);
      void adjust_velocity_blocks(const uint popID,
                                  bool doDeleteEmptyBlocks=true);
      vmesh::LocalID adjust_velocity_blocks_caller(const uint popID);
      // Templated function for storing a v-space read from a file or generated elsewhere
      template <typename fileReal> void add_velocity_blocks(const uint popID,const std::vector<vmesh::GlobalID>& blocks,fileReal* initBuffer);

      void update_velocity_block_content_lists(const uint popID);
      bool checkMesh(const uint popID);
      bool checkSizes(const uint popID);
      void clear(const uint popID, bool shrink=false);
      void setNewSizeClear(const uint popID, const vmesh::LocalID& newSize);
      void setNewSizeClear(const uint popID);

      uint64_t get_cell_memory_capacity();
      uint64_t get_cell_memory_size();
      void prepare_to_receive_blocks(const uint popID);
      bool shrink_to_fit();
      size_t size(const uint popID) const;
      void dev_resize_vmesh(const uint popID, const uint nBlocks);
      vmesh::VelocityMesh* get_velocity_mesh(const size_t& popID);
      vmesh::VelocityBlockContainer* get_velocity_blocks(const size_t& popID);
      const vmesh::VelocityBlockContainer* get_velocity_blocks(const size_t& popID) const;
      void dev_upload_population(const uint popID);
      vmesh::VelocityMesh* dev_get_velocity_mesh(const size_t& popID);
      vmesh::VelocityBlockContainer* dev_get_velocity_blocks(const size_t& popID);
      const vmesh::VelocityBlockContainer* dev_get_velocity_blocks(const size_t& popID) const;
      // Prefetches for both blockContainers and vmeshes, all populations
      void prefetchDevice();
      void prefetchHost();

      void set_max_r_dt(const uint popID,const Real& value);
      void set_max_v_dt(const uint popID,const Real& value);

      // Following functions are related to MPI //
      std::tuple<void*, int, MPI_Datatype> get_mpi_datatype(const CellID cellID,const int sender_rank,const int receiver_rank,
                                                            const bool receiving,const int neighborhood);
      static uint64_t get_mpi_transfer_type(void);
      static void set_mpi_transfer_type(const uint64_t type,bool atSysBoundaries=false);
      static void set_mpi_transfer_direction(const int dimension);
      void set_mpi_transfer_enabled(bool transferEnabled);
      void updateSparseMinValue(const uint popID);
      Real getVelocityBlockMinValue(const uint popID) const;

      // Member variables //
      std::array<Real, vderivatives::N_V_DERIVATIVES> derivativesV;           /**< Derivatives of V for vorticity AMR.*/
      std::array<Real, bvolderivatives::N_BVOL_DERIVATIVES> derivativesBVOL;  /**< Derivatives of BVOL needed by the acceleration.
                                                                               * Separate array because it does not need to be communicated.*/
      std::array<Real, CellParams::N_SPATIAL_CELL_PARAMS> parameters;         /**< Bulk variables in this spatial cell.*/
      std::array<Realf, WID3> null_block_data;

      uint64_t ioLocalCellId;                                                 /**< Local cell ID used for IO, not needed elsewhere
                                                                               * and thus not being kept up-to-date.*/
      std::array<Realf*,MAX_NEIGHBORS_PER_DIM> neighbor_block_data;           /**< Pointers for translation operator. We can point to neighbor
                                                                               * cell block data. We do not allocate memory for the pointer.*/
      std::array<vmesh::LocalID,MAX_NEIGHBORS_PER_DIM> neighbor_number_of_blocks;
      std::map<int,std::set<int>> face_neighbor_ranks;
      uint sysBoundaryFlag;                                                   /**< What type of system boundary does the cell belong to.
                                                                               * Enumerated in the sysboundarytype namespace's enum.*/
      uint sysBoundaryLayer;                                                  /**< Layers counted from closest systemBoundary. If 0 then it has not
                                                                               * been computed. First sysboundary layer is layer 1.*/
      int sysBoundaryLayerNew;                                                /** needed (by DCCRG?), do not remove. */
      split::SplitVector<vmesh::GlobalID> *velocity_block_with_content_list=0;      /**< List of existing cells with content (updated by update_has_content()).*/
      split::SplitVector<vmesh::GlobalID> *dev_velocity_block_with_content_list=0;  /**< Device pointer to list of existing cells with content.*/
      vmesh::LocalID velocity_block_with_content_list_size=0;                       /**< Size of vector. Needed for MPI communication of size before actual list transfer.*/
      vmesh::LocalID velocity_block_with_content_list_capacity=0;                   /**< Capacity of vector. Cached value.*/
      Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID> *velocity_block_with_content_map=0, *velocity_block_with_no_content_map=0;
      Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID> *dev_velocity_block_with_content_map=0, *dev_velocity_block_with_no_content_map=0;
      vmesh::LocalID vbwcl_sizePower, vbwncl_sizePower;

      split::SplitVector<vmesh::GlobalID> *list_with_replace_new=0, *dev_list_with_replace_new=0;
      split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>> *list_delete=0, *dev_list_delete=0;
      split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>> *list_to_replace=0, *dev_list_to_replace=0;
      split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>> *list_with_replace_old=0, *dev_list_with_replace_old=0;
      vmesh::LocalID list_with_replace_new_capacity=0, list_delete_capacity=0, list_to_replace_capacity=0, list_with_replace_old_capacity=0;

      Realf* gpu_rhoLossAdjust;
      Real density_pre_adjust, density_post_adjust;
      uint64_t largestvmesh = 0;

      static uint64_t mpi_transfer_type;                                      /**< Which data is transferred by the mpi datatype given by spatial cells.*/
      static bool mpiTransferAtSysBoundaries;                                 /**< Do we only transfer data at boundaries (true), or in the whole system (false).*/

   private:
      static int activePopID;
      bool initialized;
      bool mpiTransferEnabled;

      std::vector<spatial_cell::Population> populations;                        /**< Particle population variables.*/
   };

   inline void SpatialCell::debug_population_check(const uint popID) const {
      #ifdef DEBUG_SPATIAL_CELL
      if (popID >= populations.size()) {
         std::cerr << "ERROR, popID " << popID << " exceeds populations.size() " << populations.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
         exit(1);
      }
      #endif
   }
   inline void SpatialCell::debug_population_check(const uint popID, const vmesh::LocalID blockLID) const {
      debug_population_check(popID);
      #ifdef DEBUG_SPATIAL_CELL
      if (blockLID >= populations[popID].blockContainer->size()) {
         std::cerr << "ERROR, block LID out of bounds, blockContainer->size() " << populations[popID].blockContainer->size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
         exit(1);
      }
      #endif
   }

   inline Realf* SpatialCell::get_data(const uint popID) {
      debug_population_check(popID);
      return populations[popID].blockContainer->getData();
   }

   inline const Realf* SpatialCell::get_data(const uint popID) const {
      debug_population_check(popID);
      return populations[popID].blockContainer->getData();
   }

   inline Realf* SpatialCell::dev_get_data(const uint popID) {
      debug_population_check(popID);
      return populations[popID].dev_blockContainer->getData();
   }

   inline const Realf* SpatialCell::dev_get_data(const uint popID) const {
      debug_population_check(popID);
      return populations[popID].dev_blockContainer->getData();
   }

   inline Realf* SpatialCell::get_data(const vmesh::LocalID& blockLID,const uint popID) {
      debug_population_check(popID,blockLID);
      if (blockLID == vmesh::VelocityMesh::invalidLocalID()) {
         return null_block_data.data();
      }
      return populations[popID].blockContainer->getData(blockLID);
   }

   inline const Realf* SpatialCell::get_data(const vmesh::LocalID& blockLID,const uint popID) const {
      debug_population_check(popID,blockLID);
      if (blockLID == vmesh::VelocityMesh::invalidLocalID()) {
         return null_block_data.data();
      }
      return populations[popID].blockContainer->getData(blockLID);
   }

   inline Real* SpatialCell::get_block_parameters(const uint popID) {
      debug_population_check(popID);
      return populations[popID].blockContainer->getParameters();
   }

   inline const Real* SpatialCell::get_block_parameters(const uint popID) const {
      debug_population_check(popID);
      return populations[popID].blockContainer->getParameters();
   }

   inline void SpatialCell::dev_upload_population(const uint popID) {
      populations[popID].Upload();
   }

   inline Real* SpatialCell::dev_get_block_parameters(const uint popID) {
      debug_population_check(popID);
      return populations[popID].dev_blockContainer->getParameters();
   }

   inline const Real* SpatialCell::dev_get_block_parameters(const uint popID) const {
      debug_population_check(popID);
      return populations[popID].dev_blockContainer->getParameters();
   }

   inline Real* SpatialCell::get_block_parameters(const vmesh::LocalID& blockLID,const uint popID) {
      debug_population_check(popID,blockLID);
      return populations[popID].blockContainer->getParameters(blockLID);
   }

   inline const Real* SpatialCell::get_block_parameters(const vmesh::LocalID& blockLID,const uint popID) const {
      debug_population_check(popID,blockLID);
      return populations[popID].blockContainer->getParameters(blockLID);
   }

   inline Real* SpatialCell::get_cell_parameters() {
      return parameters.data();
   }

   inline const Real* SpatialCell::get_cell_parameters() const {
      return parameters.data();
   }

   inline vmesh::LocalID SpatialCell::get_number_of_velocity_blocks(const uint popID) const {
      debug_population_check(popID);
      //return populations[popID].blockContainer->size();
      // Return size from vmesh instead of VBC to allow use of host-cached value
      return populations[popID].vmesh->size();
   }

    /** Get the total number of velocity blocks in this cell, summed over
     * all existing particle populations.
     * @return Total number of velocity blocks in the cell.*/
    inline vmesh::LocalID SpatialCell::get_number_of_all_velocity_blocks() const {
        vmesh::LocalID N_blocks = 0;
        for (size_t p=0; p<populations.size(); ++p)
            N_blocks += populations[p].blockContainer->size();
        return N_blocks;
    }

   inline int SpatialCell::get_number_of_populations() const {
      return populations.size();
   }

   inline Population & SpatialCell::get_population(const uint popID) {
      return populations[popID];
   }

   inline const Population & SpatialCell::get_population(const uint popID) const {
      return populations[popID];
   }

   inline void SpatialCell::set_population(const Population& pop, cuint popID) {
      // (this->populations[popID].vmesh)->gpu_prefetchDevice();
      // (this->populations[popID].blockContainer)->gpu_prefetchDevice();
      // (pop.vmesh)->gpu_prefetchDevice();
      // (pop.blockContainer)->gpu_prefetchDevice();
      this->populations[popID] = pop;
      // Copy assign includes dev_vmesh upload
   }
   inline void SpatialCell::scale_population(creal factor, cuint popID) {
      // (this->populations[popID].vmesh)->gpu_prefetchDevice();
      // (this->populations[popID].blockContainer)->gpu_prefetchDevice();
      (this->populations[popID]).Scale(factor);
   }
   inline void SpatialCell::increment_population(const Population& pop, creal factor, cuint popID) {
      // (this->populations[popID].vmesh)->gpu_prefetchDevice();
      // (this->populations[popID].blockContainer)->gpu_prefetchDevice();
      // (pop.vmesh)->gpu_prefetchDevice();
      // (pop.blockContainer)->gpu_prefetchDevice();
      (this->populations[popID]).Increment(pop, factor);
   }
   inline void SpatialCell::increment_mass_loss(cuint popID, const Real increment) {
      (this->populations[popID]).RHOLOSSADJUST += increment;
   }

   inline std::vector<Population>& SpatialCell::get_populations() {
      return populations;
   }
   inline const std::vector<Population>& SpatialCell::get_populations() const {
      return populations;
   }

   inline const vmesh::LocalID* SpatialCell::get_velocity_grid_length(const uint popID) {
      return populations[popID].vmesh->getGridLength();
   }

   inline const vmesh::GlobalID* SpatialCell::get_velocity_grid(const uint popID) {
      return (populations[popID].vmesh->getGrid())->data();
   }

   inline const Real* SpatialCell::get_velocity_grid_block_size(const uint popID) {
      return populations[popID].vmesh->getBlockSize();
   }

   inline const Real* SpatialCell::get_velocity_grid_cell_size(const uint popID) {
      return populations[popID].vmesh->getCellSize();
   }

   inline void SpatialCell::get_velocity_block_coordinates(const uint popID,const vmesh::GlobalID& globalID,Real* coords) {
      populations[popID].vmesh->getBlockCoordinates(globalID,coords);
   }

   /*!
    Returns the indices of given velocity block
    */
   inline velocity_block_indices_t SpatialCell::get_velocity_block_indices(const uint popID,const vmesh::GlobalID block) {
      velocity_block_indices_t indices;
      populations[popID].vmesh->getIndices(block,indices[0],indices[1],indices[2]);
      return indices;
   }

   /*!
    Returns the velocity block at given indices or error_velocity_block
    */
   inline vmesh::GlobalID SpatialCell::get_velocity_block(const uint popID,const velocity_block_indices_t indices) const {
      return populations[popID].vmesh->getGlobalID(indices[0],indices[1],indices[2]);
   }

   inline vmesh::GlobalID SpatialCell::get_velocity_block(const uint popID,vmesh::GlobalID blockIndices[3]) const {
      return populations[popID].vmesh->getGlobalID(blockIndices[0],blockIndices[1],blockIndices[2]);
   }

   /*!
    Returns the velocity block at given location or
    error_velocity_block if outside of the velocity grid
    */
   inline vmesh::GlobalID SpatialCell::get_velocity_block(const uint popID,const Real vx,const Real vy,const Real vz) const {
      Real coords[3] = {vx,vy,vz};
      return populations[popID].vmesh->getGlobalID(coords);
   }

   inline vmesh::GlobalID SpatialCell::get_velocity_block(const uint popID,const Real* coords) const {
      return populations[popID].vmesh->getGlobalID(coords);
   }

   inline vmesh::GlobalID SpatialCell::get_velocity_block_global_id(const vmesh::LocalID& blockLID,const uint popID) const {
      debug_population_check(popID);
      return populations[popID].vmesh->getGlobalID(blockLID);
   }

   inline vmesh::LocalID SpatialCell::get_velocity_block_local_id(const vmesh::GlobalID& blockGID,const uint popID) const {
      debug_population_check(popID);
      return populations[popID].vmesh->getLocalID(blockGID);
   }

   inline void SpatialCell::get_velocity_block_size(const uint popID,const vmesh::GlobalID block,Real blockSize[3]) {
      populations[popID].vmesh->getBlockSize(block,blockSize);
   }

   /*!
    Returns the edge where given velocity block starts.
    */
   inline Real SpatialCell::get_velocity_block_vx_min(const uint popID,const vmesh::GlobalID block) const {
      Real coords[3];
      populations[popID].vmesh->getBlockCoordinates(block,coords);
      return coords[0];
   }

   /*!
    Returns the edge where given velocity block ends.
    */
   inline Real SpatialCell::get_velocity_block_vx_max(const uint popID,const vmesh::GlobalID block) const {
      Real coords[3];
      populations[popID].vmesh->getBlockCoordinates(block,coords);

      Real size[3];
      populations[popID].vmesh->getBlockSize(block,size);
      return coords[0]+size[0];
   }

   /*!
    Returns the edge where given velocity block starts.
    */
   inline Real SpatialCell::get_velocity_block_vy_min(const uint popID,const vmesh::GlobalID block) const {
      Real coords[3];
      populations[popID].vmesh->getBlockCoordinates(block,coords);
      return coords[1];
   }

   /*!
    Returns the edge where given velocity block ends.
    */
   inline Real SpatialCell::get_velocity_block_vy_max(const uint popID,const vmesh::GlobalID block) const {
      Real coords[3];
      populations[popID].vmesh->getBlockCoordinates(block,coords);

      Real size[3];
      populations[popID].vmesh->getBlockSize(block,size);
      return coords[1]+size[1];
   }

   /*!
    Returns the edge where given velocity block starts.
    */
   inline Real SpatialCell::get_velocity_block_vz_min(const uint popID,const vmesh::GlobalID block) const {
      Real coords[3];
      populations[popID].vmesh->getBlockCoordinates(block,coords);
      return coords[2];
   }

   /*!
    Returns the edge where given velocity block ends.
    */
   inline Real SpatialCell::get_velocity_block_vz_max(const uint popID,const vmesh::GlobalID block) const {
      Real coords[3];
      populations[popID].vmesh->getBlockCoordinates(block,coords);

      Real size[3];
      populations[popID].vmesh->getBlockSize(block,size);
      return coords[2]+size[2];
   }

   inline unsigned int SpatialCell::invalid_block_index() {
      return vmesh::VelocityMesh::invalidBlockIndex();
   }

   inline vmesh::GlobalID SpatialCell::invalid_global_id() {
      return vmesh::VelocityMesh::invalidGlobalID();
   }

   inline vmesh::GlobalID SpatialCell::invalid_local_id() {
      return vmesh::VelocityMesh::invalidLocalID();
   }

   /*!
    Returns the number of given velocity blocks that exist.
    */
   inline size_t SpatialCell::count(const vmesh::GlobalID& block,const uint popID) const {
      debug_population_check(popID);
      return populations[popID].vmesh->count(block);
   }

   /*!
    Returns the number of existing velocity blocks.
    */
   inline size_t SpatialCell::size(const uint popID) const {
      debug_population_check(popID);
      return populations[popID].vmesh->size();
   }

   /*!
    Resizes the chosen velocity mesh on-device.
    */
   inline void SpatialCell::dev_resize_vmesh(const uint popID, const uint nBlocks) {
      gpuStream_t stream = gpu_getStream();
      const bool reupload = populations[popID].vmesh->setNewCapacity(nBlocks);
      populations[popID].vmesh->setNewCachedSize(nBlocks);
      if (reupload) {
         CHK_ERR( gpuMemcpyAsync(populations[popID].dev_vmesh, populations[popID].vmesh, sizeof(vmesh::VelocityMesh), gpuMemcpyHostToDevice, stream) );
         //CHK_ERR( gpuStreamSynchronize(stream) );
      }
      spatial_cell::resize_vmesh_ondevice_kernel<<<1, 1, 0, stream>>> (
         populations[popID].dev_vmesh,
         nBlocks
         );
      CHK_ERR( gpuPeekAtLastError() );
      // leaving this sync out is a potential cause for issues during MPI communication, but a device-synchronize may do the trick.
      // Hewever, inplementing some device synchronizes in grid.cpp balanceLoad() seems to do the trick.
      //CHK_ERR( gpuStreamSynchronize(stream) );
   }

   inline vmesh::VelocityMesh* SpatialCell::get_velocity_mesh(const size_t& popID) {
      debug_population_check(popID);
      return populations[popID].vmesh;
   }
   inline vmesh::VelocityMesh* SpatialCell::dev_get_velocity_mesh(const size_t& popID) {
      debug_population_check(popID);
      return populations[popID].dev_vmesh;
   }

   inline vmesh::VelocityBlockContainer* SpatialCell::get_velocity_blocks(const size_t& popID) {
      debug_population_check(popID);
      return populations[popID].blockContainer;
   }
   inline const vmesh::VelocityBlockContainer* SpatialCell::get_velocity_blocks(const size_t& popID) const {
      debug_population_check(popID);
      return populations[popID].blockContainer;
   }
   inline vmesh::VelocityBlockContainer* SpatialCell::dev_get_velocity_blocks(const size_t& popID) {
      debug_population_check(popID);
      return populations[popID].dev_blockContainer;
   }
   inline const vmesh::VelocityBlockContainer* SpatialCell::dev_get_velocity_blocks(const size_t& popID) const {
      debug_population_check(popID);
      return populations[popID].dev_blockContainer;
   }

   inline bool SpatialCell::checkMesh(const uint popID) {
      debug_population_check(popID);
      const size_t vmeshSize = (populations[popID].vmesh)->size();
      const size_t vbcSize = (populations[popID].blockContainer)->size();
      if (vmeshSize != vbcSize) {
         printf("checkMesh ERROR: population vmesh %zu and blockcontainer %zu sizes do not match!\n",vmeshSize,vbcSize);
      }
      return populations[popID].vmesh->check();
   }
   inline bool SpatialCell::checkSizes(const uint popID) {
      debug_population_check(popID);
      const size_t vmeshSize = (populations[popID].vmesh)->size();
      const size_t vbcSize = (populations[popID].blockContainer)->size();
      if (vmeshSize != vbcSize) {
         printf("checkSizes ERROR: population vmesh %zu and blockcontainer %zu sizes do not match!\n",vmeshSize,vbcSize);
         return false;
      }
      return true;
   }

   /*!
    Removes all velocity blocks from this spatial cell and frees memory in the cell
    */
   inline void SpatialCell::clear(const uint popID, bool shrink) {
       debug_population_check(popID);
       populations[popID].vmesh->clear(shrink);
       populations[popID].blockContainer->clear(shrink);
    }

   /*!
     Ensures the selected population vmesh localToGlobalMap has sufficient capacity and is of correct size. Does not alter
     its contents. Ensures the selected population vmesh globalToLocalMap has sufficient capacity and is empty.
     Ensures the selected population VBC has sufficient capacity and is of correct size.
   */
   inline void SpatialCell::setNewSizeClear(const uint popID, const vmesh::LocalID& newSize) {
      populations[popID].ResizeClear(newSize);
   }
   inline void SpatialCell::setNewSizeClear(const uint popID) {
      populations[popID].ResizeClear(populations[popID].N_blocks);
   }

   /*!
    Return the memory consumption in bytes as reported using the size()
    functions of the containers in spatial cell
    */
   // GPUTODO Update this for GPU memory as well, same for capacity
   inline uint64_t SpatialCell::get_cell_memory_size() {
      uint64_t size = 0;
      size += WID3 * sizeof(Realf);
      size += velocity_block_with_content_list_size * sizeof(vmesh::GlobalID);
      size += velocity_block_with_content_map->size() * sizeof(Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>);
      size += velocity_block_with_no_content_map->size() * sizeof(Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>);
      size += list_with_replace_new->size() * sizeof(vmesh::GlobalID);
      size += list_delete->size() * sizeof(Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>);
      size += list_to_replace->size() * sizeof(Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>);
      size += list_with_replace_old->size() * sizeof(Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>);
      size += CellParams::N_SPATIAL_CELL_PARAMS * sizeof(Real);
      size += bvolderivatives::N_BVOL_DERIVATIVES * sizeof(Real);

      size += 2 * sizeof(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>);
      size += sizeof(split::SplitVector<vmesh::GlobalID>);
      size += 3 * sizeof(split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>);

      for (size_t popID=0; popID<populations.size(); ++popID) {
          size += populations[popID].vmesh->sizeInBytes();
          size += populations[popID].blockContainer->sizeInBytes();
          size += sizeof(vmesh::VelocityMesh);
          size += sizeof(vmesh::VelocityBlockContainer);
      }
      return size;
   }

   /*!
    Return the memory consumption in bytes as reported using
    the size() functions of the containers in spatial cell
    */
   inline uint64_t SpatialCell::get_cell_memory_capacity() {
      uint64_t capacity = 0;
      capacity += WID3 * sizeof(Realf);
      // capacity += velocity_block_with_content_list->capacity() * sizeof(vmesh::GlobalID);
      // capacity += velocity_block_with_content_map->bucket_count() * sizeof(Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>);
      // capacity += velocity_block_with_no_content_map->bucket_count() * sizeof(Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>);
      // capacity += list_with_replace_new->capacity() * sizeof(vmesh::GlobalID);
      // capacity += list_delete->capacity() * sizeof(Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>);
      // capacity += list_to_replace->capacity() * sizeof(Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>);
      // capacity += list_with_replace_old->capacity() * sizeof(Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>);
      // *** Cached versions:
      capacity += velocity_block_with_content_list_capacity * sizeof(vmesh::GlobalID);
      capacity += pow(2,vbwcl_sizePower) * sizeof(Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>);
      capacity += pow(2,vbwncl_sizePower) * sizeof(Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>);
      capacity += list_with_replace_new_capacity * sizeof(vmesh::GlobalID);
      capacity += list_delete_capacity * sizeof(Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>);
      capacity += list_to_replace_capacity * sizeof(Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>);
      capacity += list_with_replace_old_capacity * sizeof(Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>);
      // *** These two are not in GPU memory so are excluded
      // capacity += CellParams::N_SPATIAL_CELL_PARAMS * sizeof(Real);
      // capacity += bvolderivatives::N_BVOL_DERIVATIVES * sizeof(Real);

      capacity += 2 * sizeof(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>);
      capacity += sizeof(split::SplitVector<vmesh::GlobalID>);
      capacity += 3 * sizeof(split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>);
      //size_t pops1 = 0, pops2=0;
      for (size_t popID=0; popID<populations.size(); ++popID) {
         capacity += populations[popID].vmesh->capacityInBytes();
         capacity += populations[popID].blockContainer->capacityInBytes();
         capacity += sizeof(vmesh::VelocityMesh);
         capacity += sizeof(vmesh::VelocityBlockContainer);
         // pops1 += populations[popID].vmesh->capacityInBytes();
         // pops1 += populations[popID].blockContainer->capacityInBytes();
         // pops1 += sizeof(vmesh::VelocityMesh);
         // pops1 += sizeof(vmesh::VelocityBlockContainer);
         // pops2 += populations[popID].vmesh->sizeInBytes();
         // pops2 += populations[popID].blockContainer->sizeInBytes();
         // pops2 += sizeof(vmesh::VelocityMesh);
         // pops2 += sizeof(vmesh::VelocityBlockContainer);
      }
      // std::cerr<<"ratio pops/total capacity "<<(float)pops1/(float)capacity<<" pops size/capacity "<<(float)pops2/(float)pops1<<std::endl;
      // std::cerr<<"ratio vmesh "<<(float)populations[0].vmesh->sizeInBytes()/(float)populations[0].vmesh->capacityInBytes();
      // std::cerr<<" VBC "<<(float)populations[0].blockContainer->sizeInBytes()/(float)populations[0].blockContainer->capacityInBytes();
      // std::cerr<<" size(kB) vmesh "<<populations[0].vmesh->sizeInBytes()/1024;
      // std::cerr<<" VBC "<<populations[0].blockContainer->sizeInBytes()/1024;
      // std::cerr<<" capacity(kB) vmesh "<<populations[0].vmesh->capacityInBytes()/1024;
      // std::cerr<<" VBC "<<populations[0].blockContainer->capacityInBytes()/1024<<std::endl;
      return capacity;
   }

   /** Adds a vector of velocity blocks to the population, sets the parameters, and fills the data
       with phase-space densities from the provided buffer (which was read from a file).
       This version calls a kernel to perform operations on-device.
   */
   template <typename fileReal> void SpatialCell::add_velocity_blocks(const uint popID,const std::vector<vmesh::GlobalID>& blocks,fileReal* initBuffer) {
      debug_population_check(popID);
      // Add blocks to velocity mesh
      gpuStream_t stream = gpu_getStream();

      if (populations[popID].vmesh->size() != 0) {
         // TODO: make methods safe to add to a non-empty vmesh
         std::cerr << "Error in adding from buffer: Vmesh not empty!" << __FILE__ << ' ' << __LINE__ << std::endl;
         exit(1);
      }

      const uint nBlocks = blocks.size();
      if (nBlocks==0) {
         // Return if empty
         populations[popID].vmesh->setNewCachedSize(0);
         populations[popID].blockContainer->setNewCachedSize(0);
         return;
      }

      populations[popID].vmesh->setNewCapacity(nBlocks*BLOCK_ALLOCATION_PADDING);
      populations[popID].blockContainer->setNewCapacity(nBlocks*BLOCK_ALLOCATION_PADDING);

      const vmesh::LocalID adds = populations[popID].vmesh->push_back(blocks);
      // Verify that we added all requested blocks
      if (adds != nBlocks) {
         std::cerr << "Failed to add blocks" << __FILE__ << ' ' << __LINE__ << std::endl; exit(1);
         return;
      }
      // populations[popID].vmesh->setNewCachedSize(nBlocks); // handled by push_back
      // populations[popID].blockContainer->setNewCachedSize(nBlocks); // handled by push_back

      const vmesh::LocalID startLID = populations[popID].blockContainer->push_back(nBlocks);
      populations[popID].Upload();

      // Copy data to GPU
      fileReal* gpuInitBuffer;
      vmesh::GlobalID* gpuInitBlocks;
      // TODO: re-use per-thread buffers here, ensuring sufficient allocation.
      CHK_ERR( gpuMallocAsync((void**)&gpuInitBuffer,WID3*nBlocks*sizeof(fileReal), stream) );
      CHK_ERR( gpuMemcpyAsync(gpuInitBuffer, initBuffer,
                              WID3*nBlocks*sizeof(fileReal), gpuMemcpyHostToDevice, stream) );
      CHK_ERR( gpuMallocAsync((void**)&gpuInitBlocks,nBlocks*sizeof(vmesh::GlobalID), stream) );
      CHK_ERR( gpuMemcpyAsync(gpuInitBlocks, blocks.data(),
                              nBlocks*sizeof(vmesh::GlobalID), gpuMemcpyHostToDevice, stream) );

      if (nBlocks>0) {
         dim3 block(WID,WID,WID);
         // Third argument specifies the number of bytes in *shared memory* that is
         // dynamically allocated per block for this call in addition to the statically allocated memory.
         CHK_ERR( gpuStreamSynchronize(stream) );
         spatial_cell::add_blocks_from_buffer_kernel<<<nBlocks, block, 0, stream>>> (
            populations[popID].dev_vmesh,
            populations[popID].dev_blockContainer,
            startLID,
            gpuInitBlocks,
            gpuInitBuffer,
            nBlocks
            );
         CHK_ERR( gpuPeekAtLastError() );
      }
      CHK_ERR( gpuStreamSynchronize(stream) );
      CHK_ERR( gpuFreeAsync(gpuInitBuffer,stream) );
      CHK_ERR( gpuFreeAsync(gpuInitBlocks,stream) );

      #ifdef DEBUG_SPATIAL_CELL
      if (populations[popID].vmesh->size() != populations[popID].blockContainer->size()) {
         std::cerr << "size mismatch in " << __FILE__ << ' ' << __LINE__ << std::endl;
         std::cerr << " velocity mesh size "<<populations[popID].vmesh->size();
         std::cerr << " VBC size "<<populations[popID].blockContainer->size();
         std::cerr << " nBlocks "<<nBlocks;
         std::cerr << " adds " << adds << std::endl;
         exit(1);
      }
      #endif
      #ifdef DEBUG_VLASIATOR
      if (!populations[popID].vmesh->check()) {
         std::cerr << "vmesh check error in " << __FILE__ << ' ' << __LINE__ << std::endl;
         std::cerr << " velocity mesh size "<<populations[popID].vmesh->size();
         std::cerr << " VBC size "<<populations[popID].blockContainer->size();
         std::cerr << " nBlocks "<<nBlocks;
         std::cerr << " adds " << adds << std::endl;
         exit(1);
      }
      #endif
   }

   /*!
    Sets the type of data to transfer by mpi_datatype.
    */
   inline void SpatialCell::set_mpi_transfer_type(const uint64_t type,bool atSysBoundaries) {
      SpatialCell::mpi_transfer_type = type;
      SpatialCell::mpiTransferAtSysBoundaries = atSysBoundaries;
   }

   /*!
    Gets the type of data that will be transferred by mpi_datatype.
    */
   inline uint64_t SpatialCell::get_mpi_transfer_type(void) {
      return SpatialCell::mpi_transfer_type;
   }

   /*!
    Set if this cell is transferred/received using MPI in the next communication phase.
    */
   inline void SpatialCell::set_mpi_transfer_enabled(bool transferEnabled) {
      this->mpiTransferEnabled=transferEnabled;
   }

} // namespaces

#endif
