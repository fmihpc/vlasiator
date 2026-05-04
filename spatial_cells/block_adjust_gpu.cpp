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

#include "block_adjust_gpu.hpp"
#include "block_adjust_gpu_kernels.hpp"
#include "../arch/gpu_base.hpp"
#include "../object_wrapper.h"
#include "../velocity_mesh_parameters.h"

namespace spatial_cell {

   /*!\brief spatial_cell::update_velocity_block_content_lists Finds blocks above the sparsity threshold
    *
    * Bulk call over listed cells of spatial grid
    * Prepares the content / no-content velocity block lists
    * for all requested cells, for the requested popID
    */
void update_velocity_block_content_lists(
   dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   const vector<CellID>& cells,
   const uint popID) {

   const uint nCells = cells.size();
   if (nCells == 0) {
      return;
   }
   if (nCells > 65535) {
      std::cerr<<"ERROR: too many cells ("<<nCells<<") passed to GPU batch operations! Please use more GPUs / MPI tasks."<<std::endl;
      abort();
   }

   // Consider mass loss evaluation?
   const bool gatherMass = getObjectWrapper().particleSpecies[popID].sparse_conserve_mass;

   const gpuStream_t baseStream = gpu_getStream();
   // Allocate buffers for GPU operations
   phiprof::Timer mallocTimer {"allocate buffers for content list analysis"};
   gpu_batch_allocate(nCells,0);

   gpuMemoryManager.startSession(0,0);
   SESSION_HOST_ALLOCATE(gpuMemoryManager, vmesh::LocalID, host_nWithContent, nCells * sizeof(vmesh::LocalID));
   SESSION_HOST_ALLOCATE(gpuMemoryManager, Real, host_mass, nCells * sizeof(Real));
   SESSION_ALLOCATE(gpuMemoryManager, vmesh::LocalID, dev_nWithContent, nCells * sizeof(vmesh::LocalID));
   SESSION_ALLOCATE(gpuMemoryManager, Real, dev_mass, nCells * sizeof(Real));

   mallocTimer.stop();

   phiprof::Timer sparsityTimer {"update Sparsity values, apply memory reservations"};
   size_t largestSizePower = 0;
   size_t largestVelMesh = 0;
#pragma omp parallel
   {
      size_t threadLargestVelMesh = 0;
      size_t threadLargestSizePower = 0;
      SpatialCell *SC;
#pragma omp for schedule(dynamic)
      for (uint i=0; i<nCells; ++i) {
         SC = mpiGrid[cells[i]];
         SC->velocity_block_with_content_list_size = 0;
         SC->updateSparseMinValue(popID);

         vmesh::VelocityMesh* vmesh = SC->get_velocity_mesh(popID);
         // Make sure local vectors are large enough
         const size_t mySize = vmesh->size();
         SC->setReservation(popID,mySize);
         SC->applyReservation(popID);

         // might be better to apply reservation *after* clearing maps, but pointers might change.
         // Store values and pointers
         (GET_POINTER(gpuMemoryManager, vmesh::VelocityMesh*, host_vmeshes))[i] = SC->dev_get_velocity_mesh(popID);
         (GET_POINTER(gpuMemoryManager, vmesh::VelocityBlockContainer*, host_VBCs))[i] = SC->dev_get_velocity_blocks(popID);
         (GET_POINTER(gpuMemoryManager, Real, host_minValues))[i] = SC->getVelocityBlockMinValue(popID);
         (GET_POINTER(gpuMemoryManager, SINGLE_ARG(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>*), host_allMaps))[2*i] = SC->dev_velocity_block_with_content_map;
         (GET_POINTER(gpuMemoryManager, SINGLE_ARG(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>*), host_allMaps))[2*i+1] = SC->dev_velocity_block_with_no_content_map;
         (GET_POINTER(gpuMemoryManager, split::SplitVector<vmesh::GlobalID>*, host_vbwcl_vec))[i] = SC->dev_velocity_block_with_content_list;
         
         // Gather largest values
         threadLargestVelMesh = std::max(threadLargestVelMesh, mySize);
         threadLargestSizePower = std::max(threadLargestSizePower, (size_t)SC->vbwcl_sizePower);
         threadLargestSizePower = std::max(threadLargestSizePower, (size_t)SC->vbwncl_sizePower);
      }
#pragma omp critical
      {
         largestVelMesh = std::max(threadLargestVelMesh, largestVelMesh);
         largestSizePower = std::max(threadLargestSizePower, largestSizePower);
      }
   }
   sparsityTimer.stop();

   phiprof::Timer copyTimer {"copy values to device"};
   // Copy pointers and counters over to device
   CHK_ERR( gpuMemcpyAsync(GET_POINTER(gpuMemoryManager, SINGLE_ARG(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>*), dev_allMaps), GET_POINTER(gpuMemoryManager, SINGLE_ARG(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>*), host_allMaps), 2*nCells*sizeof(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>*), gpuMemcpyHostToDevice, baseStream) );
   CHK_ERR( gpuMemcpyAsync(GET_POINTER(gpuMemoryManager, split::SplitVector<vmesh::GlobalID>*, dev_vbwcl_vec), GET_POINTER(gpuMemoryManager, split::SplitVector<vmesh::GlobalID>*, host_vbwcl_vec), nCells*sizeof(split::SplitVector<vmesh::GlobalID>*), gpuMemcpyHostToDevice, baseStream) );
   CHK_ERR( gpuMemcpyAsync(GET_POINTER(gpuMemoryManager, Real, dev_minValues), GET_POINTER(gpuMemoryManager, Real, host_minValues), nCells*sizeof(Real), gpuMemcpyHostToDevice, baseStream) );
   CHK_ERR( gpuMemcpyAsync(GET_POINTER(gpuMemoryManager, vmesh::VelocityMesh*, dev_vmeshes), GET_POINTER(gpuMemoryManager, vmesh::VelocityMesh*, host_vmeshes), nCells*sizeof(vmesh::VelocityMesh*), gpuMemcpyHostToDevice, baseStream) );
   CHK_ERR( gpuMemcpyAsync(GET_POINTER(gpuMemoryManager, vmesh::VelocityBlockContainer*, dev_VBCs), GET_POINTER(gpuMemoryManager, vmesh::VelocityBlockContainer*, host_VBCs), nCells*sizeof(vmesh::VelocityBlockContainer*), gpuMemcpyHostToDevice, baseStream) );
   if (gatherMass) {
      CHK_ERR( gpuMemsetAsync(GET_SESSION_HOST_POINTER(gpuMemoryManager, Real, dev_mass), 0, nCells*sizeof(Real), baseStream) );
   }
   CHK_ERR( gpuStreamSynchronize(baseStream) );
   copyTimer.stop();

   // Batch clear all hash maps
   phiprof::Timer clearTimer {"clear all content maps"};
   clear_maps_caller(nCells,largestSizePower, baseStream);
   CHK_ERR( gpuStreamSynchronize(baseStream) );
   clearTimer.stop();

   // Batch gather GID-LID-pairs into two maps (one with content, one without)
   phiprof::Timer blockKernelTimer {"update content lists kernel"};
   const dim3 grid2(largestVelMesh,nCells,1);
   batch_update_velocity_block_content_lists_kernel<<<grid2, WID3, 0, baseStream>>> (
      GET_POINTER(gpuMemoryManager, vmesh::VelocityMesh*, dev_vmeshes),
      GET_POINTER(gpuMemoryManager, vmesh::VelocityBlockContainer*, dev_VBCs),
      GET_POINTER(gpuMemoryManager, SINGLE_ARG(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>*), dev_allMaps),
      GET_POINTER(gpuMemoryManager, Real, dev_minValues),
      gatherMass, // Also gathers total mass?
      GET_SESSION_HOST_POINTER(gpuMemoryManager, Real, dev_mass)
      );
   CHK_ERR( gpuPeekAtLastError() );
   CHK_ERR( gpuStreamSynchronize(baseStream) );
   blockKernelTimer.stop();

   // Extract all keys from content maps into content list
   phiprof::Timer extractKeysTimer {"extract content keys"};
   auto rule = []
      __device__(const Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID> *map,
                 const Hashinator::hash_pair<vmesh::GlobalID, vmesh::LocalID>& kval,
                 const vmesh::LocalID threshold,
                 const vmesh::LocalID invalidLID,
                 const vmesh::GlobalID invalidGID) -> bool {
                  // This rule does not use the threshold value
                  const vmesh::GlobalID emptybucket = map->get_emptybucket();
                  const vmesh::GlobalID tombstone   = map->get_tombstone();
                  return ( (kval.first != emptybucket) &&( kval.first != tombstone) );
               };
   // Go via launcher due to templating
   extract_GIDs_kernel_launcher<decltype(rule),vmesh::GlobalID,true>(
      GET_POINTER(gpuMemoryManager, SINGLE_ARG(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>*), dev_allMaps), // points to has_content maps
      GET_POINTER(gpuMemoryManager, split::SplitVector<vmesh::GlobalID>*, dev_vbwcl_vec), // content list vectors, output value
      GET_SESSION_POINTER(gpuMemoryManager, vmesh::LocalID, dev_nWithContent), // content list vector sizes, output value
      rule,
      GET_POINTER(gpuMemoryManager, vmesh::VelocityMesh*, dev_vmeshes), // rule_meshes, not used in this call
      GET_POINTER(gpuMemoryManager, SINGLE_ARG(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>*), dev_allMaps), // rule_maps, not used in this call
      GET_POINTER(gpuMemoryManager, split::SplitVector<vmesh::GlobalID>*, dev_vbwcl_vec), // rule_vectors, not used in this call
      nCells,
      baseStream
      );
   CHK_ERR( gpuStreamSynchronize(baseStream) );
   extractKeysTimer.stop();

   // Update host-side size values
   phiprof::Timer blocklistTimer {"update content lists extract"};
   CHK_ERR( gpuMemcpy(GET_SESSION_HOST_POINTER(gpuMemoryManager, vmesh::LocalID, host_nWithContent), GET_SESSION_POINTER(gpuMemoryManager, vmesh::LocalID, dev_nWithContent), nCells*sizeof(vmesh::LocalID), gpuMemcpyDeviceToHost) );
   CHK_ERR( gpuMemcpy(GET_SESSION_HOST_POINTER(gpuMemoryManager, Real, host_mass), GET_SESSION_HOST_POINTER(gpuMemoryManager, Real, dev_mass), nCells*sizeof(Real), gpuMemcpyDeviceToHost) );
   #pragma omp parallel for schedule(static)
   for (uint i=0; i<nCells; ++i) {
      mpiGrid[cells[i]]->velocity_block_with_content_list_size = (GET_SESSION_HOST_POINTER(gpuMemoryManager, vmesh::LocalID, host_nWithContent))[i];
      mpiGrid[cells[i]]->density_pre_adjust = (GET_SESSION_HOST_POINTER(gpuMemoryManager, Real, host_mass))[i]; // Only one counter per cell, but both this and adjustment are done per-pop before moving to next population.
   }
   blocklistTimer.stop();
   gpuMemoryManager.endSession();
}

void adjust_velocity_blocks_in_cells(
   dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   const vector<CellID>& cellsToAdjust,
   const uint popID
   ) {

   int adjustPreId {phiprof::initializeTimer("Adjusting blocks Pre")};
   int adjustId {phiprof::initializeTimer("Adjusting blocks")};
   int cleanupId {phiprof::initializeTimer("Hashmap cleanup")};
   int adjustPostId {phiprof::initializeTimer("Adjusting blocks Post")};
   const gpuStream_t baseStream = gpu_getStream();
   const gpuStream_t priorityStream = gpu_getPriorityStream();
   const uint nCells = cellsToAdjust.size();

   if (nCells > 65535) {
      std::cerr<<"ERROR: too many cells ("<<nCells<<") passed to GPU batch operations! Please use more GPUs / MPI tasks."<<std::endl;
      abort();
   }

   //GPUTODO: make nCells last dimension of grid in dim3(*,*,nCells)?
   // Allocate buffers for GPU operations
   phiprof::Timer mallocTimer {"allocate buffers for content list analysis"};
   gpu_batch_allocate(nCells,0);

   size_t maxNeighbors = 0;
   size_t largestContentList = 0;
   size_t largestContentListNeighbors = 0;
   // Count maximum number of neighbors, largest size of content blocks
#pragma omp parallel
   {
      size_t threadMaxNeighbors = 0;
      size_t threadLargestContentList = 0;
      size_t threadLargestContentListNeighbors = 0;
#pragma omp for schedule(dynamic)
      for (size_t i=0; i<nCells; ++i) {
         CellID cell_id = cellsToAdjust[i];
         SpatialCell* SC = mpiGrid[cell_id];
         if (SC->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
            continue;
         }
         threadLargestContentList = std::max(threadLargestContentList, (size_t)SC->velocity_block_with_content_list_size);
         size_t cellLargestContentListNeighbors = 0;
         std::unordered_set<CellID> uniqueNeighbors;
         const auto* neighbors = mpiGrid.get_neighbors_of(cell_id, Neighborhoods::NEAREST);
         // find only unique neighbor cells
         for ( const auto& [neighbor_id, dir] : *neighbors) {
            cellLargestContentListNeighbors = std::max(cellLargestContentListNeighbors, (size_t)(mpiGrid[neighbor_id]->velocity_block_with_content_list_size));
            if (neighbor_id != cell_id) {
               uniqueNeighbors.insert(neighbor_id);
            }
         }
         size_t reservationSize = SC->getReservation(popID);
         reservationSize = std::max(cellLargestContentListNeighbors, reservationSize);
         SC->setReservation(popID,reservationSize);
         SC->applyReservation(popID);
         size_t nNeighbors = uniqueNeighbors.size();
         threadMaxNeighbors = std::max(threadMaxNeighbors, nNeighbors);
         threadLargestContentListNeighbors = std::max(threadLargestContentListNeighbors, cellLargestContentListNeighbors);
      }
#pragma omp critical
      {
         maxNeighbors = std::max(maxNeighbors, threadMaxNeighbors);
         largestContentList = std::max(threadLargestContentList, largestContentList);
         largestContentListNeighbors = std::max(threadLargestContentListNeighbors, largestContentListNeighbors);
      }
   } // end parallel region

   // Early return if empty region for this population (
   // GPUTODO FIX: BREAKS VLASOV SUBSTEPPING
   // if (largestContentList==largestContentListNeighbors==0) {
   //    return;
   // }
   gpu_batch_allocate(nCells,maxNeighbors);

   gpuMemoryManager.startSession(0,0);
   SESSION_HOST_ALLOCATE(gpuMemoryManager, split::SplitVector<vmesh::GlobalID>*, host_vbwcl_neigh, maxNeighbors * nCells * sizeof(split::SplitVector<vmesh::GlobalID>*));
   SESSION_ALLOCATE(gpuMemoryManager, split::SplitVector<vmesh::GlobalID>*, dev_vbwcl_neigh, maxNeighbors * nCells * sizeof(split::SplitVector<vmesh::GlobalID>*));

   mallocTimer.stop();

   size_t largestVelMesh = 0;
#pragma omp parallel
   {
      phiprof::Timer timer {adjustPreId};
      size_t threadLargestVelMesh = 0;
#pragma omp for schedule(dynamic)
      for (size_t i=0; i<nCells; ++i) {
         CellID cell_id=cellsToAdjust[i];
         SpatialCell* SC = mpiGrid[cell_id];
         if (SC->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
            (GET_POINTER(gpuMemoryManager, vmesh::VelocityMesh*, host_vmeshes))[i]=0;
            (GET_POINTER(gpuMemoryManager, SINGLE_ARG(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>*), host_allMaps))[2*i]=0;
            (GET_POINTER(gpuMemoryManager, SINGLE_ARG(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>*), host_allMaps))[2*i+1]=0;
            (GET_POINTER(gpuMemoryManager, split::SplitVector<vmesh::GlobalID>*, host_vbwcl_vec))[i]=0;
            (GET_POINTER(gpuMemoryManager, split::SplitVector<vmesh::GlobalID>*, host_lists_with_replace_new))[i]=0;
            continue;
         }

         // Gather largest mesh size for launch parameters
         vmesh::VelocityMesh* vmesh = SC->get_velocity_mesh(popID);
         threadLargestVelMesh = std::max(threadLargestVelMesh, vmesh->size());

         // gather vector with pointers to spatial neighbor lists
         const auto* neighbors = mpiGrid.get_neighbors_of(cell_id, Neighborhoods::NEAREST);
         // Note: at AMR refinement boundaries this can cause blocks to propagate further
         // than absolutely required. Face neighbors, however, are not enough as we must
         // account for diagonal propagation.

         // find only unique neighbor cells
         std::unordered_set<CellID> uniqueNeighbors;
         for ( const auto& [neighbor_id, dir] : *neighbors) {
            if (neighbor_id != cell_id) {
               uniqueNeighbors.insert(neighbor_id);
            }
         }
         std::vector<CellID> reducedNeighbors;
         reducedNeighbors.insert(reducedNeighbors.end(), uniqueNeighbors.begin(), uniqueNeighbors.end());
         const uint nNeighbors = reducedNeighbors.size();
         for (uint iN = 0; iN < maxNeighbors; ++iN) {
            if (iN >= nNeighbors) {
               (GET_SESSION_HOST_POINTER(gpuMemoryManager, split::SplitVector<vmesh::GlobalID>*, host_vbwcl_neigh))[i*maxNeighbors + iN] = 0; // no neighbor at this index
               continue;
            }
            CellID neighbor_id = reducedNeighbors.at(iN);
            // store pointer to neighbor content list
            SpatialCell* NC = mpiGrid[neighbor_id];
            if (NC->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
               (GET_SESSION_HOST_POINTER(gpuMemoryManager, split::SplitVector<vmesh::GlobalID>*, host_vbwcl_neigh))[i*maxNeighbors + iN] = 0;
            } else {
               (GET_SESSION_HOST_POINTER(gpuMemoryManager, split::SplitVector<vmesh::GlobalID>*, host_vbwcl_neigh))[i*maxNeighbors + iN] = mpiGrid[neighbor_id]->dev_velocity_block_with_content_list;
            }
         }

         // Store values and pointers
         (GET_POINTER(gpuMemoryManager, vmesh::VelocityMesh*, host_vmeshes))[i] = SC->dev_get_velocity_mesh(popID);
         (GET_POINTER(gpuMemoryManager, vmesh::VelocityBlockContainer*, host_VBCs))[i] = SC->dev_get_velocity_blocks(popID);
         (GET_POINTER(gpuMemoryManager, SINGLE_ARG(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>*), host_allMaps))[2*i] = SC->dev_velocity_block_with_content_map;
         (GET_POINTER(gpuMemoryManager, SINGLE_ARG(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>*), host_allMaps))[2*i+1] = SC->dev_velocity_block_with_no_content_map;
         (GET_POINTER(gpuMemoryManager, split::SplitVector<vmesh::GlobalID>*, host_vbwcl_vec))[i] = SC->dev_velocity_block_with_content_list;
         (GET_POINTER(gpuMemoryManager, split::SplitVector<vmesh::GlobalID>*, host_lists_with_replace_new))[i] = SC->dev_list_with_replace_new;
         (GET_POINTER(gpuMemoryManager, SINGLE_ARG(split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>*), host_lists_delete))[i] = SC->dev_list_delete;
         (GET_POINTER(gpuMemoryManager, SINGLE_ARG(split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>*), host_lists_to_replace))[i] = SC->dev_list_to_replace;
         (GET_POINTER(gpuMemoryManager, SINGLE_ARG(split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>*), host_lists_with_replace_old))[i] = SC->dev_list_with_replace_old;
      }
      timer.stop();
#pragma omp critical
      {
         largestVelMesh = std::max(threadLargestVelMesh, largestVelMesh);
      }
   } // end parallel region

   /*
    * Perform block adjustment via batch operations
    * */
   phiprof::Timer copyTimer {"copy values to device"};
   // Copy pointers and counters over to device
   CHK_ERR( gpuMemcpyAsync(GET_POINTER(gpuMemoryManager, SINGLE_ARG(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>*), dev_allMaps), GET_POINTER(gpuMemoryManager, SINGLE_ARG(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>*), host_allMaps), 2*nCells*sizeof(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>*), gpuMemcpyHostToDevice, baseStream) );
   CHK_ERR( gpuMemcpyAsync(GET_POINTER(gpuMemoryManager, split::SplitVector<vmesh::GlobalID>*, dev_vbwcl_vec), GET_POINTER(gpuMemoryManager, split::SplitVector<vmesh::GlobalID>*, host_vbwcl_vec), nCells*sizeof(split::SplitVector<vmesh::GlobalID>*), gpuMemcpyHostToDevice, baseStream) );
   CHK_ERR( gpuMemcpyAsync(GET_POINTER(gpuMemoryManager, vmesh::VelocityMesh*, dev_vmeshes), GET_POINTER(gpuMemoryManager, vmesh::VelocityMesh*, host_vmeshes), nCells*sizeof(vmesh::VelocityMesh*), gpuMemcpyHostToDevice, baseStream) );
   if (maxNeighbors>0) {
      CHK_ERR( gpuMemcpyAsync(GET_SESSION_POINTER(gpuMemoryManager, split::SplitVector<vmesh::GlobalID>*, dev_vbwcl_neigh), GET_SESSION_HOST_POINTER(gpuMemoryManager, split::SplitVector<vmesh::GlobalID>*, host_vbwcl_neigh), nCells*maxNeighbors*sizeof(split::SplitVector<vmesh::GlobalID>*), gpuMemcpyHostToDevice, baseStream) );
   }
   CHK_ERR( gpuMemsetAsync(GET_POINTER(gpuMemoryManager, vmesh::LocalID, dev_nBefore), 0, nCells*sizeof(vmesh::LocalID), baseStream) );
   CHK_ERR( gpuMemsetAsync(GET_POINTER(gpuMemoryManager, vmesh::LocalID, dev_nAfter), 0, nCells*sizeof(vmesh::LocalID), baseStream) );
   CHK_ERR( gpuMemsetAsync(GET_POINTER(gpuMemoryManager, vmesh::LocalID, dev_nBlocksToChange), 0, nCells*sizeof(vmesh::LocalID), baseStream) );
   CHK_ERR( gpuMemsetAsync(GET_POINTER(gpuMemoryManager, vmesh::LocalID, dev_resizeSuccess), 0, nCells*sizeof(vmesh::LocalID), baseStream) );
   CHK_ERR( gpuMemsetAsync(GET_POINTER(gpuMemoryManager, vmesh::LocalID, dev_overflownElements), 0, nCells*sizeof(vmesh::LocalID), baseStream) );
   CHK_ERR( gpuMemcpyAsync(GET_POINTER(gpuMemoryManager, split::SplitVector<vmesh::GlobalID>*, dev_lists_with_replace_new), GET_POINTER(gpuMemoryManager, split::SplitVector<vmesh::GlobalID>*, host_lists_with_replace_new), nCells*sizeof(split::SplitVector<vmesh::GlobalID>*), gpuMemcpyHostToDevice, baseStream) );
   CHK_ERR( gpuMemcpyAsync(GET_POINTER(gpuMemoryManager, SINGLE_ARG(split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>*), dev_lists_delete), GET_POINTER(gpuMemoryManager, SINGLE_ARG(split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>*), host_lists_delete), nCells*sizeof(split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>*), gpuMemcpyHostToDevice, baseStream) );
   CHK_ERR( gpuMemcpyAsync(GET_POINTER(gpuMemoryManager, SINGLE_ARG(split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>*), dev_lists_to_replace), GET_POINTER(gpuMemoryManager, SINGLE_ARG(split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>*), host_lists_to_replace), nCells*sizeof(split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>*), gpuMemcpyHostToDevice, baseStream) );
   CHK_ERR( gpuMemcpyAsync(GET_POINTER(gpuMemoryManager, SINGLE_ARG(split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>*), dev_lists_with_replace_old), GET_POINTER(gpuMemoryManager, SINGLE_ARG(split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>*), host_lists_with_replace_old), nCells*sizeof(split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>*), gpuMemcpyHostToDevice, baseStream) );
   CHK_ERR( gpuMemcpyAsync(GET_POINTER(gpuMemoryManager, vmesh::VelocityBlockContainer*, dev_VBCs), GET_POINTER(gpuMemoryManager, vmesh::VelocityBlockContainer*, host_VBCs), nCells*sizeof(vmesh::VelocityBlockContainer*), gpuMemcpyHostToDevice, baseStream) );
   CHK_ERR( gpuStreamSynchronize(baseStream) );
   copyTimer.stop();

   // Note: Velocity halo and spatial neighbor halo can both be evaluated simultaneously.
   // Thus, we launch one into the prioritystream, the other into baseStream.

   // Evaluate velocity halo for local content blocks
   phiprof::Timer blockHaloTimer {"Block halo batch kernels"};
   const int addWidthV = getObjectWrapper().particleSpecies[popID].sparseBlockAddWidthV;
   if (addWidthV!=1) {
      std::cerr<<"Error! "<<__FILE__<<":"<<__LINE__<<" Halo extent is not 1, unsupported size."<<std::endl;
      abort();
   }
   // Halo of 1 in each direction adds up to 26 velocity neighbors.

   if (largestContentList > 0) {
      #ifdef USE_BATCH_WARPACCESSORS
      // For NVIDIA/CUDA, we can do 26 neighbors and 32 threads per warp in a single block.
      // For AMD/HIP, we can do 13 neighbors and 64 threads per warp in a single block, meaning two loops per cell.
      // In either case, we launch blocks equal to largest found velocity_block_with_content_list_size, which was stored
      // into largestContentList
      dim3 grid_vel_halo(largestContentList,nCells,1);
      batch_update_velocity_halo_kernel<<<grid_vel_halo, 26*32, 0, priorityStream>>> (
         GET_POINTER(gpuMemoryManager, vmesh::VelocityMesh*, dev_vmeshes),
         GET_POINTER(gpuMemoryManager, split::SplitVector<vmesh::GlobalID>*, dev_vbwcl_vec),
         GET_POINTER(gpuMemoryManager, SINGLE_ARG(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>*), dev_allMaps) // Needs both content and no content maps
         );
      CHK_ERR( gpuPeekAtLastError() );
      #else
      const uint warpsPerBlockBatchHalo = (threadsPerMP/GPUTHREADS + blocksPerMP - 1)/blocksPerMP;
      dim3 grid_vel_halo((largestContentList + warpsPerBlockBatchHalo - 1)/warpsPerBlockBatchHalo,nCells,1);
      dim3 block_vel_halo(GPUTHREADS, warpsPerBlockBatchHalo, 1);
      // We do 26 (launch with GPUTHREADS) neighbors in a single block at a time.
      batch_update_velocity_halo_kernel<<<grid_vel_halo, block_vel_halo, 0, priorityStream>>> (
         GET_POINTER(gpuMemoryManager, vmesh::VelocityMesh*, dev_vmeshes),
         GET_POINTER(gpuMemoryManager, split::SplitVector<vmesh::GlobalID>*, dev_vbwcl_vec),
         GET_POINTER(gpuMemoryManager, SINGLE_ARG(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>*), dev_allMaps), // Needs both content and no content maps
         warpsPerBlockBatchHalo
         );
      CHK_ERR( gpuPeekAtLastError() );
      #endif
      // CHK_ERR( gpuStreamSynchronize(priorityStream) );
   }

   if (maxNeighbors>0 && largestContentListNeighbors>0) {
      // largestContentListNeighbors accounts for remote (ghost neighbor) content list sizes as well
      #ifdef USE_BATCH_WARPACCESSORS
      // ceil int division
      const size_t blocksNeeded_neigh = 1 + ((largestContentListNeighbors - 1) / (WARPSPERBLOCK));
      dim3 grid_neigh_halo(blocksNeeded_neigh,nCells,maxNeighbors);
      // For NVIDIA/CUDA, we can do 32 neighbor GIDs and 32 threads per warp in a single block.
      // For AMD/HIP, we can do 16 neighbor GIDs and 64 threads per warp in a single block
      // This is handled in-kernel.
      batch_update_neighbour_halo_kernel<<<grid_neigh_halo, WARPSPERBLOCK*GPUTHREADS, 0, baseStream>>> (
         GET_POINTER(gpuMemoryManager, vmesh::VelocityMesh*, dev_vmeshes),
         GET_POINTER(gpuMemoryManager, SINGLE_ARG(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>*), dev_allMaps), // Needs both has_content and has_no_content maps
         GET_SESSION_POINTER(gpuMemoryManager, split::SplitVector<vmesh::GlobalID>*, dev_vbwcl_neigh)
         );
      CHK_ERR( gpuPeekAtLastError() );
      #else
      // Try smaller launch for more spatial cell -parallelism
      const size_t blocksNeeded_neigh = 1 + ((largestContentListNeighbors - 1) / (WARPSPERBLOCK*GPUTHREADS));
      dim3 grid_neigh_halo(blocksNeeded_neigh,nCells,maxNeighbors);
      // Each threads manages a single GID from the neighbour at hand
      batch_update_neighbour_halo_kernel<<<grid_neigh_halo, WARPSPERBLOCK*GPUTHREADS, 0, baseStream>>> (
         GET_POINTER(gpuMemoryManager, vmesh::VelocityMesh*, dev_vmeshes),
         GET_POINTER(gpuMemoryManager, SINGLE_ARG(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>*), dev_allMaps), // Needs both has_content and has_no_content maps
         GET_SESSION_POINTER(gpuMemoryManager, split::SplitVector<vmesh::GlobalID>*, dev_vbwcl_neigh)
         );
      CHK_ERR( gpuPeekAtLastError() );
      #endif
   }
   // Sync both streams
   CHK_ERR( gpuStreamSynchronize(priorityStream) );
   CHK_ERR( gpuStreamSynchronize(baseStream) );
   //CHK_ERR( gpuDeviceSynchronize() );
   blockHaloTimer.stop();
   gpuMemoryManager.endSession();

   // Ensure vectors in dev_lists_with_replace_new have sufficient capacity for has_content_maps.
   // Launch kernel which accesses the vector capacities with the map sizes and stores the required capacity
   // for vectors in a buffer (or 0 to indicate no need to recapacitate). After that, copy that buffer to host,
   // go through it, recapcitate as necessary, and if any recapacitiations happened, update the
   // dev_lists_with_replace_new buffer with new vector addresses and upload it to device again.
   // (this re-uploading is probably not needed, would need verifying that splitvector device handles don't get
   // reallocated)
   check_vector_capacities<<<nCells,1,0,baseStream>>>(
      GET_POINTER(gpuMemoryManager, SINGLE_ARG(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>*), dev_allMaps),
      GET_POINTER(gpuMemoryManager, split::SplitVector<vmesh::GlobalID>*, dev_lists_with_replace_new),
      GET_POINTER(gpuMemoryManager, vmesh::LocalID, dev_overflownElements)
      );
   CHK_ERR( gpuPeekAtLastError() );
   CHK_ERR( gpuMemcpyAsync(GET_POINTER(gpuMemoryManager, vmesh::LocalID, host_overflownElements), GET_POINTER(gpuMemoryManager, vmesh::LocalID, dev_overflownElements), nCells*sizeof(vmesh::LocalID), gpuMemcpyDeviceToHost, baseStream) );
   CHK_ERR( gpuStreamSynchronize(baseStream) );
   bool reUpload = false;
   for (size_t i=0; i<nCells; ++i) {
      if ((GET_POINTER(gpuMemoryManager, vmesh::LocalID, host_overflownElements))[i] != 0) {
         reUpload = true;
         CellID cell_id = cellsToAdjust[i];
         SpatialCell* SC = mpiGrid[cell_id];
         SC->setReservation(popID,(GET_POINTER(gpuMemoryManager, vmesh::LocalID, host_overflownElements))[i]*BLOCK_ALLOCATION_PADDING);
         SC->applyReservation(popID);
         (GET_POINTER(gpuMemoryManager, split::SplitVector<vmesh::GlobalID>*, host_lists_with_replace_new))[i] = SC->dev_list_with_replace_new;
      }
   }
   CHK_ERR( gpuDeviceSynchronize() );
   CHK_ERR( gpuMemsetAsync(GET_POINTER(gpuMemoryManager, vmesh::LocalID, dev_overflownElements), 0, nCells*sizeof(vmesh::LocalID), baseStream) );
   if (reUpload) {
      CHK_ERR( gpuMemcpyAsync(GET_POINTER(gpuMemoryManager, split::SplitVector<vmesh::GlobalID>*, dev_lists_with_replace_new), GET_POINTER(gpuMemoryManager, split::SplitVector<vmesh::GlobalID>*, host_lists_with_replace_new), nCells*sizeof(split::SplitVector<vmesh::GlobalID>*), gpuMemcpyHostToDevice, baseStream) );
      CHK_ERR( gpuStreamSynchronize(baseStream) );
   }

   /**
       Now extract vectors to be used in actual block adjustment.
       Previous kernels may have added dummy (to be added) entries to
       velocity_block_with_content_map with LID=vmesh->invalidLocalID()
       Or non-content blocks which should be retained (with correct LIDs).

       Note: These batch operations always include deletion of not-needed blocks.

       Rules used in extracting keys or elements from hashmaps:
       These are provided with the value of nBlocksAfterAdjust as the
       threshold argument by the kernel. To this end, rule_meshes, rule_maps,
       and rule_vectors pointer buffers are provided to the kernels.
   */
   phiprof::Timer extractKeysTimer {"extract content keys"};
   // Go via caller, then launcher due to templating. Templating manages rule lambda type,
   // output vector type, as well as a flag whether the output vector should take the whole
   // element from the map, or just the first of the pair.

   // Finds new Blocks (GID,LID) needing to be added
   // Note:list_with_replace_new then contains both new GIDs to use for replacements and new GIDs to place at end of vmesh
   extract_to_add_caller(
      GET_POINTER(gpuMemoryManager, SINGLE_ARG(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>*), dev_allMaps), // input maps: this is has_content_maps
      GET_POINTER(gpuMemoryManager, split::SplitVector<vmesh::GlobalID>*, dev_lists_with_replace_new), // output vecs
      NULL, // pass null to not store vector lengths
      GET_POINTER(gpuMemoryManager, vmesh::VelocityMesh*, dev_vmeshes), // rule_meshes, not used in this call
      GET_POINTER(gpuMemoryManager, SINGLE_ARG(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>*), dev_allMaps)+1, // rule_maps, not used in this call
      GET_POINTER(gpuMemoryManager, split::SplitVector<vmesh::GlobalID>*, dev_vbwcl_vec), // rule_vectors, not used in this call
      nCells,
      baseStream
      ); // This needs to complete before the next 3 extractions
   // Finds Blocks (GID,LID) to be rescued from end of v-space
   extract_to_delete_or_move_caller(
      GET_POINTER(gpuMemoryManager, SINGLE_ARG(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>*), dev_allMaps), // input maps: this is has_content_maps
      GET_POINTER(gpuMemoryManager, SINGLE_ARG(split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>*), dev_lists_with_replace_old), // output vecs
      NULL, // pass null to not store vector lengths
      GET_POINTER(gpuMemoryManager, vmesh::VelocityMesh*, dev_vmeshes), // rule_meshes
      GET_POINTER(gpuMemoryManager, SINGLE_ARG(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>*), dev_allMaps)+1, // rule_maps: this is has_no_content_maps
      GET_POINTER(gpuMemoryManager, split::SplitVector<vmesh::GlobalID>*, dev_lists_with_replace_new), // rule_vectors
      nCells,
      baseStream
      );
   // Find Blocks (GID,LID) to be outright deleted
   extract_to_delete_or_move_caller(
      GET_POINTER(gpuMemoryManager, SINGLE_ARG(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>*), dev_allMaps)+1, // input maps: this is has_no_content_maps
      GET_POINTER(gpuMemoryManager, SINGLE_ARG(split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>*), dev_lists_delete), // output vecs
      NULL, // pass null to not store vector lengths
      GET_POINTER(gpuMemoryManager, vmesh::VelocityMesh*, dev_vmeshes), // rule_meshes
      GET_POINTER(gpuMemoryManager, SINGLE_ARG(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>*), dev_allMaps)+1, // rule_maps: this is has_no_content_maps
      GET_POINTER(gpuMemoryManager, split::SplitVector<vmesh::GlobalID>*, dev_lists_with_replace_new), // rule_vectors
      nCells,
      baseStream
      );
   // Find Blocks (GID,LID) to be replaced with new ones
   extract_to_replace_caller(
      GET_POINTER(gpuMemoryManager, SINGLE_ARG(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>*), dev_allMaps)+1, // input maps: this is has_no_content_maps
      GET_POINTER(gpuMemoryManager, SINGLE_ARG(split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>*), dev_lists_to_replace), // output vecs
      NULL, // pass null to not store vector lengths
      GET_POINTER(gpuMemoryManager, vmesh::VelocityMesh*, dev_vmeshes), // rule_meshes
      GET_POINTER(gpuMemoryManager, SINGLE_ARG(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>*), dev_allMaps)+1, // rule_maps: this is has_no_content_maps
      GET_POINTER(gpuMemoryManager, split::SplitVector<vmesh::GlobalID>*, dev_lists_with_replace_new), // rule_vectors
      nCells,
      baseStream
      );
   CHK_ERR( gpuStreamSynchronize(baseStream) );
   extractKeysTimer.stop();

   // Call sub-function for actual block adjustment (including resizing vmeshes)
   // This same sub-function is also called from acceleration (TODO).
   // We don't need to give host or device arrays as parameters as they are universal.
   uint largestBlocksToChange = 0;
   uint largestBlocksBeforeOrAfter = 0;
   batch_adjust_blocks_caller(mpiGrid,
                              cellsToAdjust,
                              0, // no offset
                              largestBlocksToChange,
                              largestBlocksBeforeOrAfter,
                              popID);

   /* Batch tombstone cleaning
    * Extract all entries (GID,LID) which are overflown (see Hashinator for further details). At same time,
    * remove tombstones and overflown elements.
    *
    * By calling a few kernels which operate over all spatial cells at once instead of launching a few kernels per cell,
    * we reduce operational time by circa 10x.
    */
   phiprof::Timer tombstoneTimer {"GPU batch clean tombstones"};
   auto rule_overflown = []
      __device__(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID> *map,
                 Hashinator::hash_pair<vmesh::GlobalID, vmesh::LocalID>& kval) -> bool {
                            const vmesh::GlobalID emptybucket = map->get_emptybucket();
                            const vmesh::GlobalID tombstone   = map->get_tombstone();
                            if (kval.first == emptybucket) {
                               return false;
                            }
                            if (kval.first == tombstone) {
                               // Note: tombstones preceding overflown are deleted, so
                               // resetting overflown elements after this cannot rely on
                               // tombstones.
                               kval.first = emptybucket;
                               return false;
                            }
                            const size_t currentSizePower = map->getSizePower();
                            Hashinator::hash_pair<vmesh::GlobalID, vmesh::LocalID> *bck_ptr = map->expose_bucketdata<false>();
                            //const size_t hashIndex = Hashinator::HashFunction::_hash(kval.first, currentSizePower);
                            const size_t hashIndex = map->hash(kval.first);
                            const int bitMask = (1 << (currentSizePower)) - 1;
                            const bool isOverflown = (bck_ptr[hashIndex & bitMask].first != kval.first);
                            return isOverflown;
                         };
   clean_tombstones_launcher<decltype(rule_overflown)>(
      GET_POINTER(gpuMemoryManager, vmesh::VelocityMesh*, dev_vmeshes), // velocity meshes which include the hash maps to clean
      GET_POINTER(gpuMemoryManager, SINGLE_ARG(split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>*), dev_lists_with_replace_old), // use this for storing overflown elements
      GET_POINTER(gpuMemoryManager, vmesh::LocalID, dev_overflownElements), // return values: n_overflown_elements
      rule_overflown,
      nCells,
      baseStream
      );
   // Re-insert overflown elements back in vmeshes. First calculate
   // Launch parameters after using blocking memcpy to get overflow counts
   CHK_ERR( gpuMemcpyAsync(GET_POINTER(gpuMemoryManager, vmesh::LocalID, host_overflownElements), GET_POINTER(gpuMemoryManager, vmesh::LocalID, dev_overflownElements), nCells*sizeof(vmesh::LocalID), gpuMemcpyDeviceToHost, baseStream) );
   CHK_ERR( gpuStreamSynchronize(baseStream) );
   uint largestOverflow = 0;
#pragma omp parallel
   {
      uint thread_largestOverflow = 0;
#pragma omp for schedule(static)
      for (size_t i=0; i<nCells; ++i) {
         thread_largestOverflow = std::max(thread_largestOverflow, (GET_POINTER(gpuMemoryManager, vmesh::LocalID, host_overflownElements))[i]);
      }
#pragma omp critical
      {
         largestOverflow = std::max(thread_largestOverflow, largestOverflow);
      }
   } // end parallel region
   if (largestOverflow > 0) {
      dim3 grid_reinsert(largestOverflow,nCells,1);
      batch_insert_kernel<<<grid_reinsert, GPUTHREADS, 0, baseStream>>>(
         GET_POINTER(gpuMemoryManager, vmesh::VelocityMesh*, dev_vmeshes), // velocity meshes which include the hash maps to clean
         GET_POINTER(gpuMemoryManager, SINGLE_ARG(split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>*), dev_lists_with_replace_old) // use this for storing overflown elements
         );
      CHK_ERR( gpuPeekAtLastError() );
      CHK_ERR( gpuStreamSynchronize(baseStream) );
   }
   tombstoneTimer.stop();

#pragma omp parallel
   {
#pragma omp for schedule(dynamic)
      for (size_t i=0; i<nCells; ++i) {
         SpatialCell* SC = mpiGrid[cellsToAdjust[i]];
         if (SC->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
            SC->get_velocity_mesh(popID)->setNewCachedSize(0);
            SC->get_velocity_blocks(popID)->setNewCachedSize(0);
            continue;
         }
         // Perform hashmap cleanup here (instead of at acceleration mid-steps)
         phiprof::Timer cleanupTimer {cleanupId};
         //SC->get_velocity_mesh(popID)->gpu_cleanHashMap(gpu_getStream());
         //SC->dev_upload_population(popID);
         cleanupTimer.stop();

         phiprof::Timer postTimer {adjustPostId};
         #ifdef DEBUG_SPATIAL_CELL
         // Not re-doing old debug here, this should be enough
         SC->checkSizes(popID);
         #endif
         #ifdef DEBUG_VLASIATOR
         // This is a bit extreme
         SC->checkMesh(popID);
         #endif

         if (getObjectWrapper().particleSpecies[popID].sparse_conserve_mass) {
            // Block adjustment can only add empty blocks or delete existing blocks,
            // So post_adjust density must be equal to pre_adjust density minus mass loss.
            SC->density_post_adjust = SC->density_pre_adjust - (GET_POINTER(gpuMemoryManager, Real, host_massLoss))[i];
            if ( (SC->density_post_adjust > 0.0) && ((GET_POINTER(gpuMemoryManager, Real, host_massLoss))[i] != 0) ) {
               //SC->scale_population(SC->density_pre_adjust/SC->density_post_adjust, popID);
               // Now use the massloss buffer for the scaling value
               const Real mass_scaling = SC->density_pre_adjust/SC->density_post_adjust;
               (GET_POINTER(gpuMemoryManager, Real, host_massLoss))[i] = mass_scaling;
            } else {
               // Skip scaling this cell
               (GET_POINTER(gpuMemoryManager, Real, host_massLoss))[i] = 0;
            }
         } // end if conserve mass
         postTimer.stop();
      } // end cell loop
   } // end parallel region

   if ( (getObjectWrapper().particleSpecies[popID].sparse_conserve_mass)
        && (largestBlocksToChange > 0) ) {
      phiprof::Timer massConservationTimer {"GPU batch conserve mass"};
      CHK_ERR( gpuMemcpyAsync(GET_POINTER(gpuMemoryManager, Real, dev_massLoss), GET_POINTER(gpuMemoryManager, Real, host_massLoss), nCells*sizeof(Real), gpuMemcpyHostToDevice, baseStream) );
      // Launch parameters: Although post-adjustment, some VBCs can have more blocks than when entering
      // block adjustment, any new blocks will be empty and thus do not need to be scaled. Thus, we can use
      // The count which is the gathered max value over all cells of a counter which is either blocksBeforeAdjust
      // or BlocksAfterAdjust, whichever is smaller.

      // Third argument specifies the number of bytes in *shared memory* that is
      // dynamically allocated per block for this call in addition to the statically allocated memory.
      dim3 grid_mass_conservation(largestBlocksBeforeOrAfter,nCells,1);
      batch_population_scale_kernel<<<grid_mass_conservation, WID3, 0, baseStream>>> (
         GET_POINTER(gpuMemoryManager, vmesh::VelocityBlockContainer*, dev_VBCs),
         GET_POINTER(gpuMemoryManager, Real, dev_massLoss) // used now for scaling parameter
         );
      CHK_ERR( gpuPeekAtLastError() );
      CHK_ERR( gpuStreamSynchronize(baseStream) );
   }
}

void clear_maps_caller(const uint nCells,
                       const size_t largestSizePower,
                       gpuStream_t stream,
                       const size_t offset
   ) {
   const size_t largestMapSize = std::pow(2,largestSizePower);
   // fast ceil for positive ints
   //const size_t blocksNeeded = 1 + ((largestMapSize - 1) / Hashinator::defaults::MAX_BLOCKSIZE);
   size_t blocksNeeded = 1 + floor(sqrt(largestMapSize / Hashinator::defaults::MAX_BLOCKSIZE)-1);
   blocksNeeded = std::max((size_t)1, blocksNeeded);
   dim3 grid1(blocksNeeded,nCells,2);
   batch_reset_all_to_empty<<<grid1, Hashinator::defaults::MAX_BLOCKSIZE, 0, stream>>>(
      GET_POINTER(gpuMemoryManager, SINGLE_ARG(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>*), dev_allMaps)+2*offset
      );
   CHK_ERR( gpuPeekAtLastError() );
   CHK_ERR( gpuStreamSynchronize(stream) );
}


void batch_adjust_blocks_caller(
   dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   const vector<CellID>& cellsToAdjust,
   const uint cellOffset,
   uint &out_largestBlocksToChange,
   uint &out_largestBlocksBeforeOrAfter,
   const uint popID
   ) {

   const uint nCells = cellsToAdjust.size();
   if (nCells == 0) {
      return;
   }
   const gpuStream_t baseStream = gpu_getStream();

   // Resizes are faster this way with larger grid and single thread per block.
   phiprof::Timer deviceResizeTimer {"GPU resize mesh on-device"};
   batch_resize_vbc_kernel_pre<<<nCells, 1, 0, baseStream>>> (
      GET_POINTER(gpuMemoryManager, vmesh::VelocityMesh*, dev_vmeshes)+cellOffset,
      GET_POINTER(gpuMemoryManager, vmesh::VelocityBlockContainer*, dev_VBCs)+cellOffset,
      GET_POINTER(gpuMemoryManager, split::SplitVector<vmesh::GlobalID>*, dev_lists_with_replace_new)+cellOffset,
      GET_POINTER(gpuMemoryManager, SINGLE_ARG(split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>*), dev_lists_delete)+cellOffset,
      GET_POINTER(gpuMemoryManager, SINGLE_ARG(split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>*), dev_lists_to_replace)+cellOffset,
      GET_POINTER(gpuMemoryManager, SINGLE_ARG(split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>*), dev_lists_with_replace_old)+cellOffset,
      GET_POINTER(gpuMemoryManager, vmesh::LocalID, dev_nBefore)+cellOffset,
      GET_POINTER(gpuMemoryManager, vmesh::LocalID, dev_nAfter)+cellOffset,
      GET_POINTER(gpuMemoryManager, vmesh::LocalID, dev_nBlocksToChange)+cellOffset,
      GET_POINTER(gpuMemoryManager, vmesh::LocalID, dev_resizeSuccess)+cellOffset,
      GET_POINTER(gpuMemoryManager, Real, dev_massLoss)+cellOffset // mass loss, set to zero
      );
   CHK_ERR( gpuPeekAtLastError() );
   CHK_ERR( gpuMemcpyAsync(GET_POINTER(gpuMemoryManager, vmesh::LocalID, host_nBefore)+cellOffset, GET_POINTER(gpuMemoryManager, vmesh::LocalID, dev_nBefore)+cellOffset, nCells*sizeof(vmesh::LocalID), gpuMemcpyDeviceToHost, baseStream) );
   CHK_ERR( gpuMemcpyAsync(GET_POINTER(gpuMemoryManager, vmesh::LocalID, host_nAfter)+cellOffset, GET_POINTER(gpuMemoryManager, vmesh::LocalID, dev_nAfter)+cellOffset, nCells*sizeof(vmesh::LocalID), gpuMemcpyDeviceToHost, baseStream) );
   CHK_ERR( gpuMemcpyAsync(GET_POINTER(gpuMemoryManager, vmesh::LocalID, host_nBlocksToChange)+cellOffset, GET_POINTER(gpuMemoryManager, vmesh::LocalID, dev_nBlocksToChange)+cellOffset, nCells*sizeof(vmesh::LocalID), gpuMemcpyDeviceToHost, baseStream) );
   CHK_ERR( gpuMemcpyAsync(GET_POINTER(gpuMemoryManager, vmesh::LocalID, host_resizeSuccess)+cellOffset, GET_POINTER(gpuMemoryManager, vmesh::LocalID, dev_resizeSuccess)+cellOffset, nCells*sizeof(vmesh::LocalID), gpuMemcpyDeviceToHost, baseStream) );
   CHK_ERR( gpuStreamSynchronize(baseStream) );
   deviceResizeTimer.stop();

   phiprof::Timer hostResizeTimer {"GPU resize mesh from host "};
   uint largestBlocksToChange = 0;
   uint largestBlocksBeforeOrAfter = 0;
   // This loop appears to be faster non-threaded!
   for (size_t i=0; i<nCells; ++i) {
      SpatialCell* SC = mpiGrid[cellsToAdjust[i]];
      if (SC->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
         continue;
      }
      // Grow mesh if necessary and on-device resize did not work??
      const vmesh::LocalID nBlocksBeforeAdjust = (GET_POINTER(gpuMemoryManager, vmesh::LocalID, host_nBefore))[i+cellOffset];
      const vmesh::LocalID nBlocksAfterAdjust  = (GET_POINTER(gpuMemoryManager, vmesh::LocalID, host_nAfter))[i+cellOffset];
      const vmesh::LocalID nBlocksToChange     = (GET_POINTER(gpuMemoryManager, vmesh::LocalID, host_nBlocksToChange))[i+cellOffset];
      const vmesh::LocalID resizeDevSuccess    = (GET_POINTER(gpuMemoryManager, vmesh::LocalID, host_resizeSuccess))[i+cellOffset];
      largestBlocksToChange = std::max(largestBlocksToChange, nBlocksToChange);
      // This is gathered for mass loss correction: for each cell, we want the smaller of either blocks before or after. Then,
      // we want to gather the largest of those values.
      const vmesh::LocalID lowBlocks = std::min(nBlocksBeforeAdjust, nBlocksAfterAdjust);
      largestBlocksBeforeOrAfter = std::max(largestBlocksBeforeOrAfter, lowBlocks);
      if ( (nBlocksAfterAdjust > nBlocksBeforeAdjust) && (resizeDevSuccess == 0)) {
         //GPUTODO is _FACTOR enough instead of _PADDING?
         SC->get_velocity_mesh(popID)->setNewCapacity(nBlocksAfterAdjust*BLOCK_ALLOCATION_PADDING);
         SC->get_velocity_mesh(popID)->setNewSize(nBlocksAfterAdjust);
         SC->get_velocity_blocks(popID)->setNewCapacity(nBlocksAfterAdjust*BLOCK_ALLOCATION_PADDING);
         SC->get_velocity_blocks(popID)->setNewSize(nBlocksAfterAdjust);
         SC->dev_upload_population(popID);
      }
      // Update cached sizes
      SC->get_velocity_mesh(popID)->setNewCachedSize(nBlocksAfterAdjust);
      SC->get_velocity_blocks(popID)->setNewCachedSize(nBlocksAfterAdjust);
   } // end cell loop
   CHK_ERR( gpuDeviceSynchronize() );
   hostResizeTimer.stop();
   // Writing directly into pass-by-reference variables from within OMP parallel region caused issues
   out_largestBlocksToChange = largestBlocksToChange;
   out_largestBlocksBeforeOrAfter = largestBlocksBeforeOrAfter;

   // Do we actually have any changes to perform?
   if (largestBlocksToChange > 0) {
      phiprof::Timer addRemoveKernelTimer {"GPU batch add and remove blocks kernel"};
      // Third argument specifies the number of bytes in *shared memory* that is
      // dynamically allocated per block for this call in addition to the statically allocated memory.
      dim3 grid_addremove(largestBlocksToChange,nCells,1);
      // Launch grid is sized so that for all spatial cells, we launch up to the maximum number of required
      // operations (add a block, delete a block, replace a block with a new one, replace a block with an existing one)
      batch_update_velocity_blocks_kernel<<<grid_addremove, WID3, 0, baseStream>>> (
         GET_POINTER(gpuMemoryManager, vmesh::VelocityMesh*, dev_vmeshes)+cellOffset,
         GET_POINTER(gpuMemoryManager, vmesh::VelocityBlockContainer*, dev_VBCs)+cellOffset,
         GET_POINTER(gpuMemoryManager, split::SplitVector<vmesh::GlobalID>*, dev_lists_with_replace_new)+cellOffset,
         GET_POINTER(gpuMemoryManager, SINGLE_ARG(split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>*), dev_lists_delete)+cellOffset,
         GET_POINTER(gpuMemoryManager, SINGLE_ARG(split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>*), dev_lists_to_replace)+cellOffset,
         GET_POINTER(gpuMemoryManager, SINGLE_ARG(split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>*), dev_lists_with_replace_old)+cellOffset,
         GET_POINTER(gpuMemoryManager, vmesh::LocalID, dev_nBefore)+cellOffset,
         GET_POINTER(gpuMemoryManager, vmesh::LocalID, dev_nAfter)+cellOffset,
         GET_POINTER(gpuMemoryManager, vmesh::LocalID, dev_nBlocksToChange)+cellOffset,
         GET_POINTER(gpuMemoryManager, Real, dev_massLoss)+cellOffset
         );
      CHK_ERR( gpuPeekAtLastError() );
      // Pull mass loss values to host
      CHK_ERR( gpuMemcpyAsync(GET_POINTER(gpuMemoryManager, Real, host_massLoss)+cellOffset, GET_POINTER(gpuMemoryManager, Real, dev_massLoss)+cellOffset, nCells*sizeof(Real), gpuMemcpyDeviceToHost, baseStream) );
      CHK_ERR( gpuStreamSynchronize(baseStream) );
      // Update mass Loss (not worth threading)
      for (size_t i=0; i<nCells; ++i) {
         SpatialCell* SC = mpiGrid[cellsToAdjust[i]];
         if (SC->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
            continue;
         }
         SC->increment_mass_loss(popID, (GET_POINTER(gpuMemoryManager, Real, host_massLoss))[i]);
      }
      addRemoveKernelTimer.stop();

      // Should not re-allocate on shrinking, so do on-device
      phiprof::Timer deviceResizePostTimer {"GPU resize mesh on-device post"};
      // Resizes are faster this way with larger grid and single thread
      batch_resize_vbc_kernel_post<<<nCells, 1, 0, baseStream>>> (
         GET_POINTER(gpuMemoryManager, vmesh::VelocityMesh*, dev_vmeshes)+cellOffset,
         GET_POINTER(gpuMemoryManager, vmesh::VelocityBlockContainer*, dev_VBCs)+cellOffset,
         GET_POINTER(gpuMemoryManager, vmesh::LocalID, dev_nAfter)+cellOffset
         );
      CHK_ERR( gpuPeekAtLastError() );
      CHK_ERR( gpuStreamSynchronize(baseStream) );
      deviceResizePostTimer.stop();
   }
}

void extract_to_replace_caller(
   Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>** input_maps,
   split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>> **output_vecs,
   vmesh::LocalID* output_sizes,
   vmesh::VelocityMesh** rule_meshes,
   Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>** rule_maps,
   split::SplitVector<vmesh::GlobalID>** rule_vectors,
   const uint nCells,
   gpuStream_t stream
   ) {
   auto rule_to_replace = [] __device__(const Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID> *map,
                                        const Hashinator::hash_pair<vmesh::GlobalID, vmesh::LocalID>& kval,
                                        const vmesh::LocalID threshold,
                                        const vmesh::LocalID invalidLID,
                                        const vmesh::GlobalID invalidGID) -> bool {
                             const vmesh::GlobalID emptybucket = map->get_emptybucket();
                             const vmesh::GlobalID tombstone   = map->get_tombstone();
                             return kval.first  != emptybucket &&
                                    kval.first  != tombstone   &&
                                    kval.first  != invalidGID  &&
                                    kval.second <  threshold   &&
                                    kval.second != invalidLID;
                          };

   // Find Blocks (GID,LID) to be replaced with new ones
   extract_GIDs_kernel_launcher<decltype(rule_to_replace),Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>,false>(
      input_maps,
      output_vecs,
      output_sizes,
      rule_to_replace,
      rule_meshes,
      rule_maps,
      rule_vectors,
      nCells,
      stream
      );
}

void extract_to_delete_or_move_caller(
   Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>** input_maps,
   split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>> **output_vecs,
   vmesh::LocalID* output_sizes,
   vmesh::VelocityMesh** rule_meshes,
   Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>** rule_maps,
   split::SplitVector<vmesh::GlobalID>** rule_vectors,
   const uint nCells,
   gpuStream_t stream
   ) {
   auto rule_delete_move = [] __device__(const Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID> *map,
                                         const Hashinator::hash_pair<vmesh::GlobalID, vmesh::LocalID>& kval,
                                         const vmesh::LocalID threshold,
                                         const vmesh::LocalID invalidLID,
                                         const vmesh::GlobalID invalidGID) -> bool {
                              const vmesh::GlobalID emptybucket = map->get_emptybucket();
                              const vmesh::GlobalID tombstone   = map->get_tombstone();
                              return kval.first  != emptybucket &&
                                     kval.first  != tombstone   &&
                                     kval.first  != invalidGID  &&
                                     kval.second >= threshold   &&
                                     kval.second != invalidLID;
                           };
   extract_GIDs_kernel_launcher<decltype(rule_delete_move),Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>,false>(
      input_maps,
      output_vecs,
      output_sizes,
      rule_delete_move,
      rule_meshes,
      rule_maps,
      rule_vectors,
      nCells,
      stream
      );
}

void extract_to_add_caller(
   Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>** input_maps,
   split::SplitVector<vmesh::GlobalID> **output_vecs,
   vmesh::LocalID* output_sizes,
   vmesh::VelocityMesh** rule_meshes,
   Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>** rule_maps,
   split::SplitVector<vmesh::GlobalID>** rule_vectors,
   const uint nCells,
   gpuStream_t stream
   ) {
   auto rule_add = [] __device__(const Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID> *map,
                                 const Hashinator::hash_pair<vmesh::GlobalID, vmesh::LocalID>& kval,
                                 const vmesh::LocalID threshold,
                                 const vmesh::LocalID invalidLID,
                                 const vmesh::GlobalID invalidGID) -> bool {
                      // This rule does not use the threshold value
                      const vmesh::GlobalID emptybucket = map->get_emptybucket();
                      const vmesh::GlobalID tombstone   = map->get_tombstone();
                      return kval.first != emptybucket &&
                             kval.first != tombstone   &&
                             kval.first != invalidGID  &&
                             // Required GIDs which do not yet exist in vmesh were stored in
                             // velocity_block_with_content_map with kval.second==invalidLID
                             kval.second == invalidLID;
                   };
   extract_GIDs_kernel_launcher<decltype(rule_add),vmesh::GlobalID,true>(
      input_maps,
      output_vecs,
      output_sizes,
      rule_add,
      rule_meshes,
      rule_maps,
      rule_vectors,
      nCells,
      stream
      );
}

} // namespace
