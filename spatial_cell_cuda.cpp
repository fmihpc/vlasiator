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
 */

#include <unordered_set>

#include "spatial_cell_cuda.hpp"
#include "cuda_context.cuh"

//#include "velocity_blocks.h" // not needed in semilag
#include "object_wrapper.h"
#include "velocity_mesh_parameters.h"

#ifndef NDEBUG
   #define DEBUG_SPATIAL_CELL
#endif

using namespace std;

/** CUDA kernel for identifying which blocks have relevant content */
__global__ void update_velocity_block_content_lists_kernel (
   vmesh::VelocityMesh *vmesh,
   vmesh::VelocityBlockContainer *blockContainer,
   split::SplitVector<vmesh::GlobalID>* velocity_block_with_content_list,
   split::SplitVector<vmesh::GlobalID>* velocity_block_with_no_content_list,
   Realf velocity_block_min_value
   ) {

   const int cudaBlocks = gridDim.x;
   const int blocki = blockIdx.x;
   const int i = threadIdx.x;
   const int j = threadIdx.y;
   const int k = threadIdx.z;
   const uint ti = k*WID2 + j*WID + i;
   __shared__ bool has_content[WID3];
   const uint nBlocks = vmesh->size();
   for (uint blockLID=blocki; blockLID<nBlocks; blockLID += cudaBlocks) {
      vmesh::GlobalID blockGID = vmesh->getGlobalID(blockLID);
      #ifdef DEBUG_SPATIAL_CELL
      if (blockGID == vmesh->invalidGlobalID()) {
         continue;
      }
      if (blockLID == vmesh->invalidLocalID()) {
         continue;
      }
      #endif
      Realf* avgs = blockContainer->getData(blockLID);
      has_content[ti] = avgs[ti] >= velocity_block_min_value ? true : false;
      __syncthreads();
      // Implemented just a simple non-optimized thread OR
      for (unsigned int s=WID3/2; s>0; s>>=1) {
         if (ti < s) {
            has_content[ti] = has_content[ti] || has_content[ti + s];
         }
         __syncthreads();
      }
      if (ti==0) {
         if (has_content[0]) {
            velocity_block_with_content_list->device_push_back(blockGID);
         } else {
            velocity_block_with_no_content_list->device_push_back(blockGID);
         }
      }
      __syncthreads();
   }
}

/** CUDA kernel for updating blocks based on generated lists */
__global__ void update_velocity_blocks_kernel(
   vmesh::VelocityMesh *vmesh,
   vmesh::VelocityBlockContainer *blockContainer,
   split::SplitVector<vmesh::GlobalID>* BlocksToAdd,
   split::SplitVector<vmesh::GlobalID>* BlocksToRemove,
   split::SplitVector<vmesh::GlobalID>* BlocksToMove,
   vmesh::LocalID nBlocksBeforeAdjust,
   vmesh::LocalID nBlocksAfterAdjust,
   vmesh::LocalID *AddVectorIndex,
   vmesh::LocalID *MoveVectorIndex,
   Realf* dev_rhoLossAdjust
   ) {
   const int cudaBlocks = gridDim.x;
   const int blocki = blockIdx.x;
   const int i = threadIdx.x;
   const int j = threadIdx.y;
   const int k = threadIdx.z;
   const uint ti = k*WID2 + j*WID + i;

   const vmesh::LocalID nToAdd = BlocksToAdd->size();
   const vmesh::LocalID nToRemove = BlocksToRemove->size();
   const vmesh::LocalID nToMove = BlocksToMove->size();
   const vmesh::LocalID nToCreate = nToAdd > nToRemove ? (nToAdd-nToRemove) : 0;

   // For tracking mass-loss
   __shared__ Realf massloss[WID3];
   __shared__ vmesh::LocalID moveIndex;
   __shared__ vmesh::LocalID addIndex;

   for (vmesh::LocalID m=blocki; m<nToRemove; m += cudaBlocks) {
      // Go through all blocks which are to be removed.
      // If there is a corresponding block to be added, place that in its stead.
      // Otherwise, take the corresponding block from the moved list instead.
      const vmesh::GlobalID rmGID = BlocksToRemove->at(m);
      const vmesh::LocalID rmLID = vmesh->getLocalID(rmGID);

      #ifdef DEBUG_SPATIAL_CELL
      if (rmGID == vmesh->invalidGlobalID()) {
         continue;
      }
      if (rmLID == vmesh->invalidLocalID()) {
         continue;
      }
      #endif

      // Track mass loss:
      Realf* rm_avgs = blockContainer->getData(rmLID);
      Real* rm_block_parameters = blockContainer->getParameters(rmLID);
      const Real rm_DV3 = rm_block_parameters[BlockParams::DVX]
         * rm_block_parameters[BlockParams::DVY]
         * rm_block_parameters[BlockParams::DVZ];

      // thread-sum for rho
      massloss[ti] = rm_avgs[ti]*rm_DV3;
      __syncthreads();
      // Implemented just a simple non-optimized thread sum
      for (unsigned int s=WID3/2; s>0; s>>=1) {
         if (ti < s) {
            massloss[ti] += massloss[ti + s];
         }
         __syncthreads();
      }

      if (ti==0) {
         // Bookkeeping only by one thread
         size_t old= atomicAdd(dev_rhoLossAdjust, massloss[0]);
      }
      __syncthreads();

      // Figure out which GID to put here instead
      if (rmLID >= nBlocksAfterAdjust) {
         // Delete without replacing
         if (ti==0) {
            vmesh->deleteBlock(rmGID,rmLID);
         }
         __syncthreads();
         continue;
      }
      if (ti==0) moveIndex = atomicAdd(MoveVectorIndex,1);
      __syncthreads();
      if (moveIndex<nToMove) {
         // Move from latter part of vmesh
         const vmesh::GlobalID replaceGID = BlocksToMove->at(moveIndex);
         const vmesh::LocalID replaceLID = vmesh->getLocalID(replaceGID);
         Realf* repl_avgs = blockContainer->getData(replaceLID);
         Real*  repl_block_parameters = blockContainer->getParameters(replaceLID);
         rm_avgs[ti] = repl_avgs[ti];
         if (ti < BlockParams::N_VELOCITY_BLOCK_PARAMS) {
            rm_block_parameters[ti] = repl_block_parameters[ti];
         }
         __syncthreads();
         if (ti==0) {
            // Remove hashmap entry for removed block, add instead created block
            vmesh->replaceBlock(rmGID,rmLID,replaceGID);
         }
         __syncthreads();
         #ifdef DEBUG_SPATIAL_CELL
         if (vmesh->getGlobalID(rmLID) == vmesh->invalidGlobalID()) {
            continue;
         }
         if (vmesh->getLocalID(replaceGID) == vmesh->invalidLocalID()) {
            continue;
         }
         #endif
         continue;
      }
      if (ti==0) addIndex = atomicAdd(AddVectorIndex,1);
      #ifdef DEBUG_SPATIAL_CELL
      if (ti==0) printf("addIndex %d m %d\n",addIndex,m);
      #endif
      __syncthreads();
      if (addIndex<nToAdd) {
         // New GID
         const vmesh::GlobalID addGID = BlocksToAdd->at(addIndex);
         #ifdef DEBUG_SPATIAL_CELL
         if (addGID == vmesh->invalidGlobalID()) {
            __syncthreads();
            continue;
         }
         #endif
         rm_avgs[ti] = 0;
         if (ti==0) {
            // Write in block parameters
            vmesh->getCellSize(addGID,&(rm_block_parameters[BlockParams::DVX]));
            vmesh->getBlockCoordinates(addGID,&(rm_block_parameters[BlockParams::VXCRD]));
            vmesh->replaceBlock(rmGID,rmLID,addGID);
         }
         __syncthreads();
         #ifdef DEBUG_SPATIAL_CELL
         if (vmesh->getGlobalID(rmLID) == vmesh->invalidGlobalID()) {
            continue;
         }
         if (vmesh->getLocalID(addGID) == vmesh->invalidLocalID()) {
            continue;
         }
         #endif
         continue;
      }
      #ifdef DEBUG_SPATIAL_CELL
      if (ti==0) {
         printf("Error! Fall through in update_velocity_blocks_kernel! \n");
      }
      #endif
      __syncthreads();
   }
   // Now, if we need to expand the size of the vmesh, let's add blocks.
   // For thread-safety,this assumes that the localToGlobalMap is already of sufficient size, as should be
   // the block_data and block_parameters vectors.
   for (vmesh::LocalID m=blocki; m<nToCreate; m += cudaBlocks) {
      // Debug check: if we are adding elements, then nToMove should be zero
      // We have already used up nToRemove entries from the addition vector.
      const vmesh::GlobalID addGID = BlocksToAdd->at(nToRemove+m);
      // We need to add the data of addGID to a new LID:
      const vmesh::LocalID addLID = nBlocksBeforeAdjust + m;
      Realf* add_avgs = blockContainer->getData(addLID);
      Real* add_block_parameters = blockContainer->getParameters(addLID);
      // Zero out blockdata
      add_avgs[ti] = 0;
      if (ti==0) {
         // Write in block parameters
         vmesh->getCellSize(addGID,&(add_block_parameters[BlockParams::DVX]));
         vmesh->getBlockCoordinates(addGID,&(add_block_parameters[BlockParams::VXCRD]));
         vmesh->placeBlock(addGID,addLID);
      }
      __syncthreads();
      #ifdef DEBUG_SPATIAL_CELL
      if (vmesh->getGlobalID(addLID) == vmesh->invalidGlobalID()) {
         continue;
      }
      if (vmesh->getLocalID(addGID) == vmesh->invalidLocalID()) {
         continue;
      }
      #endif
   }

}


namespace spatial_cell {
   int SpatialCell::activePopID = 0;
   uint64_t SpatialCell::mpi_transfer_type = 0;
   bool SpatialCell::mpiTransferAtSysBoundaries = false;
   bool SpatialCell::mpiTransferInAMRTranslation = false;
   int SpatialCell::mpiTransferXYZTranslation = 0;

   SpatialCell::SpatialCell() {
      // Block list and cache always have room for all blocks
      this->sysBoundaryLayer=0; // Default value, layer not yet initialized
      for (unsigned int i=0; i<WID3; ++i) {
         null_block_data[i] = 0.0;
      }

      // reset spatial cell parameters
      for (unsigned int i = 0; i < CellParams::N_SPATIAL_CELL_PARAMS; i++) {
         this->parameters[i]=0.0;
      }

      // reset BVOL derivatives
      for (unsigned int i = 0; i < bvolderivatives::N_BVOL_DERIVATIVES; i++) {
         this->derivativesBVOL[i]=0;
      }

      for (unsigned int i = 0; i < MAX_NEIGHBORS_PER_DIM; ++i) {
         this->neighbor_number_of_blocks[i] = 0;
         this->neighbor_block_data[i] = NULL;
      }

      //is transferred by default
      this->mpiTransferEnabled=true;

      // Set correct number of populations
      populations.resize(getObjectWrapper().particleSpecies.size());

      // Set velocity meshes
      for (uint popID=0; popID<populations.size(); ++popID) {
         const species::Species& spec = getObjectWrapper().particleSpecies[popID];
         populations[popID].vmesh->initialize(spec.velocityMesh);
         populations[popID].velocityBlockMinValue = spec.sparseMinValue;
         populations[popID].N_blocks = 0;
      }

      // SplitVectors via pointers for unified memory
      velocity_block_with_content_list = new split::SplitVector<vmesh::GlobalID>(1);
      velocity_block_with_no_content_list = new split::SplitVector<vmesh::GlobalID>(1);
      BlocksToAdd = new split::SplitVector<vmesh::GlobalID>(1);
      BlocksToRemove = new split::SplitVector<vmesh::GlobalID>(1);
      BlocksToMove = new split::SplitVector<vmesh::GlobalID>(1);
      velocity_block_with_content_list->clear();
      velocity_block_with_no_content_list->clear();
      BlocksToAdd->clear();
      BlocksToRemove->clear();
      BlocksToMove->clear();
   }

   SpatialCell::~SpatialCell() {
      delete velocity_block_with_content_list;
      delete velocity_block_with_no_content_list;
      delete BlocksToAdd;
      delete BlocksToRemove;
      delete BlocksToMove;
   }

   SpatialCell::SpatialCell(const SpatialCell& other) {
      // These should be empty when created, but let's play safe.
      velocity_block_with_content_list = new split::SplitVector<vmesh::GlobalID>(*(other.velocity_block_with_content_list));
      velocity_block_with_no_content_list = new split::SplitVector<vmesh::GlobalID>(*(other.velocity_block_with_no_content_list));
      BlocksToAdd = new split::SplitVector<vmesh::GlobalID>(1);
      BlocksToRemove = new split::SplitVector<vmesh::GlobalID>(1);
      BlocksToMove = new split::SplitVector<vmesh::GlobalID>(1);
      BlocksToAdd->clear();
      BlocksToRemove->clear();
      BlocksToMove->clear();

      // Member variables
      ioLocalCellId = other.ioLocalCellId;
      sysBoundaryFlag = other.sysBoundaryFlag;
      sysBoundaryLayer = other.sysBoundaryLayer;
      sysBoundaryLayerNew = other.sysBoundaryLayerNew;
      velocity_block_with_content_list_size = other.velocity_block_with_content_list_size;
      initialized = other.initialized;
      mpiTransferEnabled = other.mpiTransferEnabled;
      for (unsigned int i=0; i<bvolderivatives::N_BVOL_DERIVATIVES; ++i) {
         derivativesBVOL[i] = other.derivativesBVOL[i];
      }
      for (unsigned int i=0; i<CellParams::N_SPATIAL_CELL_PARAMS; ++i) {
         parameters[i] = other.parameters[i];
      }
      for (unsigned int i=0; i<WID3; ++i) {
         null_block_data[i] = 0.0;
      }
      for (unsigned int i=0; i<MAX_NEIGHBORS_PER_DIM; ++i) {
         neighbor_block_data[i] = other.neighbor_block_data[i];
         neighbor_number_of_blocks[i] = other.neighbor_number_of_blocks[i];
      }

      if (other.face_neighbor_ranks.size()>0) {
         face_neighbor_ranks = std::map<int,std::set<int>>(other.face_neighbor_ranks);
      }
      if (other.populations.size()>0) {
         populations = std::vector<spatial_cell::Population>(other.populations);
      }
   }
   const SpatialCell& SpatialCell::operator=(const SpatialCell& other) {
      // Delete old vectors
      delete velocity_block_with_content_list;
      delete velocity_block_with_no_content_list;
      delete BlocksToAdd;
      delete BlocksToRemove;
      delete BlocksToMove;

      // These should be empty when created, but let's play safe.
      velocity_block_with_content_list = new split::SplitVector<vmesh::GlobalID>(*(other.velocity_block_with_content_list));
      velocity_block_with_no_content_list = new split::SplitVector<vmesh::GlobalID>(*(other.velocity_block_with_no_content_list));
      // These start empty
      BlocksToAdd = new split::SplitVector<vmesh::GlobalID>(1);
      BlocksToRemove = new split::SplitVector<vmesh::GlobalID>(1);
      BlocksToMove = new split::SplitVector<vmesh::GlobalID>(1);
      BlocksToAdd->clear();
      BlocksToRemove->clear();
      BlocksToMove->clear();
      // Member variables
      ioLocalCellId = other.ioLocalCellId;
      sysBoundaryFlag = other.sysBoundaryFlag;
      sysBoundaryLayer = other.sysBoundaryLayer;
      sysBoundaryLayerNew = other.sysBoundaryLayerNew;
      velocity_block_with_content_list_size = other.velocity_block_with_content_list_size;
      initialized = other.initialized;
      mpiTransferEnabled = other.mpiTransferEnabled;
      for (unsigned int i=0; i<bvolderivatives::N_BVOL_DERIVATIVES; ++i) {
         derivativesBVOL[i] = other.derivativesBVOL[i];
      }
      for (unsigned int i=0; i<CellParams::N_SPATIAL_CELL_PARAMS; ++i) {
         parameters[i] = other.parameters[i];
      }
      for (unsigned int i=0; i<WID3; ++i) {
         null_block_data[i] = 0.0;
      }
      for (unsigned int i=0; i<MAX_NEIGHBORS_PER_DIM; ++i) {
         neighbor_block_data[i] = other.neighbor_block_data[i];
         neighbor_number_of_blocks[i] = other.neighbor_number_of_blocks[i];
      }

      face_neighbor_ranks = std::map<int,std::set<int>>(other.face_neighbor_ranks);
      populations = std::vector<spatial_cell::Population>(other.populations);

      return *this;
   }

   /** Adds "important" and removes "unimportant" velocity blocks
    * to/from this cell.
    *
    * velocity_block_with_content_list needs to be up to date in local and remote cells.
    * velocity_block_with_no_content_list needs to be up to date in local cells.
    *
    * update_velocity_block_with_content_lists() should have
    * been called with the current distribution function values, and then the contetn list transferred.
    *
    * Removes all velocity blocks from this spatial cell which don't
    * have content and don't have spatial or velocity neighbors with
    * content.  Adds neighbors for all velocity blocks which do have
    * content (including spatial neighbors).  All cells in
    * spatial_neighbors are assumed to be neighbors of this cell.
    *
    * This function is thread-safe when called for different cells
    * per thread. We need the block_has_content vector from
    * neighbouring cells, but these are not written to here. We only
    * modify local cell.*/

   void SpatialCell::adjust_velocity_blocks(const std::vector<SpatialCell*>& spatial_neighbors,
                                            const uint popID,bool doDeleteEmptyBlocks) {
      #ifdef DEBUG_SPATIAL_CELL
      if (popID >= populations.size()) {
         std::cerr << "ERROR, popID " << popID << " exceeds populations.size() " << populations.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
         exit(1);
      }
      #endif

      // Need the vmesh on CPU for checking against existing blocks
      populations[popID].vmesh->dev_prefetchHost();

      //  This set contains all those cellids which have neighbors in any
      //  of the 6-dimensions Actually, we would only need to add
      //  local blocks with no content here, as blocks with content
      //  do not need to be created and also will not be removed as
      //  we only check for removal for blocks with no content
      std::unordered_set<vmesh::GlobalID> neighbors_have_content;

      //add neighbor content info for velocity space neighbors to map. We loop over blocks
      //with content and raise the neighbors_have_content for
      //itself, and for all its neighbors
      for (vmesh::LocalID block_index=0; block_index<velocity_block_with_content_list->size(); ++block_index) {
         vmesh::GlobalID block = velocity_block_with_content_list->at(block_index);

         const velocity_block_indices_t indices = SpatialCell::get_velocity_block_indices(popID,block);
         neighbors_have_content.insert(block); //also add the cell itself

         int addWidthV = getObjectWrapper().particleSpecies[popID].sparseBlockAddWidthV;
         for (int offset_vx=-addWidthV;offset_vx<=addWidthV;offset_vx++) {
            for (int offset_vy=-addWidthV;offset_vy<=addWidthV;offset_vy++) {
               for (int offset_vz=-addWidthV;offset_vz<=addWidthV;offset_vz++) {
                  const vmesh::GlobalID neighbor_block
                     = get_velocity_block(popID,{{indices[0]+offset_vx,indices[1]+offset_vy,indices[2]+offset_vz}});
                  neighbors_have_content.insert(neighbor_block); //add all potential ngbrs of this block with content
               }
            }
         }
      }

      //add neighbor content info for spatial space neighbors to map. We loop over
      //neighbor cell lists with existing blocks, and raise the
      //flag for the local block with same block id
      for (std::vector<SpatialCell*>::const_iterator neighbor=spatial_neighbors.begin();
           neighbor != spatial_neighbors.end(); ++neighbor) {
         for (vmesh::LocalID block_index=0; block_index<(*neighbor)->velocity_block_with_content_list->size(); ++block_index) {
            vmesh::GlobalID block = (*neighbor)->velocity_block_with_content_list->at(block_index);
            neighbors_have_content.insert(block);
         }
      }

      BlocksToAdd->optimizeCPU();
      BlocksToRemove->optimizeCPU();
      BlocksToMove->optimizeCPU();
      BlocksToAdd->clear();
      BlocksToRemove->clear();
      BlocksToMove->clear();
      BlocksToAdd->reserve(populations[popID].vmesh->size()*BLOCK_ALLOCATION_PADDING);
      BlocksToRemove->reserve(populations[popID].vmesh->size());
      BlocksToMove->reserve(populations[popID].vmesh->size());

      // REMOVE all blocks in this cell without content + without neighbors with content
      // better to do it in the reverse order, as then blocks at the
      // end are removed first, and we may avoid copying extra data.
      if (doDeleteEmptyBlocks) {
         for (int block_index= this->velocity_block_with_no_content_list->size()-1; block_index>=0; --block_index) {
            const vmesh::GlobalID blockGID = velocity_block_with_no_content_list->at(block_index);
            #ifdef DEBUG_SPATIAL_CELL
            if (blockGID == invalidGlobalID()) {
               cerr << "Got invalid block at " << __FILE__ << ' ' << __LINE__ << endl;
               exit(1);
            }
            const vmesh::LocalID blockLID = get_velocity_block_local_id(blockGID,popID);
            if (blockLID == invalidLocalID()) {
               cerr << "Could not find block in " << __FILE__ << ' ' << __LINE__ << endl;
               exit(1);
            }
            #endif

            std::unordered_set<vmesh::GlobalID>::iterator it = neighbors_have_content.find(blockGID);
            if (it == neighbors_have_content.end()) {
               BlocksToRemove->push_back(blockGID);
            }
         }
      }

      // ADD all blocks with neighbors in spatial or velocity space (if it exists then the block is unchanged)
      for (std::unordered_set<vmesh::GlobalID>::iterator it=neighbors_have_content.begin(); it != neighbors_have_content.end(); ++it) {

         // Only add blocks which don't yet exist to optimize cuda parallel memory management
         // Should this check against a separate list? (now we prefetch the vmesh to CPU)
         // if ( populations[popID].vmesh->count(*it) == 0 ) {
         //    BlocksToAdd->push_back(*it);
         // }
         if ( populations[popID].vmesh->getLocalID(*it) == populations[popID].vmesh->invalidLocalID()  ) {
            BlocksToAdd->push_back(*it);
         }
      }

      // On-device adjustment calling happens in separate function as it is also called from within acceleration
      adjust_velocity_blocks_caller(popID);

      // Perform hashmap cleanup here (instead of at acceleration mid-steps)
      phiprof::start("Hashinator cleanup");
      populations[popID].vmesh->dev_prefetchDevice();
      populations[popID].vmesh->dev_cleanHashMap();
      phiprof::stop("Hashinator cleanup");

   }

   void SpatialCell::adjust_velocity_blocks_caller(const uint popID) {
      /**
          Call CUDA kernel with all necessary information for creation and deletion of blocks
      **/

      phiprof::start("CUDA add and remove blocks");
      const vmesh::LocalID nBlocksBeforeAdjust = populations[popID].vmesh->size();
      const vmesh::LocalID nBlocksAfterAdjust = nBlocksBeforeAdjust + BlocksToAdd->size() - BlocksToRemove->size();

      // Figure out which blocks must be kept, but need to be moved around inside blockContainer
      phiprof::start("Vmesh prefetch");
      populations[popID].vmesh->dev_prefetchHost();
      BlocksToMove->optimizeCPU();
      BlocksToMove->clear();
      phiprof::stop("Vmesh prefetch");
      phiprof::start("Gather moved blocks");
      for (vmesh::LocalID block_index = nBlocksAfterAdjust; block_index < nBlocksBeforeAdjust; ++block_index) {
         const vmesh::GlobalID moveCandidate = populations[popID].vmesh->getGlobalID(block_index);
         auto it = std::find(BlocksToRemove->begin(), BlocksToRemove->end(), moveCandidate);
         if (it == BlocksToRemove->end()) {
            BlocksToMove->push_back(moveCandidate);
         }
      }
      phiprof::stop("Gather moved blocks");

      // Grow the vectors, if necessary
      if (nBlocksAfterAdjust > nBlocksBeforeAdjust) {
         phiprof::start("CUDA modify vmesh and VBC size (pre)");
         populations[popID].vmesh->setNewSize(nBlocksAfterAdjust);
         populations[popID].blockContainer->setSize(nBlocksAfterAdjust);
         phiprof::stop("CUDA modify vmesh and VBC size (pre)");
      }

      const int nToAdd = BlocksToAdd->size();
      const int nToRemove = BlocksToRemove->size();
      const int nBlocks = nToAdd > nToRemove ? nToAdd : nToRemove;
      const int nCudaBlocks = nBlocks > CUDABLOCKS ? CUDABLOCKS : nBlocks;

      phiprof::start("Block lists prefetch");
      BlocksToAdd->optimizeGPU();
      BlocksToRemove->optimizeGPU();
      BlocksToMove->optimizeGPU();
      phiprof::stop("Block lists prefetch");
      phiprof::start("Vmesh and VBC lists prefetch");
      populations[popID].vmesh->dev_prefetchDevice();
      populations[popID].blockContainer->dev_prefetchDevice();
      phiprof::stop("Vmesh and VBC lists prefetch");

      const uint thread_id = omp_get_thread_num();
      cudaStream_t stream = cuda_getStream();

      Realf host_rhoLossAdjust = 0;
      vmesh::LocalID *dev_AddVectorIndex;
      vmesh::LocalID *dev_MoveVectorIndex;
      vmesh::LocalID host_AddVectorIndex = 0;
      vmesh::LocalID host_MoveVectorIndex = 0;
      HANDLE_ERROR( cudaMallocAsync((void**)&dev_AddVectorIndex, sizeof(vmesh::LocalID), stream) );
      HANDLE_ERROR( cudaMallocAsync((void**)&dev_MoveVectorIndex, sizeof(vmesh::LocalID), stream) );
      HANDLE_ERROR( cudaMemcpyAsync(returnRealf[thread_id], &host_rhoLossAdjust, sizeof(Realf), cudaMemcpyHostToDevice, stream) );
      HANDLE_ERROR( cudaMemcpyAsync(dev_AddVectorIndex, &host_AddVectorIndex, sizeof(vmesh::LocalID), cudaMemcpyHostToDevice, stream) );
      HANDLE_ERROR( cudaMemcpyAsync(dev_MoveVectorIndex, &host_MoveVectorIndex, sizeof(vmesh::LocalID), cudaMemcpyHostToDevice, stream) );

      phiprof::start("CUDA add and remove blocks kernel");

      const vmesh::LocalID nToMove = BlocksToMove->size();

      if (nCudaBlocks>0) {
         //HANDLE_ERROR( cudaStreamSynchronize(stream) );
         dim3 block(WID,WID,WID);
         // Third argument specifies the number of bytes in *shared memory* that is
         // dynamically allocated per block for this call in addition to the statically allocated memory.
         //update_velocity_blocks_kernel<<<nCudaBlocks, block, WID3*sizeof(Realf), stream>>> (
         update_velocity_blocks_kernel<<<nCudaBlocks, block, 0, stream>>> (
            populations[popID].vmesh,
            populations[popID].blockContainer,
            BlocksToAdd,
            BlocksToRemove,
            BlocksToMove,
            nBlocksBeforeAdjust,
            nBlocksAfterAdjust,
            dev_AddVectorIndex,
            dev_MoveVectorIndex,
            returnRealf[thread_id]
            );
      }
      HANDLE_ERROR( cudaStreamSynchronize(stream) );
      phiprof::stop("CUDA add and remove blocks kernel");

      HANDLE_ERROR( cudaFreeAsync(dev_AddVectorIndex, stream) );
      HANDLE_ERROR( cudaFreeAsync(dev_MoveVectorIndex, stream) );
      HANDLE_ERROR( cudaMemcpyAsync(&host_rhoLossAdjust, returnRealf[thread_id], sizeof(Realf), cudaMemcpyDeviceToHost, stream) );
      HANDLE_ERROR( cudaStreamSynchronize(stream) );
      this->populations[popID].RHOLOSSADJUST += host_rhoLossAdjust;

      // Shrink the vectors, if necessary
      if (nBlocksAfterAdjust < nBlocksBeforeAdjust) {
         phiprof::start("CUDA modify vmesh and VBC size (post)");
         // do we need prefetches here?
         populations[popID].vmesh->setNewSize(nBlocksAfterAdjust);
         populations[popID].blockContainer->setSize(nBlocksAfterAdjust);
         populations[popID].blockContainer->dev_Allocate(nBlocksAfterAdjust);
         phiprof::stop("CUDA modify vmesh and VBC size (post)");
      }

      // DEBUG output after kernel
      #ifdef DEBUG_SPATIAL_CELL
      phiprof::start("Vmesh and VBC debug output");
      populations[popID].vmesh->dev_prefetchHost();
      const vmesh::LocalID nAll = populations[popID].vmesh->size();
      printf("after kernel, size is %d should be %d\n",nAll,nBlocksAfterAdjust);
      for (vmesh::LocalID m=0; m<nAll; ++m) {
         const vmesh::GlobalID GIDs = populations[popID].vmesh->getGlobalID(m);
         const vmesh::LocalID LIDs = populations[popID].vmesh->getLocalID(GIDs);
         printf("LID %d GID-solved %d LID-solved %d\n",m,GIDs,LIDs);
      }
      populations[popID].vmesh->dev_prefetchDevice();
      phiprof::start("Vmesh and VBC debug output");
      #endif
   phiprof::stop("CUDA add and remove blocks");
   }

   void SpatialCell::adjustSingleCellVelocityBlocks(const uint popID, bool doDeleteEmpty) {
      #ifdef DEBUG_SPATIAL_CELL
      if (popID >= populations.size()) {
         std::cerr << "ERROR, popID " << popID << " exceeds populations.size() " << populations.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
         exit(1);
      }
      #endif

      //neighbor_ptrs is empty as we do not have any consistent
      //data in neighbours yet, adjustments done only based on velocity
      //space. TODO: should this delete blocks or not? Now not
      std::vector<SpatialCell*> neighbor_ptrs;
      update_velocity_block_content_lists(popID);
      adjust_velocity_blocks(neighbor_ptrs,popID,doDeleteEmpty);
   }

   /** Get maximum translation timestep for the given species.
    * @param popID ID of the particle species.
    * @return Maximum timestep calculated by the Vlasov translation.*/
   const Real& SpatialCell::get_max_r_dt(const uint popID) const {
      #ifdef DEBUG_SPATIAL_CELL
      if (popID >= populations.size()) {
         std::cerr << "ERROR, popID " << popID << " exceeds populations.size() " << populations.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
         exit(1);
      }
      #endif

      return populations[popID].max_dt[species::MAXRDT];
   }

   /** Get maximum acceleration timestep for the given species.
    * @param popID ID of the particle species.
    * @return Maximum timestep calculated by Vlasov acceleration.*/
   const Real& SpatialCell::get_max_v_dt(const uint popID) const {
      #ifdef DEBUG_SPATIAL_CELL
      if (popID >= populations.size()) {
         std::cerr << "ERROR, popID " << popID << " exceeds populations.size() " << populations.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
         exit(1);
      }
      #endif

      return populations[popID].max_dt[species::MAXVDT];
   }

   /** Get MPI datatype for sending the cell data.
    * @param cellID Spatial cell (dccrg) ID.
    * @param sender_rank Rank of the MPI process sending data from this cell.
    * @param receiver_rank Rank of the MPI process receiving data to this cell.
    * @param receiving If true, this process is receiving data.
    * @param neighborhood Neighborhood ID.
    * @return MPI datatype that transfers the requested data.*/
   std::tuple<void*, int, MPI_Datatype> SpatialCell::get_mpi_datatype(
                                                                      const CellID cellID,
                                                                      const int sender_rank,
                                                                      const int receiver_rank,
                                                                      const bool receiving,
                                                                      const int neighborhood
      ) {

      std::vector<MPI_Aint> displacements;
      std::vector<int> block_lengths;

      // create datatype for actual data if we are in the first two
      // layers around a boundary, or if we send for the whole system
      // in AMR translation, only send the necessary cells
      if (this->mpiTransferEnabled && ((SpatialCell::mpiTransferAtSysBoundaries==false && SpatialCell::mpiTransferInAMRTranslation==false) ||
                                       (SpatialCell::mpiTransferAtSysBoundaries==true && (this->sysBoundaryLayer ==1 || this->sysBoundaryLayer ==2)) ||
                                       (SpatialCell::mpiTransferInAMRTranslation==true &&
                                        this->parameters[CellParams::AMR_TRANSLATE_COMM_X+SpatialCell::mpiTransferXYZTranslation]==true ))) {

         //add data to send/recv to displacement and block length lists
         if ((SpatialCell::mpi_transfer_type & Transfer::VEL_BLOCK_LIST_STAGE1) != 0) {
            //first copy values in case this is the send operation
            populations[activePopID].N_blocks = populations[activePopID].blockContainer->size();

            // send velocity block list size
            displacements.push_back((uint8_t*) &(populations[activePopID].N_blocks) - (uint8_t*) this);
            block_lengths.push_back(sizeof(vmesh::LocalID));
         }

         if ((SpatialCell::mpi_transfer_type & Transfer::VEL_BLOCK_LIST_STAGE2) != 0) {
            // STAGE1 should have been done, otherwise we have problems...
            if (receiving) {
               //mpi_number_of_blocks transferred earlier
               populations[activePopID].vmesh->setNewSize(populations[activePopID].N_blocks);
            } else {
                //resize to correct size (it will avoid reallocation if it is big enough, I assume)
                populations[activePopID].N_blocks = populations[activePopID].blockContainer->size();
            }

            // send velocity block list
            displacements.push_back((uint8_t*) &(populations[activePopID].vmesh->getGrid()[0]) - (uint8_t*) this);
            block_lengths.push_back(sizeof(vmesh::GlobalID) * populations[activePopID].vmesh->size());
         }

         if ((SpatialCell::mpi_transfer_type & Transfer::VEL_BLOCK_WITH_CONTENT_STAGE1) !=0) {
            //Communicate size of list so that buffers can be allocated on receiving side
            if (!receiving) this->velocity_block_with_content_list_size = this->velocity_block_with_content_list->size();
            displacements.push_back((uint8_t*) &(this->velocity_block_with_content_list_size) - (uint8_t*) this);
            block_lengths.push_back(sizeof(vmesh::LocalID));
         }
         if ((SpatialCell::mpi_transfer_type & Transfer::VEL_BLOCK_WITH_CONTENT_STAGE2) !=0) {
            if (receiving) {
               this->velocity_block_with_content_list->resize(this->velocity_block_with_content_list_size);
            }

            //velocity_block_with_content_list_size should first be updated, before this can be done (STAGE1)
            //displacements.push_back((uint8_t*) &(this->velocity_block_with_content_list->at(0)) - (uint8_t*) this);
            displacements.push_back((uint8_t*) this->velocity_block_with_content_list->data() - (uint8_t*) this);
            block_lengths.push_back(sizeof(vmesh::GlobalID)*this->velocity_block_with_content_list_size);
         }

         if ((SpatialCell::mpi_transfer_type & Transfer::VEL_BLOCK_DATA) !=0) {
            displacements.push_back((uint8_t*) get_data(activePopID) - (uint8_t*) this);
            block_lengths.push_back(sizeof(Realf) * WID3 * populations[activePopID].blockContainer->size());
         }

         if ((SpatialCell::mpi_transfer_type & Transfer::NEIGHBOR_VEL_BLOCK_DATA) != 0) {
            /*We are actually transferring the data of a
            * neighbor. The values of neighbor_block_data
            * and neighbor_number_of_blocks should be set in
            * solver.*/

            // Send this data only to ranks that contain face neighbors
            // this->neighbor_number_of_blocks has been initialized to 0, on other ranks it can stay that way.
            const set<int>& ranks = this->face_neighbor_ranks[neighborhood];
            if ( P::amrMaxSpatialRefLevel == 0 || receiving || ranks.find(receiver_rank) != ranks.end()) {

               for ( int i = 0; i < MAX_NEIGHBORS_PER_DIM; ++i) {
                  displacements.push_back((uint8_t*) this->neighbor_block_data[i] - (uint8_t*) this);
                  block_lengths.push_back(sizeof(Realf) * WID3 * this->neighbor_number_of_blocks[i]);
               }

            }
         }

         // send  spatial cell parameters
         if ((SpatialCell::mpi_transfer_type & Transfer::CELL_PARAMETERS)!=0){
            displacements.push_back((uint8_t*) &(this->parameters[0]) - (uint8_t*) this);
            block_lengths.push_back(sizeof(Real) * CellParams::N_SPATIAL_CELL_PARAMS);
         }

         // send spatial cell dimensions and coordinates
         if ((SpatialCell::mpi_transfer_type & Transfer::CELL_DIMENSIONS)!=0){
            displacements.push_back((uint8_t*) &(this->parameters[CellParams::XCRD]) - (uint8_t*) this);
            block_lengths.push_back(sizeof(Real) * 6);
         }

         // send  BGBXVOL BGBYVOL BGBZVOL PERBXVOL PERBYVOL PERBZVOL
         if ((SpatialCell::mpi_transfer_type & Transfer::CELL_BVOL)!=0){
            displacements.push_back((uint8_t*) &(this->parameters[CellParams::BGBXVOL]) - (uint8_t*) this);
            block_lengths.push_back(sizeof(Real) * 6);
         }

         // send RHOM, VX, VY, VZ
         if ((SpatialCell::mpi_transfer_type & Transfer::CELL_RHOM_V)!=0){
            displacements.push_back((uint8_t*) &(this->parameters[CellParams::RHOM]) - (uint8_t*) this);
            block_lengths.push_back(sizeof(Real) * 4);
         }

         // send RHOM_DT2, VX_DT2, VY_DT2, VZ_DT2
         if ((SpatialCell::mpi_transfer_type & Transfer::CELL_RHOMDT2_VDT2)!=0){
            displacements.push_back((uint8_t*) &(this->parameters[CellParams::RHOM_DT2]) - (uint8_t*) this);
            block_lengths.push_back(sizeof(Real) * 4);
         }

         // send RHOQ
         if ((SpatialCell::mpi_transfer_type & Transfer::CELL_RHOQ)!=0){
            displacements.push_back((uint8_t*) &(this->parameters[CellParams::RHOQ]) - (uint8_t*) this);
            block_lengths.push_back(sizeof(Real));
         }

         // send RHOQ_DT2
         if ((SpatialCell::mpi_transfer_type & Transfer::CELL_RHOQDT2)!=0){
            displacements.push_back((uint8_t*) &(this->parameters[CellParams::RHOQ_DT2]) - (uint8_t*) this);
            block_lengths.push_back(sizeof(Real));
         }

         // send  spatial cell BVOL derivatives
         if ((SpatialCell::mpi_transfer_type & Transfer::CELL_BVOL_DERIVATIVES)!=0){
            displacements.push_back((uint8_t*) &(this->derivativesBVOL[0]) - (uint8_t*) this);
            block_lengths.push_back(sizeof(Real) * bvolderivatives::N_BVOL_DERIVATIVES);
         }

         if ((SpatialCell::mpi_transfer_type & Transfer::CELL_IOLOCALCELLID)!=0){
            displacements.push_back((uint8_t*) &(this->ioLocalCellId) - (uint8_t*) this);
            block_lengths.push_back(sizeof(uint64_t));
         }

         // send electron pressure gradient term components
         if ((SpatialCell::mpi_transfer_type & Transfer::CELL_GRADPE_TERM)!=0){
            displacements.push_back((uint8_t*) &(this->parameters[CellParams::EXGRADPE]) - (uint8_t*) this);
            block_lengths.push_back(sizeof(Real) * 3);
         }


         // send P tensor diagonal components
         if ((SpatialCell::mpi_transfer_type & Transfer::CELL_P)!=0){
            displacements.push_back((uint8_t*) &(this->parameters[CellParams::P_11]) - (uint8_t*) this);
            block_lengths.push_back(sizeof(Real) * 3);
         }

         if ((SpatialCell::mpi_transfer_type & Transfer::CELL_PDT2)!=0){
            displacements.push_back((uint8_t*) &(this->parameters[CellParams::P_11_DT2]) - (uint8_t*) this);
            block_lengths.push_back(sizeof(Real) * 3);
         }

         // send  sysBoundaryFlag
         if ((SpatialCell::mpi_transfer_type & Transfer::CELL_SYSBOUNDARYFLAG)!=0){
            displacements.push_back((uint8_t*) &(this->sysBoundaryFlag) - (uint8_t*) this);
            block_lengths.push_back(sizeof(uint));
            displacements.push_back((uint8_t*) &(this->sysBoundaryLayer) - (uint8_t*) this);
            block_lengths.push_back(sizeof(uint));
         }

         if ((SpatialCell::mpi_transfer_type & Transfer::VEL_BLOCK_PARAMETERS) !=0) {
            displacements.push_back((uint8_t*) get_block_parameters(activePopID) - (uint8_t*) this);
            block_lengths.push_back(sizeof(Real) * size(activePopID) * BlockParams::N_VELOCITY_BLOCK_PARAMS);
         }
         // Copy particle species metadata
         if ((SpatialCell::mpi_transfer_type & Transfer::POP_METADATA) != 0) {
            for (uint popID=0; popID<populations.size(); ++popID) {
               displacements.push_back((uint8_t*) &(populations[popID].RHO) - (uint8_t*)this);
               block_lengths.push_back(offsetof(spatial_cell::Population, N_blocks));
            }
         }
         // Copy random number generator state variables
         //if ((SpatialCell::mpi_transfer_type & Transfer::RANDOMGEN) != 0) {
         //   displacements.push_back((uint8_t*)get_rng_state_buffer() - (uint8_t*)this);
         //   block_lengths.push_back(256/8);
         //   displacements.push_back((uint8_t*)get_rng_data_buffer() - (uint8_t*)this);
         //   block_lengths.push_back(sizeof(random_data));
         //}
      }

      void* address = this;
      int count;
      MPI_Datatype datatype;

      if (displacements.size() > 0) {
         count = 1;
         MPI_Type_create_hindexed(
            displacements.size(),
            &block_lengths[0],
            &displacements[0],
            MPI_BYTE,
            &datatype
         );
      } else {
         count = 0;
         datatype = MPI_BYTE;
      }

      const bool printMpiDatatype = false;
      if(printMpiDatatype) {
         int mpiSize;
         int myRank;
         MPI_Type_size(datatype,&mpiSize);
         MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
         cout << myRank << " get_mpi_datatype: " << cellID << " " << sender_rank << " " << receiver_rank << " " << mpiSize << ", Nblocks = " << populations[activePopID].N_blocks << ", nbr Nblocks =";
         for (uint i = 0; i < MAX_NEIGHBORS_PER_DIM; ++i) {
            const set<int>& ranks = this->face_neighbor_ranks[neighborhood];
            if ( receiving || ranks.find(receiver_rank) != ranks.end()) {
               cout << " " << this->neighbor_number_of_blocks[i];
            } else {
               cout << " " << 0;
            }
         }
         cout << " face_neighbor_ranks =";
         for (const auto& rank : this->face_neighbor_ranks[neighborhood]) {
            cout << " " << rank;
         }
         cout << endl;
      }

      return std::make_tuple(address,count,datatype);
   }

  /**< Minimum value of distribution function in any phase space cell
    * of a velocity block for the block to be considered to have content.
    * @param popID ID of the particle species.
    * @return Sparse min value for this species.*/
   Real SpatialCell::getVelocityBlockMinValue(const uint popID) const {
      return populations[popID].velocityBlockMinValue;
   }

   /** Prepares this spatial cell to receive the velocity grid over MPI.
    * At this stage we have received a new block list over MPI into
    * mpi_velocity_block_list, but the rest of the cell structures
    * have not been adapted to this new list. Here we re-initialize
    * the cell with empty blocks based on the new list.*/
   void SpatialCell::prepare_to_receive_blocks(const uint popID) {
      populations[popID].vmesh->setGrid();
      populations[popID].blockContainer->setSize(populations[popID].vmesh->size());

      Real* parameters = get_block_parameters(popID);

      // Set velocity block parameters:
      for (vmesh::LocalID blockLID=0; blockLID<size(popID); ++blockLID) {
         const vmesh::GlobalID blockGID = get_velocity_block_global_id(blockLID,popID);
         parameters[BlockParams::VXCRD] = get_velocity_block_vx_min(popID,blockGID);
         parameters[BlockParams::VYCRD] = get_velocity_block_vy_min(popID,blockGID);
         parameters[BlockParams::VZCRD] = get_velocity_block_vz_min(popID,blockGID);
         populations[popID].vmesh->getCellSize(blockGID,&(parameters[BlockParams::DVX]));
         parameters += BlockParams::N_VELOCITY_BLOCK_PARAMS;
      }
   }

   /** Set the particle species SpatialCell should use in functions that
    * use the velocity mesh.
    * @param popID Population ID.
    * @return If true, the new species is in use.*/
   bool SpatialCell::setCommunicatedSpecies(const uint popID) {
      #ifdef DEBUG_SPATIAL_CELL
      if (popID >= getObjectWrapper().particleSpecies.size()) {
         std::cerr << "ERROR, popID " << popID << " exceeds species.size() " << getObjectWrapper().particleSpecies.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
         exit(1);
      }
      #endif

      activePopID = popID;
      return true;
   }

   /** Set maximum translation timestep for a particle species.
    * This function is called during Vlasov translation.
    * @param popID ID of the particle species.
    * @param value New maximum timestep.*/
   void SpatialCell::set_max_r_dt(const uint popID,const Real& value) {
      #ifdef DEBUG_SPATIAL_CELL
      if (popID >= populations.size()) {
         std::cerr << "ERROR, popID " << popID << " exceeds populations.size() " << populations.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
         exit(1);
      }
      #endif

      populations[popID].max_dt[species::MAXRDT] = value;
   }

   /** Set maximum acceleration timestep for a particle species.
    * This function is called during Vlasov acceleration.
    * @param popID ID of the particle species.
    * @param value New maximum timestep.*/
   void SpatialCell::set_max_v_dt(const uint popID,const Real& value) {
      #ifdef DEBUG_SPATIAL_CELL
      if (popID >= populations.size()) {
         std::cerr << "ERROR, popID " << popID << " exceeds populations.size() " << populations.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
         exit(1);
      }
      #endif

      populations[popID].max_dt[species::MAXVDT] = value;
   }

   /**  Purges extra capacity from block vectors. It sets size to
    * num_blocks * block_allocation_factor (if capacity greater than this),
    * and also forces capacity to this new smaller value.
    * @return True on success.*/
   bool SpatialCell::shrink_to_fit() {
      bool success = true;
      return success;

      for (size_t p=0; p<populations.size(); ++p) {
         const uint64_t amount
            = 2 + populations[p].blockContainer->size()
            * populations[p].blockContainer->getBlockAllocationFactor();

         // Allow capacity to be a bit large than needed by number of blocks, shrink otherwise
         if (populations[p].blockContainer->capacity() > amount )
            if (populations[p].blockContainer->recapacitate(amount) == false) success = false;

      }
      return success;
   }

   /** Update the two lists containing blocks with content, and blocks without content.
    * @see adjustVelocityBlocks */
   void SpatialCell::update_velocity_block_content_lists(const uint popID) {
      #ifdef DEBUG_SPATIAL_CELL
      if (popID >= populations.size()) {
         std::cerr << "ERROR, popID " << popID << " exceeds populations.size() " << populations.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
         exit(1);
      }
      #endif

      phiprof::start("CUDA update spatial cell block lists");

      cudaStream_t stream = cuda_getStream();
      phiprof::start("VB content list prefetches and allocations");
      velocity_block_with_content_list->clear();
      velocity_block_with_no_content_list->clear();
      velocity_block_with_content_list->reserve(populations[popID].vmesh->size());
      velocity_block_with_no_content_list->reserve(populations[popID].vmesh->size());
      velocity_block_with_content_list->optimizeGPU(stream);
      velocity_block_with_no_content_list->optimizeGPU(stream);
      phiprof::stop("VB content list prefetches and allocations");

      const Real velocity_block_min_value = getVelocityBlockMinValue(popID);

      phiprof::start("CUDA update spatial cell block lists kernel");
      dim3 block(WID,WID,WID);
      // Third argument specifies the number of bytes in *shared memory* that is
      // dynamically allocated per block for this call in addition to the statically allocated memory.
      update_velocity_block_content_lists_kernel<<<CUDABLOCKS, block, WID3*sizeof(bool), stream>>> (
         populations[popID].vmesh,
         populations[popID].blockContainer,
         velocity_block_with_content_list,
         velocity_block_with_no_content_list,
         velocity_block_min_value
         );

      HANDLE_ERROR( cudaStreamSynchronize(stream) );
      phiprof::stop("CUDA update spatial cell block lists kernel");
      phiprof::start("VB content list prefetches");
      velocity_block_with_content_list->optimizeCPU(stream);
      velocity_block_with_no_content_list->optimizeCPU(stream);
      phiprof::stop("VB content list prefetches");
      phiprof::stop("CUDA update spatial cell block lists");
   }

   void SpatialCell::printMeshSizes() {
      cerr << "SC::printMeshSizes:" << endl;
      for (size_t p=0; p<populations.size(); ++p) {
         cerr << "\t pop " << p << " " << populations[p].vmesh->size() << ' ' << populations[p].blockContainer->size() << endl;
      }
   }

   /** Updates minValue based on algorithm value from parameters (see parameters.cpp).
    * @param popID ID of the particle species.*/
   void SpatialCell::updateSparseMinValue(const uint popID) {

      species::Species& population = getObjectWrapper().particleSpecies[popID];

      if ( population.sparseDynamicAlgorithm == 1 || population.sparseDynamicAlgorithm == 2 ) {
         // Linear algorithm for the minValue: y=kx+b
         const Real k = (population.sparseDynamicMinValue2 - population.sparseDynamicMinValue1) / (population.sparseDynamicBulkValue2 - population.sparseDynamicBulkValue1);
         const Real b = population.sparseDynamicMinValue1 - k * population.sparseDynamicBulkValue1;
         Real x;
         if ( population.sparseDynamicAlgorithm == 1 ) {
            x = this->populations[popID].RHO;
         } else {
            x = this->get_number_of_velocity_blocks(popID);
         }
         const Real newMinValue = k*x+b;
         if( newMinValue < population.sparseDynamicMinValue1 ) { // Compare against the min minValue
            populations[popID].velocityBlockMinValue = population.sparseDynamicMinValue1;
         } else if( newMinValue > population.sparseDynamicMinValue2 ) { // Compare against the max minValue
            populations[popID].velocityBlockMinValue = population.sparseDynamicMinValue2;
         } else {
            populations[popID].velocityBlockMinValue = newMinValue;
         }
         return;
      } else {
         populations[popID].velocityBlockMinValue = getObjectWrapper().particleSpecies[popID].sparseMinValue;
         return;
      }
   }

} // namespace spatial_cell
