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

#include <unordered_set>

#include "spatial_cell_wrapper.hpp"
#include "../arch/gpu_base.hpp"
#include "../object_wrapper.h"
#include "../velocity_mesh_parameters.h"

#include "spatial_cell_gpu_kernels.hpp"

// INIT_VMESH_SIZE and INIT_MAP_SIZE defined in arch/gpu_base.hpp

using namespace std;

// Certain addition lists are used in acceleration, and can require a larger allocation
// in case very many blocks are added at once.
const static uint acc_reserve_multiplier = 3;

namespace spatial_cell {
   int SpatialCell::activePopID = 0;
   uint64_t SpatialCell::mpi_transfer_type = 0;
   bool SpatialCell::mpiTransferAtSysBoundaries = false;

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
         populations[popID].vmesh->gpu_prefetchDevice();
         populations[popID].blockContainer->gpu_prefetchDevice();
         populations[popID].Upload();
         populations[popID].velocityBlockMinValue = spec.sparseMinValue;
         populations[popID].N_blocks = 0;
      }

      // SplitVectors and hashmaps via pointers for unified memory

      // create in host instead of unified memory, upload device copy
      void *buf0 = malloc(sizeof(split::SplitVector<vmesh::GlobalID>));
      velocity_block_with_content_list = ::new (buf0) split::SplitVector<vmesh::GlobalID>(INIT_VMESH_SIZE);
      //velocity_block_with_content_list = new split::SplitVector<vmesh::GlobalID>(INIT_VMESH_SIZE);
      velocity_block_with_content_list->clear();
      velocity_block_with_content_list_size=0;
      velocity_block_with_content_list_capacity=INIT_VMESH_SIZE;
      dev_velocity_block_with_content_list = velocity_block_with_content_list->upload<true>();

      // create in host instead of unified memory, upload device copy
      // velocity_block_with_content_map = new Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>(7);
      // velocity_block_with_no_content_map = new Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>(7);
      void *buf1 = malloc(sizeof(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>));
      void *buf2 = malloc(sizeof(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>));
      velocity_block_with_content_map = ::new (buf1) Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>(INIT_MAP_SIZE);
      velocity_block_with_no_content_map = ::new (buf2) Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>(INIT_MAP_SIZE);
      dev_velocity_block_with_content_map = velocity_block_with_content_map->upload<true>();
      dev_velocity_block_with_no_content_map = velocity_block_with_no_content_map->upload<true>();
      vbwcl_sizePower = INIT_MAP_SIZE;
      vbwncl_sizePower = INIT_MAP_SIZE;

      // Lists used in block adjustment
      void *buf11 = malloc(sizeof(split::SplitVector<vmesh::GlobalID>));
      void *buf12 = malloc(sizeof(split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>));
      void *buf13 = malloc(sizeof(split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>));
      void *buf14 = malloc(sizeof(split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>));
      list_with_replace_new = ::new (buf11) split::SplitVector<vmesh::GlobalID>(INIT_VMESH_SIZE*acc_reserve_multiplier);
      list_delete = ::new (buf12) split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>(INIT_VMESH_SIZE);
      list_to_replace = ::new (buf13) split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>(INIT_VMESH_SIZE);
      list_with_replace_old = ::new (buf14) split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>(INIT_VMESH_SIZE);
      // list_with_replace_new = new split::SplitVector<vmesh::GlobalID>(INIT_VMESH_SIZE*acc_reserve_multiplier);
      // list_delete = new split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>(INIT_VMESH_SIZE);
      // list_to_replace = new split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>(INIT_VMESH_SIZE);
      // list_with_replace_old = new split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>(INIT_VMESH_SIZE);
      dev_list_with_replace_new = list_with_replace_new->upload<true>();
      dev_list_delete = list_delete->upload<true>();
      dev_list_to_replace = list_to_replace->upload<true>();
      dev_list_with_replace_old = list_with_replace_old->upload<true>();
      list_with_replace_new_capacity = INIT_VMESH_SIZE*acc_reserve_multiplier;
      list_delete_capacity = INIT_VMESH_SIZE;
      list_to_replace_capacity = INIT_VMESH_SIZE;
      list_with_replace_old_capacity = INIT_VMESH_SIZE;
   }

   SpatialCell::~SpatialCell() {
      if (velocity_block_with_content_list) {
         ::delete velocity_block_with_content_list;
         velocity_block_with_content_list = 0;
      }
      if (velocity_block_with_content_map) {
         ::delete velocity_block_with_content_map;
         velocity_block_with_content_map = 0;
      }
      if (velocity_block_with_no_content_map) {
         ::delete velocity_block_with_no_content_map;
         velocity_block_with_no_content_map = 0;
      }
      velocity_block_with_content_list_size=0;
      velocity_block_with_content_list_capacity=0;
      vbwcl_sizePower = 0;
      vbwncl_sizePower = 0;

      if (list_with_replace_new) {
         ::delete list_with_replace_new;
         list_with_replace_new = 0;
      }
      if (list_delete) {
         ::delete list_delete;
         list_delete = 0;
      }
      if (list_to_replace) {
         ::delete list_to_replace;
         list_to_replace = 0;
      }
      if (list_with_replace_old) {
         ::delete list_with_replace_old;
         list_with_replace_old = 0;
      }
      list_with_replace_new_capacity = 0;
      list_delete_capacity = 0;
      list_to_replace_capacity = 0;
      list_with_replace_old_capacity = 0;
   }

   SpatialCell::SpatialCell(const SpatialCell& other) {
      std::cerr<<"Warning! Spatial Cell GPU copy constructor called. Performance may degrade."<<std::endl;
      // Note: DCCRG should not call this method ever, as far as we know. Thus untested.
      const uint reserveSize = other.velocity_block_with_content_list_capacity;

      // create in host instead of unified memory, upload device copy
      void *buf0 = malloc(sizeof(split::SplitVector<vmesh::GlobalID>));
      velocity_block_with_content_list = ::new (buf0) split::SplitVector<vmesh::GlobalID>(reserveSize);
      //velocity_block_with_content_list = new split::SplitVector<vmesh::GlobalID>(reserveSize);
      velocity_block_with_content_list->clear();
      velocity_block_with_content_list_size = 0;
      velocity_block_with_content_list_capacity = reserveSize;
      dev_velocity_block_with_content_list = velocity_block_with_content_list->upload<true>();

      // create in host instead of unified memory, upload device copy
      void *buf1 = malloc(sizeof(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>));
      void *buf2 = malloc(sizeof(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>));
      velocity_block_with_content_map = ::new (buf1) Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>(other.vbwcl_sizePower);
      velocity_block_with_no_content_map = ::new (buf2) Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>(other.vbwncl_sizePower);
      dev_velocity_block_with_content_map = velocity_block_with_content_map->upload<true>();
      dev_velocity_block_with_no_content_map = velocity_block_with_no_content_map->upload<true>();
      vbwcl_sizePower = other.vbwcl_sizePower;
      vbwncl_sizePower = other.vbwncl_sizePower;

      // Lists used in block adjustment
      void *buf11 = malloc(sizeof(split::SplitVector<vmesh::GlobalID>));
      void *buf12 = malloc(sizeof(split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>));
      void *buf13 = malloc(sizeof(split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>));
      void *buf14 = malloc(sizeof(split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>));
      list_with_replace_new = ::new (buf11) split::SplitVector<vmesh::GlobalID>(other.list_with_replace_new_capacity);
      list_delete = ::new (buf12) split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>(other.list_delete_capacity);
      list_to_replace = ::new (buf13) split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>(other.list_to_replace_capacity);
      list_with_replace_old = ::new (buf14) split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>(other.list_with_replace_old_capacity);
      // list_with_replace_new = new split::SplitVector<vmesh::GlobalID>(other.list_with_replace_new_capacity);
      // list_delete = new split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>(other.list_delete_capacity);
      // list_to_replace = new split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>(other.list_to_replace_capacity);
      // list_with_replace_old = new split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>(other.list_with_replace_old_capacity);
      dev_list_with_replace_new = list_with_replace_new->upload<true>();
      dev_list_delete = list_delete->upload<true>();
      dev_list_to_replace = list_to_replace->upload<true>();
      dev_list_with_replace_old = list_with_replace_old->upload<true>();
      list_with_replace_new_capacity = other.list_with_replace_new_capacity;
      list_delete_capacity = other.list_delete_capacity;
      list_to_replace_capacity = other.list_to_replace_capacity;
      list_with_replace_old_capacity = other.list_with_replace_old_capacity;

      // Member variables
      ioLocalCellId = other.ioLocalCellId;
      sysBoundaryFlag = other.sysBoundaryFlag;
      sysBoundaryLayer = other.sysBoundaryLayer;
      sysBoundaryLayerNew = other.sysBoundaryLayerNew;
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
         neighbor_block_data[i] = 0;
         neighbor_number_of_blocks[i] = 0;
      }
      face_neighbor_ranks.clear();
      // for (unsigned int i=0; i<MAX_NEIGHBORS_PER_DIM; ++i) {
      //    neighbor_block_data[i] = other.neighbor_block_data[i];
      //    neighbor_number_of_blocks[i] = other.neighbor_number_of_blocks[i];
      // }
      // if (other.face_neighbor_ranks.size()>0) {
      //    face_neighbor_ranks = std::map<int,std::set<int>>(other.face_neighbor_ranks);
      // }
      if (other.populations.size()>0) {
         populations = std::vector<spatial_cell::Population>(other.populations);
      }
   }

   const SpatialCell& SpatialCell::operator=(const SpatialCell& other) {
      // Used for refining spatial cells
      velocity_block_with_content_list_capacity = other.velocity_block_with_content_list_capacity;
      velocity_block_with_content_list->clear();
      velocity_block_with_content_list->reserve(velocity_block_with_content_list_capacity);
      velocity_block_with_content_list_size = 0;
      dev_velocity_block_with_content_list = velocity_block_with_content_list->upload<false>();

      gpuStream_t stream = gpu_getStream();

      if (vbwcl_sizePower < other.vbwcl_sizePower) {
         vbwcl_sizePower = other.vbwcl_sizePower;
         ::delete velocity_block_with_content_map;
         void *buf1 = malloc(sizeof(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>));
         velocity_block_with_content_map = ::new (buf1)Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>(vbwcl_sizePower);
         dev_velocity_block_with_content_map = velocity_block_with_content_map->upload<true>(stream);
      } else {
         velocity_block_with_content_map->clear<false>(Hashinator::targets::device,stream,std::pow(2,vbwcl_sizePower));
      }
      if (vbwncl_sizePower < other.vbwncl_sizePower) {
         vbwncl_sizePower = other.vbwncl_sizePower;
         ::delete velocity_block_with_no_content_map;
         void *buf2 = malloc(sizeof(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>));
         velocity_block_with_no_content_map = ::new (buf2) Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>(vbwncl_sizePower);
         dev_velocity_block_with_no_content_map = velocity_block_with_no_content_map->upload<true>(stream);
      } else {
         velocity_block_with_no_content_map->clear<false>(Hashinator::targets::device,stream,std::pow(2,vbwncl_sizePower));
      }

      list_with_replace_new_capacity = other.list_with_replace_new_capacity;
      list_delete_capacity = other.list_delete_capacity;
      list_to_replace_capacity = other.list_to_replace_capacity;
      list_with_replace_old_capacity = other.list_with_replace_old_capacity;
      list_with_replace_new->clear();
      list_delete->clear();
      list_to_replace->clear();
      list_with_replace_old->clear();

      list_with_replace_new->reserve(list_with_replace_new_capacity);
      list_delete->reserve(list_delete_capacity);
      list_to_replace->reserve(list_to_replace_capacity);
      list_with_replace_old->reserve(list_with_replace_old_capacity);
      dev_list_with_replace_new = list_with_replace_new->upload<false>();
      dev_list_delete = list_delete->upload<false>();
      dev_list_to_replace = list_to_replace->upload<false>();
      dev_list_with_replace_old = list_with_replace_old->upload<false>();

      // Member variables
      //ioLocalCellId = other.ioLocalCellId;
      sysBoundaryFlag = other.sysBoundaryFlag;
      sysBoundaryLayer = other.sysBoundaryLayer;
      sysBoundaryLayerNew = other.sysBoundaryLayerNew;
      initialized = other.initialized;
      mpiTransferEnabled = other.mpiTransferEnabled;
      for (unsigned int i=0; i<bvolderivatives::N_BVOL_DERIVATIVES; ++i) {
         derivativesBVOL[i] = other.derivativesBVOL[i]; // Will be re-calculated
      }
      for (unsigned int i=0; i<CellParams::N_SPATIAL_CELL_PARAMS; ++i) {
         parameters[i] = other.parameters[i]; // Need to be re-defined
      }
      for (unsigned int i=0; i<WID3; ++i) {
         null_block_data[i] = 0.0;
      }
      for (unsigned int i=0; i<MAX_NEIGHBORS_PER_DIM; ++i) {
         neighbor_block_data[i] = 0;
         neighbor_number_of_blocks[i] = 0;
      }
      face_neighbor_ranks.clear(); // Needs re-building after refinement
      populations = std::vector<spatial_cell::Population>(other.populations);

      return *this;
   }

   /** Sets a guidance counter so that vmesh adjustment vectors have sufficient size
    */
   void SpatialCell::setReservation(const uint popID, const vmesh::LocalID reservationsize, bool force) {
      if (force || (reservationsize > populations[popID].reservation)) {
         populations[popID].reservation = reservationsize;
      }
   }
   vmesh::LocalID SpatialCell::getReservation(const uint popID) const {
      return populations[popID].reservation;
   }
   /** Recapacitates local temporary vectors based on guidance counter
    */
   void SpatialCell::applyReservation(const uint popID) {
      const size_t reserveSize = populations[popID].reservation;
      size_t newReserve = populations[popID].reservation * BLOCK_ALLOCATION_PADDING;
      const vmesh::LocalID HashmapReqSize = ceil(log2(reserveSize))+1;
      gpuStream_t stream = gpu_getStream();
      // Now uses host-cached values
      // upload() calls include an optimizeGPU() already.

      if (velocity_block_with_content_list_capacity < reserveSize) {
         velocity_block_with_content_list->reserve(newReserve,true);
         velocity_block_with_content_list_capacity = newReserve;
         dev_velocity_block_with_content_list = velocity_block_with_content_list->upload<true>(stream);
      }
      // This one is also used in acceleration for adding new blocks to the mesh, so should have more room. 
      if (vbwcl_sizePower < HashmapReqSize+1) {
         vbwcl_sizePower = HashmapReqSize+1;
         velocity_block_with_content_map->resize(vbwcl_sizePower, Hashinator::targets::device, stream);
         // ::delete velocity_block_with_content_map;
         // void *buf = malloc(sizeof(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>));
         // velocity_block_with_content_map = ::new (buf) Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>(vbwcl_sizePower);
         // dev_velocity_block_with_content_map = velocity_block_with_content_map->upload<true>(stream);
      }
      // Here the regular size estimate should be enough.
      if (vbwncl_sizePower < HashmapReqSize) {
         vbwncl_sizePower = HashmapReqSize;
         velocity_block_with_no_content_map->resize(vbwncl_sizePower, Hashinator::targets::device, stream);
         // ::delete velocity_block_with_no_content_map;
         // void *buf = malloc(sizeof(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>));
         // velocity_block_with_no_content_map = ::new (buf) Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>(vbwncl_sizePower);
         // dev_velocity_block_with_no_content_map = velocity_block_with_no_content_map->upload<true>(stream);
      }
      // These lists are also used in acceleration, where sometimes, very many blocks may be added.
      // (Maximum possible is all existing blocks moved to a new location + 2 per column)
      // Thus, this one list needs to have larger capacity than the others..
      if (list_with_replace_new_capacity < reserveSize * acc_reserve_multiplier + 2*gpu_largest_columnCount) {
         list_with_replace_new->reserve(newReserve * acc_reserve_multiplier,true);
         list_with_replace_new_capacity = newReserve * acc_reserve_multiplier;
         dev_list_with_replace_new = list_with_replace_new->upload<true>(stream);
      }
      if (list_delete_capacity < reserveSize) {
         list_delete->reserve(newReserve,true);
         list_delete_capacity = newReserve;
         dev_list_delete = list_delete->upload<true>(stream);
      }
      if (list_to_replace_capacity < reserveSize) {
         list_to_replace->reserve(newReserve,true);
         list_to_replace_capacity = newReserve;
         dev_list_to_replace = list_to_replace->upload<true>(stream);
      }
      if (list_with_replace_old_capacity < reserveSize) {
         list_with_replace_old->reserve(newReserve,true);
         list_with_replace_old_capacity = newReserve;
         dev_list_with_replace_old = list_with_replace_old->upload<true>(stream);
      }
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

   void SpatialCell::adjust_velocity_blocks(const uint popID, bool doDeleteEmptyBlocks) {
      debug_population_check(popID);

      const uint cpuThreadID = gpu_getThread();
      const gpuStream_t stream = gpu_getStream();

      vmesh::VelocityMesh* host_vmesh    = populations[popID].vmesh;
      vmesh::VelocityMesh* dev_vmesh    = populations[popID].dev_vmesh;
      vmesh::GlobalID* _withContentData = velocity_block_with_content_list->data();

      // Evaluate velocity halo for local content blocks
      if (velocity_block_with_content_list_size>0) {
         const int addWidthV = getObjectWrapper().particleSpecies[popID].sparseBlockAddWidthV;
         if (addWidthV!=1) {
            std::cerr<<"Warning! "<<__FILE__<<":"<<__LINE__<<" Halo extent is not 1, unsupported size."<<std::endl;
         }
         // Halo of 1 in each direction adds up to 26 velocity neighbours.
         // For NVIDIA/CUDA, we dan do 26 neighbours and 32 threads per warp in a single block.
         // For AMD/HIP, we dan do 13 neighbours and 64 threads per warp in a single block, meaning two loops per cell.
         // In either case, we launch blocks equal to velocity_block_with_content_list_size
         update_velocity_halo_kernel<<<velocity_block_with_content_list_size, 26*32, 0, stream>>> (
            dev_vmesh,
            velocity_block_with_content_list_size,
            _withContentData,
            dev_velocity_block_with_content_map,
            dev_velocity_block_with_no_content_map
            );
         CHK_ERR( gpuPeekAtLastError() );
         //CHK_ERR( gpuStreamSynchronize(stream) );
      }

      /** Neighbour inclusion: GPU block adjusment now happens via batch calls.
          Single-cell block adjustments will not want to account for neighbours, so this
          code is no longer needed.
      */

      /****
      // Gather pointers and counts from neighbours
      uint neighbours_count = neighbor_ptrs.size();
      uint neighbours_blocks_count = 0;
      std::vector<vmesh::GlobalID*> neigh_vbwcls;
      std::vector<vmesh::LocalID> neigh_Nvbwcls;

      if (neighbours_count > 0) {
         for (std::vector<SpatialCell*>::const_iterator neighbor=neighbor_ptrs.begin();
              neighbor != neighbor_ptrs.end(); ++neighbor) {
            if ((*neighbor)->velocity_block_with_content_list_size > 0) {
               neigh_Nvbwcls.push_back((*neighbor)->velocity_block_with_content_list_size);
               neigh_vbwcls.push_back((*neighbor)->velocity_block_with_content_list->data());
               neighbours_blocks_count += (*neighbor)->velocity_block_with_content_list_size;
            }
         }
         neighbours_count = neigh_Nvbwcls.size(); // Only include neighbours with content.
      }

      if (neighbours_count > 0) {
         // Upload pointers and counters for neighbours
         vmesh::GlobalID** dev_neigh_vbwcls;
         vmesh::GlobalID* dev_neigh_Nvbwcls;
         CHK_ERR( gpuMallocAsync((void**)&dev_neigh_vbwcls, neighbours_count*sizeof(vmesh::GlobalID*), stream) ); // where would these be cleared?
         CHK_ERR( gpuMallocAsync((void**)&dev_neigh_Nvbwcls, neighbours_count*sizeof(vmesh::LocalID), stream) );
         CHK_ERR( gpuMemcpyAsync(dev_neigh_vbwcls, neigh_vbwcls.data(), neighbours_count*sizeof(vmesh::GlobalID*), gpuMemcpyHostToDevice, stream) );
         CHK_ERR( gpuMemcpyAsync(dev_neigh_Nvbwcls, neigh_Nvbwcls.data(), neighbours_count*sizeof(vmesh::LocalID), gpuMemcpyHostToDevice, stream) );
         // For NVIDIA/CUDA, we dan do 32 neighbour GIDs and 32 threads per warp in a single block.
         // For AMD/HIP, we dan do 16 neighbour GIDs and 64 threads per warp in a single block
         // This is handled in-kernel.
         // ceil int division
         uint launchBlocks = 1 + ((neighbours_blocks_count - 1) / WARPSPERBLOCK);
         if (launchBlocks < std::pow(2,31)) {
            update_neighbour_halo_kernel<<<launchBlocks, WARPSPERBLOCK*GPUTHREADS, 0, stream>>> (
               dev_vmesh,
               neighbours_count,
               dev_neigh_vbwcls,
               dev_neigh_Nvbwcls,
               dev_velocity_block_with_content_map,
               dev_velocity_block_with_no_content_map
               );
            CHK_ERR( gpuPeekAtLastError() );
         } else {
            // Too many launch blocks, call one by one (unlikely)
            for (uint neigh_i = 0; neigh_i < neighbours_count; neigh_i++) {
               uint launchBlocks = 1 + ((neigh_Nvbwcls[neigh_i] - 1) / WARPSPERBLOCK);
               update_neighbour_halo_kernel<<<launchBlocks, WARPSPERBLOCK*GPUTHREADS, 0, stream>>> (
                  dev_vmesh,
                  1,
                  dev_neigh_vbwcls+neigh_i,
                  dev_neigh_Nvbwcls+neigh_i,
                  dev_velocity_block_with_content_map,
                  dev_velocity_block_with_no_content_map
                  );
               CHK_ERR( gpuPeekAtLastError() );
            }
         }
         //CHK_ERR( gpuStreamSynchronize(stream) );
      }
      *****/

      /**
          Now extract vectors to be used in actual block adjustment.
          Previous kernels may have added dummy (to be added) entries to
          velocity_block_with_content_map with LID=vmesh->invalidLocalID()
          Or non-content blocks which should be retained (with correct LIDs).

          Rules used in extracting keys or elements from hashmaps
          Now these include passing pointers to GPU memory in order to evaluate
          nBlocksAfterAdjust without going via host. Pointers are copied by value.
      */
      const vmesh::GlobalID emptybucket = velocity_block_with_content_map->get_emptybucket();
      const vmesh::GlobalID tombstone   = velocity_block_with_content_map->get_tombstone();
      const vmesh::GlobalID invalidGID  = host_vmesh->invalidGlobalID();
      const vmesh::LocalID  invalidLID  = host_vmesh->invalidLocalID();

      auto rule_add = [emptybucket, tombstone, invalidGID, invalidLID]
         __device__(const Hashinator::hash_pair<vmesh::GlobalID, vmesh::LocalID>& kval) -> bool {
                         return kval.first != emptybucket &&
                            kval.first != tombstone &&
                            kval.first != invalidGID &&
                            // Required GIDs which do not yet exist in vmesh were stored in
                            // velocity_block_with_content_map with kval.second==invalidLID
                            kval.second == invalidLID; };
      velocity_block_with_content_map->extractKeysByPatternLoop(*dev_list_with_replace_new, rule_add, stream);

      if (doDeleteEmptyBlocks) {
         Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID> *vbwncm = dev_velocity_block_with_no_content_map;
         split::SplitVector<vmesh::GlobalID> *d_list_add = dev_list_with_replace_new;

         auto rule_delete_move = [emptybucket, tombstone, vbwncm, d_list_add, dev_vmesh, invalidGID, invalidLID]
            __device__(const Hashinator::hash_pair<vmesh::GlobalID, vmesh::LocalID>& kval) -> bool {
                                    const vmesh::LocalID nBlocksAfterAdjust1 = dev_vmesh->size()
                                       + d_list_add->size() - vbwncm->size();
                                    return kval.first != emptybucket &&
                                       kval.first != tombstone &&
                                       kval.first != invalidGID &&
                                       kval.second >= nBlocksAfterAdjust1 &&
                                       kval.second != invalidLID; };
         auto rule_to_replace = [emptybucket, tombstone, vbwncm, d_list_add, dev_vmesh, invalidGID, invalidLID]
            __device__(const Hashinator::hash_pair<vmesh::GlobalID, vmesh::LocalID>& kval) -> bool {
                                   const vmesh::LocalID nBlocksAfterAdjust2 = dev_vmesh->size()
                                      + d_list_add->size() - vbwncm->size();
                                   return kval.first != emptybucket &&
                                      kval.first != tombstone &&
                                      kval.first != invalidGID &&
                                      kval.second < nBlocksAfterAdjust2 &&
                                                    kval.second != invalidGID; };

         velocity_block_with_content_map->extractPatternLoop(*dev_list_with_replace_old, rule_delete_move, stream);
         velocity_block_with_no_content_map->extractPatternLoop(*dev_list_delete, rule_delete_move, stream);
         velocity_block_with_no_content_map->extractPatternLoop(*dev_list_to_replace, rule_to_replace, stream);
      } else {
         list_with_replace_old->clear();
         list_delete->clear();
         list_to_replace->clear();
      }
      // Note:list_with_replace_new contains both new GIDs to use for replacements and new GIDs to place at end of vmesh

      // Actual adjustment calling happens in separate function as it is also called from within acceleration
      vmesh::LocalID nBlocksAfterAdjust = adjust_velocity_blocks_caller(popID);

      // Perform hashmap cleanup here (instead of at acceleration mid-steps)
      populations[popID].vmesh->gpu_cleanHashMap(stream);
      populations[popID].Upload();

      #ifdef DEBUG_SPATIAL_CELL
      const size_t vmeshSize = (populations[popID].vmesh)->size();
      const size_t vbcSize = (populations[popID].blockContainer)->size();
      if (vmeshSize != vbcSize) {
         printf("ERROR: population vmesh %zu and blockcontainer %zu sizes do not match!\n",vmeshSize,vbcSize);
      }
      #endif
      #ifdef DEBUG_VLASIATOR
      // This is a bit extreme
      if (!populations[popID].vmesh->check()) {
         printf("ERROR in vmesh check: %s at %d\n",__FILE__,__LINE__);
      }
      #endif
   }

   /**
      Call GPU kernel with all necessary information for creation and deletion of blocks.
   **/
   vmesh::LocalID SpatialCell::adjust_velocity_blocks_caller(const uint popID) {
      const uint cpuThreadID = gpu_getThread();
      const gpuStream_t stream = gpu_getStream();
      (GET_SUBPOINTER(gpuMemoryManager, Realf, host_returnRealf, cpuThreadID))[0] = 0; // host_rhoLossAdjust
      // populations[popID].vmesh->print();
      // Grow the vmesh and block container, if necessary. Try performing this on-device, if possible.
      resize_vbc_kernel_pre<<<1, 1, 0, stream>>> (
         populations[popID].dev_vmesh,
         populations[popID].dev_blockContainer,
         dev_list_with_replace_new,
         dev_list_delete,
         dev_list_to_replace,
         dev_list_with_replace_old,
         GET_SUBPOINTER(gpuMemoryManager, vmesh::LocalID, returnLID, cpuThreadID), // return values: nbefore, nafter, nblockstochange, resize success
         GET_SUBPOINTER(gpuMemoryManager, Realf, returnRealf, cpuThreadID) // mass loss, set to zero
         );
      CHK_ERR( gpuPeekAtLastError() );
      CHK_ERR( gpuMemcpyAsync(GET_SUBPOINTER(gpuMemoryManager, vmesh::LocalID, host_returnLID, cpuThreadID), GET_SUBPOINTER(gpuMemoryManager, vmesh::LocalID, returnLID, cpuThreadID), 4*sizeof(vmesh::LocalID), gpuMemcpyDeviceToHost, stream) );
      CHK_ERR( gpuStreamSynchronize(stream) );
      // Grow mesh if necessary and on-device resize did not work??
      const vmesh::LocalID nBlocksBeforeAdjust = (GET_SUBPOINTER(gpuMemoryManager, vmesh::LocalID, host_returnLID, cpuThreadID))[0];
      const vmesh::LocalID nBlocksAfterAdjust = (GET_SUBPOINTER(gpuMemoryManager, vmesh::LocalID, host_returnLID, cpuThreadID))[1];
      const vmesh::LocalID nBlocksToChange = (GET_SUBPOINTER(gpuMemoryManager, vmesh::LocalID, host_returnLID, cpuThreadID))[2];
      const vmesh::LocalID resizeDevSuccess = (GET_SUBPOINTER(gpuMemoryManager, vmesh::LocalID, host_returnLID, cpuThreadID))[3];
      if ( (nBlocksAfterAdjust > nBlocksBeforeAdjust) && (resizeDevSuccess == 0)) {
         //GPUTODO is _FACTOR enough instead of _PADDING?
         populations[popID].vmesh->setNewCapacity(nBlocksAfterAdjust*BLOCK_ALLOCATION_PADDING);
         populations[popID].vmesh->setNewSize(nBlocksAfterAdjust);
         populations[popID].blockContainer->setNewCapacity(nBlocksAfterAdjust*BLOCK_ALLOCATION_PADDING);
         populations[popID].blockContainer->setNewSize(nBlocksAfterAdjust);
         populations[popID].Upload();
      }

      if (nBlocksToChange==0) {
         return nBlocksAfterAdjust;
      }

      // Each GPU block / workunit could handle several Vlasiator velocity blocks at once.
      // However, thread syncs inside the kernel prevent this.
      // const uint vlasiBlocksPerWorkUnit = WARPSPERBLOCK * GPUTHREADS / WID3;
      const uint vlasiBlocksPerWorkUnit = 1;
      // ceil int division
      const uint launchBlocks = 1 + ((nBlocksToChange - 1) / vlasiBlocksPerWorkUnit);

      // Third argument specifies the number of bytes in *shared memory* that is
      // dynamically allocated per block for this call in addition to the statically allocated memory.
      //CHK_ERR( gpuStreamSynchronize(stream) );
      update_velocity_blocks_kernel<<<launchBlocks, vlasiBlocksPerWorkUnit * WID3, 0, stream>>> (
         populations[popID].dev_vmesh,
         populations[popID].dev_blockContainer,
         dev_list_with_replace_new,
         dev_list_delete,
         dev_list_to_replace,
         dev_list_with_replace_old,
         nBlocksBeforeAdjust,
         nBlocksToChange,
         nBlocksAfterAdjust,
         GET_SUBPOINTER(gpuMemoryManager, Realf, returnRealf, cpuThreadID) // mass loss
         );
      CHK_ERR( gpuPeekAtLastError() );
      CHK_ERR( gpuMemcpyAsync(GET_SUBPOINTER(gpuMemoryManager, Realf, host_returnRealf, cpuThreadID), GET_SUBPOINTER(gpuMemoryManager, Realf, returnRealf, cpuThreadID), sizeof(Realf), gpuMemcpyDeviceToHost, stream) );

      // Shrink the vmesh and block container, if necessary
      if (nBlocksAfterAdjust < nBlocksBeforeAdjust) {
         // Should not re-allocate on shrinking, so do on-device
         resize_vbc_kernel_post<<<1, 1, 0, stream>>> (
            populations[popID].dev_vmesh,
            populations[popID].dev_blockContainer,
            nBlocksAfterAdjust
            );
         CHK_ERR( gpuPeekAtLastError() );
      }
      // Update vmesh cached size
      populations[popID].vmesh->setNewCachedSize(nBlocksAfterAdjust);
      populations[popID].blockContainer->setNewCachedSize(nBlocksAfterAdjust);

      CHK_ERR( gpuStreamSynchronize(stream) );
      this->populations[popID].RHOLOSSADJUST += (GET_SUBPOINTER(gpuMemoryManager, Realf, host_returnRealf, cpuThreadID))[0];

      // DEBUG output after kernel
      #ifdef DEBUG_SPATIAL_CELL
      const vmesh::LocalID nAll = populations[popID].vmesh->size();
      if (nAll!=nBlocksAfterAdjust) {
         populations[popID].vmesh->gpu_prefetchHost();
         CHK_ERR( gpuStreamSynchronize(stream) );
         printf("after kernel, size is %d should be %d\n",nAll,nBlocksAfterAdjust);
         for (vmesh::LocalID m=0; m<nAll; ++m) {
            const vmesh::GlobalID GIDs = populations[popID].vmesh->getGlobalID(m);
            const vmesh::LocalID LIDs = populations[popID].vmesh->getLocalID(GIDs);
            printf("LID %d GID-solved %d LID-solved %d\n",m,GIDs,LIDs);
         }
         populations[popID].vmesh->gpu_prefetchDevice();
      }
      #endif
      return nBlocksAfterAdjust;
   }

   void SpatialCell::adjustSingleCellVelocityBlocks(const uint popID, bool doDeleteEmpty) {
      debug_population_check(popID);
      update_velocity_block_content_lists(popID);
      adjust_velocity_blocks(popID,doDeleteEmpty);
   }

   /** Update the two lists containing blocks with content, and blocks without content.
       This updates lists for a single cell, unlike the batch operations.
    * @see adjustVelocityBlocks */
   void SpatialCell::update_velocity_block_content_lists(const uint popID) {
      debug_population_check(popID);
      const gpuStream_t stream = gpu_getStream();

      applyReservation(popID);

      velocity_block_with_content_list_size = 0;
      velocity_block_with_content_map->clear<false>(Hashinator::targets::device,stream,std::pow(2,vbwcl_sizePower));
      velocity_block_with_no_content_map->clear<false>(Hashinator::targets::device,stream,std::pow(2,vbwncl_sizePower));
      CHK_ERR( gpuStreamSynchronize(stream) );

      const uint nBlocks = populations[popID].vmesh->size();
      if (nBlocks==0) {
         return;
      }
      CHK_ERR( gpuStreamSynchronize(stream) );

      const Real velocity_block_min_value = getVelocityBlockMinValue(popID);
      // Each GPU block / workunit can handle several Vlasiator velocity blocks at once. (TODO FIX)
      //const uint vlasiBlocksPerWorkUnit = WARPSPERBLOCK * GPUTHREADS / WID3;
      const uint vlasiBlocksPerWorkUnit = 1;
      // ceil int division
      const uint launchBlocks = 1 + ((nBlocks - 1) / vlasiBlocksPerWorkUnit);

      // Third argument specifies the number of bytes in *shared memory* that is
      // dynamically allocated per block for this call in addition to the statically allocated memory.
      //update_velocity_block_content_lists_kernel<<<launchBlocks, WID3, WID3*sizeof(int), stream>>> (
      update_velocity_block_content_lists_kernel<<<launchBlocks, (vlasiBlocksPerWorkUnit * WID3), 0, stream>>> (
         populations[popID].dev_vmesh,
         populations[popID].dev_blockContainer,
         dev_velocity_block_with_content_map,
         dev_velocity_block_with_no_content_map,
         velocity_block_min_value
         );
      CHK_ERR( gpuPeekAtLastError() );

      // Now extract values from the map
      velocity_block_with_content_map->extractAllKeysLoop(*dev_velocity_block_with_content_list,stream);
      split::SplitInfo info;
      velocity_block_with_content_list->copyMetadata(&info, stream);
      CHK_ERR( gpuStreamSynchronize(stream) );
      velocity_block_with_content_list_size = info.size;
   }

   void SpatialCell::prefetchDevice() {
      for (size_t p=0; p<populations.size(); ++p) {
         populations[p].vmesh->gpu_prefetchDevice();
         populations[p].blockContainer->gpu_prefetchDevice();
      }
   }
   void SpatialCell::prefetchHost() {
      for (size_t p=0; p<populations.size(); ++p) {
         populations[p].vmesh->gpu_prefetchHost();
         populations[p].blockContainer->gpu_prefetchHost();
      }
   }

   /** Get maximum translation timestep for the given species.
    * @param popID ID of the particle species.
    * @return Maximum timestep calculated by the Vlasov translation.*/
   const Real& SpatialCell::get_max_r_dt(const uint popID) const {
      debug_population_check(popID);
      return populations[popID].max_dt[species::MAXRDT];
   }

   /** Get maximum acceleration timestep for the given species.
    * @param popID ID of the particle species.
    * @return Maximum timestep calculated by Vlasov acceleration.*/
   const Real& SpatialCell::get_max_v_dt(const uint popID) const {
      debug_population_check(popID);
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
      if (this->mpiTransferEnabled && (SpatialCell::mpiTransferAtSysBoundaries==false ||
                                       this->sysBoundaryLayer ==1 || this->sysBoundaryLayer ==2 )) {

         //add data to send/recv to displacement and block length lists
         if ((SpatialCell::mpi_transfer_type & Transfer::VEL_BLOCK_LIST_STAGE1) != 0) {
            // first copy values in case this is the send operation
            populations[activePopID].N_blocks = populations[activePopID].vmesh->size();

            // send velocity block list size
            displacements.push_back((uint8_t*) &(populations[activePopID].N_blocks) - (uint8_t*) this);
            block_lengths.push_back(sizeof(vmesh::LocalID));
         }

         if ((SpatialCell::mpi_transfer_type & Transfer::VEL_BLOCK_LIST_STAGE2) != 0) {
            // STAGE1 should have been done, otherwise we have problems...
            if (receiving) {
               // Set population size based on mpi_number_of_blocks transferred earlier.
               // Does not need to be cleared. Vmesh map and VBC will be prepared in prepare_to_receive_blocks.
               this->dev_resize_vmesh(activePopID,populations[activePopID].N_blocks);
            } else {
               // Ensure N_blocks is still correct
               populations[activePopID].N_blocks = populations[activePopID].vmesh->size();
            }

            // send velocity block list
            if (populations[activePopID].N_blocks > 0) {
               displacements.push_back((uint8_t*) populations[activePopID].vmesh->getGrid()->data() - (uint8_t*) this);
               block_lengths.push_back(sizeof(vmesh::GlobalID) * populations[activePopID].N_blocks);
            } else {
               displacements.push_back(0);
               block_lengths.push_back(0);
            }
         }

         if ((SpatialCell::mpi_transfer_type & Transfer::VEL_BLOCK_WITH_CONTENT_STAGE1) !=0) {
            //Communicate size of list so that buffers can be allocated on receiving side
            // if (!receiving) { // Already done during block evaluation
            //    this->velocity_block_with_content_list_size = velocity_block_with_content_list->size();
            // }
            displacements.push_back((uint8_t*) &(this->velocity_block_with_content_list_size) - (uint8_t*) this);
            block_lengths.push_back(sizeof(vmesh::LocalID));
         }
         if ((SpatialCell::mpi_transfer_type & Transfer::VEL_BLOCK_WITH_CONTENT_STAGE2) !=0) {
            const gpuStream_t stream = gpu_getStream();
            if (receiving) {
               if (this->velocity_block_with_content_list_capacity < this->velocity_block_with_content_list_size) {
                  this->velocity_block_with_content_list->reserve(this->velocity_block_with_content_list_size);
                  this->velocity_block_with_content_list_capacity = this->velocity_block_with_content_list->capacity();
                  this->velocity_block_with_content_list->resize(this->velocity_block_with_content_list_size,true);
                  this->dev_velocity_block_with_content_list = this->velocity_block_with_content_list->upload<true>();
                  //this->velocity_block_with_content_list->optimizeGPU(stream); // included in upload<true>()
               } else {
                  this->velocity_block_with_content_list->resize(this->velocity_block_with_content_list_size,true);
                  this->dev_velocity_block_with_content_list = this->velocity_block_with_content_list->upload<false>();
               }
             }
            //velocity_block_with_content_list_size should first be updated, before this can be done (STAGE1)
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

         // Refinement parameters
         if ((SpatialCell::mpi_transfer_type & Transfer::REFINEMENT_PARAMETERS)){
            displacements.push_back(reinterpret_cast<uint8_t*>(this->parameters.data() + CellParams::AMR_ALPHA1) - reinterpret_cast<uint8_t*>(this));
            block_lengths.push_back(sizeof(Real) * (CellParams::AMR_ALPHA2 - CellParams::AMR_ALPHA1 + 1)); // This is just 2, but let's be explicit
         }

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
    * At this stage we have received a new block list (over MPI or from
    * an initializatiom function), but the rest of the cell structures
    * have not been adapted to this new list. Here we re-initialize
    * the cell with empty blocks based on the new list.*/
   void SpatialCell::prepare_to_receive_blocks(const uint popID) {
      setNewSizeClear(popID);
      // As the globalToLocalMap is empty, instead of calling
      // vmesh->setGrid() we can update both that and the block
      // parameters with a single kernel launch.

      //populations[popID].vmesh->setGrid(); // Based on localToGlobalMap
      const gpuStream_t stream = gpu_getStream();
      const uint newSize = populations[popID].N_blocks;
      // Set velocity block parameters:
      if (newSize>0) {
         // ceil int division
         #ifdef USE_WARPACCESSORS
         const uint launchBlocks = 1 + ((newSize - 1) / (WARPSPERBLOCK));
         #else
         const uint launchBlocks = 1 + ((newSize - 1) / (WARPSPERBLOCK*GPUTHREADS));
         #endif
         update_vmesh_and_blockparameters_kernel<<<launchBlocks, (WARPSPERBLOCK*GPUTHREADS), 0, stream>>> (
            populations[popID].dev_vmesh,
            populations[popID].dev_blockContainer,
            newSize
            );
         CHK_ERR( gpuPeekAtLastError() );
         CHK_ERR( gpuStreamSynchronize(stream) );
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
      debug_population_check(popID);
      populations[popID].max_dt[species::MAXRDT] = value;
   }

   /** Set maximum acceleration timestep for a particle species.
    * This function is called during Vlasov acceleration.
    * @param popID ID of the particle species.
    * @param value New maximum timestep.*/
   void SpatialCell::set_max_v_dt(const uint popID,const Real& value) {
      debug_population_check(popID);
      populations[popID].max_dt[species::MAXVDT] = value;
   }

   /**  Purges extra capacity from block vectors. It sets size to
    * num_blocks * block_allocation_factor (if capacity greater than this),
    * and also forces capacity to this new smaller value.
    * @return True on success.*/
   bool SpatialCell::shrink_to_fit() {
      bool success = true;
      return true; // on AMD, shrink_to_fit appears broken.
      size_t largestAmount = 0;
      for (size_t popID=0; popID<populations.size(); ++popID) {
         const vmesh::LocalID amount
            = 2 + populations[popID].blockContainer->size()
            * populations[popID].blockContainer->getBlockAllocationFactor();
         largestAmount = std::max(largestAmount,(size_t)populations[popID].blockContainer->size());
         // Allow capacity to be a bit larger than needed by number of blocks, shrink otherwise
         if (populations[popID].blockContainer->capacity() > amount ) {
            if (populations[popID].blockContainer->setNewCapacityShrink(amount) == false) {
               success = false;
            }
            populations[popID].Upload();
         }
      }
      largestvmesh = largestAmount;
      return success;
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
