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

#ifndef VLASIATOR_BLOCK_ADJUST_GPU_HPP
#define VLASIATOR_BLOCK_ADJUST_GPU_HPP

#include "spatial_cell_gpu.hpp"

#include "../definitions.h"
#include <dccrg.hpp>
#include <dccrg_cartesian_geometry.hpp>

#ifdef DEBUG_VLASIATOR
// Re-use spatial cell debug flag also for batch operations
   #ifndef DEBUG_SPATIAL_CELL
   #define DEBUG_SPATIAL_CELL
   #endif
#endif

namespace spatial_cell {
   // Following functions act on all requested cells

   void update_velocity_block_content_lists(
      dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const vector<CellID>& cells,
      const uint popID=0);

   void adjust_velocity_blocks_in_cells(
      dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const vector<CellID>& cells,
      const uint popID=0);

   void clear_maps_caller(
      const uint nCells,
      const size_t largestSizePower,
      gpuStream_t stream=0,
      const size_t offset=0);

   void batch_adjust_blocks_caller(
      dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const vector<CellID>& cellsToAdjust,
      const uint cellOffset,
      uint &largestBlocksToChange,
      uint &largestBlocksBeforeOrAfter,
      const uint popID=0
      );

   // Non-templated caller functions due to dual use from both block adjustment and acceleration
   void extract_to_replace_caller(
      Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>** input_maps,
      split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>> **output_vecs,
      vmesh::LocalID* output_sizes,
      vmesh::VelocityMesh** rule_meshes,
      Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>** rule_maps,
      split::SplitVector<vmesh::GlobalID>** rule_vectors,
      const uint nCells,
      gpuStream_t stream
      );

   void extract_to_delete_or_move_caller(
      Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>** input_maps,
      split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>> **output_vecs,
      vmesh::LocalID* output_sizes,
      vmesh::VelocityMesh** rule_meshes,
      Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>** rule_maps,
      split::SplitVector<vmesh::GlobalID>** rule_vectors,
      const uint nCells,
      gpuStream_t stream
      );

   void extract_to_add_caller(
      Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>** input_maps,
      split::SplitVector<vmesh::GlobalID> **output_vecs,
      vmesh::LocalID* output_sizes,
      vmesh::VelocityMesh** rule_meshes,
      Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>** rule_maps,
      split::SplitVector<vmesh::GlobalID>** rule_vectors,
      const uint nCells,
      gpuStream_t stream
      );

} // namespaces

#endif
