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

} // namespaces

extern vmesh::VelocityMesh** host_vmeshes, **dev_vmeshes;
extern vmesh::VelocityBlockContainer** host_VBCs, **dev_VBCs;
extern Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>** host_allMaps, **dev_allMaps;
extern split::SplitVector<vmesh::GlobalID> ** host_vbwcl_vec, **dev_vbwcl_vec;
extern split::SplitVector<vmesh::GlobalID> ** host_lists_with_replace_new, **dev_lists_with_replace_new;
extern split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>> **host_lists_delete, **dev_lists_delete;
extern split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>> **host_lists_to_replace, **dev_lists_to_replace;
extern split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>> **host_lists_with_replace_old, **dev_lists_with_replace_old;
extern split::SplitVector<vmesh::GlobalID> ** host_vbwcl_neigh, **dev_vbwcl_neigh;
extern vmesh::LocalID* host_contentSizes, *dev_contentSizes;
extern Real* host_minValues, *dev_minValues;
extern Real* host_massLoss, *dev_massLoss;
extern Real* host_mass, *dev_mass;

#endif
