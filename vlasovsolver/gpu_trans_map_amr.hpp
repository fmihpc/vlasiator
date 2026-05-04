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
#ifndef GPU_TRANS_MAP_AMR_H
#define GPU_TRANS_MAP_AMR_H

#include <vector>
#include "vec.h"
#include "../common.h"
#include "../spatial_cells/spatial_cell_wrapper.hpp"

// Buffers, vectors and set for use in translation. Cannot be declared in arch/gpu_base.hpp, because of
// compilation errors.
extern split::SplitVector<vmesh::GlobalID> *unionOfBlocks, *dev_unionOfBlocks;
extern Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID> *unionOfBlocksSet, *dev_unionOfBlocksSet;

bool trans_map_1d_amr(const dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                  const std::vector<CellID>& localPropagatedCells,
                  const std::vector<CellID>& remoteTargetCells,
                  std::vector<uint>& nPencils,
                  const uint dimension,
                  const Realf dt,
                  const uint popID);

void update_remote_mapping_contribution_amr(dccrg::Dccrg<spatial_cell::SpatialCell,
                                            dccrg::Cartesian_Geometry>& mpiGrid,
                                            const uint dimension,
                                            int direction,
                                            const uint popID);
#endif
