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


#ifndef GPU_SORT_BLOCKS_FOR_ACC_H
#define GPU_SORT_BLOCKS_FOR_ACC_H

#include "../common.h"
#include "../spatial_cells/spatial_cell_wrapper.hpp"

void sortBlocklistByDimension( vmesh::VelocityMesh* vmesh,
                               const vmesh::LocalID nBlocks,
                               const uint dimension,
                               vmesh::GlobalID *blocksID_mapped,
                               vmesh::GlobalID *blocksID_mapped_sorted,
                               vmesh::GlobalID *blocksGID,
                               vmesh::LocalID *blocksLID_unsorted,
                               vmesh::LocalID *blocksLID,
                               vmesh::LocalID *gpu_columnNBlocks,
                               ColumnOffsets* columnData,
                               const uint gpu_async_queue_id,
                               gpuStream_t stream
   );

extern void *gpu_RadixSortTemp[]; // Declared in gpu_acc_map.cpp
extern size_t gpu_acc_RadixSortTempSize[];

#endif
