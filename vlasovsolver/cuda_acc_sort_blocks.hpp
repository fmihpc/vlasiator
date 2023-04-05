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


#ifndef CUDA_SORT_BLOCKS_FOR_ACC_H
#define CUDA_SORT_BLOCKS_FOR_ACC_H

#include "../common.h"
#include "../spatial_cell.hpp"
//inherited include of cuda_context.cuh

//#include "include/splitvector/splitvec.h"

void sortBlocklistByDimension( //const spatial_cell::SpatialCell* spatial_cell,
                               vmesh::VelocityMesh* vmesh,
                               const vmesh::LocalID nBlocks,
                               const uint dimension,
                               vmesh::GlobalID *blocksID_mapped,
                               vmesh::GlobalID *blocksID_mapped_sorted,
                               vmesh::GlobalID *blocksGID,
                               vmesh::LocalID *blocksLID_unsorted,
                               vmesh::LocalID *blocksLID,
                               vmesh::LocalID *dev_columnNBlocks,
                               ColumnOffsets* columnData,
                               // split::SplitVector<uint> & columnBlockOffsets,
                               // split::SplitVector<uint> & columnNumBlocks,
                               // split::SplitVector<uint> & setColumnOffsets,
                               // split::SplitVector<uint> & setNumColumns
                               const uint cuda_async_queue_id,
                               cudaStream_t stream
   );

#endif
