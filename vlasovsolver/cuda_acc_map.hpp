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
#ifndef CUDA_ACC_MAP_H
#define CUDA_ACC_MAP_H

#include "../spatial_cell.hpp"
#include "vec.h"
#include "../common.h"

#include "device_launch_parameters.h"
#include "cuda.h"
#include "cuda_runtime.h"
#include "cuda_acc_sort_blocks.hpp"
#include "../cuda_context.cuh"

//using namespace spatial_cell;

struct Column {
   int valuesOffset;                              // Source data values
   size_t targetBlockOffsets[MAX_BLOCKS_PER_DIM]; // Target data array offsets
   int nblocks;                                   // Number of blocks in this column
   int minBlockK,maxBlockK;                       // Column parallel coordinate limits
   int kBegin;                                    // Actual un-sheared starting block index
   int i,j;                                       // Blocks' perpendicular coordinates
};


bool cuda_acc_map_1d(spatial_cell::SpatialCell* spatial_cell,
                     const uint popID,
                     Realv intersection,
                     Realv intersection_di,
                     Realv intersection_dj,
                     Realv intersection_dk,
                     const uint dimension,
                     cudaStream_t stream
   );

extern void cuda_acc_allocate (
   uint maxBlockCount
   );
extern void cuda_acc_allocate_memory (
   uint cpuThreadID,
   uint maxBlockCount
   );
extern void cuda_acc_deallocate_memory (
   uint cpuThreadID
   );

// Device data variables, to be allocated in good time. Made into an array so that each thread has their own pointer.
extern Vec *dev_blockDataOrdered[];
//extern Column *dev_columns[];
extern uint *dev_cell_indices_to_id[];
extern uint *dev_columnNumBlocks[];
extern uint *dev_columnBlockOffsets[];
extern vmesh::GlobalID *dev_BlocksID_mapped[];
extern vmesh::GlobalID *dev_BlocksID_mapped_sorted[];
extern vmesh::GlobalID *dev_LIDlist_unsorted[];

// Unified (managed) memory variables
extern Column *unif_columns[];
extern uint *unif_GIDlist[];
extern uint *unif_LIDlist[];
extern ColumnOffsets *unif_columndata[];

extern uint cuda_acc_allocatedSize;
extern uint cuda_acc_allocatedColumns;

#endif
