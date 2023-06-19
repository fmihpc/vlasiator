/*
 * This file is part of Vlasiator.
 * Copyright 2010-2022 Finnish Meteorological Institute and University of Helsinki
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cuda_acc_map.hpp"
#include "vec.h"
#include "../definitions.h"
#include "../object_wrapper.h"
#include "../cuda_context.cuh"
#include "../spatial_cell_cuda.hpp"

#include "cpu_face_estimates.hpp"
#include "cpu_1d_pqm.hpp"
#include "cpu_1d_ppm.hpp"
#include "cpu_1d_plm.hpp"

#include "cuda_acc_sort_blocks.hpp"

// Extra profiling stream synchronizations?
#define SSYNC HANDLE_ERROR( cudaStreamSynchronize(stream) )
//#define SSYNC

#define i_pcolumnv_cuda(j, k, k_block, num_k_blocks) ( ((j) / ( VECL / WID)) * WID * ( num_k_blocks + 2) + (k) + ( k_block + 1 ) * WID )
#define i_pcolumnv_cuda_b(planeVectorIndex, k, k_block, num_k_blocks) ( planeVectorIndex * WID * ( num_k_blocks + 2) + (k) + ( k_block + 1 ) * WID )

using namespace std;
using namespace spatial_cell;

// Allocate pointers for per-thread memory regions
// #define MAXCPUTHREADS 64 now in cuda_context.hpp
Vec *dev_blockDataOrdered[MAXCPUTHREADS];
uint *dev_cell_indices_to_id[MAXCPUTHREADS];

uint *dev_block_indices_to_id[MAXCPUTHREADS];
ColumnOffsets *unif_columnOffsetData[MAXCPUTHREADS];
Column *dev_columns[MAXCPUTHREADS];

vmesh::LocalID *dev_GIDlist[MAXCPUTHREADS];
vmesh::LocalID *dev_LIDlist[MAXCPUTHREADS];
vmesh::GlobalID *dev_BlocksID_mapped[MAXCPUTHREADS];
vmesh::GlobalID *dev_BlocksID_mapped_sorted[MAXCPUTHREADS];
vmesh::GlobalID *dev_LIDlist_unsorted[MAXCPUTHREADS];
vmesh::LocalID *dev_columnNBlocks[MAXCPUTHREADS];
void *dev_RadixSortTemp[MAXCPUTHREADS];

// Memory allocation flags and values.
uint cuda_acc_allocatedSize = 0;
uint cuda_acc_allocatedColumns = 0;
uint cuda_acc_columnContainerSize = 0;
uint cuda_acc_foundColumnsCount = 0;
uint cuda_acc_RadixSortTempSize[MAXCPUTHREADS] = {0};

__host__ void cuda_acc_allocate (
   uint maxBlockCount
   ) {
   // Always prepare for at least 500 blocks
   const uint maxBlocksPerCell = maxBlockCount > 500 ? maxBlockCount : 500;
   // Check if we already have allocated enough memory?
   if ( (cuda_acc_allocatedSize > maxBlocksPerCell * BLOCK_ALLOCATION_FACTOR) &&
        (cuda_acc_allocatedColumns > cuda_acc_foundColumnsCount * BLOCK_ALLOCATION_FACTOR) &&
        (cuda_acc_foundColumnsCount > 500)) {
      return;
   }
   // Deallocate before allocating new memory
#ifdef _OPENMP
   const uint maxNThreads = omp_get_max_threads();
#else
   const uint maxNThreads = 1;
#endif
   for (uint i=0; i<maxNThreads; ++i) {
      if (cuda_acc_allocatedSize > 0) {
         cuda_acc_deallocate_memory(i);
      }
      cuda_acc_allocate_memory(i, maxBlocksPerCell);
   }
}

__host__ void cuda_acc_allocate_memory (
   uint cpuThreadID,
   uint maxBlockCount
   ) {
   // Mallocs should be in increments of 256 bytes. WID3 is at least 64, and len(Realf) is at least 4, so this is true.
   const uint blockAllocationCount = maxBlockCount * BLOCK_ALLOCATION_PADDING;
   // The worst case scenario for columns is with every block having content but no neighbours, creating up
   // to maxBlockCount columns with each needing three blocks (one value plus two for padding).
   // This value here is a reasonable estimate for most cases, though. Also checks against a counter.
   uint maxColumnsPerCell = 3 * std::pow(
      (*vmesh::getMeshWrapper()->velocityMeshes)[0].gridLength[0]
      * (*vmesh::getMeshWrapper()->velocityMeshes)[0].gridLength[1]
      * (*vmesh::getMeshWrapper()->velocityMeshes)[0].gridLength[2], 0.667);
   if (cuda_acc_foundColumnsCount * BLOCK_ALLOCATION_FACTOR > maxColumnsPerCell) {
      maxColumnsPerCell = cuda_acc_foundColumnsCount * BLOCK_ALLOCATION_PADDING;
   }
   cuda_acc_allocatedSize = blockAllocationCount;
   cuda_acc_allocatedColumns = maxColumnsPerCell;

   HANDLE_ERROR( cudaMalloc((void**)&dev_cell_indices_to_id[cpuThreadID], 3*sizeof(uint)) );
   HANDLE_ERROR( cudaMalloc((void**)&dev_block_indices_to_id[cpuThreadID], 3*sizeof(uint)) );
   HANDLE_ERROR( cudaMalloc((void**)&dev_blockDataOrdered[cpuThreadID], blockAllocationCount * (WID3 / VECL) * sizeof(Vec)) );
   HANDLE_ERROR( cudaMalloc((void**)&dev_BlocksID_mapped[cpuThreadID], blockAllocationCount*sizeof(vmesh::GlobalID)) );
   HANDLE_ERROR( cudaMalloc((void**)&dev_BlocksID_mapped_sorted[cpuThreadID], blockAllocationCount*sizeof(vmesh::GlobalID)) );
   HANDLE_ERROR( cudaMalloc((void**)&dev_LIDlist_unsorted[cpuThreadID], blockAllocationCount*sizeof(vmesh::GlobalID)) );
   HANDLE_ERROR( cudaMalloc((void**)&dev_LIDlist[cpuThreadID], blockAllocationCount*sizeof(vmesh::LocalID)) );
   HANDLE_ERROR( cudaMalloc((void**)&dev_GIDlist[cpuThreadID], blockAllocationCount*sizeof(vmesh::LocalID)) );
   HANDLE_ERROR( cudaMalloc((void**)&dev_columns[cpuThreadID], maxColumnsPerCell*sizeof(Column)) );

   // Unified memory; columndata contains several splitvectors.
   unif_columnOffsetData[cpuThreadID] = new ColumnOffsets(maxColumnsPerCell); // inherits managed

   // Potential ColumnSet block count container
   const uint c0 = (*vmesh::getMeshWrapper()->velocityMeshes)[0].gridLength[0];
   const uint c1 = (*vmesh::getMeshWrapper()->velocityMeshes)[0].gridLength[1];
   const uint c2 = (*vmesh::getMeshWrapper()->velocityMeshes)[0].gridLength[2];
   std::array<uint, 3> s = {c0,c1,c2};
   std::sort(s.begin(), s.end());
   cuda_acc_columnContainerSize = c2*c1;
   HANDLE_ERROR( cudaMalloc((void**)&dev_columnNBlocks[cpuThreadID], cuda_acc_columnContainerSize*sizeof(vmesh::LocalID)) );
}

__host__ void cuda_acc_deallocate_memory (
   uint cpuThreadID
   ) {
   HANDLE_ERROR( cudaFree(dev_cell_indices_to_id[cpuThreadID]) );
   HANDLE_ERROR( cudaFree(dev_block_indices_to_id[cpuThreadID]) );
   HANDLE_ERROR( cudaFree(dev_blockDataOrdered[cpuThreadID]) );
   HANDLE_ERROR( cudaFree(dev_BlocksID_mapped[cpuThreadID]) );
   HANDLE_ERROR( cudaFree(dev_BlocksID_mapped_sorted[cpuThreadID]) );
   HANDLE_ERROR( cudaFree(dev_LIDlist_unsorted[cpuThreadID]) );
   HANDLE_ERROR( cudaFree(dev_columnNBlocks[cpuThreadID]) );
   HANDLE_ERROR( cudaFree(dev_LIDlist[cpuThreadID]) );
   HANDLE_ERROR( cudaFree(dev_GIDlist[cpuThreadID]) );
   HANDLE_ERROR( cudaFree(dev_columns[cpuThreadID]) );

   delete unif_columnOffsetData[cpuThreadID];

   cuda_acc_allocatedSize = 0;
   cuda_acc_allocatedColumns = 0;
   cuda_acc_columnContainerSize = 0;
}

__host__ void inline swapBlockIndices(std::array<uint32_t,3> &blockIndices, const uint dimension){
   uint temp;
   // Switch block indices according to dimensions, the algorithm has
   // been written for integrating along z.
   switch (dimension){
   case 0:
      /*i and k coordinates have been swapped*/
      temp=blockIndices[2];
      blockIndices[2]=blockIndices[0];
      blockIndices[0]=temp;
      break;
   case 1:
      /*in values j and k coordinates have been swapped*/
      temp=blockIndices[2];
      blockIndices[2]=blockIndices[1];
      blockIndices[1]=temp;
      break;
   case 2:
      break;
   }
}

__device__ void inline swapBlockIndices(vmesh::LocalID &blockIndices0,vmesh::LocalID &blockIndices1,vmesh::LocalID &blockIndices2, const uint dimension){
   vmesh::LocalID temp;
   // Switch block indices according to dimensions, the algorithm has
   // been written for integrating along z.
   switch (dimension){
   case 0:
      /*i and k coordinates have been swapped*/
      temp=blockIndices2;
      blockIndices2=blockIndices0;
      blockIndices0=temp;
      break;
   case 1:
      /*in values j and k coordinates have been swapped*/
      temp=blockIndices2;
      blockIndices2=blockIndices1;
      blockIndices1=temp;
      break;
   case 2:
      break;
   }
}

__global__ void reorder_blocks_by_dimension_kernel(
   Realf *dev_blockData,
   Vec *dev_blockDataOrdered,
   uint *dev_cell_indices_to_id,
   uint totalColumns,
   vmesh::LocalID *dev_LIDlist,
   ColumnOffsets* columnData,
   // These are just cleared
   split::SplitVector<vmesh::GlobalID>* BlocksRequired,
   split::SplitVector<vmesh::GlobalID>* BlocksToRemove,
   split::SplitVector<vmesh::GlobalID>* BlocksToAdd,
   split::SplitVector<vmesh::GlobalID>* BlocksToMove
) {
   // Takes the contents of blockData, sorts it into blockDataOrdered,
   // performing transposes as necessary
   // Works column-per-column and adds the necessary one empty block at each end
   const int nThreads = blockDim.x; // should be equal to VECL
   const int ti = threadIdx.x;
   const int blocki = blockIdx.x;
   const int cudaBlocks = gridDim.x;
   if (nThreads != VECL) {
      if (ti==0) printf("Warning! VECL not matching thread count for CUDA kernel!\n");
   }
   // Loop over columns in steps of cudaBlocks. Each cudaBlock deals with one column.
   for (uint iColumn = blocki; iColumn < totalColumns; iColumn += cudaBlocks) {
      uint inputOffset = columnData->columnBlockOffsets[iColumn];
      uint outputOffset = (inputOffset + 2 * iColumn) * (WID3/VECL);
      uint columnLength = columnData->columnNumBlocks[iColumn];

      // Loop over column blocks
      for (uint b = 0; b < columnLength; b++) {
         // Slices
         for (uint k=0; k<WID; ++k) {
            // Each block slice can span multiple VECLs (equal to cudathreads per block)
            for (uint j = 0; j < WID; j += VECL/WID) {
               // full-block index
               int input = k*WID2 + j*VECL + ti;
               // directional indices
               int input_2 = input / WID2; // last (slowest) index
               int input_1 = (input - input_2 * WID2) / WID; // medium index
               int input_0 = input - input_2 * WID2 - input_1 * WID; // first (fastest) index
               // slice vector index
               int jk = j / (VECL/WID);
               int sourceindex = input_0 * dev_cell_indices_to_id[0]
                  + input_1 * dev_cell_indices_to_id[1]
                  + input_2 * dev_cell_indices_to_id[2];

               dev_blockDataOrdered[outputOffset + i_pcolumnv_cuda_b(jk, k, b, columnLength)][ti]
                  = dev_blockData[ dev_LIDlist[inputOffset + b] * WID3
                                   + sourceindex ];

            } // end loop k (layers per block)
         } // end loop b (blocks per column)
      } // end loop j (vecs per layer)

      // Set first and last blocks to zero
      for (uint k=0; k<WID; ++k) {
         for (uint j = 0; j < WID; j += VECL/WID){
               int jk = j / (VECL/WID);
               dev_blockDataOrdered[outputOffset + i_pcolumnv_cuda_b(jk, k, -1, columnLength)][ti] = 0.0;
               dev_blockDataOrdered[outputOffset + i_pcolumnv_cuda_b(jk, k, columnLength, columnLength)][ti] = 0.0;
         }
      }

   } // end loop iColumn
   // Clear vectors
   if (blockIdx.x == blockIdx.y == blockIdx.z == threadIdx.x == threadIdx.y == threadIdx.z == 0) {
      BlocksRequired->clear();
      BlocksToRemove->clear();
      BlocksToAdd->clear();
      BlocksToMove->clear();
   }

   // Note: this kernel does not memset dev_blockData to zero.
   // A separate memsetasync call is required for that.
}


__global__ void identify_block_offsets_kernel(
   const vmesh::VelocityMesh* vmesh,
   Column *columns,
   const int totalColumns,
   const uint blockDataSize,
   uint *dev_block_indices_to_id
   ) {
   const uint cudaBlocks = gridDim.x * gridDim.y * gridDim.z;
   const uint warpSize = blockDim.x * blockDim.y * blockDim.z;
   const uint blocki = blockIdx.z*gridDim.x*gridDim.y + blockIdx.y*gridDim.x + blockIdx.x;
   const uint ti = threadIdx.z*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;

   for (uint iColumn = blocki; iColumn < totalColumns; iColumn += cudaBlocks) {
      for (int blockK = columns[iColumn].minBlockK; blockK <= columns[iColumn].maxBlockK; blockK += warpSize) {
         if (blockK+ti <= columns[iColumn].maxBlockK) {
            const int targetBlock =
               columns[iColumn].i * dev_block_indices_to_id[0] +
               columns[iColumn].j * dev_block_indices_to_id[1] +
               (blockK+ti)       * dev_block_indices_to_id[2];
            const vmesh::LocalID tblockLID = vmesh->getLocalID(targetBlock);
            if (tblockLID >= blockDataSize) {
               printf("Error: block for index [%d,%d,%d] has invalid blockLID %d %d \n",columns[iColumn].i,columns[iColumn].j,blockK+ti,tblockLID,vmesh->invalidGlobalID());
            }
            columns[iColumn].targetBlockOffsets[blockK+ti] = tblockLID*WID3;
         }
      }
   }
}

// Serial kernel only to avoid page faults or prefetches
__global__ void count_columns_kernel (
   ColumnOffsets* dev_columnData,
   vmesh::LocalID* returnLID // dev_totalColumns, dev_valuesSizeRequired
   ) {
   // const int cudaBlocks = gridDim.x * gridDim.y * gridDim.z;
   // const int warpSize = blockDim.x * blockDim.y * blockDim.z;
   const int blocki = blockIdx.z*gridDim.x*gridDim.y + blockIdx.y*gridDim.x + blockIdx.x;
   const int ti = threadIdx.z*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;
   if ((blocki==0)&&(ti==0)) {
      for(uint setIndex=0; setIndex< dev_columnData->setColumnOffsets.size(); ++setIndex) {
         returnLID[0] += dev_columnData->setNumColumns[setIndex];
         for(uint columnIndex = dev_columnData->setColumnOffsets[setIndex]; columnIndex < dev_columnData->setColumnOffsets[setIndex] + dev_columnData->setNumColumns[setIndex] ; columnIndex ++){
            returnLID[1] += (dev_columnData->columnNumBlocks[columnIndex] + 2) * WID3 / VECL;
         }
      }
   }
}

// Serial kernel only to avoid page faults or prefetches
__global__ void offsets_into_columns_kernel(
   ColumnOffsets* dev_columnData,
   Column *dev_columns,
   const int valuesSizeRequired
   ) {
   // const int cudaBlocks = gridDim.x * gridDim.y * gridDim.z;
   // const int warpSize = blockDim.x * blockDim.y * blockDim.z;
   const int blocki = blockIdx.z*gridDim.x*gridDim.y + blockIdx.y*gridDim.x + blockIdx.x;
   const int ti = threadIdx.z*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;
   if ((blocki==0)&&(ti==0)) {
      uint valuesColumnOffset = 0;
      for( uint setIndex=0; setIndex< dev_columnData->setColumnOffsets.size(); ++setIndex) {
         for (uint columnIndex = dev_columnData->setColumnOffsets[setIndex]; columnIndex < dev_columnData->setColumnOffsets[setIndex] + dev_columnData->setNumColumns[setIndex] ; columnIndex ++){
            dev_columns[columnIndex].nblocks = dev_columnData->columnNumBlocks[columnIndex];
            dev_columns[columnIndex].valuesOffset = valuesColumnOffset;
            if (valuesColumnOffset >= valuesSizeRequired) {
               printf("(ERROR: Overflowing the values array (%d > %d) with column %d\n",valuesColumnOffset,valuesSizeRequired,columnIndex);
            }
            valuesColumnOffset += (dev_columnData->columnNumBlocks[columnIndex] + 2) * (WID3/VECL); // there are WID3/VECL elements of type Vec per block
         }
      }
   }
}

// Using columns, evaluate which blocks are target or source blocks
__global__ void evaluate_column_extents_kernel(
   const uint dimension,
   const vmesh::VelocityMesh* vmesh,
   ColumnOffsets* dev_columnData,
   Column *dev_columns,
   //split::SplitVector<cuda::std::pair<vmesh::GlobalID,vmesh::LocalID>>* BlocksRequired,
   split::SplitVector<vmesh::GlobalID> *BlocksRequired,
   split::SplitVector<vmesh::GlobalID> *BlocksToAdd,
   split::SplitVector<vmesh::GlobalID> *BlocksToRemove,
   vmesh::GlobalID *GIDlist,
   uint *dev_block_indices_to_id,
   Realv intersection,
   Realv intersection_di,
   Realv intersection_dj,
   Realv intersection_dk,
   uint bailout_velocity_space_wall_margin,
   uint max_v_length,
   Realv v_min,
   Realv dv,
   uint *wallspace_margin_bailout_flag
   ) {
   const uint cudaBlocks = gridDim.x * gridDim.y * gridDim.z;
   const uint warpSize = blockDim.x * blockDim.y * blockDim.z;
   const uint blocki = blockIdx.z*gridDim.x*gridDim.y + blockIdx.y*gridDim.x + blockIdx.x;
   const uint ti = threadIdx.z*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;

   __shared__ int isTargetBlock[MAX_BLOCKS_PER_DIM];
   __shared__ int isSourceBlock[MAX_BLOCKS_PER_DIM];
   for( uint setIndex=blocki; setIndex < dev_columnData->setColumnOffsets.size(); setIndex += cudaBlocks) {
      if (setIndex < dev_columnData->setColumnOffsets.size()) {

         // Clear flags used for this columnSet
         for(uint tti = 0; tti < MAX_BLOCKS_PER_DIM; tti += warpSize ) {
            uint index = tti + ti;
            if (index < MAX_BLOCKS_PER_DIM) {
               isTargetBlock[index] = 0;
               isSourceBlock[index] = 0;
            }
         }
         __syncthreads();

         /*need x,y coordinate of this column set of blocks, take it from first
           block in first column*/
         vmesh::LocalID setFirstBlockIndices0,setFirstBlockIndices1,setFirstBlockIndices2;
         uint8_t refLevel=0;
         vmesh->getIndices(GIDlist[dev_columnData->columnBlockOffsets[dev_columnData->setColumnOffsets[setIndex]]],
                           refLevel,
                           setFirstBlockIndices0, setFirstBlockIndices1, setFirstBlockIndices2);
         swapBlockIndices(setFirstBlockIndices0,setFirstBlockIndices1,setFirstBlockIndices2,dimension);
         /*compute the maximum starting point of the lagrangian (target) grid
           (base level) within the 4 corner cells in this
           block. Needed for computig maximum extent of target column*/

         Realv max_intersectionMin = intersection +
            (setFirstBlockIndices0 * WID + 0) * intersection_di +
            (setFirstBlockIndices1 * WID + 0) * intersection_dj;
         max_intersectionMin =  std::max(max_intersectionMin,
                                         intersection +
                                         (setFirstBlockIndices0 * WID + 0) * intersection_di +
                                         (setFirstBlockIndices1 * WID + WID - 1) * intersection_dj);
         max_intersectionMin =  std::max(max_intersectionMin,
                                         intersection +
                                         (setFirstBlockIndices0 * WID + WID - 1) * intersection_di +
                                         (setFirstBlockIndices1 * WID + 0) * intersection_dj);
         max_intersectionMin =  std::max(max_intersectionMin,
                                         intersection +
                                         (setFirstBlockIndices0 * WID + WID - 1) * intersection_di +
                                         (setFirstBlockIndices1 * WID + WID - 1) * intersection_dj);

         Realv min_intersectionMin = intersection +
            (setFirstBlockIndices0 * WID + 0) * intersection_di +
            (setFirstBlockIndices1 * WID + 0) * intersection_dj;
         min_intersectionMin =  std::min(min_intersectionMin,
                                         intersection +
                                         (setFirstBlockIndices0 * WID + 0) * intersection_di +
                                         (setFirstBlockIndices1 * WID + WID - 1) * intersection_dj);
         min_intersectionMin =  std::min(min_intersectionMin,
                                         intersection +
                                         (setFirstBlockIndices0 * WID + WID - 1) * intersection_di +
                                         (setFirstBlockIndices1 * WID + 0) * intersection_dj);
         min_intersectionMin =  std::min(min_intersectionMin,
                                         intersection +
                                         (setFirstBlockIndices0 * WID + WID - 1) * intersection_di +
                                         (setFirstBlockIndices1 * WID + WID - 1) * intersection_dj);

         //now, record which blocks are target blocks
         for (uint columnIndex = dev_columnData->setColumnOffsets[setIndex];
              columnIndex < dev_columnData->setColumnOffsets[setIndex] + dev_columnData->setNumColumns[setIndex] ;
              ++columnIndex) {
            // Not parallelizing this at this level; not going to be many columns within a set

            const vmesh::LocalID n_cblocks = dev_columnData->columnNumBlocks[columnIndex];
            vmesh::GlobalID* cblocks = GIDlist + dev_columnData->columnBlockOffsets[columnIndex]; //column blocks
            vmesh::LocalID firstBlockIndices0,firstBlockIndices1,firstBlockIndices2;
            vmesh::LocalID lastBlockIndices0,lastBlockIndices1,lastBlockIndices2;
            vmesh->getIndices(cblocks[0],
                              refLevel,
                              firstBlockIndices0, firstBlockIndices1, firstBlockIndices2);
            vmesh->getIndices(cblocks[n_cblocks -1],
                              refLevel,
                              lastBlockIndices0, lastBlockIndices1, lastBlockIndices2);
            swapBlockIndices(firstBlockIndices0,firstBlockIndices1,firstBlockIndices2, dimension);
            swapBlockIndices(lastBlockIndices0,lastBlockIndices1,lastBlockIndices2, dimension);

            /* firstBlockV is in z the minimum velocity value of the lower
             *  edge in source grid.
             * lastBlockV is in z the maximum velocity value of the upper
             *  edge in source grid. */
            Realv firstBlockMinV = (WID * firstBlockIndices2) * dv + v_min;
            Realv lastBlockMaxV = (WID * (lastBlockIndices2 + 1)) * dv + v_min;

            /* gk is now the k value in terms of cells in target
               grid. This distance between max_intersectionMin (so lagrangian
               plan, well max value here) and V of source grid, divided by
               intersection_dk to find out how many grid cells that is*/
            const int firstBlock_gk = (int)((firstBlockMinV - max_intersectionMin)/intersection_dk);
            const int lastBlock_gk = (int)((lastBlockMaxV - min_intersectionMin)/intersection_dk);

            int firstBlockIndexK = firstBlock_gk/WID;
            int lastBlockIndexK = lastBlock_gk/WID;

            // now enforce mesh limits for target column blocks (and check if we are
            // too close to the velocity space boundaries)
            firstBlockIndexK = (firstBlockIndexK >= 0)            ? firstBlockIndexK : 0;
            firstBlockIndexK = (firstBlockIndexK < max_v_length ) ? firstBlockIndexK : max_v_length - 1;
            lastBlockIndexK  = (lastBlockIndexK  >= 0)            ? lastBlockIndexK  : 0;
            lastBlockIndexK  = (lastBlockIndexK  < max_v_length ) ? lastBlockIndexK  : max_v_length - 1;
            if(firstBlockIndexK < bailout_velocity_space_wall_margin
               || firstBlockIndexK >= max_v_length - bailout_velocity_space_wall_margin
               || lastBlockIndexK < bailout_velocity_space_wall_margin
               || lastBlockIndexK >= max_v_length - bailout_velocity_space_wall_margin
               ) {
               // Pass bailout flag back to host
               if (ti==0) {
                  *wallspace_margin_bailout_flag = 1;
               }
            }

            //store source blocks
            for (uint blockK = firstBlockIndices2; blockK <= lastBlockIndices2; blockK +=warpSize){
               if ((blockK+ti) <= lastBlockIndices2) {
                  //const int old  = atomicAdd(&isSourceBlock[blockK+ti],1);
                  isSourceBlock[blockK+ti] = 1; // Does not need to be atomic, as long as it's no longer zero
               }
            }
            __syncthreads();

            //store target blocks
            for (uint blockK = (uint)firstBlockIndexK; blockK <= (uint)lastBlockIndexK; blockK+=warpSize){
               if ((blockK+ti) <= (uint)lastBlockIndexK) {
                  isTargetBlock[blockK+ti] = 1; // Does not need to be atomic, as long as it's no longer zero
                  //const int old  = atomicAdd(&isTargetBlock[blockK+ti],1);
               }
            }
            __syncthreads();

            if (ti==0) {
               // Set columns' transverse coordinates
               dev_columns[columnIndex].i = setFirstBlockIndices0;
               dev_columns[columnIndex].j = setFirstBlockIndices1;
               dev_columns[columnIndex].kBegin = firstBlockIndices2;

               //store also for each column firstBlockIndexK, and lastBlockIndexK
               dev_columns[columnIndex].minBlockK = firstBlockIndexK;
               dev_columns[columnIndex].maxBlockK = lastBlockIndexK;
            }
         } // end loop over columns in set
         __syncthreads();

         for (uint blockT = 0; blockT < MAX_BLOCKS_PER_DIM; blockT +=warpSize) {
            uint blockK = blockT + ti;
            if (blockK < MAX_BLOCKS_PER_DIM) {
               if(isTargetBlock[blockK]!=0)  {
                  const int targetBlock =
                     setFirstBlockIndices0 * dev_block_indices_to_id[0] +
                     setFirstBlockIndices1 * dev_block_indices_to_id[1] +
                     blockK                * dev_block_indices_to_id[2];
                  BlocksRequired->device_push_back(targetBlock);
                  //BlocksRequired->device_push_back(cuda::std::pair<vmesh::GlobalID,vmesh::LocalID>(targetBlock,targetBlock));
               }
               if(isTargetBlock[blockK]!=0 && isSourceBlock[blockK]==0 )  {
                  const int targetBlock =
                     setFirstBlockIndices0 * dev_block_indices_to_id[0] +
                     setFirstBlockIndices1 * dev_block_indices_to_id[1] +
                     blockK                * dev_block_indices_to_id[2];
                  BlocksToAdd->device_push_back(targetBlock);

               }
               if(isTargetBlock[blockK]==0 && isSourceBlock[blockK]!=0 )  {
                  const int targetBlock =
                     setFirstBlockIndices0 * dev_block_indices_to_id[0] +
                     setFirstBlockIndices1 * dev_block_indices_to_id[1] +
                     blockK                * dev_block_indices_to_id[2];
                  BlocksToRemove->device_push_back(targetBlock);
               }
            } // block within MAX_BLOCKS_PER_DIM
         } // loop over all potential blocks
      } // if valid setIndez
   } // for setIndexB
}

__global__ void acceleration_kernel(
  Realf *dev_blockData,
  Vec *dev_blockDataOrdered,
  uint *dev_cell_indices_to_id,
  Column *dev_columns,
  int totalColumns,
  Realv intersection,
  Realv intersection_di,
  Realv intersection_dj,
  Realv intersection_dk,
  Realv v_min,
  Realv i_dv,
  Realv dv,
  Realv minValue,
  const uint bdsw3
) {
   const uint cudaBlocks = gridDim.x * gridDim.y * gridDim.z;
   //const uint warpSize = blockDim.x * blockDim.y * blockDim.z;
   const uint blocki = blockIdx.z*gridDim.x*gridDim.y + blockIdx.y*gridDim.x + blockIdx.x;
   const uint index = threadIdx.z*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;

   for (uint column = blocki; column < totalColumns; column += cudaBlocks) {
      /* New threading with each warp/wavefront working on one vector */
      Realf v_r0 = ( (WID * dev_columns[column].kBegin) * dv + v_min);

      // i,j,k are relative to the order in which we copied data to the values array.
      // After this point in the k,j,i loops there should be no branches based on dimensions
      // Note that the i dimension is vectorized, and thus there are no loops over i
      // Iterate through the perpendicular directions of the column
      for (uint j = 0; j < WID; j += VECL/WID) {
         // If VECL=WID2 (WID=4, VECL=16, or WID=8, VECL=64, then j==0)
         // This loop is still needed for e.g. Warp=VECL=32, WID2=64 (then j==0 or 4)
         const vmesh::LocalID nblocks = dev_columns[column].nblocks;

         uint i_indices = index % WID;
         uint j_indices = j + index/WID;
         //int jk = j / (VECL/WID);

         int target_cell_index_common =
            i_indices * dev_cell_indices_to_id[0] +
            j_indices * dev_cell_indices_to_id[1];
         const Realf intersection_min =
            intersection +
            (dev_columns[column].i * WID + (Realv)i_indices) * intersection_di +
            (dev_columns[column].j * WID + (Realv)j_indices) * intersection_dj;

         const Realf gk_intersection_min =
            intersection +
            (dev_columns[column].i * WID + (Realv)( intersection_di > 0 ? 0 : WID-1 )) * intersection_di +
            (dev_columns[column].j * WID + (Realv)( intersection_dj > 0 ? j : j+VECL/WID-1 )) * intersection_dj;
         const Realf gk_intersection_max =
            intersection +
            (dev_columns[column].i * WID + (Realv)( intersection_di < 0 ? 0 : WID-1 )) * intersection_di +
            (dev_columns[column].j * WID + (Realv)( intersection_dj < 0 ? j : j+VECL/WID-1 )) * intersection_dj;

         // loop through all perpendicular slices in column and compute the mapping as integrals.
         for (uint k=0; k < WID * nblocks; ++k) {
            // Compute reconstructions
            // Checked on 21.01.2022: Realv a[length] goes on the register despite being an array. Explicitly declaring it
            // as __shared__ had no impact on performance.
#ifdef ACC_SEMILAG_PLM
            Realv a[2];
            compute_plm_coeff(dev_blockDataOrdered + dev_columns[column].valuesOffset + i_pcolumnv_cuda(j, 0, -1, nblocks), (k + WID), a, minValue, index);
#endif
#ifdef ACC_SEMILAG_PPM
            Realv a[3];
            compute_ppm_coeff(dev_blockDataOrdered + dev_columns[column].valuesOffset + i_pcolumnv_cuda(j, 0, -1, nblocks), h4, (k + WID), a, minValue, index);
#endif
#ifdef ACC_SEMILAG_PQM
            Realv a[5];
            compute_pqm_coeff(dev_blockDataOrdered + dev_columns[column].valuesOffset + i_pcolumnv_cuda(j, 0, -1, nblocks), h8, (k + WID), a, minValue, index);
#endif

            // set the initial value for the integrand at the boundary at v = 0
            // (in reduced cell units), this will be shifted to target_density_1, see below.
            Realf target_density_r = 0.0;

            const Realv v_r = v_r0  + (k+1)* dv;
            const Realv v_l = v_r0  + k* dv;
            const int lagrangian_gk_l = trunc((v_l-gk_intersection_max)/intersection_dk);
            const int lagrangian_gk_r = trunc((v_r-gk_intersection_min)/intersection_dk);

            //limits in lagrangian k for target column. Also take into
            //account limits of target column
            // Now all indexes in the warp should have the same gk loop extents
            const int minGk = max(lagrangian_gk_l, int(dev_columns[column].minBlockK * WID));
            const int maxGk = min(lagrangian_gk_r, int((dev_columns[column].maxBlockK + 1) * WID - 1));
            // Run along the column and perform the polynomial reconstruction
            for(int gk = minGk; gk <= maxGk; gk++) {
               const int blockK = gk/WID;
               const int gk_mod_WID = (gk - blockK * WID);

               //the block of the Lagrangian cell to which we map
               //const int target_block(target_block_index_common + blockK * block_indices_to_id[2]);
               // This already contains the value index via target_cell_index_commom
               const int tcell(target_cell_index_common + gk_mod_WID * dev_cell_indices_to_id[2]);
               //the velocity between which we will integrate to put mass
               //in the targe cell. If both v_r and v_l are in same cell
               //then v_1,v_2 should be between v_l and v_r.
               //v_1 and v_2 normalized to be between 0 and 1 in the cell.
               //For vector elements where gk is already larger than needed (lagrangian_gk_r), v_2=v_1=v_r and thus the value is zero.
               const Realf v_norm_r = (  min(  max( (gk + 1) * intersection_dk + intersection_min, v_l), v_r) - v_l) * i_dv;

               /*shift, old right is new left*/
               const Realf target_density_l = target_density_r;

               // compute right integrand
#ifdef ACC_SEMILAG_PLM
               target_density_r = v_norm_r * ( a[0] + v_norm_r * a[1] );
#endif
#ifdef ACC_SEMILAG_PPM
               target_density_r = v_norm_r * ( a[0] + v_norm_r * ( a[1] + v_norm_r * a[2] ) );
#endif
#ifdef ACC_SEMILAG_PQM
               target_density_r = v_norm_r * ( a[0] + v_norm_r * ( a[1] + v_norm_r * ( a[2] + v_norm_r * ( a[3] + v_norm_r * a[4] ) ) ) );
#endif

               //store value
               Realf tval = target_density_r - target_density_l;
               if (isfinite(tval) && tval>0) {
                  (&dev_blockData[dev_columns[column].targetBlockOffsets[blockK]])[tcell] += tval;
               }
            } // for loop over target k-indices of current source block
         } // for-loop over source blocks
      } //for loop over j index
   } // End loop over columns
} // end semilag acc kernel

/*
   Here we map from the current time step grid, to a target grid which
   is the lagrangian departure grid (so th grid at timestep +dt,
   tracked backwards by -dt)
*/
__host__ bool cuda_acc_map_1d(spatial_cell::SpatialCell* spatial_cell,
                              const uint popID,
                              Realv intersection,
                              Realv intersection_di,
                              Realv intersection_dj,
                              Realv intersection_dk,
                              const uint dimension,
                              cudaStream_t stream
   ) {
   SSYNC;
   // Thread id used for persistent device memory pointers
#ifdef _OPENMP
      const uint cpuThreadID = omp_get_thread_num();
#else
      const uint cpuThreadID = 0;
#endif

   phiprof::start("Get acc parameters");
   vmesh::VelocityMesh* vmesh = spatial_cell->get_velocity_mesh(popID);
   vmesh->dev_prefetchDevice();
   vmesh::VelocityBlockContainer* blockContainer = spatial_cell->get_velocity_blocks(popID);
   Realf *blockData = blockContainer->getData();

   const uint nBlocks = vmesh->size();
   const vmesh::LocalID D0 = vmesh->getGridLength()[0];
   const vmesh::LocalID D1 = vmesh->getGridLength()[1];
   //const vmesh::LocalID D2 = vmesh->getGridLength()[2];
   const Realv dv    = vmesh->getCellSize()[dimension];
   const Realv v_min = vmesh->getMeshMinLimits()[dimension];
   const Realv max_v_length  = vmesh->getGridLength()[dimension];
   auto minValue = spatial_cell->getVelocityBlockMinValue(popID);
   const Realv i_dv = 1.0/dv;

   // For use later
   Real host_returnReal[8];
   Real *dev_returnReal = returnReal[cpuThreadID];
   vmesh::LocalID host_returnLID[8];
   vmesh::LocalID *dev_returnLID = returnLID[cpuThreadID];
   phiprof::stop("Get acc parameters");

   //nothing to do if no blocks
   if(nBlocks == 0) {
      return true;
   }

   phiprof::start("stream Attach, prefetch");
   if (needAttachedStreams) {
      blockContainer->dev_attachToStream(stream);
      vmesh->dev_attachToStream(stream);
   }
   if (doPrefetches) {
      blockContainer->dev_prefetchDevice();
      vmesh->dev_prefetchDevice();
   }
   phiprof::stop("stream Attach, prefetch");

   phiprof::start("create columnOffsets in unified memory, attach, prefetch");
   // lists in unified memory
   ColumnOffsets *columnData = unif_columnOffsetData[cpuThreadID];
   // Verify unified memory stream attach
   if (needAttachedStreams) {
      columnData->dev_attachToStream(stream);
   }
   if (doPrefetches) {
      columnData->columnBlockOffsets.optimizeGPU(stream);
      columnData->columnNumBlocks.optimizeGPU(stream);
      columnData->setColumnOffsets.optimizeGPU(stream);
      columnData->setNumColumns.optimizeGPU(stream);
      HANDLE_ERROR( cudaMemPrefetchAsync(columnData,sizeof(ColumnOffsets),cuda_getDevice(),stream) );
   }
   SSYNC;
   phiprof::stop("create columnOffsets in unified memory, attach, prefetch");
   // Some kernels in here require this to be equal to VECL.
   // Future improvements would be to allow setting it directly to WID3.
   // Other kernels (not handling block data) can use CUDATHREADS which
   // is equal to NVIDIA: 32 or AMD: 64.
   phiprof::start("Bookkeeping");
   int cudathreadsVECL = VECL;

   /*< used when computing id of target block, 0 for compiler */
   uint block_indices_to_id[3] = {0, 0, 0};
   uint cell_indices_to_id[3] = {0, 0, 0};

   Realv is_temp;
   switch (dimension) {
    case 0: /* i and k coordinates have been swapped*/
      /*swap intersection i and k coordinates*/
      is_temp=intersection_di;
      intersection_di=intersection_dk;
      intersection_dk=is_temp;

      /*set values in array that is used to convert block indices to id using a dot product*/
      block_indices_to_id[0] = D0*D1;
      block_indices_to_id[1] = D0;
      block_indices_to_id[2] = 1;

      /*set values in array that is used to convert block indices to id using a dot product*/
      cell_indices_to_id[0]=WID2;
      cell_indices_to_id[1]=WID;
      cell_indices_to_id[2]=1;
      break;
    case 1: /* j and k coordinates have been swapped*/
      /*swap intersection j and k coordinates*/
      is_temp=intersection_dj;
      intersection_dj=intersection_dk;
      intersection_dk=is_temp;

      /*set values in array that is used to convert block indices to id using a dot product*/
      block_indices_to_id[0]=1;
      block_indices_to_id[1] = D0*D1;
      block_indices_to_id[2] = D0;

      /*set values in array that is used to convert block indices to id using a dot product*/
      cell_indices_to_id[0]=1;
      cell_indices_to_id[1]=WID2;
      cell_indices_to_id[2]=WID;
      break;
    case 2:
      /*set values in array that is used to convert block indices to id using a dot product*/
      block_indices_to_id[0]=1;
      block_indices_to_id[1] = D0;
      block_indices_to_id[2] = D0*D1;

      // set values in array that is used to convert block indices to id using a dot product.
      cell_indices_to_id[0]=1;
      cell_indices_to_id[1]=WID;
      cell_indices_to_id[2]=WID2;
      break;
   }
   // Copy indexing information to device (async)
   HANDLE_ERROR( cudaMemcpyAsync(dev_cell_indices_to_id[cpuThreadID], cell_indices_to_id, 3*sizeof(uint), cudaMemcpyHostToDevice, stream) );
   HANDLE_ERROR( cudaMemcpyAsync(dev_block_indices_to_id[cpuThreadID], block_indices_to_id, 3*sizeof(uint), cudaMemcpyHostToDevice, stream) );

   // or device memory
   vmesh::GlobalID *GIDlist = dev_GIDlist[cpuThreadID];
   vmesh::LocalID *LIDlist = dev_LIDlist[cpuThreadID];
   vmesh::GlobalID *BlocksID_mapped = dev_BlocksID_mapped[cpuThreadID];
   vmesh::GlobalID *BlocksID_mapped_sorted = dev_BlocksID_mapped[cpuThreadID];
   vmesh::LocalID *LIDlist_unsorted = dev_LIDlist_unsorted[cpuThreadID];
   vmesh::LocalID *columnNBlocks = dev_columnNBlocks[cpuThreadID];

   // Call function for sorting block list and building columns from it.
   // Can probably be further optimized.
   phiprof::stop("Bookkeeping");
   phiprof::start("sortBlockList");
   cudaStream_t priorityStream = cuda_getPriorityStream();
   HANDLE_ERROR( cudaMemsetAsync(columnNBlocks, 0, cuda_acc_columnContainerSize*sizeof(vmesh::LocalID), priorityStream) );
   HANDLE_ERROR( cudaStreamSynchronize(stream) ); // Yes needed because we use priority stream for block list sorting
   sortBlocklistByDimension(vmesh,
                            nBlocks,
                            dimension,
                            BlocksID_mapped,
                            BlocksID_mapped_sorted,
                            GIDlist,
                            LIDlist_unsorted,
                            LIDlist,
                            columnNBlocks,
                            columnData,
                            cpuThreadID,
                            priorityStream
      );
   HANDLE_ERROR( cudaStreamSynchronize(priorityStream) ); // Yes needed to get column data back to regular stream
   phiprof::stop("sortBlockList");

   // Calculate total sum of columns and total values size
   phiprof::start("count columns");
   HANDLE_ERROR( cudaMemsetAsync(dev_returnLID, 0, 2*sizeof(vmesh::LocalID), stream) );

   // this needs to be serial, but is fast.
   SSYNC;
   count_columns_kernel<<<1, 1, 0, stream>>> (
      columnData,
      dev_returnLID //dev_totalColumns,dev_valuesSizeRequired
      );
   SSYNC;
   HANDLE_ERROR( cudaMemcpyAsync(host_returnLID, dev_returnLID, 2*sizeof(vmesh::LocalID), cudaMemcpyDeviceToHost, stream) );
   HANDLE_ERROR( cudaStreamSynchronize(stream) );
   const vmesh::LocalID host_totalColumns = host_returnLID[0];
   const vmesh::LocalID host_valuesSizeRequired = host_returnLID[1];
   // Update tracker of maximum encountered column count
   if (cuda_acc_foundColumnsCount < host_totalColumns) {
      cuda_acc_foundColumnsCount = host_totalColumns;
   }
   phiprof::stop("count columns");

   phiprof::start("Host columns");
   // Create array of column objects
   Column host_columns[host_totalColumns];
   // and copy it into device memory
   Column *columns = dev_columns[cpuThreadID];
   HANDLE_ERROR( cudaMemcpyAsync(columns, &host_columns, host_totalColumns*sizeof(Column), cudaMemcpyHostToDevice, stream) );
   SSYNC;
   phiprof::stop("Host columns");

   phiprof::start("Store offsets into columns");
   // this needs to be serial, but is fast.
   offsets_into_columns_kernel<<<1, 1, 0, stream>>> (
      columnData,
      columns,
      host_valuesSizeRequired
      );
   SSYNC;
   phiprof::stop("Store offsets into columns");

   phiprof::start("Reorder blocks by dimension");
   uint cudablocks = host_totalColumns > CUDABLOCKS ? CUDABLOCKS : host_totalColumns;
   // Launch kernels for transposing and ordering velocity space data into columns
   reorder_blocks_by_dimension_kernel<<<cudablocks, cudathreadsVECL, 0, stream>>> (
      blockData, // unified memory, incoming
      dev_blockDataOrdered[cpuThreadID],
      dev_cell_indices_to_id[cpuThreadID],
      host_totalColumns,
      LIDlist,
      columnData,
      // Also clears these vectors
      spatial_cell->BlocksRequired,
      spatial_cell->BlocksToAdd,
      spatial_cell->BlocksToRemove,
      spatial_cell->BlocksToMove
      );
   SSYNC;
   phiprof::stop("Reorder blocks by dimension");

   // Calculate target column extents
   // phiprof::start("Clear block lists");
   // spatial_cell->BlocksRequired->clear();
   // spatial_cell->BlocksToAdd->clear();
   // spatial_cell->BlocksToRemove->clear();
   // spatial_cell->BlocksToMove->clear();
   // if (doPrefetches) {
   //    spatial_cell->BlocksRequired->optimizeGPU(stream);
   //    spatial_cell->BlocksToAdd->optimizeGPU(stream);
   //    spatial_cell->BlocksToRemove->optimizeGPU(stream);
   // }
   // phiprof::stop("Clear block lists");

   phiprof::start("Evaluate column extents kernel");
   HANDLE_ERROR( cudaMemsetAsync(dev_returnLID, 0, sizeof(vmesh::LocalID), stream) );
   evaluate_column_extents_kernel<<<cudablocks, CUDATHREADS, 0, stream>>> (
      dimension,
      vmesh,
      columnData,
      columns,
      spatial_cell->BlocksRequired,
      spatial_cell->BlocksToAdd,
      spatial_cell->BlocksToRemove,
      GIDlist,
      dev_block_indices_to_id[cpuThreadID],
      intersection,
      intersection_di,
      intersection_dj,
      intersection_dk,
      Parameters::bailout_velocity_space_wall_margin,
      max_v_length,
      v_min,
      dv,
      dev_returnLID //dev_wallspace_margin_bailout_flag
      );
   SSYNC;
   // Check if we need to bailout due to hitting v-space edge
   HANDLE_ERROR( cudaMemcpyAsync(host_returnLID, dev_returnLID, sizeof(vmesh::LocalID), cudaMemcpyDeviceToHost, stream) );
   HANDLE_ERROR( cudaStreamSynchronize(stream) );
   if (host_returnLID[0] != 0) { //host_wallspace_margin_bailout_flag
      string message = "Some target blocks in acceleration are going to be less than ";
      message += std::to_string(Parameters::bailout_velocity_space_wall_margin);
      message += " blocks away from the current velocity space walls for population ";
      message += getObjectWrapper().particleSpecies[popID].name;
      message += " at CellID ";
      message += std::to_string(spatial_cell->parameters[CellParams::CELLID]);
      message += ". Consider expanding velocity space for that population.";
      bailout(true, message, __FILE__, __LINE__);
   }
   phiprof::stop("Evaluate column extents kernel");

   phiprof::start("Add and delete blocks");
   // Note: in this call, BlocksToMove is empty as we only grow the vspace size.
   spatial_cell->adjust_velocity_blocks_caller(popID);
   // Velocity space has all extra blocks added and/or removed for the transform target
   // and will not change shape anymore.
   const uint newNBlocks = vmesh->size();
   const uint bdsw3 = newNBlocks * WID3;
   phiprof::stop("Add and delete blocks");

   // Zero out target data on device (unified)
   // Put into second (high-priority) stream for concurrency
   phiprof::start("Memset ACC blocks to zero");
   HANDLE_ERROR( cudaMemsetAsync(blockData, 0, bdsw3*sizeof(Realf), priorityStream) );
   SSYNC;
   phiprof::stop("Memset ACC blocks to zero");

   phiprof::start("identify new block offsets kernel");
   identify_block_offsets_kernel<<<cudablocks, CUDATHREADS, 0, stream>>> (
      vmesh,
      columns,
      host_totalColumns,
      newNBlocks,
      dev_block_indices_to_id[cpuThreadID]
      );
   SSYNC;
   phiprof::stop("identify new block offsets kernel");

   HANDLE_ERROR( cudaStreamSynchronize(priorityStream) ); // Yes needed to ensure block data was zeroed
   phiprof::start("Semi-Lagrangian acceleration kernel");
   acceleration_kernel<<<cudablocks, cudathreadsVECL, 0, stream>>> (
      blockData,
      dev_blockDataOrdered[cpuThreadID],
      dev_cell_indices_to_id[cpuThreadID],
      columns,
      host_totalColumns,
      intersection,
      intersection_di,
      intersection_dj,
      intersection_dk,
      v_min,
      i_dv,
      dv,
      minValue,
      bdsw3
      );
   SSYNC;

   phiprof::stop("Semi-Lagrangian acceleration kernel");

   return true;
}
