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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gpu_acc_map.hpp"
#include "gpu_acc_sort_blocks.hpp"
#include "vec.h"
#include "../definitions.h"
#include "../object_wrapper.h"
#include "../arch/gpu_base.hpp"
#include "../spatial_cells/spatial_cell_gpu.hpp"

#include "cpu_face_estimates.hpp"
#include "cpu_1d_pqm.hpp"
#include "cpu_1d_ppm.hpp"
#include "cpu_1d_plm.hpp"

#ifdef DEBUG_VLASIATOR
   #ifndef DEBUG_ACC
   #define DEBUG_ACC
   #endif
#endif

#define i_pcolumnv_gpu(j, k, k_block, num_k_blocks) ( ((j) / ( VECL / WID)) * WID * ( num_k_blocks + 2) + (k) + ( k_block + 1 ) * WID )
#define i_pcolumnv_gpu_b(planeVectorIndex, k, k_block, num_k_blocks) ( planeVectorIndex * WID * ( num_k_blocks + 2) + (k) + ( k_block + 1 ) * WID )

using namespace std;
using namespace spatial_cell;

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

__global__ void __launch_bounds__(VECL,4) reorder_blocks_by_dimension_kernel(
   const vmesh::VelocityBlockContainer* __restrict__ blockContainer,
   Vec *gpu_blockDataOrdered,
   const uint* __restrict__ gpu_cell_indices_to_id,
   const vmesh::LocalID* __restrict__ gpu_LIDlist,
   const ColumnOffsets* __restrict__ columnData,
   const vmesh::LocalID valuesSizeRequired
) {
   // Takes the contents of blockData, sorts it into blockDataOrdered,
   // performing transposes as necessary
   // Works column-per-column and adds the necessary one empty block at each end
   const int ti = threadIdx.x;
   const uint iColumn = blockIdx.x;
   #ifdef DEBUG_ACC
   const int nThreads = blockDim.x; // should be equal to VECL
   if (nThreads != VECL) {
      if (ti==0) {
         printf("Warning! VECL not matching thread count for GPU kernel!\n");
      }
   }
   #endif
   // Each gpuBlock deals with one column.
   {
      const uint inputOffset = columnData->columnBlockOffsets.at(iColumn);
      const uint outputOffset = (inputOffset + 2 * iColumn) * (WID3/VECL);
      const uint columnLength = columnData->columnNumBlocks.at(iColumn);

      // Loop over column blocks
      for (uint b = 0; b < columnLength; b++) {
         // Slices
         for (uint k=0; k<WID; ++k) {
            // Each block slice can span multiple VECLs (equal to gputhreads per block)
            for (uint j = 0; j < WID; j += VECL/WID) {
               // full-block index
               const int input = k*WID2 + j*VECL + ti;
               // directional indices
               const int input_2 = input / WID2; // last (slowest) index
               const int input_1 = (input - input_2 * WID2) / WID; // medium index
               const int input_0 = input - input_2 * WID2 - input_1 * WID; // first (fastest) index
               // slice vector index
               const int jk = j / (VECL/WID);
               const int sourceindex = input_0 * gpu_cell_indices_to_id[0]
                  + input_1 * gpu_cell_indices_to_id[1]
                  + input_2 * gpu_cell_indices_to_id[2];

               #ifdef DEBUG_ACC
               assert((inputOffset + b) < blockContainer->size() && "reorder_blocks_by_dimension_kernel too large LID");
               assert((outputOffset + i_pcolumnv_gpu_b(jk, k, b, columnLength)) < valuesSizeRequired && "output error");
               #endif
               const vmesh::LocalID LID = gpu_LIDlist[inputOffset + b];
               const Realf* __restrict__ gpu_blockData = blockContainer->getData(LID);
               gpu_blockDataOrdered[outputOffset + i_pcolumnv_gpu_b(jk, k, b, columnLength)][ti]
                  = gpu_blockData[sourceindex ];

            } // end loop k (layers per block)
         } // end loop b (blocks per column)
      } // end loop j (vecs per layer)

      // Set first and last blocks to zero
      for (uint k=0; k<WID; ++k) {
         for (uint j = 0; j < WID; j += VECL/WID){
               int jk = j / (VECL/WID);
               #ifdef DEBUG_ACC
               assert((outputOffset + i_pcolumnv_gpu_b(jk, k, columnLength, columnLength)) < valuesSizeRequired && "output error");
               #endif
               gpu_blockDataOrdered[outputOffset + i_pcolumnv_gpu_b(jk, k, -1, columnLength)][ti] = 0.0;
               gpu_blockDataOrdered[outputOffset + i_pcolumnv_gpu_b(jk, k, columnLength, columnLength)][ti] = 0.0;
         }
      }
   } // end iColumn
   // Note: this kernel does not memset gpu_blockData to zero.
   // A separate memsetasync call is required for that.
}

// Serial kernel only to avoid page faults or prefetches
__global__ void __launch_bounds__(1,4) count_columns_kernel (
   const ColumnOffsets* __restrict__ gpu_columnData,
   vmesh::LocalID* returnLID, // gpu_totalColumns, gpu_valuesSizeRequired
   // Pass vectors for clearing
   split::SplitVector<vmesh::GlobalID> *list_with_replace_new,
   split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>* list_delete,
   split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>* list_to_replace,
   split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>* list_with_replace_old
   ) {
   // const int gpuBlocks = gridDim.x * gridDim.y * gridDim.z;
   // const int warpSize = blockDim.x * blockDim.y * blockDim.z;
   const int blocki = blockIdx.z*gridDim.x*gridDim.y + blockIdx.y*gridDim.x + blockIdx.x;
   const int ti = threadIdx.z*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;
   if ((blocki==0)&&(ti==0)) {
      for(uint setIndex=0; setIndex< gpu_columnData->setColumnOffsets.size(); ++setIndex) {
         returnLID[0] += gpu_columnData->setNumColumns[setIndex];
         for(uint columnIndex = gpu_columnData->setColumnOffsets[setIndex]; columnIndex < gpu_columnData->setColumnOffsets[setIndex] + gpu_columnData->setNumColumns[setIndex] ; columnIndex ++){
            returnLID[1] += (gpu_columnData->columnNumBlocks[columnIndex] + 2) * WID3 / VECL;
         }
      }
   }
   // Also clear these vectors
   if ((blocki==0)&&(ti==0)) {
      list_with_replace_new->clear();
      list_delete->clear();
      list_to_replace->clear();
      list_with_replace_old->clear();
   }
}

// Serial kernel only to avoid page faults or prefetches
__global__ void __launch_bounds__(1,4) offsets_into_columns_kernel(
   const ColumnOffsets* __restrict__ gpu_columnData,
   Column *gpu_columns,
   const uint valuesSizeRequired
   ) {
   // const int gpuBlocks = gridDim.x * gridDim.y * gridDim.z;
   // const int warpSize = blockDim.x * blockDim.y * blockDim.z;
   const int blocki = blockIdx.z*gridDim.x*gridDim.y + blockIdx.y*gridDim.x + blockIdx.x;
   const int ti = threadIdx.z*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;
   if ((blocki==0)&&(ti==0)) {
      uint valuesColumnOffset = 0;
      for( uint setIndex=0; setIndex< gpu_columnData->setColumnOffsets.size(); ++setIndex) {
         for (uint columnIndex = gpu_columnData->setColumnOffsets[setIndex]; columnIndex < gpu_columnData->setColumnOffsets[setIndex] + gpu_columnData->setNumColumns[setIndex] ; columnIndex ++){
            gpu_columns[columnIndex].nblocks = gpu_columnData->columnNumBlocks[columnIndex];
            gpu_columns[columnIndex].valuesOffset = valuesColumnOffset;
            if (valuesColumnOffset >= valuesSizeRequired) {
               printf("(ERROR: Overflowing the values array (%d > %d) with column %d\n",valuesColumnOffset,valuesSizeRequired,columnIndex);
            }
            valuesColumnOffset += (gpu_columnData->columnNumBlocks[columnIndex] + 2) * (WID3/VECL); // there are WID3/VECL elements of type Vec per block
         }
      }
   }
}

// Using columns, evaluate which blocks are target or source blocks
__global__ void __launch_bounds__(GPUTHREADS,4) evaluate_column_extents_kernel(
   const uint dimension,
   const vmesh::VelocityMesh* __restrict__ vmesh,
   const ColumnOffsets* __restrict__ gpu_columnData,
   Column *gpu_columns,
   split::SplitVector<vmesh::GlobalID> *list_with_replace_new,
   Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID> *dev_map_require,
   Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID> *dev_map_remove,
   const vmesh::GlobalID* __restrict__ GIDlist,
   const uint* __restrict__ gpu_block_indices_to_id,
   const Realf intersection,
   const Realf intersection_di,
   const Realf intersection_dj,
   const Realf intersection_dk,
   const int bailout_velocity_space_wall_margin,
   const int max_v_length,
   const Realf v_min,
   const Realf dv,
   uint *bailout_flag
   ) {
   const uint warpSize = blockDim.x * blockDim.y * blockDim.z;
   const uint blocki = blockIdx.z*gridDim.x*gridDim.y + blockIdx.y*gridDim.x + blockIdx.x;
   const uint ti = threadIdx.z*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;

   // Shared within all threads in one block
   __shared__ int isTargetBlock[MAX_BLOCKS_PER_DIM];
   __shared__ int isSourceBlock[MAX_BLOCKS_PER_DIM];
   const uint setIndex=blocki;
   if (setIndex < gpu_columnData->setColumnOffsets.size()) {

      // Clear flags used for this columnSet
      for(uint tti = 0; tti < MAX_BLOCKS_PER_DIM; tti += warpSize ) {
         const uint index = tti + ti;
         if (index < MAX_BLOCKS_PER_DIM) {
            isTargetBlock[index] = 0;
            isSourceBlock[index] = 0;
         }
      }
      __syncthreads();

      /*need x,y coordinate of this column set of blocks, take it from first
        block in first column*/
      vmesh::LocalID setFirstBlockIndices0,setFirstBlockIndices1,setFirstBlockIndices2;
      vmesh->getIndices(GIDlist[gpu_columnData->columnBlockOffsets[gpu_columnData->setColumnOffsets[setIndex]]],
                        setFirstBlockIndices0, setFirstBlockIndices1, setFirstBlockIndices2);
      swapBlockIndices(setFirstBlockIndices0,setFirstBlockIndices1,setFirstBlockIndices2,dimension);
      /*compute the maximum starting point of the lagrangian (target) grid
        (base level) within the 4 corner cells in this
        block. Needed for computig maximum extent of target column*/

      Realf max_intersectionMin = intersection +
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

      Realf min_intersectionMin = intersection +
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
      for (uint columnIndex = gpu_columnData->setColumnOffsets[setIndex];
           columnIndex < gpu_columnData->setColumnOffsets[setIndex] + gpu_columnData->setNumColumns[setIndex] ;
           ++columnIndex) {
         // Not parallelizing this at this level; not going to be many columns within a set

         // Abort all threads if vector capacity bailout
         if (bailout_flag[1] ) {
            return;
         }

         const vmesh::LocalID n_cblocks = gpu_columnData->columnNumBlocks[columnIndex];
         const vmesh::GlobalID* cblocks = GIDlist + gpu_columnData->columnBlockOffsets[columnIndex]; //column blocks
         vmesh::LocalID firstBlockIndices0,firstBlockIndices1,firstBlockIndices2;
         vmesh::LocalID lastBlockIndices0,lastBlockIndices1,lastBlockIndices2;
         vmesh->getIndices(cblocks[0],
                           firstBlockIndices0, firstBlockIndices1, firstBlockIndices2);
         vmesh->getIndices(cblocks[n_cblocks -1],
                           lastBlockIndices0, lastBlockIndices1, lastBlockIndices2);
         swapBlockIndices(firstBlockIndices0,firstBlockIndices1,firstBlockIndices2, dimension);
         swapBlockIndices(lastBlockIndices0,lastBlockIndices1,lastBlockIndices2, dimension);

         /* firstBlockV is in z the minimum velocity value of the lower
          *  edge in source grid.
          * lastBlockV is in z the maximum velocity value of the upper
          *  edge in source grid. */
         const Realf firstBlockMinV = (WID * firstBlockIndices2) * dv + v_min;
         const Realf lastBlockMaxV = (WID * (lastBlockIndices2 + 1)) * dv + v_min;

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
            // Pass bailout (hitting the wall) flag back to host
            if (ti==0) {
               bailout_flag[0] = 1;
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
            gpu_columns[columnIndex].i = setFirstBlockIndices0;
            gpu_columns[columnIndex].j = setFirstBlockIndices1;
            gpu_columns[columnIndex].kBegin = firstBlockIndices2;

            //store also for each column firstBlockIndexK, and lastBlockIndexK
            gpu_columns[columnIndex].minBlockK = firstBlockIndexK;
            gpu_columns[columnIndex].maxBlockK = lastBlockIndexK;
         }
      } // end loop over columns in set
      __syncthreads();

      for (uint blockT = 0; blockT < MAX_BLOCKS_PER_DIM; blockT +=warpSize) {
         const uint blockK = blockT + ti;
         // Not using warp accessors, as each thread has different block
         if (blockK < MAX_BLOCKS_PER_DIM) {
            if(isTargetBlock[blockK]!=0)  {
               const int targetBlock =
                  setFirstBlockIndices0 * gpu_block_indices_to_id[0] +
                  setFirstBlockIndices1 * gpu_block_indices_to_id[1] +
                  blockK                * gpu_block_indices_to_id[2];
               dev_map_require->set_element(targetBlock,vmesh->getLocalID(targetBlock));
               // if(!BlocksRequired->device_push_back(targetBlock)) {
               //    bailout_flag[1]=1;
               //    return;
               // }
            }
            if(isTargetBlock[blockK]!=0 && isSourceBlock[blockK]==0 )  {
               const int targetBlock =
                  setFirstBlockIndices0 * gpu_block_indices_to_id[0] +
                  setFirstBlockIndices1 * gpu_block_indices_to_id[1] +
                  blockK                * gpu_block_indices_to_id[2];
               if(!list_with_replace_new->device_push_back(targetBlock)) {
                  bailout_flag[1]=1; // out of capacity
               }

            }
            if(isTargetBlock[blockK]==0 && isSourceBlock[blockK]!=0 )  {
               const int targetBlock =
                  setFirstBlockIndices0 * gpu_block_indices_to_id[0] +
                  setFirstBlockIndices1 * gpu_block_indices_to_id[1] +
                  blockK                * gpu_block_indices_to_id[2];
               dev_map_remove->set_element(targetBlock,vmesh->getLocalID(targetBlock));
               // GPUTODO: could use device_insert to verify insertion, but not worth it
               // if(!list_to_replace->device_push_back(
               //       Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>(targetBlock,vmesh->getLocalID(targetBlock)))) {
               //    bailout_flag[2]=1; // out of capacity
               //    return;
               // }
            }
         } // block within MAX_BLOCKS_PER_DIM
      } // loop over all potential blocks
   } // if valid setIndex
}

__global__ void __launch_bounds__(VECL,4) acceleration_kernel(
   const vmesh::VelocityMesh* __restrict__ vmesh,
   vmesh::VelocityBlockContainer *blockContainer,
   const Vec* __restrict__ gpu_blockDataOrdered,
   const uint* __restrict__ gpu_cell_indices_to_id,
   const uint* __restrict__ gpu_block_indices_to_id,
   const Column* __restrict__ gpu_columns,
   const uint totalColumns, // not used
   const Realf intersection,
   const Realf intersection_di,
   const Realf intersection_dj,
   const Realf intersection_dk,
   const Realf v_min,
   const Realf i_dv,
   const Realf dv,
   const Realf minValue,
   const size_t invalidLID
) {
   //const uint gpuBlocks = gridDim.x * gridDim.y * gridDim.z;
   //const uint warpSize = blockDim.x * blockDim.y * blockDim.z;
   const uint blocki = blockIdx.z*gridDim.x*gridDim.y + blockIdx.y*gridDim.x + blockIdx.x;
   const uint w_tid = threadIdx.z*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;

   Realf *gpu_blockData = blockContainer->getData();
   const uint column = blocki;
   {
      /* New threading with each warp/wavefront working on one vector */
      const Realf v_r0 = ( (WID * gpu_columns[column].kBegin) * dv + v_min);

      // i,j,k are relative to the order in which we copied data to the values array.
      // After this point in the k,j,i loops there should be no branches based on dimensions
      // Note that the i dimension is vectorized, and thus there are no loops over i
      // Iterate through the perpendicular directions of the column
      for (uint j = 0; j < WID; j += VECL/WID) {
         // If VECL=WID2 (WID=4, VECL=16, or WID=8, VECL=64, then j==0)
         // This loop is still needed for e.g. Warp=VECL=32, WID2=64 (then j==0 or 4)
         const vmesh::LocalID nblocks = gpu_columns[column].nblocks;

         const uint i_indices = w_tid % WID;
         const uint j_indices = j + w_tid/WID;
         //int jk = j / (VECL/WID);

         const int target_cell_index_common =
            i_indices * gpu_cell_indices_to_id[0] +
            j_indices * gpu_cell_indices_to_id[1];
         const Realf intersection_min =
            intersection +
            (gpu_columns[column].i * WID + (Realf)i_indices) * intersection_di +
            (gpu_columns[column].j * WID + (Realf)j_indices) * intersection_dj;

         const Realf gk_intersection_min =
            intersection +
            (gpu_columns[column].i * WID + (Realf)( intersection_di > 0 ? 0 : WID-1 )) * intersection_di +
            (gpu_columns[column].j * WID + (Realf)( intersection_dj > 0 ? j : j+VECL/WID-1 )) * intersection_dj;
         const Realf gk_intersection_max =
            intersection +
            (gpu_columns[column].i * WID + (Realf)( intersection_di < 0 ? 0 : WID-1 )) * intersection_di +
            (gpu_columns[column].j * WID + (Realf)( intersection_dj < 0 ? j : j+VECL/WID-1 )) * intersection_dj;

         // loop through all perpendicular slices in column and compute the mapping as integrals.
         for (uint k=0; k < WID * nblocks; ++k) {
            // Compute reconstructions
            // Checked on 21.01.2022: Realf a[length] goes on the register despite being an array. Explicitly declaring it
            // as __shared__ had no impact on performance.
#ifdef ACC_SEMILAG_PLM
            Realf a[2];
            compute_plm_coeff(gpu_blockDataOrdered + gpu_columns[column].valuesOffset + i_pcolumnv_gpu(j, 0, -1, nblocks), (k + WID), a, minValue, w_tid);
#endif
#ifdef ACC_SEMILAG_PPM
            Realf a[3];
            compute_ppm_coeff(gpu_blockDataOrdered + gpu_columns[column].valuesOffset + i_pcolumnv_gpu(j, 0, -1, nblocks), h4, (k + WID), a, minValue, w_tid);
#endif
#ifdef ACC_SEMILAG_PQM
            Realf a[5];
            compute_pqm_coeff(gpu_blockDataOrdered + gpu_columns[column].valuesOffset + i_pcolumnv_gpu(j, 0, -1, nblocks), h8, (k + WID), a, minValue, w_tid);
#endif

            // set the initial value for the integrand at the boundary at v = 0
            // (in reduced cell units), this will be shifted to target_density_1, see below.
            Realf target_density_r = 0.0;

            const Realf v_r = v_r0  + (k+1)* dv;
            const Realf v_l = v_r0  + k* dv;
            const int lagrangian_gk_l = trunc((v_l-gk_intersection_max)/intersection_dk);
            const int lagrangian_gk_r = trunc((v_r-gk_intersection_min)/intersection_dk);

            //limits in lagrangian k for target column. Also take into
            //account limits of target column
            // Now all w_tids in the warp should have the same gk loop extents
            const int minGk = max(lagrangian_gk_l, int(gpu_columns[column].minBlockK * WID));
            const int maxGk = min(lagrangian_gk_r, int((gpu_columns[column].maxBlockK + 1) * WID - 1));
            // Run along the column and perform the polynomial reconstruction
            for(int gk = minGk; gk <= maxGk; gk++) {
               const int blockK = gk/WID;
               const int gk_mod_WID = (gk - blockK * WID);

               //the block of the Lagrangian cell to which we map
               //const int target_block(target_block_index_common + blockK * block_indices_to_id[2]);
               // This already contains the value index via target_cell_index_commom
               const int tcell(target_cell_index_common + gk_mod_WID * gpu_cell_indices_to_id[2]);
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

               const int targetBlock =
                  gpu_columns[column].i * gpu_block_indices_to_id[0] +
                  gpu_columns[column].j * gpu_block_indices_to_id[1] +
                  blockK                * gpu_block_indices_to_id[2];
               const vmesh::LocalID tblockLID = vmesh->getLocalID(targetBlock);
               // Using a warp search here seems to get only partial warp masks, resulting in an error
               //const vmesh::LocalID tblockLID = vmesh->warpGetLocalID(targetBlock, w_tid);
               if (isfinite(tval) && (tval>0) && (tblockLID != invalidLID) ) {
                  (&gpu_blockData[tblockLID*WID3])[tcell] += tval;
               }
            } // for loop over target k-indices of current source block
         } // for-loop over source blocks
      } //for loop over j index
   } // End this column
} // end semilag acc kernel

/*
   Here we map from the current time step grid, to a target grid which
   is the lagrangian departure grid (so th grid at timestep +dt,
   tracked backwards by -dt)
*/
__host__ bool gpu_acc_map_1d(spatial_cell::SpatialCell* spatial_cell,
                              const uint popID,
                              Real in_intersection,
                              Real in_intersection_di,
                              Real in_intersection_dj,
                              Real in_intersection_dk,
                              const uint dimension
   ) {
   // Ensure previous actions have completed?
   //CHK_ERR( gpuStreamSynchronize(stream) );
   gpuStream_t stream = gpu_getStream();

   // Conversion here:
   Realf intersection = (Realf)in_intersection;
   Realf intersection_di = (Realf)in_intersection_di;
   Realf intersection_dj = (Realf)in_intersection_dj;
   Realf intersection_dk = (Realf)in_intersection_dk;

   //spatial_cell->dev_upload_population(popID); // Should not be necessary.
   vmesh::VelocityMesh* vmesh    = spatial_cell->get_velocity_mesh(popID);
   vmesh::VelocityBlockContainer* blockContainer = spatial_cell->get_velocity_blocks(popID);
   vmesh::VelocityMesh* dev_vmesh    = spatial_cell->dev_get_velocity_mesh(popID);
   vmesh::VelocityBlockContainer* dev_blockContainer = spatial_cell->dev_get_velocity_blocks(popID);

   //nothing to do if no blocks
   vmesh::LocalID nBlocksBeforeAdjust = vmesh->size();
   if (nBlocksBeforeAdjust == 0) {
      return true;
   }

   auto minValue = spatial_cell->getVelocityBlockMinValue(popID);
   // These query velocity mesh parameters which are duplicated for both host and device
   const vmesh::LocalID D0 = vmesh->getGridLength()[0];
   const vmesh::LocalID D1 = vmesh->getGridLength()[1];
   //const vmesh::LocalID D2 = vmesh->getGridLength()[2];
   const Realf dv    = vmesh->getCellSize()[dimension];
   const Realf v_min = vmesh->getMeshMinLimits()[dimension];
   const int max_v_length  = (int)vmesh->getGridLength()[dimension];
   const Realf i_dv = 1.0/dv;

   // Thread id used for persistent device memory pointers
   const uint cpuThreadID = gpu_getThread();

   // Some kernels in here require the number of threads to be equal to VECL.
   // Future improvements would be to allow setting it directly to WID3.
   // Other kernels (not handling block data) can use GPUTHREADS which
   // is equal to NVIDIA: 32 or AMD: 64.

   /*< used when computing id of target block, 0 for compiler */
   uint block_indices_to_id[3] = {0, 0, 0};
   uint cell_indices_to_id[3] = {0, 0, 0};
   // 13.11.2023: for some reason these hostRegister calls say the memory is already registered.
   // CHK_ERR(gpuHostRegister(block_indices_to_id, 3*sizeof(uint),gpuHostRegisterPortable));
   // CHK_ERR(gpuHostRegister(cell_indices_to_id, 3*sizeof(uint),gpuHostRegisterPortable));

   Realf is_temp;
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
   // Ensure allocations
   spatial_cell->setReservation(popID, nBlocksBeforeAdjust);
   spatial_cell->applyReservation(popID);
   gpu_vlasov_allocate_perthread(cpuThreadID, nBlocksBeforeAdjust);

   // Copy indexing information to device (async)
   CHK_ERR( gpuMemcpyAsync(gpu_cell_indices_to_id[cpuThreadID], cell_indices_to_id, 3*sizeof(uint), gpuMemcpyHostToDevice, stream) );
   CHK_ERR( gpuMemcpyAsync(gpu_block_indices_to_id[cpuThreadID], block_indices_to_id, 3*sizeof(uint), gpuMemcpyHostToDevice, stream) );

   // Re-use maps from cell itself
   Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID> *map_require = spatial_cell->velocity_block_with_content_map;
   Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID> *map_remove = spatial_cell->velocity_block_with_no_content_map;
   Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID> *dev_map_require = spatial_cell->dev_velocity_block_with_content_map;
   Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID> *dev_map_remove = spatial_cell->dev_velocity_block_with_no_content_map;

   // pointers to device memory buffers
   vmesh::GlobalID *GIDlist = gpu_GIDlist[cpuThreadID];
   vmesh::LocalID *LIDlist = gpu_LIDlist[cpuThreadID];
   vmesh::GlobalID *BlocksID_mapped = gpu_BlocksID_mapped[cpuThreadID];
   vmesh::GlobalID *BlocksID_mapped_sorted = gpu_BlocksID_mapped[cpuThreadID];
   vmesh::LocalID *LIDlist_unsorted = gpu_LIDlist_unsorted[cpuThreadID];
   vmesh::LocalID *columnNBlocks = gpu_columnNBlocks[cpuThreadID];

   // Columndata is device construct but contains splitvectors
   ColumnOffsets *columnData = gpu_columnOffsetData[cpuThreadID];
   //columnData->prefetchDevice(stream);

   // These splitvectors are in host memory (only used for a clear call)
   split::SplitVector<vmesh::GlobalID> *list_with_replace_new = spatial_cell->list_with_replace_new;
   // These splitvectors are in device memory
   split::SplitVector<vmesh::GlobalID> *dev_list_with_replace_new = spatial_cell->dev_list_with_replace_new;
   split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>* dev_list_delete = spatial_cell->dev_list_delete;
   split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>* dev_list_to_replace = spatial_cell->dev_list_to_replace;
   split::SplitVector<Hashinator::hash_pair<vmesh::GlobalID,vmesh::LocalID>>* dev_list_with_replace_old = spatial_cell->dev_list_with_replace_old;

   // Call function for sorting block list and building columns from it.
   // Can probably be further optimized.
   //gpuStream_t priorityStream = gpu_getPriorityStream();
   CHK_ERR( gpuMemsetAsync(columnNBlocks, 0, gpu_acc_columnContainerSize*sizeof(vmesh::LocalID), stream) );
   //CHK_ERR( gpuStreamSynchronize(stream) ); // Yes needed because we use priority stream for block list sorting
   sortBlocklistByDimension(dev_vmesh,
                            nBlocksBeforeAdjust,
                            dimension,
                            BlocksID_mapped,
                            BlocksID_mapped_sorted,
                            GIDlist,
                            LIDlist_unsorted,
                            LIDlist,
                            columnNBlocks,
                            columnData,
                            cpuThreadID,
                            stream
      );
   CHK_ERR( gpuStreamSynchronize(stream) ); // Yes needed to get column data back to regular stream

   // Calculate total sum of columns and total values size
   CHK_ERR( gpuMemsetAsync(returnLID[cpuThreadID], 0, 2*sizeof(vmesh::LocalID), stream) );
   // this needs to be serial, but is fast.
   count_columns_kernel<<<1, 1, 0, stream>>> (
      columnData,
      returnLID[cpuThreadID], //gpu_totalColumns,gpu_valuesSizeRequired
      // Pass vectors for clearing
      dev_list_with_replace_new,
      dev_list_delete,
      dev_list_to_replace,
      dev_list_with_replace_old
      );
   CHK_ERR( gpuPeekAtLastError() );
   CHK_ERR( gpuMemcpyAsync(host_returnLID[cpuThreadID], returnLID[cpuThreadID], 2*sizeof(vmesh::LocalID), gpuMemcpyDeviceToHost, stream) );
   CHK_ERR( gpuStreamSynchronize(stream) );
   const vmesh::LocalID host_totalColumns = host_returnLID[cpuThreadID][0];
   const vmesh::LocalID host_valuesSizeRequired = host_returnLID[cpuThreadID][1];
   // Update tracker of maximum encountered column count
   if (gpu_acc_foundColumnsCount < host_totalColumns) {
      gpu_acc_foundColumnsCount = host_totalColumns;
   }

   // Create array of column objects
   Column host_columns[host_totalColumns];
   // and copy it into device memory
   Column *columns = gpu_columns[cpuThreadID];
   CHK_ERR( gpuMemcpyAsync(columns, &host_columns, host_totalColumns*sizeof(Column), gpuMemcpyHostToDevice, stream) );
   //SSYNC;

   // this needs to be serial, but is fast.
   offsets_into_columns_kernel<<<1, 1, 0, stream>>> (
      columnData,
      columns,
      host_valuesSizeRequired
      );
   CHK_ERR( gpuPeekAtLastError() );
   //CHK_ERR( gpuStreamSynchronize(stream) );

   // Launch kernels for transposing and ordering velocity space data into columns
   reorder_blocks_by_dimension_kernel<<<host_totalColumns, VECL, 0, stream>>> (
      dev_blockContainer,
      gpu_blockDataOrdered[cpuThreadID],
      gpu_cell_indices_to_id[cpuThreadID],
      LIDlist,
      columnData,
      host_valuesSizeRequired
      );
   CHK_ERR( gpuPeekAtLastError() );
   //CHK_ERR( gpuStreamSynchronize(stream) );

   // Calculate target column extents
   do {
      CHK_ERR( gpuMemsetAsync(returnLID[cpuThreadID], 0, 2*sizeof(vmesh::LocalID), stream) );
      map_require->clear<false>(Hashinator::targets::device,stream,std::pow(2,spatial_cell->vbwcl_sizePower));
      map_remove->clear<false>(Hashinator::targets::device,stream,std::pow(2,spatial_cell->vbwncl_sizePower));
      // Hashmap clear includes a stream sync
      //CHK_ERR( gpuStreamSynchronize(stream) );
      evaluate_column_extents_kernel<<<host_totalColumns, GPUTHREADS, 0, stream>>> (
         dimension,
         dev_vmesh,
         columnData,
         columns,
         dev_list_with_replace_new,
         dev_map_require,
         dev_map_remove,
         GIDlist,
         gpu_block_indices_to_id[cpuThreadID],
         intersection,
         intersection_di,
         intersection_dj,
         intersection_dk,
         Parameters::bailout_velocity_space_wall_margin,
         max_v_length,
         v_min,
         dv,
         returnLID[cpuThreadID] //gpu_bailout_flag:
                                // - element[0]: touching velspace wall
                                // - element[1]: splitvector list_with_replace_new capacity error
         );
      CHK_ERR( gpuPeekAtLastError() );
      // Check if we need to bailout due to hitting v-space edge
      CHK_ERR( gpuMemcpyAsync(host_returnLID[cpuThreadID], returnLID[cpuThreadID], 2*sizeof(vmesh::LocalID), gpuMemcpyDeviceToHost, stream) );
      CHK_ERR( gpuStreamSynchronize(stream) );
      if (host_returnLID[cpuThreadID][0] != 0) { //host_wallspace_margin_bailout_flag
         string message = "Some target blocks in acceleration are going to be less than ";
         message += std::to_string(Parameters::bailout_velocity_space_wall_margin);
         message += " blocks away from the current velocity space walls for population ";
         message += getObjectWrapper().particleSpecies[popID].name;
         message += " at CellID ";
         message += std::to_string(spatial_cell->parameters[CellParams::CELLID]);
         message += ". Consider expanding velocity space for that population.";
         bailout(true, message, __FILE__, __LINE__);
      }

      // Check whether we exceeded the column data splitVectors on the way
      if (host_returnLID[cpuThreadID][1] != 0) {
         // If so, recapacitate and try again.
         // We'll take at least our current velspace size (plus safety factor), or, if that wasn't enough,
         // twice what we had before.
         size_t newCapacity = (size_t)(spatial_cell->getReservation(popID)*BLOCK_ALLOCATION_FACTOR);
         //printf("column data recapacitate! %lu newCapacity\n",(long unsigned)newCapacity);
         list_with_replace_new->clear();
         spatial_cell->setReservation(popID, newCapacity);
         spatial_cell->applyReservation(popID);
      }
      // Loop until we return without an out-of-capacity error
   } while (host_returnLID[cpuThreadID][1] != 0);

   /** Rules used in extracting keys or elements from hashmaps
       Now these include passing pointers to GPU memory in order to evaluate
       nBlocksAfterAdjust without going via host. Pointers are copied by value.
   */
   const vmesh::GlobalID emptybucket = map_require->get_emptybucket();
   const vmesh::GlobalID tombstone   = map_require->get_tombstone();

   auto rule_delete_move = [emptybucket, tombstone, dev_map_remove, dev_list_with_replace_new, dev_vmesh]
      __host__ __device__(const Hashinator::hash_pair<vmesh::GlobalID, vmesh::LocalID>& kval) -> bool {
                              const vmesh::LocalID nBlocksAfterAdjust1 = dev_vmesh->size()
                                 + dev_list_with_replace_new->size() - dev_map_remove->size();
                              return kval.first != emptybucket &&
                                 kval.first != tombstone &&
                                 kval.second >= nBlocksAfterAdjust1 &&
                                 kval.second != vmesh::INVALID_LOCALID; };
   auto rule_to_replace = [emptybucket, tombstone, dev_map_remove, dev_list_with_replace_new, dev_vmesh]
      __host__ __device__(const Hashinator::hash_pair<vmesh::GlobalID, vmesh::LocalID>& kval) -> bool {
                             const vmesh::LocalID nBlocksAfterAdjust2 = dev_vmesh->size()
                                + dev_list_with_replace_new->size() - dev_map_remove->size();
                             return kval.first != emptybucket &&
                                kval.first != tombstone &&
                                kval.second < nBlocksAfterAdjust2; };

   // Additions are gathered directly into list instead of a map/set
   map_require->extractPatternLoop(*dev_list_with_replace_old, rule_delete_move, stream);
   map_remove->extractPatternLoop(*dev_list_delete, rule_delete_move, stream);
   map_remove->extractPatternLoop(*dev_list_to_replace, rule_to_replace, stream);
   //CHK_ERR( gpuStreamSynchronize(stream) );

   // Note: in this call, unless hitting v-space walls, we only grow the vspace size
   // and thus do not delete blocks or replace with old blocks.
   vmesh::LocalID nBlocksAfterAdjust = spatial_cell->adjust_velocity_blocks_caller(popID);
   // Velocity space has now all extra blocks added and/or removed for the transform target
   // and will not change shape anymore.
   spatial_cell->largestvmesh = spatial_cell->largestvmesh > nBlocksAfterAdjust ? spatial_cell->largestvmesh : nBlocksAfterAdjust;

   // Zero out target data on device (unified) (note, pointer needs to be re-fetched
   // here in case VBC size was increased)
   //GPUTODO: direct access to blockContainer getData causes page fault
   Realf *blockData = blockContainer->getData();
   CHK_ERR( gpuMemsetAsync(blockData, 0, nBlocksAfterAdjust*WID3*sizeof(Realf), stream) );
   //CHK_ERR( gpuStreamSynchronize(stream) );

   // GPUTODO: Adapt to work as VECL=WID3 instead of VECL=WID2
   acceleration_kernel<<<host_totalColumns, VECL, 0, stream>>> (
      dev_vmesh,
      dev_blockContainer,
      gpu_blockDataOrdered[cpuThreadID],
      gpu_cell_indices_to_id[cpuThreadID],
      gpu_block_indices_to_id[cpuThreadID],
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
      vmesh->invalidLocalID()
      );
   CHK_ERR( gpuPeekAtLastError() );
   CHK_ERR( gpuStreamSynchronize(stream) );

   return true;
}
