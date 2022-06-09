/*
 * This file is part of Vlasiator.
 * Copyright 2010-2022 Finnish Meteorological Institute
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

#include "cuda_acc_map_kernel.cuh"
#include "vec.h"
#include "../definitions.h"
#include "../cuda_context.cuh"
#include "cpu_face_estimates.hpp"
#include "cpu_1d_pqm.hpp"
#include "cpu_1d_ppm.hpp"
#include "cpu_1d_plm.hpp"

#include "device_launch_parameters.h"
#include "cuda.h"
#include "cuda_runtime.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#define NPP_MAXABS_32F ( 3.402823466e+38f )
#define NPP_MINABS_32F ( 1.175494351e-38f )
#define NPP_MAXABS_64F ( 1.7976931348623158e+308 )
#define NPP_MINABS_64F ( 2.2250738585072014e-308 )

#define i_pcolumnv_cuda(j, k, k_block, num_k_blocks) ( ((j) / ( VECL / WID)) * WID * ( num_k_blocks + 2) + (k) + ( k_block + 1 ) * WID )
#define i_pcolumnv_cuda_b(planeVectorIndex, k, k_block, num_k_blocks) ( planeVectorIndex * WID * ( num_k_blocks + 2) + (k) + ( k_block + 1 ) * WID )

// Allocate pointers for per-thread memory regions
//#define MAXCPUTHREADS 64 now in cuda_context.hpp
//Realf *dev_blockData[MAXCPUTHREADS];
Vec *dev_blockDataOrdered[MAXCPUTHREADS];
Column *dev_columns[MAXCPUTHREADS];
uint *dev_cell_indices_to_id[MAXCPUTHREADS];
//vmesh::GlobalID *dev_GIDlist[MAXCPUTHREADS];
vmesh::LocalID *dev_LIDlist[MAXCPUTHREADS];
uint *dev_columnNumBlocks[MAXCPUTHREADS];
uint *dev_columnBlockOffsets[MAXCPUTHREADS];

Column *host_columns[MAXCPUTHREADS];
vmesh::GlobalID *host_GIDlist[MAXCPUTHREADS];
vmesh::LocalID *host_LIDlist[MAXCPUTHREADS];


// Memory allocation flags and values.
// The allocation multiplier can be adjusted upwards if necessary.
Real cudaAllocationMultiplier = 2.0;
uint cuda_acc_allocatedSize = 0;
uint cuda_acc_allocatedColumns = 0;

__host__ void cuda_acc_allocate (
   uint maxBlockCount
   ) {
   // Always prepare for at least 500 blocks
   const uint maxBlocksPerCell = maxBlockCount > 500 ? maxBlockCount : 500;
   // Check if we already have allocated enough memory?
   if (cuda_acc_allocatedSize > maxBlocksPerCell * cudaAllocationMultiplier * CUDA_ACC_SAFECTY_FACTOR) {
      return;
   }
   // Deallocate before allocating new memory
   for (uint i=0; i<omp_get_max_threads(); ++i) {
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

   // The worst case scenario is with every block having content but no neighbours, creating up
   // to maxBlockCount columns with each needing three blocks (one value plus two for padding).
   // Here we make an educated guess of  up to two symmetric smooth populations
   const uint blockAllocationCount = maxBlockCount * cudaAllocationMultiplier;
   const uint maxColumnsPerCell = 2 * std::pow(maxBlockCount, 0.667) * cudaAllocationMultiplier;
   cuda_acc_allocatedSize = blockAllocationCount;
   cuda_acc_allocatedColumns = maxColumnsPerCell;

   HANDLE_ERROR( cudaMalloc((void**)&dev_cell_indices_to_id[cpuThreadID], 3*sizeof(uint)) );
   HANDLE_ERROR( cudaMalloc((void**)&dev_columns[cpuThreadID], maxColumnsPerCell*sizeof(Column)) );
   HANDLE_ERROR( cudaMalloc((void**)&dev_columnNumBlocks[cpuThreadID], maxColumnsPerCell*sizeof(uint)) );
   HANDLE_ERROR( cudaMalloc((void**)&dev_columnBlockOffsets[cpuThreadID], maxColumnsPerCell*sizeof(uint)) );
   HANDLE_ERROR( cudaMalloc((void**)&dev_blockDataOrdered[cpuThreadID], blockAllocationCount * (WID3 / VECL) * sizeof(Vec)) );
   HANDLE_ERROR( cudaMalloc((void**)&dev_LIDlist[cpuThreadID], blockAllocationCount*sizeof(vmesh::LocalID)) );

   // Old version without checked max block count (uses too much memory to be feasible)
   // const uint maxColumnsPerCell = ( MAX_BLOCKS_PER_DIM / 2 + 1) * MAX_BLOCKS_PER_DIM * MAX_BLOCKS_PER_DIM;
   // const uint maxTargetBlocksPerColumn = 3 * ( MAX_BLOCKS_PER_DIM / 2 + 1);
   // const uint maxTargetBlocksPerCell = maxTargetBlocksPerColumn * MAX_BLOCKS_PER_DIM * MAX_BLOCKS_PER_DIM;
   // const uint maxSourceBlocksPerCell = MAX_BLOCKS_PER_DIM * MAX_BLOCKS_PER_DIM * MAX_BLOCKS_PER_DIM;
   //HANDLE_ERROR( cudaMalloc((void**)&dev_blockData[cpuThreadID], blockAllocationCount * WID3 * sizeof(Realf) ) );
   //HANDLE_ERROR( cudaMalloc((void**)&dev_GIDlist[cpuThreadID], blockAllocationCount*sizeof(vmesh::GlobalID)) );

   // Also allocate and pin memory on host for faster transfers
   HANDLE_ERROR( cudaHostAlloc((void**)&host_columns[cpuThreadID], maxColumnsPerCell*sizeof(Column), cudaHostAllocPortable) );
   HANDLE_ERROR( cudaHostAlloc((void**)&host_GIDlist[cpuThreadID], blockAllocationCount*sizeof(vmesh::LocalID), cudaHostAllocPortable) );
   HANDLE_ERROR( cudaHostAlloc((void**)&host_LIDlist[cpuThreadID], blockAllocationCount*sizeof(vmesh::LocalID), cudaHostAllocPortable) );
   // Blockdata is pinned inside cuda_acc_map_1d() in cuda_acc_map.cu
   // printf("AA addrD %d -- %lu %lu %lu %lu\n",cpuThreadID,dev_cell_indices_to_id[cpuThreadID],dev_columns[cpuThreadID],dev_blockData[cpuThreadID],dev_blockDataOrdered[cpuThreadID]);
   // printf("AA addrH %d -- %lu %lu %lu %lu\n",cpuThreadID,&dev_cell_indices_to_id[cpuThreadID],&dev_columns[cpuThreadID],&dev_blockData[cpuThreadID],&dev_blockDataOrdered[cpuThreadID]);
 }

__host__ void cuda_acc_deallocate_memory (
   uint cpuThreadID
   ) {
   // printf("DD addrD %d -- %lu %lu %lu %lu\n",cpuThreadID,dev_cell_indices_to_id[cpuThreadID],dev_columns[cpuThreadID],dev_blockData[cpuThreadID],dev_blockDataOrdered[cpuThreadID]);
   // printf("DD addrH %d -- %lu %lu %lu %lu\n",cpuThreadID,&dev_cell_indices_to_id[cpuThreadID],&dev_columns[cpuThreadID],&dev_blockData[cpuThreadID],&dev_blockDataOrdered[cpuThreadID]);
   HANDLE_ERROR( cudaFree(dev_cell_indices_to_id[cpuThreadID]) );
   HANDLE_ERROR( cudaFree(dev_columns[cpuThreadID]) );
   // HANDLE_ERROR( cudaFree(dev_blockData[cpuThreadID]) );
   HANDLE_ERROR( cudaFree(dev_blockDataOrdered[cpuThreadID]) );
   //HANDLE_ERROR( cudaFree(dev_GIDlist[cpuThreadID]) );
   HANDLE_ERROR( cudaFree(dev_LIDlist[cpuThreadID]) );
   HANDLE_ERROR( cudaFree(dev_columnNumBlocks[cpuThreadID]) );
   HANDLE_ERROR( cudaFree(dev_columnBlockOffsets[cpuThreadID]) );

   // Also de-allocate and unpin memory on host
   HANDLE_ERROR( cudaFreeHost(host_columns[cpuThreadID]) );
   HANDLE_ERROR( cudaFreeHost(host_GIDlist[cpuThreadID]) );
   HANDLE_ERROR( cudaFreeHost(host_LIDlist[cpuThreadID]) );
   cuda_acc_allocatedSize = 0;
   cuda_acc_allocatedColumns = 0;
}


__global__ void reorder_blocks_by_dimension_kernel(
   Realf *dev_blockData,
   Vec *dev_blockDataOrdered,
   uint *dev_cell_indices_to_id,
   uint totalColumns,
   vmesh::LocalID *dev_LIDlist,
   uint *dev_columnNumBlocks,
   uint *dev_columnBlockOffsets
) {
   // Takes the contents of blockData, sorts it into blockDataOrdered,
   // performing transposes as necessary
   // Works column-per-column and adds the necessary one empty block at each end
   const int nThreads = blockDim.x; // should be equal to VECL
   const int ti = threadIdx.x;
   const int start = blockIdx.x;
   const int cudaBlocks = gridDim.x;
   // if (nThreads != VECL) {
   //    printf("Warning! VECL not matching thread count for CUDA code!\n");
   // }

   // Loop over columns in steps of cudaBlocks. Each cudaBlock deals with one column.
   for (uint iColumn = start; iColumn < totalColumns; iColumn += cudaBlocks) {
      if (iColumn >= totalColumns) break;
      uint inputOffset = dev_columnBlockOffsets[iColumn];
      uint outputOffset = (inputOffset + 2 * iColumn) * (WID3/VECL);
      uint columnLength = dev_columnNumBlocks[iColumn];

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

   // Note: this kernel does not memset dev_blockData to zero.
   // A separate memsetasync call is required for that.
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
  int bdsw3
) {
   const int index = threadIdx.x;
   const int column = blockIdx.x;
   // int nThreads = blockDim.x;
   // int cudaBlocks = gridDim.x; // = totalColums
   if (column >= totalColumns) return;
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

} // end semilag acc kernel


void acceleration_1_glue(
   Realf* dev_blockData,
   Vec* dev_blockDataOrdered,
   uint dev_cell_indices_to_id[],
   Column* dev_columns,
   const int totalColumns,
   const Realv intersection,
   const Realv intersection_di,
   const Realv intersection_dj,
   const Realv intersection_dk,
   const Realv v_min,
   const Realv i_dv,
   const Realv dv,
   const Realv minValue,
   const int bdsw3,
   const int cudablocks,
   const int cudathreads,
   cudaStream_t stream
) {
   // NVIDIA: a100 64 stream multiprocessors? Blocks should be larger than this value.
   // Launch acceleration kernels
   acceleration_kernel<<<cudablocks, cudathreads, 0, stream>>> (
      dev_blockData,
      dev_blockDataOrdered,
      dev_cell_indices_to_id,
      dev_columns,
      totalColumns,
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
   return;
}

void reorder_blocks_by_dimension_glue(
   Realf* dev_blockData,
   Vec* dev_blockDataOrdered,
   uint dev_cell_indices_to_id[],
   const uint totalColumns,
   uint* dev_LIDlist,
   uint* dev_columnNumBlocks,
   uint* dev_columnBlockOffsets,
   const int cudablocks,
   const int cudathreads,
   cudaStream_t stream
) {
   reorder_blocks_by_dimension_kernel<<<cudablocks, cudathreads, 0, stream>>> (
      dev_blockData,
      dev_blockDataOrdered,
      dev_cell_indices_to_id,
      totalColumns,
      dev_LIDlist,
      dev_columnNumBlocks,
      dev_columnBlockOffsets
      );
   return;
}
