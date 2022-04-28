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

#define NPP_MAXABS_32F ( 3.402823466e+38f )
#define NPP_MINABS_32F ( 1.175494351e-38f )
#define NPP_MAXABS_64F ( 1.7976931348623158e+308 )
#define NPP_MINABS_64F ( 2.2250738585072014e-308 )

#define i_pcolumnv_cuda(j, k, k_block, num_k_blocks) ( ((j) / ( VECL / WID)) * WID * ( num_k_blocks + 2) + (k) + ( k_block + 1 ) * WID )

#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ));
static void HandleError( cudaError_t err, const char *file, int line )
{
    if (err != cudaSuccess)
    {
        printf( "%s in %s at line %d\n", cudaGetErrorString( err ), file, line );
        exit( EXIT_FAILURE );
    }
}

#define MAXCPUTHREADS 64
Realf *dev_blockData[MAXCPUTHREADS];
Column *dev_columns[MAXCPUTHREADS];
int *dev_cell_indices_to_id[MAXCPUTHREADS];
Vec *dev_values[MAXCPUTHREADS];
Column *host_columns[MAXCPUTHREADS];
Vec *host_values[MAXCPUTHREADS];
bool isCudaAllocated = false;
float cudaAllocationMultiplier = 2.0; // This can be adjusted upwards in map_1d() based on what's needed
uint cudaMaxBlockCount = 0;

__host__ void cuda_acc_allocate_memory (
   uint cpuThreadID,
   uint maxBlockCount
   ) {
   // Mallocs should be in increments of 256 bytes. WID3 is at least 64, and len(Realf) is at least 4, so this is true.

   // The worst case scenario is with every block having content but no neighbours, creating up
   // to maxBlockCount columns with each needing three blocks (one value plus two for padding).
   // We also need to pad extra (x1.5?) due to the VDF deforming between acceleration Cartesian directions.
   // Here we make an 
   const uint maxColumnsPerCell = std::pow(maxBlockCount, 0.667) * cudaAllocationMultiplier; 
        // assumes symmetric smooth population2
   const uint maxTargetBlocksPerCell = maxBlockCount * cudaAllocationMultiplier;
   const uint maxSourceBlocksPerCell = maxBlockCount * cudaAllocationMultiplier;

   // Old version without checked max block count
   // const uint maxColumnsPerCell = ( MAX_BLOCKS_PER_DIM / 2 + 1) * MAX_BLOCKS_PER_DIM * MAX_BLOCKS_PER_DIM;
   // const uint maxTargetBlocksPerColumn = 3 * ( MAX_BLOCKS_PER_DIM / 2 + 1);
   // const uint maxTargetBlocksPerCell = maxTargetBlocksPerColumn * MAX_BLOCKS_PER_DIM * MAX_BLOCKS_PER_DIM;
   // const uint maxSourceBlocksPerCell = MAX_BLOCKS_PER_DIM * MAX_BLOCKS_PER_DIM * MAX_BLOCKS_PER_DIM;

   HANDLE_ERROR( cudaMalloc((void**)&dev_blockData[cpuThreadID], maxSourceBlocksPerCell * WID3 * sizeof(Realf) ) );
   HANDLE_ERROR( cudaMalloc((void**)&dev_columns[cpuThreadID], maxColumnsPerCell*sizeof(Column)) );
   HANDLE_ERROR( cudaMalloc((void**)&dev_cell_indices_to_id[cpuThreadID], 3*sizeof(int)) );
   HANDLE_ERROR( cudaMalloc((void**)&dev_values[cpuThreadID], maxTargetBlocksPerCell * (WID3 / VECL) * sizeof(Vec)) );

   // Also allocate and pin memory on host for faster transfers
   HANDLE_ERROR( cudaHostAlloc((void**)&host_columns[cpuThreadID], maxColumnsPerCell*sizeof(Column), cudaHostAllocPortable) );
   HANDLE_ERROR( cudaHostAlloc((void**)&host_values[cpuThreadID], maxTargetBlocksPerCell * (WID3 / VECL) * sizeof(Vec), cudaHostAllocPortable) );
   // Blockdata is pinned inside map_1d() in cpu_acc_map.cpp
}

__host__ void cuda_acc_deallocate_memory (
   uint cpuThreadID
   ) {
   HANDLE_ERROR( cudaFree(dev_blockData[cpuThreadID]) );
   HANDLE_ERROR( cudaFree(dev_cell_indices_to_id[cpuThreadID]) );
   HANDLE_ERROR( cudaFree(dev_columns[cpuThreadID]) );
   HANDLE_ERROR( cudaFree(dev_values[cpuThreadID]) );
   // Also de-allocate and unpin memory on host
   HANDLE_ERROR( cudaFreeHost(host_columns[cpuThreadID]) );
   HANDLE_ERROR( cudaFreeHost(host_values[cpuThreadID]) );
}

__global__ void acceleration_kernel
(
  Realf *dev_blockData,
  Column *dev_columns,
  Vec *dev_values,
  int *dev_cell_indices_to_id,
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
)
{
   int index = threadIdx.x;
   int column = blockIdx.x;
   // int nThreads = blockDim.x;
   // int cudaBlocks = gridDim.x; // = totalColums
   // //printf("totalcolumns %d blocks %d maxcolumns %d\n",totalColumns,cudaBlocks,maxcolumns);
   // if (nThreads != VECL) {
   //    printf("Warning! VECL not matching thread count for CUDA code!\n");
   // }

   if (column >= totalColumns) return;
   /* New threading with each warp/wavefront working on one vector */
   Realf v_r0 = ( (WID * dev_columns[column].kBegin) * dv + v_min);

   // i,j,k are relative to the order in which we copied data to the values array.
   // After this point in the k,j,i loops there should be no branches based on dimensions
   // Note that the i dimension is vectorized, and thus there are no loops over i
   // Iterate through the perpendicular directions of the column
   for (uint j = 0; j < WID; j += VECL/WID) {
      // If VECL=WID2 (WID=4, VECL=16, or WID=8, VECL=64, j==0)
      // This loop is still needed for e.g. Warp=VECL=32, WID2=64
      const vmesh::LocalID nblocks = dev_columns[column].nblocks;

      uint i_indices = index % WID;
      uint j_indices = j + index/WID;

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
      for (uint k=0; k < WID * nblocks; ++k)
      {
         // Compute reconstructions
         // Checked on 21.01.2022: Realv a[length] goes on the register despite being an array. Explicitly declaring it
         // as __shared__ had no impact on performance.
#ifdef ACC_SEMILAG_PLM
         Realv a[2];
         compute_plm_coeff(dev_values + dev_columns[column].valuesOffset + i_pcolumnv_cuda(j, 0, -1, nblocks), (k + WID), a, minValue, index);
#endif
#ifdef ACC_SEMILAG_PPM
         Realv a[3];
         compute_ppm_coeff(dev_values + dev_columns[column].valuesOffset + i_pcolumnv_cuda(j, 0, -1, nblocks), h4, (k + WID), a, minValue, index);
#endif
#ifdef ACC_SEMILAG_PQM
         Realv a[5];
         compute_pqm_coeff(dev_values + dev_columns[column].valuesOffset + i_pcolumnv_cuda(j, 0, -1, nblocks), h8, (k + WID), a, minValue, index);
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
}

void acceleration_1_glue
(
  Realf *blockData,
  Column *columns,
  Vec *values,
  uint cell_indices_to_id[],
  int totalColumns,
  int valuesSizeRequired,
  int bdsw3,
  Realv intersection,
  Realv intersection_di,
  Realv intersection_dj,
  Realv intersection_dk,
  Realv v_min,
  Realv i_dv,
  Realv dv,
  Realv minValue,
  const uint cuda_async_queue_id
)
{
   // cudaEvent_t start, stop;
   // cudaEventCreate(&start);
   // cudaEventCreate(&stop);
   // cudaEventRecord(start, 0);

   // Page lock (pin) host memory for faster async transfers
   //cudaHostRegister(h_ptr,bytes,cudaHostRegisterDefault);
   //cudaHostUnregister(h_ptr);
   //cudaHostRegister(cell_indices_to_id,3*sizeof(int),cudaHostRegisterDefault);
   // cudaHostRegister(columns,totalColumns*sizeof(Column),cudaHostRegisterDefault);
   // cudaHostRegister(values,valuesSizeRequired*sizeof(Vec),cudaHostRegisterDefault);
   // cudaHostRegister(blockData, bdsw3*sizeof(Realf),cudaHostRegisterDefault);

   cudaStream_t stream; //, stream_columns, stream_memset;
   cudaStreamCreate(&stream);

   // Now done in separate call:
   // Realf *dev_blockData;
   // Column *dev_columns;
   // int *dev_cell_indices_to_id;
   // Vec *dev_values;
   // Mallocs should be in increments of 256 bytes. WID3 is at least 64, and len(Realf) is at least 4, so this is true.
   // HANDLE_ERROR( cudaMalloc((void**)&dev_blockData, bdsw3*sizeof(Realf)) );
   // HANDLE_ERROR( cudaMalloc((void**)&dev_columns, totalColumns*sizeof(Column)) );
   // HANDLE_ERROR( cudaMalloc((void**)&dev_cell_indices_to_id, 3*sizeof(int)) );
   // HANDLE_ERROR( cudaMalloc((void**)&dev_values, valuesSizeRequired*sizeof(Vec)) );

   //HANDLE_ERROR( cudaMemcpy(dev_blockData, blockData, bdsw3*sizeof(Realf), cudaMemcpyHostToDevice) );
   // This target data is initialized to zero
   HANDLE_ERROR( cudaMemsetAsync(dev_blockData[cuda_async_queue_id], 0, bdsw3*sizeof(Realf), stream) );
   // Copy source data to device (async)
   HANDLE_ERROR( cudaMemcpyAsync(dev_cell_indices_to_id[cuda_async_queue_id], cell_indices_to_id, 3*sizeof(int), cudaMemcpyHostToDevice, stream) );
   HANDLE_ERROR( cudaMemcpyAsync(dev_columns[cuda_async_queue_id], columns, totalColumns*sizeof(Column), cudaMemcpyHostToDevice, stream) );
   HANDLE_ERROR( cudaMemcpyAsync(dev_values[cuda_async_queue_id], values, valuesSizeRequired*sizeof(Vec), cudaMemcpyHostToDevice, stream) );

   int threads = VECL; // equal to CUDATHREADS; NVIDIA: 32 AMD: 64
   int blocks = totalColumns; //CUDABLOCKS; 
   // NVIDIA: a100 64 stream multiprocessors? Blocks should be larger than this value.

   // Launch acceleration kernels
   acceleration_kernel<<<blocks, threads, 0, stream>>> (
         dev_blockData[cuda_async_queue_id],
         dev_columns[cuda_async_queue_id],
         dev_values[cuda_async_queue_id],
         dev_cell_indices_to_id[cuda_async_queue_id],
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

   // Copy data back to host
   HANDLE_ERROR( cudaMemcpyAsync(blockData, dev_blockData[cuda_async_queue_id], bdsw3*sizeof(Realf), cudaMemcpyDeviceToHost, stream) );

   cudaStreamSynchronize(stream);
   
   // cudaEventRecord(stop, stream);
   // cudaEventSynchronize(stop);
   // float elapsedTime;
   // cudaEventElapsedTime(&elapsedTime, start, stop);
   //printf("%.3f ms\n", elapsedTime);

   // Now done in separate call:
   // HANDLE_ERROR( cudaFree(dev_blockData) );
   // HANDLE_ERROR( cudaFree(dev_cell_indices_to_id) );
   // HANDLE_ERROR( cudaFree(dev_columns) );
   // HANDLE_ERROR( cudaFree(dev_values) );

   cudaStreamDestroy(stream);

   // Free page locks on host memory
   //cudaHostUnregister(cell_indices_to_id);
   // cudaHostUnregister(columns);
   // cudaHostUnregister(values);
   // cudaHostUnregister(blockData);
  
   return;
}
