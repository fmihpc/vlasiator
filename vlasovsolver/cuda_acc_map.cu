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

//#include "cpu_acc_map.hpp"
#include "cuda_acc_map.cuh"
#include "vec.h"
#include "../definitions.h"
#include "../object_wrapper.h"

#include "cpu_face_estimates.hpp"
#include "cpu_1d_pqm.hpp"
#include "cpu_1d_ppm.hpp"
#include "cpu_1d_plm.hpp"

#include "cuda_acc_sort_blocks.hpp"

#ifdef USE_CUDA
#include "device_launch_parameters.h"
#include "cuda.h"
#include "cuda_runtime.h"
#endif

#define NPP_MAXABS_32F ( 3.402823466e+38f )
#define NPP_MINABS_32F ( 1.175494351e-38f )
#define NPP_MAXABS_64F ( 1.7976931348623158e+308 )
#define NPP_MINABS_64F ( 2.2250738585072014e-308 )

#define i_pcolumnv_cuda(j, k, k_block, num_k_blocks) ( ((j) / ( VECL / WID)) * WID * ( num_k_blocks + 2) + (k) + ( k_block + 1 ) * WID )

#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ));
static void HandleError( cudaError_t err, const char *file, int line ) {
    if (err != cudaSuccess) {
        printf( "%s in %s at line %d\n", cudaGetErrorString( err ), file, line );
        exit( EXIT_FAILURE );
    }
}

// Allocate pointers for per-thread memory regions
#define MAXCPUTHREADS 64
Realf *dev_blockData[MAXCPUTHREADS];
Vec *dev_blockDataOrdered[MAXCPUTHREADS];
Column *dev_columns[MAXCPUTHREADS];
int *dev_cell_indices_to_id[MAXCPUTHREADS];
//vmesh::GlobalID *dev_GIDlist[MAXCPUTHREADS];
vmesh::LocalID *dev_LIDlist[MAXCPUTHREADS];

Column *host_columns[MAXCPUTHREADS];
vmesh::GlobalID *host_GIDlist[MAXCPUTHREADS];
vmesh::LocalID *host_LIDlist[MAXCPUTHREADS];


// Memory allocation flags and values.
// The allocation multiplier can be adjusted upwards if necessary.
float cudaAllocationMultiplier = 2.0;
bool isCudaAllocated = false;
uint cudaMaxBlockCount = 0;

//using namespace std;
//using namespace spatial_cell;

/** Attempt to add the given velocity block to the given velocity mesh.
 * If the block was added to the mesh, its data is set to zero values and
 * velocity block parameters are calculated.
 * @param blockGID Global ID of the added velocity block.
 * @param vmesh Velocity mesh where the block is added.
 * @param blockContainer Velocity block data container.
 * @return Local ID of the added block. If the block was not added, the
 * local ID of the null velocity block is returned instead.*/
__host__ vmesh::LocalID addVelocityBlock(const vmesh::GlobalID& blockGID,
        vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& vmesh,
        vmesh::VelocityBlockContainer<vmesh::LocalID>& blockContainer) {
    // Block insert will fail if the block already exists, or if
    // there are too many blocks in the velocity mesh.
    if (vmesh.push_back(blockGID) == false)
        return vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>::invalidLocalID();

    // Insert velocity block data, this will set values to 0.
    const vmesh::LocalID newBlockLID = blockContainer.push_back();

    #ifdef DEBUG_ACC
        bool ok = true;
        if (vmesh.size() != blockContainer.size()) ok = false;
        if (vmesh.getLocalID(blockGID) != newBlockLID) ok = false;
        if (ok == false) {
            stringstream ss;
            ss << "ERROR in acc: sizes " << vmesh.size() << ' ' << blockContainer.size() << endl;
            ss << "\t local IDs " << vmesh.getLocalID(blockGID) << " vs " << newBlockLID << endl;
            cerr << ss.str();
            exit(1);
        }
    #endif

    // Set block parameters:
    Real* parameters = blockContainer.getParameters(newBlockLID);
    vmesh.getBlockCoordinates(blockGID,parameters+BlockParams::VXCRD);
    vmesh.getCellSize(blockGID,parameters+BlockParams::DVX);
    return newBlockLID;
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

   // Old version without checked max block count (uses too much memory to be feasible)
   // const uint maxColumnsPerCell = ( MAX_BLOCKS_PER_DIM / 2 + 1) * MAX_BLOCKS_PER_DIM * MAX_BLOCKS_PER_DIM;
   // const uint maxTargetBlocksPerColumn = 3 * ( MAX_BLOCKS_PER_DIM / 2 + 1);
   // const uint maxTargetBlocksPerCell = maxTargetBlocksPerColumn * MAX_BLOCKS_PER_DIM * MAX_BLOCKS_PER_DIM;
   // const uint maxSourceBlocksPerCell = MAX_BLOCKS_PER_DIM * MAX_BLOCKS_PER_DIM * MAX_BLOCKS_PER_DIM;

   HANDLE_ERROR( cudaMalloc((void**)&dev_cell_indices_to_id[cpuThreadID], 3*sizeof(int)) );
   HANDLE_ERROR( cudaMalloc((void**)&dev_columns[cpuThreadID], maxColumnsPerCell*sizeof(Column)) );
   HANDLE_ERROR( cudaMalloc((void**)&dev_blockData[cpuThreadID], maxSourceBlocksPerCell * WID3 * sizeof(Realf) ) );
   HANDLE_ERROR( cudaMalloc((void**)&dev_blockDataOrdered[cpuThreadID], maxTargetBlocksPerCell * (WID3 / VECL) * sizeof(Vec)) );
   //HANDLE_ERROR( cudaMalloc((void**)&dev_GIDlist[cpuThreadID], maxSourceBlocksPerCell*sizeof(vmesh::GlobalID)) );
   HANDLE_ERROR( cudaMalloc((void**)&dev_LIDlist[cpuThreadID], maxSourceBlocksPerCell*sizeof(vmesh::LocalID)) );

   // Also allocate and pin memory on host for faster transfers
   HANDLE_ERROR( cudaHostAlloc((void**)&host_columns[cpuThreadID], maxColumnsPerCell*sizeof(Column), cudaHostAllocPortable) );
   HANDLE_ERROR( cudaHostAlloc((void**)&host_GIDlist[cpuThreadID], maxSourceBlocksPerCell*sizeof(vmesh::LocalID), cudaHostAllocPortable) );
   HANDLE_ERROR( cudaHostAlloc((void**)&host_LIDlist[cpuThreadID], maxSourceBlocksPerCell*sizeof(vmesh::LocalID), cudaHostAllocPortable) );
   // Blockdata is pinned inside cuda_acc_map_1d() in cuda_acc_map.cu
}

__host__ void cuda_acc_deallocate_memory (
   uint cpuThreadID
   ) {
   HANDLE_ERROR( cudaFree(dev_cell_indices_to_id[cpuThreadID]) );
   HANDLE_ERROR( cudaFree(dev_columns[cpuThreadID]) );
   HANDLE_ERROR( cudaFree(dev_blockData[cpuThreadID]) );
   HANDLE_ERROR( cudaFree(dev_blockDataOrdered[cpuThreadID]) );
   //HANDLE_ERROR( cudaFree(dev_GIDlist[cpuThreadID]) );
   HANDLE_ERROR( cudaFree(dev_LIDlist[cpuThreadID]) );

   // Also de-allocate and unpin memory on host
   HANDLE_ERROR( cudaFreeHost(host_columns[cpuThreadID]) );
   HANDLE_ERROR( cudaFreeHost(host_GIDlist[cpuThreadID]) );
   HANDLE_ERROR( cudaFreeHost(host_LIDlist[cpuThreadID]) );
}


__global__ void reorder_blocks_by_dimension_kernel(
   Realf *dev_blockData,
   Vec *dev_blockDataOrdered,
   int *dev_cell_indices_to_id,
   size_t blockDataN,
   vmesh::LocalID *dev_LIDlist
) {
   // Takes the contents of blockData, sorts it into blockDataOrdered, performing transposes as necessary
   const int nThreads = blockDim.x; // should be equal to VECL
   const int ti = threadIdx.x;
   const int start = blockIdx.x;
   const int cudaBlocks = gridDim.x;
   // if (nThreads != VECL) {
   //    printf("Warning! VECL not matching thread count for CUDA code!\n");
   // }

   for (uint blockIndex = start; blockIndex < blockDataN; blockIndex += cudaBlocks) {
      // Each block can span multiple VECLs (equal to cudathreads per block)
      for (uint j = 0; j < WID3; j += nThreads) {
         int input = ti + j;
         
         int input_2 = input / WID2; // last (slowest) index
         int input_1 = (input - input_2 * WID2) / WID;
         int input_0 = input - input_2 * WID2 - input_1 * WID; // first (fastest) index
         
         dev_blockDataOrdered[ blockIndex * WID3
                               + input_0 * dev_cell_indices_to_id[0]
                               + input_1 * dev_cell_indices_to_id[1]
                               + input_2 * dev_cell_indices_to_id[2] ]
            = dev_blockData[ dev_LIDlist[blockIndex] * WID3 + ti ];
      } // end loop j (vecs per block)
   } // end loop blockIndex
   // Note: this kernel does not memset blockData to zero.
   // A separate memsetasync call is required for that.
}

__global__ void acceleration_kernel(
  Realf *dev_blockData,
  Vec *dev_blockDataOrdered,
  int *dev_cell_indices_to_id,
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
                     const uint dimension) {

   //nothing to do if no blocks
   vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& vmesh    = spatial_cell->get_velocity_mesh(popID);
   vmesh::VelocityBlockContainer<vmesh::LocalID>& blockContainer = spatial_cell->get_velocity_blocks(popID);
   Realf *blockData = blockContainer.getData();
   uint blockDataN = blockContainer.size();
   if(vmesh.size() == 0) {
      return true;
   }

   // init GPU stream
   cudaStream_t stream; //, stream_columns, stream_memset;
   cudaStreamCreate(&stream);
   // Here should check which rank we have on the node, and activate the correct GPU
   // cudaSetDevice(gpuid);
   // The MPI standard does not have a common standard for this information.

   // Thread id used for persistent device memory pointers
   const uint cuda_async_queue_id = omp_get_thread_num();
   int cudathreads = VECL; // equal to CUDATHREADS; NVIDIA: 32 AMD: 64

   Realv dv,v_min;
   Realv is_temp;
   uint max_v_length;
   /*< used when computing id of target block, 0 for compiler */
   uint block_indices_to_id[3] = {0, 0, 0};
   uint cell_indices_to_id[3] = {0, 0, 0};


   // Page lock (pin) host memory for faster async transfers
   HANDLE_ERROR( cudaHostRegister(blockData, blockDataN*WID3*sizeof(Realf),cudaHostRegisterDefault) );
   // Launch async copy of velocity mesh contents to GPU
   HANDLE_ERROR( cudaMemcpyAsync(dev_blockData[cuda_async_queue_id], blockData, blockDataN*WID3*sizeof(Realf), cudaMemcpyHostToDevice, stream) ); 
   
   // Velocity grid refinement level, has no effect but is
   // needed in some vmesh::VelocityMesh function calls.
   const uint8_t REFLEVEL = 0;
   dv            = vmesh.getCellSize(REFLEVEL)[dimension];
   v_min         = vmesh.getMeshMinLimits()[dimension];
   max_v_length  = vmesh.getGridLength(REFLEVEL)[dimension];
   auto minValue = spatial_cell->getVelocityBlockMinValue(popID);

   switch (dimension) {
    case 0:
      /* i and k coordinates have been swapped*/

      /*swap intersection i and k coordinates*/
      is_temp=intersection_di;
      intersection_di=intersection_dk;
      intersection_dk=is_temp;

      /*set values in array that is used to convert block indices to id using a dot product*/
      block_indices_to_id[0] = vmesh.getGridLength(REFLEVEL)[0]*vmesh.getGridLength(REFLEVEL)[1];
      block_indices_to_id[1] = vmesh.getGridLength(REFLEVEL)[0];
      block_indices_to_id[2] = 1;

      /*set values in array that is used to convert block indices to id using a dot product*/
      cell_indices_to_id[0]=WID2;
      cell_indices_to_id[1]=WID;
      cell_indices_to_id[2]=1;
      break;
    case 1:
      /* j and k coordinates have been swapped*/

      /*swap intersection j and k coordinates*/
      is_temp=intersection_dj;
      intersection_dj=intersection_dk;
      intersection_dk=is_temp;

      /*set values in array that is used to convert block indices to id using a dot product*/
      block_indices_to_id[0]=1;
      block_indices_to_id[1] = vmesh.getGridLength(REFLEVEL)[0]*vmesh.getGridLength(REFLEVEL)[1];
      block_indices_to_id[2] = vmesh.getGridLength(REFLEVEL)[0];

      /*set values in array that is used to convert block indices to id using a dot product*/
      cell_indices_to_id[0]=1;
      cell_indices_to_id[1]=WID2;
      cell_indices_to_id[2]=WID;
      break;
    case 2:
      /*set values in array that is used to convert block indices to id using a dot product*/
      block_indices_to_id[0]=1;
      block_indices_to_id[1] = vmesh.getGridLength(REFLEVEL)[0];
      block_indices_to_id[2] = vmesh.getGridLength(REFLEVEL)[0]*vmesh.getGridLength(REFLEVEL)[1];

      // set values in array that is used to convert block indices to id using a dot product.
      cell_indices_to_id[0]=1;
      cell_indices_to_id[1]=WID;
      cell_indices_to_id[2]=WID2;
      break;
   }
   // Copy source data to device (async)
   HANDLE_ERROR( cudaMemcpyAsync(dev_cell_indices_to_id[cuda_async_queue_id], cell_indices_to_id, 3*sizeof(int), cudaMemcpyHostToDevice, stream) );

   const Realv i_dv=1.0/dv;

   // Declare these in pre-pinned (Page lock) host memory for faster async transfers
   vmesh::GlobalID *GIDlist = host_GIDlist[cuda_async_queue_id];
   vmesh::LocalID *LIDlist = host_LIDlist[cuda_async_queue_id];
   GIDlist = new vmesh::GlobalID[blockDataN]; // GIDs in dimension-order (length nBlocks)
   LIDlist = new vmesh::LocalID[blockDataN]; // LIDs in dimension-order (length nBlocks)

   // sort blocks according to dimension, and divide them into columns
   std::vector<uint> columnBlockOffsets; // indexes where columns start (in blocks, length totalColumns)
   std::vector<uint> columnNumBlocks; // length of column (in blocks, length totalColumns)
   std::vector<uint> setColumnOffsets; // index from columnBlockOffsets where new set of columns starts (length nColumnSets)
   std::vector<uint> setNumColumns; // how many columns in set of columns (length nColumnSets) 

   // CPU call for now (but dedicated CUDA version which also returns LIDlist)
   // This version actually could be run on the GPU as well as
   // it does not depend on hashmap calls, just calculating new
   // indices, sorting them, and going through the list.
   // Probably won't parallelize very well, though?
   sortBlocklistByDimension(vmesh, dimension, GIDlist, LIDlist,
                            columnBlockOffsets, columnNumBlocks,
                            setColumnOffsets, setNumColumns);
   
   // Calculate total sum of columns and total values size
   uint totalColumns = 0;
   uint valuesSizeRequired = 0;
   for(uint setIndex=0; setIndex< setColumnOffsets.size(); ++setIndex) {
      totalColumns += setNumColumns[setIndex];
      for(uint columnIndex = setColumnOffsets[setIndex]; columnIndex < setColumnOffsets[setIndex] + setNumColumns[setIndex] ; columnIndex ++){
         valuesSizeRequired += (columnNumBlocks[columnIndex] + 2) * WID3 / VECL;
      }
   }
   uint cudablocks = totalColumns;

   // memcopy LIDlist to device (GIDlist isn't needed here)
   //HANDLE_ERROR( cudaMemcpyAsync(dev_GIDlist[cuda_async_queue_id], GIDlist, blockDataN*sizeof(vmesh::GlobalID), cudaMemcpyHostToDevice, stream) );
   HANDLE_ERROR( cudaMemcpyAsync(dev_LIDlist[cuda_async_queue_id], LIDlist, blockDataN*sizeof(vmesh::LocalID), cudaMemcpyHostToDevice, stream) );

   // Launch kernels for transposing and ordering velocity space data into columns
   reorder_blocks_by_dimension_kernel<<<cudablocks, cudathreads, 0, stream>>> (
      dev_blockData[cuda_async_queue_id],
      dev_blockDataOrdered[cuda_async_queue_id],
      dev_cell_indices_to_id[cuda_async_queue_id],      
      blockDataN,
      dev_LIDlist[cuda_async_queue_id]
      );
   // Unregister blockdata so that CPU can edit v-space to match requirements
   HANDLE_ERROR( cudaHostUnregister(blockData) );
   
   // pointer to columns in memory
   Column *columns = host_columns[cuda_async_queue_id];
   columns = new Column[totalColumns];

   // Store offsets into columns
   uint valuesColumnOffset = 0;
   for( uint setIndex=0; setIndex< setColumnOffsets.size(); ++setIndex) {
      for(uint columnIndex = setColumnOffsets[setIndex]; columnIndex < setColumnOffsets[setIndex] + setNumColumns[setIndex] ; columnIndex ++){
         columns[columnIndex].nblocks = columnNumBlocks[columnIndex];
         columns[columnIndex].valuesOffset = valuesColumnOffset;
         if(valuesColumnOffset >= valuesSizeRequired) {
            cerr << "ERROR: Overflowing the values array (" << valuesColumnOffset << "> " << valuesSizeRequired << ") with column " << columnIndex << std::endl;
         }
         valuesColumnOffset += (columnNumBlocks[columnIndex] + 2) * (WID3/VECL); // there are WID3/VECL elements of type Vec per block
      }
   }
   
   // Calculate target column extents
   for( uint setIndex=0; setIndex< setColumnOffsets.size(); ++setIndex) {

      bool isTargetBlock[MAX_BLOCKS_PER_DIM]= {false};
      bool isSourceBlock[MAX_BLOCKS_PER_DIM]= {false};

      /*need x,y coordinate of this column set of blocks, take it from first
        block in first column*/
      //spatial_cell::velocity_block_indices_t setFirstBlockIndices;
      std::array<uint32_t,3> setFirstBlockIndices;
      uint8_t refLevel=0;
      vmesh.getIndices(GIDlist[columnBlockOffsets[setColumnOffsets[setIndex]]],
                       refLevel,
                       setFirstBlockIndices[0], setFirstBlockIndices[1], setFirstBlockIndices[2]);
      swapBlockIndices(setFirstBlockIndices, dimension);
      /*compute the maximum starting point of the lagrangian (target) grid
        (base level) within the 4 corner cells in this
        block. Needed for computig maximum extent of target column*/

      Realv max_intersectionMin = intersection +
                                      (setFirstBlockIndices[0] * WID + 0) * intersection_di +
                                      (setFirstBlockIndices[1] * WID + 0) * intersection_dj;
      max_intersectionMin =  std::max(max_intersectionMin,
                                      intersection +
                                      (setFirstBlockIndices[0] * WID + 0) * intersection_di +
                                      (setFirstBlockIndices[1] * WID + WID - 1) * intersection_dj);
      max_intersectionMin =  std::max(max_intersectionMin,
                                      intersection +
                                      (setFirstBlockIndices[0] * WID + WID - 1) * intersection_di +
                                      (setFirstBlockIndices[1] * WID + 0) * intersection_dj);
      max_intersectionMin =  std::max(max_intersectionMin,
                                      intersection +
                                      (setFirstBlockIndices[0] * WID + WID - 1) * intersection_di +
                                      (setFirstBlockIndices[1] * WID + WID - 1) * intersection_dj);

      Realv min_intersectionMin = intersection +
                                      (setFirstBlockIndices[0] * WID + 0) * intersection_di +
                                      (setFirstBlockIndices[1] * WID + 0) * intersection_dj;
      min_intersectionMin =  std::min(min_intersectionMin,
                                      intersection +
                                      (setFirstBlockIndices[0] * WID + 0) * intersection_di +
                                      (setFirstBlockIndices[1] * WID + WID - 1) * intersection_dj);
      min_intersectionMin =  std::min(min_intersectionMin,
                                      intersection +
                                      (setFirstBlockIndices[0] * WID + WID - 1) * intersection_di +
                                      (setFirstBlockIndices[1] * WID + 0) * intersection_dj);
      min_intersectionMin =  std::min(min_intersectionMin,
                                      intersection +
                                      (setFirstBlockIndices[0] * WID + WID - 1) * intersection_di +
                                      (setFirstBlockIndices[1] * WID + WID - 1) * intersection_dj);

      //now, record which blocks are target blocks
      for(uint columnIndex = setColumnOffsets[setIndex]; columnIndex < setColumnOffsets[setIndex] + setNumColumns[setIndex] ; columnIndex ++){
         const vmesh::LocalID n_cblocks = columnNumBlocks[columnIndex];
         vmesh::GlobalID* cblocks = GIDlist + columnBlockOffsets[columnIndex]; //column blocks
         //spatial_cell::velocity_block_indices_t firstBlockIndices;
         //spatial_cell::velocity_block_indices_t lastBlockIndices;
         std::array<uint32_t,3>  firstBlockIndices;
         std::array<uint32_t,3>  lastBlockIndices;
         vmesh.getIndices(cblocks[0],
                          refLevel,
                          firstBlockIndices[0], firstBlockIndices[1], firstBlockIndices[2]);
         vmesh.getIndices(cblocks[n_cblocks -1],
                          refLevel,
                          lastBlockIndices[0], lastBlockIndices[1], lastBlockIndices[2]);
         swapBlockIndices(firstBlockIndices, dimension);
         swapBlockIndices(lastBlockIndices, dimension);

         /*firstBlockV is in z the minimum velocity value of the lower
          * edge in source grid.
           *lastBlockV is in z the maximum velocity value of the upper
          * edge in source grid. */
         Realv firstBlockMinV = (WID * firstBlockIndices[2]) * dv + v_min;
         Realv lastBlockMaxV = (WID * (lastBlockIndices[2] + 1)) * dv + v_min;

         /*gk is now the k value in terms of cells in target
         grid. This distance between max_intersectionMin (so lagrangian
         plan, well max value here) and V of source grid, divided by
         intersection_dk to find out how many grid cells that is*/
         const int firstBlock_gk = (int)((firstBlockMinV - max_intersectionMin)/intersection_dk);
         const int lastBlock_gk = (int)((lastBlockMaxV - min_intersectionMin)/intersection_dk);

         int firstBlockIndexK = firstBlock_gk/WID;
         int lastBlockIndexK = lastBlock_gk/WID;

         //now enforce mesh limits for target column blocks
         firstBlockIndexK = (firstBlockIndexK >= 0)            ? firstBlockIndexK : 0;
         firstBlockIndexK = (firstBlockIndexK < max_v_length ) ? firstBlockIndexK : max_v_length - 1;
         lastBlockIndexK  = (lastBlockIndexK  >= 0)            ? lastBlockIndexK  : 0;
         lastBlockIndexK  = (lastBlockIndexK  < max_v_length ) ? lastBlockIndexK  : max_v_length - 1;
         if(firstBlockIndexK < Parameters::bailout_velocity_space_wall_margin
            || firstBlockIndexK >= max_v_length - Parameters::bailout_velocity_space_wall_margin
            || lastBlockIndexK < Parameters::bailout_velocity_space_wall_margin
            || lastBlockIndexK >= max_v_length - Parameters::bailout_velocity_space_wall_margin
         ) {
            string message = "Some target blocks in acceleration are going to be less than ";
            message += std::to_string(Parameters::bailout_velocity_space_wall_margin);
            message += " blocks away from the current velocity space walls for population ";
            message += getObjectWrapper().particleSpecies[popID].name;
            message += " at CellID ";
            message += std::to_string(spatial_cell->parameters[CellParams::CELLID]);
            message += ". Consider expanding velocity space for that population.";
            bailout(true, message, __FILE__, __LINE__);
         }
         
         //store source blocks
         for (uint blockK = firstBlockIndices[2]; blockK <= lastBlockIndices[2]; blockK++){
            isSourceBlock[blockK] = true;
         }

         //store target blocks
         for (uint blockK = (uint)firstBlockIndexK; blockK <= (uint)lastBlockIndexK; blockK++){
            isTargetBlock[blockK]=true;
         }

         // Set columns' transverse coordinates
         columns[columnIndex].i = setFirstBlockIndices[0];
         columns[columnIndex].j = setFirstBlockIndices[1];
         columns[columnIndex].kBegin = firstBlockIndices[2];

         //store also for each column firstBlockIndexK, and lastBlockIndexK
         columns[columnIndex].minBlockK = firstBlockIndexK;
         columns[columnIndex].maxBlockK = lastBlockIndexK;
      }

      //now add target blocks that do not yet exist and remove source blocks
      //that are not target blocks
      for (uint blockK = 0; blockK < MAX_BLOCKS_PER_DIM; blockK++){
         if(isTargetBlock[blockK] && !isSourceBlock[blockK] )  {
            const int targetBlock =
               setFirstBlockIndices[0] * block_indices_to_id[0] +
               setFirstBlockIndices[1] * block_indices_to_id[1] +
               blockK                  * block_indices_to_id[2];
            addVelocityBlock(targetBlock, vmesh, blockContainer);

         }
         if(!isTargetBlock[blockK] && isSourceBlock[blockK] )  {
            const int targetBlock =
               setFirstBlockIndices[0] * block_indices_to_id[0] +
               setFirstBlockIndices[1] * block_indices_to_id[1] +
               blockK                  * block_indices_to_id[2];

            spatial_cell->remove_velocity_block(targetBlock, popID);
         }
      }
   }

   // Velocity space has all extra blocks added and/or removed for the transform target
   // and will not change shape anymore.
   // Create empty velocity space on the GPU and fill it with zeros
   size_t blockDataSize = blockContainer.size();
   size_t bdsw3 = blockDataSize * WID3;

   // Page lock (pin) again host memory for faster async transfers after kernel has run
   HANDLE_ERROR( cudaHostRegister(blockData, bdsw3*sizeof(Realf),cudaHostRegisterDefault) );
   // Zero out target data on device
   HANDLE_ERROR( cudaMemsetAsync(dev_blockData[cuda_async_queue_id], 0, bdsw3*sizeof(Realf), stream) );
   
   // Update value of cudaAllocationMultiplier if necessary
   float ratio1 = (blockDataSize / cudaMaxBlockCount);
   float ratio2 = (totalColumns / std::pow(cudaMaxBlockCount, 0.667));
   float ratio3 = ( (valuesSizeRequired*VECL/WID3) / cudaMaxBlockCount);
   if ( (ratio1 > 0.75*cudaAllocationMultiplier)
        || (ratio2 > 0.75*cudaAllocationMultiplier)
        || (ratio3 > 0.75*cudaAllocationMultiplier) ) {
      // Need to increase ratio
      cudaAllocationMultiplier *= 1.25;
      std::cerr<<"Increasing cudaAllocationMultiplier to "<<cudaAllocationMultiplier<<" with ratios "<<ratio1<<" "<<ratio2<<" "<<ratio3<<std::endl;
   }
   // Decreasing the allocation multiplier is not safe unless we know the maximum value for bdsw3 over all local cells.
   
   // Now we iterate through target columns again, identifying their block offsets
   for( uint column=0; column < totalColumns; column++) {
      //cout << "Velocity column no. " << column << ": [";
      for (int blockK = columns[column].minBlockK; blockK <= columns[column].maxBlockK; blockK++) {
         const int targetBlock =
            columns[column].i * block_indices_to_id[0] +
            columns[column].j * block_indices_to_id[1] +
            blockK            * block_indices_to_id[2];
         // The below call accesses the hashmap (CPU only for now)
         const vmesh::LocalID tblockLID = vmesh.getLocalID(targetBlock);
         // Get pointer to target block data.
         if(tblockLID >= blockContainer.size()) {
            cerr << "Error: block for index [ " << columns[column].i << ", " << columns[column].j << ", " << blockK << "] has invalid blockID " << tblockLID << std::endl;
         }
         columns[column].targetBlockOffsets[blockK] = tblockLID*WID3;
         //cout << tblockLID*WID3 << " (" << columns[column].i << "," << columns[column].j << ", " << blockK << "), ";
      }
      //cout << "]" << std::endl;
   }

   // cudaEvent_t start, stop;
   // cudaEventCreate(&start);
   // cudaEventCreate(&stop);
   // cudaEventRecord(start, 0);

   // Copy column information to device (async)
   HANDLE_ERROR( cudaMemcpyAsync(dev_columns[cuda_async_queue_id], columns, totalColumns*sizeof(Column), cudaMemcpyHostToDevice, stream) );
   
   // NVIDIA: a100 64 stream multiprocessors? Blocks should be larger than this value.
   // Launch acceleration kernels
   acceleration_kernel<<<cudablocks, cudathreads, 0, stream>>> (
         dev_blockData[cuda_async_queue_id],
         dev_blockDataOrdered[cuda_async_queue_id],
         dev_cell_indices_to_id[cuda_async_queue_id],
         dev_columns[cuda_async_queue_id],
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
   cudaStreamDestroy(stream);

   // Free page locks on host memory
   HANDLE_ERROR( cudaHostUnregister(blockData) );
   
   return true;
}
