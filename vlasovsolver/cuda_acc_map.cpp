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
#include "cuda_acc_map_kernel.cuh"
#include "vec.h"
#include "../definitions.h"
#include "../object_wrapper.h"
#include "../cuda_context.cuh"

#include "cpu_face_estimates.hpp"
#include "cpu_1d_pqm.hpp"
#include "cpu_1d_ppm.hpp"
#include "cpu_1d_plm.hpp"

#include "cuda_acc_sort_blocks.hpp"

using namespace std;
using namespace spatial_cell;

/** Attempt to add the given velocity block to the given velocity mesh.
 * If the block was added to the mesh, its data is set to zero values and
 * velocity block parameters are calculated.
 * @param blockGID Global ID of the added velocity block.
 * @param vmesh Velocity mesh where the block is added.
 * @param blockContainer Velocity block data container.
 * @return Local ID of the added block. If the block was not added, the
 * local ID of the null velocity block is returned instead.*/
__host__ vmesh::LocalID addVelocityBlock2(const vmesh::GlobalID& blockGID,
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

   //nothing to do if no blocks
   vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& vmesh    = spatial_cell->get_velocity_mesh(popID);
   vmesh::VelocityBlockContainer<vmesh::LocalID>& blockContainer = spatial_cell->get_velocity_blocks(popID);
   Realf *blockData = blockContainer.getData();
   Realf *dev_blockData = blockContainer.dev_getData();
   uint blockDataN = blockContainer.size();
   if(vmesh.size() == 0) {
      return true;
   }

   // Thread id used for persistent device memory pointers
   const uint cuda_async_queue_id = omp_get_thread_num();
   int cudathreads = VECL; // equal to CUDATHREADS; NVIDIA: 32 AMD: 64

   Realv dv,v_min;
   Realv is_temp;
   uint max_v_length;
   /*< used when computing id of target block, 0 for compiler */
   uint block_indices_to_id[3] = {0, 0, 0};
   uint cell_indices_to_id[3] = {0, 0, 0};

   // Block data copying now handled from inside cuda_acc_semilag.cpp

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
   // Copy indexing information to device (async)
   HANDLE_ERROR( cudaMemcpyAsync(dev_cell_indices_to_id[cuda_async_queue_id], cell_indices_to_id, 3*sizeof(uint), cudaMemcpyHostToDevice, stream) );

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
   phiprof::start("sortBlockList");
   sortBlocklistByDimension(vmesh, dimension, GIDlist, LIDlist,
                            columnBlockOffsets, columnNumBlocks,
                            setColumnOffsets, setNumColumns);
   phiprof::stop("sortBlockList");

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
   //cudaMemcpyAsync(dev_GIDlist[cuda_async_queue_id], GIDlist, blockDataN*sizeof(vmesh::GlobalID), cudaMemcpyHostToDevice, stream);
   HANDLE_ERROR( cudaMemcpyAsync(dev_LIDlist[cuda_async_queue_id], LIDlist, blockDataN*sizeof(vmesh::LocalID), cudaMemcpyHostToDevice, stream) );
   HANDLE_ERROR( cudaMemcpyAsync(dev_columnNumBlocks[cuda_async_queue_id], columnNumBlocks.data(), totalColumns*sizeof(uint), cudaMemcpyHostToDevice, stream) );
   HANDLE_ERROR( cudaMemcpyAsync(dev_columnBlockOffsets[cuda_async_queue_id], columnBlockOffsets.data(), totalColumns*sizeof(uint), cudaMemcpyHostToDevice, stream) );

   // Launch kernels for transposing and ordering velocity space data into columns
   reorder_blocks_by_dimension_glue(
      dev_blockData,
      dev_blockDataOrdered[cuda_async_queue_id],
      dev_cell_indices_to_id[cuda_async_queue_id],
      totalColumns,
      dev_LIDlist[cuda_async_queue_id],
      dev_columnNumBlocks[cuda_async_queue_id],
      dev_columnBlockOffsets[cuda_async_queue_id],
      cudablocks,
      cudathreads,
      stream);
   // Unregister blockdata so that CPU can edit v-space to match requirements
   //cudaHostUnregister(blockData);

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
   phiprof::start("columnExtents");
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
            addVelocityBlock2(targetBlock, vmesh, blockContainer);

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
   phiprof::stop("columnExtents");

   // Velocity space has all extra blocks added and/or removed for the transform target
   // and will not change shape anymore.
   // Create empty velocity space on the GPU and fill it with zeros
   size_t blockDataSize = blockContainer.size();
   size_t bdsw3 = blockDataSize * WID3;
   // Page lock (pin) again host memory for faster async transfers after kernel has run
   //cudaHostRegister(blockData, bdsw3*sizeof(Realf),cudaHostRegisterDefault);
   // Zero out target data on device
   HANDLE_ERROR( cudaMemsetAsync(dev_blockData, 0, bdsw3*sizeof(Realf), stream) );

   // Update value of cudaAllocationMultiplier if necessary
   float ratio1 = (blockDataSize / cuda_acc_allocatedSize);
   float ratio2 = (totalColumns / cuda_acc_allocatedColumns);
   if ( (ratio1 > CUDA_ACC_SAFECTY_FACTOR)
        || (ratio2 > CUDA_ACC_SAFECTY_FACTOR) ) {
      // Need to increase ratio
      cudaAllocationMultiplier /= CUDA_ACC_SAFECTY_FACTOR;
      std::cerr<<"Increasing cudaAllocationMultiplier to "<<cudaAllocationMultiplier<<" with ratios "<<ratio1<<" "<<ratio2<<std::endl;
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

   // CALL CUDA FUNCTION WRAPPER/GLUE
   acceleration_1_glue(
      dev_blockData,
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
      bdsw3,
      cudablocks,
      cudathreads,
      stream
      );

   // Now done in subroutine called from cuda_acc_semilag.cpp
   // Copy data back to host
   //cudaMemcpyAsync(blockData, dev_blockData[cuda_async_queue_id], bdsw3*sizeof(Realf), cudaMemcpyDeviceToHost, stream);

   // Now done in cuda_acc_semilag.cpp
   //cudaStreamSynchronize(stream);
   // cudaEventRecord(stop, stream);
   // cudaEventSynchronize(stop);
   // float elapsedTime;
   // cudaEventElapsedTime(&elapsedTime, start, stop);
   //printf("%.3f ms\n", elapsedTime);
   //cudaStreamDestroy(stream);

   // Free page locks on host memory
   //cudaHostUnregister(blockData);

   delete[] GIDlist;
   delete[] LIDlist;
   delete[] columns;

   return true;
}
