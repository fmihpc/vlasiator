/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 * 2017 CSC - IT center for Science
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

#include <math.h>
#include <algorithm>
#include <utility>

#include "vec.h"
#include "cpu_acc_sort_blocks.hpp"
#include "cpu_acc_load_blocks.hpp"
#include "cpu_1d_pqm.hpp"
#include "cpu_1d_ppm.hpp"
#include "cpu_1d_plm.hpp"
#include "cpu_acc_map.hpp"

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
vmesh::LocalID addVelocityBlock(const vmesh::GlobalID& blockGID,
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




void inline swapBlockIndices(velocity_block_indices_t &blockIndices, const uint dimension){

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

   TODO: parallelize with openMP over block-columns. If one also
   pre-creates new blocks in a separate loop first (serial operation),
   then the openmp parallization would scale well (better than over
   spatial cells), and would not need synchronization.

*/
bool map_1d(SpatialCell* spatial_cell,
            const uint popID,
            Realv intersection, Realv intersection_di, Realv intersection_dj,Realv intersection_dk,
            const uint dimension) {
   no_subnormals();

   Realv dv,v_min;
   Realv is_temp;
   uint max_v_length;
   uint block_indices_to_id[3] = {0, 0, 0}; /*< used when computing id of target block, 0 for compiler */
   uint cell_indices_to_id[3] = {0, 0, 0}; /*< used when computing id of target cell in block, 0 for compiler */

   vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& vmesh    = spatial_cell->get_velocity_mesh(popID);
   vmesh::VelocityBlockContainer<vmesh::LocalID>& blockContainer = spatial_cell->get_velocity_blocks(popID);

   //nothing to do if no blocks
   if(vmesh.size() == 0 )
      return true;


   // Velocity grid refinement level, has no effect but is
   // needed in some vmesh::VelocityMesh function calls.
   const uint8_t REFLEVEL = 0;

   dv            = vmesh.getCellSize(REFLEVEL)[dimension];
   v_min         = vmesh.getMeshMinLimits()[dimension];
   max_v_length  = vmesh.getGridLength(REFLEVEL)[dimension];

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

   const Realv i_dv=1.0/dv;

   // sort blocks according to dimension, and divide them into columns
   vmesh::LocalID* blocks = new vmesh::LocalID[vmesh.size()];
   std::vector<uint> columnBlockOffsets;
   std::vector<uint> columnNumBlocks;
   std::vector<uint> setColumnOffsets;
   std::vector<uint> setNumColumns;

   sortBlocklistByDimension(vmesh, dimension, blocks,
                            columnBlockOffsets, columnNumBlocks,
                            setColumnOffsets, setNumColumns);

   // Calculate total sum of columns and total values size
   int totalColumns = 0;
   int valuesSizeRequired = 0;

   for(uint setIndex=0; setIndex< setColumnOffsets.size(); ++setIndex) {
      totalColumns += setNumColumns[setIndex];
      for(uint columnIndex = setColumnOffsets[setIndex]; columnIndex < setColumnOffsets[setIndex] + setNumColumns[setIndex] ; columnIndex ++){
         valuesSizeRequired += (columnNumBlocks[columnIndex] + 2) * WID3 / VECL;
         }
   }

   // values array used to store column data. The actual columns index into this.
   Vec values[valuesSizeRequired];
   struct Column {
      int valuesOffset;              // Source data values
      size_t targetBlockOffsets[MAX_BLOCKS_PER_DIM]; // Target data array offsets
      int nblocks;                // Number of blocks in this column
      int minBlockK,maxBlockK;    // Column parallel coordinate limits
      int i,j;                    // Blocks' perpendicular coordinates
   } columns[totalColumns];

   // Iterate through all identified columns and shovel them into the values array.
   uint valuesColumnOffset = 0; //offset to values array for data in a column in this set
   for( uint setIndex=0; setIndex< setColumnOffsets.size(); ++setIndex) {
      //init
      //for (uint blockK = 0; blockK < MAX_BLOCKS_PER_DIM; blockK++){
      //   blockIndexToBlockDataIndex[blockK] = -1;
      //   isTargetBlock[blockK] = false;
      //   isSourceBlock[blockK] = false;
      //}

      //Load data into values array (this also zeroes the original data)
      for(uint columnIndex = setColumnOffsets[setIndex]; columnIndex < setColumnOffsets[setIndex] + setNumColumns[setIndex] ; columnIndex ++){
         const vmesh::LocalID n_cblocks = columnNumBlocks[columnIndex];
         vmesh::GlobalID* cblocks = blocks + columnBlockOffsets[columnIndex]; //column blocks
         loadColumnBlockData(vmesh, blockContainer, cblocks, n_cblocks, dimension, values + valuesColumnOffset);

         columns[columnIndex].valuesOffset = valuesColumnOffset;
         columns[columnIndex].nblocks = n_cblocks;

         valuesColumnOffset += (n_cblocks + 2) * (WID3/VECL); // there are WID3/VECL elements of type Vec per block
      }
   }

   // Now the velocity space has been emptied and all data is lying in the value array
   // (TODO: Double-check this for consistency?)
   // Upload that to the GPU.
   #pragma acc enter data copyin(values[:valuesSizeRequired]) async

   // Calculate target column extents
   for( uint setIndex=0; setIndex< setColumnOffsets.size(); ++setIndex) {

      bool isTargetBlock[MAX_BLOCKS_PER_DIM];
      bool isSourceBlock[MAX_BLOCKS_PER_DIM];
      for(int i=0; i< MAX_BLOCKS_PER_DIM; i++) {
         isTargetBlock[i] = false;
         isSourceBlock[i] = false;
      }

      /*need x,y coordinate of this column set of blocks, take it from first
        block in first column*/
      velocity_block_indices_t setFirstBlockIndices;
      uint8_t refLevel=0;
      vmesh.getIndices(blocks[columnBlockOffsets[setColumnOffsets[setIndex]]],
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
         vmesh::GlobalID* cblocks = blocks + columnBlockOffsets[columnIndex]; //column blocks
         velocity_block_indices_t firstBlockIndices;
         velocity_block_indices_t lastBlockIndices;
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
          * edge in source grid. Added 1.01*dv to account for unexpected issues*/
         double firstBlockMinV = (WID * firstBlockIndices[2]) * dv + v_min;
         double lastBlockMaxV = (WID * (lastBlockIndices[2] + 1)) * dv + v_min;

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

         //store source blocks
         for (uint blockK = firstBlockIndices[2]; blockK <= lastBlockIndices[2]; blockK++){
            isSourceBlock[blockK] = true;
         }

         //store target blocks
         for (int blockK = firstBlockIndexK; blockK <= lastBlockIndexK; blockK++){
            isTargetBlock[blockK]=true;
         }

         // Set columns' transverse coordinates
         columns[columnIndex].i = setFirstBlockIndices[0];
         columns[columnIndex].j = setFirstBlockIndices[1];

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
            //cerr << "Adding block for index [ " << setFirstBlockIndices[0] << ", " << setFirstBlockIndices[1] << ", " << blockK << "]" << std::endl;

         }
         if(!isTargetBlock[blockK] && isSourceBlock[blockK] )  {
            const int targetBlock =
               setFirstBlockIndices[0] * block_indices_to_id[0] +
               setFirstBlockIndices[1] * block_indices_to_id[1] +
               blockK                  * block_indices_to_id[2];

            spatial_cell->remove_velocity_block(targetBlock, popID);
            //cerr << "Remove block for index [ " << setFirstBlockIndices[0] << ", " << setFirstBlockIndices[1] << ", " << blockK << "]" << std::endl;
         }
      }
   }

   // Velocity space has all extra blocks added and/or removed for the transform target
   // and will not change shape anymore.
   // Create empty velocity space on the GPU and fill it with zeros
   Realf* blockData = blockContainer.getData();
   #pragma acc enter data create(blockData[:blockContainer.size()*WID3]) async
   #pragma acc parallel loop async
   for( int block=0; block < blockContainer.size()*WID3; block++) {
      blockData[block] = 0;
   }

   // Now we iterate through target columns again, identifying their block offsets
   for( uint column=0; column < totalColumns; column++) {
      for (int blockK = columns[column].minBlockK; blockK <= columns[column].maxBlockK; blockK++){
         const int targetBlock =
            columns[column].i * block_indices_to_id[0] +
            columns[column].j * block_indices_to_id[1] +
            blockK            * block_indices_to_id[2];
         const vmesh::LocalID tblockLID = vmesh.getLocalID(targetBlock);
         // Get pointer to target block data.
         if(tblockLID >= blockContainer.size()) {
            cerr << "Error: block for index [ " << columns[column].i << ", " << columns[column].j << ", " << blockK << "] has invalid blockID " << tblockLID << std::endl;
         }
         columns[column].targetBlockOffsets[blockK] = tblockLID*WID3;
      }
   }

   #pragma acc wait
   // loop over columns in set and do the mapping
   #pragma acc parallel loop copyin(columns[:totalColumns]) present(blockData[:blockContainer.size()*WID3], values[:valuesSizeRequired])
   for( uint column=0; column < totalColumns; column++) {
      const vmesh::LocalID nblocks = columns[column].nblocks;

      // i,j,k are relative to the order in which we copied data to the values array.
      // After this point in the k,j,i loops there should be no branches based on dimensions
      //
      // Note that the i dimension is vectorized, and thus there are no loops over i

      // Iterate through the perpendicular directions of the column
      #pragma acc loop
      for (uint j = 0; j < WID; j += VECL/WID){
         // create vectors with the i and j indices in the vector position on the plane.
         #if VECL == 4
         const Veci i_indices = Veci(0, 1, 2, 3);
         const Veci j_indices = Veci(j, j, j, j);
         #elif VECL == 8
         const Veci i_indices = Veci(0, 1, 2, 3,
                                     0, 1, 2, 3);
         const Veci j_indices = Veci(j, j, j, j,
                                     j + 1, j + 1, j + 1, j + 1);
         #elif VECL == 16
         const Veci i_indices = Veci(0, 1, 2, 3,
                                     0, 1, 2, 3,
                                     0, 1, 2, 3,
                                     0, 1, 2, 3);
         const Veci j_indices = Veci(j, j, j, j,
                                     j + 1, j + 1, j + 1, j + 1,
                                     j + 2, j + 2, j + 2, j + 2,
                                     j + 3, j + 3, j + 3, j + 3);
         #endif

         const Veci  target_cell_index_common =
            i_indices * cell_indices_to_id[0] +
            j_indices * cell_indices_to_id[1];

         //const int target_block_index_common =
         //   columns[column].i * block_indices_to_id[0] +
         //   columns[column].j * block_indices_to_id[1];

         // intersection_min is the intersection z coordinate (z after
         // swaps that is) of the lowest possible z plane for each i,j
         // index (i in vector)
         const Vec intersection_min =
            intersection +
            (columns[column].i * WID + to_realv(i_indices)) * intersection_di +
            (columns[column].j * WID + to_realv(j_indices)) * intersection_dj;

         /*compute some initial values, that are used to set up the
          * shifting of values as we go through all blocks in
          * order. See comments where they are shifted for
          * explanations of their meaning*/
         Vec v_r0((WID * columns[column].minBlockK) * dv + v_min);
         Vec lagrangian_v_r((v_r0-intersection_min)/intersection_dk);
         #if VECTORCLASS_H >= 20000
         Veci lagrangian_gk_r=truncatei(lagrangian_v_r);
         #else
         Veci lagrangian_gk_r=truncate_to_int(lagrangian_v_r);
         #endif

         // compute location of min and max, this does not change for one
         //  column (or even for this set of intersections, and can be used
         //  to quickly compute max and min later on*/
         //TODO, these can be computed much earlier, since they are
         //identiacal for each set of intersections
         int minGkIndex=0, maxGkIndex=0; // 0 for compiler
         {
            Realv maxV = std::numeric_limits<Realv>::min();
            Realv minV = std::numeric_limits<Realv>::max();
            for(int i = 0; i < VECL; i++) {
               if ( lagrangian_v_r[i] > maxV) {
                  maxV = lagrangian_v_r[i];
                  maxGkIndex = i;
               }
               if ( lagrangian_v_r[i] < minV) {
                  minV = lagrangian_v_r[i];
                  minGkIndex = i;
               }
            }
         }

            // loop through all blocks in column and compute the mapping as integrals.
            #pragma acc loop
            for (uint k=0; k < WID * nblocks; ++k ){
               // Compute reconstructions
               // values + i_pcolumnv(n_cblocks, -1, j, 0) is the starting point of the column data for fixed j
               // k + WID is the index where we have stored k index, WID amount of padding.
#ifdef ACC_SEMILAG_PLM
               Vec a[2];
               compute_plm_coeff(values + columns[column].valuesOffset + i_pcolumnv(j, 0, -1, nblocks), k + WID , a);
#endif
#ifdef ACC_SEMILAG_PPM
               Vec a[3];
               compute_ppm_coeff(values + columns[column].valuesOffset  + i_pcolumnv(j, 0, -1, nblocks), h4, k + WID, a);
#endif
#ifdef ACC_SEMILAG_PQM
               Vec a[5];
               compute_pqm_coeff(values + columns[column].valuesOffset  + i_pcolumnv(j, 0, -1, nblocks), h8, k + WID, a);
#endif

               // set the initial value for the integrand at the boundary at v = 0
               // (in reduced cell units), this will be shifted to target_density_1, see below.
               Vec target_density_r(0.0);
               // v_l, v_r are the left and right velocity coordinates of source cell.
               Vec v_r = v_r0  + (k+1)* dv;
               Vec v_l = v_r0  + k* dv;

               // left(l) and right(r) k values (global index) in the target
               // Lagrangian grid, the intersecting cells. Again old right is new left.
               const Veci lagrangian_gk_l = lagrangian_gk_r;
#if VECTORCLASS_H >= 20000
               lagrangian_gk_r = truncatei((v_r-intersection_min)/intersection_dk);
#else
               lagrangian_gk_r = truncate_to_int((v_r-intersection_min)/intersection_dk);

#endif

               //limits in lagrangian k for target column. Also take into
               //account limits of target column
               int minGk = max(int(lagrangian_gk_l[minGkIndex]), int(columns[column].minBlockK * WID));
               int maxGk = min(int(lagrangian_gk_r[maxGkIndex]), int((columns[column].maxBlockK + 1) * WID - 1));

               #pragma acc loop
               // Run along the column and perform the polynomial reconstruction
               for(int gk = minGk; gk <= maxGk; gk++){
                  const int blockK = gk/WID;
                  const int gk_mod_WID = (gk - blockK * WID);
                  //the block of the Lagrangian cell to which we map
                  //const int target_block(target_block_index_common + blockK * block_indices_to_id[2]);

                  //cell indices in the target block  (TODO: to be replaced by
                  //compile time generated scatter write operation)
                  const Veci target_cell(target_cell_index_common + gk_mod_WID * cell_indices_to_id[2]);

                  //the velocity between which we will integrate to put mass
                  //in the targe cell. If both v_r and v_l are in same cell
                  //then v_1,v_2 should be between v_l and v_r.
                  //v_1 and v_2 normalized to be between 0 and 1 in the cell.
                  //For vector elements where gk is already larger than needed (lagrangian_gk_r), v_2=v_1=v_r and thus the value is zero.
                  const Vec v_norm_r = (  min(  max( (gk + 1) * intersection_dk + intersection_min, v_l), v_r) - v_l) * i_dv;
                  /*shift, old right is new left*/
                  const Vec target_density_l = target_density_r;

                  // compute right integrand
                  #ifdef ACC_SEMILAG_PLM
                  target_density_r =
                     v_norm_r * ( a[0] + v_norm_r * a[1] );
                  #endif
                  #ifdef ACC_SEMILAG_PPM
                  target_density_r =
                     v_norm_r * ( a[0] + v_norm_r * ( a[1] + v_norm_r * a[2] ) );

                  #endif
                  #ifdef ACC_SEMILAG_PQM
                  target_density_r =
                     v_norm_r * ( a[0] + v_norm_r * ( a[1] + v_norm_r * ( a[2] + v_norm_r * ( a[3] + v_norm_r * a[4] ) ) ) );
                  #endif

                  //store values, one element at a time. All blocks
                  //have been created by now.
                  //TODO replace by vector version & scatter & gather operation


                  if(dimension == 2) {
                     Realf* targetDataPointer = blockData + columns[column].targetBlockOffsets[blockK] + j * cell_indices_to_id[1] + gk_mod_WID * cell_indices_to_id[2];
                     Vec targetData;
                     targetData.load_a(targetDataPointer);
                     targetData += target_density_r - target_density_l;
                     targetData.store_a(targetDataPointer);
                  }
                  else{
                     // total value of integrand
                     const Vec target_density = target_density_r - target_density_l;
   #pragma ivdep
   #pragma GCC ivdep
                     #pragma acc loop seq
                     for (int target_i=0; target_i < VECL; ++target_i) {
                        // do the conversion from Realv to Realf here, faster than doing it in accumulation
                        const Realf tval = target_density[target_i];
                        const uint tcell = target_cell[target_i];
                        (&blockData[columns[column].targetBlockOffsets[blockK]])[tcell] += tval;
                     }  // for-loop over vector elements
                  }
               } // for loop over target k-indices of current source block
            } // for-loop over source blocks
         } //for loop over j index
   } //for loop over columns

   #pragma acc exit data copyout(blockData[:blockContainer.size()*WID3]) delete(columns[:totalColumns]) delete(values[:valuesSizeRequired])

   delete [] blocks;
   return true;
}

