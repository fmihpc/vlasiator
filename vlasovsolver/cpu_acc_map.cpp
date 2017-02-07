/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://vlasiator.fmi.fi/
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

#include <cmath>
#include <algorithm>
#include <utility>

#include "cpu_acc_sort_blocks.hpp"
#include "cpu_1d_pqm.hpp"
#include "cpu_1d_ppm.hpp"
#include "cpu_1d_plm.hpp"
#include "cpu_acc_map.hpp"

using namespace std;
using namespace spatial_cell;

#define MAX_BLOCKS_PER_DIM 100

//index in the temporary and padded column data values array. Each
//column has an empty block in ether end.
#define i_pcolumnv(j, k, k_block, num_k_blocks) ( ((j) / ( VECL / WID)) * WID * ( num_k_blocks + 2) + (k) + ( k_block + 1 ) * WID )
#define i_pcolumnv_b(planeVectorIndex, k, k_block, num_k_blocks) ( planeVectorIndex * WID * ( num_k_blocks + 2) + (k) + ( k_block + 1 ) * WID )

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

void loadColumnBlockData(SpatialCell* spatial_cell,vmesh::GlobalID* blocks,
                         vmesh::LocalID n_blocks,Vec* __restrict__ values,
                         const unsigned char * const cellid_transpose);

/*!
  Copies the data of a column into the values array, and zeroes the data.
  
  
  For dimension=0 data copy  we have rotated data
  i -> k
  j -> j
  k -> i
  For dimension=1 data copy  we have rotated data
  i -> i
  j -> k
  k -> j

 * @param blocks Array containing block global IDs.
 * @param n_blocks Number of blocks in array blocks.
*/
inline void loadColumnBlockData(//SpatialCell* spatial_cell,
        const vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& vmesh,
        vmesh::VelocityBlockContainer<vmesh::LocalID>& blockContainer,
        vmesh::GlobalID* blocks,
        vmesh::LocalID n_blocks,
        Vec* __restrict__ values,
        const unsigned char * const cellid_transpose) {

   Realv blockValues[WID3];   
   // first set the 0 values for the two empty blocks 
   // we store above and below the existing blocks

   for (uint k=0; k<WID; ++k) {
      for (uint j = 0; j < WID; j += VECL/WID){ 
         values[i_pcolumnv(j, k, -1, n_blocks)] = Vec(0);
         values[i_pcolumnv(j, k, n_blocks, n_blocks)] = Vec(0);
      }
   }

   // copy block data for all blocks
   for (vmesh::LocalID block_k=0; block_k<n_blocks; ++block_k) {
      Realf* __restrict__ data = blockContainer.getData(vmesh.getLocalID(blocks[block_k]));

      //  Copy volume averages of this block, taking into account the dimension shifting
      for (uint i=0; i<WID3; ++i) {
         blockValues[i] = data[cellid_transpose[i]];
      }
      for (uint i=0; i<WID3; ++i) {
         data[i]=0;
      }

      // now load values into the actual values table..
      uint offset =0;
      for (uint k=0; k<WID; ++k) {
         for(uint planeVector = 0; planeVector < VEC_PER_PLANE; planeVector++){
            // load data from blockValues which already has shifted dimensions
            values[i_pcolumnv_b(planeVector, k, block_k, n_blocks)].load(blockValues + offset);
            offset += VECL;
         }
      }
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
bool map_1d(vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& vmesh,
            vmesh::VelocityBlockContainer<vmesh::LocalID>& blockContainer,
            Realv intersection, Realv intersection_di, Realv intersection_dj,Realv intersection_dk,
            uint dimension) {
   no_subnormals();

   Realv dv,v_min;
   Realv is_temp;
   uint max_v_length;
   uint block_indices_to_id[3]; /*< used when computing id of target block */
   uint cell_indices_to_id[3]; /*< used when computing id of target cell in block*/
   unsigned char cellid_transpose[WID3]; /*< defines the transpose for the solver internal (transposed) id: i + j*WID + k*WID2 to actual one*/

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

   // init plane_index_to_id.
   for (uint k=0; k<WID; ++k) {
      for (uint j=0; j<WID; ++j) {
         for (uint i=0; i<WID; ++i) {
            const uint cell =
               i * cell_indices_to_id[0] +
               j * cell_indices_to_id[1] +
               k * cell_indices_to_id[2];
            cellid_transpose[ i + j * WID + k * WID2] = cell;
         }
      }
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
   
/*   
        values array used to store column data The max size is the worst
        case scenario with every second block having content, creating up
        to ( MAX_BLOCKS_PER_DIM / 2 + 1) columns with each needing three
        blocks (two for padding)
   */
   Vec values[(3 * ( MAX_BLOCKS_PER_DIM / 2 + 1)) * WID3 / VECL];
   // these two temporary variables are used to optimize access to target cells
   vmesh::LocalID previous_target_block = vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>::invalidLocalID();
   Realf *target_block_data = NULL;

   // loop over block column sets  (all columns along the dimension with the other dimensions being equal )
   for( uint setIndex=0; setIndex< setColumnOffsets.size(); ++setIndex) {

      //Load data into values array (this also zeroes the original data)
      uint valuesColumnOffset = 0; //offset to values array for data in a column in this set
      for(uint columnIndex = setColumnOffsets[setIndex]; columnIndex < setColumnOffsets[setIndex] + setNumColumns[setIndex] ; columnIndex ++){
         const vmesh::LocalID n_cblocks = columnNumBlocks[columnIndex];
         vmesh::GlobalID* cblocks = blocks + columnBlockOffsets[columnIndex]; //column blocks
         loadColumnBlockData(vmesh, blockContainer, cblocks, n_cblocks, values + valuesColumnOffset, cellid_transpose);
         valuesColumnOffset += (n_cblocks + 2) * (WID3/VECL); // there are WID3/VECL elements of type Vec per block
      }

      // loop over columns in set and do the mapping
      valuesColumnOffset = 0; //offset to values array for data in a column in this set
      for(uint columnIndex = setColumnOffsets[setIndex]; columnIndex < setColumnOffsets[setIndex] + setNumColumns[setIndex] ; columnIndex ++){
         const vmesh::LocalID n_cblocks = columnNumBlocks[columnIndex];
         vmesh::GlobalID* cblocks = blocks + columnBlockOffsets[columnIndex]; //column blocks
      
         // compute the common indices for this block column set
         //First block in column
         velocity_block_indices_t block_indices_begin;
         uint8_t refLevel;
         vmesh.getIndices(cblocks[0],refLevel,block_indices_begin[0],block_indices_begin[1],block_indices_begin[2]);
         uint temp;
         // Switch block indices according to dimensions, the algorithm has
         // been written for integrating along z.
         switch (dimension){
             case 0:
                /*i and k coordinates have been swapped*/
                temp=block_indices_begin[2];
                block_indices_begin[2]=block_indices_begin[0];
                block_indices_begin[0]=temp;
                break;
             case 1:
                /*in values j and k coordinates have been swapped*/
                temp=block_indices_begin[2];
                block_indices_begin[2]=block_indices_begin[1];
                block_indices_begin[1]=temp;
                break;
             case 2:
            break;
         }

         /*  i,j,k are now relative to the order in which we copied data to the values array. 
             After this point in the k,j,i loops there should be no branches based on dimensions
          
             Note that the i dimension is vectorized, and thus there are no loops over i
         */
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
       
            const int target_block_index_common =
               block_indices_begin[0] * block_indices_to_id[0] +
               block_indices_begin[1] * block_indices_to_id[1];
       
            /* 
               intersection_min is the intersection z coordinate (z after
               swaps that is) of the lowest possible z plane for each i,j
               index (i in vector)
            */
       
            const Vec intersection_min =
               intersection +
               (block_indices_begin[0] * WID + to_realv(i_indices)) * intersection_di + 
               (block_indices_begin[1] * WID + to_realv(j_indices)) * intersection_dj;

            /*compute some initial values, that are used to set up the
             * shifting of values as we go through all blocks in
             * order. See comments where they are shifted for
             * explanations of their meaning*/
            Vec v_r((WID * block_indices_begin[2]) * dv + v_min);
            Veci lagrangian_gk_r=truncate_to_int((v_r-intersection_min)/intersection_dk);
         
            // loop through all blocks in column and compute the mapping as integrals.
            for (uint k=0; k < WID * n_cblocks; ++k ){
               // Compute reconstructions 
               // values + i_pcolumnv(n_cblocks, -1, j, 0) is the starting point of the column data for fixed j
               // k + WID is the index where we have stored k index, WID amount of padding.
               #ifdef ACC_SEMILAG_PLM
               Vec a[2];
               compute_plm_coeff(values + valuesColumnOffset + i_pcolumnv(j, 0, -1, n_cblocks), k + WID , a);
               #endif
               #ifdef ACC_SEMILAG_PPM
               Vec a[3];
               compute_ppm_coeff(values + valuesColumnOffset + i_pcolumnv(j, 0, -1, n_cblocks), h4, k + WID, a);
               #endif
               #ifdef ACC_SEMILAG_PQM
               Vec a[5];
               compute_pqm_coeff(values + valuesColumnOffset + i_pcolumnv(j, 0, -1, n_cblocks), h8, k + WID, a);
               #endif
               
               // set the initial value for the integrand at the boundary at v = 0 
               // (in reduced cell units), this will be shifted to target_density_1, see below.
               Vec target_density_r(0.0);
               // v_l, v_r are the left and right velocity coordinates of source cell. Left is the old right.
               Vec v_l = v_r; 
               v_r += dv;
               // left(l) and right(r) k values (global index) in the target
               // Lagrangian grid, the intersecting cells. Again old right is new left.
               const Veci lagrangian_gk_l = lagrangian_gk_r;
               lagrangian_gk_r = truncate_to_int((v_r-intersection_min)/intersection_dk);
            
               Veci gk(lagrangian_gk_l);
               while (horizontal_or(gk <= lagrangian_gk_r)){
                  const Veci gk_div_WID = gk/WID;
                  const Veci gk_mod_WID = (gk - gk_div_WID * WID);
                  //the block of the Lagrangian cell to which we map
                  const Veci target_block(target_block_index_common + gk_div_WID * block_indices_to_id[2]);

                  //cell index in the target block 
                  const Veci target_cell(target_cell_index_common + gk_mod_WID * cell_indices_to_id[2]);
               
                  //the velocity between which we will integrate to put mass
                  //in the targe cell. If both v_r and v_l are in same cell
                  //then v_1,v_2 should be between v_l and v_r.
                  //v_1 and v_2 normalized to be between 0 and 1 in the cell.
                  //For vector elements where gk is already larger than needed (lagrangian_gk_r), v_2=v_1=v_r and thus the value is zero.

                  const Vec v_norm_r = (min(to_realv(gk + 1) * intersection_dk + intersection_min, v_r) - v_l) * i_dv;
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

                  // total value of integrand
                  const Vec target_density = target_density_r - target_density_l;

                  //store values, one element at a time
                  for (int target_i=0; target_i < VECL; ++target_i) {
                     const vmesh::GlobalID tblock = target_block[target_i];

                     // check that we are within sane limits. If gk is negative,
                     // or above blocks_per_dim * blockcells_per_dim then we
                     // are outside of the (possible) target grid.
                     // TODO, count losses if these are not fulfilled.
                     if (gk[target_i] >=0 && gk[target_i] < max_v_length * WID) {
                        if (previous_target_block != tblock) {
                           // The code inside this block is slower with the new AMR-related interface
                           previous_target_block = tblock;

                           // Get target block local ID. Attempt to create target
                           // block if it does not exist.
                           vmesh::LocalID tblockLID = vmesh.getLocalID(tblock);
                           if (tblockLID == vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>::invalidLocalID()) {
                              tblockLID = addVelocityBlock(tblock,vmesh,blockContainer);
                           }

                           // Get pointer to target block data. If target block does 
                           // not exist, values are added to null_block_data.
                           if (tblockLID == vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>::invalidLocalID()) {
                               target_block_data = blockContainer.getNullData();
                           } else {
                               target_block_data = blockContainer.getData(tblockLID);
                           }
                        }

                        // do the conversion from Realv to Realf here, faster than doing it in accumulation
                        const Realf tval = target_density[target_i];
                        const uint tcell = target_cell[target_i];
                        target_block_data[tcell] += tval;
                     }
                  } // for-loop over vector elements
                  gk++; //next iteration in while loop
               }
            } // for-loop over blocks
         }
         valuesColumnOffset += (n_cblocks + 2) * (WID3/VECL) ;// there are WID3/VECL elements of type Vec per block    
      }
   }
   delete [] blocks;
   return true;
}
