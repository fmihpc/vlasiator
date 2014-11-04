/*
This file is part of Vlasiator.
Copyright 2013, 2014 Finnish Meteorological Institute
*/

#ifndef CPU_ACC_MAP_H
#define CPU_ACC_MAP_H

#include  "vec4.h"
#include "algorithm"
#include "cmath"
#include "utility"
#include "common.h"
#include "spatial_cell.hpp"
#include "cpu_acc_sort_blocks.hpp"
#include "cpu_1d_pqm.hpp"
#include "cpu_1d_ppm.hpp"
#include "cpu_1d_plm.hpp"

#define MAX_BLOCKS_PER_DIM 100

//index in the temporary and padded column data values array. Each
//column has an empty block in ether end.
#define i_pcolumnv(num_k_blocks, k_block, j, k) ( (j) * WID * (num_k_blocks + 2) +  (k) + ( k_block + 1 ) * WID )



using namespace std;
using namespace spatial_cell;

/*!
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
inline void load_column_block_data(SpatialCell* spatial_cell,vmesh::GlobalID* blocks,vmesh::LocalID n_blocks, Vec4 * __restrict__ values, int dimension){
   uint cell_indices_to_id[3];
   switch (dimension){
       case 0:
          /* i and k coordinates have been swapped*/
          cell_indices_to_id[0]=WID2;
          cell_indices_to_id[1]=WID;
          cell_indices_to_id[2]=1;
          break;
       case 1:
          /* i and k coordinates have been swapped*/
          cell_indices_to_id[0]=1;
          cell_indices_to_id[1]=WID2;
          cell_indices_to_id[2]=WID;
          break;
       case 2:
          cell_indices_to_id[0]=1;
          cell_indices_to_id[1]=WID;
          cell_indices_to_id[2]=WID2;
          break;
       default:
          //same as for dimension 2, mostly here to get rid of compiler warning
          cell_indices_to_id[0]=1;
          cell_indices_to_id[1]=1;
          cell_indices_to_id[2]=1;
          cerr << "Dimension argument wrong: " << dimension << " at " << __FILE__ << ":" << __LINE__ << endl;
          exit(1);
          break;
   }
      
   
   /*first set the 0 values fot the two empty blocks we store above and below the existing blosk*/

   for (uint j=0; j<WID; ++j) {
      for (uint k=0; k<WID; ++k) {
         values[i_pcolumnv(n_blocks,-1,j,k)] = Vec4(0.0);
         values[i_pcolumnv(n_blocks,n_blocks,j,k)] = Vec4(0.0);
      }
   }

   /*copy block data for all blocks*/
   for (vmesh::LocalID block_k=0; block_k<n_blocks; ++block_k) {
      const vmesh::LocalID blockLID = spatial_cell->get_velocity_block_local_id(blocks[block_k]);
      Realf* __restrict__ fx = spatial_cell->get_fx(blockLID);

      //  Copy volume averages of this block, taking into account the dimension shifting
      for (uint j=0; j<WID; ++j) {
         for (uint k=0; k<WID; ++k) {
            for (uint i=0; i<WID; ++i) {
               const uint cell =
                  i * cell_indices_to_id[0] +
                  j * cell_indices_to_id[1] +
                  k * cell_indices_to_id[2];
               values[i_pcolumnv(n_blocks,block_k,j,k)].insert(i,(Real)fx[cell]);
            }
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

bool map_1d(SpatialCell* spatial_cell,
                    Real intersection, Real intersection_di, Real intersection_dj,Real intersection_dk,
                    uint dimension ) {
   /*
     Move densities from data to fx and clear data, to prepare for mapping
   */
   for (unsigned int cell=0; cell < VELOCITY_BLOCK_LENGTH*spatial_cell->get_number_of_velocity_blocks(); ++cell) {
      //copy data to fx for solvers, and set data to zero as we will map new values there
      spatial_cell->get_fx()[cell]   = spatial_cell->get_data()[cell];
      spatial_cell->get_data()[cell] = 0.0;
   }

   Real dv,v_min;
   Real is_temp;
   uint max_v_length;
   uint block_indices_to_id[3]; /*< used when computing id of target block */
   uint cell_indices_to_id[3]; /*< used when computing id of target cell in block*/
   switch (dimension) {
    case 0:
      /* i and k coordinates have been swapped*/
      /*set cell size in dimension direction*/
      dv = SpatialCell::get_velocity_base_grid_cell_size()[0]; 
      v_min = SpatialCell::get_velocity_grid_min_limits()[0];
      max_v_length = SpatialCell::get_velocity_base_grid_length()[0];

      /*swap intersection i and k coordinates*/
      is_temp=intersection_di;
      intersection_di=intersection_dk;
      intersection_dk=is_temp;
      
      /*set values in array that is used to transfer blockindices to id using a dot product*/
      block_indices_to_id[0]=SpatialCell::get_velocity_base_grid_length()[0] * SpatialCell::get_velocity_base_grid_length()[1];
      block_indices_to_id[1]=SpatialCell::get_velocity_base_grid_length()[0];
      block_indices_to_id[2]=1;
      
      /*set values in array that is used to transfer blockindices to id using a dot product*/
      cell_indices_to_id[0]=WID2;
      cell_indices_to_id[1]=WID;
      cell_indices_to_id[2]=1;
      break;
    case 1:
      /* j and k coordinates have been swapped*/
      /*set cell size in dimension direction*/
      dv=SpatialCell::get_velocity_base_grid_cell_size()[1];
      v_min = SpatialCell::get_velocity_grid_min_limits()[1];
      max_v_length = SpatialCell::get_velocity_base_grid_length()[1];
      
      /*swap intersection j and k coordinates*/
      is_temp=intersection_dj;
      intersection_dj=intersection_dk;
      intersection_dk=is_temp;
      
      /*set values in array that is used to transfer blockindices to id using a dot product*/
      block_indices_to_id[0]=1;
      block_indices_to_id[1]=SpatialCell::get_velocity_base_grid_length()[0] * SpatialCell::get_velocity_base_grid_length()[1];
      block_indices_to_id[2]=SpatialCell::get_velocity_base_grid_length()[0];
      
      /*set values in array that is used to transfer blockindices to id using a dot product*/
      cell_indices_to_id[0]=1;
      cell_indices_to_id[1]=WID2;
      cell_indices_to_id[2]=WID;
      break;
    case 2:
      /*set cell size in dimension direction*/
      dv=SpatialCell::get_velocity_base_grid_cell_size()[2];
      v_min = SpatialCell::get_velocity_grid_min_limits()[2];
      max_v_length = SpatialCell::get_velocity_base_grid_length()[2];
      
      /*set values in array that is used to transfer blockindices to id using a dot product*/
      block_indices_to_id[0]=1;
      block_indices_to_id[1]=SpatialCell::get_velocity_base_grid_length()[0];
      block_indices_to_id[2]=SpatialCell::get_velocity_base_grid_length()[0] * SpatialCell::get_velocity_base_grid_length()[1];
      
      /*set values in array that is used to transfer blockindices to id using a dot product*/
      cell_indices_to_id[0]=1;
      cell_indices_to_id[1]=WID;
      cell_indices_to_id[2]=WID2;
      break;
   }
   const Real i_dv=1.0/dv;
   
   /*sort blocks according to dimension, and divide them into columns*/
   uint* blocks=new uint[spatial_cell->get_number_of_velocity_blocks()];
   std::vector<uint> block_column_offsets;
   std::vector<uint> block_column_lengths;
   sort_blocklist_by_dimension( spatial_cell, dimension, blocks, block_column_offsets, block_column_lengths);

   /*values array used to store column data*/
   Vec4 values[(MAX_BLOCKS_PER_DIM + 2) * WID2];
   /*these two temporary variables are used to optimize access to target cells*/
   uint previous_target_block = SpatialCell::invalid_global_id();
   Realf *target_block_data = NULL;
   
   /*loop over block columns*/
   for (vmesh::LocalID block_column_i=0; block_column_i<block_column_offsets.size(); ++block_column_i) {
      const vmesh::LocalID n_cblocks = block_column_lengths[block_column_i];
      vmesh::GlobalID* cblocks = blocks + block_column_offsets[block_column_i]; /*column blocks*/
      load_column_block_data(spatial_cell,cblocks,n_cblocks,values,dimension);
      /*compute the common indices for this block column*/

      velocity_block_indices_t block_indices_begin=SpatialCell::get_velocity_block_indices(cblocks[0]); /*First block in column*/
      uint temp;
      //Switch block indices according to dimensions, the alogirthm has
      //  been written for integrating along z.
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

      for (uint j = 0; j < WID; ++j){
         /*Now it is time to compute the actual mapping*/
         /*target cell/block index contribution not dependent on k index*/
         const Vec4i target_cell_index_common = j*cell_indices_to_id[1] + Vec4i(0, cell_indices_to_id[0], 2 * cell_indices_to_id[0], 3 * cell_indices_to_id[0]);
         const int target_block_index_common(block_indices_begin[0] * block_indices_to_id[0] + block_indices_begin[1] * block_indices_to_id[1]);
         /* intersection_min is the intersection z coordinate (z after
            swaps that is) of the lowest possible z plane for each i,j
            index (i in vector)
         */
         const Real intersection_min_base = intersection +
         (block_indices_begin[0]*WID)*intersection_di +
         (block_indices_begin[1]*WID+j)*intersection_dj;
         const Vec4 intersection_min(intersection_min_base,
                     intersection_min_base + intersection_di,
                     intersection_min_base + 2.0 * intersection_di,
                     intersection_min_base + 3.0 * intersection_di);
         
         /*compute some initial values, that are used to set up the
          * shifting of values as we go through all blocks in
          * order. See comments where they are shifted for
          * explanations of their meening*/
         Vec4 v_r((WID * block_indices_begin[2]) * dv + v_min);
         Vec4i lagrangian_gk_r=truncate_to_int((v_r-intersection_min)/intersection_dk);
         
         /*loop through all blocks in column and compute the mapping as integrals*/
         for (uint k=0; k < WID * n_cblocks; ++k ){
            /*Compute reconstructions 
            values + i_pcolumnv(n_cblocks, -1, j, 0) is the starting point of the column data for fixed j
            k + WID is the index where we have stored k index, WID amount of padding
            */
#ifdef ACC_SEMILAG_PLM
            Vec4 a[2];
            compute_plm_coeff(values + i_pcolumnv(n_cblocks, -1, j, 0), k + WID , a);
#endif
#ifdef ACC_SEMILAG_PPM
            Vec4 a[3];
            compute_ppm_coeff(values + i_pcolumnv(n_cblocks, -1, j, 0), h4, k + WID, a);
#endif
#ifdef ACC_SEMILAG_PQM
            Vec4 a[5];
            compute_pqm_coeff(values + i_pcolumnv(n_cblocks, -1, j, 0), h8, k + WID, a);
#endif
            
            /*set the initial value for the integrand at the boundary at v = 0 (in reduced cell units), this will be shifted to target_density_1, see below*/
            Vec4 target_density_r(0.0);
            /*v_l, v_r are the left and right velocity coordinates of source cell. Left is the old right*/
            Vec4 v_l = v_r; 
            v_r += dv;
            /*left(l) and right(r) k values (global index) in the target
            lagrangian grid, the intersecting cells. Again old right is new left*/
            const Vec4i lagrangian_gk_l = lagrangian_gk_r;
            lagrangian_gk_r = truncate_to_int((v_r-intersection_min)/intersection_dk);
            
            Vec4i gk(lagrangian_gk_l);
            while (horizontal_or(gk <= lagrangian_gk_r)){
               const Vec4i gk_div_WID = gk/WID;
               const Vec4i gk_mod_WID = (gk - gk_div_WID * WID);
               //the block of the lagrangian cell to which we map
               const Vec4i target_block(target_block_index_common + gk_div_WID * block_indices_to_id[2]);

               //cell index in the target block 
               const Vec4i target_cell(target_cell_index_common + gk_mod_WID * cell_indices_to_id[2]);
               
               //the velocity between which we will integrate to put mass
               //in the targe cell. If both v_r and v_l are in same cell
               //then v_1,v_2 should be between v_l and v_r.
               //v_1 and v_2 normalized to be between 0 and 1 in the cell.
               //For vector elements where gk is already larger than needed (lagrangian_gk_r), v_2=v_1=v_r and thus the value is zero.
#ifdef DP
               const Vec4 v_norm_r = (min(to_double(gk + 1) * intersection_dk + intersection_min,       v_r) - v_l) * i_dv;
#else
               const Vec4 v_norm_r = (min(to_float(gk + 1) * intersection_dk + intersection_min,       v_r) - v_l) * i_dv;
#endif
               /*shift, old right is new left*/
               const Vec4 target_density_l = target_density_r;

               /*compute right integrand*/
#ifdef ACC_SEMILAG_PLM
               target_density_r =
                  v_norm_r * a[0] +
                  v_norm_r * v_norm_r * a[1];
#endif
#ifdef ACC_SEMILAG_PPM
               target_density_r =
                  v_norm_r * a[0] +
                  v_norm_r * v_norm_r * a[1] +
                  v_norm_r * v_norm_r * v_norm_r * a[2];
#endif
#ifdef ACC_SEMILAG_PQM
               target_density_r =
                  v_norm_r * a[0] +
                  v_norm_r * v_norm_r * a[1] +
                  v_norm_r * v_norm_r * v_norm_r * a[2] +
                  v_norm_r * v_norm_r * v_norm_r * v_norm_r * a[3] +
                  v_norm_r * v_norm_r * v_norm_r * v_norm_r * v_norm_r * a[4];
#endif

               /*total value of integrand*/
               const Vec4 target_density = target_density_r - target_density_l;

               //store values, one element at a time
               for (int target_i=0; target_i<4; ++target_i) {
                  const vmesh::GlobalID tblock = target_block[target_i];

                  /*check that we are within sane limits. If gk is negative,
                  * or above blocks_per_dim * blockcells_per_dim then we
                  * are outside of the target grid.*/
                  /*TODO, count losses if these are not fulfilled*/
                  if (gk[target_i] >=0 && gk[target_i] < max_v_length * WID) {
                     if (previous_target_block != tblock) {

                        // BEGIN NOTE
                        // The code inside this block is slower with the new AMR-related interface
                        previous_target_block = tblock;

                        //not the same block as last time, lets create it if we
                        //need to and fetch its data array pointer and store it in target_block_data.
                        if (spatial_cell->count(tblock) == 0) {
                           // count is faster here since the same checks in 
                           // add_velocity_block call are more expensive
                           spatial_cell->add_velocity_block(tblock);
                           phiprof_assert(spatial_cell->count(tblock) != 0);
                        }

                        target_block_data = spatial_cell->get_data( spatial_cell->get_velocity_block_local_id(tblock) );
                        // END NOTE
                     }
                     
                     // do the conversion from Real to Realf here, faster than doing it in accumulation
                     const Realf tval = target_density[target_i];
                     const uint tcell = target_cell[target_i];
                     phiprof_assert(tcell < WID3);
                     target_block_data[tcell] += tval;
                  }
               }
               gk++; //next iteration in while loop
            }
         }
      }
   }
   delete [] blocks;
   return true;
}


#endif
