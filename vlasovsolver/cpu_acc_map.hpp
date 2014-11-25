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

using namespace spatial_cell;

struct PropagParams {
   Real intersection;
   Real intersection_di;
   Real intersection_dj;
   Real intersection_dk;
   
<<<<<<< HEAD
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
=======
   Real dv;
   Real v_min;
   Real inv_dv;

   int dimension;
   int i_mapped;
   int j_mapped;
   int k_mapped;
   int refMul;
   uint8_t refLevel;
   uint8_t maxRefLevel;
   vmesh::GlobalID Nx;
   vmesh::GlobalID Ny;
   vmesh::GlobalID k_cell_global_target_max;
};

void map_1d(SpatialCell* spatial_cell,PropagParams& params,
	    vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& vmesh,
	    vmesh::VelocityBlockContainer<vmesh::LocalID>& blockContainer);

void generateTargetMesh(SpatialCell* spatial_cell,const std::vector<vmesh::LocalID>& blocks,PropagParams& params,
			const uint8_t& targetRefLevel,const vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& vmesh);
>>>>>>> aaba9b5ef3aa49cd57a3c4bce747f441805dc757

/* 
   Here we map from the current time step grid, to a target grid which
   is the lagrangian departure grid (so th grid at timestep +dt,
   tracked backwards by -dt)

   TODO: parallelize with openMP over block-columns. If one also
   pre-creates new blocks in a separate loop first (serial operation),
   then the openmp parallization would scale well (better than over
   spatial cells), and would not need synchronization.
   
*/

<<<<<<< HEAD
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
=======
bool map_1d(SpatialCell* spatial_cell,Transform<Real,3,Affine>& fwd_transform,Transform<Real,3,Affine>& bwd_transform,int dimension,int propag) {
   // Move the old velocity mesh and data to the variables below (very fast)
   vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID> vmesh;
   vmesh::VelocityBlockContainer<vmesh::LocalID> blockContainer;
   spatial_cell->swap(vmesh,blockContainer);

   // Sort the blocks according to their refinement levels (very fast)
   const uint8_t maxRefLevel = vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>::getMaxAllowedRefinementLevel();
   std::vector<std::vector<vmesh::LocalID> > blocks(maxRefLevel+1);
   for (vmesh::LocalID block=0; block<vmesh.size(); ++block) {
      uint8_t refLevel = vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>::getRefinementLevel(vmesh.getGlobalID(block));
      blocks[refLevel].push_back(block);
>>>>>>> aaba9b5ef3aa49cd57a3c4bce747f441805dc757
   }

   // Computer intersections etc.
   PropagParams propagParams;
   propagParams.maxRefLevel = vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>::getMaxAllowedRefinementLevel();
   propagParams.Nx = SpatialCell::get_velocity_grid_length(propagParams.maxRefLevel)[0];
   propagParams.Ny = SpatialCell::get_velocity_grid_length(propagParams.maxRefLevel)[1];
   propagParams.dimension = dimension;
   switch (propag) {
    case 0:
      compute_intersections_1st(spatial_cell, bwd_transform, fwd_transform, dimension, propagParams.maxRefLevel,
				propagParams.intersection,propagParams.intersection_di,propagParams.intersection_dj,propagParams.intersection_dk);
      break;
    case 1:
      compute_intersections_2nd(spatial_cell, bwd_transform, fwd_transform, dimension, propagParams.maxRefLevel,
				propagParams.intersection,propagParams.intersection_di,propagParams.intersection_dj,propagParams.intersection_dk);
      break;
    case 2:
      compute_intersections_3rd(spatial_cell, bwd_transform, fwd_transform, dimension, propagParams.maxRefLevel,
				propagParams.intersection,propagParams.intersection_di,propagParams.intersection_dj,propagParams.intersection_dk);
      break;
    default:
      std::cerr << "error in map1d" << std::endl; exit(1);
      break;
   }

   propagParams.dimension = dimension;
   propagParams.dv    = SpatialCell::get_velocity_grid_cell_size(propagParams.maxRefLevel)[dimension];
   propagParams.v_min = SpatialCell::get_velocity_grid_min_limits()[dimension];
   propagParams.inv_dv = 1.0/propagParams.dv;
   propagParams.k_cell_global_target_max = SpatialCell::get_velocity_grid_length(propagParams.maxRefLevel)[dimension]*WID;
   switch (dimension) {
    case 0: {
      propagParams.i_mapped = 2;
      propagParams.j_mapped = 1;
      propagParams.k_mapped = 0;
      
      // Swap intersection di and dk
      const Real tmp=propagParams.intersection_di;
      propagParams.intersection_di=propagParams.intersection_dk;
      propagParams.intersection_dk=tmp;}
      break;
    case 1: {
      propagParams.i_mapped = 0;
      propagParams.j_mapped = 2;
      propagParams.k_mapped = 1;
      
      // Swap intersection dj and dk
      const Real tmp=propagParams.intersection_dj;
      propagParams.intersection_dj=propagParams.intersection_dk;
      propagParams.intersection_dk=tmp;}
      break;
    case 2:
      propagParams.i_mapped = 0;
      propagParams.j_mapped = 1;
      propagParams.k_mapped = 2;
      break;
    default:
      std::cerr << "error in map1d" << std::endl; exit(1);
      break;
   }
      
   bool rvalue = true;

   phiprof::start("mesh generation");
   for (uint8_t r=0; r<blocks.size(); ++r) {
      for (uint8_t rr=r; rr<blocks.size(); ++rr) {
	 propagParams.refLevel = rr;
	 generateTargetMesh(spatial_cell,blocks[rr],propagParams,r,vmesh);
      }
   }
   phiprof::stop("mesh generation");

   phiprof::start("mapping");
   map_1d(spatial_cell,propagParams,vmesh,blockContainer);
   phiprof::stop("mapping");

   // Merge values from coarse blocks to refined blocks wherever the same domain 
   // is covered by overlapping blocks (at different refinement levels)
   //spatial_cell->merge_values();
   
<<<<<<< HEAD
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
=======
   //if (spatial_cell->checkMesh() == false) {
      //std::cerr << "error(s) in mesh, exiting" << std::endl; exit(1);
   //}

   rvalue = false;
}

void generateTargetMesh(SpatialCell* spatial_cell,const std::vector<vmesh::LocalID>& blocks,PropagParams& params,
			const uint8_t& targetRefLevel,const vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& vmesh) {
   params.refMul = pow(2,(params.maxRefLevel-params.refLevel));
   const int baseMul = pow(2,params.maxRefLevel-targetRefLevel);

   for (size_t b=0; b<blocks.size(); ++b) {
      const vmesh::LocalID blockLID = blocks[b];
      const vmesh::GlobalID blockGID = vmesh.getGlobalID(blockLID);
      vmesh::LocalID sourceIndex[3];
      vmesh.getIndices(blockGID,params.refLevel,sourceIndex[0],sourceIndex[1],sourceIndex[2]);
      
      switch (params.dimension) {
       case 0: {
	  const vmesh::LocalID temp=sourceIndex[2];
	  sourceIndex[2]=sourceIndex[0];
	  sourceIndex[0]=temp;}
	 break;
       case 1: {
	  vmesh::LocalID temp=sourceIndex[2];
	  sourceIndex[2]=sourceIndex[1];
	  sourceIndex[1]=temp;}
	 break;
>>>>>>> aaba9b5ef3aa49cd57a3c4bce747f441805dc757
       case 2:
         break;
      }
      sourceIndex[0] *= WID*params.refMul;
      sourceIndex[1] *= WID*params.refMul;
      sourceIndex[2] *= WID*params.refMul;

      vmesh::GlobalID k_trgt_min = std::numeric_limits<vmesh::GlobalID>::max();
      vmesh::GlobalID k_trgt_max = 0;
      for (int k=0; k<2; ++k) for (int j=0; j<2; ++j) for (int i=0; i<2; ++i) {
	 const int i_cell = i*WID*params.refMul;
	 const int j_cell = j*WID*params.refMul;
	 const int k_cell = k*WID*params.refMul;

	 Realf intersection_min_base = params.intersection +
	   sourceIndex[0]*params.intersection_di +
	   (sourceIndex[1]+j_cell)*params.intersection_dj;

	 Realf intersection_min = intersection_min_base + i_cell*params.intersection_di;
	 Realf v_top = params.v_min + (sourceIndex[2]+k_cell)*params.dv;
	 vmesh::GlobalID lagrangian_gk_top = static_cast<vmesh::GlobalID>((v_top-intersection_min)/params.intersection_dk);

	 if (lagrangian_gk_top < k_trgt_min) k_trgt_min = lagrangian_gk_top;
	 if (lagrangian_gk_top > k_trgt_max) k_trgt_max = lagrangian_gk_top;
      }

<<<<<<< HEAD
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
=======
      vmesh::GlobalID targetBlockIndex[3];
      //targetBlockIndex[params.i_mapped] = sourceIndex[params.i_mapped] / (WID*baseMul);
      //targetBlockIndex[params.j_mapped] = sourceIndex[params.j_mapped] / (WID*baseMul);
      targetBlockIndex[params.i_mapped] = sourceIndex[0] / (WID*baseMul);
      targetBlockIndex[params.j_mapped] = sourceIndex[1] / (WID*baseMul);
      k_trgt_min /= (WID*baseMul);
      k_trgt_max /= (WID*baseMul);
      for (vmesh::GlobalID k=k_trgt_min; k<=k_trgt_max; ++k) {
	 targetBlockIndex[params.k_mapped] = k;
	 vmesh::GlobalID targetBlock = vmesh.getGlobalID(targetRefLevel,targetBlockIndex);

	 if (targetRefLevel == 0) {
	    spatial_cell->add_velocity_block(targetBlock);
	 } else {
	    targetBlock = vmesh.getParent(targetBlock);
	    std::map<vmesh::GlobalID,vmesh::LocalID> insertedBlocks;
	    spatial_cell->refine_block(targetBlock,insertedBlocks);
	 }
>>>>>>> aaba9b5ef3aa49cd57a3c4bce747f441805dc757
      }
   }
}

/**
 * @param spatial_cell Propagated spatial cell, contains a valid (unsorted) target mesh.
 * @param params Parameters needed in the propagation.
 * @param vmesh The source mesh.
 * @param blockContainer Source mesh data.*/
void map_1d(SpatialCell* spatial_cell,PropagParams& params,
	   vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& vmesh,
	   vmesh::VelocityBlockContainer<vmesh::LocalID>& blockContainer) {

   std::vector<vmesh::GlobalID> removeList;
   
   for (vmesh::LocalID targetLID=0; targetLID<spatial_cell->get_number_of_velocity_blocks(); ++targetLID) {
      vmesh::GlobalID targetGID = spatial_cell->get_velocity_block_global_id(targetLID);
      params.refLevel = vmesh.getRefinementLevel(targetGID);
      params.refMul = std::pow(2,params.maxRefLevel-params.refLevel);
      vmesh::LocalID targetIndex[3];
      vmesh.getIndices(targetGID,params.refLevel,targetIndex[0],targetIndex[1],targetIndex[2]);

      // Pointer to target data
      Realf* data_trgt = spatial_cell->get_data(targetLID);
      
      switch (params.dimension) {
       case 0: {
	  const vmesh::LocalID temp=targetIndex[2];
	  targetIndex[2]=targetIndex[0];
	  targetIndex[0]=temp;
       }
	 break;
       case 1: {
	  vmesh::LocalID temp=targetIndex[2];
	  targetIndex[2]=targetIndex[1];
	  targetIndex[1]=temp;
       }
	 break;
       case 2:
	 break;
      }
      targetIndex[0] *= WID*params.refMul;
      targetIndex[1] *= WID*params.refMul;
      targetIndex[2] *= WID*params.refMul;

      Real accum=0;

      // Note: the i,j,k indices below are transposed indices:
      for (int k=0; k<WID; ++k) for (int j=0; j<WID; ++j) for (int i=0; i<WID; ++i) {
	 // Target cell global i/j/k indices (transposed):
	 vmesh::GlobalID i_cell_trgt = targetIndex[0] + i*params.refMul;
	 vmesh::GlobalID j_cell_trgt = targetIndex[1] + j*params.refMul;
	 vmesh::GlobalID k_cell_trgt = targetIndex[2] + k*params.refMul;

	 Realf intersection_min = params.intersection + i_cell_trgt*params.intersection_di + j_cell_trgt*params.intersection_dj;

	 // Velocities where the bottom and top z-faces of the target cell maps to:
	 Realf v_src_bot = k_cell_trgt*params.intersection_dk + intersection_min;
	 Realf v_src_top = (k_cell_trgt+params.refMul)*params.intersection_dk + intersection_min;

	 Realf v_src_bot_l = (v_src_bot-params.v_min)/params.dv;
	 Realf v_src_top_l = (v_src_top-params.v_min)/params.dv;

	 // Source cell global k-index:
	 vmesh::GlobalID k_cell_src_bot = static_cast<vmesh::GlobalID>((v_src_bot-params.v_min)/params.dv);
	 vmesh::GlobalID k_cell_src_top = static_cast<vmesh::GlobalID>((v_src_top-params.v_min)/params.dv);
	 if (k_cell_src_top < k_cell_src_bot) {
	    vmesh::GlobalID tmp=k_cell_src_bot;
	    k_cell_src_bot = k_cell_src_top;
	    k_cell_src_top = tmp;
	 }

	 // Source cell global i/j/k indices (untransposed):
	 vmesh::GlobalID srcIndex[3];
	 srcIndex[params.i_mapped] = i_cell_trgt;
	 srcIndex[params.j_mapped] = j_cell_trgt;

	 vmesh::GlobalID k_cell_src = k_cell_src_bot;
	 while (k_cell_src <= k_cell_src_top) {
	    // Find the source block
	    srcIndex[params.k_mapped] = k_cell_src;
	    uint8_t srcRefLevel = params.refLevel;
	    vmesh::GlobalID sourceGID = vmesh.findBlockDown(srcRefLevel,srcIndex);
	    int srcRefMul = std::pow(2,params.maxRefLevel-srcRefLevel);
	    
	    if (sourceGID == vmesh.invalidGlobalID()) {
	       ++k_cell_src; 
	       continue;
	    }
	    
	    Realf* data_src = blockContainer.getData(vmesh.getLocalID(sourceGID));

	    vmesh::GlobalID i_block_src = srcIndex[params.i_mapped] / (WID*srcRefMul);
	    vmesh::GlobalID j_block_src = srcIndex[params.j_mapped] / (WID*srcRefMul);
	    vmesh::GlobalID k_block_src = k_cell_src/(WID*srcRefMul);
	    int k_cell_bot = (k_cell_src - k_block_src*(WID*srcRefMul)) / srcRefMul;

	    int srcCellIndex[3];
	    srcCellIndex[params.i_mapped] = (i_cell_trgt - i_block_src*(WID*srcRefMul)) / srcRefMul;
	    srcCellIndex[params.j_mapped] = (j_cell_trgt - j_block_src*(WID*srcRefMul)) / srcRefMul;
	    srcCellIndex[params.k_mapped] = k_cell_bot;

	    Realf v_bot = std::max(v_src_bot_l,(Realf)(k_block_src*(WID*srcRefMul)+(k_cell_bot  )*srcRefMul));
	    Realf v_top = std::min(v_src_top_l,(Realf)(k_block_src*(WID*srcRefMul)+(k_cell_bot+1)*srcRefMul));

	    Realf v_int_min = (v_bot - (k_block_src*(WID*srcRefMul)+k_cell_bot*srcRefMul)) / srcRefMul;
	    Realf v_int_max = (v_top - (k_block_src*(WID*srcRefMul)+k_cell_bot*srcRefMul)) / srcRefMul;

	    Realf factor = std::pow(2,params.refLevel-srcRefLevel);
	    
	    int trgtCellIndex[3];
	    trgtCellIndex[params.i_mapped] = i;
	    trgtCellIndex[params.j_mapped] = j;
	    trgtCellIndex[params.k_mapped] = k;
	    const int trgtCell = vblock::index(trgtCellIndex[0],trgtCellIndex[1],trgtCellIndex[2]);
	    
	    data_trgt[trgtCell] += data_src[vblock::index(srcCellIndex[0],srcCellIndex[1],srcCellIndex[2])] * (v_int_max-v_int_min) * factor;
	    accum += data_src[vblock::index(srcCellIndex[0],srcCellIndex[1],srcCellIndex[2])] * (v_int_max-v_int_min) * factor;
	    k_cell_src = k_block_src*(WID*srcRefMul)+(k_cell_bot+1)*srcRefMul;
	 }
      }

      if (accum < 1e-30) {
	 removeList.push_back(targetGID);
      }
   }

   for (size_t b=0; b<removeList.size(); ++b) spatial_cell->remove_velocity_block(removeList[b]);
}

#endif
