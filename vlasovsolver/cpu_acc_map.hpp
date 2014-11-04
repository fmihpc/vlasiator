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
      std::cerr << "Dimension argument wrong: " << dimension << " at " << __FILE__ << ":" << __LINE__ << std::endl;
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
   // Move densities from data to fx and clear data, to prepare for mapping
   for (unsigned int cell=0; cell < VELOCITY_BLOCK_LENGTH*spatial_cell->get_number_of_velocity_blocks(); ++cell) {
      //copy data to fx for solvers, and set data to zero as we will map new values there
      spatial_cell->get_fx()[cell]   = spatial_cell->get_data()[cell];
      spatial_cell->get_data()[cell] = 0.0;
   }

bool map_1d(SpatialCell* spatial_cell,Transform<Real,3,Affine>& fwd_transform,Transform<Real,3,Affine>& bwd_transform,int dimension,int propag) {
   // Set values in an array that is used to convert cell i/j/k indices to cell id
   uint cell_indices_to_id[3];
   switch (dimension) {
    case 0:
      cell_indices_to_id[0]=WID2;
      cell_indices_to_id[1]=WID;
      cell_indices_to_id[2]=1;
      break;
    case 1:
      cell_indices_to_id[0]=1;
      cell_indices_to_id[1]=WID2;
      cell_indices_to_id[2]=WID;
      break;
    case 2:
      cell_indices_to_id[0]=1;
      cell_indices_to_id[1]=WID;
      cell_indices_to_id[2]=WID2;
      break;
   }

   // Move the old velocity mesh and data to the variables below (this is very fast)
   vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID> vmesh;
   vmesh::VelocityBlockContainer<vmesh::LocalID> blockContainer;
   spatial_cell->swap(vmesh,blockContainer);

   // Sort the blocks according to their refinement levels (very fast)
   const uint8_t maxRefLevel = vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>::getMaxAllowedRefinementLevel();
   std::vector<std::vector<vmesh::LocalID> > blocks(maxRefLevel+1);
   std::vector<double> f_sums(maxRefLevel+1);
   for (int i=0; i<maxRefLevel; ++i) f_sums[i] = 0.0;

   for (vmesh::LocalID block=0; block<vmesh.size(); ++block) {
      uint8_t refLevel = vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>::getRefinementLevel(vmesh.getGlobalID(block));
      blocks[refLevel].push_back(block);
   }

   // ***** DEFINITIONS ***** //
   // The velocity mesh (at given refinement level) is a Cartesian mesh 
   // consisting of (Nx*WID) * (Ny*WID) * (Nz*WID) cells. The mesh is 
   // grouped into blocks, each consisting of WID*WID*WID cells. There are 
   // Nx*Ny*Nz blocks in total. Only some of the blocks actually exist.
   // 
   // The simulation domain is covered by several Cartesian meshes, one mesh 
   // for each refinement level. A mesh at a higher refinement level has twice 
   // as many cells per coordinate.
   // 
   // The velocity mesh can be indexed in multiple ways:
   // 1) Calculate the i/j/k indices of cells in the "global" (Nx*WID)*(Ny*WID)*(Nz*WID) mesh.
   //    IT IS NOT SAFE TO DO THIS IN PRACTICE! The global cell index will not fit in 64bit unsigned int.
   //    However, it is safe to calculate the i/j/k cell indices separately, and these 
   //    are called "cell global indices". Note that the i/j/k cell indices depend on the refinement level.
   // 
   // 2) 
   
   // PAD is the number of values copied from face neighbors (4 = all neighbor values)
   const int PAD=4;
   Realf* f_source_array = new Realf[WID*WID*(WID+2*PAD)];
   Vec4 f_tmp[WID*(WID+2*PAD)];

   // Accelerate all blocks, starting from the highest refinement level
   for (uint8_t refLevel=blocks.size()-1; refLevel<blocks.size(); --refLevel) {
      phiprof::start("compute-intersections");
      // Compute intersections at this refinement level:
      Real intersection,intersection_di,intersection_dj,intersection_dk;
      switch (propag) {
       case 0:
<<<<<<< HEAD
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
=======
	 compute_intersections_1st(spatial_cell, bwd_transform, fwd_transform, dimension, refLevel,
				   intersection,intersection_di,intersection_dj,intersection_dk);
	 break;
       case 1:
	 compute_intersections_2nd(spatial_cell, bwd_transform, fwd_transform, dimension, refLevel,
				   intersection,intersection_di,intersection_dj,intersection_dk);
	 break;
       case 2:
	 compute_intersections_3rd(spatial_cell, bwd_transform, fwd_transform, dimension, refLevel,
				   intersection,intersection_di,intersection_dj,intersection_dk);
	 break;
       default:
	 std::cerr << "error in map1d" << std::endl; exit(1);
	 break;
      }
      phiprof::stop("compute-intersections");

      phiprof::start("compute-mapping");

      // Set values in an array that is used to convert block i/j/k indices to global ID
      Real dv,v_min,inv_dv;
      uint block_indices_to_id[3];
      vmesh::GlobalID k_cell_global_target_max;
      switch (dimension) {
       case 0: {
	  dv    = SpatialCell::get_velocity_grid_cell_size(refLevel)[0];
	  v_min = SpatialCell::get_velocity_grid_min_limits()[0];
	  inv_dv = 1.0/dv;
	  k_cell_global_target_max = SpatialCell::get_velocity_grid_length(refLevel)[0]*WID;
	  
	  block_indices_to_id[0]=SpatialCell::get_velocity_grid_length(refLevel)[0]*SpatialCell::get_velocity_grid_length(refLevel)[1];
	  block_indices_to_id[1]=SpatialCell::get_velocity_grid_length(refLevel)[0];
	  block_indices_to_id[2]=1;
	  
	  // Swap intersection di and dk
	  const Real tmp=intersection_di;
	  intersection_di=intersection_dk;
	  intersection_dk=tmp;}
	 break;
       case 1: {
	  dv    = SpatialCell::get_velocity_grid_cell_size(refLevel)[1];
	  v_min = SpatialCell::get_velocity_grid_min_limits()[1];
	  inv_dv = 1.0/dv;
	  k_cell_global_target_max = SpatialCell::get_velocity_grid_length(refLevel)[1]*WID;
	  
	  block_indices_to_id[0]=1;
	  block_indices_to_id[1]=SpatialCell::get_velocity_grid_length(refLevel)[0]*SpatialCell::get_velocity_grid_length(refLevel)[1];
	  block_indices_to_id[2]=SpatialCell::get_velocity_grid_length(refLevel)[0];
	  
	  // Swap intersection dj and dk
	  const Real tmp=intersection_dj;
	  intersection_dj=intersection_dk;
	  intersection_dk=tmp;}
	 break;
       case 2:
	 dv    = SpatialCell::get_velocity_grid_cell_size(refLevel)[2];
	 v_min = SpatialCell::get_velocity_grid_min_limits()[2];
	 inv_dv = 1.0/dv;
	 k_cell_global_target_max = SpatialCell::get_velocity_grid_length(refLevel)[2]*WID;
	 
	 block_indices_to_id[0]=1;
	 block_indices_to_id[1]=SpatialCell::get_velocity_grid_length(refLevel)[0];
	 block_indices_to_id[2]=SpatialCell::get_velocity_grid_length(refLevel)[0]*SpatialCell::get_velocity_grid_length(refLevel)[1];
	 break;
       default:
	 std::cerr << "Incorrect dimension in " << __FILE__ << ' ' << __LINE__ << std::endl; exit(1);
	 break;
      }

      // Accelerate all blocks at this refinement level, blocks are in random order
      for (size_t b=0; b<blocks[refLevel].size(); ++b) {
	 const vmesh::LocalID blockLID = blocks[refLevel][b];
	 const vmesh::GlobalID blockGID = vmesh.getGlobalID(blockLID);

	 vmesh::GlobalID prevTargetBlockGID = vmesh.invalidGlobalID();
	 Realf* targetBlockData = NULL;

	 // Fetch pseudo-density of the propagated data.
	 // Note: although we call SpatialCell function here, it uses the velocity mesh 
	 // vmesh we pass as a parameter instead of the velocity mesh currently in the cell
	 phiprof::start("fetch-data");
	 spatial_cell->fetch_acc_data<PAD>(blockGID,dimension,vmesh,blockContainer.getData(),f_source_array);
	 phiprof::stop("fetch-data");

	 #warning FIXME Data transpose can be merged to fetch_acc_data function
	 phiprof::start("dummy-transpose");
	 // TEMP (swap j and k indices, load to vector)
	 for (int k=0; k<WID+2*PAD; ++k) for (int j=0; j<WID; ++j) for (int i=0; i<WID; ++i) {
	    const int tgtIndex = j*(WID+2*PAD) + k;
	    f_tmp[tgtIndex].insert(i,0);
	 }

	 for (int k=0; k<WID+2*PAD; ++k) for (int j=0; j<WID; ++j) for (int i=0; i<WID; ++i) {
	    const int srcIndex = k*WID*WID + j*WID + i;
	    const int tgtIndex = j*(WID+2*PAD) + k;
	    f_tmp[tgtIndex].insert(i,f_source_array[srcIndex]);
	 }
	 // END TEMP
	 phiprof::stop("dummy-transpose");

	 // block_indices_begin are just the (i,j,k) indices of the block
	 // (at correct refinement level) we're working with:
	 velocity_block_indices_t block_indices_begin = SpatialCell::get_velocity_block_indices(blockGID);

	 // Switch block indices according to dimensions, the alogirthm has
	 // been written for integrating along z
	 switch (dimension) {
	  case 0: {
	     // i and k coordinates have been swapped
	     const uint temp=block_indices_begin[2];
	     block_indices_begin[2]=block_indices_begin[0];
	     block_indices_begin[0]=temp; }
	    break;
	  case 1: {
	     // in values j and k coordinates have been swapped
	     const uint temp=block_indices_begin[2];
	     block_indices_begin[2]=block_indices_begin[1];
	     block_indices_begin[1]=temp; }
	    break;
	  case 2:
	    break;
	 }

	 for (int j=0; j<WID; ++j) {
	    // Now it is time to compute the actual mapping

	    // target cell/block index contribution not dependent on k index

	    // Below target_cell_index_common = j*WID + i for z-translation, and i = 0/1/2/3 (vectorized).
	    // Multiplying by cell_indices_to_id takes transposing into account:
	    const Vec4i target_cell_index_common = j*cell_indices_to_id[1] 
	                                         + Vec4i(0,cell_indices_to_id[0],2*cell_indices_to_id[0],3*cell_indices_to_id[0]);

	    // Below target_block_index_common = j_block*Nx + i_block is 
	    // calculated at current refinement level. Multiplying by block_indices_to_id
	    // takes transposing into account:
	    const vmesh::GlobalID target_block_index_common = block_indices_begin[0]*block_indices_to_id[0] + block_indices_begin[1]*block_indices_to_id[1]
	    + vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>::getGlobalIndexOffset(refLevel);
	    
	    const Real intersection_min_base = intersection +
	      (block_indices_begin[0]*WID)*intersection_di +
	      (block_indices_begin[1]*WID+j)*intersection_dj;

	    // intersection_min is the intersection z coordinate (taking transpose
	    // into account) of the lowest possible z plane for each cell i/j index
	    const Vec4 intersection_min(intersection_min_base,
					intersection_min_base + intersection_di,
					intersection_min_base + 2.0 * intersection_di,
					intersection_min_base + 3.0 * intersection_di);

	    // Compute some initial values that are used to set up the
	    // shifting of values as we go through all blocks in
	    // order. See the comments to see where the the values are shifted to.

	    // Calculate the velocity z-coordinate corresponding to the bottom z-face of the source block
	    Vec4 v_top((WID * block_indices_begin[2]) * dv + v_min);
	    Vec4i lagrangian_gk_top = truncate_to_int((v_top-intersection_min)/intersection_dk);

	    for (uint k=0; k<WID; ++k) {
	       phiprof::start("reconstruction");
	       // Compute reconstructions
	       #warning FIXME reconstructions are incorrect for AMR, cell size changes not taken into account!
	       #ifdef ACC_SEMILAG_PLM
	          Vec4 a[2];
	          compute_plm_coeff(f_tmp + j*(WID+2*PAD),k+PAD,a);
	       #endif
	       #ifdef ACC_SEMILAG_PPM
	          Vec4 a[3];
	          compute_ppm_coeff(f_tmp + j*(WID+2*PAD),h8,k+PAD,a);
	       #endif
	       #ifdef ACC_SEMILAG_PQM
	          Vec4 a[5];
	          compute_pqm_coeff(f_tmp + j*(WID+2*PAD),h8,k+PAD,a);
               #endif
	       phiprof::stop("reconstruction");

	       // Initialize cumulative pseudo-density to zero value. Cumulative 
	       // pseudo-density is the density integrated from an arbitrary origin
	       // up to the current velocity z-coordinate
	       Vec4 target_density_top(0.0);

	       // v_bot and v_top are the velocity z-coordinates of the bottom 
	       // and top faces of the source cells
	       Vec4 v_bot = v_top;
	       v_top += dv;

	       // The left(l) and right(r) k-indices (global index) in the target
	       // lagrangian grid, the intersecting cells.
	       const Vec4i lagrangian_gk_bot = lagrangian_gk_top;
	       lagrangian_gk_top = truncate_to_int((v_top-intersection_min)/intersection_dk);

	       phiprof::start("deposition");
	       // Global k-indices of the propagated cells in the target
	       // (Lagrangian) mesh after mapping. 
	       Vec4i k_cell_global_target(lagrangian_gk_bot); // old gk
	       while (horizontal_or(k_cell_global_target <= lagrangian_gk_top)) {
		  // k-indices of the blocks that own target cells
		  const Vec4i k_block_target = k_cell_global_target/WID;
	       
		  // k-indices of the target cells in the target blocks
		  const Vec4i k_cell_target = (k_cell_global_target - k_block_target*WID);

		  // Global ID of the target block
		  const Vec4i target_block(target_block_index_common + k_block_target*block_indices_to_id[2]);

		  // Global ID of the cell in the target (Lagrangian) block, 
		  // i.e., cells' 3D indices in WID*WID*WID size array
		  const Vec4i target_cell(target_cell_index_common + k_cell_target*cell_indices_to_id[2]);

		  // Calculate the integration limit in velocity z-direction. 
		  // The phase-space density is integrated from the bottom cell 
		  // z-face to the z-coordinate calculated below. The coordinates 
		  // are normalized by cell widths, i.e., values are between 0 and 1.
		  // 
		  // If k_cell_global_target is larger than what needed then the
		  // integration limits should be equal, and zero value is stored 
		  // to target cell.
		  #ifdef DP
		     const Vec4 v_norm_top = (min(to_double(k_cell_global_target+1)*intersection_dk + intersection_min,v_top) - v_bot)*inv_dv;
		  #else
		     const Vec4 v_norm_top = (min(to_float(k_cell_global_target+1)*intersection_dk + intersection_min, v_top) - v_bot)*inv_dv;
		  #endif

		  // Shift pseudo-densities, new bottom value will be the old top value
		  const Vec4 target_density_bot = target_density_top;
		  
		  // Calculate cumulative pseudo-density at vz=v_norm_top, i.e., 
		  // integrate pseudo-density from bottom vz-face to vz=v_norm_top
		  #warning FIXME Integrations may not be correct for AMR mesh!
		  #ifdef ACC_SEMILAG_PLM
		  target_density_top =
		    v_norm_top * a[0] +
		    v_norm_top * v_norm_top * a[1];
		  #endif
		  #ifdef ACC_SEMILAG_PPM
		  target_density_top =
		    v_norm_top * a[0] +
		    v_norm_top * v_norm_top * a[1] +
		    v_norm_top * v_norm_top * v_norm_top * a[2];
		  #endif
		  #ifdef ACC_SEMILAG_PQM
		  target_density_top =
		    v_norm_top * a[0] +
		    v_norm_top * v_norm_top * a[1] +
		    v_norm_top * v_norm_top * v_norm_top * a[2] +
		    v_norm_top * v_norm_top * v_norm_top * v_norm_top * a[3] +
		    v_norm_top * v_norm_top * v_norm_top * v_norm_top * v_norm_top * a[4];
		  #endif

		  // Pseudo-density deposited to target cell is the difference of the 
		  // cumulative pseudo-densities:
		  const Vec4 target_density = target_density_top - target_density_bot;

		  // There doesn't seem to be much difference between this version and the one below.
		  // It might be possible to optimize this further by only creating the target block, 
		  // not the whole octant. The blocks missing from the octant have to be created 
		  // during the data merge then.
		  const vmesh::GlobalID targetGID = target_block[0];

		  // Add the target block of each source cell
		  for (uint i=0; i<WID; ++i) {
		     // Check that target k-index is inside the target grid:
		     if (k_cell_global_target[i] < 0 || k_cell_global_target[i] >= k_cell_global_target_max) {
			continue;
		     }

		     const vmesh::GlobalID targetGID = target_block[i];
		     if (targetGID != prevTargetBlockGID) {
			spatial_cell->add_velocity_block_octant(targetGID);
			vmesh::LocalID targetLID = spatial_cell->get_velocity_block_local_id(targetGID);
			targetBlockData = spatial_cell->get_data(targetLID);
			prevTargetBlockGID = targetGID;
		     }

		     // Store values...
		     const Realf f_target = target_density[i];
		     const int index_cell_target = target_cell[i];
		     targetBlockData[index_cell_target] += f_target;
		  }

		  ++k_cell_global_target;
	       } // while (horizontal_or(k_cell_global_target <= lagrangian_gk_r))
	       phiprof::stop("deposition");
	    }
	 }
      }
      phiprof::stop("compute-mapping");
>>>>>>> Vlasiator now accelerates with AMR mesh. Reconstructions are still not
   }
   delete [] f_source_array; f_source_array = NULL;

   // Merge values from coarse blocks to refined blocks wherever the same domain 
   // is covered by overlapping blocks (at different refinement levels)
   phiprof::start("merge values");
   spatial_cell->merge_values();
   phiprof::stop("merge values");
   
   /*if (spatial_cell->checkMesh() == false) {
      //std::cerr << "error(s) in mesh, exiting" << std::endl; exit(1);
   }*/
   return true;
}

<<<<<<< HEAD

#endif
=======
#endif   
>>>>>>> Vlasiator now accelerates with AMR mesh. Reconstructions are still not
