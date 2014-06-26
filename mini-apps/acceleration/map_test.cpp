#include "vlasovsolver/vec4.h"

#include "vlasovsolver/cpu_1d_column_interpolations.hpp"
#include "definitions.h"
#include "common.h"

#define MAX_BLOCKS_PER_DIM 100
#ifdef ACC_SEMILAG_PLM
#define RECONSTRUCTION_ORDER 1
#endif
#ifdef ACC_SEMILAG_PPM
#define RECONSTRUCTION_ORDER 2
#endif
#ifdef ACC_SEMILAG_PQM
#define RECONSTRUCTION_ORDER 4
#endif


int main(void) {
  const int dv=20000;
  const Real min_v=-4e6;
  const int blocks_per_dim=100;
  const int j_block = 0; //in a real case it would be something else!
  const int i_block = 0; //in a real case it would be something else!
  const int j = 0;

  Vec4 values[(blocks_per_dim+2)*WID];
  Vec4 target[(blocks_per_dim+2)*WID];
  Vec4 a[blocks_per_dim*WID][RECONSTRUCTION_ORDER + 1];  
  
  /*init values*/
  
  Real intersection=min_v;
  Real intersection_di=dv/4.0;;
  const int iterations=1;
  
  /*square wave*/
  for(i=0;i<n_blocks_per_dim*WID;i++){
    Real v=min_v + i*dv;
    if(v > min_v*0.25 && v < 0.25*(min_v + blocks_per_dim*WID*dv))
      values[i]=1.0;
    else
      values[i]=9.0;
  }
  
  /*loop over propagations*/
  
  for(step = 0; step < iterations; step++){
#ifdef ACC_SEMILAG_PLM
         compute_plm_coeff_explicit_column(values, blocks_per_dim, a);
#endif
#ifdef ACC_SEMILAG_PPM
         compute_ppm_coeff_explicit_column(values, blocks_per_dim, a);
#endif
#ifdef ACC_SEMILAG_PQM
         compute_pqm_coeff_explicit_column(values, blocks_per_dim, a);
#endif
         /* intersection_min is the intersection z coordinate (z after
            swaps that is) of the lowest possible z plane for each i,j
            index (i in vector)
         */	 
         const Real intersection_min_base =  intersection +
	   (block_i * WID) * intersection_di + 
	   (block_j * WID + j) * intersection_dj;
         const Vec4 intersection_min(intersection_min_base + 0 * intersection_di,
                                     intersection_min_base + 1 * intersection_di,
                                     intersection_min_base + 2 * intersection_di,
                                     intersection_min_base + 3 * intersection_di);
         

         /*compute some initial values, that are used to set up the
          * shifting of values as we go through all blocks in
          * order. See comments where they are shifted for
          * explanations of their meening*/
         Vec4 v_r(v_min);
         Vec4i lagrangian_gk_r=truncate_to_int(( v_r - intersection_min)/intersection_dk);
	 
         /*loop through all blocks in column and compute the mapping as integrals*/
         for (unsigned int block_i = 0; block_i<blocks_per_dim;block_i++){
	   for (uint k=0; k<WID; ++k){ 
	     /*set the initial value for the integrand at the boundary at v = 0 (in reduced cell units), this will be shifted to target_density_1, see below*/
	     Vec4 target_density_r(0.0);
	     /*v_l, v_r are the left and right velocity coordinates of source cell. Left is the old right*/
	     Vec4 v_l = v_r; 
	     v_r += dv;
	     /*left(l) and right(r) k values (global index) in the target
	       lagrangian grid, the intersecting cells. Again old right is new left*/               
	     const Vec4i lagrangian_gk_l = lagrangian_gk_r;
	     lagrangian_gk_r = truncate_to_int((v_r - intersection_min)/intersection_dk);
	     
	     Vec4i gk(lagrangian_gk_l);	
	     while (horizontal_or(gk <= lagrangian_gk_r)){
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
		 v_norm_r * a[block_i * WID + k][0] +
		 v_norm_r * v_norm_r * a[block_i * WID + k][1];
#endif
#ifdef ACC_SEMILAG_PPM
	       target_density_r =
		 v_norm_r * a[block_i * WID + k][0] +
		 v_norm_r * v_norm_r * a[block_i * WID + k][1] +
		 v_norm_r * v_norm_r * v_norm_r * a[block_i * WID + k][2];
#endif

	       /*total value of integrand*/
	       const Vec4 target_density = target_density_r - target_density_l;
	       
	       //store values, one element at elema time
	       for(uint elem = 0; elem < 4;elem ++ ){
		 /*TODO, count losses if these are not fulfilled*/
		 int k_in_target = gk[elem];
		 if (k_in_target >=0 &&
		     k_in_target < blocks_per_dim * WID) {
		   target[k_in_target].insert(elem, target_cell[elem] + target[k_in_target]);
		 }
	       }		   
	       gk++; //next iteration in while loop
	     }
	   }
         }
  }
}
