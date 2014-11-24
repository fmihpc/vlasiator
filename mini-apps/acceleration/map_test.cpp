#include <stdio.h>
#include "common.h"
#include "cpu_1d_column_interpolations.hpp"

/*print all values in the vector valued values array. In this array
  there are blocks_per_dim blocks with a width of WID*/
void print_values(int step, Real *values, uint blocks_per_dim, Real v_min, Real dv){
  char name[256];
  sprintf(name,"dist_%03d.dat",step);

  FILE* fp=fopen(name,"w");
  for(int i=0; i < blocks_per_dim * WID; i++){
    Real v = v_min + i*dv;
    fprintf(fp,"%20.12g %20.12g\n", v, values[i + WID]);
  }
  fclose(fp);
}


void propagate(Real values[], uint  blocks_per_dim, Real v_min, Real dv,
               uint i_block, uint i_cell, uint j_block, uint j_cell,
               Real intersection, Real intersection_di, Real intersection_dj, Real intersection_dk){
  Real a[MAX_BLOCKS_PER_DIM*WID][RECONSTRUCTION_ORDER + 1];  
  Real target[(MAX_BLOCKS_PER_DIM+2)*WID]; 
  
  
#ifdef ACC_SEMILAG_PLM
  compute_plm_coeff_explicit_column(values, blocks_per_dim, a);
#endif
#ifdef ACC_SEMILAG_PPM
  compute_ppm_coeff_explicit_column(values, blocks_per_dim, a);
#endif

  /*clear temporary taret*/
  for (uint k=0; k<WID* (blocks_per_dim + 2); ++k){ 
       target[k] = 0.0;
  }
   
   /* intersection_min is the intersection z coordinate (z after
      swaps that is) of the lowest possible z plane for each i,j
      index 
   */
  const Real intersection_min = intersection +
     (i_block * WID + i_cell) * intersection_di + 
     (j_block * WID + j_cell) * intersection_dj;
  

  /*compute some initial values, that are used to set up the
   * shifting of values as we go through all blocks in
   * order. See comments where they are shifted for
   * explanations of their meening*/

  /*loop through all blocks in column and compute the mapping as integrals*/
  for (unsigned int k_block = 0; k_block<blocks_per_dim;k_block++){
    for (uint k_cell=0; k_cell<WID; ++k_cell){ 
      /*v_l, v_r are the left and right velocity coordinates of source cell*/
      Real v_l = v_min + (k_block * WID + k_cell) * dv;
      Real v_r = v_l + dv;
      /*left(l) and right(r) k values (global index) in the target
        lagrangian grid, the intersecting cells. Again old right is new left*/               
      const int target_gk_l = (int)((v_l - intersection_min)/intersection_dk);
      const int target_gk_r = (int)((v_r - intersection_min)/intersection_dk);

      for(int gk = target_gk_l; gk <= target_gk_r; gk++){
         //the velocity limits  for the integration  to put mass
         //in the targe cell. If both v_r and v_l are in same target cell
         //then v_int_l,v_int_r should be between v_l and v_r.
         //v_int_norm_l and v_int_norm_r normalized to be between 0 and 1 in the cell.
         const Real v_int_l = min( max((Real)(gk) * intersection_dk + intersection_min, v_l), v_r);
         const Real v_int_norm_l = (v_int_l - v_l)/dv;
         const Real v_int_r = min((Real)(gk + 1) * intersection_dk + intersection_min, v_r);
         const Real v_int_norm_r = (v_int_r - v_l)/dv;
         
         /*compute left and right integrand*/
#ifdef ACC_SEMILAG_PLM
         Real target_density_l =
            v_int_norm_l * a[k_block * WID + k_cell][0] +
            v_int_norm_l * v_int_norm_l * a[k_block * WID + k_cell][1];
         Real target_density_r =
            v_int_norm_r * a[k_block * WID + k_cell][0] +
            v_int_norm_r * v_int_norm_r * a[k_block * WID + k_cell][1];
#endif
#ifdef ACC_SEMILAG_PPM
         Real target_density_l =
            v_int_norm_l * a[k_block * WID + k_cell][0] +
            v_int_norm_l * v_int_norm_l * a[k_block * WID + k_cell][1] +
            v_int_norm_l * v_int_norm_l * v_int_norm_l * a[k_block * WID + k_cell][2];
         Real target_density_r =
            v_int_norm_r * a[k_block * WID + k_cell][0] +
            v_int_norm_r * v_int_norm_r * a[k_block * WID + k_cell][1] +
            v_int_norm_r * v_int_norm_r * v_int_norm_r * a[k_block * WID + k_cell][2];
#endif
         /*total value of integrand*/
         target[gk + WID] +=  target_density_r - target_density_l;
      }
    }
  }
  /*copy target to values, and clear target array*/
  for (unsigned int k_block = 0; k_block<blocks_per_dim;k_block++){
     for (uint k=0; k<WID; ++k){ 
        values[k_block * WID + k + WID] = target[k_block * WID + k + WID];
     }
  }
}

int main(void) {
  const int dv = 20000;
  const Real v_min = -4e6;
  const int blocks_per_dim = 100;
  const int i_block = 0; //x index of block, fixed in this simple test
  const int j_block = 0; //y index of block, fixed in this simple test
  const int i_cell = 0; // z index of cell within block (0..WID-1)
  const int j_cell = 0; // y index of cell within block (0..WID-1)    
  

  Real values[(blocks_per_dim+2)*WID];
   
  /*initial values*/
  
  Real intersection = v_min + 0.6*dv;
  Real intersection_di = dv/4.0;
  Real intersection_dk = dv; 
  Real intersection_dj = dv; //does not matter here, fixed j.
  
  const int iterations=1000;
  
  /*clear target & values array*/
  for (uint k=0; k<WID* (blocks_per_dim + 2); ++k){ 
     values[k] = 0.0;
  }

 /*Add square wave*/
 for(int i=0; i < blocks_per_dim * WID; i++){
   Real v=v_min + i * dv;
   if (v > v_min +  0.8 * (blocks_per_dim * WID * dv) &&
       v < v_min +  0.9 * (blocks_per_dim * WID * dv))
     values[i + WID] = 1.0;
 }
  

/*loop over propagations*/
 for(int step = 0; step < iterations; step++){
   if(step % 10 ==0)
      print_values(step,values,blocks_per_dim, v_min, dv);
   propagate(values, blocks_per_dim, v_min, dv,
             i_block, i_cell, j_block, j_cell,
             intersection, intersection_di, intersection_dj, intersection_dk);
 }
}
