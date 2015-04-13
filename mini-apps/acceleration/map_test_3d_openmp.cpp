#include <stdio.h>
#include "common.h"
#include "cpu_1d_column_interpolations.hpp"

#define index(i,j,k)   ( k + WID + j * (blocks_per_dim_z + 2) * WID + i * (blocks_per_dim_z + 2) * blocks_per_dim_y * WID2 )
#define colindex(i,j)   ( j * (blocks_per_dim_z + 2) * WID + i * (blocks_per_dim_z + 2) * blocks_per_dim_y * WID2 )

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


void propagate(const Real * const values_in, Real *values_out, 
               uint  blocks_per_dim_x, uint blocks_per_dim_y, uint blocks_per_dim_z, 
               Real v_min, Real dv,
               Real intersection, Real intersection_di, Real intersection_dj, Real intersection_dk){
#pragma omp parallel for    
  for (uint k=0; k< (blocks_per_dim_z+2) * blocks_per_dim_x * blocks_per_dim_z * WID3; ++k){ 
    values_out[k] = 0.0;
  }

#pragma omp parallel for collapse(2)  
  for(int i = 0; i < blocks_per_dim_x * WID; i++){
    for(int j = 0; j < blocks_per_dim_y * WID; j++){
      for (uint k = 0; k < blocks_per_dim_z * WID; k++){   
         const int i_block = i / WID;
         const int i_cell = i % WID;
         const int j_block = j / WID;
         const int j_cell = j % WID;
         const Real * const values = values_in + colindex(i,j);
         
         Real a[RECONSTRUCTION_ORDER + 1];
#ifdef ACC_SEMILAG_PLM
         const Real d_cv=slope_limiter(values[k - 1 + WID], values[k + WID], values[k + 1 + WID]);
         a[0] = values[k + WID] - d_cv * 0.5;
         a[1] = d_cv * 0.5;
#endif
#ifdef ACC_SEMILAG_PPM
         // TODO!
         cerr << "PPM not done yet"<<endl;
         exit(1);
#endif
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
         /*v_l, v_r are the left and right velocity coordinates of source cell*/
         const Real v_l = v_min + k * dv;
         const Real v_r = v_l + dv;
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
               v_int_norm_l * a[0] +
               v_int_norm_l * v_int_norm_l * a[1];
            Real target_density_r =
               v_int_norm_r * a[0] +
               v_int_norm_r * v_int_norm_r * a[1];
#endif
#ifdef ACC_SEMILAG_PPM
            Real target_density_l =
               v_int_norm_l * a[0] +
               v_int_norm_l * v_int_norm_l * a[1] +
               v_int_norm_l * v_int_norm_l * v_int_norm_l * a[2];
            Real target_density_r =
               v_int_norm_r * a[0] +
               v_int_norm_r * v_int_norm_r * a[1] +
               v_int_norm_r * v_int_norm_r * v_int_norm_r * a[2];
#endif
            /*total value of integrand, if it is wihtin bounds*/
            if ( gk >= 0 && gk <= blocks_per_dim_z * WID )
            //atomic not needed if k index is not threaded
            //#pragma omp atomic update
               values_out[colindex(i,j) + gk + WID] +=  target_density_r - target_density_l;
            }
      }
    }
  }
}


int main(void) {
  /*define grid size*/
  const int dv = 20000;
  const Real v_min = -2e6;
  const int blocks_per_dim_x = 10;
  const int blocks_per_dim_y = 10;
  const int blocks_per_dim_z = 50;
  

  Real *values_a = new Real[(blocks_per_dim_z+2) * blocks_per_dim_x * blocks_per_dim_z * WID3];
  Real *values_b = new Real[(blocks_per_dim_z+2) * blocks_per_dim_x * blocks_per_dim_z * WID3];
   
  /*intersection values define the acceleration transformation. These would be obtained from other routines, but are here fixed*/
  Real intersection = v_min + 0.6*dv;
  Real intersection_di = dv/4.0;
  Real intersection_dk = dv; 
  Real intersection_dj = dv; //does not matter here, fixed j.
  
  const int iterations=1000;
  
  /*clear target & values array*/
  for (uint k=0; k< (blocks_per_dim_z+2) * blocks_per_dim_x * blocks_per_dim_z * WID3; ++k){ 
     values_a[k] = 0.0;
  }

 /*Add square wave*/
  for(int i=0; i < blocks_per_dim_x * WID; i++){
    for(int j=0; j < blocks_per_dim_y * WID; j++){
      for(int k=0; k < blocks_per_dim_z * WID; k++){
         Real v = v_min + k * dv;
         if (v > v_min +  0.8 * (blocks_per_dim_z * WID * dv) &&
             v < v_min +  0.9 * (blocks_per_dim_z * WID * dv))
            values_a[index(i,j,k)] = 1.0;
      }
    }
  }

  /*loop over propagations*/
  for(int step = 0; step < iterations; step+=2){
    if(step % 10 ==0)
      print_values(step, values_a + colindex(0,0), blocks_per_dim_z, v_min, dv);
    propagate(values_a, values_b,
              blocks_per_dim_x, blocks_per_dim_y, blocks_per_dim_z, 
              v_min, dv,
              intersection, intersection_di, intersection_dj, intersection_dk);
    propagate(values_b, values_a,
              blocks_per_dim_x, blocks_per_dim_y, blocks_per_dim_z,
              v_min, dv,
              intersection, intersection_di, intersection_dj, intersection_dk);

  }
}

