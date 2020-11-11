#include <stdio.h>
#include "common.h"
#include "vlasovsolver/vec.h" 
#include "vlasovsolver/cpu_1d_pqm.hpp"
#include "vlasovsolver/cpu_1d_ppm.hpp"
#include "vlasovsolver/cpu_1d_plm.hpp"

const int fluxlimiterscalingfactor=1.e-15;
// Used for better calculation of flux limiters at extreme values.
// In vlasiator, the value of spatial_cell->getVelocityBlockMinValue(popID)
// is used here.

/*print all values in the vector valued values array. In this array
  there are blocks_per_dim blocks with a width of WID*/
void print_values(int step, Vec *values, uint blocks_per_dim, Real v_min, Real dv){
  char name[256];
  sprintf(name,"dist_%03d.dat",step);

  FILE* fp=fopen(name,"w");
  for(int i=0; i < blocks_per_dim * WID; i++){
     Real v=v_min + (i + 0.5)*dv;
     fprintf(fp,"%20.12g %20.12g %20.12g %20.12g %20.12g\n", v, values[i + WID][0], values[i + WID][1], values[i + WID][2], values[i + WID][3]);
  }
  fclose(fp);
}


void propagate(Vec values[], uint  blocks_per_dim, Real v_min, Real dv,
               uint i_block, uint j_block, uint j_cell,
               Real intersection, Real intersection_di, Real intersection_dj, Real intersection_dk){
  Vec target[(MAX_BLOCKS_PER_DIM+2)*WID]; 

  
  /*clear temporary taret*/
  for (uint k=0; k<WID* (blocks_per_dim + 2); ++k){ 
       target[k] = 0.0;
  }
   
   /* intersection_min is the intersection z coordinate (z after
      swaps that is) of the lowest possible z plane for each i,j
      index (i in vector)
   */
  const Real intersection_min_base =  
    intersection +
    (i_block * WID) * intersection_di + 
    (j_block * WID + j_cell) * intersection_dj;
  
  //const Vec intersection_min(intersection_min_base);
  const Vec intersection_min(intersection_min_base + 0 * intersection_di,
                              intersection_min_base + 1 * intersection_di,
                              intersection_min_base + 2 * intersection_di,
                              intersection_min_base + 3 * intersection_di);
  /*compute some initial values, that are used to set up the
   * shifting of values as we go through all blocks in
   * order. See comments where they are shifted for
   * explanations of their meening*/

  /*loop through all blocks in column and compute the mapping as integrals*/
  for (unsigned int k_block = 0; k_block<blocks_per_dim;k_block++){
    for (uint k_cell=0; k_cell<WID; ++k_cell){ 

#ifdef ACC_SEMILAG_PLM
      Vec a[2];
      compute_plm_coeff(values, (k_block + 1) * WID + k_cell , a, fluxlimiterscalingfactor);
#endif
#ifdef ACC_SEMILAG_PPM
      Vec a[3];
      compute_ppm_coeff(values, h6, (k_block + 1) * WID + k_cell , a, fluxlimiterscalingfactor);
#endif
#ifdef ACC_SEMILAG_PQM
      Vec a[5];
      compute_pqm_coeff(values, h8,  (k_block + 1) * WID + k_cell , a, fluxlimiterscalingfactor);
#endif



      /*v_l, v_r are the left and right velocity coordinates of source cell*/
      Vec v_l = v_min + (k_block * WID + k_cell) * dv;
      Vec v_r = v_l + dv;
      /*left(l) and right(r) k values (global index) in the target
        lagrangian grid, the intersecting cells. Again old right is new left*/               
      const Veci target_gk_l = truncate_to_int((v_l - intersection_min)/intersection_dk);
      const Veci target_gk_r = truncate_to_int((v_r - intersection_min)/intersection_dk);
      
      

      Veci gk(target_gk_l);
      while (horizontal_or(gk <= target_gk_r)){
         //the velocity limits  for the integration  to put mass
         //in the targe cell. If both v_r and v_l are in same target cell
         //then v_int_l,v_int_r should be between v_l and v_r.
         //v_int_norm_l and v_int_norm_r normalized to be between 0 and 1 in the cell.

#ifdef DP
         const Vec v_int_l = min( max(to_double(gk) * intersection_dk + intersection_min, v_l), v_r);
         const Vec v_int_norm_l = (v_int_l - v_l)/dv;
         const Vec v_int_r = min(to_double(gk + 1) * intersection_dk + intersection_min, v_r);
         const Vec v_int_norm_r = (v_int_r - v_l)/dv;
#else
         const Vec v_int_l = min( max(to_float(gk) * intersection_dk + intersection_min, v_l), v_r);
         const Vec v_int_norm_l = (v_int_l - v_l)/dv;
         const Vec v_int_r = min(to_float(gk + 1) * intersection_dk + intersection_min, v_r);
         const Vec v_int_norm_r = (v_int_r - v_l)/dv;
#endif

         /*compute left and right integrand*/
#ifdef ACC_SEMILAG_PLM
         Vec target_density_l =
            v_int_norm_l * a[0] +
            v_int_norm_l * v_int_norm_l * a[1];
         Vec target_density_r =
            v_int_norm_r * a[0] +
            v_int_norm_r * v_int_norm_r * a[1];
#endif
#ifdef ACC_SEMILAG_PPM
         Vec target_density_l =
            v_int_norm_l * a[0] +
            v_int_norm_l * v_int_norm_l * a[1] +
            v_int_norm_l * v_int_norm_l * v_int_norm_l * a[2];
         Vec target_density_r =
            v_int_norm_r * a[0] +
            v_int_norm_r * v_int_norm_r * a[1] +
            v_int_norm_r * v_int_norm_r * v_int_norm_r * a[2];
#endif
#ifdef ACC_SEMILAG_PQM
         Vec target_density_l =
            v_int_norm_l * a[0] +
            v_int_norm_l * v_int_norm_l * a[1] +
            v_int_norm_l * v_int_norm_l * v_int_norm_l * a[2] +
            v_int_norm_l * v_int_norm_l * v_int_norm_l * v_int_norm_l * a[3] +
            v_int_norm_l * v_int_norm_l * v_int_norm_l * v_int_norm_l * v_int_norm_l * a[4];

         Vec target_density_r =
            v_int_norm_r * a[0] +
            v_int_norm_r * v_int_norm_r * a[1] +
            v_int_norm_r * v_int_norm_r * v_int_norm_r * a[2] +
            v_int_norm_r * v_int_norm_r * v_int_norm_r * v_int_norm_r * a[3] +
            v_int_norm_r * v_int_norm_r * v_int_norm_r * v_int_norm_r * v_int_norm_r * a[4];

#endif
         
         
         /*total value of integrand*/
         const Vec target_density = target_density_r - target_density_l;
         
         //store values, one element at elema time
         for(uint elem = 0; elem < 4;elem ++ ){
            int k_in_target = gk[elem];
            if (k_in_target >=0 &&
               k_in_target < blocks_per_dim * WID) {
               const Real new_density = target[k_in_target + WID][elem] + target_density[elem];
               target[k_in_target + WID].insert(elem, new_density);
            }
         }
         gk++; //next iteration in while loop
      }
     }
  }
  

  /*copy target to values, and clear target array*/
  for (unsigned int k_block = 0; k_block<blocks_per_dim;k_block++){
    for (uint k=0; k<WID; ++k){ 
      values[k_block * WID + k + WID] = target[k_block * WID + k + WID];
      target[k_block * WID + k + WID] = 0.0;
     }
  }
}

void print_reconstruction(int step, Vec values[], uint  blocks_per_dim, Real v_min, Real dv,
                          uint i_block, uint j_block, uint j_cell,
                          Real intersection, Real intersection_di, Real intersection_dj, Real intersection_dk){
  char name[256];
  sprintf(name,"reconstructions_%05d.dat",step);
  FILE* fp=fopen(name,"w");

  
   
  /* intersection_min is the intersection z coordinate (z after
     swaps that is) of the lowest possible z plane for each i,j
     index (i in vector)
  */
  const Real intersection_min_base =  
    intersection +
    (i_block * WID) * intersection_di + 
    (j_block * WID + j_cell) * intersection_dj;
  
  //const Vec intersection_min(intersection_min_base);
  const Vec intersection_min(intersection_min_base + 0 * intersection_di,
                              intersection_min_base + 1 * intersection_di,
                              intersection_min_base + 2 * intersection_di,
                              intersection_min_base + 3 * intersection_di);
  /*compute some initial values, that are used to set up the
   * shifting of values as we go through all blocks in
   * order. See comments where they are shifted for
   * explanations of their meening*/
  const int subcells = 50;
  /*loop through all blocks in column and divide into subcells. Print value of reconstruction*/
  for (unsigned int k_block = 0; k_block<blocks_per_dim;k_block++){
    for (uint k_cell=0; k_cell<WID; ++k_cell){ 
#ifdef ACC_SEMILAG_PLM
      Vec a[2];
      compute_plm_coeff(values, (k_block + 1) * WID + k_cell , a, fluxlimiterscalingfactor);
#endif
#ifdef ACC_SEMILAG_PPM
      Vec a[3];
      compute_ppm_coeff(values, h6, (k_block + 1) * WID + k_cell , a, fluxlimiterscalingfactor);
#endif
#ifdef ACC_SEMILAG_PQM
      Vec a[5];
      compute_pqm_coeff(values, h8, (k_block + 1) * WID + k_cell , a, fluxlimiterscalingfactor);
#endif

      Vec v_l = v_min + (k_block * WID + k_cell) * dv;
      for (uint k_subcell=0; k_subcell< subcells; ++k_subcell){ 
           Vec v_norm = (Real)(k_subcell + 0.5)/subcells; //normalized v of subcell in source cell
           Vec v = v_l + v_norm * dv;

#ifdef ACC_SEMILAG_PLM
         Vec target = 
            a[0] +
            2.0 * v_norm * a[1];
#endif
#ifdef ACC_SEMILAG_PPM
         Vec target = 
            a[0] +
            2.0 * v_norm * a[1] +
            3.0 * v_norm * v_norm * a[2];
#endif
#ifdef ACC_SEMILAG_PQM
         Vec target = 
            a[0] +
            2.0 * v_norm * a[1] +
            3.0 * v_norm * v_norm * a[2] +
            4.0 * v_norm * v_norm * v_norm * a[3] +
            5.0 * v_norm * v_norm * v_norm * v_norm * a[4];
#endif
         fprintf(fp,"%20.12g %20.12g %20.12g\n", v[0], values[k_block * WID + k_cell + WID][0], target[0]);
      }
      fprintf(fp,"\n"); //empty line to deay wgments in gnuplot
    }
  }
  
  fclose(fp);
}






int main(void) {
  const int dv = 20000;
  const Real v_min = -4e6;
  const int blocks_per_dim = 100;
  const int i_block = 0; //x index of block, fixed in this simple test
  const int j_block = 0; //y index of block, fixed in this simple test
  const int j_cell = 0; // y index of cell within block (0..WID-1)
  

  Vec values[(blocks_per_dim+2)*WID];
   
  /*initial values*/
  
  Real intersection = v_min - 0.1*dv;
  Real intersection_di = 0.025 * dv;
  Real intersection_dk = dv;
  Real intersection_dj = 0.0 * dv; //does not matter here, fixed j.

  
  const int iterations = 1000;

   /*clear target & values array*/
  for (uint k=0; k<WID* (blocks_per_dim + 2); ++k){ 
       values[k] = 0.0;
  }

 /*Add square wave*/
  /*
  for(int i=0; i < blocks_per_dim * WID; i++){
     Real v=v_min + i*dv;
     if (v > v_min +  0.8 * (blocks_per_dim * WID * dv) & 
         v < v_min +  0.9 * (blocks_per_dim * WID * dv)) {
        values[i + WID] = Vec(1.0);
     }
  }
*/
  Real T = 500000;
  Real rho = 1.0e18;
  for(int i=0; i < blocks_per_dim * WID; i++){
     Real v=v_min + i*dv;
     values[i + WID] = rho * pow(physicalconstants::MASS_PROTON / (2.0 * M_PI * physicalconstants::K_B * T), 1.5) *
        exp(- physicalconstants::MASS_PROTON * v * v / (2.0 * physicalconstants::K_B * T));
  }
     
  
// print_values(0,values,blocks_per_dim, v_min, dv);
 print_reconstruction(0, values, blocks_per_dim, v_min, dv,
                      i_block, j_block, j_cell,
                      intersection, intersection_di, intersection_dj, intersection_dk);

 clock_t t = clock();
/*loop over propagations*/
 for(int step = 0; step <= iterations; step++){
   propagate(values, blocks_per_dim, v_min, dv,
             i_block, j_block, j_cell,
             intersection, intersection_di, intersection_dj, intersection_dk);
   if (step % 10 == 0)
     print_reconstruction(step, values, blocks_per_dim, v_min, dv,
			  i_block, j_block, j_cell,
			  intersection, intersection_di, intersection_dj, intersection_dk);
   
 }
 printf("\nTime per iteration: %12.15g\n", ((double)(clock() - t)/CLOCKS_PER_SEC)/iterations);

// print_values(iterations,values,blocks_per_dim, v_min, dv);
}
