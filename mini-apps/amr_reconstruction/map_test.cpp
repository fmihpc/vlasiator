#include <stdio.h>
#include "common.h"
#include "vlasovsolver/vec.h"
//#include "vlasovsolver/cpu_1d_ppm.hpp"
//#include "vlasovsolver/cpu_1d_ppm_nonuniform.hpp"
#include "vlasovsolver/cpu_1d_ppm_nonuniform_conserving.hpp"
#include <random>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <iostream>

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
   for(uint i=0; i < blocks_per_dim * WID; i++){
      Real v=v_min + (i + 0.5)*dv;
      fprintf(fp,"%20.12g %20.12g %20.12g %20.12g %20.12g\n", v, values[i + WID][0], values[i + WID][1], values[i + WID][2], values[i + WID][3]);
   }
   fclose(fp);
}

void propagate(Vec dr[], Vec values[], Real z_translation, uint blocks_per_dim ) {
   
   // Determine direction of translation
   // part of density goes here (cell index change along spatial direcion)
   const int target_scell_index = (z_translation > 0) ? 1: -1; 
  
   // Vector buffer where we write data, initialized to 0*/
   Vec targetValues[(blocks_per_dim + 2) * WID];
   
   for (uint k_block = 0; k_block < blocks_per_dim; k_block++){
      
      for (uint k_cell=0; k_cell < WID; ++k_cell) {
         
         uint gid = k_block * WID + k_cell + WID;
         // init target_values
         targetValues[gid] = 0.0;
         
      }
   }
   for (uint k_block = 0; k_block < blocks_per_dim; k_block++){
      
      for (uint k_cell=0; k_cell < WID; ++k_cell){
         
         uint gid = k_block * WID + k_cell + WID;
         //uint gid = (blocks_per_dim + 2) * WID - (k_block * WID + k_cell + WID);
         
         // Calculate normalized coordinates in current cell.
         // The coordinates (scaled units from 0 to 1) between which we will
         // integrate to put mass in the target  neighboring cell.
         // Normalize the coordinates to the origin cell. Then we scale with the difference
         // in volume between target and origin later when adding the integrated value.
         Realv z_1,z_2;
         if ( z_translation < 0 ) {
            z_1 = 0;
            z_2 = -z_translation / dr[gid][0]; 
         } else {
            z_1 = 1.0 - z_translation / dr[gid][0]; 
            z_2 = 1.0;
         }
         
         if( abs(z_1) > 1.0 || abs(z_2) > 1.0 ) {
            std::cout << "Error, CFL condition violated\n";
            std::cout << "Exiting\n";
            std::exit(1);
         }
         
         // Compute polynomial coefficients
         Vec a[3];
         //compute_ppm_coeff_nonuniform(dr, values, h4, gid + target_scell_index, a);
         compute_ppm_coeff_nonuniform(dr, values, h4, gid, a, fluxlimiterscalingfactor);
         
         // Compute integral
         const Vec ngbr_target_density =
            z_2 * ( a[0] + z_2 * ( a[1] + z_2 * a[2] ) ) -
            z_1 * ( a[0] + z_1 * ( a[1] + z_1 * a[2] ) );

         // Store mapped density in two target cells
         // in the neighbor cell we will put this density        
         targetValues[gid + target_scell_index] +=  ngbr_target_density * dr[gid] / dr[gid + target_scell_index];
         // in the current original cells we will put the rest of the original density
         targetValues[gid]                      +=  values[gid] - ngbr_target_density;
      }
   }

   // Store target data into source data
   for (uint k_block = 0; k_block<blocks_per_dim;k_block++){
    
      for (uint k_cell=0; k_cell<WID; ++k_cell){

         uint gid = k_block * WID + k_cell + WID;
         //uint gid = (blocks_per_dim + 2) * WID - (k_block * WID + k_cell + WID);
         values[gid] = targetValues[gid];
      
      }
    
   }  

}

void print_reconstruction(int step, Vec dr[], Vec values[], uint  blocks_per_dim, Real r_min){
   char name[256];
   sprintf(name,"reconstructions_%05d.dat",step);
   FILE* fp=fopen(name,"w");

   Vec r0 = r_min;
   const int subcells = 50;
   /*loop through all blocks in column and divide into subcells. Print value of reconstruction*/
   for (unsigned int k_block = 0; k_block<blocks_per_dim;k_block++){
      for (uint k_cell=0; k_cell<WID; ++k_cell){ 
#ifdef ACC_SEMILAG_PPM
         Vec a[3];
         //compute_ppm_coeff(               values, h4, (k_block + 1) * WID + k_cell, a);
         compute_ppm_coeff_nonuniform(dr, values, h4, (k_block + 1) * WID + k_cell, a, fluxlimiterscalingfactor);
#endif     
      
         int iend = k_block * WID + k_cell;
         if (iend > 0)
            r0 += dr[iend-1+WID];
	
      
         for (uint k_subcell=0; k_subcell< subcells; ++k_subcell){ 
            Vec r_norm = (Real)(k_subcell + 0.5)/subcells; //normalized r of subcell in source cell
            Vec r = r0 + r_norm * dr[k_block * WID + k_cell + WID];
	
#ifdef ACC_SEMILAG_PPM
            Vec target = 
               a[0] +
               2.0 * r_norm * a[1] +
               3.0 * r_norm * r_norm * a[2];
#endif

            fprintf(fp,"%20.12g %20.12e %20.12e\n", r[0], values[k_block * WID + k_cell + WID][0], target[0]);
         }
         //fprintf(fp,"\n"); //empty line to deay wgments in gnuplot
      }
   }
  
   fclose(fp);
}

void refine(Vec dr[], int ir, int max_refinement, int cells_per_level) {

   for (uint k=0; k < max_refinement * cells_per_level; ++k) {
      dr[ir + k] = dr[ir + k]/pow(2,(max_refinement - k / cells_per_level));
      if (k > 0)
         dr[ir - k] = dr[ir - k]/pow(2,(max_refinement - k / cells_per_level));
   }

}

int main(void) {
  
   const Real dr0 = 20000;
   const int blocks_per_dim = 100;
   const int i_block = 0; //x index of block, fixed in this simple test
   const int j_block = 0; //y index of block, fixed in this simple test
   const int j_cell = 0; // y index of cell within block (0..WID-1)

   Vec dr[(blocks_per_dim+2)*WID];
   Vec values[(blocks_per_dim+2)*WID];

   boost::mt19937 rng;
   boost::uniform_real<Real> u(0.0, 2.0 * dr0);
   boost::variate_generator<boost::mt19937&, boost::uniform_real<Real> > gen(rng, u);
   gen.distribution().reset();
   gen.engine().seed(12345);
  
   /*initial values*/  
   /*clear target & values array*/
   for (uint k=0; k<WID* (blocks_per_dim + 2); ++k){ 
      values[k] = 0.0;
      dr[k] = dr0;
      //dr[k] = gen();
   }

   int ir = (blocks_per_dim + 2) * WID / 2;
   int ir2 = (blocks_per_dim + 2) * WID / 3;
   int max_refinement = 5;
   int cells_per_level = 2;

   refine(dr,ir,max_refinement,cells_per_level);
   // refine(dr,ir2,max_refinement,cells_per_level);
 
   Real r_min = 0.0;
   for (uint k=WID;k < (blocks_per_dim + 2) * WID / 2; ++k) {    
      r_min -= dr[k][0];
   }
  
   Real T = 500000;
   Real rho = 1.0e18;
   Real r = r_min;
   Real r1 = -10.0 * dr0;
  
   for(uint i=0; i < blocks_per_dim * WID; i++){

      // Evaluate the function at the middle of the cell
      r = r + 0.5 * dr[i + WID][0];
      values[i + WID] = rho * pow(physicalconstants::MASS_PROTON / (2.0 * M_PI * physicalconstants::K_B * T), 1.5) *
         exp(- physicalconstants::MASS_PROTON * (r - r1) * (r - r1) / (2.0 * physicalconstants::K_B * T));    

      // if (r < 0.0 && r_min - 10.0 * r < 0.0) {
      //   values[i + WID] = abs(r_min - 10.0 * r);
      // } else {
      //   values[i + WID] = 0.0;
      // }
    
      // Move r to the end of the cell for the next iteration
      r = r + 0.5 * dr[i + WID][0];
   }
  
   print_reconstruction(0, dr, values, blocks_per_dim, r_min);

   uint nstep = 1000;
   Real step = 500.0;
  
   for (uint istep=0; istep < nstep; ++istep) {
      propagate(dr, values, step, blocks_per_dim);
      if ((istep+1) % 10 == 0)
         print_reconstruction(istep+1, dr, values, blocks_per_dim, r_min);
   }
  
}
