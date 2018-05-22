#include <stdio.h>
#include "common.h"
#include "vlasovsolver/vec.h"
#include "vlasovsolver/cpu_1d_ppm.hpp"
//#include "vlasovsolver/cpu_1d_ppm_nonuniform.hpp"
#include <random>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <iostream>

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

void print_reconstruction(int step, Vec dv[], Vec values[], uint  blocks_per_dim, Real v_min, Real dv0){
  char name[256];
  sprintf(name,"reconstructions_%03d.dat",step);
  FILE* fp=fopen(name,"w");
  FILE* fp2=fopen("a.dat","w");

  Vec v_r = v_min;
  const int subcells = 50;
  /*loop through all blocks in column and divide into subcells. Print value of reconstruction*/
  for (unsigned int k_block = 0; k_block<blocks_per_dim;k_block++){
    for (uint k_cell=0; k_cell<WID; ++k_cell){ 
#ifdef ACC_SEMILAG_PPM
      Vec a[3];
      //Vec a2[3];
      compute_ppm_coeff(               values, h6, (k_block + 1) * WID + k_cell, a);
      //compute_ppm_coeff_nonuniform(dv, values, h6, (k_block + 1) * WID + k_cell, a);
      fprintf(fp2,"%12.8g %12.8g %12.8g\n",a[0][0],a[1][0],a[2][0]);
#endif
      
      Vec v_l = v_min + (k_block * WID + k_cell) * dv0;
      
      int iend = k_block * WID + k_cell;
      if (iend > 0)
	v_r += dv[iend-1+WID];
      
      for (uint k_subcell=0; k_subcell< subcells; ++k_subcell){ 
	Vec v_norm = (Real)(k_subcell + 0.5)/subcells; //normalized v of subcell in source cell
	Vec v0 = v_l + v_norm * dv0;       

	Vec v = v_r + v_norm * dv[k_block * WID + k_cell + WID];

#ifdef ACC_SEMILAG_PPM
	Vec target = 
	  a[0] +
	  2.0 * v_norm * a[1] +
	  3.0 * v_norm * v_norm * a[2];
#endif

	fprintf(fp,"%20.12g %20.12g %20.12g\n", v[0], values[k_block * WID + k_cell + WID][0], target[0]);
      }
      //fprintf(fp,"\n"); //empty line to deay wgments in gnuplot
    }
  }
  
  fclose(fp);
  fclose(fp2);
}

int main(void) {
  
  const Real dv0 = 20000;
  const Real v_min = -4e6;
  const int blocks_per_dim = 100;
  const int i_block = 0; //x index of block, fixed in this simple test
  const int j_block = 0; //y index of block, fixed in this simple test
  const int j_cell = 0; // y index of cell within block (0..WID-1)

  Vec dv[(blocks_per_dim+2)*WID];
  Real dvTmp[(blocks_per_dim+2)*WID];
  Vec values[(blocks_per_dim+2)*WID];
  Real rnd;
  //std::random_device rd;
  //std::mt19937 e2(rd());
  //std::uniform_real_distribution<> dist(0, 10);

  boost::mt19937 rng;
  boost::uniform_real<Real> u(0.0, 2.0 * dv0);
  boost::variate_generator<boost::mt19937&, boost::uniform_real<Real> > gen(rng, u);
  
  /*initial values*/  
  /*clear target & values array*/
  for (uint k=0; k<WID* (blocks_per_dim + 2); ++k){ 
    values[k] = 0.0;
    //rnd = gen();
    rnd = dv0;
    dv[k] = rnd;
  }

  Real T = 500000;
  Real rho = 1.0e6;
  Real v = v_min;
  
  for(uint i=0; i < blocks_per_dim * WID; i++){
    // Real v=v_min + i*dv;
    values[i + WID] = rho * pow(physicalconstants::MASS_PROTON / (2.0 * M_PI * physicalconstants::K_B * T), 1.5) *
      exp(- physicalconstants::MASS_PROTON * v * v / (2.0 * physicalconstants::K_B * T));
    v = v + dv[i + WID][0];
  }
  
  // print_values(0,values,blocks_per_dim, v_min, dv);
  print_reconstruction(0, dv, values, blocks_per_dim, v_min, dv0);
}
