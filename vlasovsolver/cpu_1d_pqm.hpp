/*
This file is part of Vlasiator.
Copyright 2013, 2014 Finnish Meteorological Institute

*/

#ifndef CPU_1D_PQM_H
#define CPU_1D_PQM_H

#include <iostream>
#include "vec4.h"
#include "algorithm"
#include "cmath"
#include "cpu_slope_limiters.hpp"
#include "cpu_face_estimates.hpp"

using namespace std;



/*make sure quartic polynomial is monotonic*/
inline void filter_pqm_monotonicity(Vec4 *values, uint k, Vec4 &fv_l, Vec4 &fv_r, Vec4 &fd_l, Vec4 &fd_r){   
  const Vec4 root_outside = Vec4(100.0); //fixed values give to roots clearly outside [0,1], or nonexisting ones*/
  /*second derivative coefficients, eq 23 in white et al.*/
  Vec4 b0 =   60.0 * values[k] - 24.0 * fv_r - 36.0 * fv_l + 3.0 * (fd_r - 3.0 * fd_l);
  Vec4 b1 = -360.0 * values[k] + 36.0 * fd_l - 24.0 * fd_r + 168.0 * fv_r + 192.0 * fv_l;
  Vec4 b2 =  360.0 * values[k] + 30.0 * (fd_r - fd_l) - 180.0 * (fv_l + fv_r);
  /*let's compute sqrt value to be used for computing roots. If we
    take sqrt of negaitve numbers, then we instead set a value that
    will make the root to be +-100 which is well outside range
    of[0,1]. We also guard the sqrt against sqrt with negative
    numbers by doing a max*/
  const Vec4 sqrt_val = select(b1 * b1 - 4 * b0 * b2 < 0.0, 
			       b1 + 200.0 * b2,
			       sqrt(max(b1 * b1- 4 * b0 * b2, 0.0))); 
  //compute roots. Division is safe with vectorclass (=inf)
  const Vec4 root1 = (-b1 + sqrt_val) / (2 * b2);
  const Vec4 root2 = (-b1 - sqrt_val) / (2 * b2);

  /*PLM slope, MC limiter*/
  Vec4 plm_slope_l = 2.0 * (values[k] - values[k - 1]);
  Vec4 plm_slope_r = 2.0 * (values[k + 1] - values[k]);
  Vec4 slope_sign = plm_slope_l + plm_slope_r; //it also has some magnitude, but we will only use its sign.
  /*first derivative coefficients*/
  const Vec4 c0 = fd_l;
  const Vec4 c1 = b0;
  const Vec4 c2 = b1 / 2.0;
  const Vec4 c3 = b2 / 3.0;
  //compute both slopes at inflexion points, at least one of these
  //is with [0..1]. If the root is not in this range, we
  //simplify later if statements by setting it to the plm slope
  //sign
  Vec4 root1_slope = select(root1 >= 0.0 && root1 <= 1.0, 
			    c0  + c1 * root1 + c2 * root1 * root1 + c3 * root1 * root1 * root1,
			    slope_sign);
  Vec4 root2_slope = select(root2 >= 0.0 && root2 <= 1.0, 
			    c0  + c1 * root2 + c2 * root2 * root2 + c3 * root2 * root2 * root2,
			    slope_sign);
  if (horizontal_or (root1_slope * slope_sign < 0.0 || root2_slope * slope_sign < 0.0 )) {
    //serialized the handling of inflexion points, these do not happen for smooth regions
    for(uint i = 0;i < WID; i++) {
      if(root1_slope[i] * slope_sign[i] < 0.0 || root2_slope[i] * slope_sign[i] < 0.0 ) {
	//need to collapse, at least one inflexion point has wrong
	//sign.
	if(fabs(plm_slope_l[i]) <= fabs(plm_slope_r[i])) {
	  //collapse to left edge (eq 21)
	  fd_l.insert( i, 1.0 / 3.0 * ( 10 * values[k][i] - 2.0 * fv_r[i] - 8.0 * fv_l[i]));
	  fd_r.insert( i, -10.0 * values[k][i] + 6.0 * fv_r[i] + 4.0 * fv_l[i]);
	  //check if PLM slope is consistent (eq 28 & 29)
	  if (slope_sign[i] * fd_l[i] < 0) {
	    fd_l.insert( i, 0);
	    fv_r.insert( i, 5 * values[k][i] - 4 * fv_l[i]);
	    fd_r.insert( i, 20 * (values[k][i] - fv_l[i]));
	  }
	  else if (slope_sign[i] * fd_r[i] < 0) {
	    fd_r.insert( i, 0);
	    fv_l.insert( i, 0.5 * (5 * values[k][i] - 3 * fv_r[i]));
	    fd_l.insert( i, 10.0 / 3.0 * (-values[k][i] + fv_r[i]));
	  }
	}
	else {
	  //collapse to right edge (eq 21)
	  fd_l.insert( i, 10.0 * values[k][i] - 6.0 * fv_l[i] - 4.0 * fv_r[i]);
	  fd_r.insert( i, 1.0 / 3.0 * ( - 10.0 * values[k][i] + 2 * fv_l[i] + 8 * fv_r[i]));
	  //check if PLM slope is consistent (eq 28 & 29)
	  if (slope_sign[i] * fd_l[i] < 0) {
	    fd_l.insert( i, 0);
	    fv_r.insert( i, 0.5 * ( 5 * values[k][i] - 3 * fv_l[i]));
	    fd_r.insert( i, 10.0 / 3.0 * (values[k][i] - fv_l[i]));
	  }
	  else if (slope_sign[i] * fd_r[i] < 0) {
	    fd_r.insert( i, 0);
	    fv_l.insert( i, 5 * values[k][i] - 4 * fv_r[i]);
	    fd_l.insert( i, 20.0 * ( - values[k][i] + fv_r[i]));
	  }
	}
      }
    }
  }
}



// /*
//   PQM reconstruction as published in:
//   White, Laurent, and Alistair Adcroft. “A High-Order Finite Volume Remapping Scheme for Nonuniform Grids: The Piecewise Quartic Method (PQM).” Journal of Computational Physics 227, no. 15 (July 2008): 7394–7422. doi:10.1016/j.jcp.2008.04.026.
// */

inline void compute_pqm_coeff(Vec4 *values, face_estimate_order order, uint k, Vec4 a[5]){
   Vec4 fv_l; /*left face value*/
   Vec4 fv_r; /*right face value*/
   Vec4 fd_l; /*left face derivative*/
   Vec4 fd_r; /*right face derivative*/

   compute_filtered_face_values_derivatives(values, k, order, fv_l, fv_r, fd_l, fd_r); 
   filter_pqm_monotonicity(values, k, fv_l, fv_r, fd_l, fd_r); 

  //Fit a second order polynomial for reconstruction see, e.g., White
  //2008 (PQM article) (note additional integration factors built in,
  //contrary to White (2008) eq. 4
  a[0] = fv_l;
  a[1] = fd_l/2.0;
  a[2] =  10.0 * values[k] - 4.0 * fv_r - 6.0 * fv_l + 0.5 * (fd_r - 3 * fd_l);
  a[3] = -15.0 * values[k]  + 1.5 * fd_l - fd_r + 7.0 * fv_r + 8 * fv_l;
  a[4] =   6.0 * values[k] +  0.5 * (fd_r - fd_l) - 3.0 * (fv_l + fv_r);
}

#endif
