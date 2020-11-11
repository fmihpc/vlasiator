/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://www.physics.helsinki.fi/vlasiator/
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#ifndef CPU_1D_PQM_H
#define CPU_1D_PQM_H

#include <iostream>
#include "vec.h"
#include "algorithm"
#include "cmath"
#include "cpu_slope_limiters.hpp"
#include "cpu_face_estimates.hpp"

using namespace std;



/*make sure quartic polynomial is monotonic*/
inline void filter_pqm_monotonicity(Vec *values, uint k, Vec &fv_l, Vec &fv_r, Vec &fd_l, Vec &fd_r){   
   const Vec root_outside = Vec(100.0); //fixed values give to roots clearly outside [0,1], or nonexisting ones*/
   /*second derivative coefficients, eq 23 in white et al.*/
   Vec b0 =   60.0 * values[k] - 24.0 * fv_r - 36.0 * fv_l + 3.0 * (fd_r - 3.0 * fd_l);
   Vec b1 = -360.0 * values[k] + 36.0 * fd_l - 24.0 * fd_r + 168.0 * fv_r + 192.0 * fv_l;
   Vec b2 =  360.0 * values[k] + 30.0 * (fd_r - fd_l) - 180.0 * (fv_l + fv_r);
   /*let's compute sqrt value to be used for computing roots. If we
    take sqrt of negaitve numbers, then we instead set a value that
    will make the root to be +-100 which is well outside range
    of[0,1]. We do not catch FP exceptions, so sqrt(negative) are okish (add
    a max(val_to_sqrt,0) if not*/
   const Vec val_to_sqrt = b1 * b1 - 4 * b0 * b2;
#ifdef VEC16F_AGNER
   //this sqrt gives 10% more perf on acceleration on KNL. Also fairly
   //accurate with AVX512ER. On Xeon it is not any faster, and less accurate.
   const Vec sqrt_val = select(val_to_sqrt < 0.0, 
                               b1 + 200.0 * b2,
                               val_to_sqrt * approx_rsqrt(val_to_sqrt));
#else
   const Vec sqrt_val = select(val_to_sqrt < 0.0, 
                               b1 + 200.0 * b2,
                               sqrt(val_to_sqrt));
#endif
   //compute roots. Division is safe with vectorclass (=inf)
   const Vec root1 = (-b1 + sqrt_val) / (2 * b2);
   const Vec root2 = (-b1 - sqrt_val) / (2 * b2);

   /*PLM slope, MC limiter*/
   Vec plm_slope_l = 2.0 * (values[k] - values[k - 1]);
   Vec plm_slope_r = 2.0 * (values[k + 1] - values[k]);
   Vec slope_sign = plm_slope_l + plm_slope_r; //it also has some magnitude, but we will only use its sign.
   /*first derivative coefficients*/
   const Vec c0 = fd_l;
   const Vec c1 = b0;
   const Vec c2 = b1 / 2.0;
   const Vec c3 = b2 / 3.0;
   //compute both slopes at inflexion points, at least one of these
   //is with [0..1]. If the root is not in this range, we
   //simplify later if statements by setting it to the plm slope
   //sign
   Vec root1_slope = select(root1 >= 0.0 && root1 <= 1.0,
                             c0  + root1 * ( c1 + root1 * (c2 + root1 * c3 ) ),
                             slope_sign);
   Vec root2_slope = select(root2 >= 0.0 && root2 <= 1.0,
                            c0  + root2 * ( c1 + root2 * (c2 + root2 * c3 ) ),
                            slope_sign);
   Vecb fixInflexion = root1_slope * slope_sign < 0.0 || root2_slope * slope_sign < 0.0;
   if (horizontal_or (fixInflexion) ){ 
      Realv valuesa[VECL];
      Realv fva_l[VECL];
      Realv fva_r[VECL];
      Realv fda_l[VECL];
      Realv fda_r[VECL];
      Realv slope_signa[VECL];
      values[k].store(valuesa);
      fv_l.store(fva_l);
      fd_l.store(fda_l);
      fv_r.store(fva_r);
      fd_r.store(fda_r);
      slope_sign.store(slope_signa);
      
      //todo store and then load data to avoid inserts (is it beneficial...?)
      
//serialized the handling of inflexion points, these do not happen for smooth regions
#pragma ivdep
      for(uint i = 0;i < VECL; i++) {
         if(fixInflexion[i]){
            //need to collapse, at least one inflexion point has wrong
            //sign.
            if(fabs(plm_slope_l[i]) <= fabs(plm_slope_r[i])) {
               //collapse to left edge (eq 21)
               fda_l[i] =  1.0 / 3.0 * ( 10 * valuesa[i] - 2.0 * fva_r[i] - 8.0 * fva_l[i]);
               fda_r[i] =  -10.0 * valuesa[i] + 6.0 * fva_r[i] + 4.0 * fva_l[i];
               //check if PLM slope is consistent (eq 28 & 29)
               if (slope_signa[i] * fda_l[i] < 0) {
                  fda_l[i] =  0;
                  fva_r[i] =  5 * valuesa[i] - 4 * fva_l[i];
                  fda_r[i] =  20 * (valuesa[i] - fva_l[i]);
               }
               else if (slope_signa[i] * fda_r[i] < 0) {
                  fda_r[i] =  0;
                  fva_l[i] =  0.5 * (5 * valuesa[i] - 3 * fva_r[i]);
                  fda_l[i] =  10.0 / 3.0 * (-valuesa[i] + fva_r[i]);
               }
            }
            else {
               //collapse to right edge (eq 21)
               fda_l[i] =  10.0 * valuesa[i] - 6.0 * fva_l[i] - 4.0 * fva_r[i];
               fda_r[i] =  1.0 / 3.0 * ( - 10.0 * valuesa[i] + 2 * fva_l[i] + 8 * fva_r[i]);
               //check if PLM slope is consistent (eq 28 & 29)
               if (slope_signa[i] * fda_l[i] < 0) {
                  fda_l[i] =  0;
                  fva_r[i] =  0.5 * ( 5 * valuesa[i] - 3 * fva_l[i]);
                  fda_r[i] =  10.0 / 3.0 * (valuesa[i] - fva_l[i]);
               }
               else if (slope_signa[i] * fda_r[i] < 0) {
                  fda_r[i] =  0;
                  fva_l[i] =  5 * valuesa[i] - 4 * fva_r[i];
                  fda_l[i] =  20.0 * ( - valuesa[i] + fva_r[i]);
               }
            }
         }
      }      
      fv_l.load(fva_l);
      fd_l.load(fda_l);
      fv_r.load(fva_r);
      fd_r.load(fda_r);
   }
}



// /*
//   PQM reconstruction as published in:
//   White, Laurent, and Alistair Adcroft. “A High-Order Finite Volume Remapping Scheme for Nonuniform Grids: The Piecewise Quartic Method (PQM).” Journal of Computational Physics 227, no. 15 (July 2008): 7394–7422. doi:10.1016/j.jcp.2008.04.026.
// */

inline void compute_pqm_coeff(Vec *values, face_estimate_order order, uint k, Vec a[5], const Realv threshold){
   Vec fv_l; /*left face value*/
   Vec fv_r; /*right face value*/
   Vec fd_l; /*left face derivative*/
   Vec fd_r; /*right face derivative*/
   
   compute_filtered_face_values_derivatives(values, k, order, fv_l, fv_r, fd_l, fd_r, threshold);
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
