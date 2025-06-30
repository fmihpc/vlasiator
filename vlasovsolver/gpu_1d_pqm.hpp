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

#ifndef GPU_1D_PQM_H
#define GPU_1D_PQM_H

#include "../arch/arch_device_api.h"
#include "gpu_slope_limiters.hpp"
#include "gpu_face_estimates.hpp"

/*
  Define functions for Realf instead of Vec
*/

/*make sure quartic polynomial is monotonic*/
static ARCH_DEV inline void filter_pqm_monotonicity(const Realf* __restrict__ values, int k, Realf &fv_l, Realf &fv_r, Realf &fd_l, Realf &fd_r, const int index, const int stride) {
   /*fixed values give to roots clearly outside [0,1], or nonexisting ones*/

   /*second derivative coefficients, eq 23 in white et al.*/
   const Realf b0 =   60.0 * values[k*stride+index] - 24.0 * fv_r - 36.0 * fv_l + 3.0 * (fd_r - 3.0 * fd_l);
   const Realf b1 = -360.0 * values[k*stride+index] + 36.0 * fd_l - 24.0 * fd_r + 168.0 * fv_r + 192.0 * fv_l;
   const Realf b2 =  360.0 * values[k*stride+index] + 30.0 * (fd_r - fd_l) - 180.0 * (fv_l + fv_r);
   /*let's compute sqrt value to be used for computing roots. If we
     take sqrt of negaitve numbers, then we instead set a value that
     will make the root to be +-100 which is well outside range
     of[0,1]. We do not catch FP exceptions, so sqrt(negative) are okish (add
     a max(val_to_sqrt,0) if not*/
   const Realf val_to_sqrt = b1 * b1 - 4 * b0 * b2;
   const Realf sqrt_val = (val_to_sqrt < 0.0) ?
      b1 + 200.0 * b2 :
      sqrt(val_to_sqrt);
   //compute roots. Division is safe with vectorclass (=inf)
   const Realf root1 = (b2 != 0) ? (-b1 + sqrt_val) / (2 * b2) : 0;
   const Realf root2 = (b2 != 0) ? (-b1 - sqrt_val) / (2 * b2) : 0;

   /*PLM slope, MC limiter*/
   const Realf plm_slope_l = 2.0 * (values[k*stride+index] - values[(k-1)*stride+index]);
   const Realf plm_slope_r = 2.0 * (values[(k+1)*stride+index] - values[k*stride+index]);
   const Realf slope_sign = plm_slope_l + plm_slope_r; //it also has some magnitude, but we will only use its sign.
   /*first derivative coefficients*/
   const Realf c0 = fd_l;
   const Realf c1 = b0;
   const Realf c2 = b1 / 2.0;
   const Realf c3 = b2 / 3.0;
   //compute both slopes at inflexion points, at least one of these
   //is with [0..1]. If the root is not in this range, we
   //simplify later if statements by setting it to the plm slope
   //sign
   const Realf root1_slope = (root1 >= 0.0 && root1 <= 1.0) ?
      c0  + root1 * ( c1 + root1 * (c2 + root1 * c3 ) ) :
      slope_sign;
   const Realf root2_slope = (root2 >= 0.0 && root2 <= 1.0) ?
      c0  + root2 * ( c1 + root2 * (c2 + root2 * c3 ) ) :
      slope_sign;
   const bool fixInflexion = root1_slope * slope_sign < 0.0 || root2_slope * slope_sign < 0.0;

   if(fixInflexion) {
      const Realf valuesa = values[k*stride+index];
      Realf fva_l = fv_l;
      Realf fva_r = fv_r;
      Realf fda_l = fd_l;
      Realf fda_r = fd_r;
      const Realf slope_signa = slope_sign;
      //need to collapse, point has wrong sign

      if(fabs(plm_slope_l) <= fabs(plm_slope_r))
      {
         //collapse to left edge (eq 21)
         fda_l =  1.0 / 3.0 * ( 10 * valuesa - 2.0 * fva_r - 8.0 * fva_l);
         fda_r =  -10.0 * valuesa + 6.0 * fva_r + 4.0 * fva_l;
         //check if PLM slope is consistent (eq 28 & 29)
         if (slope_signa * fda_l < 0)
         {
            fda_l =  0;
            fva_r =  5 * valuesa - 4 * fva_l;
            fda_r =  20 * (valuesa - fva_l);
         }
         else if (slope_signa * fda_r < 0)
         {
            fda_r =  0;
            fva_l =  0.5 * (5 * valuesa - 3 * fva_r);
            fda_l =  10.0 / 3.0 * (-valuesa + fva_r);
         }
      }
      else
      {
         //collapse to right edge (eq 21)
         fda_l =  10.0 * valuesa - 6.0 * fva_l - 4.0 * fva_r;
         fda_r =  1.0 / 3.0 * ( - 10.0 * valuesa + 2 * fva_l + 8 * fva_r);
         //check if PLM slope is consistent (eq 28 & 29)
         if (slope_signa * fda_l < 0)
         {
            fda_l =  0;
            fva_r =  0.5 * ( 5 * valuesa - 3 * fva_l);
            fda_r =  10.0 / 3.0 * (valuesa - fva_l);
         }
         else if (slope_signa * fda_r < 0)
         {
            fda_r =  0;
            fva_l =  5 * valuesa - 4 * fva_r;
            fda_l =  20.0 * ( - valuesa + fva_r);
         }
      }
      fv_l = (Realf)fva_l;
      fd_l = (Realf)fda_l;
      fv_r = (Realf)fva_r;
      fd_r = (Realf)fda_r;
   }
}

static ARCH_DEV inline void compute_pqm_coeff(const Realf* __restrict__ values, face_estimate_order order, int k, Realf a[5], const Realf threshold, const int index, const int stride)
{
   Realf fv_l; /*left face value*/
   Realf fv_r; /*right face value*/
   Realf fd_l; /*left face derivative*/
   Realf fd_r; /*right face derivative*/

   compute_filtered_face_values_derivatives(values, k, order, fv_l, fv_r, fd_l, fd_r, threshold, index, stride);
   filter_pqm_monotonicity(values, k, fv_l, fv_r, fd_l, fd_r, index, stride);
   //Fit a second order polynomial for reconstruction see, e.g., White
   //2008 (PQM article) (note additional integration factors built in,
   //contrary to White (2008) eq. 4
   a[0] = fv_l;
   a[1] = fd_l/2.0;
   a[2] =  10.0 * values[k*stride+index] - 4.0 * fv_r - 6.0 * fv_l + 0.5 * (fd_r - 3 * fd_l);
   a[3] = -15.0 * values[k*stride+index]  + 1.5 * fd_l - fd_r + 7.0 * fv_r + 8 * fv_l;
   a[4] =   6.0 * values[k*stride+index] +  0.5 * (fd_r - fd_l) - 3.0 * (fv_l + fv_r);
}


#endif
