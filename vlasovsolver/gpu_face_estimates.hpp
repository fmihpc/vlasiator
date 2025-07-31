/*
 * This file is part of Vlasiator.
 * Copyright 2010-2021 Finnish Meteorological Institute
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

#ifndef GPU_FACE_ESTIMATES_H
#define GPU_FACE_ESTIMATES_H

#include "gpu_slope_limiters.hpp"
#include "../definitions.h"

#include "../arch/arch_device_api.h"

/*enum for setting face value and derivative estimates. Implicit ones
  not supported in the solver, so they are now not listed*/
enum face_estimate_order {h4, h5, h6, h8};
/**
   Define functions for Realf instead of Vec
*/


/*!
  Compute left face value based on the explicit h8 estimate.

  Right face value can be obtained as left face value of cell i + 1.

  \param values Array with volume averages. It is assumed a large enough stencil is defined around i.
  \param i Index of cell in values for which the left face is computed
  \param fv_l Face value on left face of cell i
*/
ARCH_DEV inline void compute_h8_left_face_value(const Realf* const values, int k, Realf &fv_l, const int index, const int stride)
{
   fv_l = (Realf)(1.0/840.0) * (
      (Realf)(- 3.0) * values[(k-4)*stride+index]
      + (Realf)(29.0) * values[(k-3)*stride+index]
      - (Realf)(139.0) * values[(k-2)*stride+index]
      + (Realf)(533.0) * values[(k-1)*stride+index]
      + (Realf)(533.0) * values[k*stride+index]
      - (Realf)(139.0) * values[(k+1)*stride+index]
      + (Realf)(29.0) * values[(k+2)*stride+index]
      - (Realf)(3.0) * values[(k+3)*stride+index]);
}
ARCH_DEV inline void compute_h7_left_face_derivative(const Realf* const values, int k, Realf &fd_l, const int index, const int stride){
   fd_l = (Realf)(1.0/5040.0) * (
      (Realf)(9.0) * values[(k-4)*stride+index]
      - (Realf)(119.0) * values[(k-3)*stride+index]
      + (Realf)(889.0) * values[(k-2)*stride+index]
      - (Realf)(7175.0) * values[(k-1)*stride+index]
      + (Realf)(7175.0) * values[k*stride+index]
      - (Realf)(889.0) * values[(k+1)*stride+index]
      + (Realf)(119.0) * values[(k+2)*stride+index]
      - (Realf)(9.0) * values[(k+3)*stride+index]);
}
ARCH_DEV inline void compute_h6_left_face_value(const Realf* const values, int k, Realf &fv_l, const int index, const int stride)
{
   //compute left value
   fv_l = (Realf)(1.0/60.0) * (values[(k-3)*stride+index]
                      - (Realf)(8.0) * values[(k-2)*stride+index]
                      + (Realf)(37.0) * values[(k-1)*stride+index]
                      + (Realf)(37.0) * values[k*stride+index]
                      - (Realf)(8.0) * values[(k+1)*stride+index]
                      + values[(k+2)*stride+index]);
}
ARCH_DEV inline void compute_h5_left_face_derivative(const Realf* const values, int k, Realf &fd_l, const int index, const int stride)
{
   fd_l = (Realf)(1.0/180.0) * ((Realf)(245.0) * (values[k*stride+index] - values[(k-1)*stride+index])
                       - (Realf)(25.0) * (values[(k+1)*stride+index] - values[(k-2)*stride+index])
                       + (Realf)(2.0) * (values[(k+2)*stride+index] - values[(k-3)*stride+index]));
}
ARCH_DEV inline void compute_h5_face_values(const Realf* const values, int k, Realf &fv_l, Realf &fv_r, const int index, const int stride)
{
   //compute left values
   fv_l = (Realf)(1.0/60.0) * ((Realf)(- 3.0) * values[(k-2)*stride+index]
                      + (Realf)(27.0) * values[(k-1)*stride+index]
                      + (Realf)(47.0) * values[k*stride+index]
                      - (Realf)(13.0) * values[(k+1)*stride+index]
                      + (Realf)(2.0) * values[(k+2)*stride+index]);
   fv_r = (Realf)(1.0/60.0) * ( (Realf)(2.0) * values[(k-2)*stride+index]
                       - (Realf)(13.0) * values[(k-1)*stride+index]
                       + (Realf)(47.0) * values[k*stride+index]
                       + (Realf)(27.0) * values[(k+1)*stride+index]
                       - (Realf)(3.0) * values[(k+2)*stride+index]);
}
ARCH_DEV inline void compute_h4_left_face_derivative(const Realf* const values, int k, Realf &fd_l, const int index, const int stride)
{
   fd_l = (Realf)(1.0/12.0) * ((Realf)(15.0) * (values[k*stride+index] - values[(k-1)*stride+index]) - (values[(k+1)*stride+index] - values[(k-2)*stride+index]));
}
ARCH_DEV inline void compute_h4_left_face_value(const Realf* const values, int k, Realf &fv_l, const int index, const int stride)
{
   //compute left value
   fv_l = (Realf)(1.0/12.0) * ( (Realf)(- 1.0) * values[(k-2)*stride+index]
                       + (Realf)(7.0) * values[(k-1)*stride+index]
                       + (Realf)(7.0) * values[k*stride+index]
                       - (Realf)(1.0) * values[(k+1)*stride+index]);
}

// h is bin width (dv or dx)
// u is values
ARCH_DEV inline void compute_h4_left_face_value_nonuniform(const Realf* const h, const Realf* const u, int k, Realf &fv_l, const int index, const int stride) {
   const Realf hkMinus2 = h[k-2];
   const Realf hkMinus1 = h[k-1];
   const Realf hk = h[k];
   const Realf hkPlus1 = h[k+1];
   fv_l = (
      (Realf)(1.0) / ( hkMinus2 + hkMinus1 + hk + hkPlus1 )
      * ( ( hkMinus2 + hkMinus1 ) * ( hk + hkPlus1 ) / ( hkMinus1 + hk )
          * ( u[(k-1)*stride+index] * hk + u[k*stride+index] * hkMinus1 )
          * ((Realf)(1.0) / ( hkMinus2 + hkMinus1 + hk ) + (Realf)(1.0) / ( hkMinus1 + hk + hkPlus1 ) )
          + ( hk * ( hk + hkPlus1 ) ) / ( ( hkMinus2 + hkMinus1 + hk ) * (hkMinus2 + hkMinus1 ) )
          * ( u[(k-1)*stride+index] * (hkMinus2 + (Realf)(2.0) * hkMinus1 ) - ( u[(k-2)*stride+index] * hkMinus1 ) )
          + hkMinus1 * ( hkMinus2 + hkMinus1 ) / ( ( hkMinus1 + hk + hkPlus1 ) * ( hk + hkPlus1 ) )
          * ( u[k*stride+index] * ( (Realf)(2.0) * hk + hkPlus1 ) - u[(k+1)*stride+index] * hk ) )
      );
}




ARCH_DEV inline void compute_h3_left_face_derivative(const Realf* const values, int k, Realf &fv_l, const int index, const int stride)
{
   /*compute left value*/
   fv_l = (Realf)(1.0/12.0) * ((Realf)(15.0) * (values[k*stride+index] - values[(k-1)*stride+index]) - (values[(k+1)*stride+index] - values[(k-2)*stride+index]));
}

ARCH_DEV inline void compute_filtered_face_values_derivatives(const Realf* const values, int k, face_estimate_order order, Realf &fv_l, Realf &fv_r, Realf &fd_l, Realf &fd_r, const Realf threshold, const int index,  const int stride)
{
   switch(order)
   {
      case h4:
         compute_h4_left_face_value(values, k, fv_l, index, stride);
         compute_h4_left_face_value(values, k + 1, fv_r, index, stride);
         compute_h3_left_face_derivative(values, k, fd_l, index, stride);
         compute_h3_left_face_derivative(values, k + 1, fd_r, index, stride);
         break;
      case h5:
         compute_h5_face_values(values, k, fv_l, fv_r, index, stride);
         compute_h4_left_face_derivative(values, k, fd_l, index, stride);
         compute_h4_left_face_derivative(values, k + 1, fd_r, index, stride);
         break;
      default:
      case h6:
         compute_h6_left_face_value(values, k, fv_l, index, stride);
         compute_h6_left_face_value(values, k + 1, fv_r, index, stride);
         compute_h5_left_face_derivative(values, k, fd_l, index, stride);
         compute_h5_left_face_derivative(values, k + 1, fd_r, index, stride);
         break;
      case h8:
         compute_h8_left_face_value(values, k, fv_l, index, stride);
         compute_h8_left_face_value(values, k + 1, fv_r, index, stride);
         compute_h7_left_face_derivative(values, k, fd_l, index, stride);
         compute_h7_left_face_derivative(values, k + 1, fd_r, index, stride);
         break;
   }
   Realf slope_abs,slope_sign;
   // scale values closer to 1 for more accurate slope limiter calculation
   const Realf scale = (Realf)(1.0)/threshold;
   slope_limiter(values[(k-1)*stride+index]*scale, values[k*stride+index]*scale, values[(k+1)*stride+index]*scale, slope_abs, slope_sign);
   slope_abs = slope_abs*threshold;
   //check for extrema, flatten if it is
   bool is_extrema = (slope_abs == (Realf)(0.0));
   if (is_extrema) {
      fv_r = (is_extrema) ? values[k*stride+index] : fv_r;
      fv_l = (is_extrema) ? values[k*stride+index] : fv_l;
      fd_l = (is_extrema) ? (Realf)(0.0) : fd_l;
      fd_r = (is_extrema) ? (Realf)(0.0) : fd_r;
   }
   //Fix left face if needed; boundary value is not bounded or slope is not consistent
   bool filter = (values[(k-1)*stride+index] - fv_l) * (fv_l - values[k*stride+index]) < (Realf)(0.0) || slope_sign * fd_l < (Realf)(0.0);
   if (filter) {
      //Go to linear (PLM) estimates if not ok (this is always ok!)
      fv_l= (filter) ? values[k*stride+index] - slope_sign * (Realf)(0.5) * slope_abs : fv_l;
      fd_l= (filter) ? slope_sign * slope_abs : fd_l;
   }
   //Fix right face if needed; boundary value is not bounded or slope is not consistent
   filter = (values[(k+1)*stride+index] - fv_r) * (fv_r - values[k*stride+index]) < (Realf)(0.0) || slope_sign * fd_r < (Realf)(0.0);
   if (filter) {
      //Go to linear (PLM) estimates if not ok (this is always ok!)
      fv_r= (filter) ? values[k*stride+index] + slope_sign * (Realf)(0.5) * slope_abs : fv_r;
      fd_r= (filter) ? slope_sign * slope_abs : fd_r;
   }
}

/*Filters in section 2.6.1 of white et al. to be used for PPM
  1) Checks for extrema and flattens them
  2) Makes face values bounded
  3) Makes sure face slopes are consistent with PLM slope
*/
ARCH_DEV inline void compute_filtered_face_values(const Realf* const values, int k, face_estimate_order order, Realf &fv_l, Realf &fv_r, const Realf threshold, const int index, const int stride)
{
   switch(order)
   {
      case h4:
         compute_h4_left_face_value(values, k, fv_l, index, stride);
         compute_h4_left_face_value(values, k + 1, fv_r, index, stride);
         break;
      case h5:
         compute_h5_face_values(values, k, fv_l, fv_r, index, stride);
         break;
      default:
      case h6:
         compute_h6_left_face_value(values, k, fv_l, index, stride);
         compute_h6_left_face_value(values, k + 1, fv_r, index, stride);
         break;
      case h8:
         compute_h8_left_face_value(values, k, fv_l, index, stride);
         compute_h8_left_face_value(values, k + 1, fv_r, index, stride);
         break;
   }
   Realf slope_abs, slope_sign;
   // scale values closer to 1 for more accurate slope limiter calculation
   const Realf scale = (Realf)(1.0)/threshold;
   slope_limiter(values[(k-1)*stride+index]*scale, values[k*stride+index]*scale, values[(k+1)*stride+index]*scale, slope_abs, slope_sign);
   slope_abs = slope_abs*threshold;

   //check for extrema, flatten if it is
   bool is_extrema = (slope_abs == (Realf)(0.0));
   if (is_extrema) {
      fv_r = (is_extrema) ? values[k*stride+index] : fv_r;
      fv_l = (is_extrema) ? values[k*stride+index] : fv_l;
   }
   //Fix left face if needed; boundary value is not bounded
   bool filter = (values[(k-1)*stride+index] - fv_l) * (fv_l - values[k*stride+index]) < (Realf)(0.0) ;
   if (filter) {
      //Go to linear (PLM) estimates if not ok (this is always ok!)
      fv_l = (filter) ? values[k*stride+index] - slope_sign * (Realf)(0.5) * slope_abs : fv_l;
   }
   //Fix  face if needed; boundary value is not bounded
   filter = (values[(k+1)*stride+index] - fv_r) * (fv_r - values[k*stride+index]) < (Realf)(0.0);
   if (filter) {
      //Go to linear (PLM) estimates if not ok (this is always ok!)
      fv_r = (filter) ? values[k*stride+index] + slope_sign * (Realf)(0.5) * slope_abs : fv_r;
   }
}

ARCH_DEV inline void compute_filtered_face_values_nonuniform(const Realf * const dv, const Realf* const values,int k, face_estimate_order order, Realf &fv_l, Realf &fv_r, const Realf threshold, const int index, const int stride){
   switch(order){
      case h4:
         compute_h4_left_face_value_nonuniform(dv, values, k, fv_l, index, stride);
         compute_h4_left_face_value_nonuniform(dv, values, k + 1, fv_r, index, stride);
         break;
         // case h5:
         //   compute_h5_face_values(dv, values, k, fv_l, fv_r, index, stride);
         //   break;
         // case h6:
         //   compute_h6_left_face_value(dv, values, k, fv_l, index, stride);
         //   compute_h6_left_face_value(dv, values, k + 1, fv_r, index, stride);
         //   break;
         // case h8:
         //   compute_h8_left_face_value(dv, values, k, fv_l, index, stride);
         //   compute_h8_left_face_value(dv, values, k + 1, fv_r, index, stride);
         //   break;
      default:
         printf("Order %d has not been implemented (yet)\n",order);
         break;
   }
   Realf slope_abs,slope_sign;
   if (threshold>(Realf)(0.0)) {
      // scale values closer to 1 for more accurate slope limiter calculation
      const Realf scale = (Realf)(1.0)/threshold;
      slope_limiter(values[(k-1)*stride+index]*scale, values[k*stride+index]*scale, values[(k+1)*stride+index]*scale, slope_abs, slope_sign);
      slope_abs = slope_abs*threshold;
   } else {
      slope_limiter(values[(k-1)*stride+index], values[k*stride+index], values[(k+1)*stride+index], slope_abs, slope_sign);
   }

   //check for extrema, flatten if it is
   if (slope_abs == (Realf)(0.0)) {
      fv_r = values[k*stride+index];
      fv_l = values[k*stride+index];
   }

   //Fix left face if needed; boundary value is not bounded
   if ((values[(k-1)*stride+index] - fv_l) * (fv_l - values[k*stride+index]) < (Realf)(0.0)) {
      //Go to linear (PLM) estimates if not ok (this is always ok!)
      fv_l=values[k*stride+index] - slope_sign * (Realf)(0.5) * slope_abs;
   }

   //Fix  face if needed; boundary value is not bounded
   if ((values[(k+1)*stride+index] - fv_r) * (fv_r - values[k*stride+index]) < (Realf)(0.0)) {
      //Go to linear (PLM) estimates if not ok (this is always ok!)
      fv_r=values[k*stride+index] + slope_sign * (Realf)(0.5) * slope_abs;
   }
}

ARCH_DEV inline Realf get_D2aLim(const Realf* h, const Realf* values, int k, const Realf C, Realf & fv, const int index, const int stride) {

   // Colella & Sekora, eq. 18
   Realf invh2 = 1.0 / (h[k] * h[k]);
   Realf d2a =  invh2 * 3.0 * (values[k*stride    +index] - 2.0 * fv                         + values[(k+1)*stride+index]);
   Realf d2aL = invh2       * (values[(k-1)*stride+index] - 2.0 * values[k*stride    +index] + values[(k+1)*stride+index]);
   Realf d2aR = invh2       * (values[k*stride    +index] - 2.0 * values[(k+1)*stride+index] + values[(k+2)*stride+index]);
   Realf d2aLim;
   if ( (d2a * d2aL >= 0) && (d2a * d2aR >= 0) && (d2a != 0) ) {
      d2aLim = d2a / abs(d2a) * min(abs(d2a),min(C*abs(d2aL),C*abs(d2aR)));
   } else {
      d2aLim = 0.0;
   }
   return d2aLim;
}

ARCH_DEV inline void constrain_face_values(const Realf* h, const Realf* values,int k,Realf & fv_l, Realf & fv_r, const int index, const int stride) {

   const Realf C = 1.25;
   Realf invh2 = 1.0 / (h[k] * h[k]);

   // Colella & Sekora, eq 19
   Realf p_face = 0.5 * (values[k*stride+index] + values[(k+1)*stride+index])
      - h[k] * h[k] / 3.0 * get_D2aLim(h,values,k  ,C,fv_r, index, stride);
   Realf m_face = 0.5 * (values[(k-1)*stride+index] + values[k*stride+index])
      - h[k-1] * h[k-1] / 3.0 * get_D2aLim(h,values,k-1,C,fv_l, index, stride);

   // Colella & Sekora, eq 21
   Realf d2a = -2.0 * invh2 * 6.0 * (values[k*stride+index] - 3.0 * (m_face + p_face)); // a6,j from eq. 7
   Realf d2aC = invh2 * (values[(k-1)*stride+index] - 2.0 * values[k*stride+index] + values[(k+1)*stride+index]);
   // Note: Corrected the index of 2nd term in d2aL to k - 1.
   //       In the paper it is k but that is almost certainly an error.
   Realf d2aL = invh2 * (values[(k-2)*stride+index] - 2.0 * values[(k-1)*stride+index] + values[k*stride+index]);
   Realf d2aR = invh2 * (values[k*stride+index] - 2.0 * values[(k+1)*stride+index] + values[(k+2)*stride+index]);
   Realf d2aLim;

   // Colella & Sekora, eq 22
   if ( (d2a * d2aL >= 0) && (d2a * d2aR >= 0) &&
        (d2a * d2aC >= 0) && (d2a != 0) ) {

      d2aLim = d2a / abs(d2a) * min(C * abs(d2aL), min(C * abs(d2aR), min(C * abs(d2aC), abs(d2a))));
   } else {
      d2aLim = 0.0;
      if (d2a == 0.0) {
         // Set a non-zero value for the denominator in eq. 23.
         // According to the paper the ratio d2aLim/d2a should be 0
         d2a = 1.0;
      }
   }

   // Colella & Sekora, eq 23
   fv_r = values[k*stride+index] + (p_face - values[k*stride+index]) * d2aLim / d2a;
   fv_l = values[k*stride+index] + (m_face - values[k*stride+index]) * d2aLim / d2a;

}

ARCH_DEV inline void compute_filtered_face_values_nonuniform_conserving(const Realf * const dv, const Realf* const values,int k, face_estimate_order order, Realf &fv_l, Realf &fv_r, const Realf threshold, const int index, const int stride){
   switch(order){
      case h4:
         compute_h4_left_face_value_nonuniform(dv, values, k, fv_l, index, stride);
         compute_h4_left_face_value_nonuniform(dv, values, k + 1, fv_r, index, stride);
         break;
         // case h5:
         //   compute_h5_face_values(dv, values, k, fv_l, fv_r, index, stride);
         //   break;
         // case h6:
         //   compute_h6_left_face_value(dv, values, k, fv_l, index, stride);
         //   compute_h6_left_face_value(dv, values, k + 1, fv_r, index, stride);
         //   break;
         // case h8:
         //   compute_h8_left_face_value(dv, values, k, fv_l, index, stride);
         //   compute_h8_left_face_value(dv, values, k + 1, fv_r, index, stride);
         //   break;
      default:
         printf("Order %d has not been implemented (yet)\n",order);
         break;
   }

   Realf slope_abs,slope_sign;
   if (threshold>0) {
      // scale values closer to 1 for more accurate slope limiter calculation
      const Realf scale = 1./threshold;
      slope_limiter(values[(k-1)*stride+index]*scale, values[k*stride+index]*scale, values[(k+1)*stride+index]*scale, slope_abs, slope_sign);
      slope_abs = slope_abs*threshold;
   } else {
      slope_limiter(values[(k-1)*stride+index], values[k*stride+index], values[(k+1)*stride+index], slope_abs, slope_sign);
   }

   //check for extrema
   //bool is_extrema = (slope_abs == 0.0);
   // bool filter_l = (values[(k-1)*stride+index] - fv_l) * (fv_l - values[k*stride+index]) < 0 ;
   // bool filter_r = (values[(k+1)*stride+index] - fv_r) * (fv_r - values[k*stride+index]) < 0;
   //  if(horizontal_or(is_extrema) || horizontal_or(filter_l) || horizontal_or(filter_r)) {
   // Colella & Sekora, eq. 20
   if (((fv_r - values[k*stride+index]) * (values[k*stride+index] - fv_l) <= 0.0)
       && ((values[(k-1)*stride+index] - values[k*stride+index]) * (values[k*stride+index] - values[(k+1)*stride+index]) <= 0.0)) {
      constrain_face_values(dv, values, k, fv_l, fv_r, index, stride);

      //    fv_r = select(is_extrema, values[k], fv_r);
      //    fv_l = select(is_extrema, values[k], fv_l);
   } else {

      //Fix left face if needed; boundary value is not bounded
      bool filter = (values[(k-1)*stride+index] - fv_l) * (fv_l - values[k*stride+index]) < 0 ;
      if (filter) {
         //Go to linear (PLM) estimates if not ok (this is always ok!)
         fv_l=values[k*stride+index] - slope_sign * 0.5 * slope_abs;
      }

      //Fix  face if needed; boundary value is not bounded
      filter = (values[(k+1)*stride+index] - fv_r) * (fv_r - values[k*stride+index]) < 0;
      if (filter) {
         //Go to linear (PLM) estimates if not ok (this is always ok!)
         fv_r=values[k*stride+index] + slope_sign * 0.5 * slope_abs;
      }
   }
}

















// ARCH_DEV inline void compute_h4_left_face_value_nonuniform(const Realf * const h, const Realf * const u, int k, Realf &fv_l, const int index, const int stride) {
//    fv_l = (
//       1.0 / ( h[k - 2] + h[k - 1] + h[k] + h[k + 1] )
//       * ( ( h[k - 2] + h[k - 1] ) * ( h[k] + h[k + 1] ) / ( h[k - 1] + h[k] )
//           * ( u[(k-1)*stride+index] * h[k] + u[k*stride+index] * h[k - 1] )
//           * (1.0 / ( h[k - 2] + h[k - 1] + h[k] ) + 1.0 / ( h[k - 1] + h[k] + h[k + 1] ) )
//           + ( h[k] * ( h[k] + h[k + 1] ) ) / ( ( h[k - 2] + h[k - 1] + h[k] ) * (h[k - 2] + h[k - 1] ) )
//           * ( u[(k-1)*stride+index] * (h[k - 2] + 2.0 * h[k - 1] ) - ( u[(k-2)*stride+index] * h[k - 1] ) )
//           + h[k - 1] * ( h[k - 2] + h[k - 1] ) / ( ( h[k - 1] + h[k] + h[k + 1] ) * ( h[k] + h[k + 1] ) )
//           * ( u[k*stride+index] * ( 2.0 * h[k] + h[k + 1] ) - u[(k+1)*stride+index] * h[k] ) )
//       );
// }


// ARCH_DEV inline void compute_filtered_face_values_nonuniform(const Realf * const dv, const Realf * const values,int k, face_estimate_order order, Realf &fv_l, Realf &fv_r, const Realf threshold, const int index, const int stride){
//    switch(order){
//       case h4:
//          compute_h4_left_face_value_nonuniform(dv, values, k, fv_l, index, stride);
//          compute_h4_left_face_value_nonuniform(dv, values, k + 1, fv_r, index, stride);
//          break;
//          // case h5:
//          //   compute_h5_face_values(dv, values, k, fv_l, fv_r, index, stride);
//          //   break;
//          // case h6:
//          //   compute_h6_left_face_value(dv, values, k, fv_l, index, stride);
//          //   compute_h6_left_face_value(dv, values, k + 1, fv_r, index, stride);
//          //   break;
//          // case h8:
//          //   compute_h8_left_face_value(dv, values, k, fv_l, index, stride);
//          //   compute_h8_left_face_value(dv, values, k + 1, fv_r, index, stride);
//          //   break;
//       default:
//          printf("Order %d has not been implemented (yet)\n",order);
//          break;
//    }
//    Realf slope_abs,slope_sign;
//    if (threshold>0) {
//       // scale values closer to 1 for more accurate slope limiter calculation
//       const Realf scale = 1./threshold;
//       slope_limiter(values[(k-1)*stride+index]*scale, values[k*stride+index]*scale, values[(k+1)*stride+index]*scale, slope_abs, slope_sign);
//       slope_abs = slope_abs*threshold;
//    } else {
//       slope_limiter(values[(k-1)*stride+index], values[k*stride+index], values[(k+1)*stride+index], slope_abs, slope_sign);
//    }

//    //check for extrema, flatten if it is
//    if (slope_abs == 0) {
//       fv_r = values[k*stride+index];
//       fv_l = values[k*stride+index];
//    }

//    //Fix left face if needed; boundary value is not bounded
//    if ((values[(k-1)*stride+index] - fv_l) * (fv_l - values[k*stride+index]) < 0) {
//       //Go to linear (PLM) estimates if not ok (this is always ok!)
//       fv_l=values[k*stride+index] - slope_sign * 0.5 * slope_abs;
//    }

//    //Fix  face if needed; boundary value is not bounded
//    if ((values[(k+1)*stride+index] - fv_r) * (fv_r - values[k*stride+index]) < 0) {
//       //Go to linear (PLM) estimates if not ok (this is always ok!)
//       fv_r=values[k*stride+index] + slope_sign * 0.5 * slope_abs;
//    }
// }



#endif
