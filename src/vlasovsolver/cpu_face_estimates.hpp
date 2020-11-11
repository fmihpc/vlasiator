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

#ifndef CPU_FACE_ESTIMATES_H
#define CPU_FACE_ESTIMATES_H

#include <iostream>
#include "vec.h"
#include "algorithm"
#include "cmath"
#include "cpu_slope_limiters.hpp"

/*enum for setting face value and derivative estimates. Implicit ones
  not supported in the solver, so they are now not listed*/
enum face_estimate_order {h4, h5, h6, h8};

/*! 
  Compute left face value based on the explicit h8 estimate.

  Right face value can be obtained as left face value of cell i + 1.
  
  \param values Array with volume averages. It is assumed a large enough stencil is defined around i.
  \param i Index of cell in values for which the left face is computed
  \param fv_l Face value on left face of cell i
*/ 
inline void compute_h8_left_face_value(const Vec * const values, uint k, Vec &fv_l){   
   fv_l = 1.0/840.0 * (
      - 3.0 * values[k - 4]  
      + 29.0 * values[k - 3]  
      - 139.0 * values[k - 2]  
      + 533.0 * values[k - 1] 
      + 533.0 * values[k] 
      - 139.0 * values[k + 1] 
      + 29.0 * values[k + 2]
      - 3.0 * values[k + 3]);
}


/*! 
  Compute left face derivative based on the explicit h7 estimate.

  Right face derivative can be obtained as left face derivative of cell i + 1.
  
  \param values Array with volume averages. It is assumed a large enough stencil is defined around i.
  \param i Index of cell in values for which the left face derivativeis computed
  \param fd_l Face derivative on left face of cell i
*/ 
inline void compute_h7_left_face_derivative(const Vec * const values, uint k, Vec &fd_l){   
    fd_l = 1.0/5040.0 * (
       + 9.0 * values[k - 4]  
       - 119.0 * values[k - 3]  
       + 889.0 * values[k - 2]  
       - 7175.0 * values[k - 1] 
       + 7175.0 * values[k] 
       - 889.0 * values[k + 1] 
       + 119.0 * values[k + 2]
       - 9.0 * values[k + 3]);
}
   
/*! 
  Compute left face value based on the explicit h6 estimate.

  Right face value can be obtained as left face value of cell i + 1.
  
  \param values Array with volume averages. It is assumed a large enough stencil is defined around i.
  \param i Index of cell in values for which the left face is computed
  \param fv_l Face value on left face of cell i
*/ 
inline void compute_h6_left_face_value(const Vec * const values, uint k, Vec &fv_l){   
  /*compute left value*/
   fv_l = 1.0/60.0 * (values[k - 3]  
                      - 8.0 * values[k - 2]  
                      + 37.0 * values[k - 1] 
                      + 37.0 * values[k ] 
                      - 8.0 * values[k + 1] 
                      + values[k + 2]);
}

/*! 
  Compute left face derivative based on the explicit h5 estimate.

  Right face derivative can be obtained as left face derivative of cell i + 1.
  
  \param values Array with volume averages. It is assumed a large enough stencil is defined around i.
  \param i Index of cell in values for which the left face derivativeis computed
  \param fd_l Face derivative on left face of cell i
*/ 
inline void compute_h5_left_face_derivative(const Vec * const values, uint k, Vec &fd_l){   
  fd_l = 1.0/180.0 * (245 * (values[k] - values[k - 1])  
                     - 25 * (values[k + 1] - values[k - 2]) 
                     + 2 * (values[k + 2] - values[k - 3]));
}




/*! 
  Compute left and right face value based on the explicit h5 estimate.

  \param values Array with volume averages. It is assumed a large enough stencil is defined around i.
  \param i Index of cell in values for which the left face is computed
  \param fv_l Face value on left face of cell i
*/ 
inline void compute_h5_face_values(const Vec * const values, uint k, Vec &fv_l, Vec &fv_r){   
  /*compute left values*/
  fv_l = 1.0/60.0 * (- 3.0 * values[k - 2]  
                     + 27.0 * values[k - 1] 
                     + 47.0 * values[k ] 
                     - 13.0 * values[k + 1] 
                     + 2.0 * values[k + 2]);
  fv_r = 1.0/60.0 * ( 2.0 * values[k - 2] 
                     - 13.0 * values[k - 1] 
                     + 47.0 * values[k]
                     + 27.0 * values[k + 1] 
                     - 3.0 * values[k + 2]);
}
/*! 
  Compute left face derivative based on the explicit h4 estimate.

  Right face derivative can be obtained as left face derivative of cell i + 1.
  
  \param values Array with volume averages. It is assumed a large enough stencil is defined around i.
  \param i Index of cell in values for which the left face derivativeis computed
  \param fd_l Face derivative on left face of cell i
*/ 
inline void compute_h4_left_face_derivative(const Vec * const values, uint k, Vec &fd_l){   
  fd_l = 1.0/12.0 * (15.0 * (values[k] - values[k - 1]) - (values[k + 1] - values[k - 2]));
}




/*! 
  Compute left face value based on the explicit h4 estimate.

  Right face value can be obtained as left face value of cell i + 1.
  
  \param values Array with volume averages. It is assumed a large enough stencil is defined around i.
  \param i Index of cell in values for which the left face is computed
  \param fv_l Face value on left face of cell i
*/ 
inline void compute_h4_left_face_value(const Vec * const values, uint k, Vec &fv_l){   
  /*compute left value*/
  fv_l = 1.0/12.0 * ( - 1.0 * values[k - 2]  
                      + 7.0 * values[k - 1] 
                      + 7.0 * values[k] 
                      - 1.0 * values[k + 1]);
}


/*!

  Compute left face value based on the explicit h4 estimate for nonuniform grid.
  Eqn 45 in White et. al. 2008

  \param u Array with volume averages. It is assumed a large enough stencil is defined around i.
  \param i Index of cell in values for which the left face is computed
  \param fv_l Face value on left face of cell i
  \param h Array with cell widths. Can be in abritrary units since they always cancel. Maybe 1/refinement ratio?
*/
inline void compute_h4_left_face_value_nonuniform(const Vec * const h, const Vec * const u, uint k, Vec &fv_l) {
   
   fv_l = (
           1.0 / ( h[k - 2] + h[k - 1] + h[k] + h[k + 1] )
           * ( ( h[k - 2] + h[k - 1] ) * ( h[k] + h[k + 1] ) / ( h[k - 1] + h[k] )
               * ( u[k - 1] * h[k] + u[k] * h[k - 1] )
               * (1.0 / ( h[k - 2] + h[k - 1] + h[k] ) + 1.0 / ( h[k - 1] + h[k] + h[k + 1] ) )
               + ( h[k] * ( h[k] + h[k + 1] ) ) / ( ( h[k - 2] + h[k - 1] + h[k] ) * (h[k - 2] + h[k - 1] ) )
               * ( u[k - 1] * (h[k - 2] + 2.0 * h[k - 1] ) - ( u[k - 2] * h[k - 1] ) )
               + h[k - 1] * ( h[k - 2] + h[k - 1] ) / ( ( h[k - 1] + h[k] + h[k + 1] ) * ( h[k] + h[k + 1] ) )
               * ( u[k] * ( 2.0 * h[k] + h[k + 1] ) - u[k + 1] * h[k] ) )
           );
}							    

/*! 
  Compute left face value based on the explicit h4 estimate.

  Right face value can be obtained as left face value of cell i + 1.
  
  \param values Array with volume averages. It is assumed a large enough stencil is defined around i.
  \param i Index of cell in values for which the left face is computed
  \param fv_l Face value on left face of cell i
*/ 
inline void compute_h3_left_face_derivative(const Vec * const values, uint k, Vec &fv_l){   
  /*compute left value*/
  fv_l = 1.0/12.0 * (15 * (values[k] - values[k - 1]) - (values[k + 1] - values[k - 2]));
}


/*Filters in section 2.6.1 of white et al. to be used for PQM
  1) Checks for extrema and flattens them
  2) Makes face values bounded
  3) Makes sure face slopes are consistent with PLM slope
*/
inline void compute_filtered_face_values_derivatives(const Vec * const values,uint k, face_estimate_order order,
						     Vec &fv_l, Vec &fv_r, Vec &fd_l, Vec &fd_r,
						     const Realv threshold){

   switch(order){
       case h4:
          compute_h4_left_face_value(values, k, fv_l);
          compute_h4_left_face_value(values, k + 1, fv_r);
          compute_h3_left_face_derivative(values, k, fd_l);
          compute_h3_left_face_derivative(values, k + 1, fd_r);
          break;
       case h5:
          compute_h5_face_values(values, k, fv_l, fv_r);
          compute_h4_left_face_derivative(values, k, fd_l);
          compute_h4_left_face_derivative(values, k + 1, fd_r);
          break;
       default:
       case h6:
          compute_h6_left_face_value(values, k, fv_l);
          compute_h6_left_face_value(values, k + 1, fv_r);   
          compute_h5_left_face_derivative(values, k, fd_l);
          compute_h5_left_face_derivative(values, k + 1, fd_r);
          break;
       case h8:
          compute_h8_left_face_value(values, k, fv_l);
          compute_h8_left_face_value(values, k + 1, fv_r);   
          compute_h7_left_face_derivative(values, k, fd_l);
          compute_h7_left_face_derivative(values, k + 1, fd_r);
          break;
   }
   
   Vec slope_abs,slope_sign;
   // scale values closer to 1 for more accurate slope limiter calculation
   const Realv scale = 1./threshold;
   slope_limiter(values[k -1]*scale, values[k]*scale, values[k + 1]*scale, slope_abs, slope_sign);
   slope_abs = slope_abs*threshold;

   //check for extrema, flatten if it is
   Vecb is_extrema = (slope_abs == Vec(0.0));
   if(horizontal_or(is_extrema)) {
      fv_r = select(is_extrema, values[k], fv_r);
      fv_l = select(is_extrema, values[k], fv_l);
      fd_l = select(is_extrema, 0.0 , fd_l);
      fd_r = select(is_extrema, 0.0 , fd_r);
   }

   //Fix left face if needed; boundary value is not bounded or slope is not consistent
   Vecb filter = (values[k -1] - fv_l) * (fv_l - values[k]) < 0 || slope_sign * fd_l < 0.0;
   if(horizontal_or (filter)) {  
      //Go to linear (PLM) estimates if not ok (this is always ok!)
      fv_l=select(filter, values[k ] - slope_sign * 0.5 * slope_abs, fv_l);
      fd_l=select(filter, slope_sign * slope_abs, fd_l);
   }
   
   //Fix right face if needed; boundary value is not bounded or slope is not consistent 
   filter = (values[k + 1] - fv_r) * (fv_r - values[k]) < 0 || slope_sign * fd_r < 0.0;
   if(horizontal_or (filter)) {  
      //Go to linear (PLM) estimates if not ok (this is always ok!)
      fv_r=select(filter, values[k] + slope_sign * 0.5 * slope_abs, fv_r);
      fd_r=select(filter, slope_sign * slope_abs, fd_r);
   }
}




/*Filters in section 2.6.1 of white et al. to be used for PPM
  1) Checks for extrema and flattens them
  2) Makes face values bounded
  3) Makes sure face slopes are consistent with PLM slope
*/
inline void compute_filtered_face_values(const Vec * const values,uint k, face_estimate_order order, Vec &fv_l, Vec &fv_r, const Realv threshold){
   switch(order){
       case h4:
          compute_h4_left_face_value(values, k, fv_l);
          compute_h4_left_face_value(values, k + 1, fv_r);
          break;
       case h5:
          compute_h5_face_values(values, k, fv_l, fv_r);
          break;
       default:
       case h6:
          compute_h6_left_face_value(values, k, fv_l);
          compute_h6_left_face_value(values, k + 1, fv_r);   
          break;
       case h8:
          compute_h8_left_face_value(values, k, fv_l);
          compute_h8_left_face_value(values, k + 1, fv_r);   
          break;
   }
   Vec slope_abs,slope_sign;
   // scale values closer to 1 for more accurate slope limiter calculation
   const Realv scale = 1./threshold;
   slope_limiter(values[k -1]*scale, values[k]*scale, values[k + 1]*scale, slope_abs, slope_sign);
   slope_abs = slope_abs*threshold;
   
   //check for extrema, flatten if it is
   Vecb is_extrema = (slope_abs == Vec(0.0));
   if(horizontal_or(is_extrema)) {
      fv_r = select(is_extrema, values[k], fv_r);
      fv_l = select(is_extrema, values[k], fv_l);
   }

   //Fix left face if needed; boundary value is not bounded
   Vecb filter = (values[k -1] - fv_l) * (fv_l - values[k]) < 0 ;
   if(horizontal_or (filter)) {  
      //Go to linear (PLM) estimates if not ok (this is always ok!)
      fv_l=select(filter, values[k ] - slope_sign * 0.5 * slope_abs, fv_l);
   }

   //Fix  face if needed; boundary value is not bounded    
   filter = (values[k + 1] - fv_r) * (fv_r - values[k]) < 0;
   if(horizontal_or (filter)) {  
      //Go to linear (PLM) estimates if not ok (this is always ok!)
      fv_r=select(filter, values[k] + slope_sign * 0.5 * slope_abs, fv_r);
   }
}


inline void compute_filtered_face_values_nonuniform(const Vec * const dv, const Vec * const values,uint k, face_estimate_order order, Vec &fv_l, Vec &fv_r, const Realv threshold){
  switch(order){
  case h4:     
     compute_h4_left_face_value_nonuniform(dv, values, k, fv_l);
     compute_h4_left_face_value_nonuniform(dv, values, k + 1, fv_r);
     break;
  // case h5:
  //   compute_h5_face_values(dv, values, k, fv_l, fv_r);
  //   break;
  // case h6:
  //   compute_h6_left_face_value(dv, values, k, fv_l);
  //   compute_h6_left_face_value(dv, values, k + 1, fv_r);   
  //   break;
  // case h8:
  //   compute_h8_left_face_value(dv, values, k, fv_l);
  //   compute_h8_left_face_value(dv, values, k + 1, fv_r);   
  //   break;
  default:
    std::cout << "Order " << order << " has not been implemented (yet)\n";
    break;
  }
   Vec slope_abs,slope_sign;
   // scale values closer to 1 for more accurate slope limiter calculation
   const Realv scale = 1./threshold;
   slope_limiter(values[k -1]*scale, values[k]*scale, values[k + 1]*scale, slope_abs, slope_sign);
   slope_abs = slope_abs*threshold;
   
   //check for extrema, flatten if it is
   Vecb is_extrema = (slope_abs == Vec(0.0));
   if(horizontal_or(is_extrema)) {
      fv_r = select(is_extrema, values[k], fv_r);
      fv_l = select(is_extrema, values[k], fv_l);
   }

   //Fix left face if needed; boundary value is not bounded
   Vecb filter = (values[k -1] - fv_l) * (fv_l - values[k]) < 0 ;
   if(horizontal_or (filter)) {  
      //Go to linear (PLM) estimates if not ok (this is always ok!)
      fv_l=select(filter, values[k ] - slope_sign * 0.5 * slope_abs, fv_l);
   }

   //Fix  face if needed; boundary value is not bounded    
   filter = (values[k + 1] - fv_r) * (fv_r - values[k]) < 0;
   if(horizontal_or (filter)) {  
      //Go to linear (PLM) estimates if not ok (this is always ok!)
      fv_r=select(filter, values[k] + slope_sign * 0.5 * slope_abs, fv_r);
   }
}

inline Vec get_D2aLim(const Vec * h, const Vec * values, uint k, const Vec C, Vec & fv) {

  // Colella & Sekora, eq. 18
  Vec invh2 = 1.0 / (h[k] * h[k]);
  Vec d2a =  invh2 * 3.0 * (values[k]     - 2.0 * fv            + values[k + 1]);
  Vec d2aL = invh2       * (values[k - 1] - 2.0 * values[k]     + values[k + 1]);
  Vec d2aR = invh2       * (values[k]     - 2.0 * values[k + 1] + values[k + 2]);
  Vec d2aLim;
  if ( (horizontal_or(d2a * d2aL >= 0)) && (horizontal_or(d2a * d2aR >= 0)) &&
       horizontal_and(d2a != 0)) {
    d2aLim = d2a / abs(d2a) * min(abs(d2a),min(C*abs(d2aL),C*abs(d2aR)));
  } else {
    d2aLim = 0.0;
  }

  return d2aLim;
  
}

inline pair<Vec,Vec> constrain_face_values(const Vec * h,const Vec * values,uint k,Vec & fv_l, Vec & fv_r) {

  const Vec C = 1.25;
  Vec invh2 = 1.0 / (h[k] * h[k]);
  
  // Colella & Sekora, eq 19
  Vec p_face = 0.5 * (values[k] + values[k + 1])
    - h[k] * h[k] / 3.0 * get_D2aLim(h,values,k  ,C,fv_r);
  Vec m_face = 0.5 * (values[k-1] + values[k])
    - h[k-1] * h[k-1] / 3.0 * get_D2aLim(h,values,k-1,C,fv_l);

  // Colella & Sekora, eq 21
  Vec d2a = -2.0 * invh2 * 6.0 * (values[k] - 3.0 * (m_face + p_face)); // a6,j from eq. 7 
  Vec d2aC = invh2 * (values[k - 1] - 2.0 * values[k    ] + values[k + 1]);
  // Note: Corrected the index of 2nd term in d2aL to k - 1.
  //       In the paper it is k but that is almost certainly an error.
  Vec d2aL = invh2 * (values[k - 2] - 2.0 * values[k - 1] + values[k    ]);
  Vec d2aR = invh2 * (values[k    ] - 2.0 * values[k + 1] + values[k + 2]);
  Vec d2aLim;

  // Colella & Sekora, eq 22
  if ( (horizontal_or(d2a * d2aL >= 0)) && (horizontal_or(d2a * d2aR >= 0)) &&
       (horizontal_or(d2a * d2aC >= 0)) && horizontal_and(d2a != 0)) {
    
    d2aLim = d2a / abs(d2a) * min(C * abs(d2aL), min(C * abs(d2aR), min(C * abs(d2aC), abs(d2a))));
  } else {
    d2aLim = 0.0;
    if( horizontal_or(d2a == 0.0)) {
      // Set a non-zero value for the denominator in eq. 23.
      // According to the paper the ratio d2aLim/d2a should be 0
      d2a = 1.0;
    }
  }
  
  // Colella & Sekora, eq 23
  Vec p_face_interpolant = values[k] + (p_face - values[k]) * d2aLim / d2a;
  Vec m_face_interpolant = values[k] + (m_face - values[k]) * d2aLim / d2a;
  
  pair<Vec,Vec> faceInterpolants;
  faceInterpolants = make_pair(m_face_interpolant,p_face_interpolant);
  
  // if(horizontal_and(values[k] > 1.0))
  //   std::cout << "k, a-, a+, <a>: " << k << ", " << p_face_interpolant[0] << ", " << m_face_interpolant[0] << ", " << values[k][0] << endl;
  
  return faceInterpolants;
}

inline void compute_filtered_face_values_nonuniform_conserving(const Vec * const dv, const Vec * const values,uint k, face_estimate_order order, Vec &fv_l, Vec &fv_r, const Realv threshold){
  switch(order){
  case h4:
    compute_h4_left_face_value_nonuniform(dv, values, k, fv_l);
    compute_h4_left_face_value_nonuniform(dv, values, k + 1, fv_r);
    break;
  // case h5:
  //   compute_h5_face_values(dv, values, k, fv_l, fv_r);
  //   break;
  // case h6:
  //   compute_h6_left_face_value(dv, values, k, fv_l);
  //   compute_h6_left_face_value(dv, values, k + 1, fv_r);   
  //   break;
  // case h8:
  //   compute_h8_left_face_value(dv, values, k, fv_l);
  //   compute_h8_left_face_value(dv, values, k + 1, fv_r);   
  //   break;
  default:
    std::cout << "Order " << order << " has not been implemented (yet)\n";
    break;
  }

  Vec slope_abs,slope_sign;
  // scale values closer to 1 for more accurate slope limiter calculation
  const Realv scale = 1./threshold;
  slope_limiter(values[k -1]*scale, values[k]*scale, values[k + 1]*scale, slope_abs, slope_sign);
  slope_abs = slope_abs*threshold;
  
  //check for extrema  
  Vecb is_extrema = (slope_abs == Vec(0.0));
  Vecb filter_l = (values[k - 1] - fv_l) * (fv_l - values[k]) < 0 ;
  Vecb filter_r = (values[k + 1] - fv_r) * (fv_r - values[k]) < 0;
  //  if(horizontal_or(is_extrema) || horizontal_or(filter_l) || horizontal_or(filter_r)) {
  // Colella & Sekora, eq. 20
  if(horizontal_or((fv_r - values[k]) * (values[k] - fv_l) <= Vec(0.0))
     && horizontal_or((values[k - 1] - values[k]) * (values[k] - values[k + 1]) <= Vec(0.0))) {
    auto faces = constrain_face_values(dv, values, k, fv_l, fv_r);

    fv_l = faces.first;
    fv_r = faces.second;

  //    fv_r = select(is_extrema, values[k], fv_r);
  //    fv_l = select(is_extrema, values[k], fv_l);
  } else {

     //Fix left face if needed; boundary value is not bounded
     Vecb filter = (values[k -1] - fv_l) * (fv_l - values[k]) < 0 ;
     if(horizontal_or (filter)) {  
        //Go to linear (PLM) estimates if not ok (this is always ok!)
        fv_l=select(filter, values[k ] - slope_sign * 0.5 * slope_abs, fv_l);
     }
     
     //Fix  face if needed; boundary value is not bounded    
     filter = (values[k + 1] - fv_r) * (fv_r - values[k]) < 0;
     if(horizontal_or (filter)) {  
        //Go to linear (PLM) estimates if not ok (this is always ok!)
        fv_r=select(filter, values[k] + slope_sign * 0.5 * slope_abs, fv_r);
     }
  }
  
  // //Fix left face if needed; boundary value is not bounded
  // Vecb filter = (values[k -1] - fv_l) * (fv_l - values[k]) < 0 ;
  // if(horizontal_or (filter)) {  
  //   //Go to linear (PLM) estimates if not ok (this is always ok!)
  //   fv_l=select(filter, values[k ] - slope_sign * 0.5 * slope_abs, fv_l);
  // }
  
  // //Fix  face if needed; boundary value is not bounded    
  // filter = (values[k + 1] - fv_r) * (fv_r - values[k]) < 0;
  // if(horizontal_or (filter)) {  
  //   //Go to linear (PLM) estimates if not ok (this is always ok!)
  //   fv_r=select(filter, values[k] + slope_sign * 0.5 * slope_abs, fv_r);
  // }

  
}


#endif
