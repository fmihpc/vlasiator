/*
This file is part of Vlasiator.
Copyright 2013, 2014 Finnish Meteorological Institute

*/

#ifndef CPU_FACE_ESTIMATES_H
#define CPU_FACE_ESTIMATES_H

#include <iostream>
#include "vec4.h"
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
inline void compute_h8_left_face_value(const Vec4 * const values, uint k, Vec4 &fv_l){   
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
inline void compute_h7_left_face_derivative(const Vec4 * const values, uint k, Vec4 &fd_l){   
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
inline void compute_h6_left_face_value(const Vec4 * const values, uint k, Vec4 &fv_l){   
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
inline void compute_h5_left_face_derivative(const Vec4 * const values, uint k, Vec4 &fd_l){   
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
inline void compute_h5_face_values(const Vec4 * const values, uint k, Vec4 &fv_l, Vec4 &fv_r){   
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
inline void compute_h4_left_face_derivative(const Vec4 * const values, uint k, Vec4 &fd_l){   
  fd_l = 1.0/12.0 * (15.0 * (values[k] - values[k - 1]) - (values[k + 1] - values[k - 2]));
}




/*! 
  Compute left face value based on the explicit h4 estimate.

  Right face value can be obtained as left face value of cell i + 1.
  
  \param values Array with volume averages. It is assumed a large enough stencil is defined around i.
  \param i Index of cell in values for which the left face is computed
  \param fv_l Face value on left face of cell i
*/ 
inline void compute_h4_left_face_value(const Vec4 * const values, uint k, Vec4 &fv_l){   
  /*compute left value*/
  fv_l = 1.0/12.0 * ( - 1.0 * values[k - 2]  
                      + 7.0 * values[k - 1] 
                      + 7.0 * values[k] 
                      - 1.0 * values[k + 1]);
}


/*! 
  Compute left face value based on the explicit h4 estimate.

  Right face value can be obtained as left face value of cell i + 1.
  
  \param values Array with volume averages. It is assumed a large enough stencil is defined around i.
  \param i Index of cell in values for which the left face is computed
  \param fv_l Face value on left face of cell i
*/ 
inline void compute_h3_left_face_derivative(const Vec4 * const values, uint k, Vec4 &fv_l){   
  /*compute left value*/
  fv_l = 1.0/12.0 * (15 * (values[k] - values[k - 1]) - (values[k + 1] - values[k - 2]));
}


/*Filters in section 2.6.1 of white et al. to be used for PQM
  1) Checks for extrema and flattens them
  2) Makes face values bounded
  3) Makes sure face slopes are consistent with PLM slope
*/
inline void compute_filtered_face_values_derivatives(const Vec4 * const values,uint k, face_estimate_order order,
                                                        Vec4 &fv_l, Vec4 &fv_r, Vec4 &fd_l, Vec4 &fd_r){   

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
   
   Vec4 slope_abs,slope_sign;
   slope_limiter(values[k -1], values[k], values[k + 1], slope_abs, slope_sign);
   
   //check for extrema, flatten if it is
   Vec4b is_extrema = (slope_abs == Vec4(0.0));
   if(horizontal_or(is_extrema)) {
      fv_r = select(is_extrema, values[k], fv_r);
      fv_l = select(is_extrema, values[k], fv_l);
      fd_l = select(is_extrema, 0.0 , fd_l);
      fd_r = select(is_extrema, 0.0 , fd_r);
   }

   //Fix left face if needed; boundary value is not bounded or slope is not consistent
   Vec4b filter = (values[k -1] - fv_l) * (fv_l - values[k]) < 0 || slope_sign * fd_l < 0.0;
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
inline void compute_filtered_face_values(const Vec4 * const values,uint k, face_estimate_order order, Vec4 &fv_l, Vec4 &fv_r){   
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
   Vec4 slope_abs,slope_sign;
   slope_limiter(values[k -1], values[k], values[k + 1], slope_abs, slope_sign);
   
   //check for extrema, flatten if it is
   Vec4b is_extrema = (slope_abs == Vec4(0.0));
   if(horizontal_or(is_extrema)) {
      fv_r = select(is_extrema, values[k], fv_r);
      fv_l = select(is_extrema, values[k], fv_l);
   }

   //Fix left face if needed; boundary value is not bounded
   Vec4b filter = (values[k -1] - fv_l) * (fv_l - values[k]) < 0 ;
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

#endif
