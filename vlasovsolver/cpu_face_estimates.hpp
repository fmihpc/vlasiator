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
enum face_estimate_order {lin, h4, h6, h8};

/*FIXME & NOTE. MAX_IH_CELLS should be in practice (MAX_BLOCKS_PER_DIM * WID), so 400
  supports up to 100 blocks. This can be too little at some point of
  time and no checks are done for bounds...*/
#define MAX_IH_CELLS 400



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
inline void compute_filtered_h8_face_values_derivatives(const Vec4 * const values,uint k, Vec4 &fv_l, Vec4 &fv_r, Vec4 &fd_l, Vec4 &fd_r){   
   //first use H8 estimates
   compute_h8_left_face_value(values, k, fv_l);
   compute_h8_left_face_value(values, k + 1, fv_r);   
   compute_h7_left_face_derivative(values, k, fd_l);
   compute_h7_left_face_derivative(values, k + 1, fd_r);
   
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

  //Fix left face if needed
  //Is boundary value bounded, that is, between its neighboring values.
  //Is slope ok (correct sign)
  Vec4b filter = (values[k -1] - fv_l) * (fv_l - values[k]) < 0 || slope_sign * fd_l < 0.0;
  if(horizontal_or (filter)) {  
     Vec4 h6_v_l,h4_v_l, h5_d_l,h3_d_l;
     compute_h6_left_face_value(values, k, h6_v_l);
     compute_h5_left_face_derivative(values, k, h5_d_l);
     compute_h4_left_face_value(values, k, h4_v_l);
     compute_h3_left_face_derivative(values, k, h3_d_l);
        
     //Go to H6 estimates if not ok
     fv_l = select(filter, h6_v_l, fv_l);
     fd_l = select(filter, h5_d_l, fd_l);

     //Go to H4 estimates if not ok
     filter = (values[k -1] - fv_l) * (fv_l - values[k]) < 0 || slope_sign * fd_l < 0.0;
     fv_l = select(filter, h4_v_l, fv_l);
     fd_l = select(filter, h3_d_l, fd_l);
        
     //Go to linear (PLM) estimates if not ok (this is always ok!)
     filter = (values[k -1] - fv_l) * (fv_l - values[k]) < 0 || slope_sign * fd_l < 0.0;
     fv_l=select(filter, values[k ] - slope_sign * 0.5 * slope_abs, fv_l);
     fd_l=select(filter, slope_sign * slope_abs, fd_l);
  }
  
  //Fix right face
  filter = (values[k + 1] - fv_r) * (fv_r - values[k]) < 0 || slope_sign * fd_r < 0.0;
  if(horizontal_or (filter)) {  
     Vec4 h6_v_r,h4_v_r, h5_d_r,h3_d_r;
     compute_h6_left_face_value(values, k + 1, h6_v_r);
     compute_h5_left_face_derivative(values, k + 1, h5_d_r);
     compute_h4_left_face_value(values, k + 1, h4_v_r);
     compute_h3_left_face_derivative(values, k + 1, h3_d_r);
     
     //Go to H6 estimates if not ok
     fv_r = select(filter, h6_v_r, fv_r);
     fd_r = select(filter, h5_d_r, fd_r);
     
     //Go to H4 estimates if not ok
     filter = (values[k + 1] - fv_r) * (fv_r - values[k]) < 0 || slope_sign * fd_r < 0.0;
     fv_r = select(filter, h4_v_r, fv_r);
     fd_r = select(filter, h3_d_r, fd_r);
     
     //Go to linear (PLM) estimates if not ok (this is always ok!)
     filter = (values[k + 1] - fv_r) * (fv_r - values[k]) < 0 || slope_sign * fd_r < 0.0;
     fv_r=select(filter, values[k] + slope_sign * 0.5 * slope_abs, fv_r);
     fd_r=select(filter, slope_sign * slope_abs, fd_r);
  }
}




/*Filters in section 2.6.1 of white et al. to be used for PQM
  1) Checks for extrema and flattens them
  2) Makes face values bounded
  3) Makes sure face slopes are consistent with PLM slope
*/
inline void compute_filtered_h6_face_values_derivatives(const Vec4 * const values,uint k, Vec4 &fv_l, Vec4 &fv_r, Vec4 &fd_l, Vec4 &fd_r){   
   //first use H& estimates
   compute_h6_left_face_value(values, k, fv_l);
   compute_h6_left_face_value(values, k + 1, fv_r);   
   compute_h5_left_face_derivative(values, k, fd_l);
   compute_h5_left_face_derivative(values, k + 1, fd_r);
   
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

  //Fix left face if needed
  //Is boundary value bounded, that is, between its neighboring values.
  //Is slope ok (correct sign)
  Vec4b filter = (values[k -1] - fv_l) * (fv_l - values[k]) < 0 || slope_sign * fd_l < 0.0;
  if(horizontal_or (filter)) {  
     Vec4 h4_v_l, h3_d_l;
     compute_h4_left_face_value(values, k, h4_v_l);
     compute_h3_left_face_derivative(values, k, h3_d_l);

     //Go to H4 if H6 is not ok
     fv_l = select(filter, h4_v_l, fv_l);
     fd_l = select(filter, h3_d_l, fd_l);
        
     //Go to linear (PLM) estimates if not ok (this is always ok!)
     filter = (values[k -1] - fv_l) * (fv_l - values[k]) < 0 || slope_sign * fd_l < 0.0;
     fv_l=select(filter, values[k ] - slope_sign * 0.5 * slope_abs, fv_l);
     fd_l=select(filter, slope_sign * slope_abs, fd_l);
  }
  
  //Fix right face
  filter = (values[k + 1] - fv_r) * (fv_r - values[k]) < 0 || slope_sign * fd_r < 0.0;
  if(horizontal_or (filter)) {  
     Vec4 h4_v_r,h3_d_r;
     compute_h4_left_face_value(values, k + 1, h4_v_r);
     compute_h3_left_face_derivative(values, k + 1, h3_d_r);
     
     //Go to H4 estimates if not ok
     fv_r = select(filter, h4_v_r, fv_r);
     fd_r = select(filter, h3_d_r, fd_r);
     
     //Go to linear (PLM) estimates if not ok (this is always ok!)
     filter = (values[k + 1] - fv_r) * (fv_r - values[k]) < 0 || slope_sign * fd_r < 0.0;
     fv_r=select(filter, values[k] + slope_sign * 0.5 * slope_abs, fv_r);
     fd_r=select(filter, slope_sign * slope_abs, fd_r);
  }
}




/*Filters in section 2.6.1 of white et al. to be used for PPM
  1) Checks for extrema and flattens them
  2) Makes face values bounded
  3) Makes sure face slopes are consistent with PLM slope
*/
inline void compute_filtered_h6_face_values(const Vec4 * const values,uint k, Vec4 &fv_l, Vec4 &fv_r){   
   //first use H8 estimates
   compute_h6_left_face_value(values, k, fv_l);
   compute_h6_left_face_value(values, k + 1, fv_r);   
   
   Vec4 slope_abs,slope_sign;
   slope_limiter(values[k -1], values[k], values[k + 1], slope_abs, slope_sign);
    
  //check for extrema, flatten if it is
  Vec4b is_extrema = (slope_abs == Vec4(0.0));
  if(horizontal_or(is_extrema)) {
    fv_r = select(is_extrema, values[k], fv_r);
    fv_l = select(is_extrema, values[k], fv_l);
  }

  //Fix left face if needed
  //Is boundary value bounded, that is, between its neighboring values.
  //Is slope ok (correct sign)
  Vec4b filter = (values[k -1] - fv_l) * (fv_l - values[k]) < 0 ;
  if(horizontal_or (filter)) {  
     Vec4 h4_v_l;
     compute_h4_left_face_value(values, k, h4_v_l);        
     //Go to H4 estimates if not ok
     fv_l = select(filter, h4_v_l, fv_l);
     //Go to linear (PLM) estimates if not ok (this is always ok!)
     filter = (values[k -1] - fv_l) * (fv_l - values[k]) < 0;
     fv_l=select(filter, values[k ] - slope_sign * 0.5 * slope_abs, fv_l);
  }
  
  //Fix right face
  filter = (values[k + 1] - fv_r) * (fv_r - values[k]) < 0;
  if(horizontal_or (filter)) {  
     Vec4 h4_v_r;
     compute_h4_left_face_value(values, k + 1, h4_v_r);
     //Go to H4 estimates if not ok
     fv_r = select(filter, h4_v_r, fv_r);
     //Go to linear (PLM) estimates if not ok (this is always ok!)
     filter = (values[k + 1] - fv_r) * (fv_r - values[k]) < 0;
     fv_r=select(filter, values[k] + slope_sign * 0.5 * slope_abs, fv_r);
  }
}




/*Filters in section 2.6.1 of white et al. to be used for PPM
  1) Checks for extrema and flattens them
  2) Makes face values bounded
  3) Makes sure face slopes are consistent with PLM slope
*/
inline void compute_filtered_h4_face_values(const Vec4 * const values,uint k, Vec4 &fv_l, Vec4 &fv_r){   
   //Use H4 estimates
   compute_h4_left_face_value(values, k, fv_l);
   compute_h4_left_face_value(values, k + 1, fv_r);   
   
   Vec4 slope_abs,slope_sign;
   slope_limiter(values[k -1], values[k], values[k + 1], slope_abs, slope_sign);
    
  //check for extrema, flatten if it is
  Vec4b is_extrema = (slope_abs == Vec4(0.0));
  if(horizontal_or(is_extrema)) {
    fv_r = select(is_extrema, values[k], fv_r);
    fv_l = select(is_extrema, values[k], fv_l);
  }

  //Fix left face if needed
  //Is boundary value bounded, that is, between its neighboring values.
  Vec4b filter = (values[k -1] - fv_l) * (fv_l - values[k]) < 0;
  if(horizontal_or (filter)) {  
     //Go to linear (PLM) estimates if not ok (this is always ok!)
     fv_l=select(filter, values[k ] - slope_sign * 0.5 * slope_abs, fv_l);
  }
  
  //Fix right face
  filter = (values[k + 1] - fv_r) * (fv_r - values[k]) < 0;
  if(horizontal_or (filter)) {  
     //Go to linear (PLM) estimates if not ok (this is always ok!)
     fv_r=select(filter, values[k] + slope_sign * 0.5 * slope_abs, fv_r);
  }
}




#endif
