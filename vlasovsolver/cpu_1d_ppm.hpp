/*
This file is part of Vlasiator.
Copyright 2013, 2014 Finnish Meteorological Institute

*/

#ifndef CPU_1D_PPM_H
#define CPU_1D_PPM_H

#include <iostream>
#include "vec4.h"
#include "algorithm"
#include "cmath"
#include "cpu_slope_limiters.hpp"

using namespace std;

/*
  Compute parabolic reconstruction with an explicit scheme
*/
inline void compute_ppm_coeff_explicit(const Vec4 * const values, uint k, Vec4 a[3]){
  Vec4 fv_l; /*left face value*/
  Vec4 fv_r; /*right face value*/
  compute_h5_face_values(values, k ,fv_l, fv_r); 
   
  /*Filter boundedness according to Eq. 19 in White et al. 2008  Eq. 19 & 20*/
  bool fix_bounds = horizontal_or((values[k - 1] - fv_l) * (fv_l - values[k]) < 0 ||
				  (values[k + 1] - fv_r) * (fv_r - values[k]) < 0);
  if(fix_bounds) {
    Vec4 slope_abs,slope_sign;
    slope_limiter(values[k -1], values[k], values[k + 1], slope_abs, slope_sign);
    //detect and  fix boundedness, as in WHITE 2008
    fv_l = select((values[k -1] - fv_l) * (fv_l - values[k]) < 0,
		  values[k] - slope_sign * min( 0.5 * slope_abs, abs(fv_l - values[k])),
		  fv_l);
    fv_r = select((values[k + 1] - fv_r) * (fv_r - values[k]) < 0,
		  values[k] + slope_sign * min( 0.5 * slope_abs, abs(fv_r - values[k])),
		  fv_r);
  }


  //Coella1984 eq. 1.10, detect extrema
  Vec4 extrema_check = ((fv_r - values[k]) * (values[k] - fv_l));
  fv_l = select(extrema_check < 0, values[k], fv_l);
  fv_r = select(extrema_check < 0, values[k], fv_r);


  //Coella et al, check for monotonicity   
  Vec4 m_face = fv_l;
  Vec4 p_face = fv_r;
  m_face = select((p_face - m_face) * (values[k] - 0.5 * (m_face + p_face)) >
		  (p_face - m_face)*(p_face - m_face) * one_sixth,
		  3 * values[k] - 2 * p_face,
		  m_face);
  p_face = select(-(p_face - m_face) * (p_face - m_face) * one_sixth >
		  (p_face - m_face) * (values[k] - 0.5 * (m_face + p_face)),
		  3 * values[k] - 2 * m_face,
		  p_face);
     
  //Fit a second order polynomial for reconstruction see, e.g., White
  //2008 (PQM article) (note additional integration factors built in,
  //contrary to White (2008) eq. 4
  a[0] = m_face;
  a[1] = 3.0 * values[k] - 2.0 * m_face - p_face;
  a[2] = (m_face + p_face - 2.0 * values[k]);
}
     
#endif
