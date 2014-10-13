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
inline void compute_ppm_coeff_explicit(const Vec4 * const values, face_estimate_order estimate, uint k, Vec4 a[3]){
  Vec4 fv_l; /*left face value*/
  Vec4 fv_r; /*right face value*/
  switch(estimate) {
      case h4:
         compute_filtered_h4_face_values(values, k, fv_l, fv_r); 
         break;
      case h6:
         compute_filtered_h6_face_values(values, k, fv_l, fv_r); 
         break;
  }
  

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
