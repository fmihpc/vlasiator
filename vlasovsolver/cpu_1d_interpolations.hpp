
/*
This file is part of Vlasiator.
Copyright 2013, 2014 Finnish Meteorological Institute

*/

#ifndef CPU_1D_INTERP_H
#define CPU_1D_INTERP_H

#include "vec4.h"
#include "algorithm"
#include "cmath"

using namespace std;

const Vec4 one_sixth(1.0/6.0);
const Vec4 one_twelfth(1.0/12.0);
const Vec4 seven_twelfth(7.0/12.0);
const Vec4 one_third(1.0/3.0);

// indices in padded z block
#define i_pblock(i,j,k) ( ((k) + STENCIL_WIDTH ) * WID + (j) * WID * (WID + 2* STENCIL_WIDTH) + (i) )
#define i_pblockv(j,k) ( ((k) + STENCIL_WIDTH ) * WID + (j) * WID * (WID + 2* STENCIL_WIDTH) )


/*!
  MC slope limiter
*/

inline Vec4 slope_limiter(const Vec4& l,const Vec4& m, const Vec4& r) {
  const Vec4 two(2.0);
  const Vec4 half(0.5);
  const Vec4 zero(0.0);
  Vec4 sign;
  Vec4 a=r-m;
  Vec4 b=m-l; 
  Vec4 minval=min(two*abs(a),two*abs(b));
  minval=min(minval,half*abs(a+b));
  
  //check for extrema
  Vec4 output = select(a*b < 0,zero,minval);
  
  //set sign
  return select(a + b < 0,-output,output);
}


/*!
 Compute PLM coefficients
f(v) = a[0] + a[1]/2.0*t 
t=(v-v_{i-0.5})/dv where v_{i-0.5} is the left face of a cell
The factor 2.0 is in the polynom to ease integration, then integral is a[0]*t + a[1]*t**2
*/

inline void compute_plm_coeff(Vec4 mv,Vec4 cv, Vec4 pv,
			      Vec4 * __restrict__ a){
  const Vec4 d_cv=slope_limiter(mv,cv,pv);
  a[0] = cv - d_cv * 0.5;
  a[1] = d_cv * 0.5;
}


/*!
 Compute PPM coefficients
f(v) = a[0] + a[1]/2.0*t + a[2]/3.0*t**2 
t=(v-v_{i-0.5})/dv where v_{i-0.5} is the left face of a cell
The factor 2.0 & 3.0 is in the polynom to ease integration, then integral is a[0]*t + a[1]*t**2 + a[2]*t**3
*/


inline void compute_ppm_coeff(Vec4 mmv,Vec4 mv, Vec4 cv, Vec4 pv,Vec4 ppv,
			      Vec4 * __restrict__ a){
   Vec4 p_face;
   Vec4 m_face;

   //Compute estimations of the face values. Limited estimates, like in Coella 1984. 
   const Vec4 d_cv=slope_limiter(mv,cv,pv);
   const Vec4 d_pv=slope_limiter(cv,pv,ppv);
   const Vec4 d_mv=slope_limiter(mmv,mv,cv);

   p_face=0.5*(pv+cv) + one_sixth * (d_cv-d_pv);
   m_face=0.5*(cv+mv) + one_sixth * (d_mv-d_cv);        

   //Coella1984 eq. 1.10, detect extrema
   m_face = select((p_face - cv) * (cv - m_face) < 0,cv, m_face);
   p_face = select((p_face - cv) * (cv - m_face) < 0,cv, p_face);
   
   //Coella et al, check for monotonicity   
   m_face = select((p_face-m_face)*(cv-0.5*(m_face+p_face))>(p_face-m_face)*(p_face-m_face)*one_sixth,
		  3*cv-2*p_face,
		  m_face);
   p_face = select(-(p_face-m_face)*(p_face-m_face)*one_sixth > (p_face-m_face)*(cv-0.5*(m_face+p_face)),
		   3*cv-2*m_face,
		   p_face);

   //Fit a second order polynomial for reconstruction see, e.g., White
   //2008 (PQM article) (note additional integration factors built in,
   //contrary to White (2008) eq. 4
   a[0]=m_face;
   a[1]=3.0*cv-2.0*m_face-p_face;
   a[2]=(m_face+p_face-2.0*cv);
}

#endif
