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

const Vec4 one(1.0);
const Vec4 minus_one(-1.0);
const Vec4 two(2.0);
const Vec4 half(0.5);
const Vec4 zero(0.0);

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


void slope_limiter(const Vec4& l,const Vec4& m, const Vec4& r, Vec4& slope_abs, Vec4& slope_sign) {
  const Vec4 two(2.0);
  const Vec4 half(0.5);
  const Vec4 zero(0.0);
  Vec4 sign;
  Vec4 a=r-m;
  Vec4 b=m-l; 
  Vec4 minval=min(two*abs(a),two*abs(b));
  minval=min(minval,half*abs(a+b));
  
  //check for extrema, set absolute value
  slope_abs = select(a*b < 0, zero, minval);
  slope_sign = select(a + b < 0, minus_one, one);

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

#ifdef PPM_COELLA84
inline void compute_ppm_coeff(Vec4 mmv,Vec4 mv, Vec4 cv, Vec4 pv,Vec4 ppv,
			      Vec4 * __restrict__ a){
   //Compute estimations of the face values. Limited estimates, like in Coella 1984. 
   const Vec4 d_cv=slope_limiter(mv,cv,pv);
   const Vec4 d_pv=slope_limiter(cv,pv,ppv);
   const Vec4 d_mv=slope_limiter(mmv,mv,cv);
   
   Vec4 p_face=0.5*(pv+cv) + one_sixth * (d_cv-d_pv);
   Vec4 m_face=0.5*(cv+mv) + one_sixth * (d_mv-d_cv);        
   
   
   //Coella1984 eq. 1.10, detect extrema and make algorithm constant if it is
   Vec4 extrema_check = ((p_face - cv) * (cv - m_face));
   m_face = select(extrema_check < 0, cv, m_face);
   p_face = select(extrema_check < 0,cv, p_face);
   
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


#ifdef PPM_COELLA84_WITH_WHITE08_H5FACEVALS
inline void compute_ppm_coeff(Vec4 mmv,Vec4 mv, Vec4 cv, Vec4 pv,Vec4 ppv,
			      Vec4 * __restrict__ a){
   Vec4 p_face;
   Vec4 m_face;
   //white 08 H5 face estimates
   m_face = 1.0/60.0 * ( -3 * mmv + 27 * mv + 47 * cv - 13 * pv + 2 * ppv);
   p_face = 1.0/60.0 * ( 2 * mmv -13 * mv + 47 * cv + 27 * pv - 3 * ppv);
   
   Vec4 slope_abs,slope_sign;
   slope_limiter(mv,cv,pv,slope_abs,slope_sign);

   //detect and fix boundedness, as in WHITE 2008
   m_face = select((mv - m_face) * (m_face - cv) < 0,
                   cv - slope_sign * min( 0.5 * slope_abs, abs(m_face - cv)),
                   m_face);
   p_face = select((pv - p_face) * (p_face - cv) < 0,
                   cv + slope_sign * min( 0.5 * slope_abs, abs(p_face - cv)),
                   p_face);
   
   //Coella1984 eq. 1.10, detect extrema and make algorithm constant if it is
   Vec4 extrema_check = ((p_face - cv) * (cv - m_face));
   m_face = select(extrema_check < 0, cv, m_face);
   p_face = select(extrema_check < 0,cv, p_face);
   
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


#endif
