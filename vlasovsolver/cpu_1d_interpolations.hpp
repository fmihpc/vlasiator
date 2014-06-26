/*
This file is part of Vlasiator.
Copyright 2013, 2014 Finnish Meteorological Institute

*/

#ifndef CPU_1D_INTERP_H
#define CPU_1D_INTERP_H

#include "vec4.h"
#include "algorithm"
#include "cmath"
#include "cpu_slope_limiters.hpp"

using namespace std;



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
   p_face = 1.0/60.0 * (  2 * mmv - 13 * mv + 47 * cv + 27 * pv - 3 * ppv);


   bool fix_bounds = horizontal_or((mv - m_face) * (m_face - cv) < 0 || (pv - p_face) * (p_face - cv) < 0);
   if(fix_bounds) {
      Vec4 slope_abs,slope_sign;
      slope_limiter(mv,cv,pv,slope_abs,slope_sign);

      //detect and fix boundedness, as in WHITE 2008
      m_face = select((mv - m_face) * (m_face - cv) < 0,
                      cv - slope_sign * min( 0.5 * slope_abs, abs(m_face - cv)),
                      m_face);
      p_face = select((pv - p_face) * (p_face - cv) < 0,
                      cv + slope_sign * min( 0.5 * slope_abs, abs(p_face - cv)),
                      p_face);
   }
   
   //Coella1984 eq. 1.10, detect extrema and make algorithm constant if it is
   Vec4 extrema_check = ((p_face - cv) * (cv - m_face));
   m_face = select(extrema_check < 0, cv, m_face);
   p_face = select(extrema_check < 0, cv, p_face);
   
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


#ifdef PPM_COELLA08
// A Limiter for PPM that Preserves Accuracy at Smooth Extrema, Coella
// et al.  Journal of Computational Physics, Volume 227, Issue 15,
// p. 7069-7076.


inline void compute_ppm_coeff(Vec4 mmv,Vec4 mv, Vec4 cv, Vec4 pv,Vec4 ppv,
			      Vec4 * __restrict__ a){
   Vec4 p_face;
   Vec4 m_face;
   //white 08 H5 face estimates (H6 estimates in paper would require wider stencil(!)
   //m_face = 1.0/60.0 * ( -3 * mmv + 27 * mv + 47 * cv - 13 * pv + 2 * ppv);
   //p_face = 1.0/60.0 * ( 2 * mmv -13 * mv + 47 * cv + 27 * pv - 3 * ppv);

   p_face = 7.0/12.0 * ( pv + cv ) - 1.0/12.0 * ( ppv + mv );
   m_face = 7.0/12.0 * ( cv + mv ) - 1.0/12.0 * ( pv + mmv );
   
   
   //TODO PERF, we compute a lot of things now just in case we need it if face values not bounded
   //the boundedness things. Let's be a bit smarter (compute it for
   //all if any cell has boundedness problem? )
   //computing D2_lim in Eq. (18)
   Vec4 D2vm_face = 3.0 *(mv - 2.0 * m_face + cv); //second derivative at m_face (why 3?)
   Vec4 D2vp_face = 3.0 *(cv - 2.0 * p_face + pv); //second derivative at p_face (why 3, why not 4??)
   Vec4 D2vm  = mmv - 2.0 * mv + cv;         //second derivative at m cell
   Vec4 D2vc  = mv  - 2.0 * cv + pv;            //second derivative at c cell
   Vec4 D2vp  = cv  - 2.0 * pv + ppv;         //second derivative at p cell

   const Real C = 1.25;
   Vec4 s;
   s = select(sign_bit(D2vm_face),-1.0, 1.0); //sign
   Vec4 D2vm_lim = max( min( min(C * s * D2vm, C * s *  D2vc), s * D2vm_face), zero);

   s = select(sign_bit(D2vp_face),-1.0, 1.0);
   Vec4 D2vp_lim = max( min( min(C * s * D2vc, C * s *  D2vp), s * D2vp_face), zero);

   /*and now modify face values if they are out of bounds*/
   m_face = select((mv - m_face) * (m_face - cv) < 0,
                   0.5*(cv + mv - one_third * D2vm_lim),
                   m_face);
   p_face = select((pv - p_face) * (p_face - cv) < 0,
                   0.5*(pv + cv - one_third * D2vp_lim),
                   p_face);
   
/*
  Vec4 D2vp_lim =  select(D2vp * D2vc > 0 && D2vp_face * D2vc > 0,
                           D2vp_sign * min( min(C * abs(D2vc), C * abs(D2vp)), abs(D2vp_face)),
                           zero);
*/ 


   //TODO PERF, we compute a lot of things now just in case we need if extrema
   //Coella2008 eq. 20
   Vec4 extrema_check1 = ((p_face - cv) * (cv - m_face));
   Vec4 extrema_check2 = ((pv - cv) * (cv - mv));

   //first one in Eq 21 (the others are already computed above)
   Vec4 D2a6 = 6.0 * (m_face - 2.0 * cv  + p_face); //why 6???
   Vec4 D2a6_sign = select(sign_bit(D2a6),-1.0, 1.0);
   //if all D2's have the same sign, then compute limited S2a6, otherwise it is set to zero
   Vec4 D2a6_lim = select( D2a6 * D2vm >= 0 && D2vm * D2vc >= 0 && D2vc * D2vp >= 0,
                           D2a6_sign * min(min(min( C * abs(D2vm), C * abs(D2vc)), C* abs(D2vp)), abs(D2a6)),
                           zero);
   
   Vec4 ratio = select(D2a6 != zero, D2a6_lim/D2a6, zero);
   m_face = select(extrema_check1 <= 0 || extrema_check2 <= 0,
                   cv + (m_face - cv) * ratio,
                   m_face);
   p_face = select(extrema_check1 <= 0 || extrema_check2 <= 0,
                   cv + (p_face - cv) * ratio,
                   p_face);

   


   //Coella1984 eq. 1.10, detect extrema and make algorithm constant if it is
   /* Vec4 extrema_check = ((p_face - cv) * (cv - m_face));
   m_face = select(extrema_check < 0, cv, m_face);
   p_face = select(extrema_check < 0,cv, p_face);
   */
   
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
