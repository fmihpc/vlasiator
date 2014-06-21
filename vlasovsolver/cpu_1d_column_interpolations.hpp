/*
This file is part of Vlasiator.
Copyright 2013, 2014 Finnish Meteorological Institute

*/

#ifndef CPU_1D_COLUMN_INTERP_H
#define CPU_1D_COLUMN_INTERP_H

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

inline void compute_plm_coeff_explicit_column(Real *values, uint n_cblocks, uint j, Vec4 a[][RECONSTRUCTION_ORDER + 1]){
   Vec4 mv,cv,pv; /* values in z-direction*/
   cv.load(values + i_pcolumnv(n_cblocks, 0, j, -1));
   pv.load(values + i_pcolumnv(n_cblocks, 0, j, 0));
   for (uint block_i = 0; block_i < n_cblocks; block_i++){
      for (uint k = 0; k < WID; ++k){
         mv = cv;
         cv = pv;
         pv.load(values + i_pcolumnv(n_cblocks, block_i, j, k + 1));
         const Vec4 d_cv=slope_limiter(mv,cv,pv);
         a[block_i * WID + k][0] = cv - d_cv * 0.5;
         a[block_i * WID + k][0] = d_cv * 0.5;
      }
   }
}




inline void compute_ppm_coeff_explicit_column(Real *values, uint n_cblocks, uint j, Vec4 a[][RECONSTRUCTION_ORDER + 1]){
   Vec4 p_face;
   Vec4 m_face;
   Vec4 p_face_unfiltered;
   Vec4 mmv,mv,cv,pv,ppv,pppv; /* values in z-direction*/
   /*set up shifting operation. These will be shifted to mmv - pv*/
   mmv.load(values + i_pcolumnv(n_cblocks, 0, j, -3));
   mv.load(values + i_pcolumnv(n_cblocks, 0, j, -2));
   cv.load(values + i_pcolumnv(n_cblocks, 0, j, -1));
   pv.load(values + i_pcolumnv(n_cblocks, 0, j, 0));
   ppv.load(values + i_pcolumnv(n_cblocks, 0, j, 1));
   pppv.load(values + i_pcolumnv(n_cblocks, 0, j, 2));
   p_face_unfiltered = 1.0/60.0 * (mmv  - 8.0 * mv  + 37.0 * cv + 37.0 * pv - 8.0 * ppv + pppv);   
   for (uint block_i = 0; block_i < n_cblocks; block_i++){
      for (uint k = 0; k < WID; ++k){
         /*shift values*/
         mmv = mv;
         mv = cv;
         cv = pv;
         pv = ppv;
         ppv = pppv;
         pppv.load(values + i_pcolumnv(n_cblocks, block_i, j, k + 3));
         //white 08 H6 face estimates, better than H5. Shift old unfilitered value at upper edge to the lower edge (identical edge)
         m_face = p_face_unfiltered;
         p_face_unfiltered = 1.0/60.0 * (mmv  - 8.0 * mv  + 37.0 * cv + 37.0 * pv - 8.0 * ppv + pppv);
         p_face = p_face_unfiltered;
         //white 08 H5 face estimates
         //m_face = 1.0/60.0 * ( -3 * mmv + 27 * mv + 47 * cv - 13 * pv + 2 * ppv);
         //p_face = 1.0/60.0 * (  2 * mmv - 13 * mv + 47 * cv + 27 * pv - 3 * ppv);
         bool fix_bounds = horizontal_or((mv - m_face) * (m_face - cv) < 0 || (pv - p_face) * (p_face - cv) < 0);
         if(fix_bounds) {
            Vec4 slope_abs,slope_sign;
            slope_limiter(mv, cv, pv, slope_abs, slope_sign);
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
         m_face = select((p_face - m_face) * (cv - 0.5 * (m_face + p_face))>(p_face - m_face)*(p_face - m_face) * one_sixth,
                         3 * cv - 2 * p_face,
                         m_face);
         p_face = select(-(p_face - m_face) * (p_face - m_face) * one_sixth > (p_face - m_face) * (cv - 0.5 * (m_face + p_face)),
                         3 * cv - 2 * m_face,
                         p_face);

         //Fit a second order polynomial for reconstruction see, e.g., White
         //2008 (PQM article) (note additional integration factors built in,
         //contrary to White (2008) eq. 4
         a[block_i * WID + k][0] = m_face;
         a[block_i * WID + k][1] = 3.0 * cv - 2.0 * m_face - p_face;
         a[block_i * WID + k][2] = (m_face + p_face - 2.0 * cv);
      }
   }
}

#endif
