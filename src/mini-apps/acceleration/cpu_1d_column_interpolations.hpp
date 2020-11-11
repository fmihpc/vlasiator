/*
This file is part of Vlasiator.
Copyright 2013, 2014 Finnish Meteorological Institute

*/

#ifndef CPU_1D_COLUMN_INTERP_H
#define CPU_1D_COLUMN_INTERP_H

#include "algorithm"
#include "cmath"
#include "cpu_slope_limiters.hpp"

using namespace std;


/*Compute all face values. For cell k (globla index), its left face
 * value is in fv_l[k] and right value in fv_r[k]. Based on explicit
 * h6 estimate*/
inline void compute_h6_face_values(Real *values, uint n_cblocks,  Real *fv_l, Real *fv_r){   

   /*we loop up to one extra cell. There is extra space in fv for the extra left value*/
  for (int k = 0; k < n_cblocks * WID + 1; k++){
      /*compute left values*/
      fv_l[k] = 1.0/60.0 * (values[k - 3 + WID]  - 8.0 * values[k - 2 + WID]  + 37.0 * values[k - 1 + WID] +
                37.0 * values[k  + WID] - 8.0 * values[k + 1 + WID] + values[k + 2 + WID]);
      /*set right value*/
      if(k>0)
         fv_r[k-1] = fv_l[k];
  }
}


inline void filter_extrema(Real *values, uint n_cblocks, Real *fv_l, Real *fv_r){   
   for (int k = 0; k < n_cblocks * WID; k++){
      //Coella1984 eq. 1.10, detect extrema and make algorithm constant if it is
      Real extrema_check = ((fv_r[k] - values[k + WID]) * (values[k + WID] - fv_l[k]));
      fv_l[k] = extrema_check < 0 ? values[k + WID]: fv_l[k];
      fv_r[k] = extrema_check < 0 ? values[k + WID]: fv_r[k];
   }
}

/*Filter according to Eq. 19 in White et al.*/
inline void filter_boundedness(Real *values, uint n_cblocks, Real *fv_l, Real *fv_r){   
   /*First Eq. 19 & 20*/
   for (int k = 0; k < n_cblocks * WID; k++){
      bool do_fix_bounds =
         (values[k - 1 + WID] - fv_l[k]) * (fv_l[k] - values[k + WID]) < 0 ||
         (values[k + 1 + WID] - fv_r[k]) * (fv_r[k] - values[k + WID]) < 0;
      if(do_fix_bounds) {
         Real slope_abs,slope_sign;
         slope_limiter(values[k -1 + WID], values[k + WID], values[k + 1 + WID], slope_abs, slope_sign);
         //detect and  fix boundedness, as in WHITE 2008
         fv_l[k] = (values[k -1 + WID] - fv_l[k]) * (fv_l[k] - values[k + WID]) < 0 ?
            values[k + WID] - slope_sign * min( 0.5 * slope_abs, abs(fv_l[k] - values[k + WID])) :
            fv_l[k];
         fv_r[k] = (values[k + 1 + WID] - fv_r[k]) * (fv_r[k] - values[k + WID]) < 0 ?
            values[k + WID] + slope_sign * min( 0.5 * slope_abs, abs(fv_r[k] - values[k + WID])) :
            fv_r[k];
      }
   }
}


/*!
 Compute PLM coefficients
 f(v) = a[0] + a[1]/2.0*t 
 t=(v-v_{i-0.5})/dv where v_{i-0.5} is the left face of a cell
 The factor 2.0 is in the polynom to ease integration, then integral is a[0]*t + a[1]*t**2
*/

inline void compute_plm_coeff_explicit_column(Real *values, uint n_cblocks, Real a[][RECONSTRUCTION_ORDER + 1]){
   for (uint k = 0; k < n_cblocks * WID; k++){   
      const Real d_cv=slope_limiter(values[k - 1 + WID], values[k + WID], values[k + 1 + WID]);
      a[k][0] = values[k + WID] - d_cv * 0.5;
      a[k][1] = d_cv * 0.5;
   }
}

/*
  Compute parabolic reconstruction with an explicit scheme
  
  Note that value array starts with an empty block, thus values[k + WID]
  corresponds to the current (centered) cell.
*/

inline void compute_ppm_coeff_explicit_column(Real *values, uint n_cblocks, Real a[][RECONSTRUCTION_ORDER + 1]){
   Real p_face;
   Real m_face;
   Real fv_l[MAX_BLOCKS_PER_DIM * WID + 1]; /*left face value, extra space for ease of implementation*/
   Real fv_r[MAX_BLOCKS_PER_DIM * WID + 1]; /*right face value*/

   compute_h6_face_values(values,n_cblocks,fv_l, fv_r); 
   filter_boundedness(values,n_cblocks,fv_l, fv_r); 
   filter_extrema(values,n_cblocks,fv_l, fv_r);

   for (uint k = 0; k < n_cblocks * WID; k++){
      m_face = fv_l[k];
      p_face = fv_r[k];
      
      //Coella et al, check for monotonicity   
      m_face = (p_face - m_face) * (values[k + WID] - 0.5 * (m_face + p_face)) > (p_face - m_face)*(p_face - m_face) / 6.0 ?
         3 * values[k + WID] - 2 * p_face : m_face;
      p_face = -(p_face - m_face) * (p_face - m_face) / 6.0 > (p_face - m_face) * (values[k + WID] - 0.5 * (m_face + p_face)) ?
         3 * values[k + WID] - 2 * m_face : p_face;

      //Fit a second order polynomial for reconstruction see, e.g., White
      //2008 (PQM article) (note additional integration factors built in,
      //contrary to White (2008) eq. 4
      a[k][0] = m_face;
      a[k][1] = 3.0 * values[k + WID] - 2.0 * m_face - p_face;
      a[k][2] = (m_face + p_face - 2.0 * values[k + WID]);
   }
}



#endif
