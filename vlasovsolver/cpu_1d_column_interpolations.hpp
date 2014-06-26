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


void tridiag(uint n,const Vec4 &a,const Vec4 &b,const Vec4 &c,Vec4 *u, Vec4 *r){
  Vec4 gam[MAX_BLOCKS_PER_DIM*WID];
  
  Vec4 bet=b;
  u[0]=r[0]/bet;
  for(int i = 1;i< n ;i++){
    gam[i] = c/bet;
    bet = b - a * gam[i];
    u[i] = (r[i] - a * u[i-1]) / bet;
  }
  
  for(int i = (n-2); i>=0; i--){
    u[i] -= gam[i+1]*u[i+1];
  }
}


/*Compute all face values. For cell k (globla index), its left face
 * value is in fv_l[k] and right value in fv_r[k]. Based on explicit
 * h6 estimate*/
inline void compute_h6_face_values(Vec4 *values, uint n_cblocks,  Vec4 *fv_l, Vec4 *fv_r){   

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

inline void compute_ih4_face_values(Vec4 *values, uint n_cblocks,  Vec4 *fv_l, Vec4 *fv_r){
   Vec4 r[MAX_BLOCKS_PER_DIM*WID + 1];
  const Vec4 a(1.0/4.0);
  const Vec4 b(1.0);
  const Vec4 c(1.0/4.0);
  
  for (uint k = 0; k < n_cblocks * WID + 1; k++){
    //Note that these are now defined for left face values, so not exactly like in white
    r[k] = ( 3.0 * values[k - 1 + WID] + 3.0 * values[k + WID])/4.0;
  }
  
  tridiag(n_cblocks*WID + 1,a,b,c,fv_l,r);

  /*copy right face values */
  for (uint k = 1; k < n_cblocks * WID + 1; k++){
    fv_r[k-1] = fv_l[k];
  }

}



/*Compute all face values. For cell k (globla index), its left face
 * value is in fv_l[k] and right value in fv_r[k]. Based on implicit
 * ih6 estimate from White - 2008. Note that the value for a and c is wrong in the paper. The original source is:

[ Lacor, Chris, Sergey Smirnov, and Martine Baelmans. “A Finite Volume Formulation of Compact Central Schemes on Arbitrary Structured Grids.” Journal of Computational Physics 198, no. 2 (August 2004): 535–66. doi:10.1016/j.jcp.2004.01.025.]

and solving for the variables there we get for the variable names in that paper (sympy): 
solve([Eq(a + b, 1 + 2*alpha), Eq(a+(2**3 - 1)*b,6*alpha), Eq(a + (2**5 -1)*b,10*alpha)],a,alpha,b)
Out[17]: {alpha: 1/3, b: 1/18, a: 29/18} 
*/

inline void compute_ih6_face_values(Vec4 *values, uint n_cblocks,  Vec4 *fv_l, Vec4 *fv_r){   
  Vec4 r[MAX_BLOCKS_PER_DIM*WID + 1];
  const Vec4 a(1.0/3.0);
  const Vec4 b(1.0);
  const Vec4 c(1.0/3.0);
  
  for (uint k = 0; k < n_cblocks * WID + 1; k++){
    //Note that these are now defined for left face values, so not exactly like in white
    r[k] = (values[k - 2 + WID] + 29.0 * values[k - 1 + WID] + 29.0 * values[k + WID] + values[k + 1 + WID])/36.0;
  }
  
  tridiag(n_cblocks*WID + 1,a,b,c,fv_l,r);

  /*copy right face values */
  for (uint k = 1; k < n_cblocks * WID + 1; k++){
    fv_r[k-1] = fv_l[k];
  }

}

/*
Out[43]: 
[a + b + c == 2*alpha + 1,
 a + 7*b + 19*c == 6*alpha,
 a + 31*b + 211*c == 10*alpha,
 a + 127*b + 2059*c == 14*alpha]

In [44]: solve(equations,a,b,c,alpha)
Out[44]: {c: -1/240, alpha: 3/8, b: 23/240, a: 199/120}
*/


inline void compute_ih8_face_values(Vec4 *values, uint n_cblocks,  Vec4 *fv_l, Vec4 *fv_r){   
  Vec4 r[MAX_BLOCKS_PER_DIM*WID + 1];
  const Vec4 a(3.0/8.0);
  const Vec4 b(1.0);
  const Vec4 c(3.0/8.0);
  
  for (uint k = 0; k < n_cblocks * WID + 1; k++){
    //Note that these are now defined for left face values, so not exactly like in white
    r[k] = ( -1.0*values[k-3+WID] + 23.0*values[k - 2 + WID] + 398.0*values[k - 1 + WID] + 398.0 * values[k + WID] + 23.0* values[k + 1 + WID] - 1.0*values[k+2+WID])/480.0;
  }
  
  tridiag(n_cblocks*WID + 1,a,b,c,fv_l,r);

  /*copy right face values */
  for (uint k = 1; k < n_cblocks * WID + 1; k++){
    fv_r[k-1] = fv_l[k];
  }

}


inline void filter_extrema(Vec4 *values, uint n_cblocks, Vec4 *fv_l, Vec4 *fv_r){   
   for (int k = 0; k < n_cblocks * WID; k++){
      //Coella1984 eq. 1.10, detect extrema and make algorithm constant if it is
      Vec4 extrema_check = ((fv_r[k] - values[k + WID]) * (values[k + WID] - fv_l[k]));
      fv_l[k] = select(extrema_check < 0, values[k + WID], fv_l[k]);
      fv_r[k] = select(extrema_check < 0, values[k + WID], fv_r[k]);
   }
}

/*Filter according to Eq. 19 in White et al.*/
inline void filter_boundedness(Vec4 *values, uint n_cblocks, Vec4 *fv_l, Vec4 *fv_r){   
   /*First Eq. 19 & 20*/
   for (int k = 0; k < n_cblocks * WID; k++){
      bool fix_bounds = horizontal_or((values[k - 1 + WID] - fv_l[k]) * (fv_l[k] - values[k + WID]) < 0 ||
                                      (values[k + 1 + WID] - fv_r[k]) * (fv_r[k] - values[k + WID]) < 0);
      if(fix_bounds) {
         Vec4 slope_abs,slope_sign;
         slope_limiter(values[k -1 + WID], values[k + WID], values[k + 1 + WID], slope_abs, slope_sign);
         //detect and  fix boundedness, as in WHITE 2008
         fv_l[k] = select((values[k -1 + WID] - fv_l[k]) * (fv_l[k] - values[k + WID]) < 0,
                         values[k + WID] - slope_sign * min( 0.5 * slope_abs, abs(fv_l[k] - values[k + WID])),
                         fv_l[k]);
         fv_r[k] = select((values[k + 1 + WID] - fv_r[k]) * (fv_r[k] - values[k + WID]) < 0,
                                values[k + WID] + slope_sign * min( 0.5 * slope_abs, abs(fv_r[k] - values[k + WID])),
                                fv_r[k]);
      }
   }
}

/*"A final check is conducted to make sure
  that discontinuous edge values are monotonic at any edge. If edge values are discontinuous and nonmono-
  tonic, they are both replaced by their average" */
inline void filter_face_monotonicity(Vec4 *values, uint n_cblocks, Vec4 *fv_l, Vec4 *fv_r){   
   for (int k = 0; k > n_cblocks * WID - 1; k++){
      /*shift to mean is how much we need to add (subract) from right
       * cell face (left cell face of next cell) to make both values have the average value. It is zero if no
       * filtering needs to be done here*/
      Vec4 shift_to_mean = select( (fv_l[k+1] - fv_r[k]) * (values[k + 1 + WID] - values[k + WID]) < 0.0,
                                   (fv_l[k+1] + fv_r[k]) * 0.5 - fv_r[k],
                                   0.0);
      fv_r[k] += shift_to_mean;
      fv_l[k+1] -= shift_to_mean;
   }
}


/*!
 Compute PLM coefficients
 f(v) = a[0] + a[1]/2.0*t 
t=(v-v_{i-0.5})/dv where v_{i-0.5} is the left face of a cell
The factor 2.0 is in the polynom to ease integration, then integral is a[0]*t + a[1]*t**2
*/

inline void compute_plm_coeff_explicit_column(Vec4 *values, uint n_cblocks, Vec4 a[][RECONSTRUCTION_ORDER + 1]){
   for (uint k = 0; k < n_cblocks * WID; k++){   
      const Vec4 d_cv=slope_limiter(values[k - 1 + WID], values[k + WID], values[k + 1 + WID]);
      a[k][0] = values[k + WID] - d_cv * 0.5;
      a[k][1] = d_cv * 0.5;
   }
}

/*
  Compute parabolic reconstruction with an explicit scheme
  
  Note that value array starts with an empty block, thus values[k + WID]
  corresponds to the current (centered) cell.
*/

inline void compute_ppm_coeff_explicit_column(Vec4 *values, uint n_cblocks, Vec4 a[][RECONSTRUCTION_ORDER + 1]){
   Vec4 p_face;
   Vec4 m_face;
   Vec4 fv_l[MAX_BLOCKS_PER_DIM * WID + 1]; /*left face value, extra space for ease of implementation*/
   Vec4 fv_r[MAX_BLOCKS_PER_DIM * WID + 1]; /*right face value*/

   //white 08 H6 face estimates, better than H5. Shift old unfilitered value at upper edge to the lower edge (identical edge)     

   compute_ih4_face_values(values,n_cblocks,fv_l, fv_r); 
   filter_boundedness(values,n_cblocks,fv_l, fv_r); 
   filter_extrema(values,n_cblocks,fv_l, fv_r);
//   filter_face_monotonicity(values,n_cblocks,fv);
   
   for (uint k = 0; k < n_cblocks * WID; k++){
      m_face = fv_l[k];
      p_face = fv_r[k];
      
      //Coella et al, check for monotonicity   
      m_face = select((p_face - m_face) * (values[k + WID] - 0.5 * (m_face + p_face)) >
                      (p_face - m_face)*(p_face - m_face) * one_sixth,
                      3 * values[k + WID] - 2 * p_face,
                      m_face);
      p_face = select(-(p_face - m_face) * (p_face - m_face) * one_sixth >
                      (p_face - m_face) * (values[k + WID] - 0.5 * (m_face + p_face)),
                      3 * values[k + WID] - 2 * m_face,
                      p_face);

      //Fit a second order polynomial for reconstruction see, e.g., White
      //2008 (PQM article) (note additional integration factors built in,
      //contrary to White (2008) eq. 4
      a[k][0] = m_face;
      a[k][1] = 3.0 * values[k + WID] - 2.0 * m_face - p_face;
      a[k][2] = (m_face + p_face - 2.0 * values[k + WID]);
   }
}



// /*
//   PQM reconstruction as published in:
//   White, Laurent, and Alistair Adcroft. “A High-Order Finite Volume Remapping Scheme for Nonuniform Grids: The Piecewise Quartic Method (PQM).” Journal of Computational Physics 227, no. 15 (July 2008): 7394–7422. doi:10.1016/j.jcp.2008.04.026.
// */

// inline void compute_pqm_coeff_explicit_column(Real *values, uint n_cblocks, uint j, Vec4 a[][RECONSTRUCTION_ORDER + 1]){ }

#endif
