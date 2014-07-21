/*
This file is part of Vlasiator.
Copyright 2013, 2014 Finnish Meteorological Institute

*/

#ifndef CPU_1D_COLUMN_INTERP_H
#define CPU_1D_COLUMN_INTERP_H

#include <iostream>
#include "vec4.h"
#include "algorithm"
#include "cmath"
#include "cpu_slope_limiters.hpp"

using namespace std;

/*
 * Simple tridiagonal solver, based on NR.
 */

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
 * h4 estimate*/
inline void compute_h4_face_values(Vec4 *values, uint n_cblocks,  Vec4 *fv_l, Vec4 *fv_r){   

   /*we loop up to one extra cell. There is extra space in fv for the extra left value*/
  for (int k = 0; k < n_cblocks * WID + 1; k++){
      /*compute left values*/
    fv_l[k] = 1.0/12.0 * ( - 1.0 * values[k - 2 + WID]  
			   + 7.0 * values[k - 1 + WID] 
			   + 7.0 * values[k  + WID] 
			   - 1.0 * values[k + 1 + WID]);
    /*set right value*/
    if(k>0)
      fv_r[k-1] = fv_l[k];
  }
}



/*Compute all face values. For cell k (globla index), its left face
 * value is in fv_l[k] and right value in fv_r[k]. Based on explicit
 *      h5 estimate*/
inline void compute_h5_face_values(Vec4 *values, uint n_cblocks,  Vec4 *fv_l, Vec4 *fv_r){   
   
   /*we loop up to one extra cell. There is extra space in fv for the extra left value*/
  for (int k = 0; k < n_cblocks * WID + 1; k++){
      /*compute left values*/
     fv_l[k] = 1.0/60.0 * (- 3.0 * values[k - 2 + WID]  
			   + 27.0 * values[k - 1 + WID] 
			   + 47.0 * values[k  + WID] 
			   - 13.0 * values[k + 1 + WID] 
			   + 2.0 * values[k + 2 + WID]);
     fv_r[k] = 1.0/60.0 * ( 2.0 * values[k - 2 + WID] 
			    - 13.0 * values[k - 1 + WID] 
                            + 47.0 * values[k  + WID]
			    + 27.0 * values[k + 1 + WID] 
			    - 3 * values[k + 2 + WID]);
  }
}



/*Compute all face values. For cell k (globla index), its left face
 * value is in fv_l[k] and right value in fv_r[k]. Based on explicit
 * h6 estimate*/
inline void compute_h6_face_values(Vec4 *values, uint n_cblocks,  Vec4 *fv_l, Vec4 *fv_r){   

   /*we loop up to one extra cell. There is extra space in fv for the extra left value*/
  for (int k = 0; k < n_cblocks * WID + 1; k++){
    /*compute left values*/
    fv_l[k] = 1.0/60.0 * (values[k - 3 + WID]  
			  - 8.0 * values[k - 2 + WID]  
			  + 37.0 * values[k - 1 + WID] 
			  + 37.0 * values[k  + WID] 
			  - 8.0 * values[k + 1 + WID] 
			  + values[k + 2 + WID]);
    /*set right value*/
    if(k>0)
      fv_r[k-1] = fv_l[k];
  }
}
/*Compute all face derivatives. For cell k (globla index), its left face
 * derivative is in fd_l[k] and right derivative in fd_r[k]. Based on explicit
 * h5 estimate*/
inline void compute_h5_face_derivatives(Vec4 *values, uint n_cblocks,  Vec4 *fd_l, Vec4 *fd_r){   

   /*we loop up to one extra cell. There is extra space in fv for the extra left value*/
  for (int k = 0; k < n_cblocks * WID + 1; k++){
    /*compute left values*/
    fd_l[k] = 1.0/180.0 * (245 * (values[k + WID] - values[k - 1 + WID])  
			   - 25 * (values[k + 1 + WID] - values[k - 2 + WID]) 
			   + 2 * (values[k + 2 + WID] - values[k - 3 + WID]));
    /*set right value*/
    if(k>0)
      fd_r[k-1] = fd_l[k];
  }
}


/*Compute all face derivatives. For cell k (globla index), its left face
 * derivative is in fd_l[k] and right derivative in fd_r[k]. Based on explicit
 * h5 estimate*/
inline void compute_h4_face_derivatives(Vec4 *values, uint n_cblocks,  Vec4 *fd_l, Vec4 *fd_r){   

   /*we loop up to one extra cell. There is extra space in fd for the extra left value*/
  for (int k = 0; k < n_cblocks * WID + 1; k++){
    /*compute left values*/
    fd_l[k] = 1.0/12.0 * (15.0 * (values[k + WID] - values[k - 1 + WID])
			  - (values[k + 1 + WID] - values[k - 2 + WID]));
    /*set right value*/
    if(k>0)
      fd_r[k-1] = fd_l[k];
  }
}


/*Compute all face values. For cell k (globla index), its left face
 * value is in fv_l[k] and right value in fv_r[k]. Based on implicit
 * ih4 estimate from White - 2008.
 */


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


  /*Compute all face values. For cell k (globla index), its left face
 * value is in fv_l[k] and right value in fv_r[k]. Based on implicit
 * ih8 estimates derived based on:

 [ Lacor, Chris, Sergey Smirnov, and Martine Baelmans. “A Finite Volume Formulation of Compact Central Schemes on Arbitrary Structured Grids.” Journal of Computational Physics 198, no. 2 (August 2004): 535–66. doi:10.1016/j.jcp.2004.01.025.]

 sympy:
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
    r[k] = ( -1.0*values[k-3+WID] + 23.0*values[k - 2 + WID] +
             398.0*values[k - 1 + WID] + 398.0 * values[k + WID] +
             23.0* values[k + 1 + WID] - 1.0*values[k+2+WID] )/480.0;
  }
  
  tridiag(n_cblocks*WID + 1,a,b,c,fv_l,r);

  /*copy right face values */
  for (uint k = 1; k < n_cblocks * WID + 1; k++){
    fv_r[k-1] = fv_l[k];
  }

}

//Coella1984 eq. 1.10, detect extrema and make values
inline void filter_extrema(Vec4 *values, uint n_cblocks, Vec4 *fv_l, Vec4 *fv_r){   
   for (int k = 0; k < n_cblocks * WID; k++){
      Vec4 extrema_check = ((fv_r[k] - values[k + WID]) * (values[k + WID] - fv_l[k]));
      fv_l[k] = select(extrema_check < 0, values[k + WID], fv_l[k]);
      fv_r[k] = select(extrema_check < 0, values[k + WID], fv_r[k]);
   }
}

/*Filter boundedness according to Eq. 19 in White et al.*/
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


/*Filters in section 2.6.1 of white et al. to be used for PQM
  1) Checks for extrema and flattens them
  2) Makes face values bounded
  3) Makes sure face slopes are consistent with PLM slope
*/
inline void filter_extrema_boundedness(Vec4 *values, uint n_cblocks, Vec4 *fv_l, Vec4 *fv_r, Vec4 *fd_l, Vec4 *fd_r){   
  /*First Eq. 19 & 20*/
  for (int k = 0; k < n_cblocks * WID; k++){
    Vec4 slope_abs,slope_sign;
    slope_limiter(values[k -1 + WID], values[k + WID], values[k + 1 + WID], slope_abs, slope_sign);

    //check for extrema, flatten if it is
    Vec4b is_extrema = (slope_abs == Vec4(0.0));
    bool fix_extrema = horizontal_or(is_extrema);
    if(fix_extrema) {
      fv_r[k] = select(is_extrema, values[k + WID], fv_r[k]);
      fv_l[k] = select(is_extrema, values[k + WID], fv_l[k]);
      fd_l[k] = select(is_extrema, 0.0 , fd_l[k]);
      fd_r[k] = select(is_extrema, 0.0 , fd_r[k]);
    }
    
    //Check that boundary values are bounded, that is, between their neighboring values
    Vec4b left_unbounded = (values[k -1 + WID] - fv_l[k]) * (fv_l[k] - values[k + WID]) < 0;
    Vec4b right_unbounded = (values[k + 1 + WID] - fv_r[k]) * (fv_r[k] - values[k + WID]) < 0;
    bool fix_bounds = horizontal_or(left_unbounded || right_unbounded );
     if(fix_bounds) {
       //detect and  fix boundedness, as in WHITE 2008
       fv_l[k] = select(left_unbounded,
			values[k + WID] - slope_sign * min( 0.5 * slope_abs, abs(fv_l[k] - values[k + WID])),
			fv_l[k]);
       fv_r[k] = select(right_unbounded,
			values[k + WID] + slope_sign * min( 0.5 * slope_abs, abs(fv_r[k] - values[k + WID])),
			fv_r[k]);
     }
     //Check for slope consistency. We set it to be the PLM limited slope, if it does not have the same sign
     fd_l[k] = select(slope_sign * fd_l[k] < 0.0, slope_sign * slope_abs, fd_l[k]);
     fd_r[k] = select(slope_sign * fd_r[k] < 0.0, slope_sign * slope_abs, fd_r[k]);
  }
}


/*Filter to make sure that discontinuous edge values are monotonic at
  any edge. If edge values are discontinuous and nonmono- tonic, they
  are both replaced by their average" */
inline void filter_face_value_monotonicity(Vec4 *values, uint n_cblocks, Vec4 *fv_l, Vec4 *fv_r){   
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

/*make sure quartic polynomial is monotonic*/
inline void filter_pqm_monotonicity(Vec4 *values, uint n_cblocks, Vec4 *fv_l, Vec4 *fv_r, Vec4 *fd_l, Vec4 *fd_r){   
  const Vec4 root_outside = Vec4(100.0); //fixed values give to roots clearly outside [0,1], or nonexisting ones*/
  for (int k = 0; k < n_cblocks * WID; k++){
    /*second derivative coefficients, eq 23 in white et al.*/
    Vec4 b0 =   60.0 * values[k + WID] - 24.0 * fv_r[k] - 36.0 * fv_l[k] + 3.0 * (fd_r[k] - 3.0 * fd_l[k]);
    Vec4 b1 = -360.0 * values[k + WID] + 36.0 * fd_l[k] - 24.0 * fd_r[k] + 168.0 * fv_r[k] + 192.0 * fv_l[k];
    Vec4 b2 =  360.0 * values[k + WID] + 30.0 * (fd_r[k] - fd_l[k]) - 180.0 * (fv_l[k] + fv_r[k]);
    /*let's compute sqrt value to be used for computing roots. If we
      take sqrt of negaitve numbers, then we instead set a value that
      will make the root to be +-100 which is well outside range
      of[0,1]. We also guard the sqrt against sqrt with negative
      numbers by doing a max*/
    const Vec4 sqrt_val = select(b1 * b1- 4 * b0 * b2 < 0.0, 
				 b1 + 200.0 * b2,
				 sqrt(max(b1 * b1- 4 * b0 * b2, 0.0))); 
    //compute roots. If b2 is zero, we set root value to be
    //outside. Division is safe with vectorclass in any case
    Vec4 root1 = select(b2 == 0.0, root_outside, (-b1 + sqrt_val) / (2 * b2));
    Vec4 root2 = select(b2 == 0.0, root_outside, (-b1 - sqrt_val) / (2 * b2));
    if( horizontal_and( (root1 < 0.0 || root1 > 1.0 ) &&
			(root2 < 0.0 || root2 > 1.0 )) ) {
      //no inflexion point for any vector component
      continue;
    }
    /*PLM slope, MC limiter*/
    Vec4 plm_slope_l = 2.0 * (values[k + WID] - values[k - 1 + WID]);
    Vec4 plm_slope_r = 2.0 * (values[k + 1 + WID] - values[k + WID]);
    Vec4 slope_sign = select(plm_slope_l + plm_slope_r < 0, -1.0, 1.0);
    /*first derivative coefficients*/
    Vec4 c0 = fd_l[k];
    Vec4 c1 = b0;
    Vec4 c2 = b1 / 2.0;
    Vec4 c3 = b2 / 3.0;
    //compute both slopes at inflexion points, at least one of these
    //is with [0..1]. If the root is not in this range, we
    //simplify later if statements by setting it to the plm slope
    //sign
    Vec4 root1_slope = select(root1 >= 0.0 && root1 <= 1.0, 
			 c0  + c1 * root1 + c2 * root1 * root1 + c3 * root1 * root1 * root1,
			 slope_sign);
    Vec4 root2_slope = select(root2 >= 0.0 && root2 <= 1.0, 
			 c0  + c1 * root2 + c2 * root2 * root2 + c3 * root2 * root2 * root2,
			 slope_sign);
    if (horizontal_or (root1_slope * slope_sign < 0.0 || root2_slope * slope_sign < 0.0 )) {
      //serialized the handling of inflexion points, these do not happen for smooth regions
      for(uint i = 0;i < WID; i++) {
	if(root1_slope[i] * slope_sign[i] < 0.0 || root2_slope[i] * slope_sign[i] < 0.0 ) {
      	//need to collapse, at least one inflexion point has wrong
	//sign.
	if(fabs(plm_slope_l[i]) <= fabs(plm_slope_r[i])) {
	  //collapse to left edge (eq 21)
	  fd_l[k].insert( i, 1.0 / 3.0 * ( 10 * values[k + WID][i] - 2.0 * fv_r[k][i] - 8.0 * fv_l[k][i]));
	  fd_r[k].insert( i, -10.0 * values[k + WID][i] + 6.0 * fv_r[k][i] + 4.0 * fv_l[k][i]);
	  //check if PLM slope is consistent (eq 28 & 29)
	  if (slope_sign[i] * fd_l[k][i] < 0) {
	    fd_l[k].insert( i, 0);
	    fv_r[k].insert( i, 5 * values[k + WID][i] - 4 * fv_l[k][i]);
	    fd_r[k].insert( i, 20 * (values[k + WID][i] - fv_l[k][i]));
	  }
	  else if (slope_sign[i] * fd_r[k][i] < 0) {
	    fd_r[k].insert( i, 0);
	    fv_l[k].insert( i, 0.5 * (5 * values[k + WID][i] - 3 * fv_r[k][i]));
	    fd_l[k].insert( i, 10.0 / 3.0 * (-values[k + WID][i] + fv_r[k][i]));
	  }
	}
	else {
	  //collapse to right edge (eq 21)
	  fd_l[k].insert( i, 10.0 * values[k + WID][i] - 6.0 * fv_l[k][i] - 4.0 * fv_r[k][i]);
	  fd_r[k].insert( i, 1.0 / 3.0 * ( - 10.0 * values[k + WID][i] + 2 * fv_l[k][i] + 8 * fv_r[k][i]));
	  //check if PLM slope is consistent (eq 28 & 29)
	  if (slope_sign[i] * fd_l[k][i] < 0) {
	    fd_l[k].insert( i, 0);
	    fv_r[k].insert( i, 0.5 * ( 5 * values[k + WID][i] - 3 * fv_l[k][i]));
	    fd_r[k].insert( i, 10.0 / 3.0 * (values[k + WID][i] - fv_l[k][i]));
	  }
	  else if (slope_sign[i] * fd_r[k][i] < 0) {
	    fd_r[k].insert( i, 0);
	    fv_l[k].insert( i, 5 * values[k + WID][i] - 4 * fv_r[k][i]);
	    fd_l[k].insert( i, 20.0 * ( - values[k + WID][i] + fv_r[k][i]));
	  }
	}
      }
    }
  }
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

   compute_h6_face_values(values,n_cblocks,fv_l, fv_r); 
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

inline void compute_pqm_coeff_explicit_column(Vec4 *values, uint n_cblocks, Vec4 a[][RECONSTRUCTION_ORDER + 1]){
   Vec4 p_value;
   Vec4 m_value;
   /*next face values and derivatives, extra space to ease
     implementation (final left is not used, but copied to right*/
   Vec4 fv_l[MAX_BLOCKS_PER_DIM * WID + 1]; /*left face value*/
   Vec4 fv_r[MAX_BLOCKS_PER_DIM * WID + 1]; /*right face value*/
   Vec4 fd_l[MAX_BLOCKS_PER_DIM * WID + 1]; /*left face derivative*/
   Vec4 fd_r[MAX_BLOCKS_PER_DIM * WID + 1]; /*right face derivative*/
   compute_h5_face_values(values,n_cblocks,fv_l, fv_r); 
   compute_h5_face_derivatives(values,n_cblocks,fd_l, fd_r); 
   filter_extrema_boundedness(values,n_cblocks, fv_l, fv_r, fd_l, fd_r); 
   filter_face_value_monotonicity(values,n_cblocks, fv_l, fv_r);
   filter_pqm_monotonicity(values,n_cblocks,fv_l, fv_r, fd_l, fd_r); 
   for (uint k = 0; k < n_cblocks * WID; k++){
      //Fit a second order polynomial for reconstruction see, e.g., White
      //2008 (PQM article) (note additional integration factors built in,
      //contrary to White (2008) eq. 4
      a[k][0] = fv_l[k];
      a[k][1] = fd_l[k]/2.0;
      a[k][2] = 10 * values[k + WID] - 4.0 * fv_r[k] - 6.0 * fv_l[k] + 0.5 * (fd_r[k] - 3 * fd_l[k]);
      a[k][3] = -15 * values[k + WID]  + 1.5 * fd_l[k] - fd_r[k] + 7.0 * fv_r[k] + 8 * fv_l[k];
      a[k][4] = 6.0* values[k + WID] +  0.5 * (fd_r[k] - fd_l[k]) - 3.0 * (fv_l[k] + fv_r[k]);
   }
}

#endif
