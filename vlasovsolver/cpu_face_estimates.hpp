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
enum face_value_estimate {h4, h5, h6};
enum face_derivative_estimate {dh4, dh5 };

/*FIXME & NOTE. MAX_IH_CELLS should be in practice (MAX_BLOCKS_PER_DIM * WID), so 400
  supports up to 100 blocks. This can be too little at some point of
  time and no checks are done for bounds...*/
#define MAX_IH_CELLS 400

/*!
 * Simple tridiagonal solver, based on NR.
 */

void tridiag(uint n,const Vec4 &a,const Vec4 &b,const Vec4 &c,Vec4 *u, Vec4 *r){
  Vec4 gam[MAX_IH_CELLS];
  
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



/*! 
  Compute left face value based on the explicit h4 estimate.

  Right face value can be obtained as left face value of cell i + 1.
  
  \param avgs Array with volume averages. It is assumed a large enough stencil is defined around i.
  \param i Index of cell in avgs for which the left face is computed
  \param fv_l Face value on left face of cell i
*/ 
inline void compute_h4_left_face_value(const Vec4 * const avgs, uint k, Vec4 &fv_l){   
  /*compute left value*/
  fv_l = 1.0/12.0 * ( - 1.0 * avgs[k - 2]  
		      + 7.0 * avgs[k - 1] 
		      + 7.0 * avgs[k] 
		      - 1.0 * avgs[k + 1]);
}


/*! 
  Compute left and right face value based on the explicit h5 estimate.

  \param avgs Array with volume averages. It is assumed a large enough stencil is defined around i.
  \param i Index of cell in avgs for which the left face is computed
  \param fv_l Face value on left face of cell i
*/ 
inline void compute_h5_face_values(const Vec4 * const avgs, uint k, Vec4 &fv_l, Vec4 &fv_r){   
  /*compute left values*/
  fv_l = 1.0/60.0 * (- 3.0 * avgs[k - 2]  
			+ 27.0 * avgs[k - 1] 
			+ 47.0 * avgs[k ] 
			- 13.0 * avgs[k + 1] 
			+ 2.0 * avgs[k + 2]);
  fv_r = 1.0/60.0 * ( 2.0 * avgs[k - 2] 
			 - 13.0 * avgs[k - 1] 
			 + 47.0 * avgs[k]
			 + 27.0 * avgs[k + 1] 
			 - 3.0 * avgs[k + 2]);
}

/*! 
  Compute left face value based on the explicit h6 estimate.

  Right face value can be obtained as left face value of cell i + 1.
  
  \param avgs Array with volume averages. It is assumed a large enough stencil is defined around i.
  \param i Index of cell in avgs for which the left face is computed
  \param fv_l Face value on left face of cell i
*/ 
inline void compute_h6_left_face_value(const Vec4 * const avgs, uint k, Vec4 &fv_l){   
  /*compute left value*/
  fv_l = 1.0/60.0 * (avgs[k - 3]  
			- 8.0 * avgs[k - 2]  
			+ 37.0 * avgs[k - 1] 
			+ 37.0 * avgs[k ] 
			- 8.0 * avgs[k + 1] 
			+ avgs[k + 2]);
}



/*! 
  Compute left face derivative based on the explicit h4 estimate.

  Right face derivative can be obtained as left face derivative of cell i + 1.
  
  \param avgs Array with volume averages. It is assumed a large enough stencil is defined around i.
  \param i Index of cell in avgs for which the left face derivativeis computed
  \param fd_l Face derivative on left face of cell i
*/ 
inline void compute_h4_left_face_derivative(const Vec4 * const avgs, uint k, Vec4 &fd_l){   
  fd_l = 1.0/12.0 * (15.0 * (avgs[k] - avgs[k - 1]) - (avgs[k + 1] - avgs[k - 2]));
}


/*! 
  Compute left face derivative based on the explicit h5 estimate.

  Right face derivative can be obtained as left face derivative of cell i + 1.
  
  \param avgs Array with volume averages. It is assumed a large enough stencil is defined around i.
  \param i Index of cell in avgs for which the left face derivativeis computed
  \param fd_l Face derivative on left face of cell i
*/ 
inline void compute_h5_left_face_derivative(const Vec4 * const avgs, uint k, Vec4 &fd_l){   
  fd_l = 1.0/180.0 * (245 * (avgs[k] - avgs[k - 1])  
		      - 25 * (avgs[k + 1] - avgs[k - 2]) 
		      + 2 * (avgs[k + 2] - avgs[k - 3]));
}



/*! 
  Compute all face values based on the implicit ih4 estimate.

  Based on White et al. 2008.
  
  \param avgs Array with volume averages (size is elements + 2 * avgs_offset).
  \param elements Number of cells
  \param avgs_offset This is the amount of extra boundary cells in avgs. Cell i is at i + avgs_offset.
  \param fv_l Face values on left face of cell (size is elements)
  \param fv_r Face values on right face of cell (size is elements)
*/ 


inline void compute_ih4_face_values(Vec4 *avgs, uint elements, uint avgs_offset,  
				    Vec4 *fv_l, Vec4 *fv_r){
  Vec4 r[MAX_IH_CELLS + 1];
  const Vec4 a(1.0/4.0);
  const Vec4 b(1.0);
  const Vec4 c(1.0/4.0);
  
  for (uint k = 0; k < elements + 1; k++){
    //Note that these are now defined for left face values, so not exactly like in white
    r[k] = ( 3.0 * avgs[k - 1 + avgs_offset] + 3.0 * avgs[k + avgs_offset])/4.0;
  }
  
  tridiag(elements + 1,a,b,c,fv_l,r);

  /*copy right face values */
  for (uint k = 1; k < elements + 1; k++){
    fv_r[k-1] = fv_l[k];
  }

}


/*! 
  Compute all face values based on the implicit ih6 estimate.

  Based on White et al. 2008. Note that the value for a and c is wrong
  in the paper. The original source is: [ Lacor, Chris, Sergey
  Smirnov, and Martine Baelmans. “A Finite Volume Formulation of
  Compact Central Schemes on Arbitrary Structured Grids.” Journal of
  Computational Physics 198, no. 2 (August 2004):
  535–66. doi:10.1016/j.jcp.2004.01.025.]

  and solving for the variables there we get for the variable names in
  that paper (sympy): 
  solve([Eq(a + b, 1 + 2*alpha), Eq(a+(2**3 - 1)*b,6*alpha), Eq(a + (2**5 -1)*b,10*alpha)],a,alpha,b)  
  {alpha: 1/3, b: 1/18, a: 29/18}

  
  \param avgs Array with volume averages (size is elements + 2 * avgs_offset).
  \param elements Number of cells
  \param avgs_offset This is the amount of extra boundary cells in avgs. Cell i is at i + avgs_offset.
  \param fv_l Face values on left face of cell (size is elements)
  \param fv_r Face values on right face of cell (size is elements)
*/ 

inline void compute_ih6_face_values(Vec4 *avgs, uint elements, uint avgs_offset,  
				    Vec4 *fv_l, Vec4 *fv_r){
  Vec4 r[MAX_IH_CELLS + 1];
  const Vec4 a(1.0/3.0);
  const Vec4 b(1.0);
  const Vec4 c(1.0/3.0);
  
  for (uint k = 0; k < elements + 1; k++){
    //Note that these are now defined for left face values, so not exactly like in white
    r[k] = (avgs[k - 2 + avgs_offset] + 
	    29.0 * avgs[k - 1 + avgs_offset] + 
	    29.0 * avgs[k + avgs_offset] + 
	    avgs[k + 1 + avgs_offset])/36.0;
  }
  
  tridiag(elements + 1, a, b, c, fv_l, r);

  /*copy right face values */
  for (uint k = 1; k < elements + 1; k++){
    fv_r[k-1] = fv_l[k];
  }
}



#endif
