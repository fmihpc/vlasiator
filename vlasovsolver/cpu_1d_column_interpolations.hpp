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

inline void compute_plm_coeff_explicit_column(Vec4 *values, uint n_cblocks, uint j, Vec4 a[][RECONSTRUCTION_ORDER + 1]){
   values += i_pcolumnv(n_cblocks, -1, j, 0);
   for (uint k = 0; k < n_cblocks * WID; k++){   
      const Vec4 d_cv=slope_limiter(values[k - 1 + WID], values[k + WID], values[k + 1 + WID]);
      a[k][0] = values[k + WID] - d_cv * 0.5;
      a[k][1] = d_cv * 0.5;
   }
}




inline void compute_ppm_coeff_explicit_column(Vec4 *values, uint n_cblocks, uint j, Vec4 a[][RECONSTRUCTION_ORDER + 1]){
   Vec4 p_face;
   Vec4 m_face;
   Vec4 p_face_unfiltered;

   /*move pointer to start where this column starts (for fixed j, all
     k values are in a consecutative order) We start from block at -1
     (this much extra space we have for stencils). Thus values[k +
     WID] corresponds to the current (centered) cell.
   */
   values += i_pcolumnv(n_cblocks, -1, j, 0);
   /*set up shifting operation for face value. Compute right face for
    * cell -1, this will then be shifted to left cell*/
   int k = -1;
   p_face_unfiltered = 1.0/60.0 * (values[k - 2 + WID]  - 8.0 * values[k -1 + WID]  + 37.0 * values[k + WID] +
                                   37.0 * values[k + 1 + WID] - 8.0 * values[k + 2 + WID] + values[k + 3 + WID]);   
   for (k = 0; k < n_cblocks * WID; k++){
      //white 08 H6 face estimates, better than H5. Shift old unfilitered value at upper edge to the lower edge (identical edge)
      m_face = p_face_unfiltered;
      p_face_unfiltered = 1.0/60.0 * (values[k - 2 + WID]  - 8.0 * values[k -1 + WID]  + 37.0 * values[k + WID] +
                                      37.0 * values[k + 1 + WID] - 8.0 * values[k + 2 + WID] + values[k + 3 + WID]);
      p_face = p_face_unfiltered;
         
      bool fix_bounds = horizontal_or((values[k -1 + WID] - m_face) * (m_face - values[k + WID]) < 0 ||
                                      (values[k + 1 + WID] - p_face) * (p_face - values[k + WID]) < 0);
      if(fix_bounds) {
         Vec4 slope_abs,slope_sign;
         slope_limiter(values[k -1 + WID], values[k + WID], values[k + 1 + WID], slope_abs, slope_sign);
         //detect and        fix boundedness, as in WHITE 2008
         m_face = select((values[k -1 + WID] - m_face) * (m_face - values[k + WID]) < 0,
                         values[k + WID] - slope_sign * min( 0.5 * slope_abs, abs(m_face - values[k + WID])),
                         m_face);
         p_face = select((values[k + 1 + WID] - p_face) * (p_face - values[k + WID]) < 0,
                         values[k + WID] + slope_sign * min( 0.5 * slope_abs, abs(p_face - values[k + WID])),
                         p_face);
      }
   
      //Coella1984 eq. 1.10, detect extrema and make algorithm constant if it is
      Vec4 extrema_check = ((p_face - values[k + WID]) * (values[k + WID] - m_face));
      m_face = select(extrema_check < 0, values[k + WID], m_face);
      p_face = select(extrema_check < 0, values[k + WID], p_face);
   
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

// inline void compute_h6_face_values(Real *values, uint n_cblocks, uint j, face_values fv[]){   
//    Vec4 p_face;
//    Vec4 m_face;
//    Vec4 p_face_deriv;
//    Vec4 m_face_deriv;
   
//    Vec4 p_face_unfiltered;
//    Vec4 p_face_deriv_unfiltered;
   
//    Vec4 values[k - 2 + WID],mv,values[k + WID],values[k + 1 + WID],values[k + 2 + WID],values[k + 3 + WID]; /* values in z-direction*/
//    /*set up shifting operation. These will be shifted to values[k - 2 + WID] - values[k + 1 + WID]*/
//    values[k - 2 + WID].load(values + i_pcolumnv(n_cblocks, 0, j, -3));
//    mv.load(values + i_pcolumnv(n_cblocks, 0, j, -2));
//    values[k + WID].load(values + i_pcolumnv(n_cblocks, 0, j, -1));
//    values[k + 1 + WID].load(values + i_pcolumnv(n_cblocks, 0, j, 0));
//    values[k + 2 + WID].load(values + i_pcolumnv(n_cblocks, 0, j, 1));
//    values[k + 1 + WID]alues[k + 2 + WID].load(values + i_pcolumnv(n_cblocks, 0, j, 2));
//    p_face_unfiltered = 1.0/60.0 * (values[k - 2 + WID]  - 8.0 * mv  + 37.0 * values[k + WID] + 37.0 * values[k + 1 + WID] - 8.0 * values[k + 2 + WID] + values[k + 1 + WID]alues[k + 2 + WID]);   
//    p_face_deriv_unfiltered = 1.0 / 180.0 * ( 245.0 * (values[k + 1 + WID] - values[k + WID]) - 25.0 * (values[k + 2 + WID] - mv) + 2.0 * (values[k + 1 + WID]alues[k + 2 + WID] - values[k - 2 + WID]));
   
//    for (uint block_i = 0; block_i < n_cblocks; block_i++){
//       for (uint k = 0; k < WID; ++k){
//          /*shift values*/
//          values[k - 2 + WID] = mv;
//          mv = values[k + WID];
//          values[k + WID] = values[k + 1 + WID];
//          values[k + 1 + WID] = values[k + 2 + WID];
//          values[k + 2 + WID] = values[k + 1 + WID]alues[k + 2 + WID];
//          values[k + 1 + WID]alues[k + 2 + WID].load(values + i_pcolumnv(n_cblocks, block_i, j, k + 3));
//          //white 08 H6 face estimates, better than H5. Shift old unfilitered value at upper edge to the lower edge (identical edge)
//          m_face = p_face_unfiltered;
//          p_face_unfiltered = 1.0/60.0 * (values[k - 2 + WID]  - 8.0 * mv  + 37.0 * values[k + WID] + 37.0 * values[k + 1 + WID] - 8.0 * values[k + 2 + WID] + values[k + 1 + WID]alues[k + 2 + WID]);
//          p_face = p_face_unfiltered;
   
//       }
//    }
// }

// /*
//   PQM reconstruction as published in:

//   White, Laurent, and Alistair Adcroft. “A High-Order Finite Volume Remapping Scheme for Nonuniform Grids: The Piecewise Quartic Method (PQM).” Journal of Computational Physics 227, no. 15 (July 2008): 7394–7422. doi:10.1016/j.jcp.2008.04.026.

// */

// inline void compute_pqm_coeff_explicit_column(Real *values, uint n_cblocks, uint j, Vec4 a[][RECONSTRUCTION_ORDER + 1]){
   
//    Vec4 p_face;
//    Vec4 m_face;
//    Vec4 p_face_deriv;
//    Vec4 m_face_deriv;
   
//    Vec4 p_face_unfiltered;
//    Vec4 p_face_deriv_unfiltered;
   
//    Vec4 values[k - 2 + WID],mv,values[k + WID],values[k + 1 + WID],values[k + 2 + WID],values[k + 1 + WID]alues[k + 2 + WID]; /* values in z-direction*/
//    /*set up shifting operation. These will be shifted to values[k - 2 + WID] - values[k + 1 + WID]*/
//    values[k - 2 + WID].load(values + i_pcolumnv(n_cblocks, 0, j, -3));
//    mv.load(values + i_pcolumnv(n_cblocks, 0, j, -2));
//    values[k + WID].load(values + i_pcolumnv(n_cblocks, 0, j, -1));
//    values[k + 1 + WID].load(values + i_pcolumnv(n_cblocks, 0, j, 0));
//    values[k + 2 + WID].load(values + i_pcolumnv(n_cblocks, 0, j, 1));
//    values[k + 1 + WID]alues[k + 2 + WID].load(values + i_pcolumnv(n_cblocks, 0, j, 2));
//    p_face_unfiltered = 1.0/60.0 * (values[k - 2 + WID]  - 8.0 * mv  + 37.0 * values[k + WID] + 37.0 * values[k + 1 + WID] - 8.0 * values[k + 2 + WID] + values[k + 1 + WID]alues[k + 2 + WID]);   
//    p_face_deriv_unfiltered = 1.0 / 180.0 * ( 245.0 * (values[k + 1 + WID] - values[k + WID]) - 25.0 * (values[k + 2 + WID] - mv) + 2.0 * (values[k + 1 + WID]alues[k + 2 + WID] - values[k - 2 + WID]));
   
//    for (uint block_i = 0; block_i < n_cblocks; block_i++){
//       for (uint k = 0; k < WID; ++k){
//          /*shift values*/
//          values[k - 2 + WID] = mv;
//          mv = values[k + WID];
//          values[k + WID] = values[k + 1 + WID];
//          values[k + 1 + WID] = values[k + 2 + WID];
//          values[k + 2 + WID] = values[k + 1 + WID]alues[k + 2 + WID];
//          values[k + 1 + WID]alues[k + 2 + WID].load(values + i_pcolumnv(n_cblocks, block_i, j, k + 3));
//          //white 08 H6 face estimates, better than H5. Shift old unfilitered value at upper edge to the lower edge (identical edge)
//          m_face = p_face_unfiltered;
//          p_face_unfiltered = 1.0/60.0 * (values[k - 2 + WID]  - 8.0 * mv  + 37.0 * values[k + WID] + 37.0 * values[k + 1 + WID] - 8.0 * values[k + 2 + WID] + values[k + 1 + WID]alues[k + 2 + WID]);
//          p_face = p_face_unfiltered;
//          //white-08 H5 estimates for derivatives at cell faces. Same shifting as for face values
//          m_face_deriv = p_face_deriv_unfiltered;
//          p_face_deriv_unfiltered = 1.0 / 180.0 * ( 245.0 * (values[k + 1 + WID] - values[k + WID]) - 25.0 * (values[k + 2 + WID] - mv) + 2.0 * (values[k + 1 + WID]alues[k + 2 + WID] - values[k - 2 + WID]));
//          p_face_deriv = p_face_deriv_unfiltered;

//          /*now we limit the face values and derivatives*/
         
//          /* Edge values boundedness - Eq. 19 in white-08*/
//          bool fix_bounds = horizontal_or((mv - m_face) * (m_face - values[k + WID]) < 0 || (values[k + 1 + WID] - p_face) * (p_face - values[k + WID]) < 0);
//          if(fix_bounds) {
//             Vec4 slope_abs,slope_sign;
//             slope_limiter(mv, values[k + WID], values[k + 1 + WID], slope_abs, slope_sign);
//             //detect and fix boundedness, as in WHITE 2008
//             m_face = select((mv - m_face) * (m_face - values[k + WID]) < 0,
//                             values[k + WID] - slope_sign * min( 0.5 * slope_abs, abs(m_face - values[k + WID])),
//                             m_face);
//             p_face = select((values[k + 1 + WID] - p_face) * (p_face - values[k + WID]) < 0,
//                             values[k + WID] + slope_sign * min( 0.5 * slope_abs, abs(p_face - values[k + WID])),
//                             p_face);
//          }
         


         
//          //Fit a second order polynomial for reconstruction see, e.g., White
//          //2008 (PQM article) (note additional integration factors built in,
//          //contrary to White (2008) eq. 4
//          a[block_i * WID + k][0] = m_face;
//          a[block_i * WID + k][1] = 3.0 * values[k + WID] - 2.0 * m_face - p_face;
//          a[block_i * WID + k][2] = (m_face + p_face - 2.0 * values[k + WID]);
//       }
//    }
// }

#endif
