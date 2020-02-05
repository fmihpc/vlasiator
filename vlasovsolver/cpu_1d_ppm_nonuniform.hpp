/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://www.physics.helsinki.fi/vlasiator/
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#ifndef CPU_1D_PPM_NU_H
#define CPU_1D_PPM_NU_H

#include <iostream>
#include "vec.h"
#include "algorithm"
#include "cmath"
#include "cpu_slope_limiters.hpp"
#include "cpu_face_estimates.hpp"

using namespace std;

/*
  Compute parabolic reconstruction with an explicit scheme
*/
inline void compute_ppm_coeff_nonuniform(const Vec * const dv, const Vec * const values, face_estimate_order order, uint k, Vec a[3], const Realv threshold){
   Vec fv_l; /*left face value*/
   Vec fv_r; /*right face value*/
   compute_filtered_face_values_nonuniform(dv, values, k, order, fv_l, fv_r, threshold); 
   
   //Coella et al, check for monotonicity   
   Vec m_face = fv_l;
   Vec p_face = fv_r;

   //std::cout << "value = " << values[k][0] << ", m_face = " << m_face[0] << ", p_face = " << p_face[0] << "\n";
   //std::cout << values[k][0] - m_face[0] << ", " << values[k][0] - p_face[0] << "\n";
   
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

   //std::cout << "value = " << values[k][0] << ", m_face = " << m_face[0] << ", p_face = " << p_face[0] << "\n";
   //std::cout << values[k][0] - m_face[0] << ", " << values[k][0] - p_face[0] << "\n";

   //std::cout << values[k][0] << " " << m_face[0] << " " << p_face[0] << "\n";
}

#endif
