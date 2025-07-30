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

#ifndef GPU_1D_PPM_NU_H
#define GPU_1D_PPM_NU_H

#include <iostream>
#include "algorithm"
#include "cmath"

#include "../arch/arch_device_api.h"
#include "gpu_slope_limiters.hpp"
#include "gpu_face_estimates.hpp"
#include "../definitions.h"


/****
     Define functions for Realf instead of Vec
***/

ARCH_DEV inline void compute_ppm_coeff_nonuniform(const Realf* __restrict__ const dv, const Realf* __restrict__ const values, face_estimate_order order, int k, Realf a[3], const Realf threshold, const int index, const int stride){
   Realf m_face; /*left face value*/
   Realf p_face; /*right face value*/
   compute_filtered_face_values_nonuniform(dv, values, k, order, m_face, p_face, threshold, index, stride);

   //Coella et al, check for monotonicity
   m_face = ((p_face - m_face) * (values[k*stride+index] - (Realf)(0.5) * (m_face + p_face)) >
             (p_face - m_face)*(p_face - m_face) * ((Realf)(1.0)/(Realf)(6.0))) ?
      (Realf)(3.0) * values[k*stride+index] - (Realf)(2.0) * p_face :
      m_face;
   p_face = (-(p_face - m_face) * (p_face - m_face) * ((Realf)(1.0)/(Realf)(6.0))) >
      (p_face - m_face) * (values[k*stride+index] - (Realf)(0.5) * (m_face + p_face)) ?
      (Realf)(3.0) * values[k*stride+index] - (Realf)(2.0) * m_face :
      p_face;

   //Fit a second order polynomial for reconstruction see, e.g., White
   //2008 (PQM article) (note additional integration factors built in,
   //contrary to White (2008) eq. 4
   a[0] = m_face;
   a[1] = (Realf)(3.0) * values[k*stride+index] - (Realf)(2.0) * m_face - p_face;
   a[2] = (m_face + p_face - (Realf)(2.0) * values[k*stride+index]);
}

#endif
