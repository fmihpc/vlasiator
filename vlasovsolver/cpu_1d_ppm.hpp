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

#ifndef HOSTDEV_1D_PPM_H
#define HOSTDEV_1D_PPM_H

#include "vec.h"
#include "../arch/arch_device_api.h"
#include "cpu_slope_limiters.hpp"
#include "cpu_face_estimates.hpp"

using namespace std;

/*
  Compute parabolic reconstruction with an explicit scheme
*/
static inline void compute_ppm_coeff(const Vec * const values, face_estimate_order order, uint k, Vec a[3], const Realf threshold)
{
   Vec m_face; /*left face value*/
   Vec p_face; /*right face value*/
   compute_filtered_face_values(values, k, order, m_face, p_face, threshold);
   //Coella et al, check for monotonicity
   const Vec one_sixth(1.0/6.0);
   m_face = select((p_face - m_face) * (values[k] - 0.5 * (m_face + p_face)) >
                   (p_face - m_face) * (p_face - m_face) * one_sixth,
                   3 * values[k] - 2 * p_face, m_face);
   p_face = select(-(p_face - m_face) * (p_face - m_face) * one_sixth >
                   (p_face - m_face) * (values[k] - 0.5 * (m_face + p_face)),
                   3 * values[k] - 2 * m_face, p_face);
   //Fit a second order polynomial for reconstruction see, e.g., White
   //2008 (PQM article) (note additional integration factors built in,
   //contrary to White (2008) eq. 4
   a[0] = m_face;
   a[1] = 3.0 * values[k] - 2.0 * m_face - p_face;
   a[2] = (m_face + p_face - 2.0 * values[k]);
}

/****
      Define functions for Realf instead of Vec
***/

static ARCH_DEV inline void compute_ppm_coeff(const Vec* __restrict__ const values, face_estimate_order order, uint k, Realf a[3], const Realf threshold, const int index)
{
   Realf m_face; /*left face value*/
   Realf p_face; /*right face value*/
   compute_filtered_face_values(values, k, order, m_face, p_face, threshold, index);
   //Coella et al, check for monotonicity
   const Realf one_sixth(1.0/6.0);
   m_face = ((p_face - m_face) * (values[k][index] - 0.5 * (m_face + p_face)) >
             (p_face - m_face) * (p_face - m_face) * one_sixth) ?
      3 * values[k][index] - 2 * p_face : m_face;
   p_face = (-(p_face - m_face) * (p_face - m_face) * one_sixth >
             (p_face - m_face) * (values[k][index] - 0.5 * (m_face + p_face))) ?
      3 * values[k][index] - 2 * m_face : p_face;
   //Fit a second order polynomial for reconstruction see, e.g., White
   //2008 (PQM article) (note additional integration factors built in,
   //contrary to White (2008) eq. 4
   a[0] = m_face;
   a[1] = 3.0 * values[k][index] - 2.0 * m_face - p_face;
   a[2] = (m_face + p_face - 2.0 * values[k][index]);
}

#endif
