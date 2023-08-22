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

#ifndef HOSTDEV_1D_PLM_H
#define HOSTDEV_1D_PLM_H

#include "vec.h"
#include "../arch/arch_device_api.h"
#include "cpu_slope_limiters.hpp"

using namespace std;

/*!
 Compute PLM coefficients
 f(v) = a[0] + a[1]/2.0*t 
t=(v-v_{i-0.5})/dv where v_{i-0.5} is the left face of a cell
The factor 2.0 is in the polynom to ease integration, then integral is a[0]*t + a[1]*t**2
*/

static ARCH_HOSTDEV inline void compute_plm_coeff(const Vec * const values, uint k, Vec a[2], const Realv threshold)
{
  // scale values closer to 1 for more accurate slope limiter calculation
  const Realv scale = 1./threshold;
  //Vec v_1 = values[k - 1] * scale;
  //Vec v_2 = values[k] * scale;
  //Vec v_3 = values[k + 1] * scale;
  //Vec d_cv = slope_limiter(v_1, v_2, v_3) * threshold;
  const Vec d_cv = slope_limiter( values[k-1]*scale, values[k]*scale, values[k+1]*scale)*threshold;
  a[0] = values[k] - d_cv * 0.5;
  a[1] = d_cv * 0.5;
}

/**** 
      Define functions for Realf instead of Vec 
***/

static ARCH_DEV inline void compute_plm_coeff(const Vec* const values, uint k, Realf a[2], const Realv threshold, const int index)
{
  // scale values closer to 1 for more accurate slope limiter calculation
  const Realv scale = 1./threshold;
  //Vec v_1 = values[k - 1] * scale;
  //Vec v_2 = values[k] * scale;
  //Vec v_3 = values[k + 1] * scale;
  //Vec d_cv = slope_limiter(v_1, v_2, v_3) * threshold;
  const Realf d_cv = slope_limiter( values[k-1][index]*scale, values[k][index]*scale, values[k+1][index]*scale)*threshold;
  a[0] = values[k][index] - d_cv * 0.5;
  a[1] = d_cv * 0.5;
}


#endif
