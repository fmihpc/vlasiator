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

#ifndef GPU_1D_PLM_H
#define GPU_1D_PLM_H

#include "../arch/arch_device_api.h"
#include "gpu_slope_limiters.hpp"

using namespace std;

/*!
 Compute PLM coefficients
 f(v) = a[0] + a[1]/2.0*t 
t=(v-v_{i-0.5})/dv where v_{i-0.5} is the left face of a cell
The factor 2.0 is in the polynom to ease integration, then integral is a[0]*t + a[1]*t**2
*/

static ARCH_DEV inline void compute_plm_coeff(const Realf* __restrict__ const values, int k, Realf a[2], const Realf threshold, const int index, const int stride)
{
  // scale values closer to 1 for more accurate slope limiter calculation
  const Realf scale = (Realf)(1.0)/threshold;
  const Realf d_cv = slope_limiter( values[(k-1)*stride+index]*scale, values[k*stride+index]*scale, values[(k+1)*stride+index]*scale)*threshold;
  a[0] = values[k*stride+index] - d_cv * (Realf)(0.5);
  a[1] = d_cv * (Realf)(0.5);
}


#endif
