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

#ifndef GPU_SLOPE_LIMITERS_H
#define GPU_SLOPE_LIMITERS_H

#include "../arch/arch_device_api.h"

/****
     Define functions for Realf instead of Vec
     These functions take direct arguments instead of references for non-const input
     so that registers and other optimizations are available.
***/

static ARCH_DEV inline Realf minmod(const Realf slope1, const Realf slope2)
{
   const Realf slope = (abs(slope1) < abs(slope2)) ? slope1 : slope2;
   return (slope1 * slope2 <= (Realf)(0.0)) ? (Realf)(0.0) : slope;
}
static ARCH_DEV inline Realf maxmod(const Realf slope1, const Realf slope2)
{
   const Realf slope = (abs(slope1) > abs(slope2)) ? slope1 : slope2;
   return (slope1 * slope2 <= (Realf)(0.0)) ? (Realf)(0.0) : slope;
}

/*!
  Superbee slope limiter
*/

static ARCH_DEV inline Realf slope_limiter_sb(const Realf l, const Realf m, const Realf r)
{
   const Realf a = r-m;
   const Realf b = m-l;
   const Realf slope1 = minmod(a, (Realf)(2.0)*b);
   const Realf slope2 = minmod((Realf)(2.0)*a, b);
   return maxmod(slope1, slope2);
}

/*!
  Minmod slope limiter
*/

static ARCH_DEV inline Realf slope_limiter_minmod(const Realf l, const Realf m, const Realf r)
{
   const Realf a=r-m;
   const Realf b=m-l;
   return minmod(a,b);
}

/*!
  MC slope limiter
*/

static ARCH_DEV inline Realf slope_limiter_mc(const Realf l, const Realf m, const Realf r)
{
   const Realf a=r-m;
   const Realf b=m-l;
   Realf minval=min((Realf)(2.0)*abs(a),(Realf)(2.0)*abs(b));
   minval=min(minval,(Realf)(0.5)*abs(a+b));

   //check for extrema
   const Realf output = (a*b < (Realf)(0.0)) ? (Realf)(0.0) : minval;
   //set sign
   return (a + b < (Realf)(0.0)) ? -output : output;
}

static ARCH_DEV inline Realf slope_limiter_minmod_amr(const Realf l,const Realf m, const Realf r,const Realf a,const Realf b)
{
   const Realf J = r-l;
   Realf f = (m-l)/J;
   f = min((Realf)(1.0),f);
   return min((Realf)(f)/(1+a),(Realf)(1.-f)/((Realf)(1.0)+b))*(Realf)(2.0)*J;
}

static ARCH_DEV inline Realf slope_limiter(const Realf l, const Realf m, const Realf r)
{
   return slope_limiter_sb(l,m,r);
   //return slope_limiter_minmod(l,m,r);
}

/*
 * @param a Cell size fraction dx[i-1]/dx[i] = 1/2, 1, or 2.
 * @param b Cell size fraction dx[i+1]/dx[i] = 1/2, 1, or 2.
 * @return Limited value of slope.*/
static ARCH_DEV inline Realf slope_limiter_amr(const Realf l,const Realf m, const Realf r,const Realf dx_left,const Realf dx_rght)
{
   return slope_limiter_minmod_amr(l,m,r,dx_left,dx_rght);
}

/* Slope limiter with abs and sign separately, uses the currently active slope limiter*/
static ARCH_DEV inline void slope_limiter(const Realf l,const Realf m, const Realf r, Realf& slope_abs, Realf& slope_sign)
{
   const Realf slope = slope_limiter(l,m,r);
   slope_abs = abs(slope);
   slope_sign = (slope > (Realf)(0.0)) ? (Realf)(1.0) : (Realf)(-1.0);
}



#endif
