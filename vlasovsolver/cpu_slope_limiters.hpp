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

#ifndef HOSTDEV_SLOPE_LIMITERS_H
#define HOSTDEV_SLOPE_LIMITERS_H

#include "vec.h"
#include "../arch/arch_device_api.h"

using namespace std;

static inline Vec minmod(const Vec &slope1, const Vec &slope2)
{
   const Vec veczero(0.0);
   const Vec slope = select(abs(slope1) < abs(slope2), slope1, slope2);
   return select(slope1 * slope2 <= 0, veczero, slope);
}
static inline Vec maxmod(const Vec &slope1, const Vec &slope2)
{
   const Vec veczero(0.0);
   const Vec slope = select(abs(slope1) > abs(slope2), slope1, slope2);
   return select(slope1 * slope2 <= 0, veczero, slope);
}

/*!
  Superbee slope limiter
*/

static inline Vec slope_limiter_sb(const Vec &l, const Vec &m, const Vec &r)
{
   const Vec a = r-m;
   const Vec b = m-l;
   const Vec slope1 = minmod(a, 2*b);
   const Vec slope2 = minmod(2*a, b);
   return maxmod(slope1, slope2);
}

/*!
  Minmod slope limiter
*/

static inline Vec slope_limiter_minmod(const Vec& l,const Vec& m, const Vec& r)
{
   //Vec sign;
   const Vec a=r-m;
   const Vec b=m-l;
   return minmod(a,b);
}

/*!
  MC slope limiter
*/

static inline Vec slope_limiter_mc(const Vec& l,const Vec& m, const Vec& r)
{
   const Vec veczero(0.0);
   const Vec vectwo(2.0);
   const Vec vechalf(0.5);
   //Vec sign;
   const Vec a=r-m;
   const Vec b=m-l;
   Vec minval=min(vectwo*abs(a),vectwo*abs(b));
   minval=min(minval,vechalf*abs(a+b));

   //check for extrema
   const Vec output = select(a*b < 0,veczero,minval);
   //set sign
   return select(a + b < 0,-output,output);
}

static inline Vec slope_limiter_minmod_amr(const Vec& l,const Vec& m, const Vec& r,const Vec& a,const Vec& b)
{
   const Vec J = r-l;
   Vec f = (m-l)/J;
   f = min(Vec(1.0),f);
   return min(f/(1+a),(Vec(1.)-f)/(1+b))*2*J;
}

static inline Vec slope_limiter(const Vec &l, const Vec &m, const Vec &r)
{
   return slope_limiter_sb(l,m,r);
   //return slope_limiter_minmod(l,m,r);
}

/*
 * @param a Cell size fraction dx[i-1]/dx[i] = 1/2, 1, or 2.
 * @param b Cell size fraction dx[i+1]/dx[i] = 1/2, 1, or 2.
 * @return Limited value of slope.*/
static inline Vec slope_limiter_amr(const Vec& l,const Vec& m, const Vec& r,const Vec& dx_left,const Vec& dx_rght)
{
   return slope_limiter_minmod_amr(l,m,r,dx_left,dx_rght);
}

/* Slope limiter with abs and sign separatelym, uses the currently active slope limiter*/
static inline void slope_limiter(const Vec& l,const Vec& m, const Vec& r, Vec& slope_abs, Vec& slope_sign)
{
   const Vec slope = slope_limiter(l,m,r);
   slope_abs = abs(slope);
   slope_sign = select(slope > 0, Vec(1.0), Vec(-1.0));
}


/****
     Define functions for Realf instead of Vec
     These functions take direct arguments instead of references for non-const input
     so that registers and other optimizations are available.
***/

static ARCH_DEV inline Realf minmod(const Realf slope1, const Realf slope2)
{
   const Realf slope = (abs(slope1) < abs(slope2)) ? slope1 : slope2;
   return (slope1 * slope2 <= 0) ? 0 : slope;
}
static ARCH_DEV inline Realf maxmod(const Realf slope1, const Realf slope2)
{
   const Realf slope = (abs(slope1) > abs(slope2)) ? slope1 : slope2;
   return (slope1 * slope2 <= 0) ? 0 : slope;
}

/*!
  Superbee slope limiter
*/

static ARCH_DEV inline Realf slope_limiter_sb(const Realf l, const Realf m, const Realf r)
{
   const Realf a = r-m;
   const Realf b = m-l;
   const Realf slope1 = minmod(a, 2*b);
   const Realf slope2 = minmod(2*a, b);
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
   Realf minval=min(2*abs(a),2*abs(b));
   minval=min(minval,(Realf)0.5*abs(a+b));

   //check for extrema
   const Realf output = (a*b < 0) ? 0 : minval;
   //set sign
   return (a + b < 0) ? -output : output;
}

static ARCH_DEV inline Realf slope_limiter_minmod_amr(const Realf l,const Realf m, const Realf r,const Realf a,const Realf b)
{
   const Realf J = r-l;
   Realf f = (m-l)/J;
   f = min((Realf)1.0,f);
   return min((Realf)f/(1+a),(Realf)(1.-f)/(1+b))*2*J;
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
   slope_sign = (slope > 0) ? 1 : -1.0;
}



#endif
