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

#ifndef CPU_SLOPE_LIMITERS_H
#define CPU_SLOPE_LIMITERS_H

#include "vec.h"

using namespace std;

inline Vec minmod(const Vec slope1, const Vec slope2){
   const Vec zero(0.0);
   Vec slope=select(abs(slope1) < abs(slope2), slope1, slope2);
   //check for extrema          
   return select(slope1 * slope2 <= 0, zero, slope);
}

inline Vec maxmod(const Vec slope1, const Vec slope2){
   const Vec zero(0.0);
   Vec slope=select(abs(slope1) > abs(slope2), slope1, slope2);
   //check for extrema          
   return select(slope1 * slope2 <= 0, zero, slope);
}

/*!
  Superbee slope limiter
*/

inline Vec slope_limiter_sb(const Vec& l,const Vec& m, const Vec& r) {
  Vec a=r-m;
  Vec b=m-l;
  const Vec slope1=minmod(a, 2*b);
  const Vec slope2=minmod(2*a, b);
  return maxmod(slope1, slope2);
}

/*!
  Minmod slope limiter
*/

inline Vec slope_limiter_minmod(const Vec& l,const Vec& m, const Vec& r) {
   //Vec sign;
   Vec a=r-m;
   Vec b=m-l; 
   return minmod(a,b);
}

/*!
  MC slope limiter
*/

inline Vec slope_limiter_mc(const Vec& l,const Vec& m, const Vec& r) {
  //Vec sign;
  Vec a=r-m;
  Vec b=m-l; 
  Vec minval=min(two*abs(a),two*abs(b));
  minval=min(minval,half*abs(a+b));
  
  //check for extrema
  Vec output = select(a*b < 0,zero,minval);
  //set sign
  return select(a + b < 0,-output,output);
}

inline Vec slope_limiter_minmod_amr(const Vec& l,const Vec& m, const Vec& r,const Vec& a,const Vec& b) {
   Vec J = r-l;
   Vec f = (m-l)/J;
   f = min(Vec(1.0),f);
   return min(f/(1+a),(Vec(1.)-f)/(1+b))*2*J;
}

inline Vec slope_limiter(const Vec& l,const Vec& m, const Vec& r) {
   return slope_limiter_sb(l,m,r);
   //return slope_limiter_minmod(l,m,r);
}

/*
 * @param a Cell size fraction dx[i-1]/dx[i] = 1/2, 1, or 2.
 * @param b Cell size fraction dx[i+1]/dx[i] = 1/2, 1, or 2.
 * @return Limited value of slope.*/
inline Vec slope_limiter_amr(const Vec& l,const Vec& m, const Vec& r,const Vec& dx_left,const Vec& dx_rght) {
   return slope_limiter_minmod_amr(l,m,r,dx_left,dx_rght);
}

/* Slope limiter with abs and sign separatelym, uses the currently active slope limiter*/
inline void slope_limiter(const Vec& l,const Vec& m, const Vec& r, Vec& slope_abs, Vec& slope_sign) {
   const Vec slope=slope_limiter(l,m,r);
   slope_abs=abs(slope);
   slope_sign=select(slope > 0, Vec(1.0), Vec(-1.0));
}

#endif
