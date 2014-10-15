/*
This file is part of Vlasiator.
Copyright 2013, 2014 Finnish Meteorological Institute

*/

#ifndef CPU_SLOPE_LIMITERS_H
#define CPU_SLOPE_LIMITERS_H

#include "vec4.h"

using namespace std;

inline Vec4 minmod(const Vec4 slope1, const Vec4 slope2){
   const Vec4 zero(0.0);
   Vec4 slope=select(abs(slope1) < abs(slope2), slope1, slope2);
   //check for extrema          
   return select(slope1 * slope2 <= 0, zero, slope);
}

inline Vec4 maxmod(const Vec4 slope1, const Vec4 slope2){
   const Vec4 zero(0.0);
   Vec4 slope=select(abs(slope1) > abs(slope2), slope1, slope2);
   //check for extrema          
   return select(slope1 * slope2 <= 0, zero, slope);
}

/*!
  Superbee slope limiter
*/

inline Vec4 slope_limiter_sb(const Vec4& l,const Vec4& m, const Vec4& r) {
  Vec4 a=r-m;
  Vec4 b=m-l;
  const Vec4 slope1=minmod(a, 2*b);
  const Vec4 slope2=minmod(2*a, b);
  return maxmod(slope1, slope2);
}

/*!
  Minmod slope limiter
*/

inline Vec4 slope_limiter_minmod(const Vec4& l,const Vec4& m, const Vec4& r) {
   Vec4 sign;
   Vec4 a=r-m;
   Vec4 b=m-l; 
   return minmod(a,b);
}

/*!
  MC slope limiter
*/

inline Vec4 slope_limiter_mc(const Vec4& l,const Vec4& m, const Vec4& r) {
  Vec4 sign;
  Vec4 a=r-m;
  Vec4 b=m-l; 
  Vec4 minval=min(two*abs(a),two*abs(b));
  minval=min(minval,half*abs(a+b));
  
  //check for extrema
  Vec4 output = select(a*b < 0,zero,minval);
  //set sign
  return select(a + b < 0,-output,output);
}



inline Vec4 slope_limiter(const Vec4& l,const Vec4& m, const Vec4& r) {
   return slope_limiter_sb(l,m,r);
}

/* Slope limiter with abs and sign separatelym, uses the currently active slope limiter*/
inline void slope_limiter(const Vec4& l,const Vec4& m, const Vec4& r, Vec4& slope_abs, Vec4& slope_sign) {
   const Vec4 slope=slope_limiter(l,m,r);
   slope_abs=abs(slope);
   slope_sign=select(slope > 0, 1.0, -1.0);
}



#endif
