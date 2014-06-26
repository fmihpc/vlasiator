/*
This file is part of Vlasiator.
Copyright 2013, 2014 Finnish Meteorological Institute

*/

#ifndef CPU_SLOPE_LIMITERS_H
#define CPU_SLOPE_LIMITERS_H

#include "vec4.h"

using namespace std;

/*!
  MC slope limiter
*/

inline Vec4 slope_limiter(const Vec4& l,const Vec4& m, const Vec4& r) {
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

/*!
  MC limiter. Give absolute value of slope and its sign separately
*/
void slope_limiter(const Vec4& l,const Vec4& m, const Vec4& r, Vec4& slope_abs, Vec4& slope_sign) {
  const Vec4 two(2.0);
  const Vec4 half(0.5);
  const Vec4 zero(0.0);
  Vec4 sign;
  Vec4 a=r-m;
  Vec4 b=m-l; 
  Vec4 minval=min(two*abs(a),two*abs(b));
  minval=min(minval,half*abs(a+b));
  
  //check for extrema, set absolute value
  slope_abs = select(a*b < 0, zero, minval);
  slope_sign = select(a + b < 0, minus_one, one);

}

#endif
