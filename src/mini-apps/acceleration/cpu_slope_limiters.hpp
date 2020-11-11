/*
This file is part of Vlasiator.
Copyright 2013, 2014 Finnish Meteorological Institute

*/

#ifndef CPU_SLOPE_LIMITERS_H
#define CPU_SLOPE_LIMITERS_H

using namespace std;

/*!
  MC slope limiter
*/

inline Real slope_limiter(const Real& l,const Real& m, const Real& r) {
  Real sign;
  Real a=r-m;
  Real b=m-l; 
  Real minval=min(2.0 * abs(a), 2.0 * abs(b));
  minval=min(minval, 0.5 * abs(a+b));
  
  //check for extrema
  Real output = a*b < 0 ? 0.0 : minval;
  //set sign
  return a + b < 0 ? -output : output;
}

/*!
  MC limiter. Give absolute value of slope and its sign separately
*/
void slope_limiter(const Real& l,const Real& m, const Real& r, Real& slope_abs, Real& slope_sign) {
  Real sign;
  Real a=r-m;
  Real b=m-l; 
  Real minval=min(2.0 * abs(a),2.0 * abs(b));
  minval=min(minval, 0.5 * abs(a+b));
  
  //check for extrema, set absolute value
  slope_abs = a*b < 0 ? 0.0: minval;
  slope_sign = a + b < 0 ? -1.0 : 1.0;

}

#endif
