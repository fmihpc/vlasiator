/*
This file is part of Vlasiator.

Copyright 2010, 2011 Finnish Meteorological Institute

Vlasiator is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3
as published by the Free Software Foundation.

Vlasiator is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef CPU_COMMON_H
#define CPU_COMMON_H

#include "cmath"

#include "../definitions.h"
#include "../common.h"

// ******************* DECLARATIONS *******************

template<typename UINT> inline UINT cellIndex(const UINT& i,const UINT& j,const UINT& k) {
   return k*WID2 + j*WID + i;
}

// ****************************************************
// ***** SLOPE LIMITERS FOR 2ND ORDER FVM SOLVERS *****
// ****************************************************

template<typename T> inline T minmod(const T& a,const T& b);
template<typename T> inline T nolimiter(const T& xl1,const T& xcc,const T& xr1);
template<typename T> inline T MClimiter(const T& theta);
template<typename T> inline T MClimiter(const T& xl1,const T& xcc,const T& xr1);
template<typename T> inline T superbee(const T& theta);
template<typename T> inline T superbee(const T& xl1,const T& xcc,const T& xr1);
template<typename T> inline T vanLeer(const T& theta);
template<typename T> inline T vanLeer(const T& xl1,const T& xcc,const T& xr1);

// ***************************************************
// ***** FUNCTIONS WHICH RECONSTRUCT FACE VALUES *****
// *****  FROM GIVEN VOL. AVERAGE & DERIVATIVES  *****
// *****         2ND ORDER FVM SOLVERS           *****
// ***************************************************

template<typename T> inline T reconstruct_neg(const T& avg,const T& d1,const T& d2);
template<typename T> inline T reconstruct_pos(const T& avg,const T& d1,const T& d2);

// **********************************************
// ***** LIMITERS AND RECONSTRUCT FUNCTIONS *****
// *****         SPECIFIC TO PPM            *****
// **********************************************

template<typename T> inline T limiter_ppm(const T& xl1,const T& xcc,const T& xr1);
template<typename T> inline T reconstruct_neg_ppm(const T& avg,const T& d1,const T& d2);
template<typename T> inline T reconstruct_pos_ppm(const T& avg,const T& d1,const T& d2);

// ******************* DEFINITIONS *******************

template<typename T> inline T minmod(const T& a,const T& b) {
   if (a*b < convert<T>(0.0)) return 0.0;
   if (fabs(a) < fabs(b)) return a;
   else return b;
}

template<typename T> inline T nolimiter(const T& xl1,const T& xcc,const T& xr1) {
   return convert<T>(0.5)*(xr1-xl1);
}
   
template<typename T> inline T MClimiter(const T& theta) {
   T val1;
   val1 = std::min(convert<T>(0.5)*(convert<T>(1.0)+theta),convert<T>(2.0));
   val1 = std::min(val1,2*theta);
   return std::max(convert<T>(0.0),val1);
}

template<typename T> inline T MClimiter(const T& xl1,const T& xcc,const T& xr1) {
   const T forw = xr1-xcc;
   const T back = xcc-xl1;
   const T cntr = convert<T>(0.5)*(xr1-xl1);
   const T slope1 = minmod(convert<T>(2.0)*forw,cntr);
   const T slope2 = minmod(convert<T>(2.0)*back,cntr);
   return minmod(slope1,slope2);
}

template<typename T> inline T modbee(const T& theta) {
   const T lim = convert<T>(0.5)*(convert<T>(1.0)-std::tanh(convert<T>(0.8)*(theta-convert<T>(6.0))));
   return std::max(convert<T>(0.0),std::max(std::min(convert<T>(2.0)*theta,convert<T>(1.0)),std::min(theta,convert<T>(2.0))*lim));
   //static const T CNST = 0.8;
   //static const T THETA0 = 6.0;
   //lim = 0.5*(1.0-tanh(CNST*(theta-THETA0)));
   //rvalue = std::max(std::min(2*theta,T1),std::min(theta,T2)*lim);
   //return std::max(T0,rvalue);
}

template<typename T> inline T modbee2(const T& theta) {
   T rvalue = superbee(theta);
   const T theta_pos = std::max(convert<T>(0.0),theta);
   const T CONST = 1000.0;
   return rvalue - (1.0 - std::exp(-theta_pos/CONST));
}

template<typename T> inline T superbee(const T& theta) {
   T value1 = std::min(convert<T>(2.0),theta);
   T value2 = std::min(convert<T>(1.0),convert<T>(2.0)*theta);
   value1 = std::max(convert<T>(0.0),value1);
   return std::max(value1,value2);
}

template<typename T> inline T superbee(const T& xl1,const T& xcc,const T& xr1) {
   const T forw = xr1-xcc;
   const T back = xcc-xl1;
   
   T tmp = std::min(fabs(back),fabs(forw));
   tmp = std::min(tmp,static_cast<T>(0.5*fabs(back)));
   tmp = std::min(tmp,static_cast<T>(fabs(forw)));
   
   T sgnfwd = convert<T>(1.0);
   if (forw < convert<T>(0.0)) sgnfwd = convert<T>(-1.0);
   T sgnbck = convert<T>(1.0);
   if (back < convert<T>(0.0)) sgnbck = convert<T>(-1.0);
   
   return (sgnfwd+sgnbck)*tmp;
}

template<typename T> inline T vanLeer(const T& theta) {
   return (theta + fabs(theta))/(convert<T>(1.0)+fabs(theta));
}

template<typename T> inline T vanLeer(const T& xl1,const T& xcc,const T& xr1) {
   const T EPS = convert<T>(1.0e-20);
   return std::max((xr1-xcc)*(xcc-xl1),convert<T>(0.0))/(EPS+xr1-xl1);
}

template<typename T> inline T limiter_ppm(const T& xl1,const T& xcc,const T& xr1) {
   const T forw = xr1 - xcc;
   const T back = xcc - xl1;
   if (forw*back < 0.0) return 0.0;
   
   const T dx = nolimiter(xl1,xcc,xr1);
   T rvalue = std::min(dx,static_cast<T>(2.0*fabs(back)));
   rvalue = std::min(rvalue,static_cast<T>(2.0*fabs(forw)));
   return rvalue;
}

template<typename T> inline T reconstruct_neg(const T& avg,const T& d1,const T& d2) {
   const T INV08 = 1.0/8.0;
   const T INV24 = 1.0/24.0;
   return avg - d2*INV24 + convert<T>(0.5)*d1 + d2*INV08;
   //return avg - d2/convert<T>(24.0) + convert<T>(0.5)*d1 + d2/convert<T>(8.0);
}

template<typename T> inline T reconstruct_pos(const T& avg,const T& d1,const T& d2) {
   const T INV08 = 1.0/8.0;
   const T INV24 = 1.0/24.0;
   return avg - d2*INV24 - convert<T>(0.5)*d1 + d2*INV08;
   //return avg - d2/convert<T>(24.0) - convert<T>(0.5)*d1 + d2/convert<T>(8.0);
}

template<typename T> inline T reconstruct_neg_ppm(const T& avg,const T& d1,const T& d2) {
   return d2;
}

template<typename T> inline T reconstruct_pos_ppm(const T& avg,const T& d1,const T& d2) {
   return d1;
}

#endif
