/*
This file is part of Vlasiator.

Copyright 2010, 2011, 2012, 2013 Finnish Meteorological Institute












*/

#ifndef CPU_COMMON_H
#define CPU_COMMON_H

#include "cmath"

#include "../definitions.h"
#include "../common.h"

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
   const T TWO = convert<T>(2.0);
   return TWO*std::max((xr1-xcc)*(xcc-xl1),convert<T>(0.0))/(EPS+xr1-xl1);
}



#endif
