#include "definitions.h"

// ******************* DECLARATIONS *******************

template<typename T> T minmod(const T& a,const T& b);
template<typename T> T MClimiter(const T& xl1,const T& xcc,const T& xr1);
template<typename T> T superbee(const T& xl1,const T& xcc,const T& xr1);
template<typename T> T vanLeer(const T& xl1,const T& xcc,const T& xr1);

template<typename T> T reconstruct_neg(const T& avg,const T& d1,const T& d2);
template<typename T> T reconstruct_pos(const T& avg,const T& d1,const T& d2);

// ******************* DEFINITIONS *******************

template<typename T> T minmod(const T& a,const T& b) {
   if (a*b < convert<T>(0.0)) return 0.0;
   if (fabs(a) < fabs(b)) return a;
   else return b;
}

template<typename T> T MClimiter(const T& xl1,const T& xcc,const T& xr1) {
   const T forw = xr1-xcc;
   const T back = xcc-xl1;
   const T cntr = convert<T>(0.5)*(xr1-xl1);
   const T slope1 = minmod(convert<T>(2.0)*forw,cntr);
   const T slope2 = minmod(convert<T>(2.0)*back,cntr);
   return minmod(slope1,slope2);
}

template<typename T> T superbee(const T& xl1,const T& xcc,const T& xr1) {
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

template<typename T> T vanLeer(const T& xl1,const T& xcc,const T& xr1) {
   const T EPS = convert<T>(1.0e-20);
   return std::max((xr1-xcc)*(xcc-xl1),convert<T>(0.0))/(EPS+xr1-xl1);
}

template<typename T> T reconstruct_neg(const T& avg,const T& d1,const T& d2) {
   const T INV08 = 1.0/8.0;
   const T INV24 = 1.0/24.0;
   return avg - d2*INV24 + convert<T>(0.5)*d1 + d2*INV08;
   //return avg - d2/convert<T>(24.0) + convert<T>(0.5)*d1 + d2/convert<T>(8.0);
}

template<typename T> T reconstruct_pos(const T& avg,const T& d1,const T& d2) {
   const T INV08 = 1.0/8.0;
   const T INV24 = 1.0/24.0;
   return avg - d2*INV24 - convert<T>(0.5)*d1 + d2*INV08;
   //return avg - d2/convert<T>(24.0) - convert<T>(0.5)*d1 + d2/convert<T>(8.0);
}





