
#ifndef LEVEQUE_LIMITER_VECTORIZED_H
#define LEVEQUE_LIMITER_VECTORIZED_H
#include  "vec4.h"


inline Vec4 superbee_vec4(const Vec4& theta) {
   const Vec4 two(2);
   const Vec4 one(1.0);
   const Vec4 zero(0.0);
   Vec4 value1 = min(two,theta);
   Vec4 value2 = min(one,two*theta);
   value1 = max(zero,value1);
   return max(value1,value2);
}

inline Vec4 MClimiter(const Vec4& theta) {
   Vec4 val1;
   const Vec4 two(2.0);
   const Vec4 one(1.0);
   const Vec4 half(0.5);
   const Vec4 zero(0.0);
   val1 = min(half*(one+theta),two);
   val1 = min(val1,two*theta);
   return max(zero,val1);
}


inline Vec4 limiter_vec4(const Vec4& THETA_UP,const Vec4& THETA_LO,const Vec4& XCC) {
   return MClimiter(THETA_UP/THETA_LO);
   //return superbee_vec4(THETA_UP/THETA_LO);
   //return vanLeer(THETA_UP/THETA_LO);
}

#endif
