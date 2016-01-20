/*
This file is part of Vlasiator.

Copyright 2010-2015 Finnish Meteorological Institute

*/

/*! \file fs_limiters.h
 * 
 * \brief Definitions of the limiter functions used in the field solver.
 * 
 * Definitions of the limiter functions used in the field solver.
 * The three-point limiter functions are minmod, MC, superbee and vanLeer. Their helper functions are the two-point minmod and the sign function.
 */

#ifndef FS_LIMITERS_H
#define FS_LIMITERS_H

#include <cstdlib>
#include <algorithm>
#include <limits>

template<typename T> inline int sign(const T& value) {
   const T ZERO = 0.0;
   if (value < ZERO) return -1;
   return 1;
}

template<typename T> inline T minmod(const T& a,const T& b) {
   const T ZERO = 0.0;
   
   if (a*b < ZERO) return ZERO;
   else if (fabs(a) < fabs(b)) return a;
   else return b;
}

template<typename T> inline T minmod(const T& left,const T& cent,const T& rght) {
   const T HALF = 0.5;
   const T val1 = cent-left;
   const T val2 = rght-cent;
   return HALF*(sign(val1)+sign(val2))*std::min(fabs(val1),fabs(val2));
}

template<typename T> inline T MClimiter(const T& left,const T& cent,const T& right) {
   const T HALF = 0.5;
   const T TWO  = 2.0;
   
   const T forw = right-cent;
   const T back = cent-left;
   const T cntr = HALF*(right-left);
   const T slope1 = minmod(TWO*forw,cntr);
   const T slope2 = minmod(TWO*back,cntr);
   return minmod(slope1,slope2);
}

template<typename T> inline T superbee(const T& left,const T& cent,const T& right) {
   const T HALF = 0.5; 
   
   const T back = cent-left;
   const T forw = right-cent;
   T tmp = std::min(fabs(back),fabs(forw));
   tmp = std::min(tmp,HALF*std::max(fabs(back),fabs(forw)));
   return (sign(back)+sign(forw))*tmp;
}

template<typename T> inline T vanLeer(const T& left,const T& cent,const T& right) {
   const T EPSILON = std::numeric_limits<T>::min();
   const T ZERO    = 0.0;
   const T TWO     = 2.0;

   const T numerator = std::max((right-cent)*(cent-left),ZERO);
   const T denumerator = (right-left)+EPSILON;
   return TWO * numerator / denumerator;
   
   //return TWO*std::max((right-cent)*(cent-left),ZERO)/(right-left+EPSILON);

   //const T denumerator = (right-left) + EPSILON;
   //const T rvalue = TWO*std::max((right-cent)*(cent-left),ZERO) / denumerator;
   /*
   if (rvalue > 1e10) {
      std::cerr << "ERROR vanLeer rvalue too large " << std::max((right-cent)*(cent-left),ZERO) << '\t' << denumerator << std::endl;
      std::cerr << "\t input values are " << left << '\t' << cent << '\t' << right << std::endl;
      std::cerr << "\t numerator   " << TWO*std::max((right-cent)*(cent-left),ZERO) << std::endl;
      std::cerr << "\t denumerator " << denumerator << std::endl;
   }
   //return TWO*std::max((right-cent)*(cent-left),ZERO) / denumerator;
    */
}

template<typename T> inline T limiter(const T& l_rho_v,const T& l_rho,
                               const T& c_rho_v,const T& c_rho,
                               const T& r_rho_v,const T& r_rho) {
   T left = l_rho_v / l_rho;
   if (l_rho <= 0) left = 0.0;
   
   T cent = c_rho_v / c_rho;
   if (c_rho <= 0) cent = 0.0;
   
   T rght = r_rho_v / r_rho;
   if (r_rho <= 0) rght = 0.0;

   return vanLeer(left,cent,rght);
}

template<typename T> inline T limiter(const T& left,const T& cent,const T& rght) {
   return vanLeer(left,cent,rght);
}

/*! Select the limiter to be used in the field solver. */
//Real limiter(creal& left,creal& cent,creal& rght);

#endif
