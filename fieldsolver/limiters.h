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

/*! \file limiters.h
 * 
 * \brief Definitions of the limiter functions used in the field solver.
 * 
 * Definitions of the limiter functions used in the field solver.
 * The three-point limiter functions are minmod, MC, superbee and vanLeer. Their helper functions are the two-point minmod and the sign function.
 */

#ifndef LIMITERS_H
#define LIMITERS_H

#include <cmath>

template<typename T> int sign(const T& value) {
   const T ZERO = 0.0;
   if (value < ZERO) return -1;
   return 1;
}

template<typename T> T minmod(const T& a,const T& b) {
   const T ZERO = 0.0;
   
   if (a*b < ZERO) return ZERO;
   else if (fabs(a) < fabs(b)) return a;
   else return b;
}

template<typename T> T minmod(const T& left,const T& cent,const T& rght) {
   const T HALF = 0.5;
   const T val1 = cent-left;
   const T val2 = rght-cent;
   return HALF*(sign(val1)+sign(val2))*std::min(fabs(val1),fabs(val2));
}

template<typename T> T MClimiter(const T& left,const T& cent,const T& right) {
   const T HALF = 0.5;
   const T TWO  = 2.0;
   
   const T forw = right-cent;
   const T back = cent-left;
   const T cntr = HALF*(right-left);
   const T slope1 = minmod(TWO*forw,cntr);
   const T slope2 = minmod(TWO*back,cntr);
   return minmod(slope1,slope2);
}

template<typename T> T superbee(const T& left,const T& cent,const T& right) {
   const T HALF = 0.5;
   
   const T back = cent-left;
   const T forw = right-cent;
   T tmp = std::min(fabs(back),fabs(forw));
   tmp = std::min(tmp,HALF*std::max(fabs(back),fabs(forw)));
   return (sign(back)+sign(forw))*tmp;
}

template<typename T> T vanLeer(const T& left,const T& cent,const T& right) {
   const T EPSILON = std::numeric_limits<T>::min();
   const T ZERO    = 0.0;
   const T TWO = convert<T>(2.0);
   return TWO*std::max((right-cent)*(cent-left),ZERO)/(right-left+EPSILON);
}


#endif
