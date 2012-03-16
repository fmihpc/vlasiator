/*
This file is part of Vlasiator.

Copyright 2010, 2011, 2012 Finnish Meteorological Institute

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

#ifndef LEVEQUE_COMMON_H
#define LEVEQUE_COMMON_H

#include "definitions.h"
#include "cpu_common.h"

creal ZERO    = 0.0;
creal HALF    = 0.5;
creal FOURTH  = 1.0/4.0;
creal SIXTH   = 1.0/6.0;
creal ONE     = 1.0;
creal TWO     = 2.0;
creal EPSILON = 1.0e-25;

template<typename T> inline T limiter(const T& THETA_UP,const T& THETA_LO,const T& XCC) {
   //return MClimiter(THETA_UP/THETA_LO);
   return superbee(THETA_UP/THETA_LO);
   //return vanLeer(THETA_UP/THETA_LO);
}

#endif
