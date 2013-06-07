/*
This file is part of Vlasiator.

Copyright 2010, 2011, 2012, 2013 Finnish Meteorological Institute












*/

#ifndef LEVEQUE_COMMON_H
#define LEVEQUE_COMMON_H

#include "definitions.h"
#include "limiters.h"

creal ZERO    = 0.0;
creal HALF    = 0.5;
creal FOURTH  = 1.0/4.0;
creal SIXTH   = 1.0/6.0;
creal ONE     = 1.0;
creal TWO     = 2.0;
creal EPSILON = 1.0e-25;

template<typename T> inline T limiter(const T& THETA_UP,const T& THETA_LO,const T& XCC) {
   return MClimiter(THETA_UP/THETA_LO);
   //return superbee(THETA_UP/THETA_LO);
   //return vanLeer(THETA_UP/THETA_LO);
}

#endif
