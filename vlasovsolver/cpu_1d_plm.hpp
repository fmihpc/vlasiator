/*
This file is part of Vlasiator.
Copyright 2013, 2014 Finnish Meteorological Institute

*/

#ifndef CPU_1D_PLM_H
#define CPU_1D_PLM_H

#include <iostream>
#include "vec4.h"
#include "algorithm"
#include "cmath"
#include "cpu_slope_limiters.hpp"

using namespace std;

/*!
 Compute PLM coefficients
 f(v) = a[0] + a[1]/2.0*t 
t=(v-v_{i-0.5})/dv where v_{i-0.5} is the left face of a cell
The factor 2.0 is in the polynom to ease integration, then integral is a[0]*t + a[1]*t**2
*/

inline void compute_plm_coeff(const Vec4 * const values, uint k, Vec4 a[2]){
   const Vec4 d_cv=slope_limiter(values[k - 1], values[k], values[k + 1]);
   a[0] = values[k] - d_cv * 0.5;
   a[1] = d_cv * 0.5;
}

inline void compute_plm_coeff_amr(const Vec4 * const values, uint k, Vec4 a[2],const Vec4& dx_left,const Vec4& dx_rght) {
   const Vec4 d_cv = slope_limiter_amr(values[k - 1],values[k],values[k + 1],dx_left,dx_rght);
   a[0] = values[k] - d_cv * 0.5;
   a[1] = d_cv * 0.5;
}

#endif
