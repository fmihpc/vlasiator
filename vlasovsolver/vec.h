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
#ifndef VECTORCLASS_INTERFACE_H
#define VECTORCLASS_INTERFACE_H


/*!

\file vec.h 
\brief An interface to a type with  floating point values

By setting suitable compile-time defines oen can set the length,
accuracy and implementation of the vector. The vector length has to be
a multiple of WID, which is 4 in Vlasiator. It also cannot be larger
than WID*WID or 16. Thus 4, 8 or 16 are now supported. Currently
implemented vector backends are:

VEC4D_AGNER
 - Double precision
 - Vector length of 4
 - Use Agner's vectorclass with AVX intrinisics

VEC4F_AGNER
 - Single precision
 - Vector length of 4
 - Use Agner's vectorclass with AVX/SSE2 intrinisics

VEC8F_AGNER
 - Single precision     
 - Vector length of 8
 - Use Agner's vectorclass with AVX intrinisics

 
*/




#ifdef VEC4D_AGNER
//user Agner's AVX2 optimized datatypes, double precision accuracy
#include "vectorclass.h"
typedef Vec4d Vec;
#if VECTORCLASS_H >= 20000
typedef Vec4q Veci;
#else
typedef Vec4i Veci;
#endif
typedef Vec4db Vecb;
typedef double Realv;
#define to_realv(v) to_double(v)
#define VECL 4
#define VPREC 8
#define VEC_PER_PLANE 4 //vectors per plane in block
#define VEC_PER_BLOCK 16
#endif

#ifdef VEC8D_AGNER
//user Agner's AVX512 optimized datatypes, double precision accuracy
#include "vectorclass.h"
typedef Vec8d Vec;
typedef Vec8i Veci;
typedef Vec8db Vecb;
typedef double Realv;
#define to_realv(v) to_double(v)
#define VECL 8
#define VPREC 8
#define VEC_PER_PLANE 2 //vectors per plane in block
#define VEC_PER_BLOCK 8
#endif

#ifdef VEC4F_AGNER
//user Agner's SSEx optimized datatypes, double precision accuracy
#include "vectorclass.h"
typedef Vec4f Vec;
typedef Vec4i Veci;
typedef Vec4fb Vecb;
typedef float Realv;
#define to_realv(v) to_float(v)
#define VECL 4
#define VPREC 4
#define VEC_PER_PLANE 4 //vectors per plane in block
#define VEC_PER_BLOCK 16
#endif

#ifdef VEC8F_AGNER
//user Agner's AVX2 optimized datatypes, single precision accuracy
#include "vectorclass.h"
typedef Vec8f Vec;
typedef Vec8i Veci;
typedef Vec8fb Vecb;
typedef float Realv;
#define to_realv(v) to_float(v)
#define VECL 8
#define VPREC 4
#define VEC_PER_PLANE 2 //vectors per plane in block
#define VEC_PER_BLOCK 8
#endif


#ifdef VEC16F_AGNER
//user Agner's AVX512 optimized datatypes, single precision accuracy
#include "vectorclass.h"
typedef Vec16f Vec;
typedef Vec16i Veci;
typedef Vec16fb Vecb;
typedef float Realv;
#define to_realv(v) to_float(v)
#define VECL 16
#define VPREC 4
#define VEC_PER_PLANE 1 //vectors per plane in block
#define VEC_PER_BLOCK 4
#endif



#ifdef VEC4D_FALLBACK
//user portable vectorclass
#include "vectorclass_fallback.h"
typedef Vec4Simple<double> Vec;
typedef Vec4Simple<bool> Vecb;
typedef Vec4Simple<int> Veci;
typedef double Realv;
#define to_realv(v) to_double(v)
#define VECL 4
#define VPREC 8
#define VEC_PER_PLANE 4 //vectors per plane in block
#define VEC_PER_BLOCK 16
#endif

#ifdef VEC4F_FALLBACK
//user portable vectorclass
#include "vectorclass_fallback.h"
typedef Vec4Simple<float> Vec;
typedef Vec4Simple<bool> Vecb;
typedef Vec4Simple<int> Veci;
typedef float Realv;
#define to_realv(v) to_float(v)
#define VECL 4
#define VPREC 4
#define VEC_PER_PLANE 4 //vectors per plane in block
#define VEC_PER_BLOCK 16
#endif

#ifdef VEC8D_FALLBACK
//user portable vecto rclass
#include "vectorclass_fallback.h"
typedef Vec8Simple<double> Vec;
typedef Vec8Simple<bool> Vecb;
typedef Vec8Simple<int> Veci;
typedef double Realv;
#define to_realv(v) to_double(v)
#define VECL 8
#define VPREC 8
#define VEC_PER_PLANE 2 //vectors per plane in block
#define VEC_PER_BLOCK 8
#endif


#ifdef VEC8F_FALLBACK
//user portable vectorclass
#include "vectorclass_fallback.h"
typedef Vec8Simple<float> Vec;
typedef Vec8Simple<bool> Vecb;
typedef Vec8Simple<int> Veci;
typedef float Realv;
#define to_realv(v) to_float(v)
#define VECL 8
#define VPREC 4
#define VEC_PER_PLANE 2 //vectors per plane in block
#define VEC_PER_BLOCK 8
#endif


const Vec one(1.0);
const Vec minus_one(-1.0);
const Vec two(2.0);
const Vec half(0.5);
const Vec zero(0.0);
const Vec one_sixth(1.0/6.0);
const Vec one_twelfth(1.0/12.0);
const Vec seven_twelfth(7.0/12.0);
const Vec one_third(1.0/3.0);





#endif
