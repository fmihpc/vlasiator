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

#include "../common.h"
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


#ifdef VEC_FALLBACK_GENERIC
// use portable vectorclass with specified vector length
// if not specified in Makefile, equivalent to VEC8F
#ifndef VECL
#define VECL (8)
#endif
#ifndef VEC_PER_PLANE
const int VEC_PER_PLANE = (WID*WID/VECL);
#endif
#ifndef VEC_PER_BLOCK
const int VEC_PER_BLOCK = (WID*VEC_PER_PLANE);
#endif
#ifndef VPREC
#define VPREC (8)
#endif

#include "vectorclass_fallback.h"

typedef VecSimple<bool> Vecb;
typedef VecSimple<int> Veci;

#ifdef DPF
typedef VecSimple<double> Vec;
#define to_realf(v) to_double(v)
#else
typedef VecSimple<float> Vec;
#define to_realf(v) to_float(v)
#endif

#endif

// SVE Definitions
#if (defined(VEC4F_SVE) || defined(VEC4D_SVE))
#define VECL 4
#define VEC_SVE
#endif

#if (defined(VEC8F_SVE) || defined(VEC8D_SVE))
#define VECL 8
#define VEC_SVE
#endif

#if (defined(VEC16F_SVE) || defined(VEC16D_SVE))
#define VECL 16
#define VEC_SVE
#endif

#if (defined(VEC32F_SVE) || defined(VEC32D_SVE))
#define VECL 32
#define VEC_SVE
#endif

#if (defined(VEC64F_SVE) || defined(VEC64D_SVE))
#define VECL 64
#define VEC_SVE
#endif

#if defined(VEC_SVE)

#include "vectorclass_sve.hpp"

#define VPREC 8
#ifndef VEC_PER_PLANE
const int VEC_PER_PLANE = (WID*WID/VECL);
#endif
#ifndef VEC_PER_BLOCK
const int VEC_PER_BLOCK = (WID*VEC_PER_PLANE);
#endif

#if defined(DPF)
using Vec  = VecSVE<double, VECL>;
using Veci = VecSVE<int64_t, VECL>;
using Vecb = VecSVEMask<double, VECL>;
#define to_realf(v) to_double(v)
#else
using Vec  = VecSVE<float, VECL>;
using Veci = VecSVE<int32_t, VECL>;
using Vecb = VecSVEMask<float, VECL>;
#define to_realf(v) to_float(v)
#endif

#endif // END VEC_SVE

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
#define to_realf(v) to_double(v)
#define VECL 4
#define VPREC 8
#ifndef VEC_PER_PLANE
const int VEC_PER_PLANE = (WID*WID/VECL);
#endif
#ifndef VEC_PER_BLOCK
const int VEC_PER_BLOCK = (WID*VEC_PER_PLANE);
#endif
#endif

#ifdef VEC8D_AGNER
//user Agner's AVX512 optimized datatypes, double precision accuracy
#include "vectorclass.h"
typedef Vec8d Vec;
typedef Vec8i Veci;
typedef Vec8db Vecb;
#define to_realf(v) to_double(v)
#define VECL 8
#define VPREC 8
#ifndef VEC_PER_PLANE
const int VEC_PER_PLANE = (WID*WID/VECL);
#endif
#ifndef VEC_PER_BLOCK
const int VEC_PER_BLOCK = (WID*VEC_PER_PLANE);
#endif
#endif

#ifdef VEC4F_AGNER
//user Agner's SSEx optimized datatypes, double precision accuracy
#include "vectorclass.h"
typedef Vec4f Vec;
typedef Vec4i Veci;
typedef Vec4fb Vecb;
#define to_realf(v) to_float(v)
#define VECL 4
#define VPREC 4
#ifndef VEC_PER_PLANE
const int VEC_PER_PLANE = (WID*WID/VECL);
#endif
#ifndef VEC_PER_BLOCK
const int VEC_PER_BLOCK = (WID*VEC_PER_PLANE);
#endif
#endif

#ifdef VEC8F_AGNER
//user Agner's AVX2 optimized datatypes, single precision accuracy
#include "vectorclass.h"
typedef Vec8f Vec;
typedef Vec8i Veci;
typedef Vec8fb Vecb;
#define to_realf(v) to_float(v)
#define VECL 8
#define VPREC 4
#ifndef VEC_PER_PLANE
const int VEC_PER_PLANE = (WID*WID/VECL);
#endif
#ifndef VEC_PER_BLOCK
const int VEC_PER_BLOCK = (WID*VEC_PER_PLANE);
#endif
#endif


#ifdef VEC16F_AGNER
//user Agner's AVX512 optimized datatypes, single precision accuracy
#include "vectorclass.h"
typedef Vec16f Vec;
typedef Vec16i Veci;
typedef Vec16fb Vecb;
#define to_realf(v) to_float(v)
#define VECL 16
#define VPREC 4
#ifndef VEC_PER_PLANE
const int VEC_PER_PLANE = (WID*WID/VECL);
#endif
#ifndef VEC_PER_BLOCK
const int VEC_PER_BLOCK = (WID*VEC_PER_PLANE);
#endif
#endif

#endif
