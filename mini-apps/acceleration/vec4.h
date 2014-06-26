#ifndef VECTORCLASS_INTERFACE_H
#define VECTORCLASS_INTERFACE_H


/*!

\file vec4.h 
\brief An interface to a type with four floating point values

Depending on if DP is defined or not, the type is double or float. If
USE_AGNER_VECTORCLASS is defined a fast library supporting x86
computers is used, othwerwise a simple but portable implementation in
vectorclass_fallback.h is used.

*/


#ifdef USE_AGNER_VECTORCLASS
//user Agner's SSE/AVX optimized datatypes
#include "vectorclass.h"
#ifdef DP
typedef Vec4d Vec4;
#else
typedef Vec4f Vec4;
#endif

#else
//user portable vectorclass
#include "vectorclass_fallback.h"
#ifdef DP
typedef Vec4Real<double> Vec4;
#else
typedef Vec4Real<float> Vec4;
#endif

#endif



const Vec4 one(1.0);
const Vec4 minus_one(-1.0);
const Vec4 two(2.0);
const Vec4 half(0.5);
const Vec4 zero(0.0);

const Vec4 one_sixth(1.0/6.0);
const Vec4 one_twelfth(1.0/12.0);
const Vec4 seven_twelfth(7.0/12.0);
const Vec4 one_third(1.0/3.0);

#endif


