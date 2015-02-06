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
//user Agner's SSE/AVX optimized datatypes, double precision accuracy
#include "vectorclass.h"
typedef Vec4d Vec;
typedef Vec4i Veci;
typedef Vec4db Vecb;
typedef double Realv;
#define to_realv(v) to_double(v)
#define VECL 4
#define VEC_PER_PLANE 4 //vectors per plane in block
#define VEC_PER_BLOCK 16
#endif

#ifdef VEC4F_AGNER
//user Agner's SSE/AVX optimized datatypes, double precision accuracy
#include "vectorclass.h"
typedef Vec4f Vec;
typedef Vec4i Veci;
typedef Vec4fb Vecb;
typedef float Realv;
#define to_realv(v) to_float(v)
#define VECL 4
#define VEC_PER_PLANE 4 //vectors per plane in block
#define VEC_PER_BLOCK 16
#endif

#ifdef VEC8F_AGNER
#include "vectorclass.h"
typedef Vec8f Vec;
typedef Vec8i Veci;
typedef Vec8fb Vecb;
typedef float Realv;
#define to_realv(v) to_float(v)
#define VECL 8
#define VEC_PER_PLANE 2 //vectors per plane in block
#define VEC_PER_BLOCK 8
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
#define VEC_PER_PLANE 4 //vectors per plane in block
#define VEC_PER_BLOCK 16
#endif

#ifdef VEC4F_FALLBACK
//user portable vectorclass
#include "vectorclass_fallback.h"
typedef Vec4Simple<float> Vec4;
typedef Vec4Simple<bool> Vec4b;
typedef Vec4Simple<int> Vec4i;
typedef float Realv;
#define to_realv(v) to_float(v)
#define VECL 4
#define VEC_PER_PLANE 4 //vectors per plane in block
#define VEC_PER_BLOCK 16
#endif


const Vec one(1);
const Vec minus_one(-1);
const Vec two(2);
const Vec zero(0);






#endif
