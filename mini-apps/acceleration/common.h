/*
This file is part of Vlasiator.

Copyright 2010, 2011, 2012, 2013 Finnish Meteorological Institute

*/

#ifndef COMMON_H
#define COMMON_H
# include <stdint.h>

typedef uint32_t uint;
typedef const uint32_t cuint;


const uint WID = 4;         /*!< Number of cells per coordinate in a velocity block. Only a value of 4 supported by vectorized Leveque solver */
const uint WID2 = WID*WID;  /*!< Number of cells per 2D slab in a velocity block. */
const uint WID3 = WID2*WID; /*!< Number of cells in a velocity block. */

//set general floating point precision here. Default is single precision, use -DDP to set double precision
#ifdef DP
typedef double Real;
typedef const double creal;
#else
typedef float Real;
typedef const float creal;
#endif


#define MAX_BLOCKS_PER_DIM 100
#ifdef ACC_SEMILAG_PLM
#define RECONSTRUCTION_ORDER 1
#endif
#ifdef ACC_SEMILAG_PPM
#define RECONSTRUCTION_ORDER 2
#endif



#endif


