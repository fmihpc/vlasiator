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


#define MAX_BLOCKS_PER_DIM 1000
#ifdef ACC_SEMILAG_PLM
#define RECONSTRUCTION_ORDER 1
#endif
#ifdef ACC_SEMILAG_PPM
#define RECONSTRUCTION_ORDER 2
#endif
#ifdef ACC_SEMILAG_PQM
#define RECONSTRUCTION_ORDER 4
#endif

// Natural constants
namespace physicalconstants {
   const Real MU_0 = 1.25663706e-6;  /*!< Permeability of vacuo, unit: (kg m) / (s^2 A^2).*/
   const Real K_B = 1.3806503e-23;   /*!< Boltzmann's constant, unit: (kg m^2) / (s^2 K).*/
   const Real CHARGE = 1.60217653e-19; /*!< Elementary charge, unit: C. */
   const Real MASS_PROTON = 1.67262158e-27; /*!< Proton rest mass.*/
   const Real R_E = 6.3712e6; /*!< radius of the Earth. */
}


#endif


