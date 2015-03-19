/*
This file is part of Vlasiator.

Copyright 2010, 2011, 2012, 2013 Finnish Meteorological Institute
*/

#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#include <stdint.h>
#include <limits>

//set floating point precision for storing the distribution function here. Default is single precision, use -DDPF to set double precision
#ifdef DPF
typedef double Realf;
#else
typedef float Realf;
#endif

//set general floating point precision here. Default is single precision, use -DDP to set double precision
#ifdef DP
typedef double Real;
typedef const double creal;
#else
typedef float Real;
typedef const float creal;
#endif

typedef const int cint;
typedef unsigned char uchar;
typedef const unsigned char cuchar;

typedef uint32_t uint;
typedef const uint32_t cuint;

typedef cuint csize;

typedef uint64_t CellID;

template<typename T> T convert(const T& number) {return number;}

namespace vmesh {
   #ifndef AMR
   typedef uint32_t GlobalID;              /**< Datatype used for velocity block global IDs.*/
   typedef uint32_t LocalID;               /**< Datatype used for velocity block local IDs.*/
   #else
   typedef uint32_t GlobalID;
   typedef uint32_t LocalID;
   #endif

   /** Global ID of a non-existing or otherwise erroneous velocity block.*/
   static const GlobalID INVALID_GLOBALID = std::numeric_limits<GlobalID>::max();

   /** Local ID of a non-existing or otherwise erroneous velocity block.*/
   static const LocalID INVALID_LOCALID  = std::numeric_limits<LocalID>::max();

   /** Block index of a non-existing or erroneous velocity block.*/
   static const LocalID INVALID_VEL_BLOCK_INDEX = INVALID_LOCALID;
}

/** Definition of a function that takes in a velocity block with neighbor data, 
 * and returns a number that is used to decide whether or not the block should 
 * be refined or coarsened.*/
typedef Realf (*AmrVelRefinement)(const Realf* velBlock);

#endif
