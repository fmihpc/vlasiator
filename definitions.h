/*
This file is part of Vlasiator.

Copyright 2010, 2011, 2012, 2013 Finnish Meteorological Institute
*/

#ifndef DEFINITIONS_H
#define DEFINITIONS_H

# include <stdint.h>

//set floating point precision here. Default is single precision, use -DDP to set double precision
#ifdef DPF
typedef double Realf;
#else
typedef float Realf;
#endif

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



#endif
