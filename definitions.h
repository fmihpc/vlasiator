/*
This file is part of Vlasiator.

Copyright 2010, 2011, 2012, 2013 Finnish Meteorological Institute
*/

#ifndef DEFINITIONS_H
#define DEFINITIONS_H

# include <stdint.h>

#ifndef VLASIATOR_NDEBUG
#include <cstdlib>
#include <stdio.h>
#include <mpi.h>
#endif

//set floating point precision here. Default is single precision, use -DDP to set double precision
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

// Assert function for debugging Vlasiator
void Assert( bool condition, const char * error, const char * file, int line ) {
        #ifndef VLASIATOR_NDEBUG
        if(!condition) {
           #pragma omp critical
           {
              int rank;
              MPI_Comm_rank(MPI_COMM_WORLD, &rank);
              fprintf( stderr, "ASSERT ERROR: %s; Rank: %i, File: %s, line: %i", error, file, line );
              exit(1);
           }
        }
        #endif
        return;
}

#ifndef VLASIATOR_NDEBUG
#define Assert(a) Assert( a, #a, __FILE__, __LINE__ )    //The a is the truth-value input, b is a char* input, aka const char* in our case, __FILE__ is recognized as the filename and __LINE__ the line number 
#else
#define Assert(a)
#endif

#endif
