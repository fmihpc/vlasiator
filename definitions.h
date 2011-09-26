/*
This file is part of Vlasiator.

Copyright 2010, 2011 Finnish Meteorological Institute

Vlasiator is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3
as published by the Free Software Foundation.

Vlasiator is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef DEFINITIONS_H
#define DEFINITIONS_H
#include "mpi.h"

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
typedef unsigned int uint;
typedef const unsigned int cuint;
typedef long unsigned int luint;
typedef const long unsigned int cluint;
typedef long long unsigned int lluint;
typedef const long long unsigned int clluint;

typedef cuint csize;

template<typename T> T convert(const T& number) {return number;}

#endif
