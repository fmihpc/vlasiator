/** This file is part of VLSV file format.
 * 
 *  Copyright 2011-2013 Finnish Meteorological Institute
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef MPI_CONVERSION_H
#define MPI_CONVERSION_H
#include <mpi.h>
#include <cstdlib>
#include <iostream>
#include <stdint.h>

// Overloaded templates which return the corresponding data type
// for C++ native data types. For example, if float has been 
// typedef'd as Real, then MPI_Type<Real>() should return MPI_FLOAT.
// If you later on change the typedef to double, MPI_Type<Real>() 
// still works.
template<typename T> inline MPI_Datatype MPI_Type() {
   std::cerr << "(mpiconversion.h): NULL datatype returned, byte size of original is " << sizeof(T) << std::endl;
   exit(1);
   return 0;
}

template<> inline MPI_Datatype MPI_Type<char>() {return MPI_CHAR;}

// Signed integer types
template<> inline MPI_Datatype MPI_Type<signed char>() {return MPI_CHAR;}
template<> inline MPI_Datatype MPI_Type<short>() {return MPI_SHORT;}
template<> inline MPI_Datatype MPI_Type<int>() {return MPI_INT;}
template<> inline MPI_Datatype MPI_Type<long int>() {return MPI_LONG;}
template<> inline MPI_Datatype MPI_Type<long long int>() {return MPI_LONG_LONG;}

// Unsigned integer types
template<> inline MPI_Datatype MPI_Type<unsigned char>() {return MPI_UNSIGNED_CHAR;}
template<> inline MPI_Datatype MPI_Type<unsigned short>() {return MPI_UNSIGNED_SHORT;}
template<> inline MPI_Datatype MPI_Type<unsigned int>() {return MPI_UNSIGNED;}
template<> inline MPI_Datatype MPI_Type<unsigned long int>() {return MPI_UNSIGNED_LONG;}
template<> inline MPI_Datatype MPI_Type<unsigned long long int>() {return MPI_UNSIGNED_LONG_LONG;}

// Floating point types
template<> inline MPI_Datatype MPI_Type<float>() {return MPI_FLOAT;}
template<> inline MPI_Datatype MPI_Type<double>() {return MPI_DOUBLE;}
template<> inline MPI_Datatype MPI_Type<long double>() {return MPI_LONG_DOUBLE;}

#endif

