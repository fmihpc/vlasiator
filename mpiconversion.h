/*
This file is part of Vlasiator.

Copyright 2010, 2011, 2012 Finnish Meteorological Institute

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

#ifndef MPI_CONVERSION_H
#define MPI_CONVERSION_H
#include <mpi.h>
#include <cstdlib>
#include <iostream>

// Overloaded templates which should return the corresponding data type
// for some C++ native data types. For example, if float has been 
// typedef'd as Real, then MPI_Type<Real>() should return MPI_FLOAT.
template<typename T> inline MPI_Datatype MPI_Type() {
   std::cerr << "(mpiconversion.h): NULL datatype returned" << std::endl;
   return 0;
}
template<> inline MPI_Datatype MPI_Type<char>() {return MPI_CHAR;}
template<> inline MPI_Datatype MPI_Type<unsigned char>() {return MPI_UNSIGNED_CHAR;}
template<> inline MPI_Datatype MPI_Type<short int>() {return MPI_SHORT;}
template<> inline MPI_Datatype MPI_Type<unsigned short int>() {return MPI_UNSIGNED_SHORT;}
template<> inline MPI_Datatype MPI_Type<int>() {return MPI_INT;}
template<> inline MPI_Datatype MPI_Type<unsigned int>() {return MPI_UNSIGNED;}
template<> inline MPI_Datatype MPI_Type<long int>() {return MPI_LONG;}
template<> inline MPI_Datatype MPI_Type<unsigned long int>() {return MPI_UNSIGNED_LONG;}
template<> inline MPI_Datatype MPI_Type<long long int>() {return MPI_LONG_LONG;}
template<> inline MPI_Datatype MPI_Type<unsigned long long int>() {return MPI_UNSIGNED_LONG_LONG;}
template<> inline MPI_Datatype MPI_Type<float>() {return MPI_FLOAT;}
template<> inline MPI_Datatype MPI_Type<double>() {return MPI_DOUBLE;}
template<> inline MPI_Datatype MPI_Type<long double>() {return MPI_LONG_DOUBLE;}

#endif
