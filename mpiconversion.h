#ifndef MPI_CONVERSION_H
#define MPI_CONVERSION_H
#include <cstdlib>
#include <iostream>
#include <mpi.h>

// Overloaded templates which should return the corresponding data type
// for some C++ native data types. For example, if float has been 
// typedef'd as Real, then MPI_Type<Real>() should return MPI_FLOAT.
template<typename T> inline MPI_Datatype MPI_Type() {
   std::cerr << "NULL datatype returned" << std::endl;
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
