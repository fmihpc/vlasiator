#ifndef DEFINITIONS_H
#define DEFINITIONS_H

typedef float Real;
//typedef float real;
typedef const float creal;
typedef const int cint;
typedef unsigned int uint;
typedef const unsigned int cuint;
typedef long unsigned int luint;
typedef const long unsigned int cluint;

typedef cuint csize;

template<typename T> T convert(const T& number) {return number;}

#endif
