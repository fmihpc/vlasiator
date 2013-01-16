/*
This file is part of Vlasiator.

Copyright 1997, 1998, 1999, 2000, 2001, 2010, 2011 Finnish Meteorological Institute
*/

#ifndef QUADR_HPP
#define QUADR_HPP

#include <iostream>
#include "functions.hpp"
using namespace std;

/*
  1D,2D,3D Romberg non-singular integration a'la Numerical Recipes.
  Integration bounds must be constants.
  Header file is quadr.H.
  Test program is tstquadr.C.
  Same in Mathematica is tstquadr.ma.
*/

/*
  The iterations are stopped when the results changes by less than absacc.
*/



double Romberg(const T1DFunction& func, double a, double b, double absacc);
double Romberg(const T2DFunction& func, double a, double b, double c, double d, double absacc);
double Romberg(const T3DFunction& func, double a, double b, double c, double d, double e, double f, double absacc);

#endif
