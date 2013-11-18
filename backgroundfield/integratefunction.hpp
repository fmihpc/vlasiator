/*
Background magnetic field class of Vlasiator.

Copyright 1997, 1998, 1999, 2000, 2001, 2010, 2011 Finnish Meteorological Institute
*/

#ifndef INTEGRATEFIELDFUNCTION_HPP
#define INTEGRATEFIELDFUNCTION_HPP

#include "ode.hpp"
#include "quadr.hpp"
#include "functions.hpp"
/*!
  Average of f1 along a coordinate-aligned line starting from r1,
  having length L (can be negative) and proceeding to line'th coordinate
*/
double lineAverage(
   const T3DFunction& f1,
   coordinate line,
   double accuracy,
   const double r1[3],
   double L
);

/*!
  Average of f1 along a rectangular coordinate-aligned surface
  which is orthogonal to face'th coordinate, has lower left corner at r1,
  and surface side lengths (positive) equal to L1,L2 (either yz, xz or xy,
  depending on d).
*/
double surfaceAverage(
   const T3DFunction& f1,
   coordinate face, double accuracy,
   const double r1[3],
   double L1,
   double L2
);

/*!
  Average of f1 over a rectangular coordinate-aligned volume
  having lower left corner at r1 and upper right corner at r2.
*/
double volumeAverage(
   const T3DFunction& f1,
   double accuracy,
   const double r1[3],
   const double r2[3]
);
#endif

