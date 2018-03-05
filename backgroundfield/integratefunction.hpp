/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://www.physics.helsinki.fi/vlasiator/
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */
/*
Background magnetic field class of Vlasiator.
*/

#ifndef INTEGRATEFIELDFUNCTION_HPP
#define INTEGRATEFIELDFUNCTION_HPP


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

