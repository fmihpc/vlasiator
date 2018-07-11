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
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */
/*
Background magnetic field class of Vlasiator.
*/

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>

#include "../common.h"
#include "integratefunction.hpp"
#include "functions.hpp"
#include "quadr.hpp"


double lineAverage(
   const T3DFunction& f1,
   coordinate line,
   double accuracy,
   const double r1[3],
   double L
) {
   double value;
   //subroutines not threadsafe 
#pragma omp critical
      {
         const double norm = 1/L;
         const double acc = accuracy*L;
         const double a = r1[line];
         const double b = r1[line] + L;
         
         switch (line) {
            case X:
            {
               T3D_fix23 f(f1,r1[1],r1[2]); 
               value= Romberg(f,a,b,acc)*norm;
            }
            break;
            case Y:
            {
               T3D_fix13 f(f1,r1[0],r1[2]); 
               value= Romberg(f,a,b,acc)*norm;
            }
            break;
            case Z: 
            {
               T3D_fix12 f(f1,r1[0],r1[1]); 
               value= Romberg(f,a,b,acc)*norm;
            }
            break;
            default:
               cerr << "*** lineAverage  is bad\n";
               value = 0.0;
            break;
         }
      }
   return value;
}


double surfaceAverage(
   const T3DFunction& f1, 
   coordinate face, double accuracy,
   const double r1[3],
   double L1,
   double L2
) {
   double value;
#pragma omp critical
   {
      const double acc = accuracy*L1*L2;
      const double norm = 1/(L1*L2);
      switch (face) {
         case X:
         {
            T3D_fix1 f(f1,r1[0]);
            value = Romberg(f, r1[1],r1[1]+L1, r1[2],r1[2]+L2, acc)*norm;
         }
         break;
         case Y:
         {
            T3D_fix2 f(f1,r1[1]);
            value = Romberg(f, r1[0],r1[0]+L1, r1[2],r1[2]+L2, acc)*norm; 
         }
         break;
         case Z:
         {
            T3D_fix3 f(f1,r1[2]);
            value = Romberg(f, r1[0],r1[0]+L1, r1[1],r1[1]+L2, acc)*norm;
         }
         break;
         default:
            cerr << "*** SurfaceAverage  is bad\n";
            exit(1);
         break;
      }
   }
   return value;
}


double volumeAverage(
   const T3DFunction& f1,
   double accuracy,
   const double r1[3],
   const double r2[3]
) {
   double value;
#pragma omp critical
   {
      const double acc = accuracy*(r2[0]-r1[0])*(r2[1]-r1[1])*(r2[2]-r1[2]);
      const double norm = 1.0/((r2[0]-r1[0])*(r2[1]-r1[1])*(r2[2]-r1[2]));
      value= Romberg(f1, r1[0],r2[0], r1[1],r2[1], r1[2],r2[2], acc)*norm;
   }
   return value;
}

