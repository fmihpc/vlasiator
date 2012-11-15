/*
This file is part of Vlasiator.

Copyright 2012 Finnish Meteorological Institute

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

#include "../definitions.h"
#include "cmath"
#include "backgroundfield.h"
#include "B0.hpp"

void dipole(
   creal x, creal y, creal z,
   Real& Bx, Real &By, Real& Bz
) {
   creal k_0 = 8.0e15; // Wb m
   Real r = sqrt(x*x + y*y + z*z); // radial
   Real theta = atan2(sqrt(x*x + y*y), z); // polar
   Real phi = atan2(y, x); // azimuthal
   
   Bx = sin(theta) * cos(theta) * cos(phi);
   By = sin(theta) * cos(theta) * sin(phi);
   Bz = cos(theta)*cos(theta) - 1.0 / 3.0;
   
   Bx *= 3.0 * k_0 / (r*r*r);
   By *= 3.0 * k_0 / (r*r*r);
   Bz *= 3.0 * k_0 / (r*r*r);
}



/*!
Assigns background B face average in -x direction.

Start_? and end_? correspond to the cell's vertices
at minimum and maximum coordinates respectively.
*/
void set_background_B_neg_x(
   TB0& background_B,
   creal start_x, creal start_y, creal start_z,
   creal end_x, creal end_y, creal end_z,
   Real& Bx, Real &By, Real& Bz
) {
   const double start[3] = {
      start_x,
      start_y,
      start_z
   };

   const double end[3] = {
      end_x,
      end_y,
      end_z
   };

   background_B.BackgroundSurfaceAverage(
      start,
      0,
      end[1] - start[1],
      end[2] - start[2],
      Bx,
      By,
      Bz
   );
}

/*!
Y version of set_background_B_neg_x(...).
*/
void set_background_B_neg_y(
   TB0& background_B,
   creal start_x, creal start_y, creal start_z,
   creal end_x, creal end_y, creal end_z,
   Real& Bx, Real &By, Real& Bz
) {
   const double start[3] = {
      start_x,
      start_y,
      start_z
   };

   const double end[3] = {
      end_x,
      end_y,
      end_z
   };

   background_B.BackgroundSurfaceAverage(
      start,
      1,
      end[0] - start[0],
      end[2] - start[2],
      Bx,
      By,
      Bz
   );
}

/*!
Z version of set_background_B_neg_x(...).
*/
void set_background_B_neg_z(
   TB0& background_B,
   creal start_x, creal start_y, creal start_z,
   creal end_x, creal end_y, creal end_z,
   Real& Bx, Real &By, Real& Bz
) {
   const double start[3] = {
      start_x,
      start_y,
      start_z
   };

   const double end[3] = {
      end_x,
      end_y,
      end_z
   };

   background_B.BackgroundSurfaceAverage(
      start,
      2,
      end[0] - start[0],
      end[1] - start[1],
      Bx,
      By,
      Bz
   );
}

