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

#include "../common.h"
#include "../definitions.h"
#include "cmath"
#include "backgroundfield.h"
#include "B0.hpp"


void setDipoleBackgroundBNegX(
   TB0& background_B,
   creal start_x, creal start_y, creal start_z,
   creal end_x, creal end_y, creal end_z,
   Real& Bx, Real &By, Real& Bz
);

void setDipoleBackgroundBNegY(
   TB0& background_B,
   creal start_x, creal start_y, creal start_z,
   creal end_x, creal end_y, creal end_z,
   Real& Bx, Real &By, Real& Bz
);

void setDipoleBackgroundBNegZ(
   TB0& background_B,
   creal start_x, creal start_y, creal start_z,
   creal end_x, creal end_y, creal end_z,
   Real& Bx, Real &By, Real& Bz
);


void setDipole(Real* cellParams)
{
   using namespace CellParams;
// the dipole from gumics is not threadsafe
#pragma omp critical
   {
      TB0 background_B;
      background_B.set_dipole_moment(8e15);
      background_B.initialize();

      creal
         start_x = cellParams[CellParams::XCRD],
         start_y = cellParams[CellParams::YCRD],
         start_z = cellParams[CellParams::ZCRD],
         end_x = start_x + cellParams[CellParams::DX],
         end_y = start_y + cellParams[CellParams::DY],
         end_z = start_z + cellParams[CellParams::DZ];
      
      Real Bx, By, Bz;

      // set Bx at negative x face
      setDipoleBackgroundBNegX(
         background_B,
         start_x, start_y, start_z,
         end_x, end_y, end_z,
         Bx, By, Bz
                               );
      cellParams[CellParams::BGBX] = Bx;

      // By at -y face
      setDipoleBackgroundBNegY(
         background_B,
         start_x, start_y, start_z,
         end_x, end_y, end_z,
         Bx, By, Bz
                               );
      cellParams[CellParams::BGBY] = By;

      // Bz at -z
      setDipoleBackgroundBNegZ(
         background_B,
         start_x, start_y, start_z,
         end_x, end_y, end_z,
         Bx, By, Bz
                               );
      cellParams[CellParams::BGBZ] = Bz;
   }
}

/*!
Assigns background B face average in -x direction.

Start_? and end_? correspond to the cell's vertices
at minimum and maximum coordinates respectively.
*/
void setDipoleBackgroundBNegX(
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
Y version of setDipoleBackgroundBNegx(...).
*/
void setDipoleBackgroundBNegY(
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
Z version of setDipoleBackgroundBNegx(...).
*/
void setDipoleBackgroundBNegZ(
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

