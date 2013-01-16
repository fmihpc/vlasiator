/*
 * This file is part of Vlasiator.
 * 
 * Copyright 2012 Finnish Meteorological Institute
 * 
 * Vlasiator is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 3
 * as published by the Free Software Foundation.
 * 
 * Vlasiator is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "../common.h"
#include "../definitions.h"
#include "cmath"
#include "backgroundfield.h"
#include "functions.hpp"
#include "integratefunction.hpp"


void setBackgroundField(const T3DFunction bgFunction,Real* cellParams)
{
   using namespace CellParams;
   // the dipole from gumics is not threadsafe
#pragma omp critical
   {
     Dipole bgFunction;
     creal accuracy = 1e-16;
     bgFunction.initialize(8e15);
     
     creal start[3];
     creal dx[3];
     start[0] = cellParams[CellParams::XCRD];
     start[1] = cellParams[CellParams::YCRD];
     start[2] = cellParams[CellParams::ZCRD];
     dx[0] = cellParams[CellParams::DX];
     dx[1] = cellParams[CellParams::DY];
     dx[2] = cellParams[CellParams::DZ];
     
     // set Bx at negative x face
     bgFunction.setComponent(X);
     cellParams[CellParams::BGBX] =surfaceAverage(bgFunction,
						  accuracy,
						  start,
						  dx[1],
						  dx[2]);
     // set By at negative y face
     bgFunction.setComponent(Y);
     cellParams[CellParams::BGBY] =surfaceAverage(bgFunction,
						  accuracy,
						  start,
						  dx[0],
						  dx[2]);
     // set Bx at negative x face
     bgFunction.setComponent(X);
     cellParams[CellParams::BGBZ] =surfaceAverage(bgFunction,
						  accuracy,
						  start,
						  dx[0],
						  dx[1]);
     
   }

}
