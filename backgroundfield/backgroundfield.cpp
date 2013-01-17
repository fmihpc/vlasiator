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
#include "fieldfunction.hpp"
#include "integratefunction.hpp"

//FieldFunction should be initialized
void setBackgroundField(FieldFunction& bgFunction,Real* cellParams)
{
   using namespace CellParams;
   //these are doubles, as the averaging functions copied from Gumics
   //use internally doubles. In any case, it should provide more
   //accurate results also for float simulations
   double accuracy = 1e-17;     
   double start[3];
   double end[3];
   double dx[3];

   double volB[3];
   double dvolBdr[3][3];

   start[0] = cellParams[CellParams::XCRD];
   start[1] = cellParams[CellParams::YCRD];
   start[2] = cellParams[CellParams::ZCRD];

   dx[0] = cellParams[CellParams::DX];
   dx[1] = cellParams[CellParams::DY];
   dx[2] = cellParams[CellParams::DZ];

   end[0]=start[0]+dx[0];
   end[1]=start[1]+dx[1];
   end[2]=start[2]+dx[2];
   
   
   //we are not computing derivatives
   bgFunction.setDerivative(0);
      
   // set Bx at negative x face
   bgFunction.setComponent(X);
   cellParams[CellParams::BGBX] =surfaceAverage(bgFunction,X,accuracy,start,dx[1],dx[2]);

// set By at negative y face 
   bgFunction.setComponent(Y);
   cellParams[CellParams::BGBY] =surfaceAverage(bgFunction,Y,accuracy,start,dx[0],dx[2]);

// set Bz at negative z face 
   bgFunction.setComponent(Z);
   cellParams[CellParams::BGBZ] =surfaceAverage(bgFunction,Z,accuracy,start,dx[0],dx[1]);
   

   //Volume averages
   for(unsigned int fComponent=0;fComponent<3;fComponent++){
      bgFunction.setDerivative(0);
      bgFunction.setComponent((coordinate)fComponent);      
      volB[fComponent] =volumeAverage(bgFunction,accuracy,start,end);
      
      bgFunction.setDerivative(1);
      for(unsigned int dComponent=0;dComponent<3;dComponent++){
         bgFunction.setDerivComponent((coordinate)dComponent);      
         dvolBdr[fComponent][dComponent]=volumeAverage(bgFunction,accuracy,start,end);
      }
   }
}


