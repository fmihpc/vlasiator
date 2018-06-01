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

#include "../common.h"
#include "../definitions.h"
#include "../parameters.h"
#include "cmath"
#include "backgroundfield.h"
#include "fieldfunction.hpp"
#include "integratefunction.hpp"

//FieldFunction should be initialized
void setBackgroundField(
   FieldFunction& bgFunction,
   Real* cellParams,
   Real* faceDerivatives,
   Real* volumeDerivatives,
   bool append) {
   using namespace CellParams;
   using namespace fieldsolver;
   using namespace bvolderivatives;
   
   //these are doubles, as the averaging functions copied from Gumics
   //use internally doubles. In any case, it should provide more
   //accurate results also for float simulations
   double accuracy = 1e-17;
   double start[3];
   double end[3];
   double dx[3];
   unsigned int faceCoord1[3];
   unsigned int faceCoord2[3];


   start[0] = cellParams[CellParams::XCRD];
   start[1] = cellParams[CellParams::YCRD];
   start[2] = cellParams[CellParams::ZCRD];

   dx[0] = cellParams[CellParams::DX];
   dx[1] = cellParams[CellParams::DY];
   dx[2] = cellParams[CellParams::DZ];

   end[0]=start[0]+dx[0];
   end[1]=start[1]+dx[1];
   end[2]=start[2]+dx[2];
   
   //the coordinates of the edges face with a normal in the third coordinate direction, stored here to enable looping
   faceCoord1[0]=1;
   faceCoord2[0]=2;
   faceCoord1[1]=0;
   faceCoord2[1]=2;
   faceCoord1[2]=0;
   faceCoord2[2]=1;

   /*if we do not add a new background to the existing one we first put everything to zero*/
   if(append==false) {
      setBackgroundFieldToZero(cellParams, faceDerivatives, volumeDerivatives);
   }
   
   //Face averages
   for(unsigned int fComponent=0;fComponent<3;fComponent++){
      bgFunction.setDerivative(0);
      bgFunction.setComponent((coordinate)fComponent);
      cellParams[CellParams::BGBX+fComponent] += 
      surfaceAverage(
         bgFunction,
         (coordinate)fComponent,
         accuracy,
         start,
         dx[faceCoord1[fComponent]],
         dx[faceCoord2[fComponent]]
      );
      
      //Compute derivatives. Note that we scale by dx[] as the arrays are assumed to contain differences, not true derivatives!
      bgFunction.setDerivative(1);
      bgFunction.setDerivComponent((coordinate)faceCoord1[fComponent]);
      faceDerivatives[fieldsolver::dBGBxdy+2*fComponent] +=
         dx[faceCoord1[fComponent]]*
         surfaceAverage(bgFunction,(coordinate)fComponent,accuracy,start,dx[faceCoord1[fComponent]],dx[faceCoord2[fComponent]]);
      bgFunction.setDerivComponent((coordinate)faceCoord2[fComponent]);
      faceDerivatives[fieldsolver::dBGBxdy+1+2*fComponent] +=
         dx[faceCoord2[fComponent]]*
         surfaceAverage(bgFunction,(coordinate)fComponent,accuracy,start,dx[faceCoord1[fComponent]],dx[faceCoord2[fComponent]]);
   }

   //Volume averages
   for(unsigned int fComponent=0;fComponent<3;fComponent++){
      bgFunction.setDerivative(0);
      bgFunction.setComponent((coordinate)fComponent);
      cellParams[CellParams::BGBXVOL+fComponent] += volumeAverage(bgFunction,accuracy,start,end);

      //Compute derivatives. Note that we scale by dx[] as the arrays are assumed to contain differences, not true derivatives!      
      bgFunction.setDerivative(1);
      bgFunction.setDerivComponent((coordinate)faceCoord1[fComponent]);
      volumeDerivatives[bvolderivatives::dBGBXVOLdy+2*fComponent] +=  dx[faceCoord1[fComponent]]*volumeAverage(bgFunction,accuracy,start,end);
      bgFunction.setDerivComponent((coordinate)faceCoord2[fComponent]);
      volumeDerivatives[bvolderivatives::dBGBXVOLdy+1+2*fComponent] += dx[faceCoord2[fComponent]]*volumeAverage(bgFunction,accuracy,start,end);
   }

   //TODO
   //COmpute divergence and curl of volume averaged field and check that both are zero. 
}

void setBackgroundFieldToZero(
   Real* cellParams,
   Real* faceDerivatives,
   Real* volumeDerivatives
) {
   using namespace CellParams;
   using namespace fieldsolver;
   using namespace bvolderivatives;
   
   //Face averages
   for(unsigned int fComponent=0;fComponent<3;fComponent++){
      cellParams[CellParams::BGBX+fComponent] = 0.0;
      faceDerivatives[fieldsolver::dBGBxdy+2*fComponent] = 0.0;
      faceDerivatives[fieldsolver::dBGBxdy+1+2*fComponent] = 0.0;
   }
   
   //Volume averages
   for(unsigned int fComponent=0;fComponent<3;fComponent++){
      cellParams[CellParams::BGBXVOL+fComponent] = 0.0;
      volumeDerivatives[bvolderivatives::dBGBXVOLdy+2*fComponent] = 0.0;
      volumeDerivatives[bvolderivatives::dBGBXVOLdy+1+2*fComponent] =0.0;
   }
   
}
