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
   FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, 2>& BgBGrid,
   bool append) {
   
   /*if we do not add a new background to the existing one we first put everything to zero*/
   if(append==false) {
      setBackgroundFieldToZero(BgBGrid);
   }
   
   //these are doubles, as the averaging functions copied from Gumics
   //use internally doubles. In any case, it should provide more
   //accurate results also for float simulations
   double accuracy = 1e-17;
   double start[3];
   double end[3];
   double dx[3];
   unsigned int faceCoord1[3];
   unsigned int faceCoord2[3];
   
   //the coordinates of the edges face with a normal in the third coordinate direction, stored here to enable looping
   faceCoord1[0]=1;
   faceCoord2[0]=2;
   faceCoord1[1]=0;
   faceCoord2[1]=2;
   faceCoord1[2]=0;
   faceCoord2[2]=1;
   
   auto localSize = BgBGrid.getLocalSize();
   
   // Do not thread this blindly, the bgFunction.set* calls below are not thread-safe at the moment.
   for (int x = 0; x < localSize[0]; ++x) {
      for (int y = 0; y < localSize[1]; ++y) {
         for (int z = 0; z < localSize[2]; ++z) {
            std::array<double, 3> start3 = BgBGrid.getPhysicalCoords(x, y, z);
            start[0] = start3[0];
            start[1] = start3[1];
            start[2] = start3[2];
            
            dx[0] = BgBGrid.DX;
            dx[1] = BgBGrid.DY;
            dx[2] = BgBGrid.DZ;
            
            end[0]=start[0]+dx[0];
            end[1]=start[1]+dx[1];
            end[2]=start[2]+dx[2];
            
            //Face averages
            for(uint fComponent=0; fComponent<3; fComponent++){
               bgFunction.setDerivative(0);
               bgFunction.setComponent((coordinate)fComponent);
               BgBGrid.get(x,y,z)->at(fsgrids::bgbfield::BGBX+fComponent) += 
                  surfaceAverage(bgFunction,
                     (coordinate)fComponent,
                                 accuracy,
                                 start,
                                 dx[faceCoord1[fComponent]],
                                 dx[faceCoord2[fComponent]]
                                );
               
               //Compute derivatives. Note that we scale by dx[] as the arrays are assumed to contain differences, not true derivatives!
               bgFunction.setDerivative(1);
               bgFunction.setDerivComponent((coordinate)faceCoord1[fComponent]);
               BgBGrid.get(x,y,z)->at(fsgrids::bgbfield::dBGBxdy+2*fComponent) +=
                  dx[faceCoord1[fComponent]] * 
                  surfaceAverage(bgFunction, 
                     (coordinate)fComponent,
                                 accuracy,
                                 start,
                                 dx[faceCoord1[fComponent]],
                                 dx[faceCoord2[fComponent]]
                                );
               bgFunction.setDerivComponent((coordinate)faceCoord2[fComponent]);
               BgBGrid.get(x,y,z)->at(fsgrids::bgbfield::dBGBxdy+1+2*fComponent) +=
                  dx[faceCoord2[fComponent]] *
                  surfaceAverage(bgFunction,
                     (coordinate)fComponent,
                                 accuracy,
                                 start,
                                 dx[faceCoord1[fComponent]],
                                 dx[faceCoord2[fComponent]]
                                );
            }
            
            //Volume averages
            for(unsigned int fComponent=0;fComponent<3;fComponent++){
               bgFunction.setDerivative(0);
               bgFunction.setComponent((coordinate)fComponent);
               BgBGrid.get(x,y,z)->at(fsgrids::bgbfield::BGBXVOL+fComponent) += volumeAverage(bgFunction,accuracy,start,end);
               
               //Compute derivatives. Note that we scale by dx[] as the arrays are assumed to contain differences, not true derivatives!      
               bgFunction.setDerivative(1);
               bgFunction.setDerivComponent((coordinate)faceCoord1[fComponent]);
               BgBGrid.get(x,y,z)->at(fsgrids::bgbfield::dBGBXVOLdy+2*fComponent) += dx[faceCoord1[fComponent]] * volumeAverage(bgFunction,accuracy,start,end);
               bgFunction.setDerivComponent((coordinate)faceCoord2[fComponent]);
               BgBGrid.get(x,y,z)->at(fsgrids::bgbfield::dBGBXVOLdy+1+2*fComponent) += dx[faceCoord2[fComponent]] * volumeAverage(bgFunction,accuracy,start,end);
            }
         }
      }
   }
   //TODO
   //COmpute divergence and curl of volume averaged field and check that both are zero. 
}

void setBackgroundFieldToZero(
   FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, 2>& BgBGrid
) {
   auto localSize = BgBGrid.getLocalSize().data();
   
   #pragma omp parallel for collapse(3)
   for (int x = 0; x < localSize[0]; ++x) {
      for (int y = 0; y < localSize[1]; ++y) {
         for (int z = 0; z < localSize[2]; ++z) {
            for (int i = 0; i < fsgrids::bgbfield::N_BGB; ++i) {
               BgBGrid.get(x,y,z)->at(i) = 0;
            }
         }
      }
   }
}


void setPerturbedField(
   FieldFunction& bfFunction,
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2>& perBGrid,
   bool append) {
   
   /*if we do not add a new background to the existing one we first put everything to zero*/
   if(append==false) {
      setPerturbedFieldToZero(perBGrid);
   }
   
   //these are doubles, as the averaging functions copied from Gumics
   //use internally doubles. In any case, it should provide more
   //accurate results also for float simulations
   double accuracy = 1e-17;
   double start[3];
   double end[3];
   double dx[3];
   unsigned int faceCoord1[3];
   unsigned int faceCoord2[3];
   
   //the coordinates of the edges face with a normal in the third coordinate direction, stored here to enable looping
   faceCoord1[0]=1;
   faceCoord2[0]=2;
   faceCoord1[1]=0;
   faceCoord2[1]=2;
   faceCoord1[2]=0;
   faceCoord2[2]=1;
   
   auto localSize = perBGrid.getLocalSize();
   
   // Do not thread this blindly, the bfFunction.set* calls below are not thread-safe at the moment.
   for (int x = 0; x < localSize[0]; ++x) {
      for (int y = 0; y < localSize[1]; ++y) {
         for (int z = 0; z < localSize[2]; ++z) {
            std::array<double, 3> start3 = perBGrid.getPhysicalCoords(x, y, z);
            start[0] = start3[0];
            start[1] = start3[1];
            start[2] = start3[2];
            
            dx[0] = perBGrid.DX;
            dx[1] = perBGrid.DY;
            dx[2] = perBGrid.DZ;
            
            end[0]=start[0]+dx[0];
            end[1]=start[1]+dx[1];
            end[2]=start[2]+dx[2];
            
            //Face averages
            for(uint fComponent=0; fComponent<3; fComponent++){
               bfFunction.setDerivative(0);
               bfFunction.setComponent((coordinate)fComponent);
               perBGrid.get(x,y,z)->at(fsgrids::bfield::PERBX+fComponent) += 
                  surfaceAverage(bfFunction,
                     (coordinate)fComponent,
                                 accuracy,
                                 start,
                                 dx[faceCoord1[fComponent]],
                                 dx[faceCoord2[fComponent]]
                                );
               
	    }
	    // Derivatives or volume averages are not calculated for the perBField
	 }
      }
   }
}

void setPerturbedFieldToZero(
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2> & perBGrid) {
   auto localSize = perBGrid.getLocalSize().data();
   
   #pragma omp parallel for collapse(3)
   for (int x = 0; x < localSize[0]; ++x) {
      for (int y = 0; y < localSize[1]; ++y) {
         for (int z = 0; z < localSize[2]; ++z) {
            for (int i = 0; i < fsgrids::bfield::N_BFIELD; ++i) {
               perBGrid.get(x,y,z)->at(i) = 0;
            }
         }
      }
   }  
}

