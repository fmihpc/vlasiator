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
#include "phiprof.hpp"

//FieldFunction should be initialized
void setBackgroundField(
   const FieldFunction& bgFunction,
   FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
   bool append
   ) {
   using namespace std::placeholders;

   /*if we do not add a new background to the existing one we first put everything to zero*/
   if(append==false) {
      setBackgroundFieldToZero(BgBGrid);
   }
   const FsGridTools::FsIndex_t* gridDims = &BgBGrid.getLocalSize()[0];
   const size_t N_cells = gridDims[0]*gridDims[1]*gridDims[2];
   phiprof::Timer bgTimer {"set Background field"};
   {
      //these are doubles, as the averaging functions copied from Gumics
      //use internally doubles. In any case, it should provide more
      //accurate results also for float simulations
      const double accuracy = 1e-17;
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

      int loopTopId {phiprof::initializeTimer("loop-top")};
      int loopFaceId {phiprof::initializeTimer("loop-face-averages")};
      int loopVolumeId {phiprof::initializeTimer("loop-volume-averages")};

      // These are threaded now that the dipole field is threadsafe
      #pragma omp parallel for collapse(2)
      for (FsGridTools::FsIndex_t z = 0; z < localSize[2]; ++z) {
         for (FsGridTools::FsIndex_t y = 0; y < localSize[1]; ++y) {
            for (FsGridTools::FsIndex_t x = 0; x < localSize[0]; ++x) {
               phiprof::Timer loopTopTimer {loopTopId};
               std::array<double, 3> start = BgBGrid.getPhysicalCoords(x, y, z);
               double dx[3];
               dx[0] = BgBGrid.DX;
               dx[1] = BgBGrid.DY;
               dx[2] = BgBGrid.DZ;
               double end[3];
               end[0]=start[0]+dx[0];
               end[1]=start[1]+dx[1];
               end[2]=start[2]+dx[2];
               loopTopTimer.stop();

               phiprof::Timer loopFaceTimer {loopFaceId};
               //Face averages
               for(uint fComponent=0; fComponent<3; fComponent++){
                  T3DFunction valueFunction = std::bind(bgFunction, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, (coordinate)fComponent, 0, (coordinate)0);
                  BgBGrid.get(x,y,z)->at(fsgrids::bgbfield::BGBX+fComponent) +=
                     surfaceAverage(valueFunction,
                        (coordinate)fComponent,
                                    accuracy,
                                    start.data(),
                                    dx[faceCoord1[fComponent]],
                                    dx[faceCoord2[fComponent]]
                                   );

                  //Compute derivatives. Note that we scale by dx[] as the arrays are assumed to contain differences, not true derivatives!
                  T3DFunction derivFunction1 = std::bind(bgFunction, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, (coordinate)fComponent, 1, (coordinate)faceCoord1[fComponent]);
                  BgBGrid.get(x,y,z)->at(fsgrids::bgbfield::dBGBxdy+2*fComponent) +=
                     dx[faceCoord1[fComponent]] *
                     surfaceAverage(derivFunction1,
                        (coordinate)fComponent,
                                    accuracy,
                                    start.data(),
                                    dx[faceCoord1[fComponent]],
                                    dx[faceCoord2[fComponent]]
                                   );

                  T3DFunction derivFunction2 = std::bind(bgFunction, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, (coordinate)fComponent, 1, (coordinate)faceCoord2[fComponent]);
                  BgBGrid.get(x,y,z)->at(fsgrids::bgbfield::dBGBxdy+1+2*fComponent) +=
                     dx[faceCoord2[fComponent]] *
                     surfaceAverage(derivFunction2,
                        (coordinate)fComponent,
                                    accuracy,
                                    start.data(),
                                    dx[faceCoord1[fComponent]],
                                    dx[faceCoord2[fComponent]]
                                   );
               }
               loopFaceTimer.stop();

               phiprof::Timer loopVolumeTimer {loopVolumeId};
               //Volume averages
               for(uint fComponent=0;fComponent<3;fComponent++){
                  T3DFunction valueFunction = std::bind(bgFunction, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, (coordinate)fComponent, 0, (coordinate)0);
                  BgBGrid.get(x,y,z)->at(fsgrids::bgbfield::BGBXVOL+fComponent) += volumeAverage(valueFunction,accuracy,start.data(),end);
                  
                  //Compute derivatives. Note that we scale by dx[] as the arrays are assumed to contain differences, not true derivatives!
                  for(uint dComponent=0;dComponent<3;dComponent++){
                     T3DFunction derivFunction = std::bind(bgFunction, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, (coordinate)fComponent, 1, (coordinate)dComponent);
                     BgBGrid.get(x,y,z)->at(fsgrids::bgbfield::dBGBXVOLdx+3*fComponent+dComponent) += dx[dComponent] * volumeAverage(derivFunction,accuracy,start.data(),end);
                  }
               }
               loopVolumeTimer.stop();
            }
         }
      }

   }
   bgTimer.stop(N_cells, "Spatial Cells");
   //TODO
   //Compute divergence and curl of volume averaged field and check that both are zero.
}

void setBackgroundFieldToZero(
   FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid
) {
   auto localSize = BgBGrid.getLocalSize().data();

   #pragma omp parallel for collapse(2)
   for (FsGridTools::FsIndex_t z = 0; z < localSize[2]; ++z) {
      for (FsGridTools::FsIndex_t y = 0; y < localSize[1]; ++y) {
         for (FsGridTools::FsIndex_t x = 0; x < localSize[0]; ++x) {
            for (FsGridTools::FsIndex_t i = 0; i < fsgrids::bgbfield::N_BGB; ++i) {
               BgBGrid.get(x,y,z)->at(i) = 0;
            }
         }
      }
   }
}

