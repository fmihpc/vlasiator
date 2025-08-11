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

// clang-format off
#include "backgroundfield.h"
#include "../common.h"
#include "../definitions.h"
#include "../parameters.h"
#include "cmath"
#include "phiprof.hpp"
// clang-format on

// FieldFunction should be initialized
void setBackgroundField(const FieldFunction& bgFunction, fsgrids::bgbspan bgb,
                        fsgrids::technicalspan technical, FieldSolverGrid &fsgrid, bool append) {
   const auto& gridSpacing = fsgrid.getGridSpacing();

   /*if we do not add a new background to the existing one we first put everything to zero*/
   if (append == false) {
      setBackgroundFieldToZero(fsgrid, technical, bgb);
   }
   const size_t numCells = fsgrid.getNumCells();
   phiprof::Timer bgTimer{"set Background field"};
   {
      // these are doubles, as the averaging functions copied from Gumics
      // use internally doubles. In any case, it should provide more
      // accurate results also for float simulations
      const double accuracy = 1e-17;
      unsigned int faceCoord1[3];
      unsigned int faceCoord2[3];

      // the coordinates of the edges face with a normal in the third coordinate direction, stored here to enable
      // looping
      faceCoord1[0] = 1;
      faceCoord2[0] = 2;
      faceCoord1[1] = 0;
      faceCoord2[1] = 2;
      faceCoord1[2] = 0;
      faceCoord2[2] = 1;

      int loopTopId{phiprof::initializeTimer("loop-top")};
      int loopFaceId{phiprof::initializeTimer("loop-face-averages")};
      int loopVolumeId{phiprof::initializeTimer("loop-volume-averages")};

// These are threaded now that the dipole field is threadsafe
      fsgrid.parallel_for([](int timerId) -> phiprof::Timer { return phiprof::Timer{timerId}; },
                          phiprof::initializeTimer("setBackgroundField-loop"), technical,
                          [& /*=,&fsgrid,&bgFunction*/](const fsgrid::Coordinates &coordinates, const fsgrid::FsStencil& stencil, cuint sysBoundaryFlag, cuint sysBoundaryLayer) {
         const std::array<Real, 3> start = coordinates.getPhysicalCoords(stencil.i, stencil.j, stencil.k);
         const std::array end = {
            start[0] + gridSpacing[0],
            start[1] + gridSpacing[1],
            start[2] + gridSpacing[2],
         };

         auto& field = bgb[stencil.ooo()];
         // Face averages
         for (uint fComponent = 0; fComponent < 3; fComponent++) {
            T3DFunction valueFunction =
               std::bind(bgFunction, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3,
                  (coordinate)fComponent, 0, (coordinate)0);
            field[fsgrids::bgbfield::BGBX + fComponent] +=
               surfaceAverage(valueFunction, (coordinate)fComponent, accuracy, start,
                  gridSpacing[faceCoord1[fComponent]], gridSpacing[faceCoord2[fComponent]]);

           // Compute derivatives. Note that we scale by gridSpacing[] as the arrays are assumed to contain
           // differences, not true derivatives!
           T3DFunction derivFunction1 =
              std::bind(bgFunction, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3,
                 (coordinate)fComponent, 1, (coordinate)faceCoord1[fComponent]);
           field[fsgrids::bgbfield::dBGBxdy + 2 * fComponent] +=
              gridSpacing[faceCoord1[fComponent]] *
                 surfaceAverage(derivFunction1, (coordinate)fComponent, accuracy, start,
                    gridSpacing[faceCoord1[fComponent]], gridSpacing[faceCoord2[fComponent]]);

           T3DFunction derivFunction2 =
              std::bind(bgFunction, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3,
                 (coordinate)fComponent, 1, (coordinate)faceCoord2[fComponent]);
           field[fsgrids::bgbfield::dBGBxdy + 1 + 2 * fComponent] +=
              gridSpacing[faceCoord2[fComponent]] *
                 surfaceAverage(derivFunction2, (coordinate)fComponent, accuracy, start,
                    gridSpacing[faceCoord1[fComponent]], gridSpacing[faceCoord2[fComponent]]);
         }

         // Volume averages
         for (uint fComponent = 0; fComponent < 3; fComponent++) {
            T3DFunction valueFunction =
               std::bind(bgFunction, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3,
                  (coordinate)fComponent, 0, (coordinate)0);
            field[fsgrids::bgbfield::BGBXVOL + fComponent] += volumeAverage(valueFunction, accuracy, start, end);

            // Compute derivatives. Note that we scale by gridSpacing[] as the arrays are assumed to contain
            // differences, not true derivatives!
            for (uint dComponent = 0; dComponent < 3; dComponent++) {
               T3DFunction derivFunction =
                  std::bind(bgFunction, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3,
                     (coordinate)fComponent, 1, (coordinate)dComponent);
               field[fsgrids::bgbfield::dBGBXVOLdx + 3 * fComponent + dComponent] +=
                  gridSpacing[dComponent] * volumeAverage(derivFunction, accuracy, start, end);
            }
         }
      });
   }
   bgTimer.stop(numCells, "Spatial Cells");
   // TODO
   // Compute divergence and curl of volume averaged field and check that both are zero.
}

void setBackgroundFieldToZero(
   FieldSolverGrid &fsgrid,
   fsgrids::technicalspan technical,
   fsgrids::bgbspan bgb
) {
   fsgrid.parallel_for([](int timerId) -> phiprof::Timer { return phiprof::Timer{timerId}; },
   phiprof::initializeTimer("setBackgroundFieldToZero"), technical,
   [=](const fsgrid::Coordinates &coordinates, const fsgrid::FsStencil& stencil, cuint sysBoundaryFlag, cuint sysBoundaryLayer) {
      for (size_t i = 0; i < bgb[stencil.ooo()].size(); i++) {
         bgb[stencil.ooo()][i] == 0.0;
      }
   });
}
