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

#ifndef BACKGROUNDFIELD_H
#define BACKGROUNDFIELD_H

#include "fieldfunction.hpp"
#include "integratefunction.hpp"
#include "../definitions.h"
#include "../common.h"
#include "fsgrid.hpp"

void setBackgroundField(
   const FieldFunction& bgFunction,
   FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
   bool append=false
);

void setBackgroundFieldToZero(
   FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid
);

/**
   Templated function for setting the perturbed B field to zero.
   Function is templated so it may, if necessary, be called to set the
   vector dipole correction terms inside the background field FSgrid
   object to zero (see setPerturbedField below).
*/
template<long unsigned int numFields> void setPerturbedFieldToZero(
   FsGrid< std::array<Real, numFields>, FS_STENCIL_WIDTH> & BGrid,
   int offset=fsgrids::bfield::PERBX
   ) {
   std::array<FsGridTools::FsIndex_t,3> localSize = BGrid.getLocalSize();

   #pragma omp parallel for collapse(2)
   for (FsGridTools::FsIndex_t z = 0; z < localSize[2]; ++z) {
      for (FsGridTools::FsIndex_t y = 0; y < localSize[1]; ++y) {
         for (FsGridTools::FsIndex_t x = 0; x < localSize[0]; ++x) {
            // This is still N_BFIELD (==3) instead of numFields
            for (int i = 0; i < fsgrids::bfield::N_BFIELD; ++i) {
               BGrid.get(x,y,z)->at(offset+i) = 0;
            }
         }
      }
   }
}


/**
    This function is used along with a field function to append the given function values to the
    perturbed B grid.
    A special other use case is, if using a Magnetosphere with dipole type 4 (vector dipole), when
    the corrective terms to vanish the dipole towards the inflow boundary are stored in 
    the backgroundfield FSgrid object at offset fsgrids::bgbfield::BGBXVDCORR
*/
template<long unsigned int numFields> void setPerturbedField(
   const FieldFunction& bfFunction,
   FsGrid< std::array<Real, numFields>, FS_STENCIL_WIDTH> & BGrid,
   int offset=fsgrids::bfield::PERBX,
   bool append=false) {

   using namespace std::placeholders;

   /*if we do not add a new background to the existing one we first put everything to zero*/
   if(append==false) {
      setPerturbedFieldToZero(BGrid,offset);
   }

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

   std::array<FsGridTools::FsIndex_t,3> localSize = BGrid.getLocalSize();

   // These are threaded now that the stuff around here is threadsafe
   #pragma omp parallel for collapse(2)
   for (FsGridTools::FsIndex_t z = 0; z < localSize[2]; ++z) {
      for (FsGridTools::FsIndex_t y = 0; y < localSize[1]; ++y) {
         for (FsGridTools::FsIndex_t x = 0; x < localSize[0]; ++x) {
            std::array<double, 3> start = BGrid.getPhysicalCoords(x, y, z);
            double dx[3];
            dx[0] = BGrid.DX;
            dx[1] = BGrid.DY;
            dx[2] = BGrid.DZ;

            //Face averages
            for(uint fComponent=0; fComponent<3; fComponent++){
               T3DFunction valueFunction = std::bind(bfFunction, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, (coordinate)fComponent, 0, (coordinate)0);
               BGrid.get(x,y,z)->at(offset+fComponent) += // offset defaults to fsgrids::bfield::PERBX
                  surfaceAverage(valueFunction,
                     (coordinate)fComponent,
                                 accuracy,
                                 start.data(),
                                 dx[faceCoord1[fComponent]],
                                 dx[faceCoord2[fComponent]]
                                );

	    }
	    // Derivatives or volume averages are not calculated for the perBField
	 }
      }
   }
}

#endif

