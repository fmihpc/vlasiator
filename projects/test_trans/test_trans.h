/*
This file is part of Vlasiator.

Copyright 2011, 2012 Finnish Meteorological Institute

Vlasiator is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3
as published by the Free Software Foundation.

Vlasiator is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Vlasiator. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef TEST_TRANS_H
#define TEST_TRANS_H

#include <stdlib.h>

#include "../../definitions.h"
#include "../project.h"

namespace projects {
   class test_trans: public Project {
   public:
      test_trans();
      virtual ~test_trans();

            virtual bool initialize(void);
      static void addParameters(void);
      virtual void getParameters(void);
      
   protected:
      Real getDistribValue(creal& vx, creal& vy, creal& vz);
      virtual void calcCellParameters(Real* cellParams,creal& t);
      virtual Real calcPhaseSpaceDensity(
         creal& x, creal& y, creal& z,
         creal& dx, creal& dy, creal& dz,
         creal& vx, creal& vy, creal& vz,
         creal& dvx, creal& dvy, creal& dvz
      );
      Real cellPosition;
   }; // class test_trans
} // namespace projects

#endif

/*      
struct test_transParameters {
   static Real cellPosition;
} ;
*/


/* !\brief Set the fields and distribution of a cell according to the default simulation settings.
 * This is used for the NOT_SYSBOUNDARY cells and some other system boundary conditions (e.g. Outflow).
 * \param cell Pointer to the cell to set.
 */
//void setProjectCell(SpatialCell* cell);

/*
// WARNING Lorentz force not fixed in this project (cf. JXB term in the acceleration)!!!
template<typename UINT,typename REAL>
void calcAccFaceX(
   REAL& ax, REAL& ay, REAL& az,
   const UINT& I, const UINT& J, const UINT& K,
   const REAL* const cellParams,
   const REAL* const blockParams,
   const REAL* const cellBVOLDerivatives
) {
   const REAL HALF = 0.5;
   const REAL VX = blockParams[BlockParams::VXCRD] + I*blockParams[BlockParams::DVX];
   const REAL VY = blockParams[BlockParams::VYCRD] + (J+HALF)*blockParams[BlockParams::DVY];
   const REAL VZ = blockParams[BlockParams::VZCRD] + (K+HALF)*blockParams[BlockParams::DVZ];
   ax = Parameters::q_per_m * (cellParams[CellParams::EX] +
      VY*(cellParams[CellParams::PERBZ]+cellParams[CellParams::BGBZ]) -
      VZ*(cellParams[CellParams::PERBY]+cellParams[CellParams::BGBY]));
   ay = Parameters::q_per_m * (cellParams[CellParams::EY] +
      VZ*(cellParams[CellParams::PERBX]+cellParams[CellParams::BGBX]) -
      VX*(cellParams[CellParams::PERBZ]+cellParams[CellParams::BGBZ]));
   az = Parameters::q_per_m * (cellParams[CellParams::EZ] +
      VX*(cellParams[CellParams::PERBY]+cellParams[CellParams::BGBY]) -
      VY*(cellParams[CellParams::PERBX]+cellParams[CellParams::BGBX]));
}
// WARNING Lorentz force not fixed in this project (cf. JXB term in the acceleration)!!!
template<typename UINT,typename REAL>
void calcAccFaceY(
   REAL& ax, REAL& ay, REAL& az,
   const UINT& I, const UINT& J, const UINT& K,
   const REAL* const cellParams,
   const REAL* const blockParams,
   const REAL* const cellBVOLDerivatives
) {
   const REAL HALF = 0.5;
   const REAL VX = blockParams[BlockParams::VXCRD] + (I+HALF)*blockParams[BlockParams::DVX];
   const REAL VY = blockParams[BlockParams::VYCRD] + J*blockParams[BlockParams::DVY];
   const REAL VZ = blockParams[BlockParams::VZCRD] + (K+HALF)*blockParams[BlockParams::DVZ];
   ax = Parameters::q_per_m * (cellParams[CellParams::EX] +
      VY*(cellParams[CellParams::PERBZ]+cellParams[CellParams::BGBZ]) -
      VZ*(cellParams[CellParams::PERBY]+cellParams[CellParams::BGBY]));
   ay = Parameters::q_per_m * (cellParams[CellParams::EY] +
      VZ*(cellParams[CellParams::PERBX]+cellParams[CellParams::BGBX]) -
      VX*(cellParams[CellParams::PERBZ]+cellParams[CellParams::BGBZ]));
   az = Parameters::q_per_m * (cellParams[CellParams::EZ] +
      VX*(cellParams[CellParams::PERBY]+cellParams[CellParams::BGBY]) -
      VY*(cellParams[CellParams::PERBX]+cellParams[CellParams::BGBX]));
}
// WARNING Lorentz force not fixed in this project (cf. JXB term in the acceleration)!!!
template<typename UINT,typename REAL>
void calcAccFaceZ(
   REAL& ax, REAL& ay, REAL& az,
   const UINT& I, const UINT& J, const UINT& K,
   const REAL* const cellParams,
   const REAL* const blockParams,
   const REAL* const cellBVOLDerivatives
) {
   const REAL HALF = 0.5;
   const REAL VX = blockParams[BlockParams::VXCRD] + (I+HALF)*blockParams[BlockParams::DVX];
   const REAL VY = blockParams[BlockParams::VYCRD] + (J+HALF)*blockParams[BlockParams::DVY];
   const REAL VZ = blockParams[BlockParams::VZCRD] + K*blockParams[BlockParams::DVZ];
   ax = Parameters::q_per_m * (cellParams[CellParams::EX] +
      VY*(cellParams[CellParams::PERBZ]+cellParams[CellParams::BGBZ]) -
      VZ*(cellParams[CellParams::PERBY]+cellParams[CellParams::BGBY]));
   ay = Parameters::q_per_m * (cellParams[CellParams::EY] +
      VZ*(cellParams[CellParams::PERBX]+cellParams[CellParams::BGBX]) -
      VX*(cellParams[CellParams::PERBZ]+cellParams[CellParams::BGBZ]));
   az = Parameters::q_per_m * (cellParams[CellParams::EZ] +
      VX*(cellParams[CellParams::PERBY]+cellParams[CellParams::BGBY]) -
      VY*(cellParams[CellParams::PERBX]+cellParams[CellParams::BGBX]));
}
*/

