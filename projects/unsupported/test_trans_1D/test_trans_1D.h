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

#ifndef TEST_TRANS_1D_H
#define TEST_TRANS_1D_H

#include "definitions.h"
#include "spatial_cell.hpp"
#include "parameters.h"
#include "projects/projects_common.h"
#include "projects/projects_vlasov_acceleration.h"

#include "dccrg.hpp"

/**
 * Initialize project. Can be used, e.g., to read in parameters from the input file
 */
bool initializeProject(void);

/** Register parameters that should be read in
 */
bool addProjectParameters(void);
/** Get the value that was read in
 */
bool getProjectParameters(void);

/*!\brief Set the fields and distribution of a cell according to the default simulation settings.
 * This is used for the NOT_SYSBOUNDARY cells and some other system boundary conditions (e.g. Outflow).
 * \param cell Pointer to the cell to set.
 */
void setProjectCell(SpatialCell* cell);

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
   ax = physicalconstants::CHARGE/physicalconstants::MASS_PROTON * (cellParams[CellParams::EX] +
      VY*(cellParams[CellParams::PERBZ]+cellParams[CellParams::BGBZ]) -
      VZ*(cellParams[CellParams::PERBY]+cellParams[CellParams::BGBY]));
   ay = physicalconstants::CHARGE/physicalconstants::MASS_PROTON * (cellParams[CellParams::EY] +
      VZ*(cellParams[CellParams::PERBX]+cellParams[CellParams::BGBX]) -
      VX*(cellParams[CellParams::PERBZ]+cellParams[CellParams::BGBZ]));
   az = physicalconstants::CHARGE/physicalconstants::MASS_PROTON * (cellParams[CellParams::EZ] +
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
   ax = physicalconstants::CHARGE/physicalconstants::MASS_PROTON * (cellParams[CellParams::EX] +
      VY*(cellParams[CellParams::PERBZ]+cellParams[CellParams::BGBZ]) -
      VZ*(cellParams[CellParams::PERBY]+cellParams[CellParams::BGBY]));
   ay = physicalconstants::CHARGE/physicalconstants::MASS_PROTON * (cellParams[CellParams::EY] +
      VZ*(cellParams[CellParams::PERBX]+cellParams[CellParams::BGBX]) -
      VX*(cellParams[CellParams::PERBZ]+cellParams[CellParams::BGBZ]));
   az = physicalconstants::CHARGE/physicalconstants::MASS_PROTON * (cellParams[CellParams::EZ] +
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
   ax = physicalconstants::CHARGE/physicalconstants::MASS_PROTON * (cellParams[CellParams::EX] +
      VY*(cellParams[CellParams::PERBZ]+cellParams[CellParams::BGBZ]) -
      VZ*(cellParams[CellParams::PERBY]+cellParams[CellParams::BGBY]));
   ay = physicalconstants::CHARGE/physicalconstants::MASS_PROTON * (cellParams[CellParams::EY] +
      VZ*(cellParams[CellParams::PERBX]+cellParams[CellParams::BGBX]) -
      VX*(cellParams[CellParams::PERBZ]+cellParams[CellParams::BGBZ]));
   az = physicalconstants::CHARGE/physicalconstants::MASS_PROTON * (cellParams[CellParams::EZ] +
      VX*(cellParams[CellParams::PERBY]+cellParams[CellParams::BGBY]) -
      VY*(cellParams[CellParams::PERBX]+cellParams[CellParams::BGBX]));
}


#endif
