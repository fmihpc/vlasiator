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

#ifndef DISPERSION_H
#define DISPERSION_H

#include "definitions.h"
#include "spatial_cell.hpp"
#include "projects/projects_common.h"
#include "projects/projects_fieldboundary.h"
#include "projects/projects_vlasov_acceleration.h"
#include "projects/projects_vlasov_boundary.h"
#include "fieldsolver.h"

#include "dccrg.hpp"

struct dispersionParameters {
   static Real BX0;
   static Real BY0;
   static Real BZ0;
   static Real DENSITY;
   static Real TEMPERATURE;
   static Real magPertAmp;
   static Real densityPertAmp;
   static Real velocityPertAmp;
   static uint seed;
   static uint sectorSize;
   static uint nSpaceSamples;
   static uint nVelocitySamples;
} ;

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

/** Query if spatial cell parameters (of any cell) have changed and need to be 
 * recalculated. If you have a completely static case, then you can always return 
 * false here. Otherwise you need to return true so that function calcCellParameters
 * gets called for each spatial cell.
 * @param t The current value of time.
 * @return If true, the spatial cell parameters need to be recalculated.
 */
bool cellParametersChanged(creal& t);

/** Calculate parameters for the given velocity grid block.
 * Here you need to set a value for array indices: BlockParams::Q_PER_M.
 * @param blockParams Array containing the block parameters
 */
void calcBlockParameters(Real* blockParams);

/** Calculate parameters for the given spatial cell at the given time.
 * Here you need to set values for the following array indices:
 * CellParams::EX, CellParams::EY, CellParams::EZ, CellParams::BX, CellParams::BY, and CellParams::BZ.
 * 
 * The following array indices contain the coordinates of the "lower left corner" of the cell: 
 * CellParams::XCRD, CellParams::YCRD, and CellParams::ZCRD.
 * The cell size is given in the following array indices: CellParams::DX, CellParams::DY, and CellParams::DZ.
 * @param cellParams Array containing cell parameters.
 * @param t The current value of time. This is passed as a convenience. If you need more detailed information 
 * of the state of the simulation, you can read it from Parameters.
 */
void calcCellParameters(Real* cellParams,creal& t);

void calcSimParameters(dccrg::Dccrg<SpatialCell>& mpiGrid, creal& t, Real& dt);

/** Integrate the distribution function over the given six-dimensional phase-space cell.
 * @param x Starting value of the x-coordinate of the cell.
 * @param y Starting value of the y-coordinate of the cell.
 * @param z Starting value of the z-coordinate of the cell.
 * @param dx The size of the cell in x-direction.
 * @param dy The size of the cell in y-direction.
 * @param dz The size of the cell in z-direction.
 * @param vx Starting value of the vx-coordinate of the cell.
 * @param vy Starting value of the vy-coordinate of the cell.
 * @param vz Starting value of the vz-coordinate of the cell.
 * @param dvx The size of the cell in vx-direction.
 * @param dvy The size of the cell in vy-direction.
 * @param dvz The size of the cell in vz-direction.
 * @return The volume average of the distribution function in the given phase space cell.
 * The physical unit of this quantity is 1 / (m^3 (m/s)^3).
 */
Real calcPhaseSpaceDensity(creal& x,creal& y,creal& z,creal& dx,creal& dy,creal& dz,
			   creal& vx,creal& vy,creal& vz,creal& dvx,creal& dvy,creal& dvz);

/** Calculate the boundary value of volume average of distribution function. This function 
 * should calculate the value of distribution function on the other side of given phase space 
 * cell face. The coordinate direction is given with parameter crd, and positive or negative 
 * face is given with parameter negSide. For example, if crd == 0 and negSide == true, then 
 * boundary value is requested at -x face of the given phase space cell.
 * 
 * The spatial coordinate values can be calculated using cellParams. This array contains the parameters 
 * of the spatial cell which is inside the simulation domain, and you should assume that the 
 * ghost cell has the same (dx,dy,dz) values. For example, if crd == 0 and negSide == true, then 
 * the spatial coordinates of the lower left corner of the ghost cell are 
 * (cellParams[XCRD]-cellParams[DX], cellParams[YCRD], cellParams[ZCRD]).
 * 
 * The velocity coordinates can be calculated using blockParams. This array contains the parameters 
 * of the velocity grid block which is just inside the simulation domain, and you should assume that 
 * the corresponding velocity grid block of the ghost cell has the same values. For example, the 
 * velocity coordinates of the lower left corner of the velocity block are 
 * (blockParams[VXCRD]+iv*blockParams[DVX], blockParams[VYCRD]+jv*blockParams[DVY], 
 * blockParams[VZCRD]+kv*blockParams[DVZ]).
 * 
 * Note that this function does not need to be a template.
 * 
 * @param iv The vx-index of the cell in velocity block.
 * @param jv The vy-index of the cell in velocity block.
 * @param kv The vz-index of the cell in velocity block.
 * @param cellParams Array containing the spatial cell parameters.
 * @param blockParams Array containing the velocity block parameters.
 * @param avg Volume average of distribution function in the velocity block cell just inside the 
 * simulation domain. 
 * @param crd The spatial coordinate direction (0 = x, 1 = y, 2 = z).
 * @param negSide If true, then the boundary value at the negative coordinate side is requested.
 * @return Volume average of distribution function at the given boundary. The physical 
 * unit of this quantity is 1 / (m^3 (m/s)^3).
 */
template<typename T>
T calcBoundVolAvg(cuint& iv,cuint& jv,cuint& kv,const T* const cellParams,
		  const T* const blockParams,const T& avg,const int& crd,const bool& negSide) {
   return avg;
}

template<typename CELLID,typename UINT,typename REAL>
void fieldSolverBoundaryCondDerivX(const CELLID& cellID,REAL* const array,const UINT& existingCells,const UINT& nonExistingCells,creal* const derivatives,const dccrg::Dccrg<SpatialCell> & mpiGrid) {
   return;
}

template<typename CELLID,typename UINT,typename REAL>
void fieldSolverBoundaryCondDerivY(const CELLID& cellID,REAL* const array,const UINT& existingCells,const UINT& nonExistingCells,creal* const derivatives,const dccrg::Dccrg<SpatialCell> & mpiGrid) {
   return;
}

template<typename CELLID,typename UINT,typename REAL>
void fieldSolverBoundaryCondDerivZ(const CELLID& cellID,REAL* const array,const UINT& existingCells,const UINT& nonExistingCells,creal* const derivatives,const dccrg::Dccrg<SpatialCell> & mpiGrid) {
   return;
}

template<typename CELLID,typename UINT,typename REAL>
REAL fieldSolverBoundaryCondBx(const CELLID& cellID,const UINT& existingCells,const UINT& nonExistingCells,const dccrg::Dccrg<SpatialCell> & mpiGrid) {
  return 0.0;
}

template<typename CELLID,typename UINT,typename REAL>
REAL fieldSolverBoundaryCondBy(const CELLID& cellID,const UINT& existingCells,const UINT& nonExistingCells,const dccrg::Dccrg<SpatialCell> & mpiGrid) {
  return 0.0;
}

template<typename CELLID,typename UINT,typename REAL>
REAL fieldSolverBoundaryCondBz(const CELLID& cellID,const UINT& existingCells,const UINT& nonExistingCells,const dccrg::Dccrg<SpatialCell> & mpiGrid) {
  return 0.0;
}

template<typename CELLID,typename UINT> 
void vlasovBoundaryCondition(const CELLID& cellID,const UINT& existingCells,const UINT& nonExistingCells,const dccrg::Dccrg<SpatialCell> & mpiGrid) {
  return;
}

template<typename UINT,typename REAL> void calcAccFaceX(REAL& ax,REAL& ay,REAL& az,const UINT& I,const UINT& J,const UINT& K,const REAL* const cellParams,const REAL* const blockParams) {
   lorentzForceFaceX(ax,ay,az,I,J,K,cellParams,blockParams);
}
   
template<typename UINT,typename REAL> void calcAccFaceY(REAL& ax,REAL& ay,REAL& az,const UINT& I,const UINT& J,const UINT& K,const REAL* const cellParams,const REAL* const blockParams) {
   lorentzForceFaceY(ax,ay,az,I,J,K,cellParams,blockParams);
}

template<typename UINT,typename REAL> void calcAccFaceZ(REAL& ax,REAL& ay,REAL& az,const UINT& I,const UINT& J,const UINT& K,const REAL* const cellParams,const REAL* const blockParams) {
   lorentzForceFaceZ(ax,ay,az,I,J,K,cellParams,blockParams);
}

#endif
