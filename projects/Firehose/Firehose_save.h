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

#ifndef RIEMANN_H
#define RIEMANN_H

#include "definitions.h"
#include "spatial_cell.hpp"
#include "projects/projects_common.h"
#include "projects/projects_vlasov_acceleration.h"

#include "dccrg.hpp"


struct FirehoseParameters {
   static Real rho[2];
   static Real Tx[2];
   static Real Ty[2];
   static Real Tz[2];
   static Real Vx[2];
   static Real Vy[2];
   static Real Vz[2];
   static Real Bx;
   static Real By;
   static Real Bz;   
   static Real lambda;
   static Real amp;
   static uint nSpaceSamples;
   static uint nVelocitySamples;
};

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

void calcSimParameters(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid, creal& t, Real& dt);

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

/*!\brief Set the fields and distribution of a cell according to the default simulation settings.
 * This is used for the NOT_SYSBOUNDARY cells and some other system boundary conditions (e.g. Outflow).
 * \param cell Pointer to the cell to set.
 */
void setProjectCell(SpatialCell* cell);

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

template<typename T> T spatialFluxX(cuint& i,const T& avg_neg,const T& avg_pos,const T* const blockParams) {
   creal VX = blockParams[BlockParams::VXCRD] + i*blockParams[BlockParams::DVX];
   return convert<T>(0.5)*VX*(avg_neg+avg_pos) - convert<T>(0.5)*fabs(VX)*(avg_pos-avg_neg);
}

template<typename T> T spatialFluxY(cuint& j,const T& avg_neg,const T& avg_pos,const T* const blockParams) {
   creal VY = blockParams[BlockParams::VYCRD] + j*blockParams[BlockParams::DVY];
   return convert<T>(0.5)*VY*(avg_neg+avg_pos) - convert<T>(0.5)*fabs(VY)*(avg_pos-avg_neg);
}

template<typename T> T spatialFluxZ(cuint& k,const T& avg_neg,const T& avg_pos,const T* const blockParams) {
   creal VZ = blockParams[BlockParams::VZCRD] + k*blockParams[BlockParams::DVZ];
   return convert<T>(0.5)*VZ*(avg_neg+avg_pos) - convert<T>(0.5)*fabs(VZ)*(avg_pos-avg_neg);
}

template<typename T> T velocityFluxX(const T& j,const T& k,const T& avg_neg,const T& avg_pos,const T* const cellParams,const T* const blockParams) {
   const T VY = blockParams[BlockParams::VYCRD] + (j+convert<T>(0.5))*blockParams[BlockParams::DVY];
   const T VZ = blockParams[BlockParams::VZCRD] + (k+convert<T>(0.5))*blockParams[BlockParams::DVZ];
   const T EX = cellParams[CellParams::EX];
   const T BY = cellParams[CellParams::BY];
   const T BZ = cellParams[CellParams::BZ];
   const T AX = physicalconstants::CHARGE/physicalconstants::MASS_PROTON*(EX + VY*BZ - VZ*BY);
   return convert<T>(0.5)*AX*(avg_neg + avg_pos) - convert<T>(0.5)*fabs(AX)*(avg_pos-avg_neg);
}
                                                                        
template<typename T> T velocityFluxY(const T& i,const T& k,const T& avg_neg,const T& avg_pos,const T* const cellParams,const T* const blockParams) {
   const T VX = blockParams[BlockParams::VXCRD] + (i+convert<T>(0.5))*blockParams[BlockParams::DVX];
   const T VZ = blockParams[BlockParams::VZCRD] + (k+convert<T>(0.5))*blockParams[BlockParams::DVZ];
   const T EY = cellParams[CellParams::EY];
   const T BX = cellParams[CellParams::BX];
   const T BZ = cellParams[CellParams::BZ];
   const T AY = physicalconstants::CHARGE/physicalconstants::MASS_PROTON*(EY + VZ*BX - VX*BZ);
   return convert<T>(0.5)*AY*(avg_neg + avg_pos) - convert<T>(0.5)*fabs(AY)*(avg_pos-avg_neg);
}

template<typename T> T velocityFluxZ(const T& i,const T& j,const T& avg_neg,const T& avg_pos,const T* const cellParams,const T* const blockParams) {
   const T VX = blockParams[BlockParams::VXCRD] + (i+convert<T>(0.5))*blockParams[BlockParams::DVX];
   const T VY = blockParams[BlockParams::VYCRD] + (j+convert<T>(0.5))*blockParams[BlockParams::DVY];
   const T EZ = cellParams[CellParams::EZ];
   const T BX = cellParams[CellParams::BX];
   const T BY = cellParams[CellParams::BY];
   const T AZ = physicalconstants::CHARGE/physicalconstants::MASS_PROTON*(EZ + VX*BY - VY*BX);
   return convert<T>(0.5)*AZ*(avg_neg + avg_pos) - convert<T>(0.5)*fabs(AZ)*(avg_pos-avg_neg);
}

template<typename CELLID,typename UINT,typename REAL>
void fieldSolverBoundaryCondDerivX(const CELLID& cellID,REAL* const array,const UINT& existingCells,const UINT& nonExistingCells,creal* const derivatives,const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid) {
   return;
}

template<typename CELLID,typename UINT,typename REAL>
void fieldSolverBoundaryCondDerivY(const CELLID& cellID,REAL* const array,const UINT& existingCells,const UINT& nonExistingCells,creal* const derivatives,const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid) {
   return;
}

template<typename CELLID,typename UINT,typename REAL>
void fieldSolverBoundaryCondDerivZ(const CELLID& cellID,REAL* const array,const UINT& existingCells,const UINT& nonExistingCells,creal* const derivatives,const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid) {
   return;
}

template<typename CELLID,typename UINT,typename REAL>
REAL fieldSolverBoundaryCondBx(const CELLID& cellID,const UINT& existingCells,const UINT& nonExistingCells,const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid) {
   return 0.0;
}

template<typename CELLID,typename UINT,typename REAL>
REAL fieldSolverBoundaryCondBy(const CELLID& cellID,const UINT& existingCells,const UINT& nonExistingCells,const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid) {
   return 0.0;
}

template<typename CELLID,typename UINT,typename REAL>
REAL fieldSolverBoundaryCondBz(const CELLID& cellID,const UINT& existingCells,const UINT& nonExistingCells,const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid) {
   return 0.0;
}

template<typename CELLID,typename UINT> 
void vlasovBoundaryCondition(const CELLID& cellID,const UINT& existingCells,const UINT& nonExistingCells,const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid) {
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

