#ifndef BY_CONST_H
#define BY_CONST_H

#include "definitions.h"
#include "cell_spatial.h"

#ifndef PARGRID
	#define DCCRG_SEND_SINGLE_CELLS
	#define DCCRG_CELL_DATA_SIZE_FROM_USER
	#define DCCRG_USER_MPI_DATA_TYPE
	#include "dccrg.hpp"
#else
	#include "pargrid.h"
#endif

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

#ifndef PARGRID
void calcSimParameters(dccrg<SpatialCell>& mpiGrid, creal& t);
#else
void calcSimParameters(ParGrid<SpatialCell>& mpiGrid, creal& t);
#endif

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
 * @return The integral of the distribution function over the given six-dimensional phase-space cell.
 */
Real calcPhaseSpaceDensity(creal& x,creal& y,creal& z,creal& dx,creal& dy,creal& dz,
			   creal& vx,creal& vy,creal& vz,creal& dvx,creal& dvy,creal& dvz);

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

template<typename T> T velocityFluxX(const T& J,const T& K,const T& avg_neg,const T& avg_pos,const T* const cellParams,const T* const blockParams) {
   const T VY = blockParams[BlockParams::VYCRD] + (J+convert<T>(0.5))*blockParams[BlockParams::DVY];
   const T VZ = blockParams[BlockParams::VZCRD] + (K+convert<T>(0.5))*blockParams[BlockParams::DVZ];
   const T BY = cellParams[CellParams::BY];
   
   const T AX = blockParams[BlockParams::Q_PER_M]*(-VZ*BY);
   return convert<T>(0.5)*AX*(avg_neg + avg_pos) - convert<T>(0.5)*fabs(AX)*(avg_pos-avg_neg);
}

template<typename T> T velocityFluxY(const T& I,const T& K,const T& avg_neg,const T& avg_pos,const T* const cellParams,const T* const blockParams) {
   return convert<T>(0.0);
}

template<typename T> T velocityFluxZ(const T& I,const T& J,const T& avg_neg,const T& avg_pos,const T* const cellParams,const T* const blockParams) {
   const T VX = blockParams[BlockParams::VXCRD] + (I+convert<T>(0.5))*blockParams[BlockParams::DVX];
   const T VY = blockParams[BlockParams::VYCRD] + (J+convert<T>(0.5))*blockParams[BlockParams::DVY];
   const T BY = cellParams[CellParams::BY];

   const T AZ = blockParams[BlockParams::Q_PER_M]*(VX*BY);
   return convert<T>(0.5)*AZ*(avg_neg + avg_pos) - convert<T>(0.5)*fabs(AZ)*(avg_pos-avg_neg);
}

#endif

