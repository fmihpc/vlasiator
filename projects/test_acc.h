#ifndef HARM1D_H
#define HARM1D_H

#include "definitions.h"
#include "cell_spatial.h"
#include "projects/projects_common.h"
#include "projects/projects_fieldboundary.h"
#include "projects/projects_vlasov_acceleration.h"
#include "projects/projects_vlasov_boundary.h"
#include "fieldsolver.h"
#include "arrayallocator.h"

#ifndef PARGRID
	#define DCCRG_SEND_SINGLE_CELLS
	#define DCCRG_CELL_DATA_SIZE_FROM_USER
	#define DCCRG_USER_MPI_DATA_TYPE
	#include "dccrg.hpp"
#else
	#include "pargrid.h"
#endif

/**
 * Initialize project. Can be used, e.g., to read in parameters from the input file
 */
bool initializeProject(void);

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

void calcSimParameters(
	#ifndef PARGRID
	dccrg::Dccrg<SpatialCell>& mpiGrid,
	#else
	ParGrid<SpatialCell>& mpiGrid
	#endif
	creal& t,
	Real& dt);

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
template<typename T> T calcBoundVolAvg(
	cuint& /*iv*/, cuint& /*jv*/, cuint& /*kv*/,
	const T* const /*cellParams*/,
	const T* const /*blockParams*/,
	const T& avg,
	const int& /*crd*/,
	const bool& /*negSide*/
) {
   return avg;
}

template<typename CELLID,typename UINT,typename REAL> void fieldSolverBoundaryCondDerivX(
	const CELLID& cellID,
	REAL* const array,
	const UINT& existingCells,
	const UINT& nonExistingCells,
	creal* const derivatives,
	#ifdef PARGRID
	const ParGrid<SpatialCell>& mpiGrid
	#else
	const dccrg::Dccrg<SpatialCell>& mpiGrid
	#endif
) {
   fieldSolverBoundarySetValueDerivX(cellID,array,existingCells,nonExistingCells,derivatives,mpiGrid,convert<REAL>(0.0));
}

template<typename CELLID,typename UINT,typename REAL> void fieldSolverBoundaryCondDerivY(
	const CELLID& cellID,
	REAL* const array,
	const UINT& existingCells,
	const UINT& nonExistingCells,
	creal* const derivatives,
	#ifdef PARGRID
	const ParGrid<SpatialCell>& mpiGrid
	#else
	const dccrg::Dccrg<SpatialCell>& mpiGrid
	#endif
) {
   fieldSolverBoundarySetValueDerivY(cellID,array,existingCells,nonExistingCells,derivatives,mpiGrid,convert<REAL>(0.0));
}

template<typename CELLID,typename UINT,typename REAL> void fieldSolverBoundaryCondDerivZ(
	const CELLID& cellID,
	REAL* const array,
	const UINT& existingCells,
	const UINT& nonExistingCells,
	creal* const derivatives,
	#ifdef PARGRID
	const ParGrid<SpatialCell>& mpiGrid
	#else
	const dccrg::Dccrg<SpatialCell>& mpiGrid
	#endif
) {
   fieldSolverBoundarySetValueDerivZ(cellID,array,existingCells,nonExistingCells,derivatives,mpiGrid,convert<REAL>(0.0));
}

template<typename CELLID,typename UINT,typename REAL> REAL fieldSolverBoundaryCondBx(
	const CELLID& cellID,
	const UINT& existingCells,
	const UINT& nonExistingCells,
	#ifdef PARGRID
	const ParGrid<SpatialCell>& mpiGrid
	#else
	const dccrg::Dccrg<SpatialCell>& mpiGrid
	#endif
) {
   return fieldBoundaryCopyFromExistingFaceNbrBx<CELLID,UINT,REAL>(cellID,existingCells,nonExistingCells,mpiGrid);
}

template<typename CELLID,typename UINT,typename REAL> REAL fieldSolverBoundaryCondBy(
	const CELLID& cellID,
	const UINT& existingCells,
	const UINT& nonExistingCells,
	#ifdef PARGRID
	const ParGrid<SpatialCell>& mpiGrid
	#else
	const dccrg::Dccrg<SpatialCell>& mpiGrid
	#endif
) {
   return fieldBoundaryCopyFromExistingFaceNbrBy<CELLID,UINT,REAL>(cellID,existingCells,nonExistingCells,mpiGrid);
}

template<typename CELLID,typename UINT,typename REAL> REAL fieldSolverBoundaryCondBz(
	const CELLID& cellID,
	const UINT& existingCells,
	const UINT& nonExistingCells,
	#ifdef PARGRID
	const ParGrid<SpatialCell>& mpiGrid
	#else
	const dccrg::Dccrg<SpatialCell>& mpiGrid
	#endif
) {
   return fieldBoundaryCopyFromExistingFaceNbrBz<CELLID,UINT,REAL>(cellID,existingCells,nonExistingCells,mpiGrid);
}

template<typename CELLID,typename UINT> void vlasovBoundaryCondition(
	const CELLID& cellID,
	const UINT& existingCells,
	const UINT& nonExistingCells,
	#ifdef PARGRID
	const ParGrid<SpatialCell>& mpiGrid
	#else
	const dccrg::Dccrg<SpatialCell>& mpiGrid
	#endif
) {
   vlasovBoundaryCopyFromExistingFaceNbr(cellID,existingCells,nonExistingCells,mpiGrid);
}

template<typename UINT,typename REAL> void calcAccFaceX(
	REAL& ax, REAL& ay, REAL& az,
	const UINT& I, const UINT& J, const UINT& K,
	const REAL* const cellParams,
	const REAL* const blockParams
) {
   ax = blockParams[BlockParams::VXCRD] + (I+convert<REAL>(0.5))*blockParams[BlockParams::DVX];
   ay = blockParams[BlockParams::VYCRD] + (J+convert<REAL>(0.5))*blockParams[BlockParams::DVY];
   az = blockParams[BlockParams::VZCRD] + (K+convert<REAL>(0.5))*blockParams[BlockParams::DVZ];
   
   const REAL SPEED = convert<REAL>(0.5) + convert<REAL>(0.5)*(convert<REAL>(1.0) / convert<REAL>(6.0));
   
   if (ax > convert<REAL>(0.0)) ax = SPEED;
   else ax = -SPEED;
   if (ay > convert<REAL>(0.0)) ay = SPEED;
   else ay = -SPEED;
   if (az > convert<REAL>(0.0)) az = SPEED;
   else az = -SPEED;
}

template<typename UINT,typename REAL> void calcAccFaceY(
	REAL& ax, REAL& ay, REAL& az,
	const UINT& I,const UINT& J,const UINT& K,
	const REAL* const cellParams,
	const REAL* const blockParams
) {
   ax = blockParams[BlockParams::VXCRD] + (I+convert<REAL>(0.5))*blockParams[BlockParams::DVX];
   ay = blockParams[BlockParams::VYCRD] + (J+convert<REAL>(0.5))*blockParams[BlockParams::DVY];
   az = blockParams[BlockParams::VZCRD] + (K+convert<REAL>(0.5))*blockParams[BlockParams::DVZ];
   
   const REAL SPEED = convert<REAL>(0.5) + convert<REAL>(0.5)*(convert<REAL>(1.0) / convert<REAL>(6.0));
   
   if (ax > convert<REAL>(0.0)) ax = SPEED;
   else ax = -SPEED;
   if (ay > convert<REAL>(0.0)) ay = SPEED;
   else ay = -SPEED;
   if (az > convert<REAL>(0.0)) az = SPEED;
   else az = -SPEED;
}

template<typename UINT,typename REAL> void calcAccFaceZ(
	REAL& ax, REAL& ay, REAL& az,
	const UINT& I,const UINT& J,const UINT& K,
	const REAL* const cellParams,
	const REAL* const blockParams
) {
   ax = blockParams[BlockParams::VXCRD] + (I+convert<REAL>(0.5))*blockParams[BlockParams::DVX];
   ay = blockParams[BlockParams::VYCRD] + (J+convert<REAL>(0.5))*blockParams[BlockParams::DVY];
   az = blockParams[BlockParams::VZCRD] + (K+convert<REAL>(0.5))*blockParams[BlockParams::DVZ];
   
   const REAL SPEED = convert<REAL>(0.5) + convert<REAL>(0.5)*(convert<REAL>(1.0) / convert<REAL>(6.0));
   
   if (ax > convert<REAL>(0.0)) ax = SPEED;
   else ax = -SPEED;
   if (ay > convert<REAL>(0.0)) ay = SPEED;
   else ay = -SPEED;
   if (az > convert<REAL>(0.0)) az = SPEED;
   else az = -SPEED;
}

#endif

