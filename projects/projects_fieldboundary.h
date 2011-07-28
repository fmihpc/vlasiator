#ifndef FIELDBOUNDARY_H
#define FIELDBOUNDARY_H

#include "../arrayallocator.h"
#include "../fieldsolver.h"
#include "projects_common.h"

template<typename CELLID,typename UINT,typename REAL>
REAL fieldBoundaryCopyFromExistingFaceNbrBx(const CELLID& cellID,const UINT& existingCells,const UINT& nonExistingCells,const ParGrid<SpatialCell>& mpiGrid) {
   namespace p = projects;
   switch ((nonExistingCells & p::FACE_NBR_BITMASK)) {
    case p::MISSING_ZNEG:
      // If +z neighbour exists, copy value from there:
      if (((existingCells >> p::ZP1_YCC_XCC) & 1) == 1) {
	 return mpiGrid[getNeighbour(mpiGrid,cellID,0,0,+1)]->cpu_cellParams[CellParams::BX];
      }
      break;
    case p::MISSING_YNEG:
      // If +y neighbour exists, copy value from there:
      if (((existingCells >> p::ZCC_YP1_XCC) & 1) == 1) {
	 return mpiGrid[getNeighbour(mpiGrid,cellID,0,+1,0)]->cpu_cellParams[CellParams::BX];
      }
      break;
    case p::MISSING_YPOS:
      // If -y neighbour exists, copy value from there:
      if (((existingCells >> p::ZCC_YM1_XCC) & 1) == 1) {
	 return mpiGrid[getNeighbour(mpiGrid,cellID,0,-1,0)]->cpu_cellParams[CellParams::BX];
      }
      break;
    case p::MISSING_ZPOS:
      // If -z neighbour exists, copy value from there:
      if (((existingCells >> p::ZM1_YCC_XCC) & 1) == 1) {
	 return mpiGrid[getNeighbour(mpiGrid,cellID,0,0,-1)]->cpu_cellParams[CellParams::BX];
      }
      break;
    default:
      return 0.0;
   }
   return 0.0;
}

template<typename CELLID,typename UINT,typename REAL>
REAL fieldBoundaryCopyFromExistingFaceNbrBy(const CELLID& cellID,const UINT& existingCells,const UINT& nonExistingCells,const ParGrid<SpatialCell>& mpiGrid) {
   namespace p = projects;
   switch ((nonExistingCells & p::FACE_NBR_BITMASK)) {
    case p::MISSING_ZNEG:
      // If +z neighbour exists, copy value from there:
      if (((existingCells >> p::ZP1_YCC_XCC) & 1) == 1) {
	 return mpiGrid[getNeighbour(mpiGrid,cellID,0,0,+1)]->cpu_cellParams[CellParams::BY];
      }
      break;
    case p::MISSING_XNEG:
      // If +x neighbour exists, copy value from there:
      if (((existingCells >> p::ZCC_YCC_XP1) & 1) == 1) {
	 return mpiGrid[getNeighbour(mpiGrid,cellID,+1,0,0)]->cpu_cellParams[CellParams::BY];
      }
      break;
    case p::MISSING_XPOS:
      // If -x neighbour exists, copy value from there:
      if (((existingCells >> p::ZCC_YCC_XM1) & 1) == 1) {
	 return mpiGrid[getNeighbour(mpiGrid,cellID,-1,0,0)]->cpu_cellParams[CellParams::BY];
      }
      break;
    case p::MISSING_ZPOS:
      // If -z neighbour exists, copy value from there:
      if (((existingCells >> p::ZM1_YCC_XCC) & 1) == 1) {
	 return mpiGrid[getNeighbour(mpiGrid,cellID,0,0,-1)]->cpu_cellParams[CellParams::BY];
      }
      break;
    default:
      return 0.0;
   }
   return 0.0;
}

template<typename CELLID,typename UINT,typename REAL>
REAL fieldBoundaryCopyFromExistingFaceNbrBz(const CELLID& cellID,const UINT& existingCells,const UINT& nonExistingCells,const ParGrid<SpatialCell>& mpiGrid) {
   namespace p = projects;
   switch ((nonExistingCells & p::FACE_NBR_BITMASK)) {
    case p::MISSING_YNEG:
      // If +y neighbour exists, copy value from there:
      if (((existingCells >> p::ZCC_YP1_XCC) & 1) == 1) {
	 return mpiGrid[getNeighbour(mpiGrid,cellID,0,+1,0)]->cpu_cellParams[CellParams::BY];
      }
      break;
    case p::MISSING_XNEG:
      // If +x neighbour exists, copy value from there:
      if (((existingCells >> p::ZCC_YCC_XP1) & 1) == 1) {
	 return mpiGrid[getNeighbour(mpiGrid,cellID,+1,0,0)]->cpu_cellParams[CellParams::BY];
      }
      break;
    case p::MISSING_XPOS:
      // If -x neighbour exists, copy value from there:
      if (((existingCells >> p::ZCC_YCC_XM1) & 1) == 1) {
	 return mpiGrid[getNeighbour(mpiGrid,cellID,-1,0,0)]->cpu_cellParams[CellParams::BY];
      }
      break;
    case p::MISSING_YPOS:
      // If -y neighbour exists, copy value from there:
      if (((existingCells >> p::ZCC_YM1_XCC) & 1) == 1) {
	 return mpiGrid[getNeighbour(mpiGrid,cellID,0,-1,0)]->cpu_cellParams[CellParams::BY];
      }
      break;
    default:
      return 0.0;
   }
   return 0.0;
}

template<typename CELLID,typename UINT,typename REAL>
void fieldSolverBoundarySetValueDerivX(const CELLID& cellID,REAL* const array,const UINT& existingCells,const UINT& nonExistingCells,
				       creal* derivatives,const ParGrid<SpatialCell>& mpiGrid,const REAL& value) {
   namespace fs = fieldsolver;
   array[fs::drhodx] = value;
   array[fs::dBydx]  = value;
   array[fs::dBzdx]  = value;
   array[fs::dVxdx]  = value;
   array[fs::dVydx]  = value;
   array[fs::dVzdx]  = value;
}

template<typename CELLID,typename UINT,typename REAL>
void fieldSolverBoundarySetValueDerivY(const CELLID& cellID,REAL* const array,const UINT& existingCells,const UINT& nonExistingCells,
				       creal* derivatives,const ParGrid<SpatialCell>& mpiGrid,const REAL& value) {
   namespace fs = fieldsolver;
   array[fs::drhody] = value;
   array[fs::dBxdy]  = value;
   array[fs::dBzdy]  = value;
   array[fs::dVxdy]  = value;
   array[fs::dVydy]  = value;
   array[fs::dVzdy]  = value;
}

template<typename CELLID,typename UINT,typename REAL>
void fieldSolverBoundarySetValueDerivZ(const CELLID& cellID,REAL* const array,const UINT& existingCells,const UINT& nonExistingCells,
				       creal* derivatives,const ParGrid<SpatialCell>& mpiGrid,const REAL& value) {
   namespace fs = fieldsolver;
   array[fs::drhodz] = value;
   array[fs::dBxdz]  = value;
   array[fs::dBydz]  = value;
   array[fs::dVxdz]  = value;
   array[fs::dVydz]  = value;
   array[fs::dVzdz]  = value;
}

#endif
