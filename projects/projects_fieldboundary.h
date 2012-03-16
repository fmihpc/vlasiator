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

#ifndef FIELDBOUNDARY_H
#define FIELDBOUNDARY_H

#include "../arrayallocator.h"
#include "../fieldsolver.h"
#include "projects_common.h"

template<typename CELLID,typename UINT,typename REAL>
REAL fieldBoundaryCopyFromExistingFaceNbrBx(
	const CELLID& cellID,
	const UINT& existingCells,
	const UINT& nonExistingCells,
	#ifdef PARGRID
	const ParGrid<SpatialCell>& mpiGrid
	#else
	const dccrg::Dccrg<SpatialCell>& mpiGrid
	#endif
) { 
   namespace p = projects;
   switch ((nonExistingCells & p::FACE_NBR_BITMASK)) {
    case p::MISSING_ZNEG:
      // If +z neighbour exists, copy value from there:
      if (((existingCells >> p::ZP1_YCC_XCC) & 1) == 1) {
         const CELLID neighbor = getNeighbour(mpiGrid, cellID, 0, 0, +1);
         if (mpiGrid[neighbor] == NULL) {
            std::cerr << __FILE__ << ":" << __LINE__
               << " No data for neighbor " << neighbor
               << " while solving cell " << cellID
               << std::endl;
            abort();
         }
	 return mpiGrid[neighbor]->parameters[CellParams::BX]
	    + mpiGrid[neighbor]->parameters[CellParams::BXFACEX0];
      }
      break;
    case p::MISSING_YNEG:
      // If +y neighbour exists, copy value from there:
      if (((existingCells >> p::ZCC_YP1_XCC) & 1) == 1) {
         const CELLID neighbor = getNeighbour(mpiGrid, cellID, 0, +1, 0);
         if (mpiGrid[neighbor] == NULL) {
            std::cerr << __FILE__ << ":" << __LINE__
               << " No data for neighbor " << neighbor
               << " while solving cell " << cellID
               << std::endl;
            abort();
         }
	 return mpiGrid[neighbor]->parameters[CellParams::BX]
	    + mpiGrid[neighbor]->parameters[CellParams::BXFACEX0];
      }
      break;
    case p::MISSING_YPOS:
      // If -y neighbour exists, copy value from there:
      if (((existingCells >> p::ZCC_YM1_XCC) & 1) == 1) {
         const CELLID neighbor = getNeighbour(mpiGrid, cellID, 0, -1, 0);
         if (mpiGrid[neighbor] == NULL) {
            std::cerr << __FILE__ << ":" << __LINE__
               << " No data for neighbor " << neighbor
               << " while solving cell " << cellID
               << std::endl;
            abort();
         }
	 return mpiGrid[neighbor]->parameters[CellParams::BX]
	    + mpiGrid[neighbor]->parameters[CellParams::BXFACEX0];
      }
      break;
    case p::MISSING_ZPOS:
      // If -z neighbour exists, copy value from there:
      if (((existingCells >> p::ZM1_YCC_XCC) & 1) == 1) {
         const CELLID neighbor = getNeighbour(mpiGrid, cellID, 0, 0, -1);
         if (mpiGrid[neighbor] == NULL) {
            std::cerr << __FILE__ << ":" << __LINE__
               << " No data for neighbor " << neighbor
               << " while solving cell " << cellID
               << std::endl;
            abort();
         }
	 return mpiGrid[neighbor]->parameters[CellParams::BX]
	    + mpiGrid[neighbor]->parameters[CellParams::BXFACEX0];
      }
      break;
    default:
      return 0.0;
   }
   return 0.0;
}

template<typename CELLID,typename UINT,typename REAL>
REAL fieldBoundaryCopyFromExistingFaceNbrBy(
	const CELLID& cellID,
	const UINT& existingCells,
	const UINT& nonExistingCells,
	#ifdef PARGRID
	const ParGrid<SpatialCell>& mpiGrid
	#else
	const dccrg::Dccrg<SpatialCell>& mpiGrid
	#endif
) {
   namespace p = projects;
   switch ((nonExistingCells & p::FACE_NBR_BITMASK)) {
    case p::MISSING_ZNEG:
      // If +z neighbour exists, copy value from there:
      if (((existingCells >> p::ZP1_YCC_XCC) & 1) == 1) {
         const CELLID neighbor = getNeighbour(mpiGrid,cellID,0,0,+1);
         if (mpiGrid[neighbor] == NULL) {
            std::cerr << __FILE__ << ":" << __LINE__
               << " No data for neighbor " << neighbor
               << " while solving cell " << cellID
               << std::endl;
            abort();
         }

         if (mpiGrid[neighbor]->parameters == NULL) {
            std::cerr << __FILE__ << ":" << __LINE__
               << " No cell parameters for cell " << neighbor
               << " while solving cell " << cellID
               << std::endl;
            abort();
         }

	 return mpiGrid[neighbor]->parameters[CellParams::BY]
	    + mpiGrid[neighbor]->parameters[CellParams::BYFACEY0];
      }
      break;
    case p::MISSING_XNEG:
      // If +x neighbour exists, copy value from there:
      if (((existingCells >> p::ZCC_YCC_XP1) & 1) == 1) {
         const CELLID neighbor = getNeighbour(mpiGrid,cellID,+1,0,0);
         if (mpiGrid[neighbor] == NULL) {
            std::cerr << __FILE__ << ":" << __LINE__
               << " No data for neighbor " << neighbor
               << " while solving cell " << cellID
               << std::endl;
            abort();
         }

         if (mpiGrid[neighbor]->parameters == NULL) {
            std::cerr << __FILE__ << ":" << __LINE__
               << " No cell parameters for cell " << neighbor
               << " while solving cell " << cellID
               << std::endl;
            abort();
         }

	 return mpiGrid[neighbor]->parameters[CellParams::BY]
	    + mpiGrid[neighbor]->parameters[CellParams::BYFACEY0];
      }
      break;
    case p::MISSING_XPOS:
      // If -x neighbour exists, copy value from there:
      if (((existingCells >> p::ZCC_YCC_XM1) & 1) == 1) {
         const CELLID neighbor = getNeighbour(mpiGrid,cellID,-1,0,0);
         if (mpiGrid[neighbor] == NULL) {
            std::cerr << __FILE__ << ":" << __LINE__
               << " No data for neighbor " << neighbor
               << " while solving cell " << cellID
               << std::endl;
            abort();
         }

         if (mpiGrid[neighbor]->parameters == NULL) {
            std::cerr << __FILE__ << ":" << __LINE__
               << " No cell parameters for cell " << neighbor
               << " while solving cell " << cellID
               << std::endl;
            abort();
         }

	 return mpiGrid[neighbor]->parameters[CellParams::BY]
	    + mpiGrid[neighbor]->parameters[CellParams::BYFACEY0];
      }
      break;
    case p::MISSING_ZPOS:
      // If -z neighbour exists, copy value from there:
      if (((existingCells >> p::ZM1_YCC_XCC) & 1) == 1) {
         const CELLID neighbor = getNeighbour(mpiGrid,cellID,0,0,-1);
         if (mpiGrid[neighbor] == NULL) {
            std::cerr << __FILE__ << ":" << __LINE__
               << " No data for neighbor " << neighbor
               << " while solving cell " << cellID
               << std::endl;
            abort();
         }

         if (mpiGrid[neighbor]->parameters == NULL) {
            std::cerr << __FILE__ << ":" << __LINE__
               << " No cell parameters for cell " << neighbor
               << " while solving cell " << cellID
               << std::endl;
            abort();
         }

	 return mpiGrid[neighbor]->parameters[CellParams::BY]
	    + mpiGrid[neighbor]->parameters[CellParams::BYFACEY0];
      }
      break;
    default:
      return 0.0;
   }
   return 0.0;
}

template<typename CELLID,typename UINT,typename REAL>
#ifdef PARGRID
REAL fieldBoundaryCopyFromExistingFaceNbrBz(const CELLID& cellID,const UINT& existingCells,const UINT& nonExistingCells,const ParGrid<SpatialCell>& mpiGrid) {
#else
REAL fieldBoundaryCopyFromExistingFaceNbrBz(const CELLID& cellID,const UINT& existingCells,const UINT& nonExistingCells,const dccrg::Dccrg<SpatialCell>& mpiGrid) {
#endif

    namespace p = projects;
   switch ((nonExistingCells & p::FACE_NBR_BITMASK)) {
    case p::MISSING_YNEG:
      // If +y neighbour exists, copy value from there:
      if (((existingCells >> p::ZCC_YP1_XCC) & 1) == 1) {
	 return mpiGrid[getNeighbour(mpiGrid,cellID,0,+1,0)]->parameters[CellParams::BY]
	    + mpiGrid[getNeighbour(mpiGrid,cellID,0,+1,0)]->parameters[CellParams::BYFACEY0];
      }
      break;
    case p::MISSING_XNEG:
      // If +x neighbour exists, copy value from there:
      if (((existingCells >> p::ZCC_YCC_XP1) & 1) == 1) {
	 return mpiGrid[getNeighbour(mpiGrid,cellID,+1,0,0)]->parameters[CellParams::BY]
	    + mpiGrid[getNeighbour(mpiGrid,cellID,+1,0,0)]->parameters[CellParams::BYFACEY0];
      }
      break;
    case p::MISSING_XPOS:
      // If -x neighbour exists, copy value from there:
      if (((existingCells >> p::ZCC_YCC_XM1) & 1) == 1) {
	 return mpiGrid[getNeighbour(mpiGrid,cellID,-1,0,0)]->parameters[CellParams::BY]
	    + mpiGrid[getNeighbour(mpiGrid,cellID,-1,0,0)]->parameters[CellParams::BYFACEY0];
      }
      break;
    case p::MISSING_YPOS:
      // If -y neighbour exists, copy value from there:
      if (((existingCells >> p::ZCC_YM1_XCC) & 1) == 1) {
	 return mpiGrid[getNeighbour(mpiGrid,cellID,0,-1,0)]->parameters[CellParams::BY]
	    + mpiGrid[getNeighbour(mpiGrid,cellID,0,-1,0)]->parameters[CellParams::BYFACEY0];
      }
      break;
    default:
      return 0.0;
   }
   return 0.0;
}

template<typename CELLID,typename UINT,typename REAL>
#ifdef PARGRID
void fieldSolverBoundarySetValueDerivX(const CELLID& cellID,REAL* const array,const UINT& existingCells,const UINT& nonExistingCells,
				       creal* derivatives,const ParGrid<SpatialCell>& mpiGrid,const REAL& value) {
#else
void fieldSolverBoundarySetValueDerivX(const CELLID& cellID,REAL* const array,const UINT& existingCells,const UINT& nonExistingCells,
                                           creal* derivatives,const dccrg::Dccrg<SpatialCell>& mpiGrid,const REAL& value) {
#endif
    namespace fs = fieldsolver;
   array[fs::drhodx] = value;
   array[fs::dBydx]  = value;
   array[fs::dBzdx]  = value;
   array[fs::dVxdx]  = value;
   array[fs::dVydx]  = value;
   array[fs::dVzdx]  = value;
}

template<typename CELLID,typename UINT,typename REAL>
#ifdef PARGRID
void fieldSolverBoundarySetValueDerivY(const CELLID& cellID,REAL* const array,const UINT& existingCells,const UINT& nonExistingCells,
				       creal* derivatives,const ParGrid<SpatialCell>& mpiGrid,const REAL& value) {
#else    
void fieldSolverBoundarySetValueDerivY(const CELLID& cellID,REAL* const array,const UINT& existingCells,const UINT& nonExistingCells,
				       creal* derivatives,const dccrg::Dccrg<SpatialCell>& mpiGrid,const REAL& value) {
#endif
   namespace fs = fieldsolver;
   array[fs::drhody] = value;
   array[fs::dBxdy]  = value;
   array[fs::dBzdy]  = value;
   array[fs::dVxdy]  = value;
   array[fs::dVydy]  = value;
   array[fs::dVzdy]  = value;
}

template<typename CELLID,typename UINT,typename REAL>
#ifdef PARGRID
void fieldSolverBoundarySetValueDerivZ(const CELLID& cellID,REAL* const array,const UINT& existingCells,const UINT& nonExistingCells,
				       creal* derivatives,const ParGrid<SpatialCell>& mpiGrid,const REAL& value) {
#else
void fieldSolverBoundarySetValueDerivZ(const CELLID& cellID,REAL* const array,const UINT& existingCells,const UINT& nonExistingCells,
				       creal* derivatives,const dccrg::Dccrg<SpatialCell>& mpiGrid,const REAL& value) {
#endif
    namespace fs = fieldsolver;
   array[fs::drhodz] = value;
   array[fs::dBxdz]  = value;
   array[fs::dBydz]  = value;
   array[fs::dVxdz]  = value;
   array[fs::dVydz]  = value;
   array[fs::dVzdz]  = value;
}

#endif
