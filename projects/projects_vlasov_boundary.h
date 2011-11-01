#ifndef PROJECTS_VLASOVBOUNDARY_H
#define PROJECTS_VLASOVBOUNDARY_H

#include "projects_common.h"

// *********************************
// ***** TEMPLATE DECLARATIONS *****
// *********************************

#ifdef PARGRID
template<typename CELLID,typename UINT>
void vlasovBoundaryCopyFromExistingFaceNbr(const CELLID& cellID,const UINT& boundaryFlags,const ParGrid<SpatialCell>& mpiGrid);
template<typename CELLID,typename UINT,typename REAL> 
void vlasovBoundarySetValue(const CELLID& cellID,const UINT& boundaryFlags,const ParGrid<SpatialCell>& mpiGrid,const REAL& value);
#else
template<typename CELLID,typename UINT>
void vlasovBoundaryCopyFromExistingFaceNbr(const CELLID& cellID,const UINT& boundaryFlags,const dccrg::Dccrg<SpatialCell>& mpiGrid);
template<typename CELLID,typename UINT,typename REAL> 
void vlasovBoundarySetValue(const CELLID& cellID,const UINT& boundaryFlags,const dccrg::Dccrg<SpatialCell>& mpiGrid,const REAL& value);
#endif

// ********************************
// ***** TEMPLATE DEFINITIONS *****
// ********************************
#ifdef PARGRID
template<typename CELLID,typename UINT>
void vlasovBoundaryCopyFromExistingFaceNbr(const CELLID& cellID,const UINT& existingCells,
                                           const UINT& nonExistingCells,const ParGrid<SpatialCell>& mpiGrid) {
#else
template<typename CELLID,typename UINT>
    void vlasovBoundaryCopyFromExistingFaceNbr(const CELLID& cellID,const UINT& existingCells,
                                               const UINT& nonExistingCells,const dccrg::Dccrg<SpatialCell>& mpiGrid) {
#endif
    const UINT missingPosX = (1 << 14);
   const UINT missingNegX = (1 << 12);
   const UINT missingPosY = (1 << 16);
   const UINT missingNegY = (1 << 10);
   const UINT missingPosZ = (1 << 22);
   const UINT missingNegZ = (1 << 4);
   const UINT mask = ((1 << 4) | (1 << 10) | (1 << 12) | (1 << 14) | (1 << 16) | (1 << 22));
   
   switch ((nonExistingCells & mask)) {
    case missingNegZ:
      // If +z neighbour exists, copy values from there:
      if (((existingCells >> 22) & 1) == 1) {
	 for (size_t i=0; i<mpiGrid[cellID]->N_blocks*SIZE_VELBLOCK; ++i) 
	   mpiGrid[cellID]->cpu_avgs[i] = mpiGrid[getNeighbour(mpiGrid,cellID,0,0,+1)]->cpu_avgs[i];
      } else {
	 for (size_t i=0; i<mpiGrid[cellID]->N_blocks*SIZE_VELBLOCK; ++i)
	   mpiGrid[cellID]->cpu_avgs[i] = 0.0;
      }
      break;
    case missingNegY:
      // If +y neighbour exists, copy values from there:
      if (((existingCells >> 16) & 1) == 1) {
	 for (size_t i=0; i<mpiGrid[cellID]->N_blocks*SIZE_VELBLOCK; ++i)
	   mpiGrid[cellID]->cpu_avgs[i] = mpiGrid[getNeighbour(mpiGrid,cellID,0,+1,0)]->cpu_avgs[i];
      } else {
	 for (size_t i=0; i<mpiGrid[cellID]->N_blocks*SIZE_VELBLOCK; ++i)
	   mpiGrid[cellID]->cpu_avgs[i] = 0.0;
      }
      break;
    case missingNegX:
      // If +x neighbour exists, copy values from there:
      if (((existingCells >> 14) & 1) == 1) {
	 for (size_t i=0; i<mpiGrid[cellID]->N_blocks*SIZE_VELBLOCK; ++i) 
	   mpiGrid[cellID]->cpu_avgs[i] = mpiGrid[getNeighbour(mpiGrid,cellID,+1,0,0)]->cpu_avgs[i];
      } else {
	 for (size_t i=0; i<mpiGrid[cellID]->N_blocks*SIZE_VELBLOCK; ++i)
	   mpiGrid[cellID]->cpu_avgs[i] = 0.0;
      }
      break;
    case missingPosX:
      // If -x neighbour exists, copy values from there:
      if (((existingCells >> 12) & 1) == 1) {
	 for (size_t i=0; i<mpiGrid[cellID]->N_blocks*SIZE_VELBLOCK; ++i) 
	   mpiGrid[cellID]->cpu_avgs[i] = mpiGrid[getNeighbour(mpiGrid,cellID,-1,0,0)]->cpu_avgs[i];
      } else {
	 for (size_t i=0; i<mpiGrid[cellID]->N_blocks*SIZE_VELBLOCK; ++i)
	   mpiGrid[cellID]->cpu_avgs[i] = 0.0;
      }
      break;
    case missingPosY:
      // If -y neighbour exists, copy values from there:
      if (((existingCells >> 10) & 1) == 1) {
	 for (size_t i=0; i<mpiGrid[cellID]->N_blocks*SIZE_VELBLOCK; ++i) 
	   mpiGrid[cellID]->cpu_avgs[i] = mpiGrid[getNeighbour(mpiGrid,cellID,0,-1,0)]->cpu_avgs[i];
      } else {
	 for (size_t i=0; i<mpiGrid[cellID]->N_blocks*SIZE_VELBLOCK; ++i)
	   mpiGrid[cellID]->cpu_avgs[i] = 0.0;
      }
      break;
    case missingPosZ:
      // If -z neighbour exists, copy values from there:
      if (((existingCells >>  4) & 1) == 1) {
	 for (size_t i=0; i<mpiGrid[cellID]->N_blocks*SIZE_VELBLOCK; ++i) 
	   mpiGrid[cellID]->cpu_avgs[i] = mpiGrid[getNeighbour(mpiGrid,cellID,0,0,-1)]->cpu_avgs[i];
      } else {
	 for (size_t i=0; i<mpiGrid[cellID]->N_blocks*SIZE_VELBLOCK; ++i)
	   mpiGrid[cellID]->cpu_avgs[i] = 0.0;
      }
      break;
    default:
      // Set boundary cells to zero value:
      for (size_t i=0; i<mpiGrid[cellID]->N_blocks*SIZE_VELBLOCK; ++i) mpiGrid[cellID]->cpu_avgs[i] = 0.0;
      break;
   }
}

#ifdef PARGRID 
template<typename CELLID,typename UINT,typename REAL> 
void vlasovBoundarySetValue(const CELLID& cellID,const UINT& existingCells,
                            const UINT& nonExistingCells,const ParGrid<SpatialCell>& mpiGrid,const REAL& value) {
#else
template<typename CELLID,typename UINT,typename REAL> 
void vlasovBoundarySetValue(const CELLID& cellID,const UINT& existingCells,
                            const UINT& nonExistingCells,const dccrg::Dccrg<SpatialCell>& mpiGrid,const REAL& value) {
#endif
    for (size_t i=0; i<mpiGrid[cellID]->N_blocks*SIZE_VELBLOCK; ++i) mpiGrid[cellID]->cpu_avgs[i] = value;
}





#endif
