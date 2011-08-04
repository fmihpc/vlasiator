#ifndef FIELDSOLVER_H
#define FIELDSOLVER_H

#include "definitions.h"
#include "common.h"
#include "cell_spatial.h"

#ifdef PARGRID
   #include "pargrid.h"
#else
   #define DCCRG_SEND_SINGLE_CELLS
   #define DCCRG_CELL_DATA_SIZE_FROM_USER
   #define DCCRG_USER_MPI_DATA_TYPE
   #include <dccrg.hpp>
#endif
/*
namespace fieldsolver {
   
   enum RecVars {drhodx,drhody,drhodz,
	dBxdy,dBxdz,dBydx,dBydz,dBzdx,dBzdy,
	dVxdx,dVxdy,dVxdz,dVydx,dVydy,dVydz,dVzdx,dVzdy,dVzdz
   };
   
} // namespace fieldsolver
*/
#ifdef PARGRID

void calculateFaceAveragedFields(ParGrid<SpatialCell>& mpiGrid);
void calculateVolumeAveragedFields(ParGrid<SpatialCell>& mpiGrid);
bool finalizeFieldPropagator(ParGrid<SpatialCell>& mpiGrid);
bool initializeFieldPropagator(ParGrid<SpatialCell>& mpiGrid);
bool propagateFields(ParGrid<SpatialCell>& mpiGrid,creal& dt);

#else

void calculateFaceAveragedFields(dccrg<SpatialCell>& mpiGrid);
void calculateVolumeAveragedFields(dccrg<SpatialCell>& mpiGrid);
bool finalizeFieldPropagator(dccrg<SpatialCell>& mpiGrid);
bool initializeFieldPropagator(dccrg<SpatialCell>& mpiGrid);
bool propagateFields(dccrg<SpatialCell>& mpiGrid,creal& dt);

#endif

#endif
