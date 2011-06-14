#ifndef FIELDSOLVER_H
#define FIELDSOLVER_H

#include "definitions.h"
#include "common.h"
#include "cell_spatial.h"

#ifdef PARGRID
   #include "pargrid.h"
#else
   #include <dccrg.hpp>
#endif

namespace fieldsolver {
   
   enum RecVars {drhodx,drhody,drhodz,
	dBxdy,dBxdz,dBydx,dBydz,dBzdx,dBzdy,
	dVxdx,dVxdy,dVxdz,dVydx,dVydy,dVydz,dVzdx,dVzdy,dVzdz
   };
   
} // namespace fieldsolver

#ifdef PARGRID

bool finalizeFieldPropagator(ParGrid<SpatialCell>& mpiGrid);
bool initializeFieldPropagator(ParGrid<SpatialCell>& mpiGrid);
bool propagateFields(ParGrid<SpatialCell>& mpiGrid,creal& dt);

#else

bool finalizeFieldPropagator(dccrg<SpatialCell>& mpiGrid);
bool initializeFieldPropagator(dccrg<SpatialCell>& mpiGrid);
bool propagateFields(dccrg<SpatialCell>& mpiGrid,creal& dt);

#endif

#endif
