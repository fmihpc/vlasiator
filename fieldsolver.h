#ifndef FIELDSOLVER_H
#define FIELDSOLVER_H

#include "definitions.h"
#include "common.h"
#include "pargrid.h"
#include "cell_spatial.h"

namespace fieldsolver {
   
   enum RecVars {drhodx,drhody,drhodz,
	dBxdy,dBxdz,dBydx,dBydz,dBzdx,dBzdy,
	dVxdx,dVxdy,dVxdz,dVydx,dVydy,dVydz,dVzdx,dVzdy,dVzdz
   };
   
} // namespace fieldsolver

bool finalizeFieldPropagator(ParGrid<SpatialCell>& mpiGrid);
bool initializeFieldPropagator(ParGrid<SpatialCell>& mpiGrid);
bool propagateFields(ParGrid<SpatialCell>& mpiGrid,creal& dt);

#endif
