/*
This file is part of Vlasiator.

Copyright 2010, 2011, 2012 Finnish Meteorological Institute

Vlasiator is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3
as published by the Free Software Foundation.

Vlasiator is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef FIELDSOLVER_H
#define FIELDSOLVER_H

#include "definitions.h"
#include "common.h"
#include "spatial_cell.hpp"         
using namespace spatial_cell;

#include <dccrg.hpp>

#define FIELD_SOLVER_NEIGHBORHOOD_ID 1

/*
namespace fieldsolver {
   
   enum RecVars {drhodx,drhody,drhodz,
	dBxdy,dBxdz,dBydx,dBydz,dBzdx,dBzdy,
	dVxdx,dVxdy,dVxdz,dVydx,dVydy,dVydz,dVzdx,dVzdy,dVzdz
   };
   
} // namespace fieldsolver
*/

void calculateVolumeAveragedFields(dccrg::Dccrg<SpatialCell>& mpiGrid);
bool finalizeFieldPropagator(dccrg::Dccrg<SpatialCell>& mpiGrid);
bool initializeFieldPropagator(dccrg::Dccrg<SpatialCell>& mpiGrid);
bool initializeFieldPropagatorAfterRebalance(dccrg::Dccrg<SpatialCell>& mpiGrid);
bool propagateFields(dccrg::Dccrg<SpatialCell>& mpiGrid,creal& dt);
#endif
