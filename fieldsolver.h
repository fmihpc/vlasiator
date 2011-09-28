/*
This file is part of Vlasiator.

Copyright 2010, 2011 Finnish Meteorological Institute

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
bool initializeFieldPropagator(ParGrid<SpatialCell>& mpiGrid,bool propagateFields);
bool propagateFields(ParGrid<SpatialCell>& mpiGrid,creal& dt);

#else

void calculateFaceAveragedFields(dccrg<SpatialCell>& mpiGrid);
void calculateVolumeAveragedFields(dccrg<SpatialCell>& mpiGrid);
bool finalizeFieldPropagator(dccrg<SpatialCell>& mpiGrid);
bool initializeFieldPropagator(dccrg<SpatialCell>& mpiGrid,bool propagateFields);
bool propagateFields(dccrg<SpatialCell>& mpiGrid,creal& dt);

#endif

#endif
