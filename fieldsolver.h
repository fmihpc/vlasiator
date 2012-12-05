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
#include "sysboundary/sysboundary.h"
using namespace spatial_cell;

#include <dccrg.hpp>


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
bool initializeFieldPropagator(dccrg::Dccrg<SpatialCell>& mpiGrid,
                               SysBoundary& sysBoundaries);
bool initializeFieldPropagatorAfterRebalance(dccrg::Dccrg<SpatialCell>& mpiGrid);
bool propagateFields(dccrg::Dccrg<SpatialCell>& mpiGrid,
                     SysBoundary& sysBoundaries,
                     creal& dt);
void calculateEdgeElectricFieldX(
   dccrg::Dccrg<SpatialCell>& mpiGrid,
   const CellID& cellID,
   cint& RKCase);
void calculateEdgeElectricFieldY(
   dccrg::Dccrg<SpatialCell>& mpiGrid,
   const CellID& cellID,
   cint& RKCase);
void calculateEdgeElectricFieldZ(
   dccrg::Dccrg<SpatialCell>& mpiGrid,
   const CellID& cellID,
   cint& RKCase);

CellID getNeighbourID(
   dccrg::Dccrg<SpatialCell>& mpiGrid,
   const CellID& cellID,
   const uchar& i,
   const uchar& j,
   const uchar& k
);
Real limiter(creal& left,creal& cent,creal& rght);
Real divideIfNonZero(creal rhoV, creal rho);
#endif
