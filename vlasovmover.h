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

#ifndef VLASOVMOVER_H
#define VLASOVMOVER_H

#include <vector>


#include "definitions.h"
#include "spatial_cell.hpp"
using namespace spatial_cell;

bool finalizeMover();

#include <stdint.h>


#include <dccrg.hpp>

bool initializeMover(dccrg::Dccrg<SpatialCell>& mpiGrid);
bool initMoverAfterBlockChange(dccrg::Dccrg<SpatialCell>& mpiGrid);
void calculateCellParameters(dccrg::Dccrg<SpatialCell>& mpiGrid,creal& t,uint64_t& cell);

void calculateAcceleration(dccrg::Dccrg<SpatialCell>& mpiGrid,Real dt);
void calculateCellAcceleration(dccrg::Dccrg<SpatialCell>& mpiGrid,uint64_t cellID,Real dt); //calculate acceleration in one single cell
void calculateSpatialFluxes(dccrg::Dccrg<SpatialCell>& mpiGrid,Real dt);
void calculateSpatialPropagation(dccrg::Dccrg<SpatialCell>& mpiGrid,const bool& acceleration,Real acceleration_dt);
void initialLoadBalance(dccrg::Dccrg<SpatialCell>& mpiGrid);
void calculateVelocityMoments(dccrg::Dccrg<SpatialCell>& mpiGrid);

#endif


