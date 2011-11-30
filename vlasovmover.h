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

#ifndef VLASOVMOVER_H
#define VLASOVMOVER_H

#include <vector>


#include "definitions.h"
#include "spatial_cell.hpp"
using namespace spatial_cell;

bool finalizeMover();

#include <stdint.h>

#define DCCRG_SEND_SINGLE_CELLS
#define DCCRG_CELL_DATA_SIZE_FROM_USER
#define DCCRG_USER_MPI_DATA_TYPE
#include <dccrg.hpp>

bool initializeMover(dccrg::Dccrg<SpatialCell>& mpiGrid);
void calculateSimParameters(dccrg::Dccrg<SpatialCell>& mpiGrid,creal& t,Real& dt);
void calculateCellParameters(dccrg::Dccrg<SpatialCell>& mpiGrid,creal& t,uint64_t& cell);
void calculateAcceleration(dccrg::Dccrg<SpatialCell>& mpiGrid);
void calculateSpatialDerivatives(dccrg::Dccrg<SpatialCell>& mpiGrid);
void calculateSpatialFluxes(dccrg::Dccrg<SpatialCell>& mpiGrid);
void calculateSpatialPropagation(dccrg::Dccrg<SpatialCell>& mpiGrid,const bool& secondStep,const bool& transferAvgs);
void initialLoadBalance(dccrg::Dccrg<SpatialCell>& mpiGrid);
void calculateVelocityMoments(dccrg::Dccrg<SpatialCell>& mpiGrid);

#endif


