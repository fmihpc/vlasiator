/*
  This file is part of Vlasiator.
  Copyright 2013,2014 Finnish Meteorological Institute
*/

#ifndef CPU_ACC_SEMILAG_H
#define CPU_ACC_SEMILAG_H

#include "../common.h"
#include "../spatial_cell.hpp"

void prepareAccelerateCell(spatial_cell::SpatialCell* spatial_cell, const int popID);
int getAccelerationSubcycles(spatial_cell::SpatialCell* spatial_cell, Real dt, const int& popID);



void cpu_accelerate_cell(
        spatial_cell::SpatialCell* spatial_cell,
        const int popID,
        uint map_order,
        const Real& dt);

#endif

