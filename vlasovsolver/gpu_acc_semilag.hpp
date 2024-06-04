/*
 * This file is part of Vlasiator.
 * Copyright 2010-2024 Finnish Meteorological Institute and University of Helsinki
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://www.physics.helsinki.fi/vlasiator/
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#ifndef GPU_ACC_SEMILAG_H
#define GPU_ACC_SEMILAG_H

#include "../common.h"
#include "../spatial_cell_wrapper.hpp"

void prepareAccelerateCell(spatial_cell::SpatialCell* spatial_cell, const uint popID);
uint getAccelerationSubcycles(spatial_cell::SpatialCell* spatial_cell, Real dt, const uint popID);

void gpu_accelerate_cell(
        spatial_cell::SpatialCell* spatial_cell,
        const uint popID,
        const uint map_order,
        const Real& dt);

#endif

