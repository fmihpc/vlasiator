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
#ifndef GPU_ACC_MAP_H
#define GPU_ACC_MAP_H

#include "../spatial_cells/spatial_cell_wrapper.hpp"
#include "vec.h"
#include "../common.h"

#include "gpu_acc_sort_blocks.hpp"
#include "../arch/gpu_base.hpp"

bool gpu_acc_map_1d(spatial_cell::SpatialCell* spatial_cell,
                     const uint popID,
                     Real intersection,
                     Real intersection_di,
                     Real intersection_dj,
                     Real intersection_dk,
                     const uint dimension
   );

#endif
