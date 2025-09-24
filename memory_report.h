/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://www.physics.helsinki.fi/vlasiator/
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */
#ifndef MEMORY_REPORT_H
#define MEMORY_REPORT_H

#include "definitions.h"
#include "spatial_cells/spatial_cell_wrapper.hpp"
#include <dccrg.hpp>
#include <dccrg_cartesian_geometry.hpp>

/*! Report spatial cell counts per refinement level as well as velocity cell counts per population into logfile
 */
void report_cell_and_block_counts(dccrg::Dccrg<spatial_cell::SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid);

/*! Measures memory consumption and writes it into logfile. Collective
 *  operation on MPI_COMM_WORLD
 *  extra_bytes is used for additional buffer for the high water mark,
 *  for example when estimating refinement memory usage
 */
void report_memory_consumption(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid, double extra_bytes = 0.0);

#endif
