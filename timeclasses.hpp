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
#ifndef TIMECLASSES_HPP
#define TIMECLASSES_HPP

#include "spatial_cells/spatial_cell_wrapper.hpp"
#include "spatial_cells/block_adjust_wrapper.hpp"
#include <dccrg.hpp>
#include <dccrg_cartesian_geometry.hpp>
#include "vlasovsolver/vlasovmover.h"

bool isDtTooLarge(Real dt, Real rdt, Real vdt, Real fsdt);

bool isDtTooSmall(Real dt, Real rdt, Real vdt, Real fsdt);

std::vector<CellID> checkCellTimeclasses(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid);

void updateTimeclassDts(Real fsdt, const bool applyModifier = true);

void increaseTimeclass(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                              const std::vector<CellID>& cellsToIncreaseTimeclass,
                              bool& additionalTimeclassCreated);

void calculateGlobalTcVariables(Real fsdt, Real globalMaxDt);

void initiateAllCellTimeclasses(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid);

void timeclassDebugAssertions(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid);

#endif