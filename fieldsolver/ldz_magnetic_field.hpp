/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
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

#ifndef LDZ_MAGNETIC_FIELD_HPP
#define LDZ_MAGNETIC_FIELD_HPP

#include <vector>

#include "../definitions.h"
#include "../common.h"
//#include "../spatial_cells/spatial_cell_wrapper.hpp"

#include "fs_common.h"

void propagateMagneticFieldSimple(std::span<std::array<Real, fsgrids::bfield::N_BFIELD>> perb,
                                  std::span<std::array<Real, fsgrids::bfield::N_BFIELD>> perbdt2,
                                  std::span<std::array<Real, fsgrids::bgbfield::N_BGB>> bgb,
                                  std::span<std::array<Real, fsgrids::efield::N_EFIELD>> e,
                                  std::span<std::array<Real, fsgrids::efield::N_EFIELD>> edt2,
                                  std::span<fsgrids::technical> technical, FieldSolverGrid &fsgrid,
                                  SysBoundary& sysBoundaries, creal& dt, cint& RKCase);

#endif
