/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 * Copyright 2025 CSC - IT Center for Science
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef LDZ_MAGNETIC_FIELD_HPP
#define LDZ_MAGNETIC_FIELD_HPP

#include <vector>

#include "../definitions.h"
#include "../common.h"
//#include "../spatial_cells/spatial_cell_wrapper.hpp"

#include "fs_common.h"

void propagateMagneticFieldSimple(fsgrids::perbspan perb,
                                  fsgrids::perbspan perbdt2,
                                  fsgrids::bgbspan bgb,
                                  fsgrids::efieldspan e,
                                  fsgrids::efieldspan edt2,
                                  fsgrids::technicalspan technical, FieldSolverGrid &fsgrid,
                                  SysBoundary& sysBoundaries, creal& dt, cint& RKCase);

#endif
