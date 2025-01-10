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

#ifndef LDZ_GRADPE_HPP
#define LDZ_GRADPE_HPP

#include "../definitions.h"

void calculateGradPeTermSimple(fsgrid::FsData<std::array<Real, fsgrids::egradpe::N_EGRADPE>>& egradpe,
                               fsgrid::FsData<std::array<Real, fsgrids::egradpe::N_EGRADPE>>& egradpedt2,
                               fsgrid::FsData<std::array<Real, fsgrids::moments::N_MOMENTS>>& moments,
                               fsgrid::FsData<std::array<Real, fsgrids::moments::N_MOMENTS>>& momentsdt2,
                               fsgrid::FsData<std::array<Real, fsgrids::dmoments::N_DMOMENTS>>& dmoments,
                               fsgrid::FsData<std::array<Real, fsgrids::dmoments::N_DMOMENTS>>& dmomentsdt2,
                               fsgrid::FsGrid<fsgrids::technical, FS_STENCIL_WIDTH>& technicalGrid,
                               SysBoundary& sysBoundaries, cint& RKCase);

#endif
