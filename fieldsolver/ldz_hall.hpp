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

#ifndef LDZ_HALL_HPP
#define LDZ_HALL_HPP

#include "../definitions.h"

void calculateHallTermSimple(std::span<std::array<Real, fsgrids::bfield::N_BFIELD>> perb,
                             std::span<std::array<Real, fsgrids::bfield::N_BFIELD>> perbdt2,
                             std::span<std::array<Real, fsgrids::ehall::N_EHALL>> ehall,
                             std::span<std::array<Real, fsgrids::moments::N_MOMENTS>> moments,
                             std::span<std::array<Real, fsgrids::moments::N_MOMENTS>> momentsdt2,
                             std::span<std::array<Real, fsgrids::dperb::N_DPERB>> dperb,
                             std::span<std::array<Real, fsgrids::dmoments::N_DMOMENTS>> dmoments,
                             std::span<std::array<Real, fsgrids::dmoments::N_DMOMENTS>> dmomentsdt2,
                             std::span<std::array<Real, fsgrids::bgbfield::N_BGB>> bgb,
                             fsgrid::FsGrid<fsgrids::technical, FS_STENCIL_WIDTH>& technicalGrid,
                             SysBoundary& sysBoundaries, int32_t RKCase, const bool communicateMomentsDerivatives);

#endif
