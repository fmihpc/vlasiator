/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
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

#ifndef LDZ_GRADPE_HPP
#define LDZ_GRADPE_HPP

#include "../definitions.h"

void calculateGradPeTermSimple(fsgrids::egradpespan egradpe,
                               fsgrids::egradpespan egradpedt2,
                               fsgrids::momentsspan moments,
                               fsgrids::momentsspan momentsdt2,
                               fsgrids::dmomentsspan dmoments,
                               fsgrids::dmomentsspan dmomentsdt2,
                               fsgrids::technicalspan technical,
                               FieldSolverGrid &fsgrid,
                               SysBoundary& sysBoundaries,
                               cint& RKCase);

#endif
