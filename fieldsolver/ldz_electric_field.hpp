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

#ifndef LDZ_ELECTRIC_FIELD_HPP
#define LDZ_ELECTRIC_FIELD_HPP

#include "fs_common.h"

void calculateUpwindedElectricFieldSimple(fsgrids::perbspan perb,
                                          fsgrids::perbspan perbdt2,
                                          fsgrids::efieldspan e,
                                          fsgrids::efieldspan edt2,
                                          fsgrids::ehallspan ehall,
                                          fsgrids::egradpespan egradpe,
                                          fsgrids::egradpespan egradpedt2,
                                          fsgrids::momentsspan moments,
                                          fsgrids::momentsspan momentsdt2,
                                          fsgrids::dperbspan dperb,
                                          fsgrids::dmomentsspan dmoments,
                                          fsgrids::dmomentsspan dmomentsdt2,
                                          fsgrids::bgbspan bgb,
                                          fsgrids::technicalspan technical, FieldSolverGrid &fsgrid,
                                          SysBoundary& sysBoundaries, int32_t RKCase,
                                          const bool communicateEGradPeOrMomentsDerivatives);

#endif
