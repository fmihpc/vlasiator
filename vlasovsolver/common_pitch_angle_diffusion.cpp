#pragma once

/*
 * This file is part of Vlasiator.
 * Copyright 2010-2025 Finnish Meteorological Institute and University of Helsinki
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

#include "../common.h"
#include "../spatial_cells/spatial_cell_wrapper.hpp"
#include <dccrg.hpp>
#include <dccrg_cartesian_geometry.hpp>
#include <zoltan.h>

#ifdef USE_GPU
#include "gpu_pitch_angle_diffusion.hpp"
#else
#include "cpu_pitch_angle_diffusion.h"
#endif

extern std::vector<Real> betaParaArray;
extern std::vector<Real> TanisoArray;
extern std::vector<Real> nu0Array;
extern size_t n_betaPara;
extern size_t n_Taniso;
extern bool nuArrayRead;

void readNuArrayFromFile();

Realf interpolateNuFromArray(const Real Taniso, const Real betaParallel);

void computePitchAngleDiffusionParameters(SpatialCell& cell, const uint popID, size_t CellIdx, bool& currentSpatialLoopComplete, Realf& sparsity, std::array<Real, 3>& b, Real& nu0);
