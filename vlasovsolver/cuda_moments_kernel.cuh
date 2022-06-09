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
#ifndef CUDA_MOMENTS_KERNELH
#define CUDA_MOMENTS_H

#include "../common.h"
#include "device_launch_parameters.h"
#include "cuda.h"
#include "cuda_runtime.h"

#define nMoments 8;

extern void cuda_allocateMomentCalculations();

void calculate_moments_glue(
   Realf* dev_meshDataPointers,
   Real* dev_parameterPointers,
   uint* dev_blockCounts,
   Real* dev_masses,
   Real* dev_charges,
   Real* dev_momentArrays,
   const int nPopulations,
   const bool computeSecond,
   cudaStream_t stream
   );

extern Realf *dev_meshDataPointers[];
extern Real *dev_parameterPointers[];
extern Real *dev_masses[];
extern Real *dev_charges[];
extern uint *dev_blockCounts[];
extern Real *dev_momentArrays[];

extern std::vector<Realf*> *meshDataPointers[];
extern std::vector<Real*> *parameterPointers[];
extern std::vector<Real> *masses[];
extern std::vector<Real> *charges[];
extern std::vector<uint> *blockCounts[];
extern std::vector<std::array<Real,nMoments> > *momentArrays[];

extern bool isCudaMomentsAllocated;

#endif
