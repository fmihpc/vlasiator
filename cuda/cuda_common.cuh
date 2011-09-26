/*
This file is part of Vlasiator.

Copyright 2010, 2011 Finnish Meteorological Institute

Vlasiator is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3
as published by the Free Software Foundation.

Vlasiator is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef CUDA_COMMON_CUH
#define CUDA_COMMON_CUH

#include <cstdlib>
#include <iostream>

#include <cuda_runtime.h>
#include <device_launch_parameters.h>

#include "logger.h"
#include "common.h"
#include "devicegrid.h"
#include "grid.h"
#include "parameters.h"

using namespace std;

extern Logger logger;

texture<real,1,cudaReadModeElementType> texRef_avgs;
texture<real,1,cudaReadModeElementType> texRef_cellParams;




__device__ uint bindex2(uint bix,uint biy);
__device__ uint tindex2(uint tix,uint tiy);
__device__ uint tindex3(uint tix,uint tiy,uint tiz);

__device__ void loadVelNbrs(uint MYBLOCK,uint* nbrs,uint* sha_nbr);

__device__ void transpose_yz_1warp(real* sha_arr);
__device__ void transpose_xz_1warp(real* sha_arr);

__device__ real reconstruct_neg(real avg,real d1,real d2);
__device__ real reconstruct_pos(real avg,real d1,real d2);

__device__ real minmod(real a,real b);
__device__ real MClimiter(real xl1,real xcc,real xr1);
__device__ real superbee(real xl1,real xcc,real xr1);
__device__ real vanLeer(real xl1,real xcc,real xr1);

__device__ void kernel_stencils(uint MYIND,real* sha_avg,real* d1,real* d2);
__global__ void derivatives_1warp(uint OFFSET,real* avgs,real* avgnbrx,real* avgnbry,real* avgnbrz,
				  real* d1x,real* d1y,real* d1z,real* d2x,real* d2y,real* d2z);

#endif
