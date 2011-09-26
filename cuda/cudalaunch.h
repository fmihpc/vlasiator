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

#ifndef CUDA_LAUNCH_H
#define CUDA_LAUNCH_H

#ifdef DEVICE_DEFAULT
// ---------- Default Parameters ----------
cuint VELCOPY_GSIZEX = 8;
cuint VELCOPY_WARPS = 1;
#define VELCOPY_1WARP
#define EXTRA_VELCOPY_1WARP

cuint VEL_DERIVS_GSIZEX = 8;
cuint VEL_DERIVS_WARPS = 1;
#define VEL_DERIVS_1WARP
#define EXTRA_VEL_DERIVS_1WARP

cuint COPY_VEL_DERIVS_GSIZEX = 8;
cuint COPY_VEL_DERIVS_WARPS = 1;
#define COPY_VEL_DERIVS_HWARP
#define EXTRA_COPY_VEL_DERIVS_HWARP

cuint VELFLUX_GSIZEX = 8;
cuint VELFLUX_WARPS = 1;
#define VELFLUX_2WARP
#define EXTRA_VELFLUX_1WARP

cuint COPY_VEL_FLUX_GSIZEX = 8;
cuint COPY_VEL_FLUX_WARPS = 1;
#define COPY_VEL_FLUX_HWARP
#define EXTRA_COPY_VEL_FLUX_HWARP

cuint VELPROP_GSIZEX = 8;
cuint VELPROP_WARPS = 1;
#define VELPROP_1WARP
#define EXTRA_VELPROP_1WARP

cuint SPAT_DERIVS_GSIZEX = 8;
cuint WARPS_SPAT_DERIVS = 1;
#define SPAT_DERIVS_2WARP
#define EXTRA_SPAT_DERIVS_2WARP

cuint FLUXES_GSIZEX = 8;
cuint FLUXES_WARPS = 1;
#define FLUXES_2WARP
#define EXTRA_FLUXES_2WARP

cuint PROPAG_GSIZEX = 8;
cuint PROPAG_WARPS = 1;
#define PROPAG_2WARP
#define EXTRA_PROPAG_1WARP
#endif // Default parameters 

#ifdef DEVICE_TESLA
// ---------- TESLA CARDS ----------
cuint VELCOPY_GSIZEX = 8;
cuint VELCOPY_WARPS = 1;
#define VELCOPY_1WARP
#define EXTRA_VELCOPY_1WARP

cuint VEL_DERIVS_GSIZEX = 8;
cuint VEL_DERIVS_WARPS = 1;
#define VEL_DERIVS_1WARP
#define EXTRA_VEL_DERIVS_1WARP

cuint COPY_VEL_DERIVS_GSIZEX = 8;
cuint COPY_VEL_DERIVS_WARPS = 1;
#define COPY_VEL_DERIVS_HWARP
#define EXTRA_COPY_VEL_DERIVS_HWARP

cuint VELFLUX_GSIZEX = 8;
cuint VELFLUX_WARPS = 2;
#define VELFLUX_N1WARP
#define EXTRA_VELFLUX_1WARP

cuint COPY_VEL_FLUX_GSIZEX = 1;
cuint COPY_VEL_FLUX_WARPS = 8;
#define COPY_VEL_FLUX_NHWARP
#define EXTRA_COPY_VEL_FLUX_HWARP

cuint VELPROP_GSIZEX = 8;
cuint VELPROP_WARPS = 1;
#define VELPROP_1WARP
#define EXTRA_VELPROP_1WARP

cuint SPAT_DERIVS_GSIZEX = 8;
cuint WARPS_SPAT_DERIVS = 1;
#define SPAT_DERIVS_2WARP
#define EXTRA_SPAT_DERIVS_2WARP

cuint FLUXES_GSIZEX = 8;
cuint FLUXES_WARPS = 1;
#define FLUXES_2WARP
#define EXTRA_FLUXES_2WARP

cuint PROPAG_GSIZEX = 8;
cuint PROPAG_WARPS = 1;
#define PROPAG_2WARP
#define EXTRA_PROPAG_1WARP
#endif // Tesla parameters

#endif
