/*
 * This file is part of Vlasiator.
 * Copyright 2010-2024 Finnish Meteorological Institute, University of Helsinki
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
/*!
  Spatial cell wrapper, maps to GPU or CPU version
*/
#ifndef BLOCK_ADJUST_WRAPPER_H
#define BLOCK_ADJUST_WRAPPER_H

#ifdef USE_GPU
#include "block_adjust_gpu.hpp"
#else
#include "block_adjust_cpu.hpp"
#endif

#endif
