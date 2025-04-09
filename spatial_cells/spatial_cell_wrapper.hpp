/*
 * This file is part of Vlasiator.
 * Copyright 2010-2023 Finnish Meteorological Institute, University of Helsinki
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
#ifndef SPATIAL_CELL_WRAPPER_H
#define SPATIAL_CELL_WRAPPER_H

#ifdef USE_GPU
#include "spatial_cell_gpu.hpp"
#else
#include "spatial_cell_cpu.hpp"
#endif

typedef Parameters P; // Heeded in numerous files which include this one

/*!
Used as an error from functions returning velocity cells or
as a cell that would be outside of the velocity block
*/
#define error_velocity_cell 0xFFFFFFFFu

/*!
Used as an error from functions returning velocity cell indices or
as an index that would be outside of the velocity block
*/
#define error_velocity_cell_index 0xFFFFFFFFu

namespace spatial_cell {

   namespace Transfer {
      static const uint64_t NONE                     = 0;
      static const uint64_t CELL_PARAMETERS          = (1ull<<0);
      static const uint64_t CELL_DERIVATIVES         = (1ull<<1);
      static const uint64_t VEL_BLOCK_LIST_STAGE1    = (1ull<<2);
      static const uint64_t VEL_BLOCK_LIST_STAGE2    = (1ull<<3);
      static const uint64_t VEL_BLOCK_DATA           = (1ull<<4);
      static const uint64_t VEL_BLOCK_PARAMETERS     = (1ull<<6);
      static const uint64_t VEL_BLOCK_WITH_CONTENT_STAGE1  = (1ull<<7);
      static const uint64_t VEL_BLOCK_WITH_CONTENT_STAGE2  = (1ull<<8);
      static const uint64_t CELL_SYSBOUNDARYFLAG     = (1ull<<9);
      static const uint64_t CELL_E                   = (1ull<<10);
      static const uint64_t CELL_EDT2                = (1ull<<11);
      static const uint64_t CELL_PERB                = (1ull<<12);
      static const uint64_t CELL_PERBDT2             = (1ull<<13);
      static const uint64_t CELL_RHOM_V              = (1ull<<14);
      static const uint64_t CELL_RHOMDT2_VDT2        = (1ull<<15);
      static const uint64_t CELL_RHOQ                = (1ull<<16);
      static const uint64_t CELL_RHOQDT2             = (1ull<<17);
      static const uint64_t CELL_BVOL                = (1ull<<18);
      static const uint64_t CELL_BVOL_DERIVATIVES    = (1ull<<19);
      static const uint64_t CELL_DIMENSIONS          = (1ull<<20);
      static const uint64_t CELL_IOLOCALCELLID       = (1ull<<21);
      static const uint64_t NEIGHBOR_VEL_BLOCK_DATA  = (1ull<<22);
      static const uint64_t CELL_HALL_TERM           = (1ull<<23);
      static const uint64_t CELL_P                   = (1ull<<24);
      static const uint64_t CELL_PDT2                = (1ull<<25);
      static const uint64_t POP_METADATA             = (1ull<<26);
      static const uint64_t RANDOMGEN                = (1ull<<27);
      static const uint64_t CELL_GRADPE_TERM         = (1ull<<28);
      static const uint64_t REFINEMENT_PARAMETERS    = (1ull<<29);
      //all data
      static const uint64_t ALL_DATA =
      CELL_PARAMETERS
      | CELL_DERIVATIVES | CELL_BVOL_DERIVATIVES
      | VEL_BLOCK_DATA
      | CELL_SYSBOUNDARYFLAG
      | POP_METADATA | RANDOMGEN;

      //all data, except the distribution function
      static const uint64_t ALL_SPATIAL_DATA =
      CELL_PARAMETERS
      | CELL_DERIVATIVES | CELL_BVOL_DERIVATIVES
      | CELL_SYSBOUNDARYFLAG
      | POP_METADATA | RANDOMGEN;
   }

}

#endif
