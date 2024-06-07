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

#ifndef PROJECTS_COMMON_H
#define PROJECTS_COMMON_H
#include "../spatial_cell_wrapper.hpp"
#include <dccrg.hpp>
#include <dccrg_cartesian_geometry.hpp>

namespace projects {
   enum Neighbours {
      ZM1_YM1_XM1,
      ZM1_YM1_XCC,
      ZM1_YM1_XP1,
      ZM1_YCC_XM1,
      ZM1_YCC_XCC,
      ZM1_YCC_XP1,
      ZM1_YP1_XM1,
      ZM1_YP1_XCC,
      ZM1_YP1_XP1,
      ZCC_YM1_XM1,
      ZCC_YM1_XCC,
      ZCC_YM1_XP1,
      ZCC_YCC_XM1,
      ZCC_YCC_XCC,
      ZCC_YCC_XP1,
      ZCC_YP1_XM1,
      ZCC_YP1_XCC,
      ZCC_YP1_XP1,
      ZP1_YM1_XM1,
      ZP1_YM1_XCC,
      ZP1_YM1_XP1,
      ZP1_YCC_XM1,
      ZP1_YCC_XCC,
      ZP1_YCC_XP1,
      ZP1_YP1_XM1,
      ZP1_YP1_XCC,
      ZP1_YP1_XP1
   };
   
   const uint MISSING_ZNEG = (1 << projects::ZM1_YCC_XCC);
   const uint MISSING_YNEG = (1 << projects::ZCC_YM1_XCC);
   const uint MISSING_XNEG = (1 << projects::ZCC_YCC_XM1);
   const uint MISSING_XPOS = (1 << projects::ZCC_YCC_XP1);
   const uint MISSING_YPOS = (1 << projects::ZCC_YP1_XCC);
   const uint MISSING_ZPOS = (1 << projects::ZP1_YCC_XCC);
   const uint FACE_NBR_BITMASK = (MISSING_ZNEG | MISSING_YNEG | MISSING_XNEG | MISSING_XPOS | MISSING_YPOS | MISSING_ZPOS);
}

#endif
