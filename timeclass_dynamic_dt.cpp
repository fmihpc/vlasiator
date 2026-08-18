/*
 * This file is part of Vlasiator.
 * Copyright 2010-2024 Finnish Meteorological Institute and University of Helsinki
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

#include "timeclass_dynamic_dt.hpp"

/*
dynamic dt working principle with timeclasses:

Each cell must check if its current timestep given by its timeclass is too large, every fractionaltimestep it gets propagated (get_timeclass_turn==true)
To correct a too large timestep, a rollback and re-initialization is needed. 
However, this rollback must be done for the entire simulation at once, to keep the cells synced. 
If fractionaltimestep is not zero, cells of different timeclasses are not synced. They are only synced at fractionaltimestep 0.
Therefore if a cell's timestep is too small, and fractionaltimestep is not zero, a rollback is not possible.

To facilitate this, we do the following:

1. If a cell has a too large timestep, and fractionaltimestep is not zero, instead of rolling back, we increase the cell's timeclass to one higher (a new highest timeclass is created if necessary).
    The cell is then propagated with the new timeclass until fractimestep is zero.

2. When fractimestep is zero again, we do a global rollback and re-initialization, where all cells get allocated into timeclasses based on their current timestep limits. Any new highest timeclass is removed. 
    - make sure to catch the edge case where the need for rollback happens at fractionaltimestep 0, then we do the rollback immediately.

*/

/*

computeNewTimestep - already exists, computes timeclass dts, assigns cells timeclasses, to be done at fractimestep 0.

increaseTimeclass - increases timeclass of a cell, to be done at fractimestep > 0, if cell's dt is too large.
*/


