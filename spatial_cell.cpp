/*!
Space for static variables of spatial cell class for Vlasiator.

Copyright 2011, 2012 Finnish Meteorological Institute

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

#include "spatial_cell.hpp"

namespace spatial_cell {
   unsigned int SpatialCell::vx_length = 0;
   unsigned int SpatialCell::vy_length = 0;
   unsigned int SpatialCell::vz_length = 0;
   Real SpatialCell::vx_min = 0;
   Real SpatialCell::vx_max = 0;
   Real SpatialCell::vy_min = 0;
   Real SpatialCell::vy_max = 0;
   Real SpatialCell::vz_min = 0;
   Real SpatialCell::vz_max = 0;
   Real SpatialCell::grid_dvx = 0;
   Real SpatialCell::grid_dvy = 0;
   Real SpatialCell::grid_dvz = 0;
   Real SpatialCell::block_dvx = 0;
   Real SpatialCell::block_dvy = 0;
   Real SpatialCell::block_dvz = 0;
   Real SpatialCell::cell_dvx = 0;
   Real SpatialCell::cell_dvy = 0;
   Real SpatialCell::cell_dvz = 0;
   Real SpatialCell::velocity_block_min_value = 0;    
   Real SpatialCell::velocity_block_min_avg_value = 0;
   unsigned int SpatialCell::max_velocity_blocks = 0;
   unsigned int SpatialCell::mpi_transfer_type = 0;
}

