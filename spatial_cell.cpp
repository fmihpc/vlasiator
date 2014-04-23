/*!
Space for static variables of spatial cell class for Vlasiator.

Copyright 2011, 2012 Finnish Meteorological Institute












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
   int SpatialCell::existing_spatial_cells = 0;
   unsigned int SpatialCell::max_velocity_blocks = 0;
   uint64_t SpatialCell::mpi_transfer_type = 0;
   bool SpatialCell::mpiTransferAtSysBoundaries = false;
}

