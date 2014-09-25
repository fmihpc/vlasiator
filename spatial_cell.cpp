/*!
Space for static variables of spatial cell class for Vlasiator.

Copyright 2011, 2012 Finnish Meteorological Institute












*/

#include "spatial_cell.hpp"

namespace spatial_cell {
   Real SpatialCell::velocity_block_min_value = 0;    
   uint64_t SpatialCell::mpi_transfer_type = 0;
   bool SpatialCell::mpiTransferAtSysBoundaries = false;
}

