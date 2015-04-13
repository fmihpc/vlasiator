/*
This file is part of Vlasiator.

Copyright 2010-2015 Finnish Meteorological Institute

*/

#ifndef LDZ_VOLUME
#define LDZ_VOLUME

#include <vector>

#include <dccrg.hpp>
#include <dccrg_cartesian_geometry.hpp>

#include "fs_common.h"
#include "../spatial_cell.hpp"

/*! \brief Top-level field averaging function.
 * 
 * Averages the electric and magnetic fields over the cell volumes.
 * 
 * \sa reconstructionCoefficients
 */
void calculateVolumeAveragedFields(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                   std::vector<fs_cache::CellCache>& cache,
                                   const std::vector<uint16_t>& cells);

#endif
