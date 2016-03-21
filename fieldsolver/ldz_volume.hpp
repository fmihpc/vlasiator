/*
This file is part of Vlasiator.

Copyright 2010-2015 Finnish Meteorological Institute

*/

#ifndef LDZ_VOLUME
#define LDZ_VOLUME

#include <vector>

#include "fs_common.h"
#include "../spatial_cell.hpp"

/*! \brief Top-level field averaging function.
 * 
 * Averages the electric and magnetic fields over the cell volumes.
 * 
 * \sa reconstructionCoefficients
 */
void calculateVolumeAveragedFields(
   const FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2> & perBGrid,
   const FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 2> & EGrid,
   const FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, 2> & dPerBGrid,
   FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, 2> & volGrid,
   const FsGrid< fsgrids::technical, 2> & technicalGrid
);

#endif
