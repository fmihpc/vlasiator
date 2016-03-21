/*
This file is part of Vlasiator.

Copyright 2010, 2011, 2012, 2013, 2014 Finnish Meteorological Institute

*/

#ifndef LDZ_MAGNETIC_FIELD_HPP
#define LDZ_MAGNETIC_FIELD_HPP

#include <vector>

#include "../definitions.h"
#include "../common.h"
#include "../spatial_cell.hpp"

#include "fs_common.h"

void propagateMagneticField(
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 3, 2> & perBGrid,
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 3, 2> & perBDt2Grid,
   FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 3, 2> & EGrid,
   FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 3, 2> & EDt2Grid,
   cint i,
   cint j,
   cint k,
   creal& dt,
   cint& RKCase,
   const bool doX=true,
   const bool doY=true,
   const bool doZ=true
);

void propagateMagneticFieldSimple(
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 3, 2> & perBGrid,
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 3, 2> & perBDt2Grid,
   FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 3, 2> & EGrid,
   FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 3, 2> & EDt2Grid,
   FsGrid< fsgrids::technical, 3, 2> & technicalGrid,
   SysBoundary& sysBoundaries,
   creal& dt,
   cint& RKCase
);

#endif
