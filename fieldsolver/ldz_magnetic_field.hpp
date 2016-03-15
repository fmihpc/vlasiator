/*
This file is part of Vlasiator.

Copyright 2010, 2011, 2012, 2013, 2014 Finnish Meteorological Institute

*/

#ifndef LDZ_MAGNETIC_FIELD_HPP
#define LDZ_MAGNETIC_FIELD_HPP

#include <vector>

#include <dccrg.hpp>
#include <dccrg_cartesian_geometry.hpp>

#include "../definitions.h"
#include "../common.h"
#include "../spatial_cell.hpp"

#include "fs_common.h"

void propagateMagneticField(
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 3, 2> & perBGrid,
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 3, 2> & perBDt2Grid,
   FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 3, 2> & EGrid,
   FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 3, 2> & EDt2Grid,
   const std::vector<uint16_t>& cells,
   creal& dt,
   cint& RKCase,
   const bool doX=true,
   const bool doY=true,
   const bool doZ=true
);

void propagateSysBoundaryMagneticField(
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 3, 2> & perBGrid,
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 3, 2> & perBDt2Grid,
   FsGrid< fsgrids::technical, 3, 2> & technicalGrid,
   const uint16_t& localID,
   SysBoundary& sysBoundaries,
   creal& dt,
   cint& RKCase
);

void propagateMagneticFieldSimple(
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 3, 2> & perBGrid,
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 3, 2> & perBDt2Grid,
   FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 3, 2> & EGrid,
   FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 3, 2> & EDt2Grid,
   FsGrid< fsgrids::technical, 3, 2> & technicalGrid,
   SysBoundary& sysBoundaries,
   creal& dt,
   const std::vector<CellID>& localCells,
   cint& RKCase
);

#endif
