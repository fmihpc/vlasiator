/*
This file is part of Vlasiator.

Copyright 2010-2015 Finnish Meteorological Institute

*/

#ifndef LDZ_ELECTRIC_FIELD_HPP
#define LDZ_ELECTRIC_FIELD_HPP

#include "fs_common.h"

void calculateUpwindedElectricFieldSimple(
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   SysBoundary& sysBoundaries,
   const vector<CellID>& localCells,
   cint& RKCase
);

#endif
