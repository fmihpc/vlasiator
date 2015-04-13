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

void propagateMagneticField(const CellID& cellID,
                            dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                            creal& dt,cint& RKCase); 

void propagateSysBoundaryMagneticField(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                       const CellID& cellID,SysBoundary& sysBoundaries,
                                       creal& dt,cint& RKCase);

void propagateMagneticFieldSimple(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                  SysBoundary& sysBoundaries,creal& dt,
                                  const std::vector<CellID>& localCells,cint& RKCase);

#endif
