/*
This file is part of Vlasiator.

Copyright 2010, 2011, 2012, 2013, 2014 Finnish Meteorological Institute

*/

#ifndef DERIVATIVES_HPP
#define DERIVATIVES_HPP

#include <vector>

#include "../spatial_cell.hpp"
#include "../sysboundary/sysboundary.h"

#include "fs_limiters.h"

void calculateDerivatives(
                          const CellID& cellID,
                          dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                          SysBoundary& sysBoundaries,
                          cint& RKCase,
                          const bool& doMoments);

void calculateDerivativesSimple(
                                dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                SysBoundary& sysBoundaries,
                                const std::vector<CellID>& localCells,
                                cint& RKCase,
                                const bool& doMoments);

void calculateBVOLDerivatives(
                              const CellID& cellID,
                              dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                              SysBoundary& sysBoundaries);

void calculateBVOLDerivativesSimple(
                                    dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                    SysBoundary& sysBoundaries,
                                    const std::vector<CellID>& localCells);

#endif
