/*
This file is part of Vlasiator.

Copyright 2010, 2011, 2012, 2013, 2014 Finnish Meteorological Institute

*/

#ifndef DERIVATIVES_HPP
#define DERIVATIVES_HPP

#include "fs_limiters.h"

void calculateDerivatives(
                          const CellID& cellID,
                          dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                          SysBoundary& sysBoundaries,
                          cint& RKCase,
                          const bool& doMoments);

void calculateDerivativesSimple(
                                dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                SysBoundary& sysBoundaries,
                                const vector<CellID>& localCells,
                                cint& RKCase,
                                const bool& doMoments);

void calculateBVOLDerivatives(
                              const CellID& cellID,
                              dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                              SysBoundary& sysBoundaries);

void calculateBVOLDerivativesSimple(
                                    dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                    SysBoundary& sysBoundaries,
                                    const vector<CellID>& localCells);

#endif
