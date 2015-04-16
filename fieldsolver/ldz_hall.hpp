/*
This file is part of Vlasiator.

Copyright 2010, 2011, 2012, 2013, 2014 Finnish Meteorological Institute

*/

#ifndef LDZ_HALL_HPP
#define LDZ_HALL_HPP

void calculateEdgeHallTermXComponents(Real* cp,Real* derivs,
                                      const Real* const perturbedCoefficients,
                                      cint& RKCase);

void calculateEdgeHallTermYComponents(Real* cp,Real* derivs,
                                      const Real* const perturbedCoefficients,
                                      cint& RKCase);

void calculateEdgeHallTermZComponents(Real* cp,Real* derivs,
                                      const Real* const perturbedCoefficients,
                                      cint& RKCase);

void calculateHallTerm(SysBoundary& sysBoundaries,std::vector<fs_cache::CellCache>& cache,
                       const std::vector<uint16_t>& cells,cint& RKCase);

void calculateHallTermSimple(
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   SysBoundary& sysBoundaries,
   const vector<CellID>& localCells,
   cint& RKCase,
   const bool communicateDerivatives
);

#endif
