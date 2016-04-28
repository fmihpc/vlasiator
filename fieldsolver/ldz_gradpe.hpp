/*
This file is part of Vlasiator.

Copyright 2010, 2011, 2012, 2013, 2014 Finnish Meteorological Institute

*/

#ifndef LDZ_GRADPE_HPP
#define LDZ_GRADPE_HPP

void calculateEdgeGradPeTermXComponents(
   Real* cp,Real* derivs,
   const Real* const perturbedCoefficients,
   cint& RKCase
);

void calculateEdgeGradPeTermYComponents(
   Real* cp,
   Real* derivs,
   const Real* const perturbedCoefficients,
   cint& RKCase
);

void calculateEdgeGradPeTermZComponents(
   Real* cp,
   Real* derivs,
   const Real* const perturbedCoefficients,
   cint& RKCase
);

void calculateGradPeTerm(
   SysBoundary& sysBoundaries,
   std::vector<fs_cache::CellCache>& cache,
   const std::vector<uint16_t>& cells,
   cint& RKCase
);

void calculateGradPeTermSimple(
   dccrg::Dccrg<SpatialCell,
   dccrg::Cartesian_Geometry>& mpiGrid,
   SysBoundary& sysBoundaries,
   const vector<CellID>& localCells,
   cint& RKCase
);

#endif
