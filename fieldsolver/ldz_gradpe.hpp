/*
This file is part of Vlasiator.

Copyright 2010, 2011, 2012, 2013, 2014 Finnish Meteorological Institute

*/

#ifndef LDZ_GRADPE_HPP
#define LDZ_GRADPE_HPP

// NOTE not needed??
// void calculateEdgeGradPeTermXComponents(
//    Real* cp,Real* derivs,
//    const Real* const perturbedCoefficients,
//    cint& RKCase
// );
// 
// void calculateEdgeGradPeTermYComponents(
//    Real* cp,
//    Real* derivs,
//    const Real* const perturbedCoefficients,
//    cint& RKCase
// );
// 
// void calculateEdgeGradPeTermZComponents(
//    Real* cp,
//    Real* derivs,
//    const Real* const perturbedCoefficients,
//    cint& RKCase
// );
// 
// void calculateGradPeTerm(
//    SysBoundary& sysBoundaries,
//    std::vector<fs_cache::CellCache>& cache,
//    const std::vector<uint16_t>& cells,
//    cint& RKCase
// );

void calculateGradPeTermSimple(
   FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, 3, 2> & EGradPeGrid,
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 3, 2> & momentsGrid,
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 3, 2> & momentsDt2Grid,
   FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 3, 2> & dMomentsGrid,
   FsGrid< fsgrids::technical, 3, 2> & technicalGrid,
   SysBoundary& sysBoundaries,
   const vector<CellID>& localCells,
   cint& RKCase
);

#endif
