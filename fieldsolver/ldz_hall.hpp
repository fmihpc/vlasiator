/*
This file is part of Vlasiator.

Copyright 2010, 2011, 2012, 2013, 2014 Finnish Meteorological Institute

*/

#ifndef LDZ_HALL_HPP
#define LDZ_HALL_HPP

// NOTE Not needed?
// void calculateEdgeHallTermXComponents(
//    Real* cp,
//    Real* derivs,
//    const Real* const perturbedCoefficients,
//    cint& RKCase
// );
// 
// void calculateEdgeHallTermYComponents(
//    Real* cp,
//    Real* derivs,
//    const Real* const perturbedCoefficients,
//    cint& RKCase
// );
// 
// void calculateEdgeHallTermZComponents(
//    Real* cp,
//    Real* derivs,
//    const Real* const perturbedCoefficients,
//    cint& RKCase
// );
// 
// void calculateHallTerm(
//    SysBoundary& sysBoundaries,
//    std::vector<fs_cache::CellCache>& cache,
//    const std::vector<uint16_t>& cells,
//    cint& RKCase
// );

void calculateHallTermSimple(
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 3, 2> & perBGrid,
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 3, 2> & perBDt2Grid,
   FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, 3, 2> & EHallGrid,
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 3, 2> & momentsGrid,
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 3, 2> & momentsDt2Grid,
   FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, 3, 2> & dPerBGrid,
   FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 3, 2> & dMomentsGrid,
   FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, 3, 2> & BgBGrid,
   FsGrid< fsgrids::technical, 3, 2> & technicalGrid,
   SysBoundary& sysBoundaries,
   const vector<CellID>& localCells,
   cint& RKCase,
   const bool communicateDerivatives
);

#endif
