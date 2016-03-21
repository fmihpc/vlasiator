/*
This file is part of Vlasiator.

Copyright 2010, 2011, 2012, 2013, 2014 Finnish Meteorological Institute

*/

#ifndef LDZ_GRADPE_HPP
#define LDZ_GRADPE_HPP

void calculateGradPeTermSimple(
   FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, 2> & EGradPeGrid,
   const FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 2> & momentsGrid,
   const FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 2> & dMomentsGrid,
   const FsGrid< fsgrids::technical, 2> & technicalGrid,
   SysBoundary& sysBoundaries,
   cint& RKCase
);

#endif
