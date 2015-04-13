/*
This file is part of Vlasiator.

Copyright 2010-2015 Finnish Meteorological Institute

*/

#ifndef LDZ_ELECTRIC_FIELD_HPP
#define LDZ_ELECTRIC_FIELD_HPP

#include "fs_common.h"

Real calculateWaveSpeedYZ(
   const Real* cp,
   const Real* derivs,
   const Real* nbr_cp,
   const Real* nbr_derivs,
   const Real& By,
   const Real& Bz,
   const Real& dBydx,
   const Real& dBydz,
   const Real& dBzdx,
   const Real& dBzdy,
   const Real& ydir,
   const Real& zdir,
   cint& RKCase
);

Real calculateWaveSpeedXZ(
   const Real* cp,
   const Real* derivs,
   const Real* nbr_cp,
   const Real* nbr_derivs,
   const Real& Bx,
   const Real& Bz,
   const Real& dBxdy,
   const Real& dBxdz,
   const Real& dBzdx,
   const Real& dBzdy,
   const Real& xdir,
   const Real& zdir,
   cint& RKCase
);

Real calculateWaveSpeedXY(
   const Real* cp,
   const Real* derivs,
   const Real* nbr_cp,
   const Real* nbr_derivs,
   const Real& Bx,
   const Real& By,
   const Real& dBxdy,
   const Real& dBxdz,
   const Real& dBydx,
   const Real& dBydz,
   const Real& xdir,
   const Real& ydir,
   cint& RKCase
);

void calculateUpwindedElectricFieldSimple(
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   SysBoundary& sysBoundaries,
   const vector<CellID>& localCells,
   cint& RKCase
);

#endif
