/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://www.physics.helsinki.fi/vlasiator/
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#ifndef FS_COMMON_H
#define FS_COMMON_H

#include <cstdlib>
#include <iostream>
#include <cmath>
#include <vector>
#include <map>
#include <list>
#include <set>
#include <stdint.h>

#include <fsgrid.hpp>
#include <phiprof.hpp>

#include "../common.h"
#include "../parameters.h"
#include "../projects/project.h"
#include "../sysboundary/sysboundary.h"
#include "../sysboundary/sysboundarycondition.h"

// Constants: not needed as such, but if field solver is implemented on GPUs 
// these force CPU to use float accuracy, which in turn helps to compare 
// CPU and GPU results.

const Real HALF    = 0.5;
const Real MINUS   = -1.0;
const Real PLUS    = +1.0;
const Real THIRD   = 1.0/3.0;
const Real FOURTH  = 1.0/4.0;
const Real SIXTH   = 1.0/6.0;
const Real EIGTH   = 1.0/8.0;
const Real TENTH   = 1.0/10.0;
const Real TWELWTH = 1.0/12.0;
const Real TWO     = 2.0;
const Real ZERO    = 0.0;

static creal EPS = 1.0e-30;

using namespace std;
using namespace fieldsolver;

bool initializeFieldPropagator(
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2> & perBGrid,
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2> & perBDt2Grid,
   FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 2> & EGrid,
   FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 2> & EDt2Grid,
   FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, 2> & EHallGrid,
   FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, 2> & EGradPeGrid,
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 2> & momentsGrid,
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 2> & momentsDt2Grid,
   FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, 2> & dPerBGrid,
   FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 2> & dMomentsGrid,
   FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, 2> & BgBGrid,
   FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, 2> & volGrid,
   FsGrid< fsgrids::technical, 2> & technicalGrid,
   SysBoundary& sysBoundaries
);

bool initializeFieldPropagatorAfterRebalance();

bool finalizeFieldPropagator();

bool propagateFields(
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2> & perBGrid,
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2> & perBDt2Grid,
   FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 2> & EGrid,
   FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 2> & EDt2Grid,
   FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, 2> & EHallGrid,
   FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, 2> & EGradPeGrid,
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 2> & momentsGrid,
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 2> & momentsDt2Grid,
   FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, 2> & dPerBGrid,
   FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 2> & dMomentsGrid,
   FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, 2> & BgBGrid,
   FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, 2> & volGrid,
   FsGrid< fsgrids::technical, 2> & technicalGrid,
   SysBoundary& sysBoundaries,
   creal& dt,
   cuint subcycles
);

Real divideIfNonZero(creal rhoV, creal rho);

/*! Namespace encompassing the enum defining the list of reconstruction coefficients used in field component reconstructions.*/
namespace Rec {
   /*! Enum defining the list of reconstruction coefficients used in field component reconstructions.*/
   enum Rec {
      a_0, a_x, a_y, a_z, a_xx, a_yy, a_zz, a_xy, a_xz, a_yz, a_xxx, a_xxy, a_xyy, a_xxz, a_xzz, a_xyz,
      b_0, b_x, b_y, b_z, b_xx, b_yy, b_zz, b_xy, b_xz, b_yz, b_xxy, b_xyy, b_yyy, b_yyz, b_yzz, b_xyz,
      c_0, c_x, c_y, c_z, c_xx, c_yy, c_zz, c_xy, c_xz, c_yz, c_xxz, c_xzz, c_yyz, c_yzz, c_xyz, c_zzz,
      N_REC_COEFFICIENTS
   };
}

void reconstructionCoefficients(
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2> & perBGrid,
   FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, 2> & dPerBGrid,
   Real* perturbedResult,
   cint i,
   cint j,
   cint k,
   creal& reconstructionOrder
);

#endif
