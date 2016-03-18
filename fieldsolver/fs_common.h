/*
This file is part of Vlasiator.

Copyright 2010-2015 Finnish Meteorological Institute

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

// FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 3, 2> & perBGrid,
// FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 3, 2> & perBDt2Grid,
// FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 3, 2> & EGrid,
// FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 3, 2> & EDt2Grid,
// FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, 3, 2> & EHallGrid,
// FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, 3, 2> & EGradPeGrid,
// FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 3, 2> & momentsGrid,
// FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 3, 2> & momentsDt2Grid,
// FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, 3, 2> & dPerBGrid,
// FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 3, 2> & dMomentsGrid,
// FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, 3, 2> & BgBGrid,
// FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, 3, 2> & volGrid,
// FsGrid< fsgrids::technical, 3, 2> & technicalGrid,

/*! Namespace containing enums and structs for the various field solver grid instances
 * 
 */
namespace fsgrids {
   enum bfield {
      PERBX,  /*!< Perturbed Magnetic field x-component, averaged over cell x-face. Propagated by field solver.*/
      PERBY,  /*!< Perturbed Magnetic field y-component, averaged over cell y-face. Propagated by field solver.*/
      PERBZ,  /*!< Perturbed Magnetic field z-component, averaged over cell z-face. Propagated by field solver.*/
      N_BFIELD
   };
   
   enum efield {
      EX,     /*!< Total electric field x-component, averaged over cell edge. Used to propagate BX,BY,BZ.*/
      EY,     /*!< Total electric field y-component, averaged over cell edge. Used to propagate BX,BY,BZ.*/
      EZ,     /*!< Total electric field z-component, averaged over cell edge. Used to propagate BX,BY,BZ.*/
      N_EFIELD
   };
   
   enum ehall {
      EXHALL_000_100,   /*!< Hall term x averaged along x on -y/-z edge of spatial cell.*/
      EYHALL_000_010,   /*!< Hall term y averaged along y on -x/-z edge of spatial cell.*/
      EZHALL_000_001,   /*!< Hall term z averaged along z on -x/-y edge of spatial cell.*/
      EYHALL_100_110,   /*!< Hall term y averaged along y on +x/-z edge of spatial cell.*/
      EZHALL_100_101,   /*!< Hall term z averaged along z on +x/-y edge of spatial cell.*/
      EXHALL_010_110,   /*!< Hall term x averaged along x on +y/-z edge of spatial cell.*/
      EZHALL_010_011,   /*!< Hall term z averaged along z on +y/-x edge of spatial cell.*/
      EZHALL_110_111,   /*!< Hall term z averaged along z on +x/+y edge of spatial cell.*/
      EXHALL_001_101,   /*!< Hall term x averaged along x on -y/+z edge of spatial cell.*/
      EYHALL_001_011,   /*!< Hall term y averaged along y on -x/+z edge of spatial cell.*/
      EYHALL_101_111,   /*!< Hall term y averaged along y on +x/+z edge of spatial cell.*/
      EXHALL_011_111,   /*!< Hall term x averaged along x on +y/+z edge of spatial cell.*/
      N_EHALL
   };
   
   enum egradpe {
      EXGRADPE,         /*!< Electron pressure gradient term x.*/
      EYGRADPE,         /*!< Electron pressure gradient term y.*/
      EZGRADPE,         /*!< Electron pressure gradient term z.*/
      N_EGRADPE
   };
   
   enum moments {
      RHO,    /*!< Number density. Calculated by Vlasov propagator, used to propagate BX,BY,BZ.*/
      RHOVX,  /*!< x-component of number density times Vx. Calculated by Vlasov propagator, used to propagate BX,BY,BZ.*/
      RHOVY,  /*!< y-component of number density times Vy. Calculated by Vlasov propagator, used to propagate BX,BY,BZ.*/
      RHOVZ,  /*!< z-component of number density times Vz. Calculated by Vlasov propagator, used to propagate BX,BY,BZ.*/
      P_11,     /*!< Pressure P_xx component, computed by Vlasov propagator. */
      P_22,     /*!< Pressure P_yy component, computed by Vlasov propagator. */
      P_33,     /*!< Pressure P_zz component, computed by Vlasov propagator. */
      N_MOMENTS
   };
   
   enum dperb {
      dPERBxdy,     /*!< Derivative of face-averaged Bx to y-direction. */
      dPERBxdz,     /*!< Derivative of face-averaged Bx to z-direction. */
      dPERBydx,     /*!< Derivative of face-averaged By to x-direction. */
      dPERBydz,     /*!< Derivative of face-averaged By to z-direction. */
      dPERBzdx,     /*!< Derivative of face-averaged Bz to x-direction. */
      dPERBzdy,     /*!< Derivative of face-averaged Bz to y-direction. */
      dPERBxdyy,     /*!< Second derivative of face-averaged Bx to yy-direction. */
      dPERBxdzz,     /*!< Second derivative of face-averaged Bx to zz-direction. */
      dPERBxdyz,     /*!< Second derivative of face-averaged Bx to yz-direction. */
      dPERBydxx,     /*!< Second derivative of face-averaged By to xx-direction. */
      dPERBydzz,     /*!< Second derivative of face-averaged By to zz-direction. */
      dPERBydxz,     /*!< Second derivative of face-averaged By to xz-direction. */
      dPERBzdxx,     /*!< Second derivative of face-averaged Bz to xx-direction. */
      dPERBzdyy,     /*!< Second derivative of face-averaged Bz to yy-direction. */
      dPERBzdxy,     /*!< Second derivative of face-averaged Bz to xy-direction. */
      N_DPERB
   }
   
   enum dmoments {
      drhodx,    /*!< Derivative of volume-averaged number density to x-direction. */
      drhody,    /*!< Derivative of volume-averaged number density to y-direction. */
      drhodz,    /*!< Derivative of volume-averaged number density to z-direction. */
      dp11dx,        /*!< Derivative of P_11 to x direction. */
      dp11dy,        /*!< Derivative of P_11 to x direction. */
      dp11dz,        /*!< Derivative of P_11 to x direction. */
      dp22dx,        /*!< Derivative of P_22 to y direction. */
      dp22dy,        /*!< Derivative of P_22 to y direction. */
      dp22dz,        /*!< Derivative of P_22 to y direction. */
      dp33dx,        /*!< Derivative of P_33 to z direction. */
      dp33dy,        /*!< Derivative of P_33 to z direction. */
      dp33dz,        /*!< Derivative of P_33 to z direction. */
      dVxdx,     /*!< Derivative of volume-averaged Vx to x-direction. */
      dVxdy,     /*!< Derivative of volume-averaged Vx to y-direction. */
      dVxdz,     /*!< Derivative of volume-averaged Vx to z-direction. */
      dVydx,     /*!< Derivative of volume-averaged Vy to x-direction. */
      dVydy,     /*!< Derivative of volume-averaged Vy to y-direction. */
      dVydz,     /*!< Derivative of volume-averaged Vy to z-direction. */
      dVzdx,     /*!< Derivative of volume-averaged Vz to x-direction. */
      dVzdy,     /*!< Derivative of volume-averaged Vz to y-direction. */
      dVzdz,     /*!< Derivative of volume-averaged Vz to z-direction. */
      N_DMOMENTS
   };
   
   // NOTE This contains the BGB derivatives as they do not change either
   enum bgbfield {
      BGBX,   /*!< Background magnetic field x-component, averaged over cell x-face.*/
      BGBY,   /*!< Background magnetic field x-component, averaged over cell x-face.*/
      BGBZ,   /*!< Background magnetic field x-component, averaged over cell x-face.*/
      BGBXVOL,   /*!< background magnetic field averaged over spatial cell.*/
      BGBYVOL,   /*!< background magnetic field averaged over spatial cell.*/
      BGBZVOL,   /*!< background magnetic field averaged over spatial cell.*/
      BGBX_000_010,   /*!< Background Bx averaged along y on -x/-z edge of spatial cell (for Hall term only).*/
      BGBX_100_110,   /*!< Background Bx averaged along y on +x/-z edge of spatial cell (for Hall term only).*/
      BGBX_001_011,   /*!< Background Bx averaged along y on -x/+z edge of spatial cell (for Hall term only).*/
      BGBX_101_111,   /*!< Background Bx averaged along y on +x/+z edge of spatial cell (for Hall term only).*/
      BGBX_000_001,   /*!< Background Bx averaged along z on -x/-y edge of spatial cell (for Hall term only).*/
      BGBX_100_101,   /*!< Background Bx averaged along z on +x/-y edge of spatial cell (for Hall term only).*/
      BGBX_010_011,   /*!< Background Bx averaged along z on +y/-x edge of spatial cell (for Hall term only).*/
      BGBX_110_111,   /*!< Background Bx averaged along z on +x/+y edge of spatial cell (for Hall term only).*/
      BGBY_000_100,   /*!< Background By averaged along x on -y/-z edge of spatial cell (for Hall term only).*/
      BGBY_010_110,   /*!< Background By averaged along x on +y/-z edge of spatial cell (for Hall term only).*/
      BGBY_001_101,   /*!< Background By averaged along x on -y/+z edge of spatial cell (for Hall term only).*/
      BGBY_011_111,   /*!< Background By averaged along x on +y/+z edge of spatial cell (for Hall term only).*/
      BGBY_000_001,   /*!< Background By averaged along z on -x/-y edge of spatial cell (for Hall term only).*/
      BGBY_100_101,   /*!< Background By averaged along z on +x/-y edge of spatial cell (for Hall term only).*/
      BGBY_010_011,   /*!< Background By averaged along z on +y/-x edge of spatial cell (for Hall term only).*/
      BGBY_110_111,   /*!< Background By averaged along z on +x/+y edge of spatial cell (for Hall term only).*/
      BGBZ_000_100,   /*!< Background Bz averaged along x on -y/-z edge of spatial cell (for Hall term only).*/
      BGBZ_010_110,   /*!< Background Bz averaged along x on +y/-z edge of spatial cell (for Hall term only).*/
      BGBZ_001_101,   /*!< Background Bz averaged along x on -y/+z edge of spatial cell (for Hall term only).*/
      BGBZ_011_111,   /*!< Background Bz averaged along x on +y/+z edge of spatial cell (for Hall term only).*/
      BGBZ_000_010,   /*!< Background Bz averaged along y on -x/-z edge of spatial cell (for Hall term only).*/
      BGBZ_100_110,   /*!< Background Bz averaged along y on +x/-z edge of spatial cell (for Hall term only).*/
      BGBZ_001_011,   /*!< Background Bz averaged along y on -x/+z edge of spatial cell (for Hall term only).*/
      BGBZ_101_111,   /*!< Background Bz averaged along y on +x/+z edge of spatial cell (for Hall term only).*/
      dBGBxdy,     /*!< Derivative of face-averaged Bx to y-direction. */
      dBGBxdz,     /*!< Derivative of face-averaged Bx to z-direction. */
      dBGBydx,     /*!< Derivative of face-averaged By to x-direction. */
      dBGBydz,     /*!< Derivative of face-averaged By to z-direction. */
      dBGBzdx,     /*!< Derivative of face-averaged Bz to x-direction. */
      dBGBzdy,     /*!< Derivative of face-averaged Bz to y-direction. */
      dBGBXVOLdy,
      dBGBXVOLdz,
      dBGBYVOLdx,
      dBGBYVOLdz,
      dBGBZVOLdx,
      dBGBZVOLdy,
      N_BGB
   };
   
   // NOTE This contains the PERBVOL derivatives
   enum volfields {
      PERBXVOL,  /*!< perturbed magnetic field  PERBX averaged over spatial cell.*/
      PERBYVOL,  /*!< perturbed magnetic field  PERBY averaged over spatial cell.*/
      PERBZVOL,  /*!< perturbed magnetic field  PERBZ averaged over spatial cell.*/
      EXVOL,     /*!< Ex averaged over spatial cell.*/
      EYVOL,     /*!< Ey averaged over spatial cell.*/
      EZVOL,     /*!< Ez averaged over spatial cell.*/
      dPERBXVOLdy,
      dPERBXVOLdz,
      dPERBYVOLdx,
      dPERBYVOLdz,
      dPERBZVOLdx,
      dPERBZVOLdy,
      N_VOL
   };
   
   struct technical {
      int sysBoundaryFlag;  /*!< System boundary flags. */
      int sysBoundaryLayer; /*!< System boundary layer index. */
      Real maxFsDt;         /*!< maximum timestep allowed in ordinary space by fieldsolver for this cell**/
   }
   
}

bool initializeFieldPropagator(
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 3, 2> & perBGrid,
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 3, 2> & perBDt2Grid,
   FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 3, 2> & EGrid,
   FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 3, 2> & EDt2Grid,
   FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, 3, 2> & EHallGrid,
   FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, 3, 2> & EGradPeGrid,
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 3, 2> & momentsGrid,
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 3, 2> & momentsDt2Grid,
   FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, 3, 2> & dPerBGrid,
   FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 3, 2> & dMomentsGrid,
   FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, 3, 2> & BgBGrid,
   FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, 3, 2> & volGrid,
   FsGrid< fsgrids::technical, 3, 2> & technicalGrid,
   SysBoundary& sysBoundaries
   SysBoundary& sysBoundaries
);

bool initializeFieldPropagatorAfterRebalance();

bool finalizeFieldPropagator();

bool propagateFields(
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 3, 2> & perBGrid,
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 3, 2> & perBDt2Grid,
   FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 3, 2> & EGrid,
   FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 3, 2> & EDt2Grid,
   FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, 3, 2> & EHallGrid,
   FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, 3, 2> & EGradPeGrid,
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 3, 2> & momentsGrid,
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 3, 2> & momentsDt2Grid,
   FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, 3, 2> & dPerBGrid,
   FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 3, 2> & dMomentsGrid,
   FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, 3, 2> & BgBGrid,
   FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, 3, 2> & volGrid,
   FsGrid< fsgrids::technical, 3, 2> & technicalGrid,
   SysBoundary& sysBoundaries,
   creal& dt,
   cint& subcycles
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
   const FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 3, 2> & perBGrid,
   FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, 3, 2> & dPerBGrid,
   Real* perturbedResult,
   i,
   j,
   k,
   creal& reconstructionOrder
);

#endif
