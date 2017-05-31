/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://vlasiator.fmi.fi/
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#ifndef COMMON_H
#define COMMON_H

#include <limits>
#include <string>
#include <vector>
#include "definitions.h"

#ifdef DEBUG_SOLVERS
#define CHECK_FLOAT(x) \
   if ((x) != (x)) {\
      std::cerr << __FILE__ << ":" << __LINE__ << " Illegal value: " << x << std::endl;\
      abort();\
   }
#else
#define CHECK_FLOAT(x) {}
#endif

void bailout(
   const bool condition,
   const std::string message
);
void bailout(
   const bool condition,
   const char * const file,
   const int line
);
void bailout(
   const bool condition,
   const std::string message,
   const char * const file,
   const int line
);

#define sqr(x) ((x)*(x))
#define pow2(x) sqr(x)
#define pow3(x) ((x)*(x)*(x))

#define MASTER_RANK 0


/* Maximum number of blocks in each dimension in velocity space. The
   size of velocity space defined in cfg can at maximum be this large
*/
#define MAX_BLOCKS_PER_DIM 200


/*! A namespace for storing indices into an array which contains 
 * neighbour list for each spatial cell. These indices refer to 
 * the CPU memory, i.e. the device does not use these.
 */
namespace NbrsSpa {
   const uint INNER = 0;      /*!< The cell is an inner cell, i.e. all its neighbours are located on the same computation node.*/
   const uint X_NEG_BND = (1 << 0);  /*!< The cell is a boundary cell in -x direction.*/
   const uint X_POS_BND = (1 << 1);  /*!< The cell is a boundary cell in +x direction.*/
   const uint Y_NEG_BND = (1 << 2);  /*!< The cell is a boundary cell in -y direction.*/
   const uint Y_POS_BND = (1 << 3);  /*!< The cell is a boundary cell in +y direction.*/
   const uint Z_NEG_BND = (1 << 4); /*!< The cell is a boundary cell in -z direction.*/
   const uint Z_POS_BND = (1 << 5); /*!< The cell is a boundary cell in +z direction.*/
   
   enum {
      STATE, /*!< Contains the neighbour information of this cell, i.e. whether it is an inner cell or a boundary cell in one or more coordinate directions.*/
      MYIND, /*!< The index of this cell.*/
      X1NEG,  /*!< The index of the -x neighbouring block, distance 1.*/
      Y1NEG,  /*!< The index of the -y neighbouring block, distance 1.*/
      Z1NEG,  /*!< The index of the -z neighbouring block, distance 1.*/
      X1POS,  /*!< The index of the +x neighbouring block, distance 1.*/
      Y1POS,  /*!< The index of the +y neighbouring block, distance 1.*/
      Z1POS,  /*!< The index of the +z neighbouring block, distance 1.*/
      X2NEG,  /*!< The index of the -x neighbouring block, distance 1.*/
      Y2NEG,  /*!< The index of the -y neighbouring block, distance 1.*/
      Z2NEG,  /*!< The index of the -z neighbouring block, distance 1.*/
      X2POS,  /*!< The index of the +x neighbouring block, distance 1.*/
      Y2POS,  /*!< The index of the +y neighbouring block, distance 1.*/
      Z2POS   /*!< The index of the +z neighbouring block, distance 1.*/
   };
}


/*! A namespace for storing indices into an array which contains 
 * the physical parameters of each velocity block.*/
namespace BlockParams {
   enum {
      VXCRD,   /*!< vx-coordinate of the bottom left corner of the block.*/
      VYCRD,   /*!< vy-coordinate of the bottom left corner of the block.*/
      VZCRD,   /*!< vz-coordinate of the bottom left corner of the block.*/
      DVX,     /*!< Grid separation in vx-coordinate for the block.*/
      DVY,     /*!< Grid separation in vy-coordinate for the block.*/
      DVZ,     /*!< Grid separation in vz-coordinate for the block.*/
      N_VELOCITY_BLOCK_PARAMS
   };
}

/*! A namespace for storing indices into an array which contains the 
 * physical parameters of each spatial cell. Do not change the order 
 * of variables unless you know what you are doing - MPI transfers in 
 * field solver are optimised for this particular ordering.
 */
namespace CellParams {
#warning RHOMV not needed, switch to V?
   enum {
      XCRD,   /*!< x-coordinate of the bottom left corner.*/
      YCRD,   /*!< y-coordinate of the bottom left corner.*/
      ZCRD,   /*!< z-coordinate of the bottom left corner.*/
      DX,     /*!< Grid separation in x-coordinate.*/
      DY,     /*!< Grid separation in y-coordinate.*/
      DZ,     /*!< Grid separation in z-coordinate.*/
      EX,     /*!< Total electric field x-component, averaged over cell edge. Used to propagate BX,BY,BZ.*/
      EY,     /*!< Total wlectric field y-component, averaged over cell edge. Used to propagate BX,BY,BZ.*/
      EZ,     /*!< Total electric field z-component, averaged over cell edge. Used to propagate BX,BY,BZ.*/
      BGBX,   /*!< Background magnetic field x-component, averaged over cell x-face.*/
      BGBY,   /*!< Background magnetic field x-component, averaged over cell x-face.*/
      BGBZ,   /*!< Background magnetic field x-component, averaged over cell x-face.*/
      PERBX,  /*!< Perturbed Magnetic field x-component, averaged over cell x-face. Propagated by field solver.*/
      PERBY,  /*!< Perturbed Magnetic field y-component, averaged over cell y-face. Propagated by field solver.*/
      PERBZ,  /*!< Perturbed Magnetic field z-component, averaged over cell z-face. Propagated by field solver.*/
      RHO,    /*!< Number density. Calculated by Vlasov propagator, used to propagate BX,BY,BZ.*/
      RHOVX,  /*!< x-component of number density times Vx. Calculated by Vlasov propagator, used to propagate BX,BY,BZ.*/
      RHOVY,  /*!< y-component of number density times Vy. Calculated by Vlasov propagator, used to propagate BX,BY,BZ.*/
      RHOVZ,  /*!< z-component of number density times Vz. Calculated by Vlasov propagator, used to propagate BX,BY,BZ.*/
      RHOM,    /*!< Total mass density. Calculated by Vlasov propagator, used to propagate BX,BY,BZ.*/
      RHOMVX,  /*!< x-component of mass density times Vx. Calculated by Vlasov propagator, used to propagate BX,BY,BZ.*/
      RHOMVY,  /*!< y-component of mass density times Vy. Calculated by Vlasov propagator, used to propagate BX,BY,BZ.*/
      RHOMVZ,  /*!< z-component of mass density times Vz. Calculated by Vlasov propagator, used to propagate BX,BY,BZ.*/
      RHOQ,    /*!< Total charge density. Calculated by Vlasov propagator, used to propagate BX,BY,BZ.*/
      EX_DT2,    /*!< Intermediate step value for RK2 time stepping in field solver.*/
      EY_DT2,    /*!< Intermediate step value for RK2 time stepping in field solver.*/
      EZ_DT2,    /*!< Intermediate step value for RK2 time stepping in field solver.*/
      PERBX_DT2, /*!< Intermediate step value for PERBX for RK2 time stepping in field solver.*/
      PERBY_DT2, /*!< Intermediate step value for PERBY for RK2 time stepping in field solver.*/
      PERBZ_DT2, /*!< Intermediate step value for PERBZ for RK2 time stepping in field solver.*/
      RHOM_DT2,    /*!< Total mass density. Calculated by Vlasov propagator, used to propagate BX,BY,BZ.*/
      RHOMVX_DT2,  /*!< x-component of mass density times Vx. Calculated by Vlasov propagator, used to propagate BX,BY,BZ.*/
      RHOMVY_DT2,  /*!< y-component of mass density times Vy. Calculated by Vlasov propagator, used to propagate BX,BY,BZ.*/
      RHOMVZ_DT2,  /*!< z-component of mass density times Vz. Calculated by Vlasov propagator, used to propagate BX,BY,BZ.*/
      RHOQ_DT2,    /*!< Total charge density. Calculated by Vlasov propagator, used to propagate BX,BY,BZ.*/
      BGBXVOL,   /*!< background magnetic field averaged over spatial cell.*/
      BGBYVOL,   /*!< background magnetic field averaged over spatial cell.*/
      BGBZVOL,   /*!< background magnetic field averaged over spatial cell.*/
      PERBXVOL,  /*!< perturbed magnetic field  PERBX averaged over spatial cell.*/
      PERBYVOL,  /*!< perturbed magnetic field  PERBY averaged over spatial cell.*/
      PERBZVOL,  /*!< perturbed magnetic field  PERBZ averaged over spatial cell.*/
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
      EXVOL,     /*!< Ex averaged over spatial cell.*/
      EYVOL,     /*!< Ey averaged over spatial cell.*/
      EZVOL,     /*!< Ez averaged over spatial cell.*/
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
      EXGRADPE,         /*!< Electron pressure gradient term x.*/
      EYGRADPE,         /*!< Electron pressure gradient term y.*/
      EZGRADPE,         /*!< Electron pressure gradient term z.*/
      RHOM_R,     /*!< RHO after propagation in ordinary space*/
      RHOMVX_R,   /*!< RHOVX after propagation in ordinary space*/
      RHOMVY_R,   /*!< RHOVX after propagation in ordinary space*/
      RHOMVZ_R,   /*!< RHOVX after propagation in ordinary space*/
      RHOQ_R,     /*!< RHO after propagation in ordinary space*/
      RHOM_V,     /*!< RHO after propagation in velocity space*/
      RHOMVX_V,   /*!< RHOVX after propagation in velocity space*/
      RHOMVY_V,   /*!< RHOVX after propagation in velocity space*/
      RHOMVZ_V,   /*!< RHOVX after propagation in velocity space*/
      RHOQ_V,     /*!< RHO after propagation in velocity space*/
      P_11,     /*!< Pressure P_xx component, computed by Vlasov propagator. */
      P_22,     /*!< Pressure P_yy component, computed by Vlasov propagator. */
      P_33,     /*!< Pressure P_zz component, computed by Vlasov propagator. */
      P_11_DT2, /*!< Intermediate step value for RK2 time stepping in field solver. Computed from P_11_R and P_11_V. */
      P_22_DT2, /*!< Intermediate step value for RK2 time stepping in field solver. Computed from P_22_R and P_22_V. */
      P_33_DT2, /*!< Intermediate step value for RK2 time stepping in field solver. Computed from P_33_R and P_33_V. */
      P_11_R,   /*!< P_xx component after propagation in ordinary space */
      P_22_R,   /*!< P_yy component after propagation in ordinary space */
      P_33_R,   /*!< P_zz component after propagation in ordinary space */
      P_11_V,   /*!< P_xx component after propagation in velocity space */
      P_22_V,   /*!< P_yy component after propagation in velocity space */
      P_33_V,   /*!< P_zz component after propagation in velocity space */
      MAXVDT,             /*!< maximum timestep allowed in velocity space for this cell, 
                           * this is the max allowed timestep over all particle species.*/
      MAXRDT,             /*!< maximum timestep allowed in ordinary space for this cell,
                           * this is the max allowed timestep over all particle species.*/
      MAXFDT,             /*!< maximum timestep allowed in ordinary space by fieldsolver for this cell**/
      LBWEIGHTCOUNTER,    /*!< Counter for storing compute time weights needed by the load balancing**/
      ISCELLSAVINGF,      /*!< Value telling whether a cell is saving its distribution function when partial f data is written out. */
      PHI,        /*!< Electrostatic potential.*/
      PHI_TMP,    /*!< Temporary electrostatic potential.*/
      RHOQ_TOT,   /*!< Total charge density, summed over all particle populations.*/
      RHOQ_EXT,   /*<! External charge density.*/
      BGEXVOL,    /*!< Background electric field averaged over spatial cell, x-component.*/
      BGEYVOL,    /*!< Background electric field averaged over spatial cell, y-component.*/
      BGEZVOL,    /*!< Background electric field averaged over spatial cell, z-component.*/
      N_SPATIAL_CELL_PARAMS
   };
}

/*! Namespace fieldsolver contains indices into arrays which store 
 * variables required by the field solver. These quantities are derivatives 
 * of variables described in namespace CellParams.
 * Do not change the order of variables unless you know what you are doing: 
 * in several places the size of cpu_derivatives array in cell_spatial is calculated 
 * as fieldsolver::dVzdz+1.
 */
namespace fieldsolver {
   enum {
      drhomdx,    /*!< Derivative of volume-averaged mass density to x-direction. */
      drhomdy,    /*!< Derivative of volume-averaged mass density to y-direction. */
      drhomdz,    /*!< Derivative of volume-averaged mass density to z-direction. */
      drhoqdx,    /*!< Derivative of volume-averaged charge density to x-direction. */
      drhoqdy,    /*!< Derivative of volume-averaged charge density to y-direction. */
      drhoqdz,    /*!< Derivative of volume-averaged charge density to z-direction. */
      dBGBxdy,     /*!< Derivative of face-averaged Bx to y-direction. */
      dBGBxdz,     /*!< Derivative of face-averaged Bx to z-direction. */
      dBGBydx,     /*!< Derivative of face-averaged By to x-direction. */
      dBGBydz,     /*!< Derivative of face-averaged By to z-direction. */
      dBGBzdx,     /*!< Derivative of face-averaged Bz to x-direction. */
      dBGBzdy,     /*!< Derivative of face-averaged Bz to y-direction. */
      dPERBxdy,     /*!< Derivative of face-averaged Bx to y-direction. */
      dPERBxdz,     /*!< Derivative of face-averaged Bx to z-direction. */
      dPERBydx,     /*!< Derivative of face-averaged By to x-direction. */
      dPERBydz,     /*!< Derivative of face-averaged By to z-direction. */
      dPERBzdx,     /*!< Derivative of face-averaged Bz to x-direction. */
      dPERBzdy,     /*!< Derivative of face-averaged Bz to y-direction. */
      // Insert for Hall term
      // NOTE 2nd derivatives of BGBn are not needed as curl(dipole) = 0.0
      // will change if BGB is not curl-free
//       dBGBxdyy,     /*!< Second derivative of face-averaged Bx to yy-direction. */
//       dBGBxdzz,     /*!< Second derivative of face-averaged Bx to zz-direction. */
//       dBGBxdyz,     /*!< Second derivative of face-averaged Bx to yz-direction. */
//       dBGBydxx,     /*!< Second derivative of face-averaged By to xx-direction. */
//       dBGBydzz,     /*!< Second derivative of face-averaged By to zz-direction. */
//       dBGBydxz,     /*!< Second derivative of face-averaged By to xz-direction. */
//       dBGBzdxx,     /*!< Second derivative of face-averaged Bz to xx-direction. */
//       dBGBzdyy,     /*!< Second derivative of face-averaged Bz to yy-direction. */
//       dBGBzdxy,     /*!< Second derivative of face-averaged Bz to xy-direction. */
      dPERBxdyy,     /*!< Second derivative of face-averaged Bx to yy-direction. */
      dPERBxdzz,     /*!< Second derivative of face-averaged Bx to zz-direction. */
      dPERBxdyz,     /*!< Second derivative of face-averaged Bx to yz-direction. */
      dPERBydxx,     /*!< Second derivative of face-averaged By to xx-direction. */
      dPERBydzz,     /*!< Second derivative of face-averaged By to zz-direction. */
      dPERBydxz,     /*!< Second derivative of face-averaged By to xz-direction. */
      dPERBzdxx,     /*!< Second derivative of face-averaged Bz to xx-direction. */
      dPERBzdyy,     /*!< Second derivative of face-averaged Bz to yy-direction. */
      dPERBzdxy,     /*!< Second derivative of face-averaged Bz to xy-direction. */
      dp11dx,        /*!< Derivative of P_11 to x direction. */
      dp11dy,        /*!< Derivative of P_11 to x direction. */
      dp11dz,        /*!< Derivative of P_11 to x direction. */
      dp22dx,        /*!< Derivative of P_22 to y direction. */
      dp22dy,        /*!< Derivative of P_22 to y direction. */
      dp22dz,        /*!< Derivative of P_22 to y direction. */
      dp33dx,        /*!< Derivative of P_33 to z direction. */
      dp33dy,        /*!< Derivative of P_33 to z direction. */
      dp33dz,        /*!< Derivative of P_33 to z direction. */
      // End of insert for Hall term
      dVxdx,     /*!< Derivative of volume-averaged Vx to x-direction. */
      dVxdy,     /*!< Derivative of volume-averaged Vx to y-direction. */
      dVxdz,     /*!< Derivative of volume-averaged Vx to z-direction. */
      dVydx,     /*!< Derivative of volume-averaged Vy to x-direction. */
      dVydy,     /*!< Derivative of volume-averaged Vy to y-direction. */
      dVydz,     /*!< Derivative of volume-averaged Vy to z-direction. */
      dVzdx,     /*!< Derivative of volume-averaged Vz to x-direction. */
      dVzdy,     /*!< Derivative of volume-averaged Vz to y-direction. */
      dVzdz,     /*!< Derivative of volume-averaged Vz to z-direction. */
      N_SPATIAL_CELL_DERIVATIVES
   };
}

/*! The namespace bvolderivatives contains the indices to an array which stores the spatial
 * derivatives of the volume-averaged magnetic field, needed for Lorentz force. 
 * TODO: Vol values may be removed if background field is curlfree
 */
namespace bvolderivatives {
   enum {
      dBGBXVOLdy,
      dBGBXVOLdz,
      dBGBYVOLdx,
      dBGBYVOLdz,
      dBGBZVOLdx,
      dBGBZVOLdy,
      dPERBXVOLdy,
      dPERBXVOLdz,
      dPERBYVOLdx,
      dPERBYVOLdz,
      dPERBZVOLdx,
      dPERBZVOLdy,
      N_BVOL_DERIVATIVES
   };
}

/*! The namespace sysboundarytype contains the identification index of the boundary condition types applied to a cell,
 * it is stored in SpatialCell::sysBoundaryFlag and used by the BoundaryCondition class' functions to determine what type of BC to apply to a cell.
 * At least for the workings of vlasovmover_leveque.cpp the order of the first two entries should not be changed.
 */
namespace sysboundarytype {
   enum {
      DO_NOT_COMPUTE,   /*!< E.g. cells within the ionospheric outer radius should not be computed at all. */
      NOT_SYSBOUNDARY,  /*!< Cells within the simulation domain are not boundary cells. */
      IONOSPHERE,       /*!< Initially a perfectly conducting sphere. */
      OUTFLOW,          /*!< No fixed conditions on the fields and distribution function. */
      SET_MAXWELLIAN,   /*!< Set Maxwellian boundary condition, i.e. set fields and distribution function. */
      ANTISYMMETRIC,    /*!< System is antisymmetric wrt. to the boundary.*/
      PROJECT,         /*!< Simulated project's setCell and setCellBackgroundField functions are used 
                        * to set the boundary conditions.*/
      N_SYSBOUNDARY_CONDITIONS
   };
}

/*! Steps in Runge-Kutta methods */
enum {RK_ORDER1,   /*!< First order method, one step (and initialisation) */
RK_ORDER2_STEP1,   /*!< Two-step second order method, first step */
RK_ORDER2_STEP2    /*!< Two-step second order method, second step */
};

const uint WID = 4;         /*!< Number of cells per coordinate in a velocity block. Only a value of 4 supported by vectorized Leveque solver */
const uint WID2 = WID*WID;  /*!< Number of cells per 2D slab in a velocity block. */
const uint WID3 = WID2*WID; /*!< Number of cells in a velocity block. */

/*!
Get the cellindex in the velocity space block
*/
template<typename UINT> inline UINT cellIndex(const UINT& i,const UINT& j,const UINT& k) {
   return k*WID2 + j*WID + i;
}

const uint SIZE_VELBLOCK    = WID3; /*!< Number of cells in a velocity block. */

/*!
 * Name space for flags needed globally, such as the bailout flag.
 */
struct globalflags {
   static int bailingOut; /*!< Global flag raised to true if a run bailout (write restart if requested/set and stop the simulation peacefully) is needed. */
   static bool writeRestart; /*!< Global flag raised to true if a restart writing is needed (without bailout). NOTE: used only by MASTER_RANK in vlasiator.cpp. */
};

/*!
 * Name space for flags going into the project hook function.
 */
namespace hook {
   enum {
      END_OF_TIME_STEP
   };
}


// Natural constants
namespace physicalconstants {
   const Real EPS_0 = 8.85418782e-12; /*!< Permittivity of value, unit: C / (V m).*/
   const Real MU_0 = 1.25663706e-6;  /*!< Permeability of vacuo, units: (kg m) / (s^2 A^2).*/
   const Real K_B = 1.3806503e-23;   /*!< Boltzmann's constant, units: (kg m^2) / (s^2 K).*/
   const Real CHARGE = 1.60217653e-19; /*!< Elementary charge, units: C. */
   const Real MASS_ELECTRON = 9.10938188e-31; /**< Electron rest mass, units: kg.*/
   const Real MASS_PROTON = 1.67262158e-27; /*!< Proton rest mass, units: kg.*/
   const Real R_E = 6.3712e6; /*!< radius of the Earth, units: m. */
}

const std::vector<CellID>& getLocalCells();
void recalculateLocalCellsCache();

#endif
