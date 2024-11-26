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
   const std::string& message
);
void bailout(
   const bool condition,
   const char * const file,
   const int line
);
void bailout(
   const bool condition,
   const std::string& message,
   const char * const file,
   const int line
);

[[ noreturn ]] void abort_mpi(const std::string str, const int err_type = 0);

#define sqr(x) ((x)*(x))
#define pow2(x) sqr(x)
#define pow3(x) ((x)*(x)*(x))

#define MASTER_RANK 0


/* Maximum number of blocks in each dimension in velocity space. The
   size of velocity space defined in cfg can at maximum be this large
*/
#define MAX_BLOCKS_PER_DIM 256


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
 * field solver are relying on this particular ordering, even though the actual
 * fsgrid data layouts might be slightly different (see below).
 *
 * Note: RHOM and RHOQ are somewhat out-of-order for backwards compatibility
 * with pre-multipop tools.
 */
namespace CellParams {
   enum {
      XCRD,   /*!< x-coordinate of the bottom left corner.*/
      YCRD,   /*!< y-coordinate of the bottom left corner.*/
      ZCRD,   /*!< z-coordinate of the bottom left corner.*/
      // DX,DY,DZ have to be consecutive.
      DX,     /*!< Grid separation in x-coordinate.*/
      DY,     /*!< Grid separation in y-coordinate.*/
      DZ,     /*!< Grid separation in z-coordinate.*/
      RHOM,    /*!< Total mass density. Calculated by Vlasov propagator, used to propagate BX,BY,BZ.*/
      VX,  /*!< Vx. Calculated by Vlasov propagator, used to propagate BX,BY,BZ.*/
      VY,  /*!< Vy. Calculated by Vlasov propagator, used to propagate BX,BY,BZ.*/
      VZ,  /*!< Vz. Calculated by Vlasov propagator, used to propagate BX,BY,BZ.*/
      RHOQ,    /*!< Total charge density. Calculated by Vlasov propagator, used to propagate BX,BY,BZ.*/
      RHOM_DT2,    /*!< Total mass density. Calculated by Vlasov propagator, used to propagate BX,BY,BZ.*/
      VX_DT2,  /*!< Vx. Calculated by Vlasov propagator, used to propagate BX,BY,BZ.*/
      VY_DT2,  /*!< Vy. Calculated by Vlasov propagator, used to propagate BX,BY,BZ.*/
      VZ_DT2,  /*!< Vz. Calculated by Vlasov propagator, used to propagate BX,BY,BZ.*/
      RHOQ_DT2,    /*!< Total charge density. Calculated by Vlasov propagator, used to propagate BX,BY,BZ.*/
      BGBXVOL,   /*!< background magnetic field averaged over spatial cell.*/
      BGBYVOL,   /*!< background magnetic field averaged over spatial cell.*/
      BGBZVOL,   /*!< background magnetic field averaged over spatial cell.*/
      PERBXVOL,  /*!< perturbed magnetic field  PERBX averaged over spatial cell.*/
      PERBYVOL,  /*!< perturbed magnetic field  PERBY averaged over spatial cell.*/
      PERBZVOL,  /*!< perturbed magnetic field  PERBZ averaged over spatial cell.*/
      EXGRADPE,         /*!< Electron pressure gradient term x.*/
      EYGRADPE,         /*!< Electron pressure gradient term y.*/
      EZGRADPE,         /*!< Electron pressure gradient term z.*/
      RHOM_R,     /*!< RHO after propagation in ordinary space*/
      VX_R,   /*!< VX after propagation in ordinary space*/
      VY_R,   /*!< VY after propagation in ordinary space*/
      VZ_R,   /*!< VZ after propagation in ordinary space*/
      RHOQ_R,     /*!< RHOQ after propagation in ordinary space*/
      RHOM_V,     /*!< RHOM after propagation in velocity space*/
      VX_V,   /*!< VX after propagation in velocity space*/
      VY_V,   /*!< VY after propagation in velocity space*/
      VZ_V,   /*!< VZ after propagation in velocity space*/
      RHOQ_V,     /*!< RHOQ after propagation in velocity space*/
      P_11,     /*!< Pressure P_xx component, computed by Vlasov propagator. */
      P_22,     /*!< Pressure P_yy component, computed by Vlasov propagator. */
      P_33,     /*!< Pressure P_zz component, computed by Vlasov propagator. */
      P_23,     /*!< Pressure P_yz component, computed by Vlasov propagator. */
      P_13,     /*!< Pressure P_xz component, computed by Vlasov propagator. */
      P_12,     /*!< Pressure P_xy component, computed by Vlasov propagator. */
      P_11_DT2, /*!< Intermediate step value for RK2 time stepping in field solver. Computed from P_11_R and P_11_V. */
      P_22_DT2, /*!< Intermediate step value for RK2 time stepping in field solver. Computed from P_22_R and P_22_V. */
      P_33_DT2, /*!< Intermediate step value for RK2 time stepping in field solver. Computed from P_33_R and P_33_V. */
      P_23_DT2, /*!< Intermediate step value for RK2 time stepping in field solver. Computed from P_23_R and P_23_V. */
      P_13_DT2, /*!< Intermediate step value for RK2 time stepping in field solver. Computed from P_13_R and P_13_V. */
      P_12_DT2, /*!< Intermediate step value for RK2 time stepping in field solver. Computed from P_12_R and P_12_V. */
      P_11_R,   /*!< P_xx component after propagation in ordinary space */
      P_22_R,   /*!< P_yy component after propagation in ordinary space */
      P_33_R,   /*!< P_zz component after propagation in ordinary space */
      P_23_R,   /*!< P_yz component after propagation in ordinary space */
      P_13_R,   /*!< P_xz component after propagation in ordinary space */
      P_12_R,   /*!< P_xy component after propagation in ordinary space */
      P_11_V,   /*!< P_xx component after propagation in velocity space */
      P_22_V,   /*!< P_yy component after propagation in velocity space */
      P_33_V,   /*!< P_zz component after propagation in velocity space */
      P_23_V,   /*!< P_yz component after propagation in velocity space */
      P_12_V,   /*!< P_xy component after propagation in velocity space */
      P_13_V,   /*!< P_xz component after propagation in velocity space */
      EXVOL,    /*!< Volume electric field averaged over spatial cell, x-component.*/
      EYVOL,    /*!< Volume electric field averaged over spatial cell, y-component.*/
      EZVOL,    /*!< Volume electric field averaged over spatial cell, z-component.*/
      MAXVDT,             /*!< maximum timestep allowed in velocity space for this cell,
                           * this is the max allowed timestep over all particle species.*/
      MAXRDT,             /*!< maximum timestep allowed in ordinary space for this cell,
                           * this is the max allowed timestep over all particle species.*/
      MAXFDT,             /*!< maximum timestep allowed in ordinary space by fieldsolver for this cell**/
      LBWEIGHTCOUNTER,    /*!< Counter for storing compute time weights needed by the load balancing**/
      ISCELLSAVINGF,      /*!< Value telling whether a cell is saving its distribution function when partial f data is written out. */
      FSGRID_RANK, /*!< Rank of this cell in the FsGrid cartesian communicator */
      FSGRID_BOUNDARYTYPE, /*!< Boundary type of this cell, as stored in the fsGrid */
      CELLID, /*! < DCCRG cell index */
      REFINEMENT_LEVEL, /*! < Refinement level */
      CONNECTION, /*!< Magnetic connection. See TracingPointConnectionType for assigned values. */
      CONNECTION_FW_X, /*!< Endpoint x (forward-propagated) for the magnetic connection tracer*/
      CONNECTION_FW_Y, /*!< Endpoint y (forward-propagated) for the magnetic connection tracer*/
      CONNECTION_FW_Z, /*!< Endpoint z (forward-propagated) for the magnetic connection tracer*/
      CONNECTION_BW_X, /*!< Endpoint x (backward-propagated) for the magnetic connection tracer*/
      CONNECTION_BW_Y, /*!< Endpoint y (backward-propagated) for the magnetic connection tracer*/
      CONNECTION_BW_Z, /*!< Endpoint z (backward-propagated) for the magnetic connection tracer*/
      CURVATUREX, /*!< Magnetic field curvature vector x component */
      CURVATUREY, /*!< Magnetic field curvature vector y component */
      CURVATUREZ, /*!< Magnetic field curvature vector z component */
      FLUXROPE,   /*!< 0 for regular cells, 1 for cells considered to be in a flux rope */
      AMR_DRHO,
      AMR_DU,
      AMR_DPSQ,
      AMR_DBSQ,
      AMR_DB,
      AMR_ALPHA1,
      AMR_ALPHA2,
      P_ANISOTROPY,
      AMR_VORTICITY,
      RECENTLY_REFINED,
      BULKV_FORCING_X, /*! Externally forced drift velocity (ex. from the ionosphere) */
      BULKV_FORCING_Y, /*! Externally forced drift velocity (ex. from the ionosphere) */
      BULKV_FORCING_Z, /*! Externally forced drift velocity (ex. from the ionosphere) */
      N_SPATIAL_CELL_PARAMS
   };
}

/*! The namespace bvolderivatives contains the indices to an array which stores the spatial
 * derivatives of the volume-averaged magnetic field, needed for Lorentz force.
 */
namespace bvolderivatives {
   enum {
      dPERBXVOLdx, /*!< Derivative of perturbed volume-averaged Bx in x-direction. */
      dPERBXVOLdy, /*!< Derivative of perturbed volume-averaged Bx in y-direction. */
      dPERBXVOLdz, /*!< Derivative of perturbed volume-averaged Bx in z-direction. */
      dPERBYVOLdx, /*!< Derivative of perturbed volume-averaged By in x-direction. */
      dPERBYVOLdy, /*!< Derivative of perturbed volume-averaged By in y-direction. */
      dPERBYVOLdz, /*!< Derivative of perturbed volume-averaged By in z-direction. */
      dPERBZVOLdx, /*!< Derivative of perturbed volume-averaged Bz in x-direction. */
      dPERBZVOLdy, /*!< Derivative of perturbed volume-averaged Bz in y-direction. */
      dPERBZVOLdz, /*!< Derivative of perturbed volume-averaged Bz in z-direction. */
      N_BVOL_DERIVATIVES
   };
}

namespace vderivatives {
   // Essentially a copy from dmoments
   enum {
      dVxdx,     /*!< Derivative of volume-averaged Vx to x-direction. */
      dVxdy,     /*!< Derivative of volume-averaged Vx to y-direction. */
      dVxdz,     /*!< Derivative of volume-averaged Vx to z-direction. */
      dVydx,     /*!< Derivative of volume-averaged Vy to x-direction. */
      dVydy,     /*!< Derivative of volume-averaged Vy to y-direction. */
      dVydz,     /*!< Derivative of volume-averaged Vy to z-direction. */
      dVzdx,     /*!< Derivative of volume-averaged Vz to x-direction. */
      dVzdy,     /*!< Derivative of volume-averaged Vz to y-direction. */
      dVzdz,     /*!< Derivative of volume-averaged Vz to z-direction. */
      N_V_DERIVATIVES
   };
}

// FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
// FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBDt2Grid,
// FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
// FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EDt2Grid,
// FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
// FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
// FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
// FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsDt2Grid,
// FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
// FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
// FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
// FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
// FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,

/*! Namespace containing enums and structs for the various field solver grid instances
 *
 * Note that in some of these, the order of members differs from the cell
 * parameter fields (see above). So double-check before blindly copying data
 * back and forth.
 */
namespace fsgrids {
   enum bfield : int {
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
      RHOM, /*!< Overall mass density. Calculated by Vlasov propagator, used to propagate fields.*/
      RHOQ, /*!< Overall charge density. Calculated by Vlasov propagator, used to propagate fields.*/
      VX,   /*!< Vx. Calculated by Vlasov propagator, used to propagate fields.*/
      VY,   /*!< Vy. Calculated by Vlasov propagator, used to propagate fields.*/
      VZ,   /*!< Vz. Calculated by Vlasov propagator, used to propagate fields.*/
      P_11, /*!< Pressure P_xx component, computed by Vlasov propagator. */
      P_22, /*!< Pressure P_yy component, computed by Vlasov propagator. */
      P_33, /*!< Pressure P_zz component, computed by Vlasov propagator. */
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
   };

   enum dmoments {
      drhomdx,    /*!< Derivative of mass density to x-direction. */
      drhomdy,    /*!< Derivative of mass density to y-direction. */
      drhomdz,    /*!< Derivative of mass density to z-direction. */
      drhoqdx,    /*!< Derivative of charge density to x-direction. */
      drhoqdy,    /*!< Derivative of charge density to y-direction. */
      drhoqdz,    /*!< Derivative of charge density to z-direction. */
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
      dPedx,    /*!< Derivative of electron pressure to x-direction. */
      dPedy,    /*!< Derivative of electron pressure to y-direction. */
      dPedz,    /*!< Derivative of electron pressure to z-direction. */
      N_DMOMENTS
   };

   // NOTE This contains the BGB derivatives as they do not change either
   enum bgbfield : int  {
      BGBX,   /*!< Background magnetic field x-component, averaged over cell x-face.*/
      BGBY,   /*!< Background magnetic field y-component, averaged over cell y-face.*/
      BGBZ,   /*!< Background magnetic field z-component, averaged over cell z-face.*/
      BGBXVOL,   /*!< background magnetic field x-component averaged over spatial cell.*/
      BGBYVOL,   /*!< background magnetic field y-component averaged over spatial cell.*/
      BGBZVOL,   /*!< background magnetic field z-component averaged over spatial cell.*/
      BGBXVDCORR, /*!< correction term for background magnetic field, used by vector dipole.*/
      BGBYVDCORR, /*!< correction term for background magnetic field, used by vector dipole.*/
      BGBZVDCORR, /*!< correction term for background magnetic field, used by vector dipole.*/
      dBGBxdy,     /*!< Derivative of background face-averaged Bx in y-direction. */
      dBGBxdz,     /*!< Derivative of background face-averaged Bx in z-direction. */
      dBGBydx,     /*!< Derivative of background face-averaged By in x-direction. */
      dBGBydz,     /*!< Derivative of background face-averaged By in z-direction. */
      dBGBzdx,     /*!< Derivative of background face-averaged Bz in x-direction. */
      dBGBzdy,     /*!< Derivative of background face-averaged Bz in y-direction. */
      dBGBXVOLdx,  /*!< Derivative of background volume-averaged Bx in x-direction. */
      dBGBXVOLdy,  /*!< Derivative of background volume-averaged Bx in y-direction. */
      dBGBXVOLdz,  /*!< Derivative of background volume-averaged Bx in z-direction. */
      dBGBYVOLdx,  /*!< Derivative of background volume-averaged By in x-direction. */
      dBGBYVOLdy,  /*!< Derivative of background volume-averaged By in y-direction. */
      dBGBYVOLdz,  /*!< Derivative of background volume-averaged By in z-direction. */
      dBGBZVOLdx,  /*!< Derivative of background volume-averaged Bz in x-direction. */
      dBGBZVOLdy,  /*!< Derivative of background volume-averaged Bz in y-direction. */
      dBGBZVOLdz,  /*!< Derivative of background volume-averaged Bz in z-direction. */
      N_BGB
   };

   // NOTE This contains the PERBVOL derivatives
   enum volfields {
      PERBXVOL,  /*!< perturbed magnetic field  PERBX averaged over spatial cell.*/
      PERBYVOL,  /*!< perturbed magnetic field  PERBY averaged over spatial cell.*/
      PERBZVOL,  /*!< perturbed magnetic field  PERBZ averaged over spatial cell.*/
      dPERBXVOLdx, /*!< Derivative of perturbed volume-averaged Bx in x-direction. */
      dPERBXVOLdy, /*!< Derivative of perturbed volume-averaged Bx in y-direction. */
      dPERBXVOLdz, /*!< Derivative of perturbed volume-averaged Bx in z-direction. */
      dPERBYVOLdx, /*!< Derivative of perturbed volume-averaged By in x-direction. */
      dPERBYVOLdy, /*!< Derivative of perturbed volume-averaged By in y-direction. */
      dPERBYVOLdz, /*!< Derivative of perturbed volume-averaged By in z-direction. */
      dPERBZVOLdx, /*!< Derivative of perturbed volume-averaged Bz in x-direction. */
      dPERBZVOLdy, /*!< Derivative of perturbed volume-averaged Bz in y-direction. */
      dPERBZVOLdz, /*!< Derivative of perturbed volume-averaged Bz in z-direction. */
      EXVOL,   /*!< volume-averaged electric field x component */
      EYVOL,   /*!< volume-averaged electric field y component */
      EZVOL,   /*!< volume-averaged electric field z component */
      CURVATUREX, /*!< Magnetic field curvature vector x component, grid-glued to DCCRG */
      CURVATUREY, /*!< Magnetic field curvature vector y component, grid-glued to DCCRG */
      CURVATUREZ, /*!< Magnetic field curvature vector z component, grid-glued to DCCRG */
      N_VOL
   };

   struct technical {
      uint sysBoundaryFlag;  /*!< System boundary flags. */
      int sysBoundaryLayer; /*!< System boundary layer index. */
      Real maxFsDt;         /*!< maximum timestep allowed in ordinary space by fieldsolver for this cell**/
      int fsGridRank;       /*!< Rank in the fsGrids cartesian coordinator */
      uint SOLVE;           /*!< Bit mask to determine whether a given cell should solve E or B components. */
      int refLevel;         /*!<AMR Refinement Level*/
   };

}

// Ionosphere node parameters
enum ionosphereParameters {
  SOURCE,    /*!< Field aligned current source term (Ampere). Note: this is current, *not* density. */
  SIGMA,   SIGMA12, SIGMA13,
  SIGMA21, SIGMA22, SIGMA23, /*!< Overall conductivity tensor */
  SIGMA31, SIGMA32, SIGMA33,

  SIGMAP,   /*!< Scalar Pedersen conductivity */
  SIGMAH,   /*!< Scalar Hall conductivity */
  SIGMAPARALLEL, /*!< Scalar parallel conductivity */
  PRECIP,   /*!< Precipitation */
  RHON,     /*!< Downmapped magnetospheric plasma number density */
  TEMPERATURE, /*!< Downmapped electron temperature */
  NODE_BX,NODE_BY,NODE_BZ, /*!< Magnetic field at the node */
  UPMAPPED_BX,UPMAPPED_BY,UPMAPPED_BZ, /*!< Magnetic field at the upper and of the mapping fieldline */
  SOLUTION, /*!< Currently considered solution potential */
  BEST_SOLUTION, /*!< Best solution found so far */
  RESIDUAL, /*!< Residual of the current solution */
  RRESIDUAL,
  ZPARAM, ZZPARAM,
  PPPARAM, PPARAM,
  N_IONOSPHERE_PARAMETERS
};

/*! The namespace sysboundarytype contains the identification index of the boundary condition types applied to a cell,
 * it is stored in SpatialCell::sysBoundaryFlag and used by the BoundaryCondition class' functions to determine what type of BC to apply to a cell.
 * At least for the workings of vlasovmover_leveque.cpp the order of the first two entries should not be changed.
 */
namespace sysboundarytype {
   enum {
      DO_NOT_COMPUTE,   /*!< E.g. cells within the ionospheric outer radius should not be computed at all. */
      NOT_SYSBOUNDARY,  /*!< Cells within the simulation domain are not boundary cells. */
      IONOSPHERE,       /*!< Ionospheric current model. */
      OUTFLOW,          /*!< No fixed conditions on the fields and distribution function. */
      MAXWELLIAN,       /*!< Set Maxwellian boundary condition, i.e. set fields and distribution function. */
      COPYSPHERE,       /*!< A sphere with copy-condition for perturbed B as the simple inner boundary */
      OUTER_BOUNDARY_PADDING, /*!< These cells only occur on FSGrid, where boundaries are not at the highest refinement level */
      N_SYSBOUNDARY_CONDITIONS
   };
}

namespace compute {
   const uint BX = (1 << 0); // 1
   const uint BY = (1 << 1); // 2
   const uint BZ = (1 << 2); // 4
   const uint EX = (1 << 3); // 8
   const uint EY = (1 << 4); // 16
   const uint EZ = (1 << 5); // 32
}

/*! Steps in Runge-Kutta methods */
enum {RK_ORDER1,   /*!< First order method, one step (and initialisation) */
RK_ORDER2_STEP1,   /*!< Two-step second order method, first step */
RK_ORDER2_STEP2    /*!< Two-step second order method, second step */
};

const int WID = 4;         /*!< Number of cells per coordinate in a velocity block. Only a value of 4 supported by vectorized Leveque solver */
const int WID2 = WID*WID;  /*!< Number of cells per 2D slab in a velocity block. */
const int WID3 = WID2*WID; /*!< Number of cells in a velocity block. */

/*!
Get the cellindex in the velocity space block
*/
template<typename INT> inline INT cellIndex(const INT& i,const INT& j,const INT& k) {
   return k*WID2 + j*WID + i;
}

const int SIZE_VELBLOCK    = WID3; /*!< Number of cells in a velocity block. */

template<typename T> inline int sign(const T& value) {
   const T ZERO = 0.0;
   if (value < ZERO) return -1;
   return 1;
}

/*!
 * Name space for flags needed globally, such as the bailout flag.
 */
struct globalflags {
   static int bailingOut; /*!< Global flag raised to true if a run bailout (write restart if requested/set and stop the simulation peacefully) is needed. */
   static bool writeRestart; /*!< Global flag raised to true if a restart writing is needed (without bailout). NOTE: used only by MASTER_RANK in vlasiator.cpp. */
   static bool balanceLoad; /*!< Global flag raised to true if a load balancing is needed. NOTE: used only by MASTER_RANK in vlasiator.cpp. */
   static bool doRefine; /*!< Global flag raised to true if a re-refine is needed. NOTE: used only by MASTER_RANK in vlasiator.cpp. */
   static bool ionosphereJustSolved; /*!< Flag used to notify that the ionosphere has been freshly solved, used to check whether the Vlasov boundary/bulk forcing need updating. */
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

#endif
