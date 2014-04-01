/*
This file is part of Vlasiator.

Copyright 2010, 2011, 2012, 2013, 2014 Finnish Meteorological Institute

*/

#ifndef COMMON_H
#define COMMON_H

#include <limits>
#include <string>
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

#define BAILOUT(condition) \
   if ((condition) && (globalflags::bailingOut == 0)) { \
      int myRank; \
      MPI_Comm_rank(MPI_COMM_WORLD,&myRank); \
      std::cerr << "Process " << myRank << " bailing out at " << __FILE__ << ":" << __LINE__ << "." << std::endl; \
      globalflags::bailingOut = 1; \
   }

#define sqr(x) ((x)*(x))
#define pow2(x) sqr(x)
#define pow3(x) ((x)*(x)*(x))

#define MASTER_RANK 0

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
      //Q_PER_M, /*!< The charge-to-mass ratio of the particle species. DEPRECATED: GOING TO BE REMOVED!.*/
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
      EX_DT2,    /*!< Intermediate step value for RK2 time stepping in field solver.*/
      EY_DT2,    /*!< Intermediate step value for RK2 time stepping in field solver.*/
      EZ_DT2,    /*!< Intermediate step value for RK2 time stepping in field solver.*/
      PERBX_DT2, /*!< Intermediate step value for PERBX for RK2 time stepping in field solver.*/
      PERBY_DT2, /*!< Intermediate step value for PERBY for RK2 time stepping in field solver.*/
      PERBZ_DT2, /*!< Intermediate step value for PERBZ for RK2 time stepping in field solver.*/
      RHO_DT2,   /*!< Intermediate step value for RK2 time stepping in field solver. Computed from RHO_R and RHO_V*/
      RHOVX_DT2, /*!< Intermediate step value for RK2 time stepping in field solver. Computed from RHOVX_R and RHOVX_V*/
      RHOVY_DT2, /*!< Intermediate step value for RK2 time stepping in field solver. Computed from RHOVY_R and RHOVY_V*/
      RHOVZ_DT2, /*!< Intermediate step value for RK2 time stepping in field solver. Computed from RHOVZ_R and RHOVZ_V*/
      BGBXVOL,   /*!< background magnetic field averaged over spatial cell.*/
      BGBYVOL,   /*!< background magnetic field averaged over spatial cell.*/
      BGBZVOL,   /*!< background magnetic field averaged over spatial cell.*/
      PERBXVOL,  /*!< perturbed magnetic field  PERBX averaged over spatial cell.*/
      PERBYVOL,  /*!< perturbed magnetic field  PERBY averaged over spatial cell.*/
      PERBZVOL,  /*!< perturbed magnetic field  PERBZ averaged over spatial cell.*/
      EXVOL,     /*!< Ex averaged over spatial cell.*/
      EYVOL,     /*!< Ey averaged over spatial cell.*/
      EZVOL,     /*!< Ez averaged over spatial cell.*/
      RHO_R,     /*!< RHO after propagation in ordinary space*/
      RHOVX_R,   /*!< RHOVX after propagation in ordinary space*/
      RHOVY_R,   /*!< RHOVX after propagation in ordinary space*/
      RHOVZ_R,   /*!< RHOVX after propagation in ordinary space*/
      RHO_V,     /*!< RHO after propagation in velocity space*/
      RHOVX_V,   /*!< RHOVX after propagation in velocity space*/
      RHOVY_V,   /*!< RHOVX after propagation in velocity space*/
      RHOVZ_V,   /*!< RHOVX after propagation in velocity space*/
      RHOLOSSADJUST,      /*!< Counter for massloss from the destroying blocks in blockadjustment*/
      RHOLOSSVELBOUNDARY, /*!< Counter for massloss through outflow boundaries in velocity space*/
      MAXVDT,             /*!< maximum timestep allowed in velocity space for this cell**/
      MAXRDT,             /*!< maximum timestep allowed in ordinary space for this cell **/
      MAXFDT,             /*!< maximum timestep allowed in ordinary space by fieldsolver for this cell**/
      LBWEIGHTCOUNTER,    /*!< Counter for storing compute time weights needed by the load balancing**/
      ISCELLSAVINGF,      /*!< Value telling whether a cell is saving its distribution function when partial f data is written out. */
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
      drhodx,    /*!< Derivative of volume-averaged number density to x-direction. */
      drhody,    /*!< Derivative of volume-averaged number density to y-direction. */
      drhodz,    /*!< Derivative of volume-averaged number density to z-direction. */
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
      DO_NOT_COMPUTE, /*!< E.g. cells within the ionospheric outer radius should not be computed at all. */
      NOT_SYSBOUNDARY, /*!< Cells within the simulation domain are not boundary cells. */
      IONOSPHERE, /*!< Initially a perfectly conducting sphere. */
      OUTFLOW, /*!< No fixed conditions on the fields and distribution function. */
      SET_MAXWELLIAN, /*!< Set Maxwellian boundary condition, i.e. set fields and distribution function. */
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
   static int bailingOut; /*!< Global flag raised to true if a run bailout (write restart and stop the simulation peacefully) is needed. */
};

// Natural constants
namespace physicalconstants {
   const Real MU_0 = 1.25663706e-6;  /*!< Permeability of vacuo, unit: (kg m) / (s^2 A^2).*/
   const Real K_B = 1.3806503e-23;   /*!< Boltzmann's constant, unit: (kg m^2) / (s^2 K).*/
   
   const Real MASS_ELECTRON = 1.60238188e-31; /*!< Electron rest mass.*/
   const Real MASS_PROTON = 1.67262158e-27; /*!< Proton rest mass.*/
   const Real R_E = 6.3712e6; /*!< radius of the Earth. */
}




#endif


