/*
This file is part of Vlasiator.

Copyright 2010, 2011, 2012 Finnish Meteorological Institute

Vlasiator is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3
as published by the Free Software Foundation.

Vlasiator is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef COMMON_H
#define COMMON_H

#include <limits>
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

/** A namespace for storing indices into an array which contains 
 * neighbour list for each spatial cell. These indices refer to 
 * the CPU memory, i.e. the device does not use these.
 */
namespace NbrsSpa {
   const uint INNER = 0;      /**< The cell is an inner cell, i.e. all its neighbours are located on the same computation node.*/
   const uint X_NEG_BND = (1 << 0);  /**< The cell is a boundary cell in -x direction.*/
   const uint X_POS_BND = (1 << 1);  /**< The cell is a boundary cell in +x direction.*/
   const uint Y_NEG_BND = (1 << 2);  /**< The cell is a boundary cell in -y direction.*/
   const uint Y_POS_BND = (1 << 3);  /**< The cell is a boundary cell in +y direction.*/
   const uint Z_NEG_BND = (1 << 4); /**< The cell is a boundary cell in -z direction.*/
   const uint Z_POS_BND = (1 << 5); /**< The cell is a boundary cell in +z direction.*/
   
   enum {
      STATE, /**< Contains the neighbour information of this cell, i.e. whether it is an inner cell or a boundary cell in one or more coordinate directions.*/
      MYIND, /**< The index of this cell.*/
      X1NEG,  /**< The index of the -x neighbouring block, distance 1.*/
      Y1NEG,  /**< The index of the -y neighbouring block, distance 1.*/
      Z1NEG,  /**< The index of the -z neighbouring block, distance 1.*/
      X1POS,  /**< The index of the +x neighbouring block, distance 1.*/
      Y1POS,  /**< The index of the +y neighbouring block, distance 1.*/
      Z1POS,  /**< The index of the +z neighbouring block, distance 1.*/
      X2NEG,  /**< The index of the -x neighbouring block, distance 1.*/
      Y2NEG,  /**< The index of the -y neighbouring block, distance 1.*/
      Z2NEG,  /**< The index of the -z neighbouring block, distance 1.*/
      X2POS,  /**< The index of the +x neighbouring block, distance 1.*/
      Y2POS,  /**< The index of the +y neighbouring block, distance 1.*/
      Z2POS   /**< The index of the +z neighbouring block, distance 1.*/
   };
}


/** A namespace for storing indices into an array which contains 
 * the physical parameters of each velocity block.*/
namespace BlockParams {
   enum {
      //Q_PER_M, /**< The charge-to-mass ratio of the particle species. DEPRECATED: GOING TO BE REMOVED!.*/
      VXCRD,   /**< vx-coordinate of the bottom left corner of the block.*/
      VYCRD,   /**< vy-coordinate of the bottom left corner of the block.*/
      VZCRD,   /**< vz-coordinate of the bottom left corner of the block.*/
      DVX,     /**< Grid separation in vx-coordinate for the block.*/
      DVY,     /**< Grid separation in vy-coordinate for the block.*/
      DVZ,     /**< Grid separation in vz-coordinate for the block.*/
      N_VELOCITY_BLOCK_PARAMS
   };
}

/** A namespace for storing indices into an array which contains the 
 * physical parameters of each spatial cell. Do not change the order 
 * of variables unless you know what you are doing - MPI transfers in 
 * field solver are optimised for this particular ordering.
 */
namespace CellParams {
   enum {
      XCRD,    /**< x-coordinate of the bottom left corner.*/
      YCRD,    /**< y-coordinate of the bottom left corner.*/
      ZCRD,    /**< z-coordinate of the bottom left corner.*/
      DX,      /**< Grid separation in x-coordinate.*/
      DY,      /**< Grid separation in y-coordinate.*/
      DZ,      /**< Grid separation in z-coordinate.*/
      EX,      /**< Electric field x-component, averaged over cell edge. Used to propagate BX,BY,BZ.*/
      EY,      /**< Electric field y-component, averaged over cell edge. Used to propagate BX,BY,BZ.*/
      EZ,      /**< Electric field z-component, averaged over cell edge. Used to propagate BX,BY,BZ.*/
      BX,      /**< Magnetic field x-component, averaged over cell x-face. Propagated by field solver.*/
      BY,      /**< Magnetic field y-component, averaged over cell y-face. Propagated by field solver.*/
      BZ,      /**< Magnetic field z-component, averaged over cell z-face. Propagated by field solver.*/
      RHO,     /**< Number density. Calculated by Vlasov propagator, used to propagate BX,BY,BZ.*/
      RHOVX,   /**< x-component of number density times Vx. Calculated by Vlasov propagator, used to propagate BX,BY,BZ.*/
      RHOVY,   /**< y-component of number density times Vy. Calculated by Vlasov propagator, used to propagate BX,BY,BZ.*/
      RHOVZ,   /**< z-component of number density times Vz. Calculated by Vlasov propagator, used to propagate BX,BY,BZ.*/
      BXFACEX, /**< Bx averaged over x-face. Used to propagate distribution function.*/
      BYFACEX, /**< By averaged over x-face. Used to propagate distribution function.*/
      BZFACEX, /**< Bz averaged over x-face. Used to propagate distribution function.*/
      BXFACEY, /**< Bx averaged over y-face. Used to propagate distribution function.*/
      BYFACEY, /**< By averaged over y-face. Used to propagate distribution function.*/
      BZFACEY, /**< Bz averaged over y-face. Used to propagate distribution function.*/
      BXFACEZ, /**< Bx averaged over z-face. Used to propagate distribution function.*/
      BYFACEZ, /**< By averaged over z-face. Used to propagate distribution function.*/
      BZFACEZ, /**< Bz averaged over z-face. Used to propagate distribution function.*/
      EXFACEX, /**< Ex averaged over x-face. Used to propagate distribution function.*/
      EYFACEX, /**< Ey averaged over x-face. Used to propagate distribution function.*/
      EZFACEX, /**< Ez averaged over x-face. Used to propagate distribution function.*/
      EXFACEY, /**< Ex averaged over y-face. Used to propagate distribution function.*/
      EYFACEY, /**< Ey averaged over y-face. Used to propagate distribution function.*/
      EZFACEY, /**< Ez averaged over y-face. Used to propagate distribution function.*/
      EXFACEZ, /**< Ex averaged over z-face. Used to propagate distribution function.*/
      EYFACEZ, /**< Ey averaged over z-face. Used to propagate distribution function.*/
      EZFACEZ, /**< Ez averaged over z-face. Used to propagate distribution function.*/
      BXVOL,   /**< Bx averaged over spatial cell.*/
      BYVOL,   /**< By averaged over spatial cell.*/
      BZVOL,   /**< Bz averaged over spatial cell.*/
      EXVOL,   /**< Ex averaged over spatial cell.*/
      EYVOL,   /**< Ey averaged over spatial cell.*/
      EZVOL,   /**< Ez averaged over spatial cell.*/
      N_SPATIAL_CELL_PARAMS
   };
}

/** Namespace fieldsolver contains indices into arrays which store 
 * variables required by the field solver. These quantities are derivatives 
 * of variables described in namespace CellParams.
 * Do not change the order of variables unless you know what you are doing: 
 * in several places the size of cpu_derivatives array in cell_spatial is calculated 
 * as fieldsolver::dVzdz+1.
 */
namespace fieldsolver {
   enum {
      drhodx,    /**< Derivative of volume-averaged number density to x-direction. */
      drhody,    /**< Derivative of volume-averaged number density to y-direction. */
      drhodz,    /**< Derivative of volume-averaged number density to z-direction. */
      dBxdy,     /**< Derivative of face-averaged Bx to y-direction. */
      dBxdz,     /**< Derivative of face-averaged Bx to z-direction. */
      dBydx,     /**< Derivative of face-averaged By to x-direction. */
      dBydz,     /**< Derivative of face-averaged By to z-direction. */
      dBzdx,     /**< Derivative of face-averaged Bz to x-direction. */
      dBzdy,     /**< Derivative of face-averaged Bz to y-direction. */
      dVxdx,     /**< Derivative of volume-averaged Vx to x-direction. */
      dVxdy,     /**< Derivative of volume-averaged Vx to y-direction. */
      dVxdz,     /**< Derivative of volume-averaged Vx to z-direction. */
      dVydx,     /**< Derivative of volume-averaged Vy to x-direction. */
      dVydy,     /**< Derivative of volume-averaged Vy to y-direction. */
      dVydz,     /**< Derivative of volume-averaged Vy to z-direction. */
      dVzdx,     /**< Derivative of volume-averaged Vz to x-direction. */
      dVzdy,     /**< Derivative of volume-averaged Vz to y-direction. */
      dVzdz,     /**< Derivative of volume-averaged Vz to z-direction. */
      N_SPATIAL_CELL_DERIVATIVES
   };
}




// TODO: change into a template parameter
const uint WID = 4;         /**< Number of cells per coordinate in a velocity block. */
const uint WID2 = WID*WID;
const uint WID3 = WID2*WID; 

//const uint SIZE_NBRS_VEL    = 8;    /**< The size of velocity grid neighbour list per velocity block. */
//const uint SIZE_NBRS_VEL    = 28;    /**< The size of velocity grid neighbour list per velocity block. */
const uint SIZE_NBRS_SPA    = 31;   /**< The size of spatial grid neighbour list per spatial cell. */
const uint SIZE_VELBLOCK    = WID3; /**< Number of cells in a velocity block. */
const uint SIZE_BLOCKPARAMS = 6;    /**< Number of parameters per velocity block. */

const uint SIZE_BOUND       = WID3;
const uint SIZE_BDERI       = WID2;
const uint SIZE_BFLUX       = WID2;
const uint SIZE_DERIV       = WID3;
const uint SIZE_FLUXS       = WID3;

// Natural constants
namespace physicalconstants {
   const Real MU_0 = 1.25663706e-6;  /**< Permeability of vacuo, unit: (kg m) / (s^2 A^2).*/
   const Real K_B = 1.3806503e-23;   /**< Boltzmann's constant, unit: (kg m^2) / (s^2 K).*/
   
   const Real MASS_ELECTRON = 1.60238188e-31; /**< Electron rest mass.*/
   const Real MASS_PROTON = 1.67262158e-27; /**< Proton rest mass.*/
}

#endif
