#ifndef COMMON_H
#define COMMON_H

#include <limits>
#include "definitions.h"

/** A namespace for storing indices into an array which contains 
 * neighbour list for each spatial cell. These indices refer to 
 * the CPU memory, i.e. the device does not use these.
 */
namespace NbrsSpa {
   cuint INNER = 0;      /**< The cell is an inner cell, i.e. all its neighbours are located on the same computation node.*/
   cuint X_NEG_BND = (1 << 0);  /**< The cell is a boundary cell in -x direction.*/
   cuint X_POS_BND = (1 << 1);  /**< The cell is a boundary cell in +x direction.*/
   cuint Y_NEG_BND = (1 << 2);  /**< The cell is a boundary cell in -y direction.*/
   cuint Y_POS_BND = (1 << 3);  /**< The cell is a boundary cell in +y direction.*/
   cuint Z_NEG_BND = (1 << 4); /**< The cell is a boundary cell in -z direction.*/
   cuint Z_POS_BND = (1 << 5); /**< The cell is a boundary cell in +z direction.*/
   
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

/** A namespace for storing indices into an array containing neighbour information 
 * of velocity grid blocks. These indices are used by the device.
 */
namespace NbrsVel {
//   cuint INNER = 0;         /**< The block is an inner block, i.e. all its neighbours are stored on the same computation node.*/
//   cuint X_NEG_BND = (1 << 0);     /**< The block is a boundary block in -x direction.*/
//   cuint X_POS_BND = (1 << 1);     /**< The block is a boundary block in +x direction.*/
//   cuint Y_NEG_BND = (1 << 2);     /**< The block is a boundary block in -y direction.*/
//   cuint Y_POS_BND = (1 << 3);     /**< The block is a boundary block in +y direction.*/
//   cuint Z_NEG_BND = (1 << 4);    /**< The block is a boundary block in -z direction.*/
//   cuint Z_POS_BND = (1 << 5);    /**< The block is a boundary block in +z direction.*/
//   cuint VX_NEG_BND = (1 << 6);   /**< The block is a boundary block in -vx direction.*/
//   cuint VX_POS_BND = (1 << 7);  /**< The block is a boundary block in +vx direction.*/
//   cuint VY_NEG_BND = (1 << 8);  /**< The block is a boundary block in -vy direction.*/
//   cuint VY_POS_BND = (1 << 9);  /**< The block is a boundary block in +vy direction.*/
//   cuint VZ_NEG_BND = (1 << 10); /**< The block is a boundary block in -vz direction.*/
//   cuint VZ_POS_BND = (1 << 11); /**< The block is a boundary block in +vz direction.*/
   
//   enum {
//      STATE, /**< Contains the neighbour information bits of a velocity block.*/
//      MYIND, /**< The index of the block.*/   
//      VXNEG, /**< The index of -vx neighbour.*/
//      VYNEG, /**< The index of -vy neighbour.*/
//      VZNEG, /**< The index of -vz neighbour.*/
//      VXPOS, /**< The index of +vx neighbour.*/
//      VYPOS, /**< The index of +vy neighbour.*/
//      VZPOS  /**< The index of +vz neighbour.*/
//   };

   cuint XM1_YM1_ZM1 = 0;  /**< Index of (x-1,y-1,z-1) neighbour.*/
   cuint XCC_YM1_ZM1 = 1;  /**< Index of (x  ,y-1,z-1) neighbour.*/
   cuint XP1_YM1_ZM1 = 2;  /**< Index of (x+1,y-1,z-1) neighbour.*/
   cuint XM1_YCC_ZM1 = 3;  /**< Index of (x-1,y  ,z-1) neighbour.*/
   cuint XCC_YCC_ZM1 = 4;  /**< Index of (x  ,y  ,z-1) neighbour.*/
   cuint XP1_YCC_ZM1 = 5;  /**< Index of (x+1,y  ,z-1) neighbour.*/
   cuint XM1_YP1_ZM1 = 6;  /**< Index of (x-1,y+1,z-1) neighbour.*/
   cuint XCC_YP1_ZM1 = 7;  /**< Index of (x  ,y+1,z-1) neighbour.*/
   cuint XP1_YP1_ZM1 = 8;  /**< Index of (x+1,y+1,z-1) neighbour.*/   
   cuint XM1_YM1_ZCC = 9;  /**< Index of (x-1,y-1,z  ) neighbour.*/
   cuint XCC_YM1_ZCC = 10; /**< Index of (x  ,y-1,z  ) neighbour.*/
   cuint XP1_YM1_ZCC = 11; /**< Index of (x+1,y-1,z  ) neighbour.*/
   cuint XM1_YCC_ZCC = 12; /**< Index of (x-1,y  ,z  ) neighbour.*/
   cuint XCC_YCC_ZCC = 13; /**< Index of (x  ,y  ,z  ) neighbour.*/
   cuint XP1_YCC_ZCC = 14; /**< Index of (x+1,y  ,z  ) neighbour.*/
   cuint XM1_YP1_ZCC = 15; /**< Index of (x-1,y+1,z  ) neighbour.*/
   cuint XCC_YP1_ZCC = 16; /**< Index of (x  ,y+1,z  ) neighbour.*/
   cuint XP1_YP1_ZCC = 17; /**< Index of (x+1,y+1,z  ) neighbour.*/   
   cuint XM1_YM1_ZP1 = 18; /**< Index of (x-1,y-1,z+1) neighbour.*/
   cuint XCC_YM1_ZP1 = 19; /**< Index of (x  ,y-1,z+1) neighbour.*/
   cuint XP1_YM1_ZP1 = 20; /**< Index of (x+1,y-1,z+1) neighbour.*/
   cuint XM1_YCC_ZP1 = 21; /**< Index of (x-1,y  ,z+1) neighbour.*/
   cuint XCC_YCC_ZP1 = 22; /**< Index of (x  ,y  ,z+1) neighbour.*/
   cuint XP1_YCC_ZP1 = 23; /**< Index of (x+1,y  ,z+1) neighbour.*/
   cuint XM1_YP1_ZP1 = 24; /**< Index of (x-1,y+1,z+1) neighbour.*/
   cuint XCC_YP1_ZP1 = 25; /**< Index of (x  ,y+1,z+1) neighbour.*/
   cuint XP1_YP1_ZP1 = 26; /**< Index of (x+1,y+1,z+1) neighbour.*/
   
   cuint NON_EXISTING = std::numeric_limits<uint>::max(); /**< Invalid block ID, indicating that the block does not exist.*/
   
   cuint NBRFLAGS = 27; /**< Flags for existing neighbours.*/
   cuint MYIND    = 13; /**< Index of the block. Required for KT solver.*/
   cuint VXNEG    = 12; /**< Index of -vx neighbour. Required for KT solver.*/
   cuint VYNEG    = 10; /**< Index of -vy neighbour. Required for KT solver.*/
   cuint VZNEG    = 4;  /**< Index of -vz neighbour. Required for KT solver.*/
   cuint VXPOS    = 14; /**< Index of +vx neighbour. Required for KT solver.*/
   cuint VYPOS    = 16; /**< Index of +vy neighbour. Required for KT solver.*/
   cuint VZPOS    = 22; /**< Index of +vz neighbour. Required for KT solver.*/
   
   cuint VX_NEG_BND = (1 << VXNEG);
   cuint VX_POS_BND = (1 << VXPOS);
   cuint VY_NEG_BND = (1 << VYNEG);
   cuint VY_POS_BND = (1 << VYPOS);
   cuint VZ_NEG_BND = (1 << VZNEG);
   cuint VZ_POS_BND = (1 << VZPOS);
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
      DVZ      /**< Grid separation in vz-coordinate for the block.*/
   };
}

/** A namespace for storing indices into an array which contains the 
 * physical parameters of each spatial cell.*/
namespace CellParams {
   enum {
      XCRD,  /**< x-coordinate of the bottom left corner.*/
      YCRD,  /**< y-coordinate of the bottom left corner.*/
      ZCRD,  /**< z-coordinate of the bottom left corner.*/
      DX,    /**< Grid separation in x-coordinate.*/
      DY,    /**< Grid separation in y-coordinate.*/
      DZ,    /**< Grid separation in z-coordinate.*/
      EX,    /**< Electric field x-component.*/
      EY,    /**< Electric field y-component.*/
      EZ,    /**< Electric field z-component.*/
      BX,    /**< Magnetic field x-component.*/
      BY,    /**< Magnetic field y-component.*/
      BZ,    /**< Magnetic field z-component.*/
      RHO,   /**< Density.*/
      RHOVX, /**< x-component of the momentum density.*/
      RHOVY, /**< y-component of the momentum density.*/
      RHOVZ  /**< z-component of the momentum density.*/
   };
}

cuint WID = 4;         /**< Number of cells per coordinate in a velocity block. */
cuint WID2 = WID*WID;
cuint WID3 = WID2*WID; 

cuint SIZE_CELLPARAMS  = 16;   /**< The number of parameters for one spatial cell. */
//cuint SIZE_NBRS_VEL    = 8;    /**< The size of velocity grid neighbour list per velocity block. */
cuint SIZE_NBRS_VEL    = 28;    /**< The size of velocity grid neighbour list per velocity block. */
cuint SIZE_NBRS_SPA    = 31;   /**< The size of spatial grid neighbour list per spatial cell. */
cuint SIZE_VELBLOCK    = WID3; /**< Number of cells in a velocity block. */
cuint SIZE_BLOCKPARAMS = 6;    /**< Number of parameters per velocity block. */

cuint SIZE_BOUND       = WID3;
cuint SIZE_BDERI       = WID2;
cuint SIZE_BFLUX       = WID2;
cuint SIZE_DERIV       = WID3;
cuint SIZE_FLUXS       = WID3;

// Natural constants
namespace physicalconstants {
   //creal MU_0 = 1.25663706e-6;
   creal MU_0 = 1.0;
}

#endif
