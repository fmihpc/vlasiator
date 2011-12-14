/*!
  Spatial cell class for Vlasiator that supports a variable number of velocity blocks.
*/

#ifndef VLASIATOR_SPATIAL_CELL_HPP
#define VLASIATOR_SPATIAL_CELL_HPP

#include "algorithm"
#include "boost/array.hpp"

#ifndef NO_SPARSE
#include "boost/unordered_map.hpp"
#include "boost/unordered_set.hpp"
#endif

#include "cmath"
#include "fstream"
#include "iostream"
#include "mpi.h"
#include "limits"
#include "stdint.h"
#include "vector"

#include "common.h"
#include "parameters.h"
#include "definitions.h"
typedef Parameters P;


namespace spatial_cell {
   //fixme namespaces in lower case
   namespace Transfer {
      const unsigned int NONE=0;
      const unsigned int CELL_PARAMETERS=(1<<0);
      const unsigned int CELL_DERIVATIVES=(1<<1);
      const unsigned int VEL_BLOCK_LIST=(1<<2);
      const unsigned int VEL_BLOCK_DATA=(1<<3);
      const unsigned int VEL_BLOCK_FLUXES=(1<<4);
      const unsigned int VEL_BLOCK_KT_DERIVATIVES=(1<<5);
      const unsigned int VEL_BLOCK_PARAMETERS=(1<<6);
      const unsigned int CELL_B_RHO_RHOV=(1<<7);
      const unsigned int CELL_E=(1<<8);
      const unsigned int CELL_GHOSTFLAG=(1<<9);
      const unsigned int ALL = CELL_PARAMETERS|CELL_DERIVATIVES|VEL_BLOCK_LIST|VEL_BLOCK_DATA|
      VEL_BLOCK_FLUXES | VEL_BLOCK_KT_DERIVATIVES|CELL_GHOSTFLAG;
   };
   
// length of a velocity block in velocity cells
   const unsigned int block_len_x = WID;
   const unsigned int block_len_y = WID;
   const unsigned int block_len_z = WID;

   const Real cell_dvx = (P::vxmax - P::vxmin) / (P::vxblocks_ini * block_len_x);
   const Real cell_dvy = (P::vymax - P::vymin) / (P::vyblocks_ini * block_len_y);
   const Real cell_dvz = (P::vzmax - P::vzmin) / (P::vzblocks_ini * block_len_z);


/** A namespace for storing indices into an array containing neighbour information 
 * of velocity grid blocks. 
 */
namespace velocity_neighbor {
   const uint WIDTH = 3; //width of neighborhood
   const uint NNGBRS = WIDTH*WIDTH*WIDTH; //number of neighbors per block in velocity space (27)           
   const uint SIZE = NNGBRS+1; //number of neighbors per block in velocity space (27) plus one integer for flags
   const uint XM1_YM1_ZM1 = 0;  /**< Index of (x-1,y-1,z-1) neighbour.*/
   const uint XCC_YM1_ZM1 = 1;  /**< Index of (x  ,y-1,z-1) neighbour.*/
   const uint XP1_YM1_ZM1 = 2;  /**< Index of (x+1,y-1,z-1) neighbour.*/
   const uint XM1_YCC_ZM1 = 3;  /**< Index of (x-1,y  ,z-1) neighbour.*/
   const uint XCC_YCC_ZM1 = 4;  /**< Index of (x  ,y  ,z-1) neighbour.*/
   const uint XP1_YCC_ZM1 = 5;  /**< Index of (x+1,y  ,z-1) neighbour.*/
   const uint XM1_YP1_ZM1 = 6;  /**< Index of (x-1,y+1,z-1) neighbour.*/
   const uint XCC_YP1_ZM1 = 7;  /**< Index of (x  ,y+1,z-1) neighbour.*/
   const uint XP1_YP1_ZM1 = 8;  /**< Index of (x+1,y+1,z-1) neighbour.*/   
   const uint XM1_YM1_ZCC = 9;  /**< Index of (x-1,y-1,z  ) neighbour.*/
   const uint XCC_YM1_ZCC = 10; /**< Index of (x  ,y-1,z  ) neighbour.*/
   const uint XP1_YM1_ZCC = 11; /**< Index of (x+1,y-1,z  ) neighbour.*/
   const uint XM1_YCC_ZCC = 12; /**< Index of (x-1,y  ,z  ) neighbour.*/
   const uint XCC_YCC_ZCC = 13; /**< Index of (x  ,y  ,z  ) neighbour.*/
   const uint XP1_YCC_ZCC = 14; /**< Index of (x+1,y  ,z  ) neighbour.*/
   const uint XM1_YP1_ZCC = 15; /**< Index of (x-1,y+1,z  ) neighbour.*/
   const uint XCC_YP1_ZCC = 16; /**< Index of (x  ,y+1,z  ) neighbour.*/
   const uint XP1_YP1_ZCC = 17; /**< Index of (x+1,y+1,z  ) neighbour.*/   
   const uint XM1_YM1_ZP1 = 18; /**< Index of (x-1,y-1,z+1) neighbour.*/
   const uint XCC_YM1_ZP1 = 19; /**< Index of (x  ,y-1,z+1) neighbour.*/
   const uint XP1_YM1_ZP1 = 20; /**< Index of (x+1,y-1,z+1) neighbour.*/
   const uint XM1_YCC_ZP1 = 21; /**< Index of (x-1,y  ,z+1) neighbour.*/
   const uint XCC_YCC_ZP1 = 22; /**< Index of (x  ,y  ,z+1) neighbour.*/
   const uint XP1_YCC_ZP1 = 23; /**< Index of (x+1,y  ,z+1) neighbour.*/
   const uint XM1_YP1_ZP1 = 24; /**< Index of (x-1,y+1,z+1) neighbour.*/
   const uint XCC_YP1_ZP1 = 25; /**< Index of (x  ,y+1,z+1) neighbour.*/
   const uint XP1_YP1_ZP1 = 26; /**< Index of (x+1,y+1,z+1) neighbour.*/
//   const uint NON_EXISTING = std::numeric_limits<uint>::max(); /**< Invalid block ID, indicating that the block does not exist.*/
   const uint NBRFLAGS = 27; /**< Index for flags for existing neighbours.*/

   const uint MYIND    = 13; /**< Index of the block. Required for KT solver.*/
   const uint VXNEG    = 12; /**< Index of -vx neighbour. Required for KT solver.*/
   const uint VYNEG    = 10; /**< Index of -vy neighbour. Required for KT solver.*/
   const uint VZNEG    = 4;  /**< Index of -vz neighbour. Required for KT solver.*/
   const uint VXPOS    = 14; /**< Index of +vx neighbour. Required for KT solver.*/
   const uint VYPOS    = 16; /**< Index of +vy neighbour. Required for KT solver.*/
   const uint VZPOS    = 22; /**< Index of +vz neighbour. Required for KT solver.*/
   const uint VX_NEG_BND = (1 << VXNEG);
   const uint VX_POS_BND = (1 << VXPOS);
   const uint VY_NEG_BND = (1 << VYNEG);
   const uint VY_POS_BND = (1 << VYPOS);
   const uint VZ_NEG_BND = (1 << VZNEG);
   const uint VZ_POS_BND = (1 << VZPOS);
}

   
/*!
  
  Used as an error from functions returning velocity cells or
  as a cell that would be outside of the velocity block
*/
   const unsigned int error_velocity_cell = std::numeric_limits<unsigned int>::max();

/*!
  Defines the indices of a velocity cell in a velocity block.
  Indices start from 0 and the first value is the index in x direction.
*/
   typedef boost::array<unsigned int, 3> velocity_cell_indices_t;

/*!
  Used as an error from functions returning velocity cell indices or
  as an index that would be outside of the velocity block
*/
   const unsigned int error_velocity_cell_index = std::numeric_limits<unsigned int>::max();

   const unsigned int velocity_block_len = block_len_x * block_len_y * block_len_z;

// All velocity neigboring velocity blocks are neighbors   
   const unsigned int n_neighbor_velocity_blocks = 28;

   class Velocity_Block {
   public:
      // value of the distribution function
      Real data[velocity_block_len];
   
      // spatial fluxes of this block
      //fixme, fx could be called flux for leveque
      Real fx[SIZE_FLUXS];
#ifdef SOLVER_KT
      Real fy[SIZE_FLUXS];
      Real fz[SIZE_FLUXS];
// spatial derivatives of the distribution functio
      Real d1x[SIZE_DERIV];
      Real d2x[SIZE_DERIV];
      Real d1y[SIZE_DERIV];
      Real d2y[SIZE_DERIV];
      Real d1z[SIZE_DERIV];
      Real d2z[SIZE_DERIV];
#endif
      
      Real parameters[BlockParams::N_VELOCITY_BLOCK_PARAMS];
      Velocity_Block* neighbors[n_neighbor_velocity_blocks];
//      Velocity_Block* neighborsSpatial[SIZE_NBRS_SPA];
      
      /*!
        Sets data, derivatives and fluxes of this block to zero.
      */
      void clear(void)
         {
            for (unsigned int i = 0; i < velocity_block_len; i++) {
               this->data[i] = 0;
            }
#ifndef SOLVER_KT       
            for (unsigned int i = 0; i < SIZE_FLUXS; i++) {
               this->fx[i] = 0;
            }
#else   
            for (unsigned int i = 0; i < SIZE_FLUXS; i++) {
               this->fx[i] = 0;
               this->fy[i] = 0;
               this->fz[i] = 0;
            }

            for (unsigned int i = 0; i < SIZE_DERIV; i++) {
               this->d1x[i] = 0;
               this->d2x[i] = 0;
               this->d1y[i] = 0;
               this->d2y[i] = 0;
               this->d1z[i] = 0;
               this->d2z[i] = 0;
            }
#endif  
         }
   };

/*!
  Defines the indices of a velocity block in the velocity grid.
  Indices start from 0 and the first value is the index in x direction.
*/
   typedef boost::array<unsigned int, 3> velocity_block_indices_t;

/*!
  Used as an error from functions returning velocity blocks indices or
  as an index that would be outside of the velocity grid in this cell
*/
   const unsigned int error_velocity_block_index = std::numeric_limits<unsigned int>::max();

/*!
  Used as an error from functions returning velocity blocks or
  as a block that would be outside of the velocity grid in this cell
*/
   const unsigned int error_velocity_block = std::numeric_limits<unsigned int>::max();

   const unsigned int max_velocity_blocks = P::vxblocks_ini * P::vyblocks_ini * P::vzblocks_ini;


   const Real cell_dx = (P::vxmax - P::vxmin) / P::vxblocks_ini;
   const Real cell_dy = (P::vymax - P::vymin) / P::vyblocks_ini;
   const Real cell_dz = (P::vzmax - P::vzmin) / P::vzblocks_ini;

// TODO: typedef unsigned int velocity_cell_t;
// TODO: typedef unsigned int velocity_block_t;

/****************************
 * Velocity block functions *
 ****************************/

/*!
  Returns the indices of given velocity block
*/
   inline velocity_block_indices_t get_velocity_block_indices(const unsigned int block) {
      velocity_block_indices_t indices;

      if (block >= max_velocity_blocks) {
         indices[0] = indices[1] = indices[2] = error_velocity_block_index;
      } else {
         indices[0] = block % P::vxblocks_ini;
         indices[1] = (block / P::vxblocks_ini) % P::vyblocks_ini;
         indices[2] = block / (P::vxblocks_ini * P::vyblocks_ini);
      }

      return indices;
   }


/*!
  Returns the velocity block at given indices or error_velocity_block
*/
   inline   unsigned int get_velocity_block(const velocity_block_indices_t indices) {
      if (indices[0] >= P::vxblocks_ini
          || indices[1] >= P::vyblocks_ini
          || indices[2] >= P::vzblocks_ini) {
         return error_velocity_block;
      }

      return indices[0] + indices[1] * P::vxblocks_ini + indices[2] * P::vxblocks_ini * P::vyblocks_ini;
   }

/*!
  Returns the velocity block at given location or
  error_velocity_block if outside of the velocity grid
*/
   inline unsigned int get_velocity_block(
      const Real vx,
      const Real vy,
      const Real vz
      )
   {
      if (vx < P::vxmin || vx >= P::vxmax
          || vy < P::vymin || vy >= P::vymax
          || vz < P::vzmin || vz >= P::vzmax) {
         return error_velocity_block;
      }

      const velocity_block_indices_t indices = {
         (unsigned int) floor((vx - P::vxmin) / cell_dx),
         (unsigned int) floor((vy - P::vymin) / cell_dy),
         (unsigned int) floor((vz - P::vzmin) / cell_dz)
      };

      return get_velocity_block(indices);
   }

   /*!
  Returns the id of a velocity beock that is neighboring given block in given direction vx,vy,vz.
  Returns error_velocity_block in case the neighboring velocity block would be outside
  of the velocity grid.
*/
   inline unsigned int get_velocity_block(
      const unsigned int block,
      const  int direction_vx,
      const  int direction_vy,
      const  int direction_vz
      )
   {
      unsigned int neighborBlock;
      const velocity_block_indices_t indices = get_velocity_block_indices(block);
      if (indices[0] == error_velocity_block_index) {
         return error_velocity_block;
      }
      
      int xyzDirection[3];
      unsigned int vxyzblocks_ini[3];
      unsigned int offset[3];
      xyzDirection[0]=direction_vx;
      xyzDirection[1]=direction_vy;
      xyzDirection[2]=direction_vz;
      
      vxyzblocks_ini[0]=P::vxblocks_ini;
      vxyzblocks_ini[1]=P::vyblocks_ini;
      vxyzblocks_ini[2]=P::vzblocks_ini;

      offset[0]=1;
      offset[1]=vxyzblocks_ini[0];
      offset[2]=vxyzblocks_ini[0]*vxyzblocks_ini[1];
      
      
      neighborBlock=block;
      //loop over vx,vy,vz
      for(unsigned int c=0;c<3;c++){
         switch (xyzDirection[c]) {
             case -1: //in negative direction
                if (indices[c] == 0) {
                   return error_velocity_block;
                } else {
                    neighborBlock -= offset[c];
                }
                break;
             case 0:
                break;
             case 1: //in positive direction
                if (indices[c] >= vxyzblocks_ini[c] - 1) {
                   return error_velocity_block;
                } else {
                   neighborBlock += offset[c];
                }
                break;
             default:
                return error_velocity_block;
                break;
         }
      }
      return neighborBlock;
   }
   
   /*!
  Returns the id of a velocity block that is neighboring given block in given direction.
  Returns error_velocity_block in case the neighboring velocity block would be outside
  of the velocity grid.
*/

   
  inline unsigned int get_velocity_block(
      const unsigned int block,
      const unsigned int direction
      )
   {
      unsigned int w=velocity_neighbor::WIDTH;
      return get_velocity_block(direction%w-1,
                                (direction/w)%w-1,
                                (direction/(w*w))-1);
   }




/*!
  Returns the edge where given velocity block starts.
*/
   inline Real get_velocity_block_vx_min(const unsigned int block) {
      if (block == error_velocity_block) {
         return std::numeric_limits<Real>::quiet_NaN();
      }

      if (block >= max_velocity_blocks) {
         return std::numeric_limits<Real>::quiet_NaN();
      }

      const velocity_block_indices_t indices = get_velocity_block_indices(block);
      if (indices[0] == error_velocity_block_index) {
         return std::numeric_limits<Real>::quiet_NaN();
      }

      return P::vxmin + cell_dx * indices[0];
   }

/*!
  Returns the edge where given velocity block ends.
*/
   inline Real get_velocity_block_vx_max(const unsigned int block) {
      return get_velocity_block_vx_min(block) + cell_dx;
   }


/*!
  Returns the edge where given velocity block starts.
*/
   inline Real get_velocity_block_vy_min(const unsigned int block) {
      if (block == error_velocity_block) {
         return std::numeric_limits<Real>::quiet_NaN();
      }

      if (block >= max_velocity_blocks) {
         return std::numeric_limits<Real>::quiet_NaN();
      }

      const velocity_block_indices_t indices = get_velocity_block_indices(block);
      if (indices[1] == error_velocity_block_index) {
         return std::numeric_limits<Real>::quiet_NaN();
      }

      return P::vymin + cell_dy * indices[1];
   }

/*!
  Returns the edge where given velocity block ends.
*/
   inline Real get_velocity_block_vy_max(const unsigned int block) {
      return get_velocity_block_vy_min(block) + cell_dy;
   }


/*!
  Returns the edge where given velocity block starts.
*/
   inline Real get_velocity_block_vz_min(const unsigned int block) {
      if (block == error_velocity_block) {
         return std::numeric_limits<Real>::quiet_NaN();
      }

      if (block >= max_velocity_blocks) {
         return std::numeric_limits<Real>::quiet_NaN();
      }

      const velocity_block_indices_t indices = get_velocity_block_indices(block);
      if (indices[2] == error_velocity_block_index) {
         return std::numeric_limits<Real>::quiet_NaN();
      }

      return P::vzmin + cell_dz * indices[2];
   }

/*!
  Returns the edge where given velocity block ends.
*/
  inline Real get_velocity_block_vz_max(const unsigned int block) {
      return get_velocity_block_vz_min(block) + cell_dz;
   }


/***************************
 * Velocity cell functions *
 ***************************/

/*!
  Returns the indices of given velocity cell
*/
  inline velocity_cell_indices_t get_velocity_cell_indices(const unsigned int cell) {
      velocity_cell_indices_t indices;

      if (cell >= velocity_block_len) {
         indices[0] = indices[1] = indices[2] = error_velocity_cell_index;
      } else {
         indices[0] = cell % block_len_x;
         indices[1] = (cell / block_len_x) % block_len_y;
         indices[2] = cell / (block_len_x * block_len_y);
      }

      return indices;
   }

/*!
  Returns the velocity cell at given indices or error_velocity_cell
*/
  inline unsigned int get_velocity_cell(const velocity_cell_indices_t indices) {
      if (indices[0] >= block_len_x
          || indices[1] >= block_len_y
          || indices[2] >= block_len_z) {
         return error_velocity_cell;
      }
      return indices[0] + indices[1] * block_len_x + indices[2] * block_len_x * block_len_y;
   }


/*!
  Returns the id of a velocity cell that is neighboring given cell in given direction.
  Returns error_velocity_cell in case the neighboring velocity cell would be outside
  of the velocity block.
*/
  inline unsigned int get_velocity_cell(
      const unsigned int cell,
      const unsigned int direction
      )
   {
      int xyzDirection[3];
      unsigned int block_len[3];
      unsigned int offset[3];
      
      unsigned int neighborCell;
      
      const velocity_cell_indices_t indices = get_velocity_cell_indices(cell);
      if (indices[0] == error_velocity_cell_index) {
         return error_velocity_cell;
      }
      
      unsigned int w=velocity_neighbor::WIDTH;
      xyzDirection[0]=direction%w-1;
      xyzDirection[1]=(direction/w)%w-1;
      xyzDirection[2]=(direction/(w*w))%w-1; 
      
      block_len[0]=block_len_x;
      block_len[1]=block_len_y;
      block_len[2]=block_len_z;

      offset[0]=1;
      offset[1]=block_len_x;
      offset[2]=block_len_x * block_len_y;

      neighborCell=cell;
      //loop over vx,vy,vz
      for(int c=0;c<3;c++){
         switch (xyzDirection[c]) {
             case -1: //in negative direction
                if (indices[c] == 0) {
                   return error_velocity_cell;
                } else {
                   neighborCell+=-offset[c];
                }
                break;
             case 0:
                break;
             case 1: //in positive direction
                if (indices[c] >= block_len[c] - 1) {
                   return error_velocity_cell;
                } else {
                   neighborCell += offset[c];
                }
                break;
             default:
                return error_velocity_cell;
                break;
         }
      }
      return neighborCell;
   }

   
/*!     
  Returns the velocity cell at given location or
  error_velocity_cell if outside of given velocity block.
*/
   inline unsigned int get_velocity_cell(
      const unsigned int velocity_block,
      const Real vx,
      const Real vy,
      const Real vz
      )
   {
      const Real block_vx_min = get_velocity_block_vx_min(velocity_block);
      const Real block_vx_max = get_velocity_block_vx_max(velocity_block);
      const Real block_vy_min = get_velocity_block_vy_min(velocity_block);
      const Real block_vy_max = get_velocity_block_vy_max(velocity_block);
      const Real block_vz_min = get_velocity_block_vz_min(velocity_block);
      const Real block_vz_max = get_velocity_block_vz_max(velocity_block);

      if (vx < block_vx_min || vx >= block_vx_max
          || vy < block_vy_min || vy >= block_vy_max
          || vz < block_vz_min || vz >= block_vz_max
          ) {
         return error_velocity_cell;
      }

      const velocity_block_indices_t indices = {
         (unsigned int) floor((vx - block_vx_min) / ((block_vx_max - block_vx_min) / block_len_x)),
         (unsigned int) floor((vy - block_vy_min) / ((block_vy_max - block_vy_min) / block_len_y)),
         (unsigned int) floor((vz - block_vz_min) / ((block_vz_max - block_vz_min) / block_len_z))
      };

      return get_velocity_cell(indices);
   }


/*!
  Returns the edge where given velocity cell in the given velocity block starts.
  TODO: move these to velocity cell class?
*/
   inline Real get_velocity_cell_vx_min(
      const unsigned int velocity_block,
      const unsigned int velocity_cell
      )
   {
      if (velocity_cell == error_velocity_cell) {
         return std::numeric_limits<Real>::quiet_NaN();
      }

      const velocity_cell_indices_t indices = get_velocity_cell_indices(velocity_cell);
      if (indices[0] == error_velocity_cell_index) {
         return std::numeric_limits<Real>::quiet_NaN();
      }

      const Real block_vx_min = get_velocity_block_vx_min(velocity_block);
      const Real block_vx_max = get_velocity_block_vx_max(velocity_block);

      return block_vx_min + (block_vx_max - block_vx_min) / block_len_x * indices[0];
   }

/*!
  Returns the edge where given velocity cell in the given velocity block ends.
*/
 inline  Real get_velocity_cell_vx_max(
      const unsigned int velocity_block,
      const unsigned int velocity_cell
      )
   {
      if (velocity_cell == error_velocity_cell) {
         return std::numeric_limits<Real>::quiet_NaN();
      }

      const velocity_cell_indices_t indices = get_velocity_cell_indices(velocity_cell);
      if (indices[0] == error_velocity_cell_index) {
         return std::numeric_limits<Real>::quiet_NaN();
      }

      const Real block_vx_min = get_velocity_block_vx_min(velocity_block);
      const Real block_vx_max = get_velocity_block_vx_max(velocity_block);

      return block_vx_min + (block_vx_max - block_vx_min) / block_len_x * (indices[0] + 1);
   }

/*!
  Returns the edge where given velocity cell in the given velocity block starts.
*/
  inline Real get_velocity_cell_vy_min(
      const unsigned int velocity_block,
      const unsigned int velocity_cell
      )
   {
      if (velocity_cell == error_velocity_cell) {
         return std::numeric_limits<Real>::quiet_NaN();
      }

      const velocity_cell_indices_t indices = get_velocity_cell_indices(velocity_cell);
      if (indices[1] == error_velocity_cell_index) {
         return std::numeric_limits<Real>::quiet_NaN();
      }

      const Real block_vy_min = get_velocity_block_vy_min(velocity_block);
      const Real block_vy_max = get_velocity_block_vy_max(velocity_block);

      return block_vy_min + (block_vy_max - block_vy_min) / block_len_y * indices[1];
   }

/*!
  Returns the edge where given velocity cell in the given velocity block ends.
*/
   inline Real get_velocity_cell_vy_max(
      const unsigned int velocity_block,
      const unsigned int velocity_cell
      )
   {
      if (velocity_cell == error_velocity_cell) {
         return std::numeric_limits<Real>::quiet_NaN();
      }

      const velocity_cell_indices_t indices = get_velocity_cell_indices(velocity_cell);
      if (indices[1] == error_velocity_cell_index) {
         return std::numeric_limits<Real>::quiet_NaN();
      }

      const Real block_vy_min = get_velocity_block_vy_min(velocity_block);
      const Real block_vy_max = get_velocity_block_vy_max(velocity_block);

      return block_vy_min + (block_vy_max - block_vy_min) / block_len_y * (indices[1] + 1);
   }

/*!
  Returns the edge where given velocity cell in the given velocity block starts.
*/
  inline  Real get_velocity_cell_vz_min(
      const unsigned int velocity_block,
      const unsigned int velocity_cell
      )
   {
      if (velocity_cell == error_velocity_cell) {
         return std::numeric_limits<Real>::quiet_NaN();
      }

      const velocity_cell_indices_t indices = get_velocity_cell_indices(velocity_cell);
      if (indices[2] == error_velocity_cell_index) {
         return std::numeric_limits<Real>::quiet_NaN();
      }

      const Real block_vz_min = get_velocity_block_vz_min(velocity_block);
      const Real block_vz_max = get_velocity_block_vz_max(velocity_block);

      return block_vz_min + (block_vz_max - block_vz_min) / block_len_z * indices[2];
   }

/*!
  Returns the edge where given velocity cell in the given velocity block ends.
*/
   inline Real get_velocity_cell_vz_max(
      const unsigned int velocity_block,
      const unsigned int velocity_cell
      )
   {
      if (velocity_cell == error_velocity_cell) {
         return std::numeric_limits<Real>::quiet_NaN();
      }

      const velocity_cell_indices_t indices = get_velocity_cell_indices(velocity_cell);
      if (indices[2] == error_velocity_cell_index) {
         return std::numeric_limits<Real>::quiet_NaN();
      }

      const Real block_vz_min = get_velocity_block_vz_min(velocity_block);
      const Real block_vz_max = get_velocity_block_vz_max(velocity_block);

      return block_vz_min + (block_vz_max - block_vz_min) / block_len_z * (indices[2] + 1);
   }


   class SpatialCell {
   public:

      SpatialCell()
         {
            /*
              Block list always has room for all blocks
            */
            this->velocity_block_list.reserve(max_velocity_blocks);

#ifdef NO_SPARSE
            this->velocity_blocks.resize(max_velocity_blocks);
            for (unsigned int block = 0; block < max_velocity_blocks; block++) {
               this->velocity_block_list.push_back(block);
            }
#else
            this->block_address_cache.reserve(max_velocity_blocks);
            for (unsigned int block = 0; block < max_velocity_blocks; block++) {
               this->velocity_block_list.push_back(error_velocity_block);
               this->block_address_cache.push_back(&(this->null_block));
            }
#endif

            this->null_block.clear();
            // zero neighbor lists of null block
            for (unsigned int i = 0; i < n_neighbor_velocity_blocks; i++) {
               this->null_block.neighbors[i] = NULL;
            }

            this->velocity_block_min_value = 0;
            this->velocity_block_min_avg_value = 0;

            // reset spatial cell parameters
            for (unsigned int i = 0; i < CellParams::N_SPATIAL_CELL_PARAMS; i++) {
               this->parameters[i]=0;
            }

            // reset spatial cell derivatives
            for (unsigned int i = 0; i < fieldsolver::N_SPATIAL_CELL_DERIVATIVES; i++) {
               this->derivatives[i]=0;
            }
         }
   

       
      /*!
        Returns a reference to the given velocity block or to
        the null block if given velocity block doesn't exist.
      */
      Velocity_Block& at(const unsigned int block)
         {
            if (block == error_velocity_block
                || block >= max_velocity_blocks) {
               return this->null_block;
            } else {
#ifdef NO_SPARSE
               return this->velocity_blocks.at(block);
#else
               return *(this->block_address_cache.at(block));
#endif
            }
         }

      /*!
        A const version of the non-const at function.
      */
      Velocity_Block const& at(const unsigned int block) const
         {
            if (block == error_velocity_block
                || block >= max_velocity_blocks) {
               return this->null_block;
            } else {
#ifdef NO_SPARSE
               return this->velocity_blocks.at(block);
#else
               return *(this->block_address_cache.at(block));
#endif
            }
         }


      /*!
        Returns the number of given velocity blocks that exist.
      */
      size_t count(const unsigned int block) const
         {
#ifdef NO_SPARSE
            if (block == error_velocity_block
                || block >= max_velocity_blocks) {
               return 0;
            } else {
               return 1;
            }
#else
            return this->velocity_blocks.count(block);
#endif
         }


      /*!
        Returns the number of existing velocity blocks.
      */
      size_t size(void) const
         {
            return this->velocity_blocks.size();
         }


      /*!
        Sets the given value to a velocity cell at given coordinates.

        Creates the velocity block at given coordinates if it doesn't exist.
      */
      void set_value(const Real vx, const Real vy, const Real vz, const Real value)
         {
            const unsigned int block = get_velocity_block(vx, vy, vz);
#ifndef NO_SPARSE
            if (this->velocity_blocks.count(block) == 0) {
               if (!this->add_velocity_block(block)) {
                  std::cerr << "Couldn't add velocity block " << block << std::endl;
                  abort();
               }
            }
#endif
            std::cout << "Added block " << block << " in set_value" << std::endl;

            Velocity_Block* block_ptr = &(this->velocity_blocks.at(block));

            const unsigned int cell = get_velocity_cell(block, vx, vy, vz);
            block_ptr->data[cell] = value;
         }


      void* at(void)
         {
            return this;
         }


      MPI_Datatype mpi_datatype(void){
            MPI_Datatype type;
            std::vector<MPI_Aint> displacements;
            std::vector<int> block_lengths;
            unsigned int block_index = 0;

            //add data to send/recv to displacement and block length lists

            if((SpatialCell::mpi_transfer_type & Transfer::VEL_BLOCK_LIST)!=0){
               // send velocity block list  
               displacements.push_back((uint8_t*) &(this->velocity_block_list[0]) - (uint8_t*) this);
               block_lengths.push_back(sizeof(unsigned int) * max_velocity_blocks);
            }

            if((SpatialCell::mpi_transfer_type & Transfer::VEL_BLOCK_DATA)!=0){
               displacements.reserve(displacements.size()+this->velocity_blocks.size());
               block_lengths.reserve(block_lengths.size()+this->velocity_blocks.size());
               
               while (block_index < max_velocity_blocks
                      && this->velocity_block_list[block_index] != error_velocity_block) {

// TODO: use cached block addresses
                      displacements.push_back((uint8_t*) this->velocity_blocks.at(this->velocity_block_list[block_index]).data - (uint8_t*) this);
                      block_lengths.push_back(sizeof(Real) * velocity_block_len);
                      
                      block_index++;
               }
            }

            if((SpatialCell::mpi_transfer_type & Transfer::VEL_BLOCK_FLUXES)!=0){
               displacements.reserve(displacements.size()+this->velocity_blocks.size());
               block_lengths.reserve(block_lengths.size()+this->velocity_blocks.size());
               //send velocity block fluxes

               while (block_index < max_velocity_blocks
                      && this->velocity_block_list[block_index] != error_velocity_block) {
                  
                  displacements.push_back((uint8_t*) this->velocity_blocks.at(this->velocity_block_list[block_index]).fx - (uint8_t*) this);
#ifdef  SOLVER_KT                 
                  block_lengths.push_back(sizeof(Real) * 3 * SIZE_FLUXS);
#else  //leveque                                        
                  block_lengths.push_back(sizeof(Real) * SIZE_FLUXS);
#endif                        
                  block_index++;
               }
            }
            
            // send  spatial cell parameters
            if((SpatialCell::mpi_transfer_type & Transfer::CELL_PARAMETERS)!=0){
                displacements.push_back((uint8_t*) &(this->velocity_block_min_value) - (uint8_t*) this);
                block_lengths.push_back(sizeof(Real));
                
                displacements.push_back((uint8_t*) &(this->velocity_block_min_avg_value) - (uint8_t*) this);
                block_lengths.push_back(sizeof(Real));
                
                displacements.push_back((uint8_t*) &(this->parameters[0]) - (uint8_t*) this);
                block_lengths.push_back(sizeof(Real) * CellParams::N_SPATIAL_CELL_PARAMS);
            }

            // send  BX, BY, BZ, RHO, RHOVx, RHOVY, RHOVZ (order in enum should never change(!)
            if((SpatialCell::mpi_transfer_type & Transfer::CELL_B_RHO_RHOV)!=0){
               displacements.push_back((uint8_t*) &(this->parameters[CellParams::BX]) - (uint8_t*) this);
               block_lengths.push_back(sizeof(Real) * 7);
            }

            // send  BX, BY, BZ, RHO, RHOVx, RHOVY, RHOVZ (order in enum should never change(!)
            if((SpatialCell::mpi_transfer_type & Transfer::CELL_E)!=0){
               displacements.push_back((uint8_t*) &(this->parameters[CellParams::EX]) - (uint8_t*) this);
               block_lengths.push_back(sizeof(Real) * 3);
            }

            // send  spatial cell parameters
            if((SpatialCell::mpi_transfer_type & Transfer::CELL_DERIVATIVES)!=0){
                displacements.push_back((uint8_t*) &(this->derivatives[0]) - (uint8_t*) this);
                block_lengths.push_back(sizeof(Real) * fieldsolver::N_SPATIAL_CELL_DERIVATIVES);
            }

            // send  isGhostCell        
            if((SpatialCell::mpi_transfer_type & Transfer::CELL_GHOSTFLAG)!=0){
                displacements.push_back((uint8_t*) &(this->isGhostCell) - (uint8_t*) this);
                block_lengths.push_back(sizeof(Real));
            }
            
            
// send velocity block derivatives
            if((SpatialCell::mpi_transfer_type & Transfer::VEL_BLOCK_KT_DERIVATIVES)!=0){
#ifdef  SOLVER_KT
               displacements.reserve(displacements.size()+this->velocity_blocks.size());
               block_lengths.reserve(block_lengths.size()+this->velocity_blocks.size());
               while (block_index < max_velocity_blocks
                      && this->velocity_block_list[block_index] != error_velocity_block) {
                  
                  displacements.push_back((uint8_t*) this->velocity_blocks.at(this->velocity_block_list[block_index]).d1x - (uint8_t*) this);
                  block_lengths.push_back(sizeof(Real) * 6 * SIZE_DERIV);
                  block_index++;
               }
#endif                  
            }

            if((SpatialCell::mpi_transfer_type & Transfer::VEL_BLOCK_PARAMETERS)!=0){
               displacements.reserve(displacements.size()+this->velocity_blocks.size());
               block_lengths.reserve(block_lengths.size()+this->velocity_blocks.size());
               while (block_index < max_velocity_blocks
                      && this->velocity_block_list[block_index] != error_velocity_block) {
                  
                      // TODO: use cached block addresses
                  displacements.push_back((uint8_t*) this->velocity_blocks.at(this->velocity_block_list[block_index]).parameters - (uint8_t*) this);
                  block_lengths.push_back(sizeof(Real) * BlockParams::N_VELOCITY_BLOCK_PARAMS);
                  block_index++;
               }
            }

            
            
            if (displacements.size() > 0) {
               MPI_Type_create_hindexed(
                  displacements.size(),
                  &block_lengths[0],
                  &displacements[0],
                  MPI_BYTE,
                  &type
                  );
            } else {
               MPI_Type_contiguous(0, MPI_BYTE, &type);
            }
            
            return type;
         }
      

      /*!
        Sets the minimum velocity cell value of a distrubution function for
        that velocity block to be considered to have contents.
      */
      void set_block_minimum(const Real value)
         {
            this->velocity_block_min_value = value;
         }

      /*!
        Sets the minimum average velocity cell value of a distrubution function
        within a block for that block to be considered to have contents.
      */
      void set_block_average_minimum(const Real value)
         {
            this->velocity_block_min_avg_value = value;
         }

      /*!
        Returns true if given velocity block has enough of a distribution function.
        Returns false if the value of the distribution function is too low in every
        sense in given block.
        Also returns false if given block doesn't exist or is an error block.
      */
      bool velocity_block_has_contents(
#ifdef NO_SPARSE
         const unsigned int /*block*/
#else
         const unsigned int block
#endif
         ) const
         {
#ifndef NO_SPARSE
            if (block == error_velocity_block
                || this->velocity_blocks.count(block) == 0) {
               return false;
            }

            bool has_content = false;

            Real total = 0;
            const Velocity_Block* block_ptr = &(this->velocity_blocks.at(block));

            for (unsigned int i = 0; i < velocity_block_len; i++) {
               total += block_ptr->data[i];
               if (block_ptr->data[i] >= this->velocity_block_min_value) {
                  has_content = true;
                  break;
               }
            }

            if (total >= this->velocity_block_min_avg_value * velocity_block_len) {
               has_content = true;
            }

            return has_content;
#else
            return true;
#endif
         }


      /*!
        Returns the total value of the distribution function within this spatial cell.
      */
      Real get_total_value(void) const
         {
            Real total = 0;

            for (auto block = this->velocity_blocks.cbegin(); block != this->velocity_blocks.cend(); block++) {
               for (unsigned int i = 0; i < velocity_block_len; i++) {
#ifdef NO_SPARSE
                  total += block->data[i];
#else
                  total += block->second.data[i];
#endif
               }
            }

            return total;
         }


      /*!
        Returns the total size of the data in this spatial cell in bytes.

        Does not include velocity block lists, the null velocity block or velocity block neighbor lists.
      */
      size_t get_data_size(void) const
         {
            const unsigned int n = this->velocity_blocks.size();

            return 2 * sizeof(Real)
               + n * sizeof(unsigned int)
               + n * 2 * velocity_block_len * sizeof(Real);
         }


      /*!
        Checks velocity blocks in the velocity block list.
      */
      void check_velocity_block_list(void) const
         {
            for (unsigned int i = 0; i < max_velocity_blocks; i++) {
               if (this->velocity_block_list[i] == error_velocity_block) {
                  for (unsigned int j = i; j < max_velocity_blocks; j++) {
                     if (this->velocity_block_list[i] != error_velocity_block) {
                        std::cerr << __FILE__ << ":" << __LINE__
                                  << "Velocity block list has holes"
                                  << std::endl;
                        abort();
                     }
                  }
                  break;
               }

#ifndef NO_SPARSE
               if (this->velocity_blocks.count(this->velocity_block_list[i]) == 0) {
                  std::cerr << __FILE__ << ":" << __LINE__
                            << " Velocity block " << this->velocity_block_list[i]
                            << " doesn't exist"
                            << std::endl;
                  abort();
               }
#endif
            }
         }

      /*!
        Prints velocity blocks in the velocity block list.
      */
      void print_velocity_block_list(void) const
         {
            std::cout << this->velocity_blocks.size() << " blocks: ";
            for (unsigned int i = 0; i < max_velocity_blocks; i++) {

               if (this->velocity_block_list[i] == error_velocity_block) {
                  // debug
                  for (unsigned int j = i; j < max_velocity_blocks; j++) {
                     if (this->velocity_block_list[i] != error_velocity_block) {
                        std::cerr << "Velocity block list has holes" << std::endl;
                        abort();
                     }
                  }
                  break;
               }
               std::cout << this->velocity_block_list[i] << " ";
            }
            std::cout << std::endl;
         }

      /*!
        Prints given velocity block's velocity neighbor list.
      */
      void print_velocity_neighbor_list(const unsigned int block) const
         {
#ifndef NO_SPARSE
            if (this->velocity_blocks.count(block) == 0) {
               return;
            }
#endif

            const Velocity_Block* block_ptr = &(this->velocity_blocks.at(block));

            std::cout << block << " neighbors: ";
            for (unsigned int neighbor = 0; neighbor < n_neighbor_velocity_blocks; neighbor++) {
               if (block_ptr->neighbors[neighbor] == NULL) {
                  std::cout << "NULL ";
               } else if (block_ptr->neighbors[neighbor] == &(this->null_block)) {
                  std::cout << "null ";
               } else {
                  std::cout << block_ptr->neighbors[neighbor] << " ";
               }
            }
            std::cout << std::endl;
         }

      /*!
        Prints all velocity blocks' velocity neighbor list.
      */
      void print_velocity_neighbor_lists(void) const
         {
            for (unsigned int i = 0; i < max_velocity_blocks; i++) {

               if (this->velocity_block_list[i] == error_velocity_block) {
                  break;
               }

               this->print_velocity_neighbor_list(this->velocity_block_list[i]);
            }
         }

      /*!
        Adds "important" and removes "unimportant" velocity blocks to/from this cell.

        Removes all velocity blocks from this spatial cell which don't have content
        and don't have spatial or velocity neighbors with content.
        Adds neighbors for all velocity blocks which do have content (including spatial neighbors).
        Assumes that only blocks with a shared face in velocity space are neighbors.
        All cells in spatial_neighbors are assumed to be neighbors of this cell.
      */
      void adjust_velocity_blocks(
#ifdef NO_SPARSE
         const std::vector<SpatialCell*>& /*spatial_neighbors*/
#else
         const std::vector<SpatialCell*>& spatial_neighbors
#endif
         ) {
         // debug
         this->check_velocity_block_list();

#ifndef NO_SPARSE
         // don't iterate over blocks created / removed by this function
         std::vector<unsigned int> original_block_list;
         for (
            unsigned int block = this->velocity_block_list[0], block_i = 0;
            block_i < max_velocity_blocks
               && this->velocity_block_list[block_i] != error_velocity_block;
            block = this->velocity_block_list[++block_i]
            ) {
            original_block_list.push_back(block);
         }

         // get all velocity blocks with content in neighboring spatial cells
         boost::unordered_set<unsigned int> neighbors_with_content;
         for (std::vector<SpatialCell*>::const_iterator
                 neighbor = spatial_neighbors.begin();
              neighbor != spatial_neighbors.end();
              neighbor++
              ) {

            for (std::vector<unsigned int>::const_iterator
                    block = (*neighbor)->velocity_block_list.begin();
                 block != (*neighbor)->velocity_block_list.end();
                 block++
                 ) {
               if (*block == error_velocity_block) {
                  break;
               }

               if ((*neighbor)->velocity_block_has_contents(*block)) {
                  neighbors_with_content.insert(*block);
               }
            }
         }

         // remove all local blocks without content and without neighbors with content
         for (std::vector<unsigned int>::const_iterator
                 original = original_block_list.begin();
              original != original_block_list.end();
              original++
              ) {
            const bool original_has_content = this->velocity_block_has_contents(*original);

            if (original_has_content) {

               // add missing neighbors in velocity space
               for (unsigned int direction = 0; direction <velocity_neighbor::NNGBRS; direction++) {
                  const unsigned int neighbor_block = get_velocity_block(*original, direction);
                  if (neighbor_block == error_velocity_block) {
                     continue;
                  }

                  if (this->velocity_blocks.count(neighbor_block) > 0) {
                     continue;
                  }

                  if (!this->add_velocity_block(neighbor_block)) {
                     std::cerr << __FILE__ << ":" << __LINE__
                               << " Failed to add neighbor block " << neighbor_block
                               << " for block " << *original
                               << std::endl;
                     abort();
                  }
               }

            } else {

               // check if any neighbour has contents
               bool velocity_neighbors_have_content = false;
               for (unsigned int direction = 0; direction <velocity_neighbor::NNGBRS; direction++) {
                  const unsigned int neighbor_block = get_velocity_block(*original, direction);
                  if (this->velocity_block_has_contents(neighbor_block)) {
                     velocity_neighbors_have_content = true;
                     break;
                  }
               }
               
               if (!velocity_neighbors_have_content
                   && neighbors_with_content.count(*original) == 0) {
                  this->remove_velocity_block(*original);
               }
            }
         }

         // add local blocks for spatial neighbors with content
         for (boost::unordered_set<unsigned int>::const_iterator
                 neighbor = neighbors_with_content.begin();
              neighbor != neighbors_with_content.end();
              neighbor++
              ) {
            this->add_velocity_block(*neighbor);
         }
#endif
      }

      /*!
        Saves this spatial cell in vtk ascii format into the given filename.
      */
      void save_vtk(const char* filename) const {

         // do nothing if one cell or block dimension is 0
         if (block_len_x == 0
             || block_len_y == 0
             || block_len_z == 0
             || P::vxblocks_ini == 0
             || P::vyblocks_ini == 0
             || P::vzblocks_ini == 0) {
            return;
         }

         std::ofstream outfile(filename);
         if (!outfile.is_open()) {
            std::cerr << "Couldn't open file " << filename << std::endl;
            // TODO: throw an exception instead
            abort();
         }

         outfile << "# vtk DataFile Version 2.0" << std::endl;
         outfile << "Vlasiator spatial cell" << std::endl;
         outfile << "ASCII" << std::endl;
         outfile << "DATASET UNSTRUCTURED_GRID" << std::endl;

         // write separate points for every velocity cells' corners
         outfile << "POINTS " << (this->velocity_blocks.size() * velocity_block_len + 2) * 8 << " Real" << std::endl;
         for (std::vector<unsigned int>::const_iterator
                 block = this->velocity_block_list.begin();
              block != this->velocity_block_list.end();
              block++
              ) {

            if (*block == error_velocity_block) {
               // assume no blocks after first error block
               break;
            }

            for (unsigned int z_index = 0; z_index < block_len_z; z_index++)
               for (unsigned int y_index = 0; y_index < block_len_y; y_index++)
                  for (unsigned int x_index = 0; x_index < block_len_x; x_index++) {

                     const velocity_cell_indices_t indices = {x_index, y_index, z_index};
                     const unsigned int velocity_cell = get_velocity_cell(indices);

                     outfile << get_velocity_cell_vx_min(*block, velocity_cell) << " "
                             << get_velocity_cell_vy_min(*block, velocity_cell) << " "
                             << get_velocity_cell_vz_min(*block, velocity_cell) << std::endl;

                     outfile << get_velocity_cell_vx_max(*block, velocity_cell) << " "
                             << get_velocity_cell_vy_min(*block, velocity_cell) << " "
                             << get_velocity_cell_vz_min(*block, velocity_cell) << std::endl;

                     outfile << get_velocity_cell_vx_min(*block, velocity_cell) << " "
                             << get_velocity_cell_vy_max(*block, velocity_cell) << " "
                             << get_velocity_cell_vz_min(*block, velocity_cell) << std::endl;

                     outfile << get_velocity_cell_vx_max(*block, velocity_cell) << " "
                             << get_velocity_cell_vy_max(*block, velocity_cell) << " "
                             << get_velocity_cell_vz_min(*block, velocity_cell) << std::endl;

                     outfile << get_velocity_cell_vx_min(*block, velocity_cell) << " "
                             << get_velocity_cell_vy_min(*block, velocity_cell) << " "
                             << get_velocity_cell_vz_max(*block, velocity_cell) << std::endl;

                     outfile << get_velocity_cell_vx_max(*block, velocity_cell) << " "
                             << get_velocity_cell_vy_min(*block, velocity_cell) << " "
                             << get_velocity_cell_vz_max(*block, velocity_cell) << std::endl;

                     outfile << get_velocity_cell_vx_min(*block, velocity_cell) << " "
                             << get_velocity_cell_vy_max(*block, velocity_cell) << " "
                             << get_velocity_cell_vz_max(*block, velocity_cell) << std::endl;

                     outfile << get_velocity_cell_vx_max(*block, velocity_cell) << " "
                             << get_velocity_cell_vy_max(*block, velocity_cell) << " "
                             << get_velocity_cell_vz_max(*block, velocity_cell) << std::endl;
                  }
         }
         /*
           Add small velocity cells to the negative and positive corners of the grid
           so VisIt knows the maximum size of the velocity grid regardless of existing cells
         */
         outfile << P::vxmin - 0.1 * cell_dvx << " "
                 << P::vymin - 0.1 * cell_dvy << " "
                 << P::vzmin - 0.1 * cell_dvz << std::endl;
         outfile << P::vxmin << " "
                 << P::vymin - 0.1 * cell_dvy << " "
                 << P::vzmin - 0.1 * cell_dvz << std::endl;
         outfile << P::vxmin - 0.1 * cell_dvx << " "
                 << P::vymin << " "
                 << P::vzmin - 0.1 * cell_dvz << std::endl;
         outfile << P::vxmin << " "
                 << P::vymin << " "
                 << P::vzmin - 0.1 * cell_dvz << std::endl;
         outfile << P::vxmin - 0.1 * cell_dvx << " "
                 << P::vymin - 0.1 * cell_dvy << " "
                 << P::vzmin << std::endl;
         outfile << P::vxmin << " "
                 << P::vymin - 0.1 * cell_dvy << " "
                 << P::vzmin << std::endl;
         outfile << P::vxmin - 0.1 * cell_dvx << " "
                 << P::vymin << " "
                 << P::vzmin << std::endl;
         outfile << P::vxmin << " "
                 << P::vymin << " "
                 << P::vzmin << std::endl;

         outfile << P::vxmax << " "
                 << P::vymax << " "
                 << P::vzmax << std::endl;
         outfile << P::vxmax + 0.1 * cell_dvx << " "
                 << P::vymax << " "
                 << P::vzmax << std::endl;
         outfile << P::vxmax << " "
                 << P::vymax + 0.1 * cell_dvy << " "
                 << P::vzmax << std::endl;
         outfile << P::vxmax + 0.1 * cell_dvx << " "
                 << P::vymax + 0.1 * cell_dvy << " "
                 << P::vzmax << std::endl;
         outfile << P::vxmax << " "
                 << P::vymax << " "
                 << P::vzmax + 0.1 * cell_dvz << std::endl;
         outfile << P::vxmax + 0.1 * cell_dvx << " "
                 << P::vymax << " "
                 << P::vzmax + 0.1 *  cell_dvz << std::endl;
         outfile << P::vxmax << " "
                 << P::vymax + 0.1 * cell_dvy << " "
                 << P::vzmax + 0.1 * cell_dvz << std::endl;
         outfile << P::vxmax + 0.1 * cell_dvx << " "
                 << P::vymax + 0.1 * cell_dvy << " "
                 << P::vzmax + 0.1 * cell_dvz << std::endl;

         // map cells to written points
         outfile << "CELLS "
                 << this->velocity_blocks.size() * velocity_block_len + 2 << " "
                 << (this->velocity_blocks.size() * velocity_block_len + 2)* 9 << std::endl;

         unsigned int j = 0;
         for (std::vector<unsigned int>::const_iterator
                 block = this->velocity_block_list.begin();
              block != this->velocity_block_list.end();
              block++
              ) {

            if (*block == error_velocity_block) {
               // assume no blocks after first error block
               break;
            }

            for (unsigned int z_index = 0; z_index < block_len_z; z_index++)
               for (unsigned int y_index = 0; y_index < block_len_y; y_index++)
                  for (unsigned int x_index = 0; x_index < block_len_x; x_index++) {

                     outfile << "8 ";
                     for (int i = 0; i < 8; i++) {
                        outfile << j * 8 + i << " ";
                     }
                     outfile << std::endl;

                     j++;
                  }
         }
         outfile << "8 ";
         for (unsigned int i = 0; i < 8; i++) {
            outfile << j * 8 + i << " ";
         }
         outfile << std::endl;
         outfile << "8 ";
         for (unsigned int i = 0; i < 8; i++) {
            outfile << j * 8 + i << " ";
         }
         outfile << std::endl;

         // cell types
         outfile << "CELL_TYPES " << this->velocity_blocks.size() * velocity_block_len + 2 << std::endl;
         for (unsigned int i = 0; i < this->velocity_blocks.size() * velocity_block_len + 2; i++) {
            outfile << 11 << std::endl;
         }

         // Put minimum value from existing blocks into two additional cells
         Real min_value = std::numeric_limits<Real>::max();

         // distribution function
         outfile << "CELL_DATA " << this->velocity_blocks.size() * velocity_block_len + 2 << std::endl;
         outfile << "SCALARS rho Real 1" << std::endl;
         outfile << "LOOKUP_TABLE default" << std::endl;
         for (std::vector<unsigned int>::const_iterator
                 block = this->velocity_block_list.begin();
              block != this->velocity_block_list.end();
              block++
              ) {

            if (*block == error_velocity_block) {
               // assume no blocks after first error block
               if (min_value == std::numeric_limits<Real>::max()) {
                  min_value = 0;
               }
               break;
            }

            const Velocity_Block* block_ptr = &(this->velocity_blocks.at(*block));

            for (unsigned int z_index = 0; z_index < block_len_z; z_index++)
               for (unsigned int y_index = 0; y_index < block_len_y; y_index++)
                  for (unsigned int x_index = 0; x_index < block_len_x; x_index++) {

                     const velocity_cell_indices_t indices = {x_index, y_index, z_index};
                     const unsigned int velocity_cell = get_velocity_cell(indices);
                     const Real value = block_ptr->data[velocity_cell];
                     outfile << value << " ";
                     min_value = std::min(min_value, value);
                  }
            outfile << std::endl;
         }
         outfile << min_value << " " << min_value << std::endl;

         if (!outfile.good()) {
            std::cerr << "Writing of vtk file probably failed" << std::endl;
            // TODO: throw an exception instead
            abort();
         }

         outfile.close();
      }


      /*!
        Adds an empty velocity block into this spatial cell.

        Returns true if given block was added or already exists.
        Returns false if given block is invalid or would be outside
        of the velocity grid.
      */
      bool add_velocity_block(const unsigned int block) {
         if (block == error_velocity_block) {
            return false;
         }

         if (block >= max_velocity_blocks) {
            return false;
         }

#ifndef NO_SPARSE
         if (this->velocity_blocks.count(block) > 0) {
            return true;
         }
#endif

#ifdef NO_SPARSE
         // assume blocks were added in default constructor
#else
         this->velocity_blocks[block];
         this->block_address_cache[block] = &(this->velocity_blocks.at(block));
#endif

#ifdef NO_SPARSE
         Velocity_Block* block_ptr = &(this->velocity_blocks.at(block));
#else
         Velocity_Block* block_ptr = this->block_address_cache[block];
#endif

         block_ptr->clear();

         // set block parameters
         block_ptr->parameters[BlockParams::VXCRD] = get_velocity_block_vx_min(block);
         block_ptr->parameters[BlockParams::VYCRD] = get_velocity_block_vy_min(block);
         block_ptr->parameters[BlockParams::VZCRD] = get_velocity_block_vz_min(block);
         block_ptr->parameters[BlockParams::DVX] =
            (get_velocity_block_vx_max(block) - get_velocity_block_vx_min(block)) / 2;
         block_ptr->parameters[BlockParams::DVY] =
            (get_velocity_block_vy_max(block) - get_velocity_block_vy_min(block)) / 2;
         block_ptr->parameters[BlockParams::DVZ] =
            (get_velocity_block_vz_max(block) - get_velocity_block_vz_min(block)) / 2;

         // set neighbour pointers
         unsigned int neighbor_block;
         int neighbor_index=0;
         for(int neighbor_index=0;neighbor_index<velocity_neighbor::NNGBRS;neighbor_index++){
            if(neighbor_index==velocity_neighbor::MYIND) continue; //self index
            neighbor_block = get_velocity_block(block, neighbor_index);
            if (neighbor_block == error_velocity_block) {
               block_ptr->neighbors[neighbor_index] = NULL;
            }
#ifndef NO_SPARSE   
            else if (this->velocity_blocks.count(neighbor_block) == 0) {
               block_ptr->neighbors[neighbor_index] = &(this->null_block);
            }
#endif      
            else {
               block_ptr->neighbors[neighbor_index] = &(this->velocity_blocks.at(neighbor_block));
               
               //  update the neighbor list of neighboring block
               unsigned int w=velocity_neighbor::WIDTH;
               int dvx=neighbor_index%w-1;
               int dvy=(neighbor_index/w)%w-1;
               int dvz=(neighbor_index/(w*w))-1;
               //index of current block in neighbors neighbor table
               int neighbor_neighbor_index=(1-dvx)+(1-dvy)*3+(1-dvz)*9;
               Velocity_Block* neighbor_ptr = &(this->velocity_blocks.at(neighbor_block));
               neighbor_ptr->neighbors[neighbor_neighbor_index] = block_ptr;
            }
         }


#ifndef NO_SPARSE
         unsigned int first_error_block = 0;
         while (first_error_block < max_velocity_blocks
                && this->velocity_block_list[first_error_block] != error_velocity_block
                // block is in the list when preparing to receive blocks
                && this->velocity_block_list[first_error_block] != block
                ) {
            first_error_block++;
         }
         
         this->velocity_block_list[first_error_block] = block;
#endif                  
         
         return true;
      }
      
      /*!
        Adds all velocity blocks than don't exist into the velocity grid.

        Returns true if all non-existing blocks were added, false otherwise.
      */
      bool add_all_velocity_blocks(void)
         {
            bool result = true;

            for (unsigned int i = 0; i < max_velocity_blocks; i++) {
#ifndef NO_SPARSE
               if (this->velocity_blocks.count(i) > 0) {
                  continue;
               }
#endif

               if (!this->add_velocity_block(i)) {
                  result = false;
               }
            }

            return result;
         }


      /*!
        Removes given block from the velocity grid.
        Does nothing if given block doesn't exist.
      */
      void remove_velocity_block(
#ifdef NO_SPARSE
         const unsigned int /*block*/
#else
         const unsigned int block
#endif
         ) {
#ifndef NO_SPARSE
         if (block == error_velocity_block) {
            return;
         }

         if (this->velocity_blocks.count(block) == 0) {
            return;
         }

         // set neighbour pointers to this block to NULL
         unsigned int neighbour_block;


         //remove block from neighbors neighbor lists
         unsigned int neighbor_block;
         int neighbor_index=0;
         for(int neighbor_index=0;neighbor_index<velocity_neighbor::NNGBRS;neighbor_index++){
            if(neighbor_index==velocity_neighbor::MYIND) continue;
            neighbor_block = get_velocity_block(block, neighbor_index);
            if (neighbour_block != error_velocity_block
                && this->velocity_blocks.count(neighbour_block) > 0) {
               // TODO use cached addresses of neighbors
               Velocity_Block* neighbour_data = &(this->velocity_blocks.at(neighbour_block));
               //  update the neighbor list of neighboring block
               unsigned int w=velocity_neighbor::WIDTH;
               int dvx=neighbor_index%w-1;
               int dvy=(neighbor_index/w)%w-1;
               int dvz=(neighbor_index/(w*w))-1;
               //index of current block in neighbors neighbor table
               int neighbor_neighbor_index=(1-dvx)+(1-dvy)*3+(1-dvz)*9;
               neighbour_data->neighbors[neighbor_neighbor_index] = &(this->null_block);
            }
         }
         

         this->velocity_blocks.erase(block);
         this->block_address_cache[block] = &(this->null_block);

         /*
           Move the last existing block in the block list
           to the removed block's position
         */
         unsigned int block_index = 0;
         while (block_index < max_velocity_blocks
                && this->velocity_block_list[block_index] != block) {
            block_index++;
         }
         //debug
         if (block_index == max_velocity_blocks) {
            std::cerr << "Velocity block " << block << " not found in list" << std::endl;
            abort();
         }

         unsigned int first_error_block = 0;
         while (first_error_block < max_velocity_blocks
                && this->velocity_block_list[first_error_block] != error_velocity_block) {
            first_error_block++;
         }

         if (block_index == first_error_block - 1) {
            this->velocity_block_list[block_index] = error_velocity_block;
         } else {
            this->velocity_block_list[block_index] = this->velocity_block_list[first_error_block - 1];
            this->velocity_block_list[first_error_block - 1] = error_velocity_block;
         }
#endif
      }


      /*!
        Removes all velocity blocks from this spatial cell.
      */
      void clear(void)
         {
#ifndef NO_SPARSE
            this->velocity_blocks.clear();
            for (unsigned int i = 0; i < max_velocity_blocks; i++) {
               this->block_address_cache[i] = &(this->null_block);
               this->velocity_block_list[i] = error_velocity_block;
            }
#endif
         }


      /*!
        Prepares this spatial cell to receive the velocity grid over MPI.
      */
      void prepare_to_receive_blocks(void)
         {
#ifndef NO_SPARSE
            unsigned int number_of_blocks = 0;
            while (number_of_blocks < max_velocity_blocks
                   && this->velocity_block_list[number_of_blocks] != error_velocity_block) {
               number_of_blocks++;
            }

            // add_velocity_block overwrites the block list every time so:
            std::vector<unsigned int> old_block_list(number_of_blocks, error_velocity_block);
            for (unsigned int block_index = 0; block_index < number_of_blocks; block_index++) {
               old_block_list[block_index] = this->velocity_block_list[block_index];
            }

            this->clear();

            // add velocity blocks that are about to be received with MPI
            for (unsigned int block_index = 0; block_index < number_of_blocks; block_index++) {
               this->add_velocity_block(old_block_list[block_index]);
            }
#endif
         }


      /*!
        Sets the type of data to transfer by mpi_datatype.
      */
      static void set_mpi_transfer_type(const unsigned int type)
         {
            SpatialCell::mpi_transfer_type = type;
         }

      /*!
        Gets the type of data that will be transferred by mpi_datatype.
      */
      static unsigned int get_mpi_transfer_type(void)
         {
            return SpatialCell::mpi_transfer_type;
         }



   private:
      /*
        Which data is transferred by the mpi datatype given by spatial cells.
      */
      static unsigned int mpi_transfer_type;


      /*
        Minimum value of distribution function
        in any cell of a velocity block for the
        block to be considered to have contents
      */
      Real velocity_block_min_value;

      /*
        Minimum value of the average of distribution
        function within a velocity block for the
        block to be considered to have contents
      */
      Real velocity_block_min_avg_value;

      /*!
        Used as a neighbour instead of blocks that don't
        exist but would be inside of the velocity grid.
        Neighbors that would be outside of the grid are always NULL.
      */
      Velocity_Block null_block;

      // data of velocity blocks that exist in this cell
#ifdef NO_SPARSE
      
      std::vector<Velocity_Block> velocity_blocks;
#else
      boost::unordered_map<unsigned int, Velocity_Block> velocity_blocks;
      /*
        Speed up search of velocity block addresses in the hash table above.
        Addresses of non-existing blocks point to the null block.
      */
      //FIXME, static array?      
      std::vector<Velocity_Block*> block_address_cache;
#endif



   public:
      /*
        List of velocity blocks in this cell, used for
        transferring spatial cells between processes using MPI.
      */
      
      std::vector<unsigned int> velocity_block_list;
      
      /*
        Bulk variables in this spatial cell.
      */
//      std::vector<Real> parameters;
      Real parameters[CellParams::N_SPATIAL_CELL_PARAMS];
      /*
        Derivatives of bulk variables in this spatial cell.
      */
      Real derivatives[fieldsolver::N_SPATIAL_CELL_DERIVATIVES];
      
      //neighbor id's. Kept up to date in solvers, not by the spatial_cell class
      std::vector<uint64_t> neighbors;
      unsigned int boundaryFlag;
      bool isGhostCell;
      
//      std::vector<Real> derivatives;

   }; // class SpatialCell
   

} // namespaces
#endif

