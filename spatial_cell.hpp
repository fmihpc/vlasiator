/*!
Spatial cell class for Vlasiator that supports a variable number of velocity blocks.

Copyright 2011 Finnish Meteorological Institute

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

#ifndef VLASIATOR_SPATIAL_CELL_HPP
#define VLASIATOR_SPATIAL_CELL_HPP

#include "algorithm"
#include "boost/array.hpp"
#include "boost/unordered_map.hpp"
#include "boost/unordered_set.hpp"
#include "boost/lexical_cast.hpp"
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

// size of velocity blocks in velocity cells
#define block_vx_length WID
#define block_vy_length WID
#define block_vz_length WID
#define VELOCITY_BLOCK_LENGTH (block_vx_length * block_vy_length * block_vz_length)

#define N_NEIGHBOR_VELOCITY_BLOCKS 28


/*!
Used as an error from functions returning velocity cells or
as a cell that would be outside of the velocity block
*/
#define error_velocity_cell 0xFFFFFFFFu

/*!
Used as an error from functions returning velocity cell indices or
as an index that would be outside of the velocity block
*/
#define error_velocity_cell_index 0xFFFFFFFFu

/*!
Used as an error from functions returning velocity blocks or
as a block that would be outside of the velocity grid in this cell
*/
#define error_velocity_block 0xFFFFFFFFu

/*!
Used as an error from functions returning velocity blocks indices or
as an index that would be outside of the velocity grid in this cell
*/
#define error_velocity_block_index 0xFFFFFFFFu



namespace spatial_cell {
   //fixme namespaces in lower case
   namespace Transfer {
      const unsigned int NONE                     = 0;
      const unsigned int CELL_PARAMETERS          = (1<<0);
      const unsigned int CELL_DERIVATIVES         = (1<<1);
      const unsigned int VEL_BLOCK_LIST           = (1<<2);
      const unsigned int VEL_BLOCK_DATA           = (1<<3);
      const unsigned int VEL_BLOCK_FLUXES         = (1<<4);
      const unsigned int VEL_BLOCK_KT_DERIVATIVES = (1<<5);
      const unsigned int VEL_BLOCK_PARAMETERS     = (1<<6);
      const unsigned int CELL_B_RHO_RHOV          = (1<<7);
      const unsigned int CELL_E                   = (1<<8);
      const unsigned int CELL_GHOSTFLAG           = (1<<9);
//does not update block-lists, need to be updated in separate stage
      const unsigned int ALL_DATA =
      CELL_PARAMETERS
      | CELL_DERIVATIVES
      | VEL_BLOCK_DATA
      | VEL_BLOCK_FLUXES
      | CELL_GHOSTFLAG
      | VEL_BLOCK_KT_DERIVATIVES;
   };
   
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
  Defines the indices of a velocity cell in a velocity block.
  Indices start from 0 and the first value is the index in x direction.
*/
   typedef boost::array<unsigned int, 3> velocity_cell_indices_t;


   /* FIXME FIXME   No copy constructor, even if we need one. At the moment it is called at initialization of spatial cells, which also calls its copy constructors. As it is empty of blocks, no problems arise because of our missing velocity block copy constructor. The problem is the neigbors member, which has direct pointers to neighboring blocks in velocity space. In a copy these would be invalid....
   */
   
   class Velocity_Block {
   public:

      // value of the distribution function
      Real data[VELOCITY_BLOCK_LENGTH];
   
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
      Velocity_Block* neighbors[N_NEIGHBOR_VELOCITY_BLOCKS];
      
      /*!
        Sets data, derivatives and fluxes of this block to zero.
      */
      void clear(void)
         {
            for (unsigned int i = 0; i < VELOCITY_BLOCK_LENGTH; i++) {
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

// TODO: typedef unsigned int velocity_cell_t;
// TODO: typedef unsigned int velocity_block_t;


   class SpatialCell {
   public:

      /****************************
       * Velocity block functions *
       ****************************/

      /*!
      Returns the indices of given velocity block
      */
      static velocity_block_indices_t get_velocity_block_indices(const unsigned int block) {
         velocity_block_indices_t indices;

         if (block >= SpatialCell::max_velocity_blocks) {
            indices[0] = indices[1] = indices[2] = error_velocity_block_index;
         } else {
            indices[0] = block % SpatialCell::vx_length;
            indices[1] = (block / SpatialCell::vx_length) % SpatialCell::vy_length;
            indices[2] = block / (SpatialCell::vx_length * SpatialCell::vy_length);
         }

         return indices;
      }


      /*!
      Returns the velocity block at given indices or error_velocity_block
      */
      static unsigned int get_velocity_block(const velocity_block_indices_t indices) {
         if (indices[0] >= SpatialCell::vx_length
             || indices[1] >= SpatialCell::vy_length
             || indices[2] >= SpatialCell::vz_length) {
            return error_velocity_block;
         }

         return indices[0]
            + indices[1] * SpatialCell::vx_length
            + indices[2] * SpatialCell::vx_length * SpatialCell::vy_length;
      }

      /*!
      Returns the velocity block at given location or
      error_velocity_block if outside of the velocity grid
      */
      static unsigned int get_velocity_block(
         const Real vx,
         const Real vy,
         const Real vz
      ) {
         if (vx < SpatialCell::vx_min || vx >= SpatialCell::vx_max
             || vy < SpatialCell::vy_min || vy >= SpatialCell::vy_max
             || vz < SpatialCell::vz_min || vz >= SpatialCell::vz_max) {
            return error_velocity_block;
         }

         const velocity_block_indices_t indices = {
            (unsigned int) floor((vx - SpatialCell::vx_min) / SpatialCell::block_dvx),
            (unsigned int) floor((vy - SpatialCell::vy_min) / SpatialCell::block_dvy),
            (unsigned int) floor((vz - SpatialCell::vz_min) / SpatialCell::block_dvz)
         };

         return SpatialCell::get_velocity_block(indices);
      }
      
      /*!
      Returns the edge where given velocity block starts.
      */
      static Real get_velocity_block_vx_min(const unsigned int block) {
         if (block == error_velocity_block) {
            return std::numeric_limits<Real>::quiet_NaN();
         }

         if (block >= SpatialCell::max_velocity_blocks) {
            return std::numeric_limits<Real>::quiet_NaN();
         }

         const velocity_block_indices_t indices = get_velocity_block_indices(block);
         if (indices[0] == error_velocity_block_index) {
            return std::numeric_limits<Real>::quiet_NaN();
         }

         return SpatialCell::vx_min + SpatialCell::block_dvx * indices[0];
      }

      /*!
      Returns the edge where given velocity block ends.
      */
      static Real get_velocity_block_vx_max(const unsigned int block) {
         return SpatialCell::get_velocity_block_vx_min(block) + SpatialCell::block_dvx;
      }

      /*!
      Returns the edge where given velocity block starts.
      */
      static Real get_velocity_block_vy_min(const unsigned int block) {
         if (block == error_velocity_block) {
            return std::numeric_limits<Real>::quiet_NaN();
         }

         if (block >= SpatialCell::max_velocity_blocks) {
            return std::numeric_limits<Real>::quiet_NaN();
         }

         const velocity_block_indices_t indices = get_velocity_block_indices(block);
         if (indices[1] == error_velocity_block_index) {
            return std::numeric_limits<Real>::quiet_NaN();
         }

         return SpatialCell::vy_min + SpatialCell::block_dvy * indices[1];
      }

      /*!
      Returns the edge where given velocity block ends.
      */
      static Real get_velocity_block_vy_max(const unsigned int block) {
         return SpatialCell::get_velocity_block_vy_min(block) + SpatialCell::block_dvy;
      }

      /*!
      Returns the edge where given velocity block starts.
      */
      static Real get_velocity_block_vz_min(const unsigned int block) {
         if (block == error_velocity_block) {
            return std::numeric_limits<Real>::quiet_NaN();
         }

         if (block >= SpatialCell::max_velocity_blocks) {
            return std::numeric_limits<Real>::quiet_NaN();
         }

         const velocity_block_indices_t indices = get_velocity_block_indices(block);
         if (indices[2] == error_velocity_block_index) {
            return std::numeric_limits<Real>::quiet_NaN();
         }

         return SpatialCell::vz_min + SpatialCell::block_dvz * indices[2];
      }

      /*!
      Returns the edge where given velocity block ends.
      */
      static Real get_velocity_block_vz_max(const unsigned int block) {
          return SpatialCell::get_velocity_block_vz_min(block) + SpatialCell::block_dvz;
      }


      /***************************
       * Velocity cell functions *
       ***************************/

      /*!
      Returns the indices of given velocity cell
      */
      static velocity_cell_indices_t get_velocity_cell_indices(const unsigned int cell) {
         velocity_cell_indices_t indices;

         if (cell >= VELOCITY_BLOCK_LENGTH) {
            indices[0] = indices[1] = indices[2] = error_velocity_cell_index;
         } else {
            indices[0] = cell % block_vx_length;
            indices[1] = (cell / block_vx_length) % block_vy_length;
            indices[2] = cell / (block_vx_length * block_vy_length);
         }

         return indices;
      }

      /*!
      Returns the velocity cell at given indices or error_velocity_cell
      */
      static unsigned int get_velocity_cell(const velocity_cell_indices_t indices) {
         if (indices[0] >= block_vx_length
             || indices[1] >= block_vy_length
             || indices[2] >= block_vz_length) {
            return error_velocity_cell;
         }
         return indices[0] + indices[1] * block_vx_length + indices[2] * block_vx_length * block_vy_length;
      }


      /*!
      Returns the id of a velocity cell that is neighboring given cell in given direction.
      Returns error_velocity_cell in case the neighboring velocity cell would be outside
      of the velocity block.
      */
      static unsigned int get_velocity_cell(
         const unsigned int cell,
         const unsigned int direction
      ) {

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
      
      block_len[0]=block_vx_length;
      block_len[1]=block_vy_length;
      block_len[2]=block_vz_length;

      offset[0]=1;
      offset[1]=block_vx_length;
      offset[2]=block_vx_length * block_vy_length;

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
      static unsigned int get_velocity_cell(
         const unsigned int velocity_block,
         const Real vx,
         const Real vy,
         const Real vz
         ) {
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
            (unsigned int) floor((vx - block_vx_min) / ((block_vx_max - block_vx_min) / block_vx_length)),
            (unsigned int) floor((vy - block_vy_min) / ((block_vy_max - block_vy_min) / block_vy_length)),
            (unsigned int) floor((vz - block_vz_min) / ((block_vz_max - block_vz_min) / block_vz_length))
         };

         return SpatialCell::get_velocity_cell(indices);
      }


      /*!
      Returns the edge where given velocity cell in the given velocity block starts.
      TODO: move these to velocity cell class?
      */
      static Real get_velocity_cell_vx_min(
         const unsigned int velocity_block,
         const unsigned int velocity_cell
      ) {
         if (velocity_cell == error_velocity_cell) {
             return std::numeric_limits<Real>::quiet_NaN();
         }

         const velocity_cell_indices_t indices = get_velocity_cell_indices(velocity_cell);
         if (indices[0] == error_velocity_cell_index) {
            return std::numeric_limits<Real>::quiet_NaN();
         }

         const Real block_vx_min = get_velocity_block_vx_min(velocity_block);
         const Real block_vx_max = get_velocity_block_vx_max(velocity_block);

         return block_vx_min + (block_vx_max - block_vx_min) / block_vx_length * indices[0];
      }

      /*!
      Returns the edge where given velocity cell in the given velocity block ends.
      */
      static Real get_velocity_cell_vx_max(
         const unsigned int velocity_block,
         const unsigned int velocity_cell
      ) {
         if (velocity_cell == error_velocity_cell) {
            return std::numeric_limits<Real>::quiet_NaN();
         }

         const velocity_cell_indices_t indices = get_velocity_cell_indices(velocity_cell);
         if (indices[0] == error_velocity_cell_index) {
            return std::numeric_limits<Real>::quiet_NaN();
         }

         const Real block_vx_min = get_velocity_block_vx_min(velocity_block);
         const Real block_vx_max = get_velocity_block_vx_max(velocity_block);

         return block_vx_min + (block_vx_max - block_vx_min) / block_vx_length * (indices[0] + 1);
      }

      /*!
      Returns the edge where given velocity cell in the given velocity block starts.
      */
      static Real get_velocity_cell_vy_min(
         const unsigned int velocity_block,
         const unsigned int velocity_cell
      ) {
         if (velocity_cell == error_velocity_cell) {
            return std::numeric_limits<Real>::quiet_NaN();
         }

         const velocity_cell_indices_t indices = get_velocity_cell_indices(velocity_cell);
         if (indices[1] == error_velocity_cell_index) {
            return std::numeric_limits<Real>::quiet_NaN();
         }

         const Real block_vy_min = get_velocity_block_vy_min(velocity_block);
         const Real block_vy_max = get_velocity_block_vy_max(velocity_block);

         return block_vy_min + (block_vy_max - block_vy_min) / block_vy_length * indices[1];
      }

      /*!
      Returns the edge where given velocity cell in the given velocity block ends.
      */
      static Real get_velocity_cell_vy_max(
         const unsigned int velocity_block,
         const unsigned int velocity_cell
      ) {
         if (velocity_cell == error_velocity_cell) {
            return std::numeric_limits<Real>::quiet_NaN();
         }

         const velocity_cell_indices_t indices = get_velocity_cell_indices(velocity_cell);
         if (indices[1] == error_velocity_cell_index) {
            return std::numeric_limits<Real>::quiet_NaN();
         }

         const Real block_vy_min = get_velocity_block_vy_min(velocity_block);
         const Real block_vy_max = get_velocity_block_vy_max(velocity_block);

         return block_vy_min + (block_vy_max - block_vy_min) / block_vy_length * (indices[1] + 1);
      }

      /*!
      Returns the edge where given velocity cell in the given velocity block starts.
      */
      static  Real get_velocity_cell_vz_min(
         const unsigned int velocity_block,
         const unsigned int velocity_cell
      ) {
         if (velocity_cell == error_velocity_cell) {
            return std::numeric_limits<Real>::quiet_NaN();
         }

         const velocity_cell_indices_t indices = get_velocity_cell_indices(velocity_cell);
         if (indices[2] == error_velocity_cell_index) {
            return std::numeric_limits<Real>::quiet_NaN();
         }

         const Real block_vz_min = get_velocity_block_vz_min(velocity_block);
         const Real block_vz_max = get_velocity_block_vz_max(velocity_block);

         return block_vz_min + (block_vz_max - block_vz_min) / block_vz_length * indices[2];
      }

      /*!
      Returns the edge where given velocity cell in the given velocity block ends.
      */
      static Real get_velocity_cell_vz_max(
         const unsigned int velocity_block,
         const unsigned int velocity_cell
      ) {
         if (velocity_cell == error_velocity_cell) {
            return std::numeric_limits<Real>::quiet_NaN();
         }

         const velocity_cell_indices_t indices = get_velocity_cell_indices(velocity_cell);
         if (indices[2] == error_velocity_cell_index) {
            return std::numeric_limits<Real>::quiet_NaN();
         }

         const Real block_vz_min = get_velocity_block_vz_min(velocity_block);
         const Real block_vz_max = get_velocity_block_vz_max(velocity_block);

         return block_vz_min + (block_vz_max - block_vz_min) / block_vz_length * (indices[2] + 1);
      }


      SpatialCell()
      {
         /*
           Block list and cache always have room for all blocks
         */
         this->velocity_block_list.reserve(SpatialCell::max_velocity_blocks);
         this->number_of_blocks=0;
         
         this->block_address_cache.reserve(SpatialCell::max_velocity_blocks);
         for (unsigned int block = 0; block < SpatialCell::max_velocity_blocks; block++) {
            this->velocity_block_list.push_back(error_velocity_block);
            this->block_address_cache.push_back(&(this->null_block));
         }
         this->null_block.clear();
         // zero neighbor lists of null block
         for (unsigned int i = 0; i < N_NEIGHBOR_VELOCITY_BLOCKS; i++) {
            this->null_block.neighbors[i] = NULL;
         }

         // reset spatial cell parameters
         for (unsigned int i = 0; i < CellParams::N_SPATIAL_CELL_PARAMS; i++) {
            this->parameters[i]=0;
         }

         // reset spatial cell derivatives
         for (unsigned int i = 0; i < fieldsolver::N_SPATIAL_CELL_DERIVATIVES; i++) {
            this->derivatives[i]=0;
         }
      }
      
      SpatialCell(const SpatialCell& other):
         number_of_blocks(other.number_of_blocks),
         velocity_block_list(other.velocity_block_list),
         velocity_blocks(other.velocity_blocks),
         initialized(other.initialized),
         isGhostCell(other.isGhostCell),
         boundaryFlag(other.boundaryFlag),
         neighbors(other.neighbors){

         for(unsigned int i=0;i< CellParams::N_SPATIAL_CELL_PARAMS;i++){
            parameters[i]=other.parameters[i];
         }
         for(unsigned int i=0;i< fieldsolver::N_SPATIAL_CELL_DERIVATIVES;i++){
            derivatives[i]=other.derivatives[i];
         }
         
         for (unsigned int block = 0; block < SpatialCell::max_velocity_blocks; block++) {
            if( other.block_address_cache[block] == &(other.null_block))
               this->block_address_cache.push_back(&(this->null_block));
            else
               this->block_address_cache.push_back(&(this->velocity_blocks.at(block)));
         }
         this->null_block.clear();
         // zer o         neighbor lists of null block
         for (unsigned int i = 0; i < N_NEIGHBOR_VELOCITY_BLOCKS; i++) {
            this->null_block.neighbors[i] = NULL;
         }         
         
         
      }
      
       
      /*!
        Returns a pointer to the given velocity block or to
        the null block if given velocity block doesn't exist.
      */
      Velocity_Block* at(const unsigned int block)
         {
            if (block == error_velocity_block
                || block >= SpatialCell::max_velocity_blocks) {
               return &(this->null_block);
            } else {
               return this->block_address_cache.at(block);
            }
         }

      /*!
        A const version of the non-const at function.
      */
      Velocity_Block const* at(const unsigned int block) const
         {
            if (block == error_velocity_block
                || block >= SpatialCell::max_velocity_blocks) {
               return &(this->null_block);
            } else {
               return this->block_address_cache.at(block);
            }
         }


      /*!
        Returns the number of given velocity blocks that exist.
      */
      size_t count(const unsigned int block) const
         {
            return this->velocity_blocks.count(block);
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
            if (this->velocity_blocks.count(block) == 0) {
               if (!this->add_velocity_block(block)) {
                  std::cerr << "Couldn't add velocity block " << block << std::endl;
                  abort();
               }
            }
            Velocity_Block* block_ptr = &(this->velocity_blocks.at(block));
            if (block_ptr == NULL) {
               std::cerr << __FILE__ << ":" << __LINE__
                  << " block_ptr == NULL" << std::endl;
               abort();
            }

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
               displacements.push_back((uint8_t*) &(this->number_of_blocks) - (uint8_t*) this);
               block_lengths.push_back(sizeof(unsigned int) * SpatialCell::max_velocity_blocks);
               block_lengths.push_back(sizeof(unsigned int));
            }

            if((SpatialCell::mpi_transfer_type & Transfer::VEL_BLOCK_DATA)!=0){
               displacements.reserve(displacements.size()+this->velocity_blocks.size());
               block_lengths.reserve(block_lengths.size()+this->velocity_blocks.size());

               for(block_index=0;block_index<this->number_of_blocks;block_index++){
// TODO: use cached block addresses
                  displacements.push_back((uint8_t*) this->velocity_blocks.at(this->velocity_block_list[block_index]).data - (uint8_t*) this);
                  block_lengths.push_back(sizeof(Real) * VELOCITY_BLOCK_LENGTH);
               }
            }

            if((SpatialCell::mpi_transfer_type & Transfer::VEL_BLOCK_FLUXES)!=0){
               displacements.reserve(displacements.size()+this->velocity_blocks.size());
               block_lengths.reserve(block_lengths.size()+this->velocity_blocks.size());
               //send velocity block fluxes

               for(block_index=0;block_index<this->number_of_blocks;block_index++){
                  displacements.push_back((uint8_t*) this->velocity_blocks.at(this->velocity_block_list[block_index]).fx - (uint8_t*) this);
#ifdef  SOLVER_KT                 
                  block_lengths.push_back(sizeof(Real) * 3 * SIZE_FLUXS);
#else  //leveque                                        
                  block_lengths.push_back(sizeof(Real) * SIZE_FLUXS);
#endif                        

               }
            }
            
            // send  spatial cell parameters
            if((SpatialCell::mpi_transfer_type & Transfer::CELL_PARAMETERS)!=0){
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
                block_lengths.push_back(sizeof(bool));
            }
            
            
// send velocity block derivatives
            if((SpatialCell::mpi_transfer_type & Transfer::VEL_BLOCK_KT_DERIVATIVES)!=0){
#ifdef  SOLVER_KT
               displacements.reserve(displacements.size()+this->velocity_blocks.size());
               block_lengths.reserve(block_lengths.size()+this->velocity_blocks.size());
               for(block_index=0;block_index<this->number_of_blocks;block_index++){
                  displacements.push_back((uint8_t*) this->velocity_blocks.at(this->velocity_block_list[block_index]).d1x - (uint8_t*) this);
                  block_lengths.push_back(sizeof(Real) * 6 * SIZE_DERIV);
               }
#endif                  
            }

            if((SpatialCell::mpi_transfer_type & Transfer::VEL_BLOCK_PARAMETERS)!=0){
               displacements.reserve(displacements.size()+this->velocity_blocks.size());
               block_lengths.reserve(block_lengths.size()+this->velocity_blocks.size());
               for(block_index=0;block_index<this->number_of_blocks;block_index++){
                  // TODO: use cached block addresses
                  displacements.push_back((uint8_t*) this->velocity_blocks.at(this->velocity_block_list[block_index]).parameters - (uint8_t*) this);
                  block_lengths.push_back(sizeof(Real) * BlockParams::N_VELOCITY_BLOCK_PARAMS);
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
        Returns true if given velocity block has enough of a distribution function.
        Returns false if the value of the distribution function is too low in every
        sense in given block.
        Also returns false if given block doesn't exist or is an error block.
      */
      bool velocity_block_has_contents(
         const unsigned int block
      ) const
      {
            if (block == error_velocity_block
                || this->velocity_blocks.count(block) == 0) {
               return false;
            }

            bool has_content = false;

            Real total = 0;
            const Velocity_Block* block_ptr = &(this->velocity_blocks.at(block));

            for (unsigned int i = 0; i < VELOCITY_BLOCK_LENGTH; i++) {
               total += block_ptr->data[i];
               if (block_ptr->data[i] >= SpatialCell::velocity_block_min_value) {
                  has_content = true;
                  break;
               }
            }

            if (total >= SpatialCell::velocity_block_min_avg_value * VELOCITY_BLOCK_LENGTH) {
               has_content = true;
            }

            return has_content;
         }


      /*!
        Returns the total value of the distribution function within this spatial cell.
      */
      Real get_total_value(void) const
         {
            Real total = 0;

            for (boost::unordered_map<unsigned int, Velocity_Block>::const_iterator
               block = this->velocity_blocks.begin();
               block != this->velocity_blocks.end();
               block++
            ) {
               for (unsigned int i = 0; i < VELOCITY_BLOCK_LENGTH; i++) {
                  total += block->second.data[i];
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
               + n * 2 * VELOCITY_BLOCK_LENGTH * sizeof(Real);
         }


      /*!
      Returns a block at given offsets from given block.

      Returns error_velocity_block if the returned block
      would be outside of the velocity grid or all given
      offsets are 0.
      */
      unsigned int get_velocity_block_from_offsets(
         const unsigned int block,
         const int x_offset, 
         const int y_offset, 
         const int z_offset
      ) {
         if (x_offset == 0  &&  y_offset == 0 && z_offset == 0) {
            return error_velocity_block;
         }

         const velocity_block_indices_t indices = get_velocity_block_indices(block);
         if (indices[0] == error_velocity_block_index) {
            return error_velocity_block;
         }

         const velocity_block_indices_t neighbor_indices = {
            indices[0] + x_offset,
            indices[1] + y_offset,
            indices[2] + z_offset
         };

         if (neighbor_indices[0] == error_velocity_block_index) {
            return error_velocity_block;
         }

         return get_velocity_block(neighbor_indices);
      }


      /*!
        Checks velocity blocks in the velocity block list.
      */
      void check_velocity_block_list(void) const
         {
            for (unsigned int i = 0; i < SpatialCell::max_velocity_blocks; i++) {
               if (this->velocity_block_list[i] == error_velocity_block) {
                  for (unsigned int j = i; j < SpatialCell::max_velocity_blocks; j++) {
                     if (this->velocity_block_list[i] != error_velocity_block) {
                        std::cerr << __FILE__ << ":" << __LINE__
                                  << "Velocity block list has holes"
                                  << std::endl;
                        abort();
                     }
                  }
                  break;
               }

               if (this->velocity_blocks.count(this->velocity_block_list[i]) == 0) {
                  std::cerr << __FILE__ << ":" << __LINE__
                            << " Velocity block " << this->velocity_block_list[i]
                            << " doesn't exist"
                            << std::endl;
                  abort();
               }
            }
         }
      
      inline void check_number_of_blocks(const std::string file,const int line) const {
         unsigned int numBlocks = 0;
         while (numBlocks < SpatialCell::max_velocity_blocks
                && this->velocity_block_list[numBlocks] != error_velocity_block) {
            numBlocks++;
         }
         if(numBlocks!=this->number_of_blocks){
            std::cerr <<  file << ":" << line << "Number of blocks " << this->number_of_blocks << " should be "<< numBlocks<<std::endl;
            abort();
         }
      }
      
      /*!
        Prints velocity blocks in the velocity block list.
      */
      void print_velocity_block_list(void) const
         {
            std::cout << this->velocity_blocks.size() << " blocks: ";
            for (unsigned int i = 0; i < SpatialCell::max_velocity_blocks; i++) {

               if (this->velocity_block_list[i] == error_velocity_block) {
                  // debug
                  for (unsigned int j = i; j < SpatialCell::max_velocity_blocks; j++) {
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
         if (this->velocity_blocks.count(block) == 0) {
            return;
         }

         const Velocity_Block* block_ptr = &(this->velocity_blocks.at(block));

         std::cout << block << " neighbors: ";
         for (unsigned int neighbor = 0; neighbor < N_NEIGHBOR_VELOCITY_BLOCKS; neighbor++) {
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
            for (unsigned int i = 0; i < SpatialCell::max_velocity_blocks; i++) {

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
        All cells in spatial_neighbors are assumed to be neighbors of this cell.
      */
      void adjust_velocity_blocks(const std::vector<SpatialCell*>& spatial_neighbors)
      {
         #ifdef DEBUG
         this->check_velocity_block_list();
         #endif
         this->check_number_of_blocks(__FILE__,__LINE__);
         // don't iterate over blocks created / removed by this function
         std::vector<unsigned int> original_block_list;
         original_block_list.reserve(SpatialCell::max_velocity_blocks);
         for(unsigned int block_index=0;block_index<this->number_of_blocks;block_index++){
            unsigned int block = this->velocity_block_list[block_index];
            original_block_list.push_back(block);
         }
         
         // remove all blocks in this cell without content + without neighbors with content
         // add neighboring blocks in velocity space for blocks with content
         for (std::vector<unsigned int>::const_iterator
              original = original_block_list.begin();
              original != original_block_list.end();
              original++
         ) {
            const bool original_has_content = this->velocity_block_has_contents(*original);
            
            if (original_has_content) {             
               // add missing neighbors in velocity space
               for(int offset_vx=-1;offset_vx<=1;offset_vx++)
               for(int offset_vy=-1;offset_vy<=1;offset_vy++)
               for(int offset_vz=-1;offset_vz<=1;offset_vz++){                  
                  const unsigned int neighbor_block = get_velocity_block_from_offsets(*original, offset_vx,offset_vy,offset_vz);
                  if (neighbor_block == error_velocity_block) {
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

               // remove local block if also no neighbor has content
               bool neighbors_have_content = false;
               
               // velocity space neighbors
               for(int offset_vx=-1;offset_vx<=1;offset_vx++)
               for(int offset_vy=-1;offset_vy<=1;offset_vy++)
               for(int offset_vz=-1;offset_vz<=1;offset_vz++) {
                  const unsigned int neighbor_block = get_velocity_block_from_offsets(*original, offset_vx,offset_vy,offset_vz);
                  if (this->velocity_block_has_contents(neighbor_block)) {
                     neighbors_have_content = true;
                     break;
                  }
               }
               
               // real space neighbors
               for (std::vector<SpatialCell*>::const_iterator
                  neighbor = spatial_neighbors.begin();
                  neighbor != spatial_neighbors.end();
                  neighbor++
               ) {
                  if ((*neighbor)->velocity_block_has_contents(*original)) {
                     neighbors_have_content = true;
                     break;
                  }
               }
               
               if (!neighbors_have_content) {
                  this->remove_velocity_block(*original);
               }
               
            }
            
         }

         // add local blocks for neighbors in real space with content
         for (std::vector<SpatialCell*>::const_iterator
            neighbor = spatial_neighbors.begin();
            neighbor != spatial_neighbors.end();
            neighbor++
              ) {
            for(unsigned int block_index=0;block_index<(*neighbor)->number_of_blocks;block_index++){
               unsigned int block = (*neighbor)->velocity_block_list[block_index];
               if ((*neighbor)->velocity_block_has_contents(block)) {
                  this->add_velocity_block(block);
               }
            }
         }
         
      }
      
      /*!
        Saves this spatial cell in vtk ascii format into the given filename.
      */
      void save_vtk(const char* filename) const {

         // do nothing if one cell or block dimension is 0
         if (block_vx_length == 0
             || block_vy_length == 0
             || block_vz_length == 0
             || SpatialCell::vx_length == 0
             || SpatialCell::vy_length == 0
             || SpatialCell::vz_length == 0) {
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
         outfile << "POINTS " << (this->velocity_blocks.size() * VELOCITY_BLOCK_LENGTH + 2) * 8 << " double" << std::endl;
         for (std::vector<unsigned int>::const_iterator
            block = this->velocity_block_list.begin();
            block != this->velocity_block_list.end();
            block++
         ) {

            if (*block == error_velocity_block) {
               // assume no blocks after first error block
               break;
            }

            for (unsigned int z_index = 0; z_index < block_vz_length; z_index++)
               for (unsigned int y_index = 0; y_index < block_vy_length; y_index++)
                  for (unsigned int x_index = 0; x_index < block_vx_length; x_index++) {

                     const velocity_cell_indices_t indices = {x_index, y_index, z_index};
                     const unsigned int velocity_cell = get_velocity_cell(indices);

                     outfile << get_velocity_cell_vx_min(*block, velocity_cell) << " "
                             << get_velocity_cell_vy_min(*block, velocity_cell) << " "
                             << get_velocity_cell_vz_min(*block, velocity_cell) << "\n";

                     outfile << get_velocity_cell_vx_max(*block, velocity_cell) << " "
                             << get_velocity_cell_vy_min(*block, velocity_cell) << " "
                             << get_velocity_cell_vz_min(*block, velocity_cell) << "\n";

                     outfile << get_velocity_cell_vx_min(*block, velocity_cell) << " "
                             << get_velocity_cell_vy_max(*block, velocity_cell) << " "
                             << get_velocity_cell_vz_min(*block, velocity_cell) << "\n";

                     outfile << get_velocity_cell_vx_max(*block, velocity_cell) << " "
                             << get_velocity_cell_vy_max(*block, velocity_cell) << " "
                             << get_velocity_cell_vz_min(*block, velocity_cell) << "\n";

                     outfile << get_velocity_cell_vx_min(*block, velocity_cell) << " "
                             << get_velocity_cell_vy_min(*block, velocity_cell) << " "
                             << get_velocity_cell_vz_max(*block, velocity_cell) << "\n";

                     outfile << get_velocity_cell_vx_max(*block, velocity_cell) << " "
                             << get_velocity_cell_vy_min(*block, velocity_cell) << " "
                             << get_velocity_cell_vz_max(*block, velocity_cell) << "\n";
 
                     outfile << get_velocity_cell_vx_min(*block, velocity_cell) << " "
                             << get_velocity_cell_vy_max(*block, velocity_cell) << " "
                             << get_velocity_cell_vz_max(*block, velocity_cell) << "\n";

                     outfile << get_velocity_cell_vx_max(*block, velocity_cell) << " "
                             << get_velocity_cell_vy_max(*block, velocity_cell) << " "
                             << get_velocity_cell_vz_max(*block, velocity_cell) << "\n";
                  }
         }
         /*
           Add small velocity cells to the negative and positive corners of the grid
           so VisIt knows the maximum size of the velocity grid regardless of existing cells
         */
         outfile << SpatialCell::vx_min - 0.1 * SpatialCell::cell_dvx << " "
                 << SpatialCell::vy_min - 0.1 * SpatialCell::cell_dvy << " "
                 << SpatialCell::vz_min - 0.1 * SpatialCell::cell_dvz << "\n";
         outfile << SpatialCell::vx_min << " "
                 << SpatialCell::vy_min - 0.1 * SpatialCell::cell_dvy << " "
                 << SpatialCell::vz_min - 0.1 * SpatialCell::cell_dvz << "\n";
         outfile << SpatialCell::vx_min - 0.1 * SpatialCell::cell_dvx << " "
                 << SpatialCell::vy_min << " "
                 << SpatialCell::vz_min - 0.1 * SpatialCell::cell_dvz << "\n";
         outfile << SpatialCell::vx_min << " "
                 << SpatialCell::vy_min << " "
                 << SpatialCell::vz_min - 0.1 * SpatialCell::cell_dvz << "\n";
         outfile << SpatialCell::vx_min - 0.1 * SpatialCell::cell_dvx << " "
                 << SpatialCell::vy_min - 0.1 * SpatialCell::cell_dvy << " "
                 << SpatialCell::vz_min << "\n";
         outfile << SpatialCell::vx_min << " "
                 << SpatialCell::vy_min - 0.1 * SpatialCell::cell_dvy << " "
                 << SpatialCell::vz_min << "\n";
         outfile << SpatialCell::vx_min - 0.1 * SpatialCell::cell_dvx << " "
                 << SpatialCell::vy_min << " "
                 << SpatialCell::vz_min << "\n";
         outfile << SpatialCell::vx_min << " "
                 << SpatialCell::vy_min << " "
                 << SpatialCell::vz_min << "\n";

         outfile << SpatialCell::vx_max << " "
                 << SpatialCell::vy_max << " "
                 << SpatialCell::vz_max << "\n";
         outfile << SpatialCell::vx_max + 0.1 * SpatialCell::cell_dvx << " "
                 << SpatialCell::vy_max << " "
                 << SpatialCell::vz_max << "\n";
         outfile << SpatialCell::vx_max << " "
                 << SpatialCell::vy_max + 0.1 * SpatialCell::cell_dvy << " "
                 << SpatialCell::vz_max << "\n";
         outfile << SpatialCell::vx_max + 0.1 * SpatialCell::cell_dvx << " "
                 << SpatialCell::vy_max + 0.1 * SpatialCell::cell_dvy << " "
                 << SpatialCell::vz_max << "\n";
         outfile << SpatialCell::vx_max << " "
                 << SpatialCell::vy_max << " "
                 << SpatialCell::vz_max + 0.1 * SpatialCell::cell_dvz << "\n";
         outfile << SpatialCell::vx_max + 0.1 * SpatialCell::cell_dvx << " "
                 << SpatialCell::vy_max << " "
                 << SpatialCell::vz_max + 0.1 * SpatialCell::cell_dvz << "\n";
         outfile << SpatialCell::vx_max << " "
                 << SpatialCell::vy_max + 0.1 * SpatialCell::cell_dvy << " "
                 << SpatialCell::vz_max + 0.1 * SpatialCell::cell_dvz << "\n";
         outfile << SpatialCell::vx_max + 0.1 * SpatialCell::cell_dvx << " "
                 << SpatialCell::vy_max + 0.1 * SpatialCell::cell_dvy << " "
                 << SpatialCell::vz_max + 0.1 * SpatialCell::cell_dvz << "\n";

         // map cells to written points
         outfile << "CELLS "
                 << this->velocity_blocks.size() * VELOCITY_BLOCK_LENGTH + 2 << " "
                 << (this->velocity_blocks.size() * VELOCITY_BLOCK_LENGTH + 2)* 9 << std::endl;

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

            for (unsigned int z_index = 0; z_index < block_vz_length; z_index++)
               for (unsigned int y_index = 0; y_index < block_vy_length; y_index++)
                  for (unsigned int x_index = 0; x_index < block_vx_length; x_index++) {

                     outfile << "8 ";
                     for (int i = 0; i < 8; i++) {
                        outfile << j * 8 + i << " ";
                     }
                     outfile << "\n";

                     j++;
                  }
         }
         outfile << "8 ";
         for (unsigned int i = 0; i < 8; i++) {
            outfile << j * 8 + i << " ";
         }
         j++;
         outfile << "\n 8 ";
         for (unsigned int i = 0; i < 8; i++) {
            outfile << j * 8 + i << " ";
         }
         outfile << "\n";

         // cell types
         outfile << "CELL_TYPES " << this->velocity_blocks.size() * VELOCITY_BLOCK_LENGTH + 2 << std::endl;
         for (unsigned int i = 0; i < this->velocity_blocks.size() * VELOCITY_BLOCK_LENGTH + 2; i++) {
            outfile << "11\n";
         }

         // Put minimum value from existing blocks into two additional cells
         Real min_value = std::numeric_limits<Real>::max();

         // distribution function
         outfile << "CELL_DATA " << this->velocity_blocks.size() * VELOCITY_BLOCK_LENGTH + 2 << std::endl;
         outfile << "SCALARS rho double 1" << std::endl;
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

            const Velocity_Block block_data = this->velocity_blocks.at(*block);

            for (unsigned int z_index = 0; z_index < block_vz_length; z_index++)
               for (unsigned int y_index = 0; y_index < block_vy_length; y_index++)
                  for (unsigned int x_index = 0; x_index < block_vx_length; x_index++) {

                     const velocity_cell_indices_t indices = {x_index, y_index, z_index};
                     const unsigned int velocity_cell = get_velocity_cell(indices);
                     const Real value = block_data.data[velocity_cell];
                     outfile << value << " ";
                     min_value = std::min(min_value, value);
                  }
            outfile << "\n";
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

         if (block >= SpatialCell::max_velocity_blocks) {
            return false;
         }

         if (this->velocity_blocks.count(block) > 0) {
            return true;
         }

         this->velocity_blocks[block];
         if (this->velocity_blocks.count(block) == 0) {
         }
         this->block_address_cache[block] = &(this->velocity_blocks.at(block));

         Velocity_Block* block_ptr = this->block_address_cache[block];

	if (block_ptr == NULL) {
	   std::cerr << __FILE__ << ":" << __LINE__
	      << " Block pointer == NULL" << std::endl;
	   abort();
	}

         block_ptr->clear();

         // set block parameters
         block_ptr->parameters[BlockParams::VXCRD] = get_velocity_block_vx_min(block);
         block_ptr->parameters[BlockParams::VYCRD] = get_velocity_block_vy_min(block);
         block_ptr->parameters[BlockParams::VZCRD] = get_velocity_block_vz_min(block);
         block_ptr->parameters[BlockParams::DVX] = SpatialCell::cell_dvx;
         block_ptr->parameters[BlockParams::DVY] = SpatialCell::cell_dvy;
         block_ptr->parameters[BlockParams::DVZ] = SpatialCell::cell_dvz;

         // set neighbour pointers
         unsigned int neighbor_block;
         int neighbor_index=0;
         
         for(int offset_vx=-1;offset_vx<=1;offset_vx++)
         for(int offset_vy=-1;offset_vy<=1;offset_vy++)
         for(int offset_vz=-1;offset_vz<=1;offset_vz++){

            // location of neighbor pointer in this block's neighbor list
            neighbor_index=(1+offset_vx)+(1+offset_vy)*3+(1+offset_vz)*9;
            if(offset_vx==0 && offset_vy==0 && offset_vz==0) continue;

            neighbor_block = get_velocity_block_from_offsets(block, offset_vx,offset_vy,offset_vz);
            if (neighbor_block == error_velocity_block) {
               block_ptr->neighbors[neighbor_index] = NULL;
            }
            else if (this->velocity_blocks.count(neighbor_block) == 0) {
               block_ptr->neighbors[neighbor_index] = &(this->null_block);
            } else {
               block_ptr->neighbors[neighbor_index] = &(this->velocity_blocks.at(neighbor_block));
               //update the neighbor list of neighboring block
               //index of current block in neighbors neighbor table
               //FIXME: check this properly, or write it in a better way
               int neighbor_neighbor_index=(1-offset_vx)+(1-offset_vy)*3+(1-offset_vz)*9;
               Velocity_Block* neighbor_ptr = &(this->velocity_blocks.at(neighbor_block));
               if (neighbor_ptr == NULL) {
                  std::cerr << __FILE__ << ":" << __LINE__
                            << " Block pointer == NULL" << std::endl;
                  abort();
               }
               neighbor_ptr->neighbors[neighbor_neighbor_index] = block_ptr;
            }
         }

         this->velocity_block_list[this->number_of_blocks] = block;
         //add stopper value  (deprecated, use number_of_blocks)
         if(this->number_of_blocks<SpatialCell::max_velocity_blocks)
            this->velocity_block_list[this->number_of_blocks+1] = error_velocity_block;
         this->number_of_blocks++;


         return true;
      }
      
      /*!
        Adds all velocity blocks than don't exist into the velocity grid.

        Returns true if all non-existing blocks were added, false otherwise.
      */
      bool add_all_velocity_blocks(void)
         {
            bool result = true;
            for (unsigned int i = 0; i < SpatialCell::max_velocity_blocks; i++) {
               if (this->velocity_blocks.count(i) > 0) {
                  continue;
               }

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
      void remove_velocity_block(const unsigned int block)
      {
         if (block == error_velocity_block) {
            return;
         }

         if (this->velocity_blocks.count(block) == 0) {
            return;
         }

         //remove block from neighbors neighbor lists
         unsigned int neighbor_block;
         for(int offset_vx=-1;offset_vx<=1;offset_vx++)
         for(int offset_vy=-1;offset_vy<=1;offset_vy++)
         for(int offset_vz=-1;offset_vz<=1;offset_vz++) {

            if(offset_vx==0 && offset_vy==0 && offset_vz==0) continue;

            neighbor_block = get_velocity_block_from_offsets(block, offset_vx,offset_vy,offset_vz);
            if (neighbor_block != error_velocity_block
                && this->velocity_blocks.count(neighbor_block) > 0) {
               // TODO use cached addresses of neighbors
               Velocity_Block* neighbor_data = &(this->velocity_blocks.at(neighbor_block));
               // update the neighbor list of neighboring block
               // index of current block in neighbors neighbor table
               int neighbor_neighbor_index=(1-offset_vx)+(1-offset_vy)*3+(1-offset_vz)*9;
               neighbor_data->neighbors[neighbor_neighbor_index] = &(this->null_block);
            }
         }

         this->velocity_blocks.erase(block);
         this->block_address_cache[block] = &(this->null_block);

         /*
           Move the last existing block in the block list
           to the removed block's position
         */
         unsigned int block_index = 0;
         while (block_index < SpatialCell::max_velocity_blocks
                && this->velocity_block_list[block_index] != block) {
            block_index++;
         }
         //debug
         if (block_index == SpatialCell::max_velocity_blocks) {
            std::cerr << "Velocity block " << block << " not found in list" << std::endl;
            abort();
         }


         
         if (block_index == this->number_of_blocks - 1) {
            this->velocity_block_list[block_index] = error_velocity_block;
         } else {
            this->velocity_block_list[block_index] = this->velocity_block_list[this->number_of_blocks - 1];
            this->velocity_block_list[this->number_of_blocks - 1] = error_velocity_block;
         }
         this->number_of_blocks--;

      }


      /*!
        Removes all velocity blocks from this spatial cell.
      */
      void clear(void)
         {
            this->velocity_blocks.clear();
            //FIXME, it shold not be necessary to fill all list values with error_velocity_block, only first.
            for (unsigned int i = 0; i < SpatialCell::max_velocity_blocks; i++) {
               this->block_address_cache[i] = &(this->null_block);
               this->velocity_block_list[i] = error_velocity_block;
            }
            this->number_of_blocks=0;
         }


      /*!
        Prepares this spatial cell to receive the velocity grid over MPI.

        At this stage we have received a new blocklist over MPI, but
        the rest of the cell structres have ntbee adapted to this new
        list. Here we re-initialize the cell wit empty blocks based on
        the new list
      */
      void prepare_to_receive_blocks(void)
         {
            unsigned int oldNumBlocks = this->number_of_blocks;
            // clear + add_velocity_block overwrites the block list so:
            std::vector<unsigned int> old_block_list(oldNumBlocks, error_velocity_block);
            for (unsigned int block_index = 0; block_index < oldNumBlocks; block_index++) {
               old_block_list[block_index] = this->velocity_block_list[block_index];
            }

            this->clear();

            // add velocity blocks that are about to be received with MPI
            for (unsigned int block_index = 0; block_index < oldNumBlocks; block_index++) {
               this->add_velocity_block(old_block_list[block_index]);
            }


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

      /*
      Primary velocity grid parameters
      */

      // size of the velocity grid in velocity blocks
      static unsigned int vx_length, vy_length, vz_length;

      // physical size of the velocity grid
      static Real vx_min, vx_max, vy_min, vy_max, vz_min, vz_max;

      /*
      Derived velocity grid parameters
      */
      static unsigned int max_velocity_blocks;

      // physical size of velocity grid
      static Real grid_dvx, grid_dvy, grid_dvz;

      // physical size of velocity blocks
      static Real block_dvx, block_dvy, block_dvz;

      // physical size of velocity cells
      static Real cell_dvx, cell_dvy, cell_dvz;

      /*
        Which data is transferred by the mpi datatype given by spatial cells.
      */
      static unsigned int mpi_transfer_type;

      /*  
        Minimum value of distribution function
        in any cell of a velocity block for the
        block to be considered to have contents
      */
      static Real velocity_block_min_value;

      /*
        Minimum value of the average of distribution
        function within a velocity block for the
        block to be considered to have contents
      */
      static Real velocity_block_min_avg_value;

   private:

      SpatialCell & operator= (const SpatialCell&);
      
      
      /*
      Per-instance parameters
      */

      bool initialized;
      /*!
        Used as a neighbour instead of blocks that don't
        exist but would be inside of the velocity grid.
        Neighbors that would be outside of the grid are always NULL.
      */
      Velocity_Block null_block;


      // data of velocity blocks that exist in this cell
      boost::unordered_map<unsigned int, Velocity_Block> velocity_blocks;

      /*
        Speed up search of velocity block addresses in the hash table above.
        Addresses of non-existing blocks point to the null block.
      */

      std::vector<Velocity_Block*> block_address_cache;



   public:
      /*
        List of velocity blocks in this cell, used for
        transferring spatial cells between processes using MPI.
      */
      
      std::vector<unsigned int> velocity_block_list;
      unsigned int number_of_blocks;

      /*
        Bulk variables in this spatial cell.
      */
      Real parameters[CellParams::N_SPATIAL_CELL_PARAMS];

      /*
        Derivatives of bulk variables in this spatial cell.
      */
      Real derivatives[fieldsolver::N_SPATIAL_CELL_DERIVATIVES];
      
      //neighbor id's. Kept up to date in solvers, not by the spatial_cell class
      std::vector<uint64_t> neighbors;
      unsigned int boundaryFlag;
      bool isGhostCell;
      
      
   }; // class SpatialCell
   

} // namespaces
#endif

