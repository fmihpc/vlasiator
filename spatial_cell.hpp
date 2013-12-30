
/*!
Spatial cell class for Vlasiator that supports a variable number of velocity blocks.

Copyright 2011 Finnish Meteorological Institute
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
#include "set"

#include "memoryallocation.h"

#include "phiprof.hpp"
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
      const uint64_t NONE                     = 0;
      const uint64_t CELL_PARAMETERS          = (1<<0);
      const uint64_t CELL_DERIVATIVES         = (1<<1);
      const uint64_t VEL_BLOCK_LIST_STAGE1    = (1<<2);
      const uint64_t VEL_BLOCK_LIST_STAGE2    = (1<<3);
      const uint64_t VEL_BLOCK_DATA           = (1<<4);
      const uint64_t VEL_BLOCK_FLUXES         = (1<<5);
      const uint64_t VEL_BLOCK_PARAMETERS     = (1<<6);
      const uint64_t VEL_BLOCK_WITH_CONTENT_STAGE1  = (1<<7); 
      const uint64_t VEL_BLOCK_WITH_CONTENT_STAGE2  = (1<<8); 
      const uint64_t CELL_SYSBOUNDARYFLAG     = (1<<9);
      const uint64_t CELL_E                   = (1<<10);
      const uint64_t CELL_EDT2                = (1<<11);
      const uint64_t CELL_PERB                = (1<<12);
      const uint64_t CELL_PERBDT2             = (1<<13);
      const uint64_t CELL_BGB                 = (1<<14);
      const uint64_t CELL_RHO_RHOV            = (1<<15);
      const uint64_t CELL_RHODT2_RHOVDT2      = (1<<16);
      const uint64_t CELL_BVOL                = (1<<17);
      const uint64_t CELL_BVOL_DERIVATIVES    = (1<<18);
      const uint64_t CELL_DIMENSIONS          = (1<<19);
      const uint64_t CELL_IOLOCALCELLID       = (1<<20);
      
      const uint64_t ALL_DATA =
      CELL_PARAMETERS
      | CELL_DERIVATIVES | CELL_BVOL_DERIVATIVES
      | VEL_BLOCK_DATA
      | VEL_BLOCK_FLUXES
      | CELL_SYSBOUNDARYFLAG;
   }

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


   
   class Velocity_Block {
   public:
      // value of the distribution function
      Real *data;   
      // spatial fluxes of this block
      //fixme, fx could be called flux for leveque
      Real *fx;
      
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
            for (unsigned int i = 0; i < SIZE_FLUXS; i++) {
               this->fx[i] = 0;
            }
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

         const velocity_block_indices_t indices = {{
            (unsigned int) floor((vx - SpatialCell::vx_min) / SpatialCell::block_dvx),
            (unsigned int) floor((vy - SpatialCell::vy_min) / SpatialCell::block_dvy),
            (unsigned int) floor((vz - SpatialCell::vz_min) / SpatialCell::block_dvz)
         }};

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

         const velocity_block_indices_t indices = {{
            (unsigned int) floor((vx - block_vx_min) / ((block_vx_max - block_vx_min) / block_vx_length)),
            (unsigned int) floor((vy - block_vy_min) / ((block_vy_max - block_vy_min) / block_vy_length)),
            (unsigned int) floor((vz - block_vz_min) / ((block_vz_max - block_vz_min) / block_vz_length))
         }};

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
         this->number_of_blocks=0;
         this->mpi_number_of_blocks=0;

         //allocate memory for null block
         this->null_block_data.resize(VELOCITY_BLOCK_LENGTH);
         this->null_block_fx.resize(SIZE_FLUXS);
         this->null_block.data=&(this->null_block_data[0]);
         this->null_block.fx=&(this->null_block_fx[0]);
         this->sysBoundaryLayer=0; /*!< Default value, layer not yet initialized*/
         this->null_block.clear();

         // zero neighbor lists of null block
         for (unsigned int i = 0; i < N_NEIGHBOR_VELOCITY_BLOCKS; i++) {
            this->null_block.neighbors[i] = NULL;
         }

         // reset spatial cell parameters
         for (unsigned int i = 0; i < CellParams::N_SPATIAL_CELL_PARAMS; i++) {
            this->parameters[i]=0.0;
         }

         // reset spatial cell derivatives
         for (unsigned int i = 0; i < fieldsolver::N_SPATIAL_CELL_DERIVATIVES; i++) {
            this->derivatives[i]=0;
         }
         
         // reset BVOL derivatives
         for (unsigned int i = 0; i < bvolderivatives::N_BVOL_DERIVATIVES; i++) {
            this->derivativesBVOL[i]=0;
         }
	 //is transferred by default
	 this->mpiTransferEnabled=true;

      }
      
      SpatialCell(const SpatialCell& other):
         initialized(other.initialized),
         velocity_blocks(other.velocity_blocks),
	 mpiTransferEnabled(other.mpiTransferEnabled),
         number_of_blocks(other.number_of_blocks),
         velocity_block_list(other.velocity_block_list),
         mpi_number_of_blocks(other.mpi_number_of_blocks),
         mpi_velocity_block_list(other.mpi_velocity_block_list),
         velocity_block_with_content_list(other.velocity_block_with_content_list),
         velocity_block_with_no_content_list(other.velocity_block_with_no_content_list),
         block_data(other.block_data),
         block_fx(other.block_fx),
         null_block_data(other.null_block_data),
         null_block_fx(other.null_block_fx),
         neighbors(other.neighbors),
         procBoundaryFlag(other.procBoundaryFlag),
         sysBoundaryFlag(other.sysBoundaryFlag),
         sysBoundaryLayer(other.sysBoundaryLayer)
         {

//       phiprof::initializeTimer("SpatialCell copy", "SpatialCell copy");
//       phiprof::start("SpatialCell copy");

         //copy parameters
         for(unsigned int i=0;i< CellParams::N_SPATIAL_CELL_PARAMS;i++){
            parameters[i]=other.parameters[i];
         }
         //copy derivatives
         for(unsigned int i=0;i< fieldsolver::N_SPATIAL_CELL_DERIVATIVES;i++){
            derivatives[i]=other.derivatives[i];
         }
         //copy BVOL derivatives
         for(unsigned int i=0;i< bvolderivatives::N_BVOL_DERIVATIVES;i++){
            derivativesBVOL[i]=other.derivativesBVOL[i];
         }

         
         // zero neighbor lists of null block
         for (unsigned int i = 0; i < N_NEIGHBOR_VELOCITY_BLOCKS; i++) {
            this->null_block.neighbors[i] = NULL;
         }
         //set null block data
         this->null_block.data=&(this->null_block_data[0]);
         this->null_block.fx=&(this->null_block_fx[0]);
         this->null_block.clear();


         for(unsigned int block_i=0;block_i<number_of_blocks;block_i++){
            unsigned int block=velocity_block_list[block_i];
            Velocity_Block* block_ptr = this->at(block);
            //fix block data pointers   
            set_block_data_pointers(block_i);
            
//fix neighbor lists of all normal blocks
            //loop through neighbors
            for(int offset_vx=-1;offset_vx<=1;offset_vx++)
               for(int offset_vy=-1;offset_vy<=1;offset_vy++)
                  for(int offset_vz=-1;offset_vz<=1;offset_vz++){
                     // location of neighbor pointer in this block's neighbor list
                     int neighbor_index=(1+offset_vx)+(1+offset_vy)*3+(1+offset_vz)*9;
                     if(offset_vx==0 && offset_vy==0 && offset_vz==0)
                        continue;

                     unsigned int neighbor_block = get_velocity_block_from_offsets(block, offset_vx,offset_vy,offset_vz);
                     //error, make neighbor null
                     if (neighbor_block == error_velocity_block) {
                        block_ptr->neighbors[neighbor_index] = NULL;
                     }
                     //neighbor not in cell, add null_block as neighbour
                     else if (this->velocity_blocks.count(neighbor_block) == 0) {
                        block_ptr->neighbors[neighbor_index] = &(this->null_block);
                     }
                     //add neighbor pointer
                     else {
                        block_ptr->neighbors[neighbor_index] = &(this->velocity_blocks.at(neighbor_block));
                     }
                  }
         }
//         phiprof::stop("SpatialCell copy");
      }
      
      /*!
        Returns a pointer to the given velocity block. Fast version, no error checking. Need to check separately for existence.
      */
      Velocity_Block* at_fast(const unsigned int block){
	return &(this->velocity_blocks.at(block));
      }
      
      /*!
        Returns a pointer to the given velocity block or to
        the null block if given velocity block doesn't exist.
      */
      Velocity_Block* at(const unsigned int block){
         // check if block cannot exist, non-allowed block index
         if(block==error_velocity_block || block >= SpatialCell::max_velocity_blocks) {
            return &(this->null_block);
         }
         
         boost::unordered_map<unsigned int, Velocity_Block>::iterator iter=this->velocity_blocks.find(block);
         //check if it exists in map
         if(iter == this->velocity_blocks.end()){
            return &(this->null_block);
         }
         else {
            return &(iter->second);
         }
      }

      /*!
        A const version of the non-const at function.
      */
      Velocity_Block const* at(const unsigned int block) const {
         //  check if block cannot exist, non-allowed block inde        x
         if(block==error_velocity_block || block >= SpatialCell::max_velocity_blocks) {
            return &(this->null_block);
         }
         
         boost::unordered_map<unsigned int, Velocity_Block>::const_iterator iter=this->velocity_blocks.find(block);
         //check if it exists in map
         if(iter == this->velocity_blocks.end()){
            return &(this->null_block);
         }
         else {
            return &(iter->second);
         }
      }
      

      bool is_null_block(Velocity_Block* block_ptr) const
         {
            return (block_ptr== &(this->null_block)|| block_ptr==NULL);
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
      void set_value(const Real vx, const Real vy, const Real vz, const Real value) {
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

//TODO - thread safe set/increment functions which do not create blocks automatically

      /*! Sets the value of a particular cell in a block. The block is
       *  created if it does not exist. This version is faster than
       *  the velocity value based version.
       *
       * This function is not thread safe due to the creation of
       * blocks.

        \param block Block index of velocity-cell
        \param cell  Cell index (0..WID3-1) of velocity-cell in block
        \param value Value that is set for velocity-cell
      */
     void set_value(const unsigned int block,const unsigned int cell, const Real value) {
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
            block_ptr->data[cell] = value;
      }
      
          
      /*!
        Increments the value of velocity cell at given coordinate-

        Creates the velocity block at given coordinates if it doesn't exist.
      */
      void increment_value(const Real vx, const Real vy, const Real vz, const Real value) {
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
            block_ptr->data[cell] += value;
      }

      /*!
        Increments the value of velocity cell at given index

        Creates the velocity block if it doesn't exist.
      */
     void increment_value(const unsigned int block,const unsigned int cell, const Real value) {
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
            block_ptr->data[cell] += value;
      }
      

      /*!
       * Gets the value of a velocity cell at given coordinates.
       * 
       * Returns 0 if it doesn't exist.
       */
      Real get_value(const Real vx, const Real vy, const Real vz) const
      {
         const unsigned int block = get_velocity_block(vx, vy, vz);
         if (this->velocity_blocks.count(block) == 0) {
            return 0.0;
         }
         const Velocity_Block* const block_ptr = &(this->velocity_blocks.at(block));
         if (block_ptr == NULL) {
            std::cerr << __FILE__ << ":" << __LINE__
            << " block_ptr == NULL" << std::endl; 
            abort();
         }
         
         const unsigned int cell = get_velocity_cell(block, vx, vy, vz);
         return block_ptr->data[cell];
      }
      
      void mpi_datatype(
         void*& address,
         int& count,
         MPI_Datatype& datatype,
         const uint64_t /*cell_id*/,
         const int /*sender*/,
         const int /*receiver*/,
         const bool receiving
      ) {
            address = this;

            std::vector<MPI_Aint> displacements;
            std::vector<int> block_lengths;
            unsigned int block_index = 0;

            /*create daatype for actual data if we are in the first two layers around a boundary, or if we send for thw whole system*/
            if(this->mpiTransferEnabled &&
	       (SpatialCell::mpiTransferAtSysBoundaries==false || this->sysBoundaryLayer ==1 || this->sysBoundaryLayer ==2 )
	       ) {
	      
               //add data to send/recv to displacement and block length lists
               if((SpatialCell::mpi_transfer_type & Transfer::VEL_BLOCK_LIST_STAGE1)!=0){
                  //first copy values in case this is the send operation
                  this->mpi_number_of_blocks=this->number_of_blocks;
                  // send velocity block list size
                  displacements.push_back((uint8_t*) &(this->mpi_number_of_blocks) - (uint8_t*) this);
                  block_lengths.push_back(sizeof(unsigned int));
               }
            
               if((SpatialCell::mpi_transfer_type & Transfer::VEL_BLOCK_LIST_STAGE2)!=0){
                  // STAGE1 should have been done, otherwise we have problems...
                  if(receiving) {
                     //mpi_number_of_blocks transferred earlier
                     this->mpi_velocity_block_list.resize(this->mpi_number_of_blocks);
                  }
                  else {
                     //resize to correct size (it will avoid reallociation if it is big enough, I assume)
                     this->mpi_number_of_blocks=this->number_of_blocks;
                     this->mpi_velocity_block_list.resize(this->number_of_blocks);
                     //copy values if this is the send operation                     
                     for(unsigned int i=0;i< this->number_of_blocks;i++)
                        this->mpi_velocity_block_list[i]=this->velocity_block_list[i];
                     
                  }
                  // send velocity block list
                  displacements.push_back((uint8_t*) &(this->mpi_velocity_block_list[0]) - (uint8_t*) this);
                  block_lengths.push_back(sizeof(unsigned int) * this->mpi_number_of_blocks);
               }

               if((SpatialCell::mpi_transfer_type & Transfer::VEL_BLOCK_WITH_CONTENT_STAGE1)!=0){
                  //Communicate size of list so that buffers can be allocated on receiving side
                  if(!receiving)
                     this->velocity_block_with_content_list_size=this->velocity_block_with_content_list.size();
                  displacements.push_back((uint8_t*) &(this->velocity_block_with_content_list_size) - (uint8_t*) this);
                  block_lengths.push_back(sizeof(unsigned int));
               }
               if((SpatialCell::mpi_transfer_type & Transfer::VEL_BLOCK_WITH_CONTENT_STAGE2)!=0){
                  if(receiving) {
                     this->velocity_block_with_content_list.resize(this->velocity_block_with_content_list_size);
                  }
                  //velocity_block_with_content_list_size should first be updated, before this can be done (STAGE1)
                  displacements.push_back((uint8_t*) &(this->velocity_block_with_content_list[0]) - (uint8_t*) this);
                  block_lengths.push_back(sizeof(unsigned int)*this->velocity_block_with_content_list_size);
               }

               if((SpatialCell::mpi_transfer_type & Transfer::VEL_BLOCK_DATA)!=0){
                  displacements.push_back((uint8_t*) &(this->block_data[0]) - (uint8_t*) this);               
                  block_lengths.push_back(sizeof(Real) * VELOCITY_BLOCK_LENGTH* this->number_of_blocks);
               }
               
               if((SpatialCell::mpi_transfer_type & Transfer::VEL_BLOCK_FLUXES)!=0){
                  displacements.push_back((uint8_t*) &(this->block_fx[0]) - (uint8_t*) this);               
                  block_lengths.push_back(sizeof(Real) * SIZE_FLUXS* this->number_of_blocks);
               }
               
               // send  spatial cell parameters
               if((SpatialCell::mpi_transfer_type & Transfer::CELL_PARAMETERS)!=0){
                  displacements.push_back((uint8_t*) &(this->parameters[0]) - (uint8_t*) this);
                  block_lengths.push_back(sizeof(Real) * CellParams::N_SPATIAL_CELL_PARAMS);
               }
               
               // send  spatial cell dimensions
               if((SpatialCell::mpi_transfer_type & Transfer::CELL_DIMENSIONS)!=0){
                  displacements.push_back((uint8_t*) &(this->parameters[CellParams::DX]) - (uint8_t*) this);
                  block_lengths.push_back(sizeof(Real) * 3);
               }
               
               // send  BGBX BGBY BGBZ
               if((SpatialCell::mpi_transfer_type & Transfer::CELL_BGB)!=0){
                  displacements.push_back((uint8_t*) &(this->parameters[CellParams::BGBX]) - (uint8_t*) this);
                  block_lengths.push_back(sizeof(Real) * 3);
               }
            
               // send  BXVOL BYVOL BZVOL
               if((SpatialCell::mpi_transfer_type & Transfer::CELL_BVOL)!=0){
                  displacements.push_back((uint8_t*) &(this->parameters[CellParams::BGBXVOL]) - (uint8_t*) this);
                  block_lengths.push_back(sizeof(Real) * 6);
               }
            
               // send  EX, EY EZ
               if((SpatialCell::mpi_transfer_type & Transfer::CELL_E)!=0){
                  displacements.push_back((uint8_t*) &(this->parameters[CellParams::EX]) - (uint8_t*) this);
                  block_lengths.push_back(sizeof(Real) * 3);
               }
            
               // send  EX_DT2, EY_DT2, EZ_DT2
               if((SpatialCell::mpi_transfer_type & Transfer::CELL_EDT2)!=0){
                  displacements.push_back((uint8_t*) &(this->parameters[CellParams::EX_DT2]) - (uint8_t*) this);
                  block_lengths.push_back(sizeof(Real) * 3);
               }
            
               // send  PERBX, PERBY, PERBZ
               if((SpatialCell::mpi_transfer_type & Transfer::CELL_PERB)!=0){
                  displacements.push_back((uint8_t*) &(this->parameters[CellParams::PERBX]) - (uint8_t*) this);
                  block_lengths.push_back(sizeof(Real) * 3);
               }
            
               // send  PERBX_DT2, PERBY_DT2, PERBZ_DT2
               if((SpatialCell::mpi_transfer_type & Transfer::CELL_PERBDT2)!=0){
                  displacements.push_back((uint8_t*) &(this->parameters[CellParams::PERBX_DT2]) - (uint8_t*) this);
                  block_lengths.push_back(sizeof(Real) * 3);
               }
            
               // send RHO, RHOVX, RHOVY, RHOVZ
               if((SpatialCell::mpi_transfer_type & Transfer::CELL_RHO_RHOV)!=0){
                  displacements.push_back((uint8_t*) &(this->parameters[CellParams::RHO]) - (uint8_t*) this);
                  block_lengths.push_back(sizeof(Real) * 4);
               }
            
               // send RHO_DT2, RHOVX_DT2, RHOVY_DT2, RHOVZ_DT2
               if((SpatialCell::mpi_transfer_type & Transfer::CELL_RHODT2_RHOVDT2)!=0){
                  displacements.push_back((uint8_t*) &(this->parameters[CellParams::RHO_DT2]) - (uint8_t*) this);
                  block_lengths.push_back(sizeof(Real) * 4);
               }
            
               // send  spatial cell derivatives
               if((SpatialCell::mpi_transfer_type & Transfer::CELL_DERIVATIVES)!=0){
                  displacements.push_back((uint8_t*) &(this->derivatives[0]) - (uint8_t*) this);
                  block_lengths.push_back(sizeof(Real) * fieldsolver::N_SPATIAL_CELL_DERIVATIVES);
                  
               }

               // send  spatial cell BVOL derivatives
               if((SpatialCell::mpi_transfer_type & Transfer::CELL_BVOL_DERIVATIVES)!=0){
                  displacements.push_back((uint8_t*) &(this->derivativesBVOL[0]) - (uint8_t*) this);
                  block_lengths.push_back(sizeof(Real) * bvolderivatives::N_BVOL_DERIVATIVES);
               }

               if((SpatialCell::mpi_transfer_type & Transfer::CELL_IOLOCALCELLID)!=0){
                  displacements.push_back((uint8_t*) &(this->ioLocalCellId) - (uint8_t*) this);
                  block_lengths.push_back(sizeof(uint64_t));
               }

               // send  sysBoundaryFlag        
               if((SpatialCell::mpi_transfer_type & Transfer::CELL_SYSBOUNDARYFLAG)!=0){
                  displacements.push_back((uint8_t*) &(this->sysBoundaryFlag) - (uint8_t*) this);
                  block_lengths.push_back(sizeof(uint));
                  displacements.push_back((uint8_t*) &(this->sysBoundaryLayer) - (uint8_t*) this);
                  block_lengths.push_back(sizeof(uint));
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
            }
            
            if (displacements.size() > 0) {
               count = 1;
               MPI_Type_create_hindexed(
                  displacements.size(),
                  &block_lengths[0],
                  &displacements[0],
                  MPI_BYTE,
                  &datatype
                  );
            } else {
               count = 0;
               datatype = MPI_BYTE;
            }
      }

      /*!
        Returns true if given velocity block has enough of a distribution function.
        Returns false if the value of the distribution function is too low in every
        sense in given block.
        Also returns false if given block doesn't exist or is an error block.
      */
      bool compute_block_has_content(const unsigned int block) const {
         if (block == error_velocity_block
             || this->velocity_blocks.count(block) == 0) {
            return false;
         }
         
         bool has_content = false;
         const Velocity_Block* block_ptr = &(this->velocity_blocks.at(block));
         
         for (unsigned int i = 0; i < VELOCITY_BLOCK_LENGTH; i++) {
            if (block_ptr->data[i] >= SpatialCell::velocity_block_min_value) {
               has_content = true;
               break;
            }
         }
         
         return has_content;
      }

      /*!
        Update the two lists containing blocks with content, and blocks without content.
      */
      void update_velocity_block_content_lists(void){
         this->velocity_block_with_content_list.clear();
         this->velocity_block_with_no_content_list.clear();
         
         for(unsigned int block_index=0;block_index<this->number_of_blocks;block_index++){
            unsigned int block = this->velocity_block_list[block_index];
            if(this->compute_block_has_content(block)){
               this->velocity_block_with_content_list.push_back(block);
            }
            else{
               this->velocity_block_with_no_content_list.push_back(block);
            }
         }
         
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

         const velocity_block_indices_t neighbor_indices = {{
            indices[0] + x_offset,
            indices[1] + y_offset,
            indices[2] + z_offset
         }};

         if (neighbor_indices[0] == error_velocity_block_index) {
            return error_velocity_block;
         }

         return get_velocity_block(neighbor_indices);
      }


      /*!  Adds "important" and removes "unimportant" velocity blocks
        to/from this cell.
        
        velocity_block_with_content_list needs to be up to date in local and remote cells.
        velocity_block_with_no_content_list needs to be up to date in local cells.
        
        update_velocity_block_with_content_lists() should have
        been called with the current distribution function values, and then the contetn list transferred.
        
        Removes all velocity blocks from this spatial cell which don't
        have content and don't have spatial or velocity neighbors with
        content.  Adds neighbors for all velocity blocks which do have
        content (including spatial neighbors).  All cells in
        spatial_neighbors are assumed to be neighbors of this cell.
        
        
        This function is thread-safe when called for different cells
        per thread. We need the block_has_content vector from
        neighbouring cells, but these are not written to here. We only
        modify local cell.
      */
      void adjust_velocity_blocks(const std::vector<SpatialCell*>& spatial_neighbors, bool doDeleteEmptyBlocks=true)
      {
         //  This set contains all those cellids which have neighbors in any
         //  of the 6-dimensions Actually, we would only need to add
         //  local blocks with no content here, as blocks with content
         //  do not need to be created and also will not be removed as
         //  we only check for removal for blocks with no content
         boost::unordered_set<unsigned int> neighbors_have_content;
         
         //add neighbor content info for velocity space neighbors to map. We loop over blocks
         //with content and raise the neighbors_have_content for
         //itself, and for all its neighbors
         for(unsigned int block_index=0;block_index< this->velocity_block_with_content_list.size();block_index++){
            unsigned int block = this->velocity_block_with_content_list[block_index];
            neighbors_have_content.insert(block); //also add the cell itself
            for(int offset_vx=-P::sparseBlockAddWidthV;offset_vx<=P::sparseBlockAddWidthV;offset_vx++)
               for(int offset_vy=-P::sparseBlockAddWidthV;offset_vy<=P::sparseBlockAddWidthV;offset_vy++)
		  for(int offset_vz=-P::sparseBlockAddWidthV;offset_vz<=P::sparseBlockAddWidthV;offset_vz++){                  
                     const unsigned int neighbor_block = get_velocity_block_from_offsets(block, offset_vx,offset_vy,offset_vz);
                     if (neighbor_block == error_velocity_block) {
                        continue;
                     }
                     neighbors_have_content.insert(neighbor_block);
                  }
         }
         //add neighbor content info for spatial space neighbors to map. We loop over
         //neighbor cell lists with existing blocks, and raise the
         //flag for the local block with same block id
         for (std::vector<SpatialCell*>::const_iterator neighbor = spatial_neighbors.begin();
              neighbor != spatial_neighbors.end(); neighbor++ ) {
            for(unsigned int block_index=0;block_index< (*neighbor)->velocity_block_with_content_list.size();block_index++){
               unsigned int block = (*neighbor)->velocity_block_with_content_list[block_index];
               neighbors_have_content.insert(block);
            }
         }

         // ADD all blocks with neighbors in spatial or velocity space (if it exists then the block is unchanged)
         for (boost::unordered_set<unsigned int>::iterator it= neighbors_have_content.begin(); it != neighbors_have_content.end();++it) {
            this->add_velocity_block(*it);
         }
         
         
         // REMOVE all blocks in this cell without content + without neighbors with content
         if(doDeleteEmptyBlocks) {
            for(unsigned int block_index=0;block_index< this->velocity_block_with_no_content_list.size();block_index++){
               unsigned int block = this->velocity_block_with_no_content_list[block_index];
               boost::unordered_set<unsigned int>::iterator it = neighbors_have_content.find(block);
               if(it==neighbors_have_content.end()) {
                  //No content, and also no neighbor have content -> remove
                  //increment rho loss counters
                  Velocity_Block* block_ptr = this->at(block);
                  const Real DV3 = block_ptr->parameters[BlockParams::DVX]*
                     block_ptr->parameters[BlockParams::DVY]*
                     block_ptr->parameters[BlockParams::DVZ];
                  Real sum=0;
                  for(unsigned int i=0;i<WID3;i++)
                     sum+=block_ptr->data[i];
                  this->parameters[CellParams::RHOLOSSADJUST]+=DV3*sum;
                  //and finally remove block
                  this->remove_velocity_block(block);
               }
            }
         }
      }
      
      void adjustSingleCellVelocityBlocks() {
         //neighbor_ptrs is empty as we do not have any consistent
         //data in neighbours yet, adjustments done only based on velocity
         //space. Do not delete blocks, because otherwise we might loose teneous flow in spatial space
         std::vector<SpatialCell*> neighbor_ptrs;
         this->update_velocity_block_content_lists();
         this->adjust_velocity_blocks(neighbor_ptrs,false);
      }
      
      //set block data pointers data and fx for block
      //velocity_block_list[block_index], so that they point to the
      //same index as in block_list in the block_data and block_fx
      //which contain all data for all blocks. It just sets the
      //pointers and does not care about what is there earlier.
      //velocity_block_list needs to be up-to-date
      void set_block_data_pointers(int block_index){
         Velocity_Block* tmp_block_ptr = this->at(this->velocity_block_list[block_index]);
         tmp_block_ptr->data=&(this->block_data[block_index*VELOCITY_BLOCK_LENGTH]);
         tmp_block_ptr->fx=&(this->block_fx[block_index*SIZE_FLUXS]);

      }
      
      //This function will resize block data if needed, resize happen
      //in steps of block_allocation_chunk. It will always preserve
      //space for one extra block, resize can happen after adding the
      //block. We need up-to-date velocity_block_list as
      //set_block_data_pointers needs it.
      //If there is only free space left for 1 additional block (extra
      //block should not be needed, but there for safety), resize it
      //so that we have free space for block_allocation_chunk blocks.
      //If there is free space for more than 2*block_allocation_chunk
      //blocks, resize it so that we have free space for an additional
      //block_allocation_chunks blocks.
      void resize_block_data(){
         const int block_allocation_chunk=100;
         
         if((this->number_of_blocks+1)*VELOCITY_BLOCK_LENGTH>=this->block_data.size() ||
            (this->number_of_blocks+2*block_allocation_chunk)*VELOCITY_BLOCK_LENGTH<this->block_data.size()){
            
            //resize so that free space is block_allocation_chunk blocks
            int new_size=this->number_of_blocks+block_allocation_chunk;
            this->block_data.resize(new_size*VELOCITY_BLOCK_LENGTH);
            this->block_fx.resize(new_size*SIZE_FLUXS);
            
            //fix block pointers if a reallocation occured
            for(unsigned int block_index=0;block_index<this->number_of_blocks;block_index++){
               set_block_data_pointers(block_index);
            }
         }
      }
      
      
      /*!
        Adds an empty velocity block into this spatial cell.
        Returns true if given block was added or already exists.
        Returns false if given block is invalid or would be outside
        of the velocity grid.
      */
      bool add_velocity_block(const unsigned int block) {
         if (this->velocity_blocks.count(block) > 0) {
            return true;
         }

         if (block == error_velocity_block) {
            std::cerr << "ERROR: trying to add a block with the key error_velocity_block!" << std::endl;
            return false;
         }

         if (block >= SpatialCell::max_velocity_blocks) {
            std::cerr << "ERROR: trying to add a block with a key > max_velocity_blocks!" << std::endl;
            return false;
         }


         //create block
         this->velocity_blocks[block];
         //update number of blocks
         this->number_of_blocks++;
                  
         //add block to list of existing block
         this->velocity_block_list.push_back(block);
         //add more space for block data 
         resize_block_data();
         //fix block data pointers   
         set_block_data_pointers(this->number_of_blocks-1);
         
         //get pointer to block
         Velocity_Block* block_ptr = this->at(block);



         //clear block
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
         

         return true;
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
         
         //Find where in the block list the removed block was (index to block list). We need to fill this hole.    
         unsigned int block_index = -1;
         for(unsigned int i=0; i< velocity_block_list.size();i++) {
            if(this->velocity_block_list[i] == block) {
               block_index=i;
               break;
            }
         }
         
         /*
           Move the last existing block in the block list
           to the removed block's position
         */
         this->velocity_block_list[block_index] = this->velocity_block_list.back();
         this->velocity_block_list.pop_back(); //remove last item
            

         //copy velocity block data to the removed blocks position in order to fill the hole
         for(unsigned int i=0;i<VELOCITY_BLOCK_LENGTH;i++){
            this->block_data[block_index*VELOCITY_BLOCK_LENGTH+i] = this->block_data[(this->number_of_blocks - 1)*VELOCITY_BLOCK_LENGTH+i];
         }
         for(unsigned int i=0;i<SIZE_FLUXS;i++){
            this->block_fx[block_index*SIZE_FLUXS+i] = this->block_fx[(this->number_of_blocks - 1)*SIZE_FLUXS+i];
         }
         //set block data pointers to the location where we copied data
         set_block_data_pointers(block_index);

         //reduce number of blocks
         this->number_of_blocks--;

         //also remove velocity block structure
         this->velocity_blocks.erase(block);

         //check if we can decrease memory consumption
         resize_block_data();         
      }

      
      /*!
        Removes all velocity blocks from this spatial cell and frees memory in the cell
      */
      void clear(void)
         {
            //use the swap trick to force c++ to release the memory held by the vectors & maps
            //FIXME: we could jsut as well just set this->vector=std::vector<> ?
            boost::unordered_map<unsigned int, Velocity_Block> ().swap(this->velocity_blocks);
            //remove block data (value and fluxes)
            std::vector<Real,aligned_allocator<Real,64> >().swap(this->block_data);
            std::vector<Real,aligned_allocator<Real,64> >().swap(this->block_fx);
            std::vector<unsigned int>().swap(this->velocity_block_list);
            std::vector<unsigned int>().swap(this->mpi_velocity_block_list);
            this->number_of_blocks=0;
         }

      /*!       
        Prepares this spatial cell to receive the velocity grid over MPI.

        At this stage we have received a new blocklist over MPI into
        mpi_velocity_block_list, but the rest of the cell structures
        have not been adapted to this new list. Here we re-initialize
        the cell with empty blocks based on the new list
      */
      void prepare_to_receive_blocks(void) {
         // take copy of existing block  list so that we can safely loop over it while removing blocks
         std::vector<unsigned int> old_block_list(&(this->velocity_block_list[0]),
                                                  &(this->velocity_block_list[this->number_of_blocks]));
         
         //make set of  mpi block list to make element finding faster, use pointers as iterators
         std::set<unsigned int> mpi_velocity_block_set(&(this->mpi_velocity_block_list[0]),
                                                  &(this->mpi_velocity_block_list[this->mpi_number_of_blocks]));

         
         
         //remove blocks that are not in the new list
         for (unsigned int block_index = 0; block_index < old_block_list.size(); block_index++) {
            if(mpi_velocity_block_set.count(old_block_list[block_index])==0){
               this->remove_velocity_block(old_block_list[block_index]);
            }
         }
         
         // add velocity blocks that are about to be received with MPI, if it exists then
         //add_velocity_blocks safely returns
         for (unsigned int block_index = 0; block_index < this->mpi_number_of_blocks; block_index++) {
            this->add_velocity_block(this->mpi_velocity_block_list[block_index]);
         }

         //now we need to fix the order of blocks, the pointers in the velocity blocks have to point to the correct places
         //we do not copy the actual data, so that will be out of sync until we have received new data 
         for (unsigned int block_index = 0; block_index < this->mpi_number_of_blocks; block_index++) {
            //re-order velocity_block list to be the same that we received
            this->velocity_block_list[block_index]=this->mpi_velocity_block_list[block_index];

            //set array pointers to correct memory segment in velocity block
            set_block_data_pointers(block_index);

         }
      }
     /*! Apply velocity boundary condition.
       
       The outermost layer of velocity cells are set to zero, and the loss counters are updated with the mass lost in this operation.
     */
      void applyVelocityBoundaryCondition(){
            //is not     computed as the other cell does not exist = no outflow).
            //x-1 face
         for(unsigned int i=0;i<number_of_blocks;i++){
            Velocity_Block* block_ptr=this->at(velocity_block_list[i]);
            const Real DV3 = block_ptr->parameters[BlockParams::DVX]*
               block_ptr->parameters[BlockParams::DVY]*
               block_ptr->parameters[BlockParams::DVZ];
            
            if(is_null_block(block_ptr->neighbors[velocity_neighbor::XM1_YCC_ZCC])){
               unsigned int i=0;
               for(unsigned int k=0;k<WID;k++)
                  for(unsigned int j=0;j<WID;j++){
                     this->parameters[CellParams::RHOLOSSVELBOUNDARY]+=DV3*block_ptr->data[i+WID*j+WID2*k];;
                     block_ptr->data[i+WID*j+WID2*k]=0.0;
                  }
            }
            //x+1 face           
            if(is_null_block(block_ptr->neighbors[velocity_neighbor::XP1_YCC_ZCC])){
               unsigned int i=WID-1;
               for(unsigned int k=0;k<WID;k++)
                  for(unsigned int j=0;j<WID;j++){
                     this->parameters[CellParams::RHOLOSSVELBOUNDARY]+=DV3*block_ptr->data[i+WID*j+WID2*k];;
                     block_ptr->data[i+WID*j+WID2*k]=0.0;
                  }
            }
            
            //y-1  face
            if(is_null_block(block_ptr->neighbors[velocity_neighbor::XCC_YM1_ZCC])){
               unsigned int j=0;
               for(unsigned int k=0;k<WID;k++)
                  for(unsigned int i=0;i<WID;i++){
                     this->parameters[CellParams::RHOLOSSVELBOUNDARY]+=DV3*block_ptr->data[i+WID*j+WID2*k];;
                     block_ptr->data[i+WID*j+WID2*k]=0.0;
                  }
            }
            //y+1 face   
            if(is_null_block(block_ptr->neighbors[velocity_neighbor::XCC_YP1_ZCC])){
               unsigned int j=WID-1;
               for(unsigned int k=0;k<WID;k++)
                  for(unsigned int i=0;i<WID;i++){
                     this->parameters[CellParams::RHOLOSSVELBOUNDARY]+=DV3*block_ptr->data[i+WID*j+WID2*k];;
                     block_ptr->data[i+WID*j+WID2*k]=0.0;
                  }
            }
            //z-1 face          
            if(is_null_block(block_ptr->neighbors[velocity_neighbor::XCC_YCC_ZM1])){
               unsigned int k=0;
               for(unsigned int j=0;j<WID;j++)
                  for(unsigned int i=0;i<WID;i++){
                     this->parameters[CellParams::RHOLOSSVELBOUNDARY]+=DV3*block_ptr->data[i+WID*j+WID2*k];;
                     block_ptr->data[i+WID*j+WID2*k]=0.0;
                  }
            }
            //z+1 face           
            if(is_null_block(block_ptr->neighbors[velocity_neighbor::XCC_YCC_ZP1])){
               unsigned int k=WID-1;
               for(unsigned int j=0;j<WID;j++)
                  for(unsigned int i=0;i<WID;i++){
                     this->parameters[CellParams::RHOLOSSVELBOUNDARY]+=DV3*block_ptr->data[i+WID*j+WID2*k];;
                     block_ptr->data[i+WID*j+WID2*k]=0.0;
                  }
            }
         }
      }
      
      /*!
        Sets the type of data to transfer by mpi_datatype.
      */
      static void set_mpi_transfer_type(const uint64_t type,bool atSysBoundaries=false)
      {
         SpatialCell::mpi_transfer_type = type;
         SpatialCell::mpiTransferAtSysBoundaries=atSysBoundaries;
      }

      /*!
        Gets the type of data that will be transferred by mpi_datatype.
      */
      static uint64_t get_mpi_transfer_type(void)
      {
         return SpatialCell::mpi_transfer_type;
      }

     /*!
       Set if this cell is transferred/received using MPI in the next communication phase.
     */

      void set_mpi_transfer_enabled(bool transferEnabled){
        this->mpiTransferEnabled=transferEnabled;
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
      static uint64_t mpi_transfer_type;

      /*
        Do we only transfer data at boundaries (true), or in the whole system (false)
      */
      static bool mpiTransferAtSysBoundaries;
      /*  
        Minimum value of distribution function
        in any phase space cell of a velocity block for the
        block to be considered to have content
      */
      static Real velocity_block_min_value;


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


      //Storage container for velocity blocks that exist in this cell
      boost::unordered_map<unsigned int, Velocity_Block> velocity_blocks;

     
     bool mpiTransferEnabled;
      
   public:
      //number of blocks in cell (should be equal to velocity_block_list.size())
      unsigned int number_of_blocks;
      // List of velocity blocks in this cell,
      std::vector<unsigned int> velocity_block_list;
      //number of blocks in mpi_velocity_block_list
      unsigned int mpi_number_of_blocks;
      //this list is used for communicating a velocity block list over MPI
      std::vector<unsigned int>  mpi_velocity_block_list;

      //List of existing cells with content, only up-to-date after
      //call to update_has_content
      std::vector<unsigned int> velocity_block_with_content_list;
      //Size of vector. Needed for MPI communication of size before actual list transfer.
      unsigned int velocity_block_with_content_list_size;
      //List of existing cells with no content, only up-to-date after
      //call to update_has_content. This is also never transferred
      //over MPI, so is invalid on remote cells
      std::vector<unsigned int> velocity_block_with_no_content_list;
      

      //vectors for storing block data. We set pointers to this
      //datastorage in set_block_date_pointers()
      std::vector<Real,aligned_allocator<Real,64> > block_data;
      std::vector<Real,aligned_allocator<Real,64> > block_fx;

      //vectors for storing null block data
      std::vector<Real,aligned_allocator<Real,64> > null_block_data;
      std::vector<Real,aligned_allocator<Real,64> > null_block_fx;
      
      /*
        Bulk variables in this spatial cell.
      */
      Real parameters[CellParams::N_SPATIAL_CELL_PARAMS];

      /*
        Derivatives of bulk variables in this spatial cell.
      */
      Real derivatives[fieldsolver::N_SPATIAL_CELL_DERIVATIVES];
      
      // Derivatives of BVOL needed by the acceleration. Separate array because it does not need to be communicated.
      Real derivativesBVOL[bvolderivatives::N_BVOL_DERIVATIVES];
      
      //neighbor id's. Kept up to date in solvers, not by the spatial_cell class
      std::vector<uint64_t> neighbors;
     
      
      unsigned int procBoundaryFlag; /*!< bitfield usied in leveque vlasov solver to see if a neighbor exists, or if it is outside the system. TODO: bad/missleading name */
      uint sysBoundaryFlag;          /*!< What type of system boundary does the cell belong to. Enumerated in the sysboundarytype namespace's enum */
      uint sysBoundaryLayer;         /*!< Layers counted from closest systemBoundary. If 0 then it has not been computed. First sysboundary layer is layer 1 */
      uint64_t ioLocalCellId;       /*!< Local cell ID used for IO, not needed elsewhere and thus not being kept up-to-date*/ 
   }; // class SpatialCell
   

} // namespaces
#endif

