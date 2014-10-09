
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

#include "velocity_mesh_old.h"
#include "velocity_block_container.h"

typedef Parameters P;

// size of velocity blocks in velocity cells
#define block_vx_length WID
#define block_vy_length WID
#define block_vz_length WID
//this is also defined in common.h as SIZE_VELBLOCK, we should remove either one
#define VELOCITY_BLOCK_LENGTH WID3
#define N_NEIGHBOR_VELOCITY_BLOCKS 28

//extra memory allocated for block data. Should be in parameters and read in
//#define block_allocation_factor 1.1


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

namespace spatial_cell {

   namespace Transfer {
      const uint64_t NONE                     = 0;
      const uint64_t CELL_PARAMETERS          = (1<<0);
      const uint64_t CELL_DERIVATIVES         = (1<<1);
      const uint64_t VEL_BLOCK_LIST_STAGE1    = (1<<2);
      const uint64_t VEL_BLOCK_LIST_STAGE2    = (1<<3);
      const uint64_t VEL_BLOCK_DATA           = (1<<4);
      const uint64_t VEL_BLOCK_DATA_TO_FLUXES         = (1<<5);
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
      const uint64_t NEIGHBOR_VEL_BLOCK_FLUXES = (1<<21);
      const uint64_t CELL_HALL_TERM           = (1<<22);
      const uint64_t CELL_P                   = (1<<23);
      const uint64_t CELL_PDT2                = (1<<24);
      
      //all data, expect for the fx table (never needed on remote cells)
      const uint64_t ALL_DATA =
      CELL_PARAMETERS
      | CELL_DERIVATIVES | CELL_BVOL_DERIVATIVES
      | VEL_BLOCK_DATA
      | CELL_SYSBOUNDARYFLAG;
      //all data, except the distribution function
      const uint64_t ALL_SPATIAL_DATA =
      CELL_PARAMETERS
      | CELL_DERIVATIVES | CELL_BVOL_DERIVATIVES
      | CELL_SYSBOUNDARYFLAG;
   }

   #warning TODO: typedef unsigned int velocity_cell_t;
   #warning TODO: typedef unsigned int velocity_block_t;
   
   typedef boost::array<unsigned int, 3> velocity_cell_indices_t;             /**< Defines the indices of a velocity cell in a velocity block.
									       * Indices start from 0 and the first value is the index in x direction.
									       * Note: these are the (i,j,k) indices of the cell within the block.
									       * Valid values are ([0,block_vx_length[,[0,block_vy_length[,[0,block_vz_length[).*/
   
   typedef boost::array<unsigned int, 3> velocity_block_indices_t;            /**< Defines the indices of a velocity block in the velocity grid.
									       * Indices start from 0 and the first value is the index in x direction.
									       * Note: these are the (i,j,k) indices of the block.
									       * Valid values are ([0,vx_length[,[0,vy_length[,[0,vz_length[).*/

   class SpatialCell {
   public:
      SpatialCell();
      SpatialCell(const SpatialCell& other);

      // Following functions return velocity grid metadata //
      Realf* get_data();
      const Realf* get_data() const;
      Realf* get_data(const vmesh::LocalID& blockLID);
      const Realf* get_data(const vmesh::LocalID& blockLID) const;
      Realf* get_fx();
      Realf* get_fx(const vmesh::LocalID& blockLID);
      Real* get_block_parameters();
      Real* get_block_parameters(const vmesh::LocalID& blockLID);
      const Real* get_block_parameters(const vmesh::LocalID& blockLID) const;
      vmesh::LocalID get_number_of_velocity_blocks() const;
      static const unsigned int* get_velocity_base_grid_length();
      static const Real* get_velocity_base_grid_block_size();
      static const Real* get_velocity_base_grid_cell_size();
      static velocity_block_indices_t get_velocity_block_indices(const vmesh::GlobalID globalID);                             // OK
      static vmesh::GlobalID get_velocity_block(const velocity_block_indices_t indices);                                      // OK
      static vmesh::GlobalID get_velocity_block(const Real vx,const Real vy,const Real vz);
      vmesh::GlobalID get_velocity_block_global_id(const vmesh::LocalID& blockLID) const;
      vmesh::LocalID get_velocity_block_local_id(const vmesh::GlobalID& blockGID) const;
      static void get_velocity_block_size(const vmesh::GlobalID block,Real size[3]);
      static Real get_velocity_block_vx_min(const vmesh::GlobalID block);
      static Real get_velocity_block_vx_max(const vmesh::GlobalID block);
      static Real get_velocity_block_vy_min(const vmesh::GlobalID block);
      static Real get_velocity_block_vy_max(const vmesh::GlobalID block);
      static Real get_velocity_block_vz_min(const vmesh::GlobalID block);
      static Real get_velocity_block_vz_max(const vmesh::GlobalID block);
      static velocity_cell_indices_t get_velocity_cell_indices(const unsigned int cell);
      static unsigned int get_velocity_cell(const velocity_cell_indices_t indices);
      static unsigned int get_velocity_cell(const vmesh::GlobalID velocity_block,const Real vx,const Real vy,const Real vz);
      static Real get_velocity_cell_vx_min(const vmesh::GlobalID velocity_block,const unsigned int velocity_cell);
      static Real get_velocity_cell_vx_max(const vmesh::GlobalID velocity_block,const unsigned int velocity_cell);
      static Real get_velocity_cell_vy_min(const vmesh::GlobalID velocity_block,const unsigned int velocity_cell);
      static Real get_velocity_cell_vy_max(const vmesh::GlobalID velocity_block,const unsigned int velocity_cell);
      static Real get_velocity_cell_vz_min(const vmesh::GlobalID velocity_block,const unsigned int velocity_cell);
      static Real get_velocity_cell_vz_max(const vmesh::GlobalID velocity_block,const unsigned int velocity_cell);
      static const Real* get_velocity_grid_min_limits();
      static void initialize_mesh(Real v_limits[6],unsigned int meshSize[3],unsigned int blockSize[3],Real f_min);
      static unsigned int invalid_block_index();
      static vmesh::GlobalID invalid_global_id();
      static vmesh::LocalID invalid_local_id();

      size_t count(const vmesh::GlobalID& block) const;

      // Following functions adjust velocity blocks stored on the cell //
      bool add_velocity_block(const vmesh::GlobalID& block);
      void adjustSingleCellVelocityBlocks();
      void adjust_velocity_blocks(const std::vector<SpatialCell*>& spatial_neighbors, bool doDeleteEmptyBlocks=true);
      void update_velocity_block_content_lists(void);
      void clear(void);
      uint64_t get_cell_memory_capacity();
      uint64_t get_cell_memory_size();
      void prepare_to_receive_blocks(void);
      bool shrink_to_fit();
      size_t size(void) const;
      void remove_velocity_block(const vmesh::GlobalID& block);

      Real get_value(const Real vx, const Real vy, const Real vz) const;
      void increment_value(const Real vx, const Real vy, const Real vz, const Realf value);
      void increment_value(const vmesh::GlobalID& block,const unsigned int cell, const Real value);
      void set_value(const Real vx, const Real vy, const Real vz, const Realf value);
      void set_value(const vmesh::GlobalID& block,const unsigned int cell, const Realf value);



      // Following functions are related to MPI //
      boost::tuple<void*, int, MPI_Datatype> get_mpi_datatype(const CellID cellID,const int sender_rank,const int receiver_rank,
							      const bool receiving,const int neighborhood);
      static uint64_t get_mpi_transfer_type(void);
      static void set_mpi_transfer_type(const uint64_t type,bool atSysBoundaries=false);
      void set_mpi_transfer_enabled(bool transferEnabled);
      
      // Member variables //
      Real derivatives[fieldsolver::N_SPATIAL_CELL_DERIVATIVES];              /**< Derivatives of bulk variables in this spatial cell.*/
      Real derivativesBVOL[bvolderivatives::N_BVOL_DERIVATIVES];              /**< Derivatives of BVOL needed by the acceleration. 
									       * Separate array because it does not need to be communicated.*/
      Real parameters[CellParams::N_SPATIAL_CELL_PARAMS];                     /**< Bulk variables in this spatial cell.*/
      Realf null_block_data[WID3];
      Realf null_block_fx[WID3];

      uint64_t ioLocalCellId;                                                 /**< Local cell ID used for IO, not needed elsewhere 
									       * and thus not being kept up-to-date.*/
      vmesh::LocalID mpi_number_of_blocks;                                    /**< Number of blocks in mpi_velocity_block_list.*/
      std::vector<vmesh::GlobalID>  mpi_velocity_block_list;                  /**< This list is used for communicating a velocity block list over MPI.*/
      Realf* neighbor_block_data;                                             /**< Pointers for translation operator. We can point to neighbor
									       * cell block data. We do not allocate memory for the pointer.*/
      vmesh::LocalID neighbor_number_of_blocks;
      uint sysBoundaryFlag;                                                   /**< What type of system boundary does the cell belong to. 
									       * Enumerated in the sysboundarytype namespace's enum.*/
      uint sysBoundaryLayer;                                                  /**< Layers counted from closest systemBoundary. If 0 then it has not 
									       * been computed. First sysboundary layer is layer 1.*/
      std::vector<vmesh::GlobalID> velocity_block_with_content_list;          /**< List of existing cells with content, only up-to-date after
									       * call to update_has_content().*/
      vmesh::LocalID velocity_block_with_content_list_size;                   /**< Size of vector. Needed for MPI communication of size before actual list transfer.*/
      std::vector<vmesh::GlobalID> velocity_block_with_no_content_list;       /**< List of existing cells with no content, only up-to-date after
									       * call to update_has_content. This is also never transferred
									       * over MPI, so is invalid on remote cells.*/
      static uint64_t mpi_transfer_type;                                      /**< Which data is transferred by the mpi datatype given by spatial cells.*/
      static bool mpiTransferAtSysBoundaries;                                 /**< Do we only transfer data at boundaries (true), or in the whole system (false).*/
      static Real velocity_block_min_value;                                   /**< Minimum value of distribution function in any phase space cell 
									       * of a velocity block for the block to be considered to have content.*/

    private:
      SpatialCell& operator=(const SpatialCell&);
      
      bool compute_block_has_content(const vmesh::GlobalID& block) const;
      vmesh::GlobalID get_velocity_block_from_offsets(const vmesh::GlobalID& block,const int x_offset,
						      const int y_offset,const int z_offset);
      void resize_block_data();
      void set_block_data_pointers(int block_index);

      bool initialized;
      bool mpiTransferEnabled;
      vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID> vmesh;
      vmesh::VelocityBlockContainer<vmesh::LocalID> blockContainer;
   };

   /****************************
    * Velocity block functions *
    ****************************/

   inline Realf* SpatialCell::get_data() {
      return blockContainer.getData();
   }
   
   inline const Realf* SpatialCell::get_data() const {
      return blockContainer.getData();
   }

   inline Realf* SpatialCell::get_data(const vmesh::LocalID& blockLID) {
      return blockContainer.getData(blockLID);
   }
   
   inline const Realf* SpatialCell::get_data(const vmesh::LocalID& blockLID) const {
      return blockContainer.getData(blockLID);
   }

   inline Realf* SpatialCell::get_fx() {
      return blockContainer.getFx();
   }

   inline Realf* SpatialCell::get_fx(const vmesh::LocalID& blockLID) {
      return blockContainer.getFx(blockLID);
   }
   
   inline Real* SpatialCell::get_block_parameters() {
      return blockContainer.getParameters();
   }
   
   inline Real* SpatialCell::get_block_parameters(const vmesh::LocalID& blockLID) {
      return blockContainer.getParameters(blockLID);
   }
   
   inline const Real* SpatialCell::get_block_parameters(const vmesh::LocalID& blockLID) const {
      return blockContainer.getParameters(blockLID);
   }

   inline vmesh::LocalID SpatialCell::get_number_of_velocity_blocks() const {
      return blockContainer.size();
   }

   inline const unsigned int* SpatialCell::get_velocity_base_grid_length() {
      return vmesh::VelocityMesh<uint32_t,uint32_t>::getBaseGridLength();
   }

   inline const Real* SpatialCell::get_velocity_base_grid_block_size() {
      return vmesh::VelocityMesh<uint32_t,uint32_t>::getBaseGridBlockSize();
   }

   inline const Real* SpatialCell::get_velocity_base_grid_cell_size() {
      return vmesh::VelocityMesh<uint32_t,uint32_t>::getBaseGridCellSize();
   }

   /*!
    Returns the indices of given velocity block
    */
   inline velocity_block_indices_t SpatialCell::get_velocity_block_indices(const vmesh::GlobalID block) {
      velocity_block_indices_t indices;
      uint32_t refLevel;
      vmesh::VelocityMesh<uint32_t,uint32_t>::getIndices(block,refLevel,indices[0],indices[1],indices[2]);
      return indices;
   }

   /*!
    Returns the velocity block at given indices or error_velocity_block
    */
   inline vmesh::GlobalID SpatialCell::get_velocity_block(const velocity_block_indices_t indices) {
      const uint32_t refLevel = 0;
      return vmesh::VelocityMesh<uint32_t,uint32_t>::getGlobalID(refLevel,indices[0],indices[1],indices[2]);
   }

   /*!
    Returns the velocity block at given location or
    error_velocity_block if outside of the velocity grid
    */
   inline vmesh::GlobalID SpatialCell::get_velocity_block(
							  const Real vx,
							  const Real vy,
							  const Real vz
							  ) {
      return vmesh::VelocityMesh<uint32_t,uint32_t>::getGlobalID(vx,vy,vz);
   }
   
   inline vmesh::GlobalID SpatialCell::get_velocity_block_global_id(const vmesh::LocalID& blockLID) const {
      return vmesh.getGlobalID(blockLID);
   }
   
   inline vmesh::LocalID SpatialCell::get_velocity_block_local_id(const vmesh::GlobalID& blockGID) const {
      return vmesh.getLocalID(blockGID);
   }

   inline void SpatialCell::get_velocity_block_size(const vmesh::GlobalID block,Real blockSize[3]) {
      vmesh::VelocityMesh<uint32_t,uint32_t>::getBlockSize(block,blockSize);
   }
   
   /*!
    Returns the edge where given velocity block starts.
    */
   inline Real SpatialCell::get_velocity_block_vx_min(const vmesh::GlobalID block) {
      Real coords[3];
      vmesh::VelocityMesh<uint32_t,uint32_t>::getBlockCoordinates(block,coords);
      return coords[0];
   }

   /*!
    Returns the edge where given velocity block ends.
    */
   inline Real SpatialCell::get_velocity_block_vx_max(const vmesh::GlobalID block) {
      Real coords[3];
      vmesh::VelocityMesh<uint32_t,uint32_t>::getBlockCoordinates(block,coords);
      
      Real size[3];
      vmesh::VelocityMesh<uint32_t,uint32_t>::getBlockSize(block,size);
      return coords[0]+size[0];
   }

   /*!
    Returns the edge where given velocity block starts.
    */
   inline Real SpatialCell::get_velocity_block_vy_min(const vmesh::GlobalID block) {
      Real coords[3];
      vmesh::VelocityMesh<uint32_t,uint32_t>::getBlockCoordinates(block,coords);
      return coords[1];
   }

   /*!
    Returns the edge where given velocity block ends.
    */
   inline Real SpatialCell::get_velocity_block_vy_max(const vmesh::GlobalID block) {
      Real coords[3];
      vmesh::VelocityMesh<uint32_t,uint32_t>::getBlockCoordinates(block,coords);
      
      Real size[3];
      vmesh::VelocityMesh<uint32_t,uint32_t>::getBlockSize(block,size);
      return coords[1]+size[1];
   }

   /*!
    Returns the edge where given velocity block starts.
    */
   inline Real SpatialCell::get_velocity_block_vz_min(const vmesh::GlobalID block) {
      Real coords[3];
      vmesh::VelocityMesh<uint32_t,uint32_t>::getBlockCoordinates(block,coords);
      return coords[2];
   }
   
   /*!
    Returns the edge where given velocity block ends.
    */
   inline Real SpatialCell::get_velocity_block_vz_max(const vmesh::GlobalID block) {
      Real coords[3];
      vmesh::VelocityMesh<uint32_t,uint32_t>::getBlockCoordinates(block,coords);
      
      Real size[3];
      vmesh::VelocityMesh<uint32_t,uint32_t>::getBlockSize(block,size);
      return coords[2]+size[2];
   }

      /***************************
       * Velocity cell functions *
       ***************************/

   /*!
    Returns the indices of given velocity cell
    */
   inline velocity_cell_indices_t SpatialCell::get_velocity_cell_indices(const unsigned int cell) {
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
   inline unsigned int SpatialCell::get_velocity_cell(const velocity_cell_indices_t indices) {
      if (indices[0] >= block_vx_length
	  || indices[1] >= block_vy_length
	  || indices[2] >= block_vz_length) {
	 return error_velocity_cell;
      }
      return indices[0] + indices[1] * block_vx_length + indices[2] * block_vx_length * block_vy_length;
   }

   /*!     
    Returns the velocity cell at given location or
    error_velocity_cell if outside of given velocity block.
    */
   inline unsigned int SpatialCell::get_velocity_cell(
						      const vmesh::GlobalID velocity_block,
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
   inline Real SpatialCell::get_velocity_cell_vx_min(
						     const vmesh::GlobalID velocity_block,
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
   inline Real SpatialCell::get_velocity_cell_vx_max(
						     const vmesh::GlobalID velocity_block,
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
   inline Real SpatialCell::get_velocity_cell_vy_min(
						     const vmesh::GlobalID velocity_block,
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
   inline Real SpatialCell::get_velocity_cell_vy_max(
						     const vmesh::GlobalID velocity_block,
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
   inline Real SpatialCell::get_velocity_cell_vz_min(
						     const vmesh::GlobalID velocity_block,
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
   inline Real SpatialCell::get_velocity_cell_vz_max(
						     const vmesh::GlobalID velocity_block,
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
   
   inline const Real* SpatialCell::get_velocity_grid_min_limits() {
      return vmesh::VelocityMesh<uint32_t,uint32_t>::getMeshMinLimits();
   }

   inline void SpatialCell::initialize_mesh(Real v_limits[6],unsigned int meshSize[3],unsigned int blockSize[3],Real f_min) {
      vmesh::VelocityMesh<uint32_t,uint32_t>::initialize(v_limits,meshSize,blockSize);
      velocity_block_min_value = f_min;
   }

   inline unsigned int SpatialCell::invalid_block_index() {
      return vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>::invalidBlockIndex();
   }

   inline vmesh::GlobalID SpatialCell::invalid_global_id() {
      return vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>::invalidGlobalID();
   }

   inline vmesh::GlobalID SpatialCell::invalid_local_id() {
      return vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>::invalidLocalID();
   }

   inline SpatialCell::SpatialCell() {
      /*
       Block list and cache always have room for all blocks
       */
      this->mpi_number_of_blocks=0;
      this->sysBoundaryLayer=0; /*!< Default value, layer not yet initialized*/
      for (unsigned int i=0; i<WID3; ++i) null_block_data[i] = 0.0;
      for (unsigned int i=0; i<WID3; ++i) null_block_fx[i] = 0.0;
      
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
   
   inline SpatialCell::SpatialCell(const SpatialCell& other):
     initialized(other.initialized),
     mpiTransferEnabled(other.mpiTransferEnabled),
     mpi_number_of_blocks(other.mpi_number_of_blocks),
     mpi_velocity_block_list(other.mpi_velocity_block_list),
     velocity_block_with_content_list(other.velocity_block_with_content_list),
     velocity_block_with_no_content_list(other.velocity_block_with_no_content_list),
     sysBoundaryFlag(other.sysBoundaryFlag),
     sysBoundaryLayer(other.sysBoundaryLayer),
     vmesh(other.vmesh), blockContainer(other.blockContainer) {
	//       phiprof::initializeTimer("SpatialCell copy", "SpatialCell copy");
	//       phiprof::start("SpatialCell copy");
	// 
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
	  
	//set null block data
	for (unsigned int i=0; i<WID3; ++i) null_block_data[i] = 0.0;
	for (unsigned int i=0; i<WID3; ++i) null_block_fx[i] = 0.0;
	//         phiprof::stop("SpatialCell copy");
     }

   /*!
    Returns the number of given velocity blocks that exist.
    */
   inline size_t SpatialCell::count(const vmesh::GlobalID& block) const {
      return vmesh.count(block);
   }

   /*!
    Returns the number of existing velocity blocks.
    */
   inline size_t SpatialCell::size(void) const {
      return vmesh.size();
   }
   
   /*!
    Sets the given value to a velocity cell at given coordinates.
    * 
    Creates the velocity block at given coordinates if it doesn't exist.
    */
   inline void SpatialCell::set_value(const Real vx, const Real vy, const Real vz, const Realf value) {
      const vmesh::GlobalID blockGID = get_velocity_block(vx, vy, vz);
      if (count(blockGID) == 0) {
	 if (!this->add_velocity_block(blockGID)) {
	    std::cerr << "Couldn't add velocity block " << blockGID << std::endl;
	    abort();
	 }
      }
      
      const vmesh::LocalID blockLID = get_velocity_block_local_id(blockGID);
      if (blockLID == invalid_local_id()) {
	 std::cerr << __FILE__ << ":" << __LINE__
	   << " block_ptr == NULL" << std::endl;
	 abort();
      }

      const unsigned int cell = get_velocity_cell(blockGID, vx, vy, vz);
      get_data(blockLID)[cell] = value;
   }

//TODO - thread safe set/increment functions which do not create blocks automatically

   /*! Sets the value of a particular cell in a block. The block is
    *  created if it does not exist. This version is faster than
    *  the velocity value based version.
    *
    * This function is not thread safe due to the creation of
    * blocks.
    * 
    \param block Block index of velocity-cell
    \param cell  Cell index (0..WID3-1) of velocity-cell in block
    \param value Value that is set for velocity-cell
    */
   inline void SpatialCell::set_value(const vmesh::GlobalID& blockGID,const unsigned int cell, const Realf value) {
      if (count(blockGID) == 0) {
	 if (!this->add_velocity_block(blockGID)) {
	    std::cerr << "Couldn't add velocity block " << blockGID << std::endl;
	    abort();
	 }
      }
      const vmesh::LocalID blockLID = get_velocity_block_local_id(blockGID);
      if (blockLID == invalid_local_id()) {
	 std::cerr << __FILE__ << ":" << __LINE__
	   << " block_ptr == NULL" << std::endl;
	 abort();
      }
      get_data(blockLID)[cell] = value;
   }
          
   /*!
    Increments the value of velocity cell at given coordinate-
    * 
    Creates the velocity block at given coordinates if it doesn't exist.
    */
   inline void SpatialCell::increment_value(const Real vx, const Real vy, const Real vz, const Realf value) {
      const vmesh::GlobalID blockGID = get_velocity_block(vx, vy, vz);
      if (count(blockGID) == 0) {
	 if (!this->add_velocity_block(blockGID)) {
	    std::cerr << "Couldn't add velocity block " << blockGID << std::endl;
	    abort();
	 }
      }
      const vmesh::LocalID blockLID = get_velocity_block_local_id(blockGID);
      if (blockLID == invalid_local_id()) {
	 std::cerr << __FILE__ << ":" << __LINE__
	   << " block_ptr == NULL" << std::endl;
	 abort();
      }

      const unsigned int cell = get_velocity_cell(blockGID, vx, vy, vz);
      get_data(blockLID)[cell] += value;
   }
   
   /*!
    Increments the value of velocity cell at given index
    * 
    Creates the velocity block if it doesn't exist.
    */
   inline void SpatialCell::increment_value(const vmesh::GlobalID& blockGID,const unsigned int cell, const Real value) {
      if (count(blockGID) == 0) {
	 if (!this->add_velocity_block(blockGID)) {
	    std::cerr << "Couldn't add velocity block " << blockGID << std::endl;
	    abort();
	 }
      }
      const vmesh::LocalID blockLID = get_velocity_block_local_id(blockGID);
      if (blockLID == invalid_local_id()) {
	 std::cerr << __FILE__ << ":" << __LINE__
	   << " block_ptr == NULL" << std::endl;
	 abort();
      }
      get_data(blockLID)[cell] += value;
   }

   /*!
    * Gets the value of a velocity cell at given coordinates.
    * 
    * Returns 0 if it doesn't exist.
    */
   inline Real SpatialCell::get_value(const Real vx, const Real vy, const Real vz) const {
      const vmesh::GlobalID blockGID = get_velocity_block(vx, vy, vz);
      if (count(blockGID) == 0) {
	 return 0.0;
      }
      const vmesh::LocalID blockLID = get_velocity_block_local_id(blockGID);
      //const Velocity_Block block_ptr = at(block);
      if (blockLID == invalid_local_id()) {
      //if (block_ptr.is_null() == true) {
	 std::cerr << __FILE__ << ":" << __LINE__
	   << " block_ptr == NULL" << std::endl; 
	 abort();
      }

      const unsigned int cell = get_velocity_cell(blockGID, vx, vy, vz);
      // Cast to real: Note block_ptr->data[cell] is Realf type
      const Real value = get_data(blockLID)[cell];
      return value;
   }
   
   /*! get mpi datatype for sending the cell data. */
   inline boost::tuple<void*, int, MPI_Datatype> SpatialCell::get_mpi_datatype(
									       const CellID cellID/*cell_id*/,
									       const int sender_rank/*sender*/,
									       const int receiver_rank/*receiver*/,
									       const bool receiving,
									       const int neighborhood) {
      std::vector<MPI_Aint> displacements;
      std::vector<int> block_lengths;
      vmesh::LocalID block_index = 0;
      
      /*create datatype for actual data if we are in the first two layers around a boundary, or if we send for thw whole system*/
      if (this->mpiTransferEnabled && (SpatialCell::mpiTransferAtSysBoundaries==false || this->sysBoundaryLayer ==1 || this->sysBoundaryLayer ==2 )) {
	 //add data to send/recv to displacement and block length lists
	 if ((SpatialCell::mpi_transfer_type & Transfer::VEL_BLOCK_LIST_STAGE1) != 0) {
	    //first copy values in case this is the send operation
	    this->mpi_number_of_blocks = blockContainer.size();

	    // send velocity block list size
	    displacements.push_back((uint8_t*) &(this->mpi_number_of_blocks) - (uint8_t*) this);
	    block_lengths.push_back(sizeof(vmesh::LocalID));
	 }

	 if ((SpatialCell::mpi_transfer_type & Transfer::VEL_BLOCK_LIST_STAGE2) != 0) {
	    // STAGE1 should have been done, otherwise we have problems...
	    if (receiving) {
	       //mpi_number_of_blocks transferred earlier
	       this->mpi_velocity_block_list.resize(this->mpi_number_of_blocks);
	    } else {
	       //resize to correct size (it will avoid reallociation if it is big enough, I assume)
	       this->mpi_number_of_blocks = blockContainer.size();
	       this->mpi_velocity_block_list.resize(blockContainer.size());
	       
	       //copy values if this is the send operation
	       for (vmesh::LocalID i=0; i<blockContainer.size(); ++i) {
		  this->mpi_velocity_block_list[i] = vmesh.getGlobalID(i);
	       }
	    }

	    // send velocity block list
	    displacements.push_back((uint8_t*) &(this->mpi_velocity_block_list[0]) - (uint8_t*) this);
	    block_lengths.push_back(sizeof(vmesh::GlobalID) * this->mpi_number_of_blocks);
	 }

	 if ((SpatialCell::mpi_transfer_type & Transfer::VEL_BLOCK_WITH_CONTENT_STAGE1) !=0) {
	    //Communicate size of list so that buffers can be allocated on receiving side
	    if (!receiving) this->velocity_block_with_content_list_size = this->velocity_block_with_content_list.size();
	    displacements.push_back((uint8_t*) &(this->velocity_block_with_content_list_size) - (uint8_t*) this);
	    block_lengths.push_back(sizeof(vmesh::LocalID));
	 }
	 if ((SpatialCell::mpi_transfer_type & Transfer::VEL_BLOCK_WITH_CONTENT_STAGE2) !=0) {
	    if (receiving) {
	       this->velocity_block_with_content_list.resize(this->velocity_block_with_content_list_size);
	    }
	    
	    //velocity_block_with_content_list_size should first be updated, before this can be done (STAGE1)
	    displacements.push_back((uint8_t*) &(this->velocity_block_with_content_list[0]) - (uint8_t*) this);
	    block_lengths.push_back(sizeof(vmesh::GlobalID)*this->velocity_block_with_content_list_size);
	 }
	 
	 if ((SpatialCell::mpi_transfer_type & Transfer::VEL_BLOCK_DATA) !=0) {
	    displacements.push_back((uint8_t*) get_data() - (uint8_t*) this);
	    block_lengths.push_back(sizeof(Realf) * VELOCITY_BLOCK_LENGTH * blockContainer.size());
	 }

	 if ((SpatialCell::mpi_transfer_type & Transfer::VEL_BLOCK_DATA_TO_FLUXES) != 0) {
	    if (receiving) {
	       displacements.push_back((uint8_t*) get_fx() - (uint8_t*) this);
	    } else {
	       displacements.push_back((uint8_t*) get_data() - (uint8_t*) this);
	    }
	    block_lengths.push_back(sizeof(Realf) * VELOCITY_BLOCK_LENGTH * blockContainer.size());
	 }

	 if ((SpatialCell::mpi_transfer_type & Transfer::NEIGHBOR_VEL_BLOCK_FLUXES) != 0) {
	    /*We are actually transfering the data of a
	     * neighbor. The values of neighbor_block_data
	     * and neighbor_number_of_blocks should be set in
	     * solver.*/               
	    displacements.push_back((uint8_t*) this->neighbor_block_data - (uint8_t*) this);               
	    block_lengths.push_back(sizeof(Realf) * VELOCITY_BLOCK_LENGTH* this->neighbor_number_of_blocks);
	 }

	 // send  spatial cell parameters
	 if ((SpatialCell::mpi_transfer_type & Transfer::CELL_PARAMETERS)!=0){
	    displacements.push_back((uint8_t*) &(this->parameters[0]) - (uint8_t*) this);
	    block_lengths.push_back(sizeof(Real) * CellParams::N_SPATIAL_CELL_PARAMS);
	 }

	 // send  spatial cell dimensions
	 if ((SpatialCell::mpi_transfer_type & Transfer::CELL_DIMENSIONS)!=0){
	    displacements.push_back((uint8_t*) &(this->parameters[CellParams::DX]) - (uint8_t*) this);
	    block_lengths.push_back(sizeof(Real) * 3);
	 }
	 
	 // send  BGBX BGBY BGBZ and all edge-averaged BGBs
	 if ((SpatialCell::mpi_transfer_type & Transfer::CELL_BGB)!=0){
	    displacements.push_back((uint8_t*) &(this->parameters[CellParams::BGBX]) - (uint8_t*) this);
	    block_lengths.push_back(sizeof(Real) * 3);
	    displacements.push_back((uint8_t*) &(this->parameters[CellParams::BGBX_000_010]) - (uint8_t*) this);
	    block_lengths.push_back(sizeof(Real) * 24);
	 }
            
	 // send  BGBXVOL BGBYVOL BGBZVOL PERBXVOL PERBYVOL PERBZVOL
	 if ((SpatialCell::mpi_transfer_type & Transfer::CELL_BVOL)!=0){
	    displacements.push_back((uint8_t*) &(this->parameters[CellParams::BGBXVOL]) - (uint8_t*) this);
	    block_lengths.push_back(sizeof(Real) * 6);
	 }
	 
	 // send  EX, EY EZ
	 if ((SpatialCell::mpi_transfer_type & Transfer::CELL_E)!=0){
	    displacements.push_back((uint8_t*) &(this->parameters[CellParams::EX]) - (uint8_t*) this);
	    block_lengths.push_back(sizeof(Real) * 3);
	 }
            
	 // send  EX_DT2, EY_DT2, EZ_DT2
	 if ((SpatialCell::mpi_transfer_type & Transfer::CELL_EDT2)!=0){
	    displacements.push_back((uint8_t*) &(this->parameters[CellParams::EX_DT2]) - (uint8_t*) this);
	    block_lengths.push_back(sizeof(Real) * 3);
	 }
	 
	 // send  PERBX, PERBY, PERBZ
	 if ((SpatialCell::mpi_transfer_type & Transfer::CELL_PERB)!=0){
	    displacements.push_back((uint8_t*) &(this->parameters[CellParams::PERBX]) - (uint8_t*) this);
	    block_lengths.push_back(sizeof(Real) * 3);
	 }
            
	 // send  PERBX_DT2, PERBY_DT2, PERBZ_DT2
	 if ((SpatialCell::mpi_transfer_type & Transfer::CELL_PERBDT2)!=0){
	    displacements.push_back((uint8_t*) &(this->parameters[CellParams::PERBX_DT2]) - (uint8_t*) this);
	    block_lengths.push_back(sizeof(Real) * 3);
	 }
	 
	 // send RHO, RHOVX, RHOVY, RHOVZ
	 if ((SpatialCell::mpi_transfer_type & Transfer::CELL_RHO_RHOV)!=0){
	    displacements.push_back((uint8_t*) &(this->parameters[CellParams::RHO]) - (uint8_t*) this);
	    block_lengths.push_back(sizeof(Real) * 4);
	 }
	 
	 // send RHO_DT2, RHOVX_DT2, RHOVY_DT2, RHOVZ_DT2
	 if ((SpatialCell::mpi_transfer_type & Transfer::CELL_RHODT2_RHOVDT2)!=0){
	    displacements.push_back((uint8_t*) &(this->parameters[CellParams::RHO_DT2]) - (uint8_t*) this);
	    block_lengths.push_back(sizeof(Real) * 4);
	 }
	 
	 // send  spatial cell derivatives
	 if ((SpatialCell::mpi_transfer_type & Transfer::CELL_DERIVATIVES)!=0){
	    displacements.push_back((uint8_t*) &(this->derivatives[0]) - (uint8_t*) this);
	    block_lengths.push_back(sizeof(Real) * fieldsolver::N_SPATIAL_CELL_DERIVATIVES);
	 }
	 
	 // send  spatial cell BVOL derivatives
	 if ((SpatialCell::mpi_transfer_type & Transfer::CELL_BVOL_DERIVATIVES)!=0){
	    displacements.push_back((uint8_t*) &(this->derivativesBVOL[0]) - (uint8_t*) this);
	    block_lengths.push_back(sizeof(Real) * bvolderivatives::N_BVOL_DERIVATIVES);
	 }
            
	 if ((SpatialCell::mpi_transfer_type & Transfer::CELL_IOLOCALCELLID)!=0){
	    displacements.push_back((uint8_t*) &(this->ioLocalCellId) - (uint8_t*) this);
	    block_lengths.push_back(sizeof(uint64_t));
	 }
	 
	 // send Hall term components
	 if ((SpatialCell::mpi_transfer_type & Transfer::CELL_HALL_TERM)!=0){
	    displacements.push_back((uint8_t*) &(this->parameters[CellParams::EXHALL_000_100]) - (uint8_t*) this);
	    block_lengths.push_back(sizeof(Real) * 12);
	 }
	 
	 // send P tensor diagonal components
	 if ((SpatialCell::mpi_transfer_type & Transfer::CELL_P)!=0){
	    displacements.push_back((uint8_t*) &(this->parameters[CellParams::P_11]) - (uint8_t*) this);
	    block_lengths.push_back(sizeof(Real) * 3);
	 }
	 
	 if ((SpatialCell::mpi_transfer_type & Transfer::CELL_PDT2)!=0){
	    displacements.push_back((uint8_t*) &(this->parameters[CellParams::P_11_DT2]) - (uint8_t*) this);
	    block_lengths.push_back(sizeof(Real) * 3);
	 }
	 
	 // send  sysBoundaryFlag
	 if ((SpatialCell::mpi_transfer_type & Transfer::CELL_SYSBOUNDARYFLAG)!=0){
	    displacements.push_back((uint8_t*) &(this->sysBoundaryFlag) - (uint8_t*) this);
	    block_lengths.push_back(sizeof(uint));
	    displacements.push_back((uint8_t*) &(this->sysBoundaryLayer) - (uint8_t*) this);
	    block_lengths.push_back(sizeof(uint));
	 }
            
	 if ((SpatialCell::mpi_transfer_type & Transfer::VEL_BLOCK_PARAMETERS) !=0) {
	    displacements.push_back((uint8_t*) get_block_parameters() - (uint8_t*) this);
	    block_lengths.push_back(sizeof(Real) * size() * BlockParams::N_VELOCITY_BLOCK_PARAMS);
	 }
      }

      void* address = this;
      int count;
      MPI_Datatype datatype;      
      
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

      return boost::make_tuple(address,count,datatype);
   }

   /*!
    Returns true if given velocity block has enough of a distribution function.
    Returns false if the value of the distribution function is too low in every
    sense in given block.
    Also returns false if given block doesn't exist or is an error block.
    */
   inline bool SpatialCell::compute_block_has_content(const vmesh::GlobalID& blockGID) const {
      if (blockGID == invalid_global_id() || count(blockGID) == 0) return false;
      
      bool has_content = false;
      const vmesh::LocalID blockLID = get_velocity_block_local_id(blockGID);
      const Realf* block_data = blockContainer.getData(blockLID);
      for (unsigned int i=0; i<VELOCITY_BLOCK_LENGTH; ++i) {
	 if (block_data[i] >= SpatialCell::velocity_block_min_value) {
	    has_content = true;
	    break;
	 }
      }

      return has_content;
   }
   
   /*!
    Update the two lists containing blocks with content, and blocks without content.
    */
   inline void SpatialCell::update_velocity_block_content_lists(void) {
      this->velocity_block_with_content_list.clear();
      this->velocity_block_with_no_content_list.clear();

      for (vmesh::LocalID block_index=0; block_index<vmesh.size(); ++block_index) {
	 const vmesh::GlobalID globalID = vmesh.getGlobalID(block_index);
	 if (this->compute_block_has_content(globalID)){
	    this->velocity_block_with_content_list.push_back(globalID);
	 } else {
	    this->velocity_block_with_no_content_list.push_back(globalID);
	 }
      }
   }

   /*!
    Removes all velocity blocks from this spatial cell and frees memory in the cell
    */
   inline void SpatialCell::clear(void) {
      vmesh.clear();
      blockContainer.clear();

      // use the swap trick to force c++ to release the memory held by the vectors & maps
      std::vector<vmesh::GlobalID>().swap(this->mpi_velocity_block_list);
   }

   /*!  Purges extra capacity from block vectors. It sets size to
    num_blocks * block_allocation_factor  (if capacity greater than this), and also forces capacity to this new smaller value.
    \return True on success
    */
   inline bool SpatialCell::shrink_to_fit() {
      const size_t amount = 2 + blockContainer.size() * blockContainer.getBlockAllocationFactor();

      // Allow capacity to be a bit large than needed by number of blocks, shrink otherwise
      if (blockContainer.capacity() > amount * VELOCITY_BLOCK_LENGTH) return blockContainer.recapacitate(amount);
      return true;
   }

         /*!
     This function will resize (increase) block data if needed, resize happen
     in steps of block_allocation_chunk. It will always preserve
     space for one extra block, resize can happen after adding the
     block. We need up-to-date velocity_block_list as
     set_block_data_pointers needs it.
     If there is only free space left for 1 additional block (extra
     block should not be needed, but there for safety), resize it
     so that we have free space for block_allocation_chunk blocks.

     This function never decreases memory space. To do that one should
     call shrink_to_fit(!)
     
   */
   inline void SpatialCell::resize_block_data() {
      std::cerr << "resize_block_data(): error not implemented" << std::endl;
      exit(1);
      /*
      if ((this->number_of_blocks+1)*VELOCITY_BLOCK_LENGTH >= this->block_data.size()){
	 // Resize so that free space is block_allocation_chunk blocks, and at least two in case of having zero blocks
	 int new_size = 2 + this->number_of_blocks * block_allocation_factor;
	 this->block_data.resize(new_size*VELOCITY_BLOCK_LENGTH);
	 this->block_fx.resize(new_size*VELOCITY_BLOCK_LENGTH);	 
	 blocks.resize(new_size);

	 // Fix block pointers if a reallocation occured
	 for (unsigned int block_index=0; block_index<this->number_of_blocks; ++block_index){
	    set_block_data_pointers(block_index);
	 }
      }*/
   }

   /*!
    Return the memory consumption in bytes as reported using the size()
    functions of the containers in spatial cell
    */
   inline uint64_t SpatialCell::get_cell_memory_size() {
      const uint64_t VEL_BLOCK_SIZE = 2*WID3*sizeof(Realf) + BlockParams::N_VELOCITY_BLOCK_PARAMS*sizeof(Real);
      uint64_t size = 0;
      size += vmesh.sizeInBytes();
      size += blockContainer.sizeInBytes();
      size += 2 * WID3 * sizeof(Realf);
      size += mpi_velocity_block_list.size() * sizeof(vmesh::GlobalID);
      size += velocity_block_with_content_list.size() * sizeof(vmesh::GlobalID);
      size += velocity_block_with_no_content_list.size() * sizeof(vmesh::GlobalID);
      size += CellParams::N_SPATIAL_CELL_PARAMS * sizeof(Real);
      size += fieldsolver::N_SPATIAL_CELL_DERIVATIVES * sizeof(Real);
      size += bvolderivatives::N_BVOL_DERIVATIVES * sizeof(Real);
      return size;
   }
   
   /*!
    Return the memory consumption in bytes as reported using
    the size() functions of the containers in spatial cell
    */
   inline uint64_t SpatialCell::get_cell_memory_capacity() {
      const uint64_t VEL_BLOCK_SIZE = 2*WID3*sizeof(Realf) + BlockParams::N_VELOCITY_BLOCK_PARAMS*sizeof(Real);
      uint64_t capacity = 0;
      capacity += vmesh.capacityInBytes();
      capacity += blockContainer.capacityInBytes();
      capacity += 2 * WID3 * sizeof(Realf);
      capacity += mpi_velocity_block_list.capacity()  * sizeof(vmesh::GlobalID);
      capacity += velocity_block_with_content_list.capacity()  * sizeof(vmesh::GlobalID);
      capacity += velocity_block_with_no_content_list.capacity()  * sizeof(vmesh::GlobalID);
      capacity += CellParams::N_SPATIAL_CELL_PARAMS * sizeof(Real);
      capacity += fieldsolver::N_SPATIAL_CELL_DERIVATIVES * sizeof(Real);
      capacity += bvolderivatives::N_BVOL_DERIVATIVES * sizeof(Real);
      return capacity;
   }

   /*!
    Returns a block at given offsets from given block.
    * 
    Returns error_velocity_block if the returned block
    would be outside of the velocity grid or all given
    offsets are 0.
    */
   inline vmesh::GlobalID SpatialCell::get_velocity_block_from_offsets(
								       const vmesh::GlobalID& block,
								       const int x_offset, 
								       const int y_offset, 
								       const int z_offset
								       ) {
      #warning This function needs to be redefined for AMR mesh
      if (x_offset == 0  &&  y_offset == 0 && z_offset == 0) {
	 return invalid_global_id();
      }
      
      const velocity_block_indices_t indices = get_velocity_block_indices(block);
      if (indices[0] == invalid_block_index()) {
	 return invalid_global_id();
      }
      
      const velocity_block_indices_t neighbor_indices = {{
	 indices[0] + x_offset,
	   indices[1] + y_offset,
	   indices[2] + z_offset
      }};
      
      if (neighbor_indices[0] == invalid_block_index()) {
	 return invalid_global_id();
      }

      return get_velocity_block(neighbor_indices);
   }

   /*!  Adds "important" and removes "unimportant" velocity blocks
    to/from this cell.
    * 
    velocity_block_with_content_list needs to be up to date in local and remote cells.
    velocity_block_with_no_content_list needs to be up to date in local cells.
        
    update_velocity_block_with_content_lists() should have
    been called with the current distribution function values, and then the contetn list transferred.
    * 
    Removes all velocity blocks from this spatial cell which don't
    have content and don't have spatial or velocity neighbors with
    content.  Adds neighbors for all velocity blocks which do have
    content (including spatial neighbors).  All cells in
    spatial_neighbors are assumed to be neighbors of this cell.
    * 
        
    This function is thread-safe when called for different cells
    per thread. We need the block_has_content vector from
    neighbouring cells, but these are not written to here. We only
    modify local cell.
    */
   inline void SpatialCell::adjust_velocity_blocks(const std::vector<SpatialCell*>& spatial_neighbors, bool doDeleteEmptyBlocks) {
      //  This set contains all those cellids which have neighbors in any
      //  of the 6-dimensions Actually, we would only need to add
      //  local blocks with no content here, as blocks with content
      //  do not need to be created and also will not be removed as
      //  we only check for removal for blocks with no content
      boost::unordered_set<vmesh::GlobalID> neighbors_have_content;
      
      //add neighbor content info for velocity space neighbors to map. We loop over blocks
      //with content and raise the neighbors_have_content for
      //itself, and for all its neighbors
      for (vmesh::LocalID block_index=0; block_index<velocity_block_with_content_list.size(); ++block_index) {
	 vmesh::GlobalID block = velocity_block_with_content_list[block_index];
	 const velocity_block_indices_t indices = get_velocity_block_indices(block);
	 neighbors_have_content.insert(block); //also add the cell itself
	 for (int offset_vx=-P::sparseBlockAddWidthV;offset_vx<=P::sparseBlockAddWidthV;offset_vx++)
	   for (int offset_vy=-P::sparseBlockAddWidthV;offset_vy<=P::sparseBlockAddWidthV;offset_vy++)
	     for (int offset_vz=-P::sparseBlockAddWidthV;offset_vz<=P::sparseBlockAddWidthV;offset_vz++){                  
		const vmesh::GlobalID neighbor_block = get_velocity_block({{indices[0] + offset_vx, indices[1] + offset_vy, indices[2] + offset_vz}});
		neighbors_have_content.insert(neighbor_block); //add all potential ngbrs of this block with content
	     }
      }
      //add neighbor content info for spatial space neighbors to map. We loop over
      //neighbor cell lists with existing blocks, and raise the
      //flag for the local block with same block id
      for (std::vector<SpatialCell*>::const_iterator neighbor = spatial_neighbors.begin();
	   neighbor != spatial_neighbors.end(); neighbor++ ) {
	 for (vmesh::LocalID block_index=0;block_index< (*neighbor)->velocity_block_with_content_list.size();block_index++){
	    vmesh::GlobalID block = (*neighbor)->velocity_block_with_content_list[block_index];
	    neighbors_have_content.insert(block);
	 }
      }
      
      // REMOVE all blocks in this cell without content + without neighbors with content
      /*better to do it in the reverse order, as then blocks at the
       * end are removed first, and we may avoid copying extra
       * data.*/
      if (doDeleteEmptyBlocks) {
	 for (int block_index= this->velocity_block_with_no_content_list.size()-1; block_index>=0; --block_index) {
	    const vmesh::GlobalID blockGID = this->velocity_block_with_no_content_list[block_index];
	    const vmesh::LocalID blockLID = get_velocity_block_local_id(blockGID);
	    boost::unordered_set<vmesh::GlobalID>::iterator it = neighbors_have_content.find(blockGID);
	    if (it == neighbors_have_content.end()) {
	       //No content, and also no neighbor have content -> remove
	       //increment rho loss counters
	       const Real* block_parameters = get_block_parameters(blockLID);
	       const Real DV3 = block_parameters[BlockParams::DVX]
		              * block_parameters[BlockParams::DVY]
		              * block_parameters[BlockParams::DVZ];
	       Real sum=0;
	       for (unsigned int i=0; i<WID3; ++i) sum += get_data(blockLID)[i];
	       this->parameters[CellParams::RHOLOSSADJUST] += DV3*sum;

	       // and finally remove block
	       this->remove_velocity_block(blockGID);
	    }
	 }
      }
      
      // ADD all blocks with neighbors in spatial or velocity space (if it exists then the block is unchanged)
      for (boost::unordered_set<vmesh::GlobalID>::iterator it=neighbors_have_content.begin(); it != neighbors_have_content.end(); ++it) {
	 this->add_velocity_block(*it);
      }
   }
   
   inline void SpatialCell::adjustSingleCellVelocityBlocks() {
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
   //velocity_block_list needs to be up-to-date.
   //We also set block_list_index here to the correct position
   inline void SpatialCell::set_block_data_pointers(int block_index) { }
      
   /*!
    Adds an empty velocity block into this spatial cell.
    Returns true if given block was added or already exists.
    Returns false if given block is invalid or would be outside
    of the velocity grid.
    */
   inline bool SpatialCell::add_velocity_block(const vmesh::GlobalID& block) {
      // Block insert will fail, if the block already exists, or if 
      // there are too many blocks in the spatial cell
      bool success = true;
      if (vmesh.push_back(block) == false) {
	 return false;
      }

      const vmesh::LocalID VBC_LID = blockContainer.push_back();

      // Set block parameters:
      Real* parameters = get_block_parameters(vmesh.getLocalID(block));
      parameters[BlockParams::VXCRD] = get_velocity_block_vx_min(block);
      parameters[BlockParams::VYCRD] = get_velocity_block_vy_min(block);
      parameters[BlockParams::VZCRD] = get_velocity_block_vz_min(block);
      vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>::getCellSize(block,&(parameters[BlockParams::DVX]));

      // The following call 'should' be the fastest, but is actually 
      // much slower that the parameter setting above
      //vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>::getBlockInfo(block,get_block_parameters( blockContainer.push_back() ));
      return success;
   }

   /*!
    Removes given block from the velocity grid.
    Does nothing if given block doesn't exist.
    */
   inline void SpatialCell::remove_velocity_block(const vmesh::GlobalID& block) {
      if (block == invalid_global_id()) {
	 return;
      }
      if (count(block) == 0) return;

      // Get local ID of the removed block, and local ID of the last block:
      const vmesh::LocalID removedLID = vmesh.getLocalID(block);
      const vmesh::LocalID lastLID = vmesh.size()-1;
      vmesh.copy(lastLID,removedLID);
      vmesh.pop();

      blockContainer.copy(lastLID,removedLID);
      blockContainer.pop();
   }

   /*!       
    Prepares this spatial cell to receive the velocity grid over MPI.
    * 
    At this stage we have received a new blocklist over MPI into
    mpi_velocity_block_list, but the rest of the cell structures
    have not been adapted to this new list. Here we re-initialize
    the cell with empty blocks based on the new list
    */
   inline void SpatialCell::prepare_to_receive_blocks(void) {
      vmesh.setGrid(mpi_velocity_block_list);
      blockContainer.setSize(mpi_velocity_block_list.size());

      // Set velocity block parameters:
      for (vmesh::LocalID blockLID=0; blockLID<size(); ++blockLID) {
	 const vmesh::GlobalID blockGID = get_velocity_block_global_id(blockLID);
	 get_block_parameters(blockLID)[BlockParams::VXCRD] = get_velocity_block_vx_min(blockGID);
	 get_block_parameters(blockLID)[BlockParams::VYCRD] = get_velocity_block_vy_min(blockGID);
	 get_block_parameters(blockLID)[BlockParams::VZCRD] = get_velocity_block_vz_min(blockGID);
	 vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>::getCellSize(blockGID,&(get_block_parameters(blockLID)[BlockParams::DVX]));
      }
   }

   /*!
    Sets the type of data to transfer by mpi_datatype.
    */
   inline void SpatialCell::set_mpi_transfer_type(const uint64_t type,bool atSysBoundaries) {
      SpatialCell::mpi_transfer_type = type;
      SpatialCell::mpiTransferAtSysBoundaries = atSysBoundaries;
   }

   /*!
    Gets the type of data that will be transferred by mpi_datatype.
    */
   inline uint64_t SpatialCell::get_mpi_transfer_type(void) {
      return SpatialCell::mpi_transfer_type;
   }
   
   /*!
    Set if this cell is transferred/received using MPI in the next communication phase.
    */
   inline void SpatialCell::set_mpi_transfer_enabled(bool transferEnabled) {
      this->mpiTransferEnabled=transferEnabled;
   }

   inline SpatialCell& SpatialCell::operator=(const SpatialCell&) { 
      return *this;
   }

} // namespaces

#endif
