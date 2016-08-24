/*!
Spatial cell class for Vlasiator that supports a variable number of velocity blocks.

Copyright 2011-2015 Finnish Meteorological Institute
*/

#ifndef VLASIATOR_SPATIAL_CELL_HPP
#define VLASIATOR_SPATIAL_CELL_HPP

#include <algorithm>
#include <boost/array.hpp>
#include <boost/unordered_map.hpp>
#include <boost/lexical_cast.hpp>
#include <cmath>
#include <fstream>
#include <iostream>
#include <mpi.h>
#include <limits>
#include <stdint.h>
#include <vector>
#include <set>
#include <phiprof.hpp>
#include <tuple>

#include "memoryallocation.h"
#include "common.h"
#include "parameters.h"
#include "definitions.h"

#ifndef AMR
   #include "velocity_mesh_old.h"
#else
   #include "velocity_mesh_amr.h"
#endif

#include "amr_refinement_criteria.h"
#include "velocity_blocks.h"
#include "velocity_block_container.h"

#include "logger.h"
extern Logger logFile;

#ifndef NDEBUG
   #define DEBUG_SPATIAL_CELL
#endif

typedef Parameters P;

// size of velocity blocks in velocity cells
#define block_vx_length WID
#define block_vy_length WID
#define block_vz_length WID
//this is also defined in common.h as SIZE_VELBLOCK, we should remove either one
#define VELOCITY_BLOCK_LENGTH WID3
//#define N_NEIGHBOR_VELOCITY_BLOCKS 28

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
      const uint64_t NEIGHBOR_VEL_BLOCK_DATA  = (1<<21);
      const uint64_t CELL_HALL_TERM           = (1<<22);
      const uint64_t CELL_P                   = (1<<23);
      const uint64_t CELL_PDT2                = (1<<24);
      const uint64_t CELL_RHOQ_TOT            = (1<<25);
      const uint64_t CELL_PHI                 = (1<<26);
      const uint64_t POP_METADATA             = (1<<27);
      const uint64_t RANDOMGEN                = (1<<28);
      const uint64_t CELL_GRADPE_TERM         = (1<<29);
      //all data
      const uint64_t ALL_DATA =
      CELL_PARAMETERS
      | CELL_DERIVATIVES | CELL_BVOL_DERIVATIVES
      | VEL_BLOCK_DATA
      | CELL_SYSBOUNDARYFLAG
      | POP_METADATA | RANDOMGEN;

      //all data, except the distribution function
      const uint64_t ALL_SPATIAL_DATA =
      CELL_PARAMETERS
      | CELL_DERIVATIVES | CELL_BVOL_DERIVATIVES
      | CELL_SYSBOUNDARYFLAG
      | POP_METADATA | RANDOMGEN;
   }

   typedef boost::array<unsigned int, 3> velocity_cell_indices_t;             /**< Defines the indices of a velocity cell in a velocity block.
                                                                               * Indices start from 0 and the first value is the index in x direction.
                                                                               * Note: these are the (i,j,k) indices of the cell within the block.
                                                                               * Valid values are ([0,block_vx_length[,[0,block_vy_length[,[0,block_vz_length[).*/

   typedef boost::array<vmesh::LocalID,3> velocity_block_indices_t;           /**< Defines the indices of a velocity block in the velocity grid.
                                                                               * Indices start from 0 and the first value is the index in x direction.
                                                                               * Note: these are the (i,j,k) indices of the block.
                                                                               * Valid values are ([0,vx_length[,[0,vy_length[,[0,vz_length[).*/

   /** Wrapper for variables needed for each particle species.*/
   struct Population {
      Real max_dt[2];                                                /**< Element[0] is max_r_dt, element[1] max_v_dt.*/
      vmesh::LocalID N_blocks;                                       /**< Number of velocity blocks, used when receiving velocity 
                                                                      * mesh from remote neighbors using MPI.*/
      Real velocityBlockMinValue;
      vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID> vmesh;     /**< Velocity mesh. Contains all velocity blocks that exist 
                                                                      * in this spatial cell. Cells are identified by their unique 
                                                                      * global IDs.*/
      vmesh::VelocityBlockContainer<vmesh::LocalID> blockContainer;  /**< Velocity block data.*/
   };

   class SpatialCell {
   public:
      SpatialCell();
      SpatialCell(const SpatialCell& other);

      // Following functions return velocity grid metadata //
      template<int PAD> void fetch_data(const vmesh::GlobalID& blockGID,const vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& vmesh,
                                        const Realf* src,Realf* array);
      template<int PAD>	void fetch_acc_data(const vmesh::GlobalID& blockGID,const int& dim,
					    vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& vmesh,
					    const Realf* src,Realf* array,Real cellSizeFractions[2]);

      vmesh::GlobalID find_velocity_block(uint8_t& refLevel,vmesh::GlobalID cellIndices[3],const int& popID);
      Realf* get_data(const int& popID);
      const Realf* get_data(const int& popID) const;
      Realf* get_data(const vmesh::LocalID& blockLID,const int& popID);
      const Realf* get_data(const vmesh::LocalID& blockLID,const int& popID) const;
      Real* get_block_parameters(const int& popID);
      const Real* get_block_parameters(const int& popID) const;
      Real* get_block_parameters(const vmesh::LocalID& blockLID,const int& popID);
      const Real* get_block_parameters(const vmesh::LocalID& blockLID,const int& popID) const;

      Real* get_cell_parameters();
      const Real* get_cell_parameters() const;

      vmesh::LocalID get_number_of_velocity_blocks(const int& popID) const;
      vmesh::LocalID get_number_of_all_velocity_blocks() const;
      int get_number_of_populations() const;

      uint8_t get_maximum_refinement_level(const int& popID);
      const Real& get_max_r_dt(const int& popID) const;
      const Real& get_max_v_dt(const int& popID) const;

      const vmesh::LocalID* get_velocity_grid_length(const int& popID,const uint8_t& refLevel=0);
      const Real* get_velocity_grid_block_size(const int& popID,const uint8_t& refLevel=0);
      const Real* get_velocity_grid_cell_size(const int& popID,const uint8_t& refLevel=0);
      void get_velocity_block_coordinates(const int& popID,const vmesh::GlobalID& globalID,Real* coords);
      velocity_block_indices_t get_velocity_block_indices(const int& popID,const vmesh::GlobalID globalID);                             // OK
      velocity_block_indices_t get_velocity_block_indices(const int& popID,const vmesh::GlobalID globalID,uint8_t& refLevel);
      vmesh::GlobalID get_velocity_block(const int& popID,vmesh::GlobalID blockIndices[3],const uint8_t& refLevel) const;
      vmesh::GlobalID get_velocity_block(const int& popID,const velocity_block_indices_t indices,const uint8_t& refLevel) const;
      vmesh::GlobalID get_velocity_block(const int& popID,const Real* coords,const uint8_t& refLevel=0) const;
      vmesh::GlobalID get_velocity_block(const int& popID,const Real vx,const Real vy,const Real vz,const uint8_t& refLevel=0) const;
      vmesh::GlobalID get_velocity_block_child(const int& popID,const vmesh::GlobalID& blockGID,const uint8_t& refLevel,
                                               const int& i_cell,const int& j_cell,const int& k_cell);
      void get_velocity_block_children_local_ids(const vmesh::GlobalID& blockGID,
                                                 std::vector<vmesh::LocalID>& childrenLIDs,
                                                 const int& popID);
      vmesh::GlobalID get_velocity_block_parent(const int& popID,const vmesh::GlobalID& blockGID);
      vmesh::GlobalID get_velocity_block_global_id(const vmesh::LocalID& blockLID,const int& popID) const;
      vmesh::LocalID get_velocity_block_local_id(const vmesh::GlobalID& blockGID,const int& popID) const;
      void get_velocity_block_size(const int& popID,const vmesh::GlobalID block,Real size[3]);
      Real get_velocity_block_vx_min(const int& popID,const vmesh::GlobalID block) const;
      Real get_velocity_block_vx_max(const int& popID,const vmesh::GlobalID block) const;
      Real get_velocity_block_vy_min(const int& popID,const vmesh::GlobalID block) const;
      Real get_velocity_block_vy_max(const int& popID,const vmesh::GlobalID block) const;
      Real get_velocity_block_vz_min(const int& popID,const vmesh::GlobalID block) const;
      Real get_velocity_block_vz_max(const int& popID,const vmesh::GlobalID block) const;
      velocity_cell_indices_t get_velocity_cell_indices(const unsigned int cell) const;
      unsigned int get_velocity_cell(const velocity_cell_indices_t indices) const;
      unsigned int get_velocity_cell(const int& popID,const vmesh::GlobalID velocity_block,const Real vx,const Real vy,const Real vz) const;
      Real get_velocity_cell_vx_min(const int& popID,const vmesh::GlobalID velocity_block,const unsigned int velocity_cell) const;
      Real get_velocity_cell_vx_max(const int& popID,const vmesh::GlobalID velocity_block,const unsigned int velocity_cell) const;
      Real get_velocity_cell_vy_min(const int& popID,const vmesh::GlobalID velocity_block,const unsigned int velocity_cell) const;
      Real get_velocity_cell_vy_max(const int& popID,const vmesh::GlobalID velocity_block,const unsigned int velocity_cell) const;
      Real get_velocity_cell_vz_min(const int& popID,const vmesh::GlobalID velocity_block,const unsigned int velocity_cell) const;
      Real get_velocity_cell_vz_max(const int& popID,const vmesh::GlobalID velocity_block,const unsigned int velocity_cell) const;
      const Real* get_velocity_grid_min_limits(const int& popID);
      const Real* get_velocity_grid_max_limits(const int& popID);
      bool initialize_mesh();

      static unsigned int invalid_block_index();
      static vmesh::GlobalID invalid_global_id();
      static vmesh::LocalID invalid_local_id();

      size_t count(const vmesh::GlobalID& block,const int& popID) const;

      void add_values(const vmesh::GlobalID& targetGID,
		      std::unordered_map<vmesh::GlobalID,Realf[(WID+2)*(WID+2)*(WID+2)]>& sourceData,
                      const int& popID);

      void printMeshSizes();
      static bool setCommunicatedSpecies(const int& popID);

      // Following functions adjust velocity blocks stored on the cell //
      bool add_velocity_block(const vmesh::GlobalID& block,const int& popID);
      void add_velocity_blocks(const std::vector<vmesh::GlobalID>& blocks,const int& popID);
      bool add_velocity_block_octant(const vmesh::GlobalID& blockGID,const int& popID);
      void adjustSingleCellVelocityBlocks(const int& popID);
      void adjust_velocity_blocks(const std::vector<SpatialCell*>& spatial_neighbors,
                                  const int& popID,
                                  bool doDeleteEmptyBlocks=true);
      void update_velocity_block_content_lists(const int& popID);
      bool checkMesh(const int& popID);
      void clear(const int& popID);
      void coarsen_block(const vmesh::GlobalID& parent,const std::vector<vmesh::GlobalID>& children,const int& popID);
      void coarsen_blocks(amr_ref_criteria::Base* evaluator,const int& popID);
      uint64_t get_cell_memory_capacity();
      uint64_t get_cell_memory_size();
      void merge_values(const int& popID);
      void prepare_to_receive_blocks(const int& popID);
      bool shrink_to_fit();
      size_t size(const int& popID) const;
      void remove_velocity_block(const vmesh::GlobalID& block,const int& popID);
      void swap(vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& vmesh,
                vmesh::VelocityBlockContainer<vmesh::LocalID>& blockContainer,const int& popID);
      vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& get_velocity_mesh(const size_t& popID);
      vmesh::VelocityBlockContainer<vmesh::LocalID>& get_velocity_blocks(const size_t& popID);
      vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& get_velocity_mesh_temporary();
      vmesh::VelocityBlockContainer<vmesh::LocalID>& get_velocity_blocks_temporary();

      Realf get_value(const Real vx,const Real vy,const Real vz,const int& popID) const;
      Realf get_value(const vmesh::GlobalID& blockGID, const unsigned int cell, const int& popID) const;
      void increment_value(const Real vx,const Real vy,const Real vz,const Realf value,const int& popID);
      void increment_value(const vmesh::GlobalID& block,const unsigned int cell,const Realf value,const int& popID);
      void set_max_r_dt(const int& popID,const Real& value);
      void set_max_v_dt(const int& popID,const Real& value);
      void set_value(const Real vx, const Real vy, const Real vz, const Realf value,const int& popID);
      void set_value(const vmesh::GlobalID& block,const unsigned int cell, const Realf value,const int& popID);
      void refine_block(const vmesh::GlobalID& block,std::map<vmesh::GlobalID,vmesh::LocalID>& insertedBlocks,
                        const int& popID);
      bool velocity_block_has_children(const vmesh::GlobalID& blockGID,const int& popID) const;
      vmesh::GlobalID velocity_block_has_grandparent(const vmesh::GlobalID& blockGID,const int& popID) const;

      // Following functions are related to MPI //
      std::tuple<void*, int, MPI_Datatype> get_mpi_datatype(const CellID cellID,const int sender_rank,const int receiver_rank,
                                                            const bool receiving,const int neighborhood);
      static uint64_t get_mpi_transfer_type(void);
      static void set_mpi_transfer_type(const uint64_t type,bool atSysBoundaries=false);
      void set_mpi_transfer_enabled(bool transferEnabled);
      void updateSparseMinValue(const int& popID);
      Real getVelocityBlockMinValue(const int& popID) const;

      // Random number generator functions
      //char* get_rng_state_buffer();
      //random_data* get_rng_data_buffer();

      // Member variables //
      Real derivatives[fieldsolver::N_SPATIAL_CELL_DERIVATIVES];              /**< Derivatives of bulk variables in this spatial cell.*/
      Real derivativesBVOL[bvolderivatives::N_BVOL_DERIVATIVES];              /**< Derivatives of BVOL needed by the acceleration. 
                                                                               * Separate array because it does not need to be communicated.*/
      Real parameters[CellParams::N_SPATIAL_CELL_PARAMS];                     /**< Bulk variables in this spatial cell.*/
      Realf null_block_data[WID3];

      uint64_t ioLocalCellId;                                                 /**< Local cell ID used for IO, not needed elsewhere 
                                                                               * and thus not being kept up-to-date.*/
      //vmesh::LocalID mpi_number_of_blocks;                                    /**< Number of blocks in mpi_velocity_block_list.*/
      Realf* neighbor_block_data;                                             /**< Pointers for translation operator. We can point to neighbor
                                                                               * cell block data. We do not allocate memory for the pointer.*/
      vmesh::LocalID neighbor_number_of_blocks;
      uint sysBoundaryFlag;                                                   /**< What type of system boundary does the cell belong to. 
                                                                               * Enumerated in the sysboundarytype namespace's enum.*/
      uint sysBoundaryLayer;                                                  /**< Layers counted from closest systemBoundary. If 0 then it has not 
                                                                               * been computed. First sysboundary layer is layer 1.*/
      int sysBoundaryLayerNew;
      std::vector<vmesh::GlobalID> velocity_block_with_content_list;          /**< List of existing cells with content, only up-to-date after
                                                                               * call to update_has_content().*/
      vmesh::LocalID velocity_block_with_content_list_size;                   /**< Size of vector. Needed for MPI communication of size before actual list transfer.*/
      std::vector<vmesh::GlobalID> velocity_block_with_no_content_list;       /**< List of existing cells with no content, only up-to-date after
                                                                               * call to update_has_content. This is also never transferred
                                                                               * over MPI, so is invalid on remote cells.*/
      static uint64_t mpi_transfer_type;                                      /**< Which data is transferred by the mpi datatype given by spatial cells.*/
      static bool mpiTransferAtSysBoundaries;                                 /**< Do we only transfer data at boundaries (true), or in the whole system (false).*/

    private:
      SpatialCell& operator=(const SpatialCell&);
      
      bool compute_block_has_content(const vmesh::GlobalID& block,const int& popID) const;
      void merge_values_recursive(const int& popID,vmesh::GlobalID parentGID,vmesh::GlobalID blockGID,uint8_t refLevel,bool recursive,const Realf* data,
				  std::set<vmesh::GlobalID>& blockRemovalList);

      static int activePopID;
      bool initialized;
      bool mpiTransferEnabled;

      // Random number generator state variables, used for running reproducible 
      // simulations that do not depend on the number of threads of MPI processes used.
      //char rngStateBuffer[256];                                                 /**< Random number generator state buffer.*/
      //random_data rngDataBuffer;                                                /**< Random number generator data buffer.*/

      // Temporary mesh used in acceleration and propagation. 
      vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID> vmeshTemp;            /**< Temporary velocity mesh that is used in Vlasov solver.
                                                                                 * NOTE: Do not call the get-functions using this mesh as object
                                                                                 * before you have set the correct meshID using setMesh function.*/
      vmesh::VelocityBlockContainer<vmesh::LocalID> blockContainerTemp;
      std::vector<spatial_cell::Population> populations;                        /**< Particle population variables.*/
   };

   /****************************
    * Velocity block functions *
    ****************************/   
   
   template<int PAD> inline
   void SpatialCell::fetch_acc_data(const vmesh::GlobalID& blockGID,const int& dim,
                                    vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& vmesh,
                                    const Realf* src,Realf* array,Real cellSizeFractions[2]) {
      const vmesh::LocalID blockLID = vmesh.getLocalID(blockGID);
      
      #ifdef DEBUG_SPATIAL_CELL
      if (blockGID == vmesh.invalidGlobalID() || blockLID == vmesh.invalidLocalID()) {
         std::cerr << "ERROR: block has invalid global or local index " << __FILE__ << ':' << __LINE__ << std::endl;
         exit(1);
      }
      #endif

      const Realf* ptr = NULL;
      uint8_t refLevel;
      vmesh::LocalID i_block,j_block,k_block;
      vmesh.getIndices(blockGID,refLevel,i_block,j_block,k_block);

      // Copy values from x face neighbors:
      std::vector<vmesh::LocalID> nbrIDs;
      int32_t refLevelDiff;
      Real crd;
      switch (dim) {
       case 0: // Transpose i->k, j->j, k->i
         ptr = src + blockLID*WID3; // Copy values from this block
         for (int k=0; k<WID; ++k) for (int j=0; j<WID; ++j) for (int i=0; i<WID; ++i) {
            array[vblock::index(k,j,i+PAD)] = ptr[vblock::index(i,j,k)];
         }

         for (int i_nbr_off=-1; i_nbr_off<2; i_nbr_off+=2) { // Copy values from x face neighbors:
            // Get local IDs of neighbor blocks
            vmesh.getNeighborsExistingAtOffset(blockGID,i_nbr_off,+0,+0,nbrIDs,refLevelDiff);
            
            // Position that is used to interpolate values from neighbor blocks
            Real pos[3];
            if (i_nbr_off < 0) crd = WID-0.5-(PAD-1);
            else crd = 0.5;
            
            // i-index to array where interpolated values are stored
            uint32_t i_trgt = 0;
            if (i_nbr_off > 0) i_trgt = WID+PAD;
            
            if (nbrIDs.size() > 0) {     // This block has at least one existing neighbor
               if (refLevelDiff == -1) { // Neighbor is one level coarser, interpolate
                  if (nbrIDs[0] == invalid_local_id()) ptr = null_block_data; // (this check might not be necessary here)
                  else ptr = src + nbrIDs[0]*WID3;
                  for (uint32_t i=0; i<PAD; ++i) for (uint32_t k=0; k<WID; ++k) for (uint32_t j=0; j<WID; ++j) {
                     pos[0] = crd + i;
                     pos[1] = 2*(j_block%2) + j/2 + 0.5;
                     pos[2] = 2*(k_block%2) + k/2 + 0.5;
                     array[vblock::index(k,j,i_trgt+i)] = vblock::interp_xy<vblock::interpmethod::NGP>(pos,ptr);
                  }
                  cellSizeFractions[(i_nbr_off+1)/2] = 2.0;
               } else if (refLevelDiff == 0) { // Neighbor at same level, copy data
                  if (nbrIDs[0] == invalid_local_id()) ptr = null_block_data; // (this check might not be necessary here)
                  else ptr = src + nbrIDs[0]*WID3;
                  uint32_t i_src = 0;
                  if (i_nbr_off < 0) i_src = WID-PAD;
                  for (uint32_t i=0; i<PAD; ++i) for (uint32_t k=0; k<WID; ++k) for (uint32_t j=0; j<WID; ++j) {
                     array[vblock::index(k,j,i_trgt+i)] = ptr[vblock::index(i_src+i,j,k)];
                  }
                  cellSizeFractions[(i_nbr_off+1)/2] = 1.0;
               } else if (refLevelDiff == +1) { // nbr one level more refined, interpolate from four neighbors
                  for (uint32_t i=0; i<PAD; ++i) for (uint32_t k=0; k<WID; ++k) for (uint32_t j=0; j<WID; ++j) {
                     
                     int index = (k/2)*2 + j/2;
                     if (nbrIDs[index] == invalid_local_id()) ptr = null_block_data;
                     else ptr = src + nbrIDs[index]*WID3;
                     
                     pos[0] = crd + i;
                     pos[1] = 2*(j%2) + 1;
                     pos[2] = 2*(k%2) + 1;
                     array[vblock::index(k,j,i_trgt+i)] = vblock::interp_xy<vblock::interpmethod::CIC>(pos,ptr);
                  }
                  cellSizeFractions[(i_nbr_off+1)/2] = 0.5;
               }
            } else { // Neighbor does not exist, return zero values
               for (uint32_t i=0; i<PAD; ++i) for (uint32_t k=0; k<WID; ++k) for (uint32_t j=0; j<WID; ++j) {
                  array[vblock::index(k,j,i_trgt+i)] = 0.0;
               }
               cellSizeFractions[(i_nbr_off+1)/2] = 1.0;
            }
         }
         break;
       case 1: // Transpose i->i, j->k, k->j
         ptr = src + blockLID*WID3; // Copy values from this block
         for (int k=0; k<WID; ++k) for (int j=0; j<WID; ++j) for (int i=0; i<WID; ++i) {
            array[vblock::index(i,k,j+PAD)] = ptr[vblock::index(i,j,k)];
         }
         
         for (int j_nbr_off=-1; j_nbr_off<2; j_nbr_off+=2) { // Copy values from y face neighbors:
            // Get local IDs of neighbor blocks
            vmesh.getNeighborsExistingAtOffset(blockGID,+0,j_nbr_off,+0,nbrIDs,refLevelDiff);
            
            // Position that is used to interpolate values from neighbor blocks
            Real pos[3];
            if (j_nbr_off < 0) crd = WID-0.5-(PAD-1);
            else crd = 0.5;
            
            // j-index to array where interpolated values are stored
            uint32_t j_trgt = 0;
            if (j_nbr_off > 0) j_trgt = WID+PAD;
            
            if (nbrIDs.size() > 0) {     // This block has at least one existing neighbor
               if (refLevelDiff == -1) { // Neighbor is one level coarser, interpolate
                  if (nbrIDs[0] == invalid_local_id()) ptr = null_block_data; // (this check might not be necessary here)
                  else ptr = src + nbrIDs[0]*WID3;
                  for (uint32_t j=0; j<PAD; ++j) for (uint32_t k=0; k<WID; ++k) for (uint32_t i=0; i<WID; ++i) {
                     pos[0] = 2*(i_block%2) + i/2 + 0.5;
                     pos[1] = crd + j;
                     pos[2] = 2*(k_block%2) + k/2 + 0.5;
                     array[vblock::index(i,k,j_trgt+j)] = vblock::interp_xy<vblock::interpmethod::NGP>(pos,ptr);
                  }
                  cellSizeFractions[(j_nbr_off+1)/2] = 2.0;
               } else if (refLevelDiff == 0) { // Neighbor at same level, copy data
                  if (nbrIDs[0] == invalid_local_id()) ptr = null_block_data; // (this check might not be necessary here)
                  else ptr = src + nbrIDs[0]*WID3;
                  uint32_t j_src = 0;
                  if (j_nbr_off < 0) j_src = WID-PAD;
                  for (uint32_t j=0; j<PAD; ++j) for (uint32_t k=0; k<WID; ++k) for (uint32_t i=0; i<WID; ++i) {
                     array[vblock::index(i,k,j_trgt+j)] = ptr[vblock::index(i,j_src+j,k)];
                  }
                  cellSizeFractions[(j_nbr_off+1)/2] = 1.0;
               } else if (refLevelDiff == +1) { // nbr one level more refined, interpolate from four neighbors
                  for (uint32_t j=0; j<PAD; ++j) for (uint32_t k=0; k<WID; ++k) for (uint32_t i=0; i<WID; ++i) {
                     // Iterate over the four neighbors. If the neighbor does not exist, 
                     // interpolate values from the null block
                     int index = (k/2)*2 + i/2;
                     if (nbrIDs[index] == invalid_local_id()) ptr = null_block_data;
                     else ptr = src + nbrIDs[index]*WID3;
                     
                     pos[0] = 2*(i%2) + 1;
                     pos[1] = crd + j;
                     pos[2] = 2*(k%2) + 1;
                     array[vblock::index(i,k,j_trgt+j)] = vblock::interp_xy<vblock::interpmethod::CIC>(pos,ptr);
                  }
                  cellSizeFractions[(j_nbr_off+1)/2] = 0.5;
               }
            } else { // Neighbor does not exist, return zero values
               for (uint32_t j=0; j<PAD; ++j) for (uint32_t k=0; k<WID; ++k) for (uint32_t i=0; i<WID; ++i) {
                  array[vblock::index(i,k,j_trgt+j)] = 0.0;
               }
               cellSizeFractions[(j_nbr_off+1)/2] = 1.0;
            }
         }
         break;
       case 2:
         ptr = src + blockLID*WID3; // Copy values from this block
         for (int k=0; k<WID; ++k) for (int j=0; j<WID; ++j) for (int i=0; i<WID; ++i) {
            array[vblock::index(i,j,k+PAD)] = ptr[vblock::index(i,j,k)];
         }
         
         for (int k_nbr_off=-1; k_nbr_off<2; k_nbr_off+=2) { // Copy values from z face neighbors:
            // Get local IDs of neighbor blocks
            vmesh.getNeighborsExistingAtOffset(blockGID,+0,+0,k_nbr_off,nbrIDs,refLevelDiff);
            
            // Position that is used to interpolate values from neighbor blocks
            Real pos[3];
            if (k_nbr_off < 0) crd = WID-0.5-(PAD-1);
            else crd = 0.5;

            // k-index to array where interpolated values are stored
            uint32_t k_trgt = 0;
            if (k_nbr_off > 0) k_trgt = WID+PAD;
            
            if (nbrIDs.size() > 0) {     // This block has at least one existing neighbor
               if (refLevelDiff == -1) { // Neighbor is one level coarser, interpolate
                  if (nbrIDs[0] == invalid_local_id()) ptr = null_block_data; // (this check might not be necessary here)
                  else ptr = src + nbrIDs[0]*WID3;
                  for (uint32_t k=0; k<PAD; ++k) for (uint32_t j=0; j<WID; ++j) for (uint32_t i=0; i<WID; ++i) {
                     pos[0] = 2*(i_block%2) + i/2 + 0.5;
                     pos[1] = 2*(j_block%2) + j/2 + 0.5;
                     pos[2] = crd + k;
                     array[vblock::index(i,j,k_trgt+k)] = vblock::interp_xy<vblock::interpmethod::NGP>(pos,ptr);
                  }
                  cellSizeFractions[(k_nbr_off+1)/2] = 2.0;
               } else if (refLevelDiff == 0) { // Neighbor at same level, copy data
                  if (nbrIDs[0] == invalid_local_id()) ptr = null_block_data; // (this check might not be necessary here)
                  else ptr = src + nbrIDs[0]*WID3;
                  uint32_t k_src = 0;
                  if (k_nbr_off < 0) k_src = WID-PAD;
                  for (uint32_t k=0; k<PAD; ++k) for (uint32_t j=0; j<WID; ++j) for (uint32_t i=0; i<WID; ++i) {
                     array[vblock::index(i,j,k_trgt+k)] = ptr[vblock::index(i,j,k_src+k)];
                  }
                  cellSizeFractions[(k_nbr_off+1)/2] = 1.0;
               } else if (refLevelDiff == +1) { // nbr one level more refined, interpolate from four neighbors
                  for (uint32_t k=0; k<PAD; ++k) for (uint32_t j=0; j<WID; ++j) for (uint32_t i=0; i<WID; ++i) {
                     // Iterate over the four neighbors. If the neighbor does not exist, 
                     // interpolate values from the null block
                     int index = (j/2)*2 + i/2;
                     if (nbrIDs[index] == invalid_local_id()) ptr = null_block_data;
                     else ptr = src + nbrIDs[index]*WID3;
                     
                     pos[0] = 2*(i%2) + 1;
                     pos[1] = 2*(j%2) + 1;
                     pos[2] = crd + k;
                     array[vblock::index(i,j,k_trgt+k)] = vblock::interp_xy<vblock::interpmethod::CIC>(pos,ptr);
                  }
                  cellSizeFractions[(k_nbr_off+1)/2] = 0.5;
               }
            } else { // Neighbor does not exist, return zero values
               for (uint32_t k=0; k<PAD; ++k) for (uint32_t j=0; j<WID; ++j) for (uint32_t i=0; i<WID; ++i) {
                  array[vblock::index(i,j,k_trgt+k)] = 0.0;
               }
               cellSizeFractions[(k_nbr_off+1)/2] = 1.0;
            }
         }
	 break;

      } // end switch
   }
   
   template<int PAD> inline
   void SpatialCell::fetch_data(const vmesh::GlobalID& blockGID,
                                const vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& vmesh,
                                const Realf* src,Realf* array) {
      //const vmesh::LocalID blockLID = get_velocity_block_local_id(blockGID);
      const vmesh::LocalID blockLID = vmesh.getLocalID(blockGID);

      if (blockLID == invalid_local_id()) {
         std::cerr << "ERROR: invalid local id in " << __FILE__ << ' ' << __LINE__ << std::endl;
         exit(1);
      }

      // Copy values from this block:
      const Realf* ptr = src + blockLID*WID3;
      for (int k=0; k<WID; ++k) for (int j=0; j<WID; ++j) for (int i=0; i<WID; ++i) {
         array[vblock::padIndex<PAD>(i+PAD,j+PAD,k+PAD)] = ptr[vblock::index(i,j,k)];
      }

      uint8_t refLevel;
      vmesh::LocalID i_block,j_block,k_block;
      vmesh.getIndices(blockGID,refLevel,i_block,j_block,k_block);

      // Copy values from x face neighbors:
      std::vector<vmesh::LocalID> nbrIDs;
      int32_t refLevelDiff;

      Real crd;
      for (int i_nbr_off=-1; i_nbr_off<2; i_nbr_off+=2) {
         vmesh.getNeighborsExistingAtOffset(blockGID,i_nbr_off,+0,+0,nbrIDs,refLevelDiff);
         Real pos[3];
         if (i_nbr_off < 0) crd = WID-0.5-(PAD-1);
         else crd = 0.5;

         uint32_t i_trgt = 0;
         if (i_nbr_off > 0) i_trgt = WID+PAD;

         if (nbrIDs.size() > 0) {
            if (refLevelDiff == -1) { // nbr one level coarser, interpolate
               if (nbrIDs[0] == invalid_local_id()) ptr = null_block_data;
               else ptr = src + nbrIDs[0]*WID3;
               for (uint32_t i=0; i<PAD; ++i) for (uint32_t k=0; k<WID; ++k) for (uint32_t j=0; j<WID; ++j) {
                  pos[0] = crd + i;
                  pos[1] = 2*(j_block%2) + j/2 + 0.5;
                  pos[2] = 2*(k_block%2) + k/2 + 0.5;
                  array[vblock::padIndex<PAD>(i_trgt+i,j+PAD,k+PAD)] = vblock::interp_yz<vblock::interpmethod::NGP>(pos,ptr);
               }
            } else if (refLevelDiff == 0) { // nbr at same level, simple data copy
               if (nbrIDs[0] == invalid_local_id()) ptr = null_block_data;
               else ptr = src + nbrIDs[0]*WID3;
               uint32_t i_src = 0;
               if (i_nbr_off < 0) i_src = WID-PAD;
               for (uint32_t i=0; i<PAD; ++i) for (uint32_t k=0; k<WID; ++k) for (uint32_t j=0; j<WID; ++j) {
                  array[vblock::padIndex<PAD>(i_trgt+i,j+PAD,k+PAD)] = ptr[vblock::index(i_src+i,j,k)];
               }
            } else if (refLevelDiff == +1) { // nbr one level more refined, interpolate from four neighbors
               for (uint32_t i=0; i<PAD; ++i) for (uint32_t k=0; k<WID; ++k) for (uint32_t j=0; j<WID; ++j) {
                  int index = (k/2)*2 + j/2;
                  if (nbrIDs[index] == invalid_local_id()) ptr = null_block_data;
                  else ptr = src + nbrIDs[index]*WID3;

                  pos[0] = crd + i;
                  pos[1] = 2*(j%2) + 1;
                  pos[2] = 2*(k%2) + 1;
                  array[vblock::padIndex<PAD>(i_trgt+i,j+PAD,k+PAD)] = vblock::interp_yz<vblock::interpmethod::CIC>(pos,ptr);
               }
            }
         } else { // Neighbor does not exist, return zero values
            for (uint32_t i=0; i<PAD; ++i) for (uint32_t k=0; k<WID; ++k) for (uint32_t j=0; j<WID; ++j) {
               array[vblock::padIndex<PAD>(i_trgt+i,j+PAD,k+PAD)] = 0.0;
            }
         }
      }

      // Copy values from y face neighbors:
      for (int j_nbr_off=-1; j_nbr_off<2; j_nbr_off+=2) {
         vmesh.getNeighborsExistingAtOffset(blockGID,+0,j_nbr_off,+0,nbrIDs,refLevelDiff);
         Real pos[3];
         if (j_nbr_off < 0) crd = WID-0.5-(PAD-1);
         else crd = 0.5;
         
         uint32_t j_trgt = 0;
         if (j_nbr_off > 0) j_trgt = WID+PAD;
         
         if (nbrIDs.size() > 0) {
            if (refLevelDiff == -1) { // nbr one level coarser, interpolate
               if (nbrIDs[0] == invalid_local_id()) ptr = null_block_data;
               else ptr = src + nbrIDs[0]*WID3;
               for (uint32_t j=0; j<PAD; ++j) for (uint32_t k=0; k<WID; ++k) for (uint32_t i=0; i<WID; ++i) {
                  pos[0] = 2*(i_block%2) + i/2 + 0.5;
                  pos[1] = crd + j;
                  pos[2] = 2*(k_block%2) + k/2 + 0.5;
                  array[vblock::padIndex<PAD>(i+PAD,j_trgt+j,k+PAD)] = vblock::interp_xz<vblock::interpmethod::NGP>(pos,ptr);
               }
            } else if (refLevelDiff == 0) { // nbr at same level, simple data copy
               if (nbrIDs[0] == invalid_local_id()) ptr = null_block_data;
               else ptr = src + nbrIDs[0]*WID3;
               uint32_t j_src = 0;
               if (j_nbr_off < 0) j_src = WID-PAD;
               for (uint32_t j=0; j<PAD; ++j) for (uint32_t k=0; k<WID; ++k) for (uint32_t i=0; i<WID; ++i) {
                  array[vblock::padIndex<PAD>(i+PAD,j_trgt+j,k+PAD)] = ptr[vblock::index(i,j_src+j,k)];
               }
            } else if (refLevelDiff == +1) { // nbr one level more refined, interpolate from four neighbors
               for (uint32_t j=0; j<PAD; ++j) for (uint32_t k=0; k<WID; ++k) for (uint32_t i=0; i<WID; ++i) {
                  int index = (k/2)*2 + i/2;
                  if (nbrIDs[index] == invalid_local_id()) ptr = null_block_data;
                  else ptr = src + nbrIDs[index]*WID3;

                  pos[0] = 2*(i%2) + 1;
                  pos[1] = crd + j;
                  pos[2] = 2*(k%2) + 1;
                  array[vblock::padIndex<PAD>(i+PAD,j_trgt+j,k+PAD)] = vblock::interp_xz<vblock::interpmethod::CIC>(pos,ptr);
               }
            }
         } else { // Neighbor does not exist, return zero values
            for (uint32_t j=0; j<PAD; ++j) for (uint32_t k=0; k<WID; ++k) for (uint32_t i=0; i<WID; ++i) {
               array[vblock::padIndex<PAD>(i+PAD,j_trgt+j,k+PAD)] = 0.0;
            }
         }
      }

      // Copy values from z face neighbors:
      for (int k_nbr_off=-1; k_nbr_off<2; k_nbr_off+=2) {
         vmesh.getNeighborsExistingAtOffset(blockGID,+0,+0,k_nbr_off,nbrIDs,refLevelDiff);
         Real pos[3];
         uint32_t k_trgt = 0;
         if (k_nbr_off > 0) k_trgt = WID+PAD;
         
         if (k_nbr_off < 0) crd = WID-0.5-(PAD-1);
         else crd = 0.5;
         
         if (nbrIDs.size() > 0) {
            if (refLevelDiff == -1) { // nbr one level coarser, interpolate
               if (nbrIDs[0] == invalid_local_id()) ptr = null_block_data;
               else ptr = src + nbrIDs[0]*WID3;
               for (uint32_t k=0; k<PAD; ++k) for (uint32_t j=0; j<WID; ++j) for (uint32_t i=0; i<WID; ++i) {
                  pos[0] = 2*(i_block%2) + i/2 + 0.5;
                  pos[1] = 2*(j_block%2) + j/2 + 0.5;
                  pos[2] = crd + k;
                  array[vblock::padIndex<PAD>(i+PAD,j+PAD,k_trgt+k)] = vblock::interp_xy<vblock::interpmethod::NGP>(pos,ptr);
               }
            } else if (refLevelDiff == 0) { // nbr at same level, simple data copy
               if (nbrIDs[0] == invalid_local_id()) ptr = null_block_data;
               else ptr = src + nbrIDs[0]*WID3;
               uint32_t k_src = 0;
               if (k_nbr_off < 0) k_src = WID-PAD;
               for (uint32_t k=0; k<PAD; ++k) for (uint32_t j=0; j<WID; ++j) for (uint32_t i=0; i<WID; ++i) {
                  array[vblock::padIndex<PAD>(i+PAD,j+PAD,k_trgt+k)] = ptr[vblock::index(i,j,k_src+k)];
               }
            } else if (refLevelDiff == +1) { // nbr one level more refined, interpolate from four neighbors
               for (uint32_t k=0; k<PAD; ++k) for (uint32_t j=0; j<WID; ++j) for (uint32_t i=0; i<WID; ++i) {
                  int index = (j/2)*2 + i/2;
                  if (nbrIDs[index] == invalid_local_id()) ptr = null_block_data;
                  else ptr = src + nbrIDs[index]*WID3;
                  
                  pos[0] = 2*(i%2) + 1;
                  pos[1] = 2*(j%2) + 1;
                  pos[2] = crd + k;
                  array[vblock::padIndex<PAD>(i+PAD,j+PAD,k_trgt+k)] = vblock::interp_xy<vblock::interpmethod::CIC>(pos,ptr);
               }
            }
         } else { // Neighbor does not exist, return zero values
            for (uint32_t k=0; k<PAD; ++k) for (uint32_t j=0; j<WID; ++j) for (uint32_t i=0; i<WID; ++i) {
               array[vblock::padIndex<PAD>(i+PAD,j+PAD,k_trgt+k)] = 0.0;
            }
         }
      }
   }
   
   inline vmesh::GlobalID SpatialCell::find_velocity_block(uint8_t& refLevel,vmesh::GlobalID cellIndices[3],const int& popID) {
      #ifdef DEBUG_SPATIAL_CELL
      if (popID >= populations.size()) {
         std::cerr << "ERROR, popID " << popID << " exceeds populations.size() " << populations.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;             
         exit(1);
      }
      #endif
      return populations[popID].vmesh.findBlock(refLevel,cellIndices);
   }

   inline Realf* SpatialCell::get_data(const int& popID) {
      #ifdef DEBUG_SPATIAL_CELL
      if (popID >= populations.size()) {
         std::cerr << "ERROR, popID " << popID << " exceeds populations.size() " << populations.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;             
         exit(1);
      }
      #endif
      return populations[popID].blockContainer.getData();
   }
   
   inline const Realf* SpatialCell::get_data(const int& popID) const {
      #ifdef DEBUG_SPATIAL_CELL
      if (popID >= populations.size()) {
         std::cerr << "ERROR, popID " << popID << " exceeds populations.size() " << populations.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;             
         exit(1);
      }
      #endif
      return populations[popID].blockContainer.getData();
   }

   inline Realf* SpatialCell::get_data(const vmesh::LocalID& blockLID,const int& popID) {
      #ifdef DEBUG_SPATIAL_CELL
      if (popID >= populations.size()) {
         std::cerr << "ERROR, popID " << popID << " exceeds populations.size() " << populations.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;             
         exit(1);
      }
      if (blockLID >= populations[popID].blockContainer.size()) {
         std::cerr << "ERROR, block LID out of bounds, blockContainer.size() " << populations[popID].blockContainer.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
         exit(1);
      }
      #endif
      if (blockLID == vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>::invalidLocalID()) return null_block_data;
      return populations[popID].blockContainer.getData(blockLID);
   }
   
   inline const Realf* SpatialCell::get_data(const vmesh::LocalID& blockLID,const int& popID) const {
      #ifdef DEBUG_SPATIAL_CELL
      if (popID >= populations.size()) {
         std::cerr << "ERROR, popID " << popID << " exceeds populations.size() " << populations.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;             
         exit(1);
      }
      if (blockLID >= populations[popID].blockContainer.size()) {
         std::cerr << "ERROR, block LID out of bounds, blockContainer.size() " << populations[popID].blockContainer.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
         exit(1);
      }
      #endif
      if (blockLID == vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>::invalidLocalID()) return null_block_data;
      return populations[popID].blockContainer.getData(blockLID);
   }

   inline Real* SpatialCell::get_block_parameters(const int& popID) {
      #ifdef DEBUG_SPATIAL_CELL
      if (popID >= populations.size()) {
         std::cerr << "ERROR, popID " << popID << " exceeds populations.size() " << populations.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;             
         exit(1);
      }
      #endif
      return populations[popID].blockContainer.getParameters();
   }
   
   inline const Real* SpatialCell::get_block_parameters(const int& popID) const {
      #ifdef DEBUG_SPATIAL_CELL
      if (popID >= populations.size()) {
         std::cerr << "ERROR, popID " << popID << " exceeds populations.size() " << populations.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;             
         exit(1);
      }
      #endif
      return populations[popID].blockContainer.getParameters();
   }
   
   inline Real* SpatialCell::get_block_parameters(const vmesh::LocalID& blockLID,const int& popID) {
      #ifdef DEBUG_SPATIAL_CELL
      if (popID >= populations.size()) {
         std::cerr << "ERROR, popID " << popID << " exceeds populations.size() " << populations.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;             
         exit(1);
      }
      if (blockLID >= populations[popID].blockContainer.size()) {
         std::cerr << "ERROR, block LID out of bounds, blockContainer.size() " << populations[popID].blockContainer.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
         exit(1);
      }
      #endif
      return populations[popID].blockContainer.getParameters(blockLID);
   }
   
   inline const Real* SpatialCell::get_block_parameters(const vmesh::LocalID& blockLID,const int& popID) const {
      #ifdef DEBUG_SPATIAL_CELL
      if (popID >= populations.size()) {
         std::cerr << "ERROR, popID " << popID << " exceeds populations.size() " << populations.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;             
         exit(1);
      }
      if (blockLID >= populations[popID].blockContainer.size()) {
         std::cerr << "ERROR, block LID out of bounds, blockContainer.size() " << populations[popID].blockContainer.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
         exit(1);
      }
      #endif
      return populations[popID].blockContainer.getParameters(blockLID);
   }
   
   inline Real* SpatialCell::get_cell_parameters() {
       return parameters;
   }

   inline const Real* SpatialCell::get_cell_parameters() const {
      return parameters;
   }

   inline uint8_t SpatialCell::get_maximum_refinement_level(const int& popID) {
      return populations[popID].vmesh.getMaxAllowedRefinementLevel();
   }

   inline vmesh::LocalID SpatialCell::get_number_of_velocity_blocks(const int& popID) const {
      #ifdef DEBUG_SPATIAL_CELL
      if (popID >= populations.size()) {
         std::cerr << "ERROR, popID " << popID << " exceeds populations.size() " << populations.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;             
         exit(1);
      }
      #endif
      return populations[popID].blockContainer.size();
   }

    /** Get the total number of velocity blocks in this cell, summed over 
     * all existing particle populations.
     * @return Total number of velocity blocks in the cell.*/
    inline vmesh::LocalID SpatialCell::get_number_of_all_velocity_blocks() const {
        vmesh::LocalID N_blocks = 0;
        for (size_t p=0; p<populations.size(); ++p)
            N_blocks += populations[p].blockContainer.size();
        return N_blocks;
    }
    
   inline int SpatialCell::get_number_of_populations() const {
      return populations.size();
   }

   inline const vmesh::LocalID* SpatialCell::get_velocity_grid_length(const int& popID,const uint8_t& refLevel) {
      return populations[popID].vmesh.getGridLength(refLevel);
   }

   inline const Real* SpatialCell::get_velocity_grid_block_size(const int& popID,const uint8_t& refLevel) {
      return populations[popID].vmesh.getBlockSize(refLevel);
   }

   inline const Real* SpatialCell::get_velocity_grid_cell_size(const int& popID,const uint8_t& refLevel) {
      return populations[popID].vmesh.getCellSize(refLevel);
   }

   inline void SpatialCell::get_velocity_block_coordinates(const int& popID,const vmesh::GlobalID& globalID,Real* coords) {
      populations[popID].vmesh.getBlockCoordinates(globalID,coords);
   }
   
   /*!
    Returns the indices of given velocity block
    */
   inline velocity_block_indices_t SpatialCell::get_velocity_block_indices(const int& popID,const vmesh::GlobalID block) {
      velocity_block_indices_t indices;
      uint8_t refLevel;
      populations[popID].vmesh.getIndices(block,refLevel,indices[0],indices[1],indices[2]);
      return indices;
   }

   inline velocity_block_indices_t SpatialCell::get_velocity_block_indices(const int& popID,const vmesh::GlobalID block,uint8_t& refLevel) {
      velocity_block_indices_t indices;
      populations[popID].vmesh.getIndices(block,refLevel,indices[0],indices[1],indices[2]);
      return indices;
   }

   /*!
    Returns the velocity block at given indices or error_velocity_block
    */
   inline vmesh::GlobalID SpatialCell::get_velocity_block(const int& popID,const velocity_block_indices_t indices,const uint8_t& refLevel) const {
      return populations[popID].vmesh.getGlobalID(refLevel,indices[0],indices[1],indices[2]);
   }
   
   inline vmesh::GlobalID SpatialCell::get_velocity_block(const int& popID,vmesh::GlobalID blockIndices[3],const uint8_t& refLevel) const {
      return populations[popID].vmesh.getGlobalID(refLevel,blockIndices[0],blockIndices[1],blockIndices[2]);
   }
   
   /*!
    Returns the velocity block at given location or
    error_velocity_block if outside of the velocity grid
    */
   inline vmesh::GlobalID SpatialCell::get_velocity_block(const int& popID,const Real vx,const Real vy,const Real vz,const uint8_t& refLevel) const {
      Real coords[3] = {vx,vy,vz};
      return populations[popID].vmesh.getGlobalID(refLevel,coords);
   }
   
   inline vmesh::GlobalID SpatialCell::get_velocity_block(const int& popID,const Real* coords,const uint8_t& refLevel) const {
      return populations[popID].vmesh.getGlobalID(refLevel,coords);
   }
   
   inline vmesh::GlobalID SpatialCell::get_velocity_block_child(const int& popID,const vmesh::GlobalID& blockGID,const uint8_t& refLevel,
                                                                const int& i_cell,const int& j_cell,const int& k_cell) {
      uint8_t ref = refLevel;

      vmesh::LocalID i_child,j_child,k_child;
      i_child = 2*i_child + i_cell/2;
      j_child = 2*j_child + j_cell/2;
      k_child = 2*k_child + k_cell/2;

      while (ref != populations[popID].vmesh.getMaxAllowedRefinementLevel()) {
         vmesh::LocalID i_child,j_child,k_child;
         populations[popID].vmesh.getIndices(blockGID,ref,i_child,j_child,k_child);
         
         return populations[popID].vmesh.getGlobalID(refLevel+1,i_child,j_child,k_child);
      }
      return populations[popID].vmesh.invalidGlobalID();
   }
   
   inline void SpatialCell::get_velocity_block_children_local_ids(
         const vmesh::GlobalID& blockGID,
         std::vector<vmesh::LocalID>& childrenLIDs,
         const int& popID) {
      #ifdef DEBUG_SPATIAL_CELL
      if (popID >= populations.size()) {
         std::cerr << "ERROR, popID " << popID << " exceeds populations.size() " << populations.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;             
         exit(1);
      }
      #endif
      
      std::vector<vmesh::GlobalID> childrenGIDs;
      populations[popID].vmesh.getChildren(blockGID,childrenGIDs);
      childrenLIDs.resize(childrenGIDs.size());
      for (size_t c=0; c<childrenGIDs.size(); ++c) 
          childrenLIDs[c] = populations[popID].vmesh.getLocalID(childrenGIDs[c]);
   }

   inline vmesh::GlobalID SpatialCell::get_velocity_block_parent(const int& popID,const vmesh::GlobalID& blockGID) {
      return populations[popID].vmesh.getParent(blockGID);
   }

   inline vmesh::GlobalID SpatialCell::get_velocity_block_global_id(const vmesh::LocalID& blockLID,const int& popID) const {
      #ifdef DEBUG_SPATIAL_CELL
      if (popID >= populations.size()) {
         std::cerr << "ERROR, popID " << popID << " exceeds populations.size() " << populations.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;             
         exit(1);
      }
      #endif

      return populations[popID].vmesh.getGlobalID(blockLID);
   }
      
   inline vmesh::LocalID SpatialCell::get_velocity_block_local_id(const vmesh::GlobalID& blockGID,const int& popID) const {
      #ifdef DEBUG_SPATIAL_CELL
      if (popID >= populations.size()) {
         std::cerr << "ERROR, popID " << popID << " exceeds populations.size() " << populations.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;             
         exit(1);
      }
      #endif
      
      return populations[popID].vmesh.getLocalID(blockGID);
   }

   inline void SpatialCell::get_velocity_block_size(const int& popID,const vmesh::GlobalID block,Real blockSize[3]) {
      populations[popID].vmesh.getBlockSize(block,blockSize);
   }
   
   /*!
    Returns the edge where given velocity block starts.
    */
   inline Real SpatialCell::get_velocity_block_vx_min(const int& popID,const vmesh::GlobalID block) const {
      Real coords[3];
      populations[popID].vmesh.getBlockCoordinates(block,coords);
      return coords[0];
   }

   /*!
    Returns the edge where given velocity block ends.
    */
   inline Real SpatialCell::get_velocity_block_vx_max(const int& popID,const vmesh::GlobalID block) const {
      Real coords[3];
      populations[popID].vmesh.getBlockCoordinates(block,coords);
      
      Real size[3];
      populations[popID].vmesh.getBlockSize(block,size);
      return coords[0]+size[0];
   }

   /*!
    Returns the edge where given velocity block starts.
    */
   inline Real SpatialCell::get_velocity_block_vy_min(const int& popID,const vmesh::GlobalID block) const {
      Real coords[3];
      populations[popID].vmesh.getBlockCoordinates(block,coords);
      return coords[1];
   }

   /*!
    Returns the edge where given velocity block ends.
    */
   inline Real SpatialCell::get_velocity_block_vy_max(const int& popID,const vmesh::GlobalID block) const {
      Real coords[3];
      populations[popID].vmesh.getBlockCoordinates(block,coords);
      
      Real size[3];
      populations[popID].vmesh.getBlockSize(block,size);
      return coords[1]+size[1];
   }

   /*!
    Returns the edge where given velocity block starts.
    */
   inline Real SpatialCell::get_velocity_block_vz_min(const int& popID,const vmesh::GlobalID block) const {
      Real coords[3];
      populations[popID].vmesh.getBlockCoordinates(block,coords);
      return coords[2];
   }
   
   /*!
    Returns the edge where given velocity block ends.
    */
   inline Real SpatialCell::get_velocity_block_vz_max(const int& popID,const vmesh::GlobalID block) const {
      Real coords[3];
      populations[popID].vmesh.getBlockCoordinates(block,coords);
      
      Real size[3];
      populations[popID].vmesh.getBlockSize(block,size);
      return coords[2]+size[2];
   }

      /***************************
       * Velocity cell functions *
       ***************************/

   /*!
    Returns the indices of given velocity cell
    */
   inline velocity_cell_indices_t SpatialCell::get_velocity_cell_indices(const unsigned int cell) const {
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
   inline unsigned int SpatialCell::get_velocity_cell(const velocity_cell_indices_t indices) const {
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
   inline unsigned int SpatialCell::get_velocity_cell(const int& popID,
                                                      const vmesh::GlobalID velocity_block,
                                                      const Real vx,
                                                      const Real vy,
                                                      const Real vz
                                                     ) const {
      const Real block_vx_min = get_velocity_block_vx_min(popID,velocity_block);
      const Real block_vx_max = get_velocity_block_vx_max(popID,velocity_block);
      const Real block_vy_min = get_velocity_block_vy_min(popID,velocity_block);
      const Real block_vy_max = get_velocity_block_vy_max(popID,velocity_block);
      const Real block_vz_min = get_velocity_block_vz_min(popID,velocity_block);
      const Real block_vz_max = get_velocity_block_vz_max(popID,velocity_block);
      
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
      
      return get_velocity_cell(indices);
   }

   /*!
    Returns the edge where given velocity cell in the given velocity block starts.
    TODO: move these to velocity cell class?
    */
   inline Real SpatialCell::get_velocity_cell_vx_min(const int& popID,
      const vmesh::GlobalID velocity_block,
      const unsigned int velocity_cell
   ) const {
      if (velocity_cell == error_velocity_cell) {
         return std::numeric_limits<Real>::quiet_NaN();
      }
      
      const velocity_cell_indices_t indices = get_velocity_cell_indices(velocity_cell);
      if (indices[0] == error_velocity_cell_index) {
         return std::numeric_limits<Real>::quiet_NaN();
      }
      
      const Real block_vx_min = get_velocity_block_vx_min(popID,velocity_block);
      const Real block_vx_max = get_velocity_block_vx_max(popID,velocity_block);
      
      return block_vx_min + (block_vx_max - block_vx_min) / block_vx_length * indices[0];
   }

   /*!
    Returns the edge where given velocity cell in the given velocity block ends.
    */
   inline Real SpatialCell::get_velocity_cell_vx_max(const int& popID,
      const vmesh::GlobalID velocity_block,
      const unsigned int velocity_cell
   ) const {
      if (velocity_cell == error_velocity_cell) {
         return std::numeric_limits<Real>::quiet_NaN();
      }
      
      const velocity_cell_indices_t indices = get_velocity_cell_indices(velocity_cell);
      if (indices[0] == error_velocity_cell_index) {
         return std::numeric_limits<Real>::quiet_NaN();
      }
      
      const Real block_vx_min = get_velocity_block_vx_min(popID,velocity_block);
      const Real block_vx_max = get_velocity_block_vx_max(popID,velocity_block);

      return block_vx_min + (block_vx_max - block_vx_min) / block_vx_length * (indices[0] + 1);
   }

   /*!
    Returns the edge where given velocity cell in the given velocity block starts.
    */
   inline Real SpatialCell::get_velocity_cell_vy_min(const int& popID,
      const vmesh::GlobalID velocity_block,
      const unsigned int velocity_cell
   ) const {
      if (velocity_cell == error_velocity_cell) {
         return std::numeric_limits<Real>::quiet_NaN();
      }
      
      const velocity_cell_indices_t indices = get_velocity_cell_indices(velocity_cell);
      if (indices[1] == error_velocity_cell_index) {
         return std::numeric_limits<Real>::quiet_NaN();
      }
      
      const Real block_vy_min = get_velocity_block_vy_min(popID,velocity_block);
      const Real block_vy_max = get_velocity_block_vy_max(popID,velocity_block);
      
      return block_vy_min + (block_vy_max - block_vy_min) / block_vy_length * indices[1];
   }

   /*!
    Returns the edge where given velocity cell in the given velocity block ends.
    */
   inline Real SpatialCell::get_velocity_cell_vy_max(const int& popID,
      const vmesh::GlobalID velocity_block,
      const unsigned int velocity_cell
   ) const {
      if (velocity_cell == error_velocity_cell) {
         return std::numeric_limits<Real>::quiet_NaN();
      }
      
      const velocity_cell_indices_t indices = get_velocity_cell_indices(velocity_cell);
      if (indices[1] == error_velocity_cell_index) {
         return std::numeric_limits<Real>::quiet_NaN();
      }
      
      const Real block_vy_min = get_velocity_block_vy_min(popID,velocity_block);
      const Real block_vy_max = get_velocity_block_vy_max(popID,velocity_block);
      
      return block_vy_min + (block_vy_max - block_vy_min) / block_vy_length * (indices[1] + 1);
   }

   /*!
    Returns the edge where given velocity cell in the given velocity block starts.
    */
   inline Real SpatialCell::get_velocity_cell_vz_min(const int& popID,
      const vmesh::GlobalID velocity_block,
      const unsigned int velocity_cell
   ) const {
      if (velocity_cell == error_velocity_cell) {
         return std::numeric_limits<Real>::quiet_NaN();
      }
      
      const velocity_cell_indices_t indices = get_velocity_cell_indices(velocity_cell);
      if (indices[2] == error_velocity_cell_index) {
         return std::numeric_limits<Real>::quiet_NaN();
      }
      
      const Real block_vz_min = get_velocity_block_vz_min(popID,velocity_block);
      const Real block_vz_max = get_velocity_block_vz_max(popID,velocity_block);
      
      return block_vz_min + (block_vz_max - block_vz_min) / block_vz_length * indices[2];
   }

   /*!
    Returns the edge where given velocity cell in the given velocity block ends.
    */
   inline Real SpatialCell::get_velocity_cell_vz_max(const int& popID,
      const vmesh::GlobalID velocity_block,
      const unsigned int velocity_cell
   ) const {
      if (velocity_cell == error_velocity_cell) {
         return std::numeric_limits<Real>::quiet_NaN();
      }
      
      const velocity_cell_indices_t indices = get_velocity_cell_indices(velocity_cell);
      if (indices[2] == error_velocity_cell_index) {
         return std::numeric_limits<Real>::quiet_NaN();
      }
      
      const Real block_vz_min = get_velocity_block_vz_min(popID,velocity_block);
      const Real block_vz_max = get_velocity_block_vz_max(popID,velocity_block);

      return block_vz_min + (block_vz_max - block_vz_min) / block_vz_length * (indices[2] + 1);
   }
   
   inline const Real* SpatialCell::get_velocity_grid_min_limits(const int& popID) {
      return populations[popID].vmesh.getMeshMinLimits();
   }
   
   inline const Real* SpatialCell::get_velocity_grid_max_limits(const int& popID) {
      return populations[popID].vmesh.getMeshMaxLimits();
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

   /*!
    Returns the number of given velocity blocks that exist.
    */
   inline size_t SpatialCell::count(const vmesh::GlobalID& block,const int& popID) const {
      #ifdef DEBUG_SPATIAL_CELL
      if (popID >= populations.size()) {
         std::cerr << "ERROR, popID " << popID << " exceeds populations.size() " << populations.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;             
         exit(1);
      }
      #endif

      return populations[popID].vmesh.count(block);
   }

   /*!
    Returns the number of existing velocity blocks.
    */
   inline size_t SpatialCell::size(const int& popID) const {
      #ifdef DEBUG_SPATIAL_CELL
      if (popID >= populations.size()) {
         std::cerr << "ERROR, popID " << popID << " exceeds populations.size() " << populations.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;             
         exit(1);
      }
      #endif

      return populations[popID].vmesh.size();
   }

   /*!
    Sets the given value to a velocity cell at given coordinates.
    * 
    Creates the velocity block at given coordinates if it doesn't exist.
    */
   inline void SpatialCell::set_value(const Real vx,const Real vy,const Real vz,
                                      const Realf value,const int& popID) {
      #ifdef DEBUG_SPATIAL_CELL
      if (popID >= populations.size()) {
         std::cerr << "ERROR, popID " << popID << " exceeds populations.size() " << populations.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;             
         exit(1);
      }
      #endif
      
      const vmesh::GlobalID blockGID = get_velocity_block(popID,vx, vy, vz);
      vmesh::LocalID blockLID = populations[popID].vmesh.getLocalID(blockGID);
      if (blockLID == vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>::invalidLocalID()) {
         if (!add_velocity_block(blockGID,popID)) {
            std::cerr << "Couldn't add velocity block " << blockGID << std::endl;
            abort();
         }
         blockLID = populations[popID].vmesh.getLocalID(blockGID);
      }

      const unsigned int cell = get_velocity_cell(popID,blockGID, vx, vy, vz);
      get_data(blockLID,popID)[cell] = value;
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
   inline void SpatialCell::set_value(const vmesh::GlobalID& blockGID,const unsigned int cell,
                                      const Realf value,const int& popID) {
      #ifdef DEBUG_SPATIAL_CELL
      if (popID >= populations.size()) {
         std::cerr << "ERROR, popID " << popID << " exceeds populations.size() " << populations.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;             
         exit(1);
      }
      #endif
      
      vmesh::LocalID blockLID = populations[popID].vmesh.getLocalID(blockGID);
      if (blockLID == vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>::invalidLocalID()) {
         if (!add_velocity_block(blockGID,popID)) {
            std::cerr << "Couldn't add velocity block " << blockGID << std::endl;
            abort();
         }
         blockLID = populations[popID].vmesh.getLocalID(blockGID);
      }

      get_data(blockLID,popID)[cell] = value;
   }

   /*!
    Increments the value of velocity cell at given coordinate-
    * 
    Creates the velocity block at given coordinates if it doesn't exist.
    */
   inline void SpatialCell::increment_value(const Real vx,const Real vy,const Real vz,
                                            const Realf value,const int& popID) {
      #ifdef DEBUG_SPATIAL_CELL
      if (popID >= populations.size()) {
         std::cerr << "ERROR, popID " << popID << " exceeds populations.size() " << populations.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;             
         exit(1);
      }
      #endif

      const vmesh::GlobalID blockGID = SpatialCell::get_velocity_block(popID,vx,vy,vz);
      vmesh::LocalID blockLID = get_velocity_block_local_id(blockGID,popID);

      if (blockLID == SpatialCell::invalid_local_id()) {
          if (!add_velocity_block(blockGID,popID)) {
             std::cerr << "Couldn't add velocity block " << blockGID << std::endl;
             abort();
          }
          blockLID = get_velocity_block_local_id(blockGID,popID);
      }

      if (blockLID == invalid_local_id()) {
         std::cerr << __FILE__ << ":" << __LINE__
           << " block_ptr == NULL" << std::endl;
         abort();
      }

      const unsigned int cell = get_velocity_cell(popID,blockGID,vx,vy,vz);
      get_data(blockLID,popID)[cell] += value;
   }

   /*!
    Increments the value of velocity cell at given index
    * 
    Creates the velocity block if it doesn't exist.
    */
   inline void SpatialCell::increment_value(const vmesh::GlobalID& blockGID,
                                            const unsigned int cell,
                                            const Realf value,const int& popID) {
      #ifdef DEBUG_SPATIAL_CELL
      if (popID >= populations.size()) {
         std::cerr << "ERROR, popID " << popID << " exceeds populations.size() " << populations.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;             
         exit(1);
      }
      #endif

      vmesh::LocalID blockLID = populations[popID].vmesh.getLocalID(blockGID);
      if (blockLID == vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>::invalidLocalID()) {
         if (!add_velocity_block(blockGID,popID)) {
            std::cerr << "Couldn't add velocity block " << blockGID << std::endl;
            abort();
         }
         blockLID = populations[popID].vmesh.getLocalID(blockGID);
      }

      get_data(blockLID,popID)[cell] += value;
   }

   inline vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& SpatialCell::get_velocity_mesh(const size_t& popID) {
      #ifdef DEBUG_SPATIAL_CELL
      if (popID >= populations.size()) {
         std::cerr << "ERROR, popID " << popID << " exceeds populations.size() " << populations.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;             
         exit(1);
      }
      #endif

      return populations[popID].vmesh;
   }

   inline vmesh::VelocityBlockContainer<vmesh::LocalID>& SpatialCell::get_velocity_blocks(const size_t& popID) {
      #ifdef DEBUG_SPATIAL_CELL
      if (popID >= populations.size()) {
         std::cerr << "ERROR, popID " << popID << " exceeds populations.size() " << populations.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;             
         exit(1);
      }
      #endif
      
      return populations[popID].blockContainer;
   }

   inline vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& SpatialCell::get_velocity_mesh_temporary() {
      return vmeshTemp;
   }
   
   inline vmesh::VelocityBlockContainer<vmesh::LocalID>& SpatialCell::get_velocity_blocks_temporary() {
      return blockContainerTemp;
   }

   /*!
    * Gets the value of a velocity cell at given coordinates.
    * 
    * Returns 0 if it doesn't exist.
    */
   inline Realf SpatialCell::get_value(const Real vx, const Real vy, const Real vz,const int& popID) const {
      #ifdef DEBUG_SPATIAL_CELL
      if (popID >= populations.size()) {
         std::cerr << "ERROR, popID " << popID << " exceeds populations.size() " << populations.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;             
         exit(1);
      }
      #endif
      const vmesh::GlobalID blockGID = get_velocity_block(popID, vx, vy, vz);
      if (count(blockGID,popID) == 0) {
         return 0.0;
      }
      const vmesh::LocalID blockLID = get_velocity_block_local_id(blockGID, popID);
      //const Velocity_Block block_ptr = at(block);
      if (blockLID == invalid_local_id()) {
      //if (block_ptr.is_null() == true) {
         std::cerr << __FILE__ << ":" << __LINE__
            << " block_ptr == NULL" << std::endl; 
         abort();
      }

      const unsigned int cell = get_velocity_cell(popID, blockGID, vx, vy, vz);
      
      return get_data(blockLID,popID)[cell];
   }

   inline Realf SpatialCell::get_value(const vmesh::GlobalID& blockGID, const unsigned int cell, const int& popID) const {
      if (count(blockGID, popID) == 0) {
         return 0.0;
      }
      const vmesh::LocalID blockLID = get_velocity_block_local_id(blockGID, popID);
      //const Velocity_Block block_ptr = at(block);
      if (blockLID == invalid_local_id()) {
      //if (block_ptr.is_null() == true) {
         std::cerr << __FILE__ << ":" << __LINE__
            << " block_ptr == NULL" << std::endl; 
         abort();
      }
      return get_data(blockLID,popID)[cell];
   }

   inline bool SpatialCell::checkMesh(const int& popID) {
      #ifdef DEBUG_SPATIAL_CELL
      if (popID >= populations.size()) {
         std::cerr << "ERROR, popID " << popID << " exceeds populations.size() " << populations.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;             
         exit(1);
      }
      #endif
      
      return populations[popID].vmesh.check();
   }

   /*!
    Removes all velocity blocks from this spatial cell and frees memory in the cell
    */
    inline void SpatialCell::clear(const int& popID) {
       #ifdef DEBUG_SPATIAL_CELL
      if (popID >= populations.size()) {
         std::cerr << "ERROR, popID " << popID << " exceeds populations.size() " << populations.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;             
         exit(1);
      }
      #endif
       
      populations[popID].vmesh.clear();
      populations[popID].blockContainer.clear();
    }

   /*!
    Return the memory consumption in bytes as reported using the size()
    functions of the containers in spatial cell
    */
   inline uint64_t SpatialCell::get_cell_memory_size() {
      const uint64_t VEL_BLOCK_SIZE = 2*WID3*sizeof(Realf) + BlockParams::N_VELOCITY_BLOCK_PARAMS*sizeof(Real);
      uint64_t size = 0;
      size += vmeshTemp.sizeInBytes();
      size += blockContainerTemp.sizeInBytes();
      size += 2 * WID3 * sizeof(Realf);
      //size += mpi_velocity_block_list.size() * sizeof(vmesh::GlobalID);
      size += velocity_block_with_content_list.size() * sizeof(vmesh::GlobalID);
      size += velocity_block_with_no_content_list.size() * sizeof(vmesh::GlobalID);
      size += CellParams::N_SPATIAL_CELL_PARAMS * sizeof(Real);
      size += fieldsolver::N_SPATIAL_CELL_DERIVATIVES * sizeof(Real);
      size += bvolderivatives::N_BVOL_DERIVATIVES * sizeof(Real);

      for (size_t p=0; p<populations.size(); ++p) {
          size += populations[p].vmesh.sizeInBytes();
          size += populations[p].blockContainer.sizeInBytes();
      }

      return size;
   }

   /*!
    Return the memory consumption in bytes as reported using
    the size() functions of the containers in spatial cell
    */
   inline uint64_t SpatialCell::get_cell_memory_capacity() {
      const uint64_t VEL_BLOCK_SIZE = 2*WID3*sizeof(Realf) + BlockParams::N_VELOCITY_BLOCK_PARAMS*sizeof(Real);
      uint64_t capacity = 0;
      
      capacity += vmeshTemp.capacityInBytes();
      capacity += blockContainerTemp.capacityInBytes();
      capacity += 2 * WID3 * sizeof(Realf);
      //capacity += mpi_velocity_block_list.capacity()  * sizeof(vmesh::GlobalID);
      capacity += velocity_block_with_content_list.capacity()  * sizeof(vmesh::GlobalID);
      capacity += velocity_block_with_no_content_list.capacity()  * sizeof(vmesh::GlobalID);
      capacity += CellParams::N_SPATIAL_CELL_PARAMS * sizeof(Real);
      capacity += fieldsolver::N_SPATIAL_CELL_DERIVATIVES * sizeof(Real);
      capacity += bvolderivatives::N_BVOL_DERIVATIVES * sizeof(Real);
      
      for (size_t p=0; p<populations.size(); ++p) {
        capacity += populations[p].vmesh.capacityInBytes();
        capacity += populations[p].blockContainer.capacityInBytes();
      }
      
      return capacity;
   }
      
   /*!
    Adds an empty velocity block into this spatial cell.
    Returns true if given block was added or already exists.
    Returns false if given block is invalid or would be outside
    of the velocity grid.
    */
   inline bool SpatialCell::add_velocity_block(const vmesh::GlobalID& block,const int& popID) {
      #ifdef DEBUG_SPATIAL_CELL
      if (popID >= populations.size()) {
         std::cerr << "ERROR, popID " << popID << " exceeds populations.size() " << populations.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;             
         exit(1);
      }
      #endif
      
      // Block insert will fail, if the block already exists, or if 
      // there are too many blocks in the spatial cell
      bool success = true;
      if (populations[popID].vmesh.push_back(block) == false) {
         return false;
      }

      const vmesh::LocalID VBC_LID = populations[popID].blockContainer.push_back();

      // Set block data to zero values:
      Realf* data = populations[popID].blockContainer.getData(VBC_LID);
      for (unsigned int i=0; i<WID*WID*WID; ++i) data[i] = 0;

      // Set block parameters:
//      Real* parameters = get_block_parameters(populations[popID].vmesh.getLocalID(block));
      Real* parameters = get_block_parameters(VBC_LID,popID);
      parameters[BlockParams::VXCRD] = get_velocity_block_vx_min(popID,block);
      parameters[BlockParams::VYCRD] = get_velocity_block_vy_min(popID,block);
      parameters[BlockParams::VZCRD] = get_velocity_block_vz_min(popID,block);
      populations[popID].vmesh.getCellSize(block,&(parameters[BlockParams::DVX]));

      // The following call 'should' be the fastest, but is actually 
      // much slower that the parameter setting above
      //vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>::getBlockInfo(block,get_block_parameters( blockContainer.push_back() ));
      return success;
   }
   
   inline void SpatialCell::add_velocity_blocks(const std::vector<vmesh::GlobalID>& blocks,const int& popID) {
      #ifdef DEBUG_SPATIAL_CELL
      if (popID >= populations.size()) {
         std::cerr << "ERROR, popID " << popID << " exceeds populations.size() " << populations.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;             
         exit(1);
      }
      #endif
      
      // Add blocks to mesh
      const uint8_t adds = populations[popID].vmesh.push_back(blocks);
      if (adds == 0) {
         std::cerr << "Failed to add blocks" << std::endl;
         return;
      }

      // Add blocks to block container
      vmesh::LocalID startLID = populations[popID].blockContainer.push_back(blocks.size());
      Real* parameters = populations[popID].blockContainer.getParameters(startLID);

      #ifdef DEBUG_SPATIAL_CELL
         if (populations[popID].vmesh.size() != populations[popID].blockContainer.size()) {
	    std::cerr << "size mismatch in " << __FILE__ << ' ' << __LINE__ << std::endl; exit(1);
	 }
      #endif

      // Set block parameters
      for (size_t b=0; b<blocks.size(); ++b) {
         parameters[BlockParams::VXCRD] = get_velocity_block_vx_min(popID,blocks[b]);
         parameters[BlockParams::VYCRD] = get_velocity_block_vy_min(popID,blocks[b]);
         parameters[BlockParams::VZCRD] = get_velocity_block_vz_min(popID,blocks[b]);
         populations[popID].vmesh.getCellSize(blocks[b],&(parameters[BlockParams::DVX]));
         parameters += BlockParams::N_VELOCITY_BLOCK_PARAMS;
      }
   }

   inline bool SpatialCell::add_velocity_block_octant(const vmesh::GlobalID& blockGID,const int& popID) {
      #ifdef DEBUG_SPATIAL_CELL
      if (popID >= populations.size()) {
         std::cerr << "ERROR, popID " << popID << " exceeds populations.size() " << populations.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;             
         exit(1);
      }
      #endif
      
      // Return immediately if the block already exists or if an 
      // invalid block is attempted to be created:
      const vmesh::LocalID invalid = vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>::invalidLocalID();
      if (populations[popID].vmesh.getLocalID(blockGID) != invalid) return false;
      if (populations[popID].vmesh.count(blockGID) > 0) return false;

      /*
      // If parent exists, refine it:
      const vmesh::GlobalID parentGID = vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>::getParent(blockGID);
      if (parentGID != blockGID) {
	 if (vmesh.getLocalID(parentGID) != invalid) {
	    std::map<vmesh::GlobalID,vmesh::LocalID> insertedBlocks;
	    refine_block(parentGID,insertedBlocks);
	    return true;
	 }
      }*/

      // Attempt to add the block and its siblings. The vector
      // siblings also includes this block.
      // 
      // Note: These functions do not check for errors, such as creation of invalid blocks.
      // Creating all siblings at once is faster than creating one sibling at a time (16% difference)
      //std::vector<vmesh::GlobalID> siblings;
      //vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>::getSiblings(blockGID,siblings);
      //add_velocity_blocks(siblings);

//      for (size_t s=0; s<siblings.size(); ++s) {
//	 add_velocity_block(siblings[s]);
	 /*if (siblings[s] != vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>::invalidGlobalID()) {
	    // Check that the children of the sibling do not exist. 
	    // If they do, the block volume is already covered by 
	    // more refined blocks:
	    bool createBlock = true;
	    std::vector<vmesh::GlobalID> children;
	    vmesh.getChildren(siblings[s],children);
	    for (size_t c=0; c<children.size(); ++c) {
	       if (vmesh.getLocalID(children[c]) != invalid) {
		  createBlock = false;
		  break;
	       }
	    }
	    if (createBlock == true) add_velocity_block(siblings[s]);
	    add_velocity_block(siblings[s]);
	 }*/
//      }
    
      int32_t refLevel = populations[popID].vmesh.getRefinementLevel(blockGID);
      vmesh::GlobalID currentBlock = blockGID;
      do {
         std::vector<vmesh::GlobalID> siblings;
         populations[popID].vmesh.getSiblings(currentBlock,siblings);
         for (size_t s=0; s<siblings.size(); ++s) add_velocity_block(siblings[s],popID);

         currentBlock = populations[popID].vmesh.getParent(currentBlock);
         --refLevel;
      } while (refLevel > 0);

      return true;
   }

   /*!
    Removes given block from the velocity grid.
    Does nothing if given block doesn't exist.
    */
   inline void SpatialCell::remove_velocity_block(const vmesh::GlobalID& block,const int& popID) {
      #ifdef DEBUG_SPATIAL_CELL
      if (popID >= populations.size()) {
         std::cerr << "ERROR, popID " << popID << " exceeds populations.size() " << populations.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;             
         exit(1);
      }
      #endif
      
      if (block == invalid_global_id()) {
         //std::cerr << "not removing, block " << block << " is invalid" << std::endl;
         return;
      }
      
      const vmesh::LocalID removedLID = populations[popID].vmesh.getLocalID(block);      
      if (removedLID == invalid_local_id()) {
         //std::cerr << "not removing since block " << block << " does not exist" << std::endl;
         return;
      }

      // Get local ID of the last block:
      const vmesh::LocalID lastLID = populations[popID].vmesh.size()-1;

      populations[popID].vmesh.copy(lastLID,removedLID);
      populations[popID].vmesh.pop();

      populations[popID].blockContainer.copy(lastLID,removedLID);
      populations[popID].blockContainer.pop();
   }

   inline void SpatialCell::swap(vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& vmesh,
                                 vmesh::VelocityBlockContainer<vmesh::LocalID>& blockContainer,
                                 const int& popID) {
      #ifdef DEBUG_SPATIAL_CELL
      if (populations[popID].vmesh.size() != populations[popID].blockContainer.size()) {
         std::cerr << "Error, velocity mesh size and block container size do not agree in " << __FILE__ << ' ' << __LINE__ << std::endl;
         exit(1);
      }
      if (popID >= populations.size()) {
         std::cerr << "ERROR, popID " << popID << " exceeds populations.size() " << populations.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;             
         exit(1);
      }
      #endif

      populations[popID].vmesh.swap(vmesh);
      populations[popID].blockContainer.swap(blockContainer);
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
   
   inline bool SpatialCell::velocity_block_has_children(const vmesh::GlobalID& blockGID,const int& popID) const {
      #ifdef DEBUG_SPATIAL_CELL
      if (popID >= populations.size()) {
         std::cerr << "ERROR, popID " << popID << " exceeds populations.size() " << populations.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;             
         exit(1);
      }
      #endif
      return populations[popID].vmesh.hasChildren(blockGID);
   }

   inline vmesh::GlobalID SpatialCell::velocity_block_has_grandparent(const vmesh::GlobalID& blockGID,const int& popID) const {
      #ifdef DEBUG_SPATIAL_CELL
      if (popID >= populations.size()) {
         std::cerr << "ERROR, popID " << popID << " exceeds populations.size() " << populations.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;             
         exit(1);
      }
      #endif
      return populations[popID].vmesh.hasGrandParent(blockGID);
   }
   
   inline SpatialCell& SpatialCell::operator=(const SpatialCell&) { 
      return *this;
   }

} // namespaces

#endif
