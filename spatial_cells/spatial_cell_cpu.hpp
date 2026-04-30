/*
 * This file is part of Vlasiator.
 * Copyright 2010-2024 Finnish Meteorological Institute and University of Helsinki
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://www.physics.helsinki.fi/vlasiator/
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

/*!
Spatial cell class for Vlasiator that supports a variable number of velocity blocks.
*/

#ifndef VLASIATOR_SPATIAL_CELL_CPU_HPP
#define VLASIATOR_SPATIAL_CELL_CPU_HPP

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <mpi.h>
#include <limits>
#include <stdint.h>
#include <vector>
#include <array>
#include <unordered_map>
#include <set>
#include <map>
#include <phiprof.hpp>
#include <tuple>

#include "../memoryallocation.h"
#include "../common.h"
#include "../parameters.h"
#include "../definitions.h"

#include "velocity_mesh_cpu.h"
#include "velocity_block_container.h"

#ifdef DEBUG_VLASIATOR
   #ifndef DEBUG_SPATIAL_CELL
   #define DEBUG_SPATIAL_CELL
   #endif
#endif

namespace spatial_cell {

   /** Wrapper for variables needed for each particle species.
    *  Change order if you know what you are doing.
    * All Real fields should be consecutive, as they are communicated as a block.
    *
    */
   struct Population {
      Real RHO;
      Real V[3];
      Real RHO_R;
      Real V_R[3];
      Real RHO_V;
      Real V_V[3];
      Real P[6];
      Real P_R[6];
      Real P_V[6];
      Real RHO_R_PREV;
      Real V_R_PREV[3];
      Real P_R_PREV[3];
      Real RHOLOSSADJUST = 0.0;      /*!< Counter for particle number loss from the destroying blocks in blockadjustment*/
      Real max_dt[2];                                                /**< Element[0] is max_r_dt, element[1] max_v_dt.*/
      Real velocityBlockMinValue;

      // enum PropagationState = {UNINITIALIZED, TRANSLATED, ACCELERATED};
      Real T, T_R, T_V;

      uint ACCSUBCYCLES;          /*!< number of subcyles for each cell*/
      vmesh::LocalID N_blocks;    /**< Number of velocity blocks, used when receiving velocity
                                   * mesh from remote neighbors using MPI.*/
      vmesh::VelocityMesh *vmesh; /**< Velocity mesh. Contains all velocity blocks that exist
                                   * in this spatial cell. Cells are identified by their unique
                                   * global IDs.*/
      vmesh::VelocityBlockContainer *blockContainer;  /**< Velocity block data.*/

      /**< Temporary storage of acceleration transform intersections and sybcycling dt.*/
      Real intersection_z,intersection_z_di,intersection_z_dj,intersection_z_dk;
      Real intersection_x,intersection_x_di,intersection_x_dj,intersection_x_dk;
      Real intersection_y,intersection_y_di,intersection_y_dj,intersection_y_dk;
      Real subcycleDt;

      // Constructor, destructor
      Population() {
         vmesh = new vmesh::VelocityMesh();
         blockContainer = new vmesh::VelocityBlockContainer();
         // Set values to zero in case of zero-block populations
         RHO = RHO_R = RHO_V = RHOLOSSADJUST = velocityBlockMinValue = ACCSUBCYCLES = N_blocks = 0;
         for (uint i=0; i<2; ++i) {
            max_dt[i] = 0;
         }
         for (uint i=0; i<3; ++i) {
            V[i] = V_R[i] = V_V[i] = 0;
         }
         for (uint i=0; i<6; i++) {
            P[i] = P_R[i] = P_V[i] = 0;
         }
      }
      ~Population() {
         delete vmesh;
         delete blockContainer;
      }
      // Copy constructor
      Population(const Population& other) {
         vmesh = new vmesh::VelocityMesh(*(other.vmesh));
         blockContainer = new vmesh::VelocityBlockContainer(*(other.blockContainer));

         RHO = other.RHO;
         RHO_R = other.RHO_R;
         RHO_V = other.RHO_V;
         RHOLOSSADJUST = other.RHOLOSSADJUST;
         velocityBlockMinValue = other.velocityBlockMinValue;
         ACCSUBCYCLES = other.ACCSUBCYCLES;
         N_blocks = other.N_blocks;
         for (uint i=0; i<2; ++i) {
            max_dt[i] = other.max_dt[i];
         }
         for (uint i=0; i<3; ++i) {
            V[i] = other.V[i];
            V_R[i] = other.V_R[i];
            V_V[i] = other.V_V[i];
         }
         for (uint i=0; i<6; i++) {
            P[i] = other.P[i];
            P_R[i] = other.P_R[i];
            P_V[i] = other.P_V[i];
         }
      }

      Population(Population&& other) = delete; // Move constructor not implemented

      const Population& operator=(const Population& other) { // Copy assignment
         delete vmesh;
         delete blockContainer;
         vmesh = new vmesh::VelocityMesh(*(other.vmesh));
         blockContainer = new vmesh::VelocityBlockContainer(*(other.blockContainer));

         RHO = other.RHO;
         RHO_R = other.RHO_R;
         RHO_V = other.RHO_V;
         RHOLOSSADJUST = other.RHOLOSSADJUST;
         velocityBlockMinValue = other.velocityBlockMinValue;
         ACCSUBCYCLES = other.ACCSUBCYCLES;
         N_blocks = other.N_blocks;
         for (uint i=0; i<2; ++i) {
            max_dt[i] = other.max_dt[i];
         }
         for (uint i=0; i<3; ++i) {
            V[i] = other.V[i];
            V_R[i] = other.V_R[i];
            V_V[i] = other.V_V[i];
         }
         for (uint i=0; i<6; i++) {
            P[i] = other.P[i];
            P_R[i] = other.P_R[i];
            P_V[i] = other.P_V[i];
         }
         return *this;
      }

      Population& operator=(const Population&& other) = delete; // Move assignment not implemented


      void ResizeClear(const uint newSize) {
         // Resizes the vmesh localToGlobalMap, clears the vmesh GlobalToLocalMap,
         // and resizes the velocity block container.
         // Contents of the localToGlobalMap or the VBC are not edited.
         vmesh->setNewSize(newSize);
         vmesh->clearMap(newSize);
         blockContainer->setNewSize(newSize);
      }
      void Scale(creal factor) {
         RHO *= factor;
         RHO_R *= factor;
         RHO_V *= factor;
         for (uint i=0; i<3; ++i) {
            P[i] *= factor;
            P_R[i] *= factor;
            P_V[i] *= factor;
         }
         // Now loop over whole velocity space and scale the values
         for (vmesh::LocalID blockLID=0; blockLID < vmesh->size(); ++blockLID) {
            // Pointer to target block data
            Realf* toData = blockContainer->getData(blockLID);
            // Scale data
            for (uint i=0; i<WID3; ++i) {
               toData[i] = toData[i] * factor;
            }
         }
      }
      void Increment(const Population& other, creal factor) {
         // Note: moments will be invalidated.
         // Loop over the whole velocity space, and add scaled values.
         for (vmesh::LocalID incBlockLID=0; incBlockLID<(other.vmesh)->size(); ++incBlockLID) {
            const Realf* fromData = (other.blockContainer)->getData(incBlockLID);

            // Global ID of the block containing incoming data
            vmesh::GlobalID incBlockGID = (other.vmesh)->getGlobalID(incBlockLID);

            // Get local ID of the target block. If the block doesn't exist, create it.
            vmesh::GlobalID toBlockLID = vmesh->getLocalID(incBlockGID);
            if (toBlockLID == vmesh->invalidLocalID()) {
               bool success = vmesh->push_back(incBlockGID);
               toBlockLID = blockContainer->push_back_and_zero();
               Real* parameters = blockContainer->getParameters(toBlockLID);
               vmesh->getBlockInfo(incBlockGID, parameters+BlockParams::VXCRD);
            }

            // Pointer to target block data
            Realf* toData = blockContainer->getData(toBlockLID);
            // Add scaled values from source cells
            for (uint i=0; i<WID3; ++i) {
               toData[i] += fromData[i] * factor;
            }
         } // for-loop over velocity blocks
      }
   };

   typedef std::array<unsigned int, 3> velocity_cell_indices_t;   /**< Defines the indices of a velocity cell in a velocity block.
                                                                   * Indices start from 0 and the first value is the index in x direction.
                                                                   * Note: these are the (i,j,k) indices of the cell within the block.
                                                                   * Valid values are ([0,WID[,[0,WID[,[0,WID[).*/

   typedef std::array<vmesh::LocalID,3> velocity_block_indices_t; /**< Defines the indices of a velocity block in the velocity grid.
                                                                   * Indices start from 0 and the first value is the index in x direction.
                                                                   * Note: these are the (i,j,k) indices of the block.
                                                                   * Valid values are ([0,vx_length[,[0,vy_length[,[0,vz_length[).*/

   class SpatialCell {
   public:
      SpatialCell();
      ~SpatialCell();
      SpatialCell(const SpatialCell& other);
      const SpatialCell& operator=(const SpatialCell& other);

      // Null vmesh to check against null data
      static vmesh::VelocityMesh null_vmesh;
      // Following functions return velocity grid metadata // TODO timeclasses
      template<int PAD> void fetch_data(const vmesh::GlobalID& blockGID,const vmesh::VelocityMesh* vmesh,
                                        const Realf* src,Realf* array);
      template<int PAD>	void fetch_acc_data(const vmesh::GlobalID& blockGID,const int& dim,
					    vmesh::VelocityMesh* vmesh,
					    const Realf* src,Realf* array,Real cellSizeFractions[2]);

      vmesh::GlobalID find_velocity_block(vmesh::GlobalID cellIndices[3],const uint popID, const int timeclass);
      Realf* get_data(const uint popID, const int timeclass);
      const Realf* get_data(const uint popID, const int timeclass) const;
      Realf* get_data(const vmesh::LocalID& blockLID,const uint popID, const int timeclass);
      const Realf* get_data(const vmesh::LocalID& blockLID,const uint popID, const int timeclass) const;
      Real* get_block_parameters(const uint popID, const int timeclass);
      const Real* get_block_parameters(const uint popID, const int timeclass) const;
      Real* get_block_parameters(const vmesh::LocalID& blockLID,const uint popID, const int timeclass);
      const Real* get_block_parameters(const vmesh::LocalID& blockLID,const uint popID, const int timeclass) const;

      Real* get_cell_parameters();
      const Real* get_cell_parameters() const;

      vmesh::LocalID get_number_of_velocity_blocks(const uint popID, const int timeclass) const;
      vmesh::LocalID get_number_of_velocity_blocks_ghost(const uint popID) const;
      vmesh::LocalID get_number_of_all_velocity_blocks() const;
      int get_number_of_populations() const;
      void debug_population_check(const uint popID) const;
      void debug_population_check(const uint popID, const int timeclass) const;
      void debug_population_check(const uint popID, const vmesh::LocalID blockLID, const int timeclass) const;

      Population & get_population(const uint popID, const int timeclass);
      const Population & get_population(const uint popID, const int timeclass) const;
      void set_population(const Population& pop, cuint popID);
      void set_ghost_population(const Population& pop, cuint popID, const int timeclass);
      void scale_population(creal factor, cuint popID, const int timeclass);
      void increment_population(const Population& pop, creal factor, cuint popID, const int timeclass);

      std::vector<Population>& get_populations();
      const std::vector<Population>& get_populations() const;

      const Real& get_max_r_dt(const uint popID) const;
      const Real& get_max_v_dt(const uint popID) const;

      const Real& get_tc_dt() const;
      int get_tc() const;
      bool has_timeclass(int) const;
      bool get_timeclass_turn_r() const;
      bool get_timeclass_turn_v() const;
      bool get_timeclass_turn_v(int tc) const;

      const vmesh::LocalID* get_velocity_grid_length(const uint popID, const int timeclass);
      const Real* get_velocity_grid_block_size(const uint popID, const int timeclass);
      const Real* get_velocity_grid_cell_size(const uint popID, const int timeclass);
      void get_velocity_block_coordinates(const uint popID,const vmesh::GlobalID& globalID,Real* coords, const int timeclass);
      velocity_block_indices_t get_velocity_block_indices(const uint popID,const vmesh::GlobalID globalID, const int timeclass);
      vmesh::GlobalID get_velocity_block(const uint popID,vmesh::GlobalID blockIndices[3], const int timeclass) const;
      vmesh::GlobalID get_velocity_block(const uint popID,const velocity_block_indices_t indices, const int timeclass) const;
      vmesh::GlobalID get_velocity_block(const uint popID,const Real* coords, const int timeclass) const;
      vmesh::GlobalID get_velocity_block(const uint popID,const Real vx,const Real vy,const Real vz, const int timeclass) const;
      vmesh::GlobalID get_velocity_block_child(const uint popID,const vmesh::GlobalID& blockGID,
                                               const int& i_cell,const int& j_cell,const int& k_cell, const int timeclass =-1);
      void get_velocity_block_children_local_ids(const vmesh::GlobalID& blockGID,
                                                 std::vector<vmesh::LocalID>& childrenLIDs,
                                                 const uint popID, const int timeclass = -1);
      vmesh::GlobalID get_velocity_block_parent(const uint popID,const vmesh::GlobalID& blockGID, const int timeclass = -1); 
      vmesh::GlobalID get_velocity_block_global_id(const vmesh::LocalID& blockLID,const uint popID, const int timeclass) const;
      vmesh::LocalID get_velocity_block_local_id(const vmesh::GlobalID& blockGID,const uint popID, const int timeclass) const;
      void get_velocity_block_size(const uint popID,const vmesh::GlobalID block,Real size[3], const int timeclass);
      Real get_velocity_block_vx_min(const uint popID,const vmesh::GlobalID block, const int timeclass) const;
      Real get_velocity_block_vx_max(const uint popID,const vmesh::GlobalID block, const int timeclass) const;
      Real get_velocity_block_vy_min(const uint popID,const vmesh::GlobalID block, const int timeclass) const;
      Real get_velocity_block_vy_max(const uint popID,const vmesh::GlobalID block, const int timeclass) const;
      Real get_velocity_block_vz_min(const uint popID,const vmesh::GlobalID block, const int timeclass) const;
      Real get_velocity_block_vz_max(const uint popID,const vmesh::GlobalID block, const int timeclass) const;

      static unsigned int invalid_block_index();
      static vmesh::GlobalID invalid_global_id();
      static vmesh::LocalID invalid_local_id();

      size_t count(const vmesh::GlobalID& block,const uint popID, const int timeclass) const;

      void add_values(const vmesh::GlobalID& targetGID,
		      std::unordered_map<vmesh::GlobalID,Realf[(WID+2)*(WID+2)*(WID+2)]>& sourceData,
                      const uint popID);

      void printMeshSizes();
      static bool setCommunicatedSpecies(const uint popID, const int timeclass = -1);

      // Following functions adjust velocity blocks stored on the cell //
      bool add_velocity_block(const vmesh::GlobalID& block,const uint popID, const int timeclass);
      void adjustSingleCellVelocityBlocks(const uint popID, bool doDeleteEmpty=false, const int timeclass=-1);
      void adjust_velocity_blocks(const std::vector<SpatialCell*>& spatial_neighbors,
                                  const uint popID,
                                  bool doDeleteEmptyBlocks=true, const int timeclass=-1);
      // Templated function for storing a v-space read from a file or generated elsewhere
      template <typename fileReal> void add_velocity_blocks(const uint popID,const std::vector<vmesh::GlobalID>& blocks,fileReal* avgBuffer, const int timeclass=-1);

      void update_velocity_block_content_lists(const uint popID, const int timeclass=-1);
      bool checkMesh(const uint popID, const int timeclass);
      void clear(const uint popID, bool shrink=false, const int timeclass=-1);
      void setNewSizeClear(const uint popID, const vmesh::LocalID& newSize, const int timeclass=-1);
      void setNewSizeClear(const uint popID, const int timeclass=-1);

      uint64_t get_cell_memory_capacity();
      uint64_t get_cell_memory_size();
      void merge_values(const uint popID, const int timeclass=-1);
      void prepare_to_receive_blocks(const uint popID, const int timeclass=-1);
      bool shrink_to_fit();
      size_t size(const uint popID, const int timeclass) const;
      void remove_velocity_block(const vmesh::GlobalID& block,const uint popID, const int timeclass);      
      vmesh::VelocityMesh* get_velocity_mesh(const size_t& popID, const int timeclass);
      const vmesh::VelocityMesh* get_velocity_mesh(const size_t& popID, const int timeclass) const;
      vmesh::VelocityBlockContainer* get_velocity_blocks(const size_t& popID, const int timeclass = -1);
      const vmesh::VelocityBlockContainer* get_velocity_blocks(const size_t& popID, const int timeclass = -1) const;

      void set_velocity_mesh_ghost(const size_t& popID, const int timeclass, const int src_timeclass);
      void set_velocity_blocks_ghost(const size_t& popID, const int timeclass, const int src_timeclass);
      vmesh::VelocityMesh* get_velocity_mesh_ghost(const size_t& popID, const int timeclass);
      vmesh::VelocityBlockContainer* get_velocity_blocks_ghost(const size_t& popID, const int timeclass);

      void set_max_r_dt(const uint popID,const Real& value);
      void set_max_v_dt(const uint popID,const Real& value);

      // Following functions are related to MPI //
      std::tuple<void*, int, MPI_Datatype> get_mpi_datatype(const CellID cellID,const int sender_rank,const int receiver_rank,
                                                            const bool receiving,const int neighborhood);
      static uint64_t get_mpi_transfer_type(void);
      static void set_mpi_transfer_type(const uint64_t type,bool atSysBoundaries=false);
      static void set_mpi_transfer_direction(const int dimension);
      void set_mpi_transfer_enabled(bool transferEnabled);
      void updateSparseMinValue(const uint popID, const int timeclass =-1);
      Real getVelocityBlockMinValue(const uint popID, const int timeclass = -1) const;

      // Member variables //
      std::array<Real, vderivatives::N_V_DERIVATIVES> derivativesV;           /**< Derivatives of V for vorticity AMR.*/
      std::array<Real, bvolderivatives::N_BVOL_DERIVATIVES> derivativesBVOL;  /**< Derivatives of BVOL needed by the acceleration.
                                                                               * Separate array because it does not need to be communicated.*/
      std::array<Real, CellParams::N_SPATIAL_CELL_PARAMS> parameters;         /**< Bulk variables in this spatial cell.*/
      std::array<Realf, WID3> null_block_data;

      CellID get_cellid(){return parameters[CellParams::CELLID];};
      uint64_t ioLocalCellId;                                                 /**< Local cell ID used for IO, not needed elsewhere
                                                                               * and thus not being kept up-to-date.*/
      std::array<Realf*,MAX_NEIGHBORS_PER_DIM> neighbor_block_data;           /**< Pointers for translation operator. We can point to neighbor
                                                                               * cell block data. We do not allocate memory for the pointer.*/
      std::array<vmesh::LocalID,MAX_NEIGHBORS_PER_DIM> neighbor_number_of_blocks;
      std::map<int,std::set<int>> face_neighbor_ranks;
      uint sysBoundaryFlag;                                                   /**< What type of system boundary does the cell belong to.
                                                                               * Enumerated in the sysboundarytype namespace's enum.*/
      uint sysBoundaryLayer;                                                  /**< Layers counted from closest systemBoundary. If 0 then it has not
                                                                               * been computed. First sysboundary layer is layer 1.*/
      int sysBoundaryLayerNew;                                                /** needed (by DCCRG?), do not remove. */
      std::vector<vmesh::GlobalID> *velocity_block_with_content_list;         /**< List of existing cells with content, only up-to-date after
                                                                               * call to update_has_content().*/
      vmesh::LocalID velocity_block_with_content_list_size;                   /**< Size of vector. Needed for MPI communication of size before actual list transfer.*/
      std::vector<vmesh::GlobalID> *velocity_block_with_no_content_list;      /**< List of existing cells with no content, only up-to-date after
                                                                               * call to update_has_content. This is also never transferred
                                                                               * over MPI, so is invalid on remote cells.*/
      static uint64_t mpi_transfer_type;                                      /**< Which data is transferred by the mpi datatype given by spatial cells.*/
      static bool mpiTransferAtSysBoundaries;                                 /**< Do we only transfer data at boundaries (true), or in the whole system (false).*/

      std::set<int> requested_timeclass_ghosts = {};                       /**< See Pencil construction. Translation stencil neighbours may want v-space values at
      *   varying timeclass synchronizations. This keeps track which levels are requested of this
      *   cell. Populations struct contains mappings of these timeclasses to ghost vmeshes. */
      std::set<int> requested_timeclass_copy_ghosts = {};                       /**< See Pencil construction. Translation stencil neighbours may want v-space values at
      *   varying timeclass synchronizations. This keeps track which levels are requested of this
      *   cell. Populations struct contains mappings of these timeclasses to ghost vmeshes. */

      inline std::set<int> get_all_ghosts();
      //SpatialCell& operator=(const SpatialCell& other);
    private:
      //SpatialCell& operator=(const SpatialCell&);

      bool compute_block_has_content(const vmesh::GlobalID& block,const uint popID, const int timeclass=-1) const;

      static int activePopID;
      static int activeTimeclass;
      bool initialized;
      bool mpiTransferEnabled;

      std::vector<spatial_cell::Population> populations;                        /**< Particle population variables.*/
      std::map<std::pair<const uint, const int>, spatial_cell::Population> ghostPopulations; // Key is {popID, timeclass}
   };

   inline void SpatialCell::debug_population_check(const uint popID) const {
      #ifdef DEBUG_SPATIAL_CELL
      if (popID >= populations.size()) {
         std::cerr << "ERROR, popID " << popID << " exceeds populations.size() " << populations.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
         exit(1);
      }
      #endif
   }

   inline void SpatialCell::debug_population_check(const uint popID, const int timeclass) const {
      debug_population_check(popID);
      #ifdef DEBUG_SPATIAL_CELL
      if (this->ghostPopulations.find({popID, timeclass}) == ghostPopulations.end()) {
         std::cerr << "ERROR, popID " << popID << " with timeclass " << timeclass << " not found in ghostPopulations in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
         exit(1);
      }
      #endif
   }

   inline void SpatialCell::debug_population_check(const uint popID, const vmesh::LocalID blockLID, const int timeclass = -1) const {
      debug_population_check(popID, timeclass);
      #ifdef DEBUG_SPATIAL_CELL
      if (blockLID >= this->get_velocity_blocks(popID, timeclass)->size()) {
         std::cerr << "ERROR, block LID out of bounds, blockContainer->size() " << populations[popID].blockContainer->size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
         exit(1);
      }
      #endif
   }

   inline vmesh::GlobalID SpatialCell::find_velocity_block(vmesh::GlobalID cellIndices[3],const uint popID, const int timeclass = -1) {
      debug_population_check(popID, timeclass);
      return this->get_velocity_mesh(popID, timeclass)->findBlock(cellIndices);
   }

   inline Realf* SpatialCell::get_data(const uint popID, const int timeclass = -1) {
      debug_population_check(popID, timeclass);
      return this->get_velocity_blocks(popID, timeclass)->getData();
   }

   inline const Realf* SpatialCell::get_data(const uint popID, const int timeclass = -1) const {
      debug_population_check(popID, timeclass);
      return this->get_velocity_blocks(popID, timeclass)->getData();
   }

   inline Realf* SpatialCell::get_data(const vmesh::LocalID& blockLID,const uint popID, const int timeclass = -1) {
      debug_population_check(popID, blockLID, timeclass);
      if (blockLID == vmesh::VelocityMesh::invalidLocalID()) {
         return null_block_data.data();
      }
      return this->get_velocity_blocks(popID, timeclass)->getData(blockLID);
   }

   inline const Realf* SpatialCell::get_data(const vmesh::LocalID& blockLID,const uint popID, const int timeclass = -1) const {
      debug_population_check(popID,blockLID, timeclass);
      if (blockLID == vmesh::VelocityMesh::invalidLocalID()) {
         return null_block_data.data();
      }
      return this->get_velocity_blocks(popID, timeclass)->getData(blockLID);
   }

   inline Real* SpatialCell::get_block_parameters(const uint popID, const int timeclass = -1) {
      debug_population_check(popID,timeclass);
      return this->get_velocity_blocks(popID, timeclass)->getParameters();
   }

   inline const Real* SpatialCell::get_block_parameters(const uint popID, const int timeclass = -1) const {
      debug_population_check(popID,timeclass);
      return this->get_velocity_blocks(popID, timeclass)->getParameters();
   }

   inline Real* SpatialCell::get_block_parameters(const vmesh::LocalID& blockLID, const uint popID, const int timeclass = -1) {
      debug_population_check(popID,blockLID, timeclass);
      return this->get_velocity_blocks(popID, timeclass)->getParameters(blockLID);
   }

   inline const Real* SpatialCell::get_block_parameters(const vmesh::LocalID& blockLID, const uint popID, const int timeclass = -1) const {
      debug_population_check(popID,blockLID, timeclass);
      return this->get_velocity_blocks(popID, timeclass)->getParameters(blockLID);
   }

   inline Real* SpatialCell::get_cell_parameters() {
      return parameters.data();
   }

   inline const Real* SpatialCell::get_cell_parameters() const {
      return parameters.data();
   }

   inline vmesh::LocalID SpatialCell::get_number_of_velocity_blocks(const uint popID, const int timeclass=-1) const {
      debug_population_check(popID, timeclass);
      return this->get_velocity_blocks(popID, timeclass)->size();
   }

   /** Get the total number of velocity blocks in the time ghosts of this cell, for a given population.
     * @return Number of blocks in ghost populations for popID*/
   inline vmesh::LocalID SpatialCell::get_number_of_velocity_blocks_ghost(const uint popID) const {
      debug_population_check(popID);
      vmesh::LocalID sum = 0;
      for (int tc : this->requested_timeclass_ghosts) {
         sum += ghostPopulations.at({popID,tc}).blockContainer->size();
      }
      for (int tc : this->requested_timeclass_copy_ghosts) {
         sum += ghostPopulations.at({popID,tc}).blockContainer->size();
      }
      return sum;
   }

    /** Get the total number of velocity blocks in this cell, summed over
     * all existing particle populations.
     * @return Total number of velocity blocks in the cell.*/
   inline vmesh::LocalID SpatialCell::get_number_of_all_velocity_blocks() const {
      vmesh::LocalID N_blocks = 0;
      for (size_t p=0; p<populations.size(); ++p)
         N_blocks += get_velocity_blocks(p)->size();
      for(auto& pop:this->ghostPopulations){
         N_blocks += pop.second.blockContainer->size();
      }
      return N_blocks;
   }

   inline int SpatialCell::get_number_of_populations() const {
      return populations.size();
   }

   inline Population & SpatialCell::get_population(const uint popID, const int timeclass=-1){
      if (timeclass < 0 || this->parameters[CellParams::TIMECLASS] == timeclass){
         return populations[popID];
      }
      else{
         return ghostPopulations.at({popID, timeclass});
      }
      
   }
   
   inline const Population & SpatialCell::get_population(const uint popID, const int timeclass=-1) const {
      if (timeclass < 0 || this->parameters[CellParams::TIMECLASS] == timeclass){
         return populations[popID];
      }
      else{
         return ghostPopulations.at({popID, timeclass}); // try-emplace?
      }
   }

   inline void SpatialCell::set_population(const Population& pop, cuint popID) {
      this->populations[popID] = pop;
   }

   inline void SpatialCell::set_ghost_population(const Population& pop, cuint popID, const int timeclass = -1) {
      ghostPopulations.erase({popID, timeclass});
      ghostPopulations.emplace(std::pair<cuint, const int>({popID, timeclass}), pop);
   }

   
   inline void SpatialCell::scale_population(creal factor, cuint popID, const int timeclass=-1) {
      (this->get_population(popID, timeclass)).Scale(factor);
   }
   inline void SpatialCell::increment_population(const Population& pop, creal factor, cuint popID, const int timeclass=-1) {
      (this->get_population(popID, timeclass)).Increment(pop, factor);
   }
   inline std::vector<Population>& SpatialCell::get_populations() {
      return populations;
   }
   inline const std::vector<Population>& SpatialCell::get_populations() const {
      return populations;
   }

   inline const vmesh::LocalID* SpatialCell::get_velocity_grid_length(const uint popID, const int timeclass=-1) {
      return this->get_velocity_mesh(popID, timeclass)->getGridLength();
   }

   inline const Real* SpatialCell::get_velocity_grid_block_size(const uint popID, const int timeclass=-1) {
      return this->get_velocity_mesh(popID, timeclass)->getBlockSize();
   }

   inline const Real* SpatialCell::get_velocity_grid_cell_size(const uint popID, const int timeclass=-1) {
      return this->get_velocity_mesh(popID, timeclass)->getCellSize();
   }

   inline void SpatialCell::get_velocity_block_coordinates(const uint popID,const vmesh::GlobalID& globalID,Real* coords, const int timeclass=-1) {
      this->get_velocity_mesh(popID, timeclass)->getBlockCoordinates(globalID,coords);
   }

   /*!
    Returns the indices of given velocity block
    */
   inline velocity_block_indices_t SpatialCell::get_velocity_block_indices(const uint popID,const vmesh::GlobalID block, const int timeclass=-1) {
      velocity_block_indices_t indices;
      this->get_velocity_mesh(popID, timeclass)->getIndices(block,indices[0],indices[1],indices[2]);
      return indices;
   }

   /*!
    Returns the velocity block at given indices or error_velocity_block
    */
   inline vmesh::GlobalID SpatialCell::get_velocity_block(const uint popID,const velocity_block_indices_t indices, const int timeclass=-1) const {
      return this->get_velocity_mesh(popID, timeclass)->getGlobalID(indices[0],indices[1],indices[2]);
   }

   inline vmesh::GlobalID SpatialCell::get_velocity_block(const uint popID,vmesh::GlobalID blockIndices[3], const int timeclass=-1) const {
      return this->get_velocity_mesh(popID, timeclass)->getGlobalID(blockIndices[0],blockIndices[1],blockIndices[2]);
   }

   /*!
    Returns the velocity block at given location or
    error_velocity_block if outside of the velocity grid
    */
   inline vmesh::GlobalID SpatialCell::get_velocity_block(const uint popID,const Real vx,const Real vy,const Real vz, const int timeclass=-1) const {
      Real coords[3] = {vx,vy,vz};
      return this->get_velocity_mesh(popID, timeclass)->getGlobalID(coords);
   }

   inline vmesh::GlobalID SpatialCell::get_velocity_block(const uint popID,const Real* coords, const int timeclass=-1) const {
      return this->get_velocity_mesh(popID, timeclass)->getGlobalID(coords);
   }

   inline vmesh::GlobalID SpatialCell::get_velocity_block_global_id(const vmesh::LocalID& blockLID,const uint popID, const int timeclass=-1) const {
      debug_population_check(popID);
      return this->get_velocity_mesh(popID, timeclass)->getGlobalID(blockLID);
   }

   inline vmesh::LocalID SpatialCell::get_velocity_block_local_id(const vmesh::GlobalID& blockGID,const uint popID, const int timeclass = -1) const {
      debug_population_check(popID);
      return this->get_velocity_mesh(popID, timeclass)->getLocalID(blockGID);
   }

   inline void SpatialCell::get_velocity_block_size(const uint popID,const vmesh::GlobalID block,Real blockSize[3], const int timeclass = -1) {
      this->get_velocity_mesh(popID, timeclass)->getBlockSize(block,blockSize);
   }

   /*!
    Returns the edge where given velocity block starts.
    */
   inline Real SpatialCell::get_velocity_block_vx_min(const uint popID,const vmesh::GlobalID block, const int timeclass = -1) const {
      Real coords[3];
      this->get_velocity_mesh(popID, timeclass)->getBlockCoordinates(block,coords);
      return coords[0];
   }

   /*!
    Returns the edge where given velocity block ends.
    */
   inline Real SpatialCell::get_velocity_block_vx_max(const uint popID,const vmesh::GlobalID block, const int timeclass = -1) const {
      Real coords[3];
      this->get_velocity_mesh(popID, timeclass)->getBlockCoordinates(block,coords);

      Real size[3];
      this->get_velocity_mesh(popID, timeclass)->getBlockSize(block,size);
      return coords[0]+size[0];
   }

   /*!
    Returns the edge where given velocity block starts.
    */
   inline Real SpatialCell::get_velocity_block_vy_min(const uint popID,const vmesh::GlobalID block, const int timeclass = -1) const {
      Real coords[3];
      this->get_velocity_mesh(popID, timeclass)->getBlockCoordinates(block,coords);
      return coords[1];
   }

   /*!
    Returns the edge where given velocity block ends.
    */
   inline Real SpatialCell::get_velocity_block_vy_max(const uint popID,const vmesh::GlobalID block, const int timeclass = -1) const {
      Real coords[3];
      this->get_velocity_mesh(popID, timeclass)->getBlockCoordinates(block,coords);

      Real size[3];
      this->get_velocity_mesh(popID, timeclass)->getBlockSize(block,size);
      return coords[1]+size[1];
   }

   /*!
    Returns the edge where given velocity block starts.
    */
   inline Real SpatialCell::get_velocity_block_vz_min(const uint popID,const vmesh::GlobalID block, const int timeclass = -1) const {
      Real coords[3];
      this->get_velocity_mesh(popID, timeclass)->getBlockCoordinates(block,coords);
      return coords[2];
   }

   /*!
    Returns the edge where given velocity block ends.
    */
   inline Real SpatialCell::get_velocity_block_vz_max(const uint popID,const vmesh::GlobalID block, const int timeclass = -1) const {
      Real coords[3];
      this->get_velocity_mesh(popID, timeclass)->getBlockCoordinates(block,coords);

      Real size[3];
      this->get_velocity_mesh(popID, timeclass)->getBlockSize(block,size);
      return coords[2]+size[2];
   }

   inline unsigned int SpatialCell::invalid_block_index() {
      return vmesh::VelocityMesh::invalidBlockIndex();
   }

   inline vmesh::GlobalID SpatialCell::invalid_global_id() {
      return vmesh::VelocityMesh::invalidGlobalID();
   }

   inline vmesh::GlobalID SpatialCell::invalid_local_id() {
      return vmesh::VelocityMesh::invalidLocalID();
   }

   /*!
    Returns the number of given velocity blocks that exist.
    */
   inline size_t SpatialCell::count(const vmesh::GlobalID& block,const uint popID, const int timeclass = -1) const {
      debug_population_check(popID);
      return this->get_velocity_mesh(popID,timeclass)->count(block);
   }

   /*!
    Returns the number of existing velocity blocks.
    */
   inline size_t SpatialCell::size(const uint popID, const int timeclass = -1) const {
      debug_population_check(popID);
      return get_velocity_mesh(popID,timeclass)->size();
   }

   inline vmesh::VelocityMesh* SpatialCell::get_velocity_mesh(const size_t& popID, const int timeclass = -1) {
      debug_population_check(popID);
      if (timeclass < 0 || this->parameters[CellParams::TIMECLASS] == timeclass) {
         // std::cout << "regular vmesh\n";
         return this->populations[popID].vmesh;
      } else {
         // std::cout << "ghost vmesh\n";
         return this->ghostPopulations.at({popID,timeclass}).vmesh;
      }
   }

   inline const vmesh::VelocityMesh* SpatialCell::get_velocity_mesh(const size_t& popID, const int timeclass = -1) const {
      debug_population_check(popID);
      if (timeclass < 0 || this->parameters[CellParams::TIMECLASS] == timeclass) {
         return this->populations[popID].vmesh;
      } else {
         return this->ghostPopulations.at({popID,timeclass}).vmesh;
      }
   }   

   inline vmesh::VelocityBlockContainer* SpatialCell::get_velocity_blocks(const size_t& popID, const int timeclass) {
      debug_population_check(popID);
      if (timeclass < 0 || this->parameters[CellParams::TIMECLASS] == timeclass) {
         return this->populations[popID].blockContainer;
      } else {
         return this->ghostPopulations.at({popID,timeclass}).blockContainer;
      }   
   }

   inline const vmesh::VelocityBlockContainer* SpatialCell::get_velocity_blocks(const size_t& popID, const int timeclass) const {
      debug_population_check(popID);
      if (timeclass < 0 || this->parameters[CellParams::TIMECLASS] == timeclass) {
         return this->populations[popID].blockContainer;
      } else {
         return this->ghostPopulations.at({popID,timeclass}).blockContainer;
      }   
   }

   inline void SpatialCell::set_velocity_mesh_ghost(const size_t& popID, const int timeclass, const int src_timeclass = -1) {
      debug_population_check(popID);
      assert(timeclass != src_timeclass);
      delete this->ghostPopulations[{popID,timeclass}].vmesh;
      // std::cout << __FILE__ <<":" <<__LINE__ << "\n";
      vmesh::VelocityMesh* newVmesh = new vmesh::VelocityMesh(*this->get_velocity_mesh(popID, src_timeclass));
      // std::cout <<"A new vmesh appeared, with meshID " << newVmesh->getMesh() << " at "<< this->get_cellid() << "\n";
      this->ghostPopulations[{popID,timeclass}].vmesh = newVmesh;
      // this->ghostPopulations[{popID,timeclass}].vmesh = vmesh::VelocityMesh(*this->populations[popID].vmesh);
      int vmsum = 0;
      for(auto it : *(this->ghostPopulations[{popID,timeclass}].vmesh)->getGrid()){
         vmsum += it;
      }

      #ifdef DEBUG_SPATIAL_CELL
      std::cout << "Copy-constructed ghostPopulations[{"<<popID<<","<<timeclass<<"}].vmesh (idsum " << vmsum <<"); " << &this->ghostPopulations[{popID,timeclass}].vmesh<< " from initial at " << &this->populations[popID].vmesh <<"\n";
      #endif
   }

   inline void SpatialCell::set_velocity_blocks_ghost(const size_t& popID, const int timeclass, const int src_timeclass = -1) {
      debug_population_check(popID);
      assert(timeclass != src_timeclass);
      delete this->ghostPopulations[{popID,timeclass}].blockContainer;
      vmesh::VelocityBlockContainer* newBlockContainer = new vmesh::VelocityBlockContainer(*this->get_velocity_blocks(popID, src_timeclass));
      // this->ghostPopulations[{popID,timeclass}].blockContainer = this->get_velocity_blocks(popID, src_timeclass);
      this->ghostPopulations[{popID,timeclass}].blockContainer = newBlockContainer;
      #ifdef DEBUG_SPATIAL_CELL
      int blocksum = 0;
      for (vmesh::LocalID block_i=0; block_i<this->ghostPopulations[{popID,timeclass}].blockContainer->size(); ++block_i) {
         const vmesh::GlobalID block = (this->ghostPopulations[{popID,timeclass}].vmesh)->getGlobalID(block_i);
         blocksum += block;
      }

      std::cout << "Copy-constructed ghostPopulations[{"<<popID<<","<<timeclass<<"}].blockContainer to " << &this->ghostPopulations[{popID,timeclass}].blockContainer<< " (blockidsum " << blocksum << ")  from initial at " << &this->populations[popID].blockContainer << ", new size " << this->ghostPopulations[{popID,timeclass}].blockContainer->size() << "; old size:" << this->ghostPopulations[{popID,timeclass}].blockContainer->size() << "\n";
      #endif
   }

   inline vmesh::VelocityMesh* SpatialCell::get_velocity_mesh_ghost(const size_t& popID, const int timeclass) {
      debug_population_check(popID);
      // std::cerr << "get_velocity_mesh_ghost with tc " << timeclass << " at " << &ghostPopulations[{popID,timeclass}].vmesh <<"\n";

      return this->ghostPopulations.at({popID,timeclass}).vmesh; // OBS try-emplace
   }

   inline vmesh::VelocityBlockContainer* SpatialCell::get_velocity_blocks_ghost(const size_t& popID, const int timeclass) {
      debug_population_check(popID);
      
      return this->ghostPopulations[{popID,timeclass}].blockContainer; // OBS try-emplace
   }

   inline bool SpatialCell::checkMesh(const uint popID, const int timeclass = -1) {
      debug_population_check(popID);
      return this->get_velocity_mesh(popID, timeclass)->check();
   }

   /*!
    * Removes all velocity blocks from this spatial cell and frees memory in the cell
   */
   inline void SpatialCell::clear(const uint popID, bool shrink, const int timeclass) {
      debug_population_check(popID);
      this->get_velocity_mesh(popID, timeclass)->clear(shrink);
      this->get_velocity_blocks(popID, timeclass)->clear(shrink);
    }

   /*!
     Ensures the selected population vmesh localToGlobalMap has sufficient capacity and is of correct size. Does not alter
     its contents. Ensures the selected population vmesh globalToLocalMap has sufficient capacity and is empty.
     Ensures the selected population VBC has sufficient capacity and is of correct size.
   */
  //!TODO timeclasses
   inline void SpatialCell::setNewSizeClear(const uint popID, const vmesh::LocalID& newSize, const int timeclass) {
      get_population(popID, timeclass).ResizeClear(newSize);
   }
   inline void SpatialCell::setNewSizeClear(const uint popID, const int timeclass) {
      get_population(popID, timeclass).ResizeClear(get_population(popID, timeclass).N_blocks);
   }

   /*!
    Return the memory consumption in bytes as reported using the size()
    functions of the containers in spatial cell
    */
   // TODO timeclasses
   inline uint64_t SpatialCell::get_cell_memory_size() {
      uint64_t size = 0;
      size += 2 * WID3 * sizeof(Realf);
      size += velocity_block_with_content_list->size() * sizeof(vmesh::GlobalID);
      size += velocity_block_with_no_content_list->size() * sizeof(vmesh::GlobalID);
      size += CellParams::N_SPATIAL_CELL_PARAMS * sizeof(Real);
      size += bvolderivatives::N_BVOL_DERIVATIVES * sizeof(Real);

      for (size_t p=0; p<populations.size(); ++p) {
          size += get_velocity_mesh(p)->sizeInBytes();
          size += get_velocity_blocks(p)->sizeInBytes();
      }

      for (size_t p=0; p<populations.size(); ++p) {
         for (int tc : this->requested_timeclass_ghosts) {
            size += get_velocity_mesh(p, tc)->sizeInBytes();
            size += get_velocity_blocks(p, tc)->sizeInBytes();
         }
      }

      return size;
   }

   /*!
    Return the memory consumption in bytes as reported using
    the size() functions of the containers in spatial cell
    */
   // TODO timeclasses
   inline uint64_t SpatialCell::get_cell_memory_capacity() {
      uint64_t capacity = 0;

      capacity += 2 * WID3 * sizeof(Realf);
      capacity += velocity_block_with_content_list->capacity()  * sizeof(vmesh::GlobalID);
      capacity += velocity_block_with_no_content_list->capacity()  * sizeof(vmesh::GlobalID);
      capacity += CellParams::N_SPATIAL_CELL_PARAMS * sizeof(Real);
      capacity += bvolderivatives::N_BVOL_DERIVATIVES * sizeof(Real);

      for (size_t p=0; p<populations.size(); ++p) {
        capacity += get_velocity_mesh(p)->capacityInBytes();
        capacity += get_velocity_blocks(p)->capacityInBytes();
      }

      for (size_t p=0; p<populations.size(); ++p) {
         for (int tc : this->requested_timeclass_ghosts) {
            capacity += get_velocity_mesh(p, tc)->capacityInBytes();
            capacity += get_velocity_blocks(p, tc)->capacityInBytes();
         }
      }

      return capacity;
   }

   /*!
    Adds an empty velocity block into this spatial cell.
    Returns true if given block was added or already exists.
    Returns false if given block is invalid or would be outside
    of the velocity grid.
    */
   inline bool SpatialCell::add_velocity_block(const vmesh::GlobalID& block,const uint popID, const int timeclass = -1) {
      debug_population_check(popID);
      // Block insert will fail, if the block already exists, or if
      // there are too many blocks in the spatial cell
      bool success = true;
      if (this->get_velocity_mesh(popID, timeclass)->push_back(block) == false) {
         return false;
      }

      const vmesh::LocalID VBC_LID = this->get_velocity_blocks(popID, timeclass)->push_back_and_zero();

      // Set block parameters:
      Real* parameters = get_block_parameters(VBC_LID, popID, timeclass);
      this->get_velocity_mesh(popID, timeclass)->getBlockInfo(block, parameters);

      // The following call 'should' be the fastest, but is actually
      // much slower that the parameter setting above
      //vmesh::VelocityMesh::getBlockInfo(block,get_block_parameters( blockContainer->push_back() ));
      return success;
   }

   /** Adds a vector of velocity blocks to the population, sets the parameters, and fills the data
       with phase-space densities from the provided buffer (which was read from a file).
   */
   template <typename fileReal> void SpatialCell::add_velocity_blocks(const uint popID, const std::vector<vmesh::GlobalID>& blocks,fileReal* avgBuffer, const int timeclass) {
      debug_population_check(popID);
      // Return if no blocks to add
      const uint nBlocks = blocks.size();
      if (nBlocks==0) {
         return;
      }
      this->get_velocity_mesh(popID, timeclass)->setNewCapacity(nBlocks);
      this->get_velocity_blocks(popID, timeclass)->setNewCapacity(nBlocks);

      // Add blocks to mesh
      const uint8_t adds = this->get_velocity_mesh(popID, timeclass)->push_back(blocks);
      if (adds == 0) {
         std::cerr << "Failed to add blocks" << std::endl;
         return;
      }

      // Add blocks to block container
      vmesh::LocalID startLID = this->get_velocity_blocks(popID, timeclass)->push_back(nBlocks);
      Real* parameters = this->get_velocity_blocks(popID, timeclass)->getParameters(startLID);

      #ifdef DEBUG_SPATIAL_CELL
      if (this->get_velocity_mesh(popID, timeclass)->size() != this->get_velocity_blocks(popID, timeclass)->size()) {
         std::cerr << "size mismatch in " << __FILE__ << ' ' << __LINE__ << std::endl; exit(1);
      }
      #endif

      // Set block parameters
      for (size_t b=0; b<nBlocks; ++b) {
         vmesh::GlobalID VBC_GID = blocks.at(b);
         this->get_velocity_mesh(popID, timeclass)->getBlockInfo(VBC_GID, parameters);
         parameters += BlockParams::N_VELOCITY_BLOCK_PARAMS;
      }

      //copy avgs data, here a conversion may happen between float and double
      Realf *cellBlockData = this->get_velocity_blocks(popID, timeclass)->getData(startLID);
      for(uint64_t i = 0; i< WID3 * nBlocks ; i++){
         cellBlockData[i] = avgBuffer[i];
      }
   }

   /*!
    Removes given block from the velocity grid.
    Does nothing if given block doesn't exist.
    */
   inline void SpatialCell::remove_velocity_block(const vmesh::GlobalID& block,const uint popID, const int timeclass = -1) {
      debug_population_check(popID);
      if (block == invalid_global_id()) {
         //std::cerr << "not removing, block " << block << " is invalid" << std::endl;
         return;
      }

      const vmesh::LocalID removedLID = this->get_velocity_mesh(popID, timeclass)->getLocalID(block);
      if (removedLID == invalid_local_id()) {
         //std::cerr << "not removing since block " << block << " does not exist" << std::endl;
         return;
      }

      // Get local ID of the last block:
      const vmesh::LocalID lastLID = this->get_velocity_mesh(popID, timeclass)->size()-1;
      // If block to remove is already last:
      if (lastLID == removedLID) {
         // Just remove the block
         this->get_velocity_mesh(popID, timeclass)->pop();
         this->get_velocity_blocks(popID, timeclass)->pop();
      } else {
         // Move the last block to the removed position
         this->get_velocity_mesh(popID, timeclass)->move(lastLID,removedLID);
         this->get_velocity_blocks(popID, timeclass)->move(lastLID,removedLID);
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

   inline std::set<int> SpatialCell::get_all_ghosts(){
      std::set<int> allghosts = this->requested_timeclass_ghosts;
      allghosts.insert(this->requested_timeclass_copy_ghosts.begin(),this->requested_timeclass_copy_ghosts.end());
      return allghosts;
   }

   // Used inside a population loop -> no pop information
   typedef struct {
      int timeclass;
      SpatialCell* cellptr;
      Real dt;
      // int step; // unused
   } AccelerationPayload;

} // namespaces

#endif
