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

#include <Eigen/Geometry>

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
      Real RHOLOSSADJUST = 0.0;      /*!< Counter for particle number loss from the destroying blocks in blockadjustment*/
      Real max_dt[2];                                                /**< Element[0] is max_r_dt, element[1] max_v_dt.*/
      Real velocityBlockMinValue;

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

      Eigen::Transform<Real,3, Eigen::Affine> fwd_transform;
      Eigen::Transform<Real,3, Eigen::Affine> bwd_transform;

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
      const Population& operator=(const Population& other) {
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

      // Following functions return velocity grid metadata //
      template<int PAD> void fetch_data(const vmesh::GlobalID& blockGID,const vmesh::VelocityMesh* vmesh,
                                        const Realf* src,Realf* array);
      template<int PAD>	void fetch_acc_data(const vmesh::GlobalID& blockGID,const int& dim,
					    vmesh::VelocityMesh* vmesh,
					    const Realf* src,Realf* array,Real cellSizeFractions[2]);

      vmesh::GlobalID find_velocity_block(vmesh::GlobalID cellIndices[3],const uint popID);
      Realf* get_data(const uint popID);
      const Realf* get_data(const uint popID) const;
      Realf* get_data(const vmesh::LocalID& blockLID,const uint popID);
      const Realf* get_data(const vmesh::LocalID& blockLID,const uint popID) const;
      Real* get_block_parameters(const uint popID);
      const Real* get_block_parameters(const uint popID) const;
      Real* get_block_parameters(const vmesh::LocalID& blockLID,const uint popID);
      const Real* get_block_parameters(const vmesh::LocalID& blockLID,const uint popID) const;

      Real* get_cell_parameters();
      const Real* get_cell_parameters() const;

      vmesh::LocalID get_number_of_velocity_blocks(const uint popID) const;
      vmesh::LocalID get_number_of_all_velocity_blocks() const;
      int get_number_of_populations() const;
      void debug_population_check(const uint popID) const;
      void debug_population_check(const uint popID, const vmesh::LocalID blockLID) const;

      Population & get_population(const uint popID);
      const Population & get_population(const uint popID) const;
      void set_population(const Population& pop, cuint popID);
      void scale_population(creal factor, cuint popID);
      void increment_population(const Population& pop, creal factor, cuint popID);

      std::vector<Population>& get_populations();
      const std::vector<Population>& get_populations() const;

      const Real& get_max_r_dt(const uint popID) const;
      const Real& get_max_v_dt(const uint popID) const;

      const std::array<uint32_t, 3>& get_velocity_grid_length(const uint popID);
      void get_velocity_block_coordinates(const uint popID,const vmesh::GlobalID& globalID,Real* coords);
      velocity_block_indices_t get_velocity_block_indices(const uint popID,const vmesh::GlobalID globalID);
      vmesh::GlobalID get_velocity_block(const uint popID,vmesh::GlobalID blockIndices[3]) const;
      vmesh::GlobalID get_velocity_block(const uint popID,const velocity_block_indices_t indices) const;
      vmesh::GlobalID get_velocity_block(const uint popID,const Real* coords) const;
      vmesh::GlobalID get_velocity_block(const uint popID,const Real vx,const Real vy,const Real vz) const;
      vmesh::GlobalID get_velocity_block_child(const uint popID,const vmesh::GlobalID& blockGID,
                                               const int& i_cell,const int& j_cell,const int& k_cell);
      void get_velocity_block_children_local_ids(const vmesh::GlobalID& blockGID,
                                                 std::vector<vmesh::LocalID>& childrenLIDs,
                                                 const uint popID);
      vmesh::GlobalID get_velocity_block_parent(const uint popID,const vmesh::GlobalID& blockGID);
      vmesh::GlobalID get_velocity_block_global_id(const vmesh::LocalID& blockLID,const uint popID) const;
      vmesh::LocalID get_velocity_block_local_id(const vmesh::GlobalID& blockGID,const uint popID) const;
      void get_velocity_block_size(const uint popID,const vmesh::GlobalID block,Real size[3]);
      Real get_velocity_block_vx_min(const uint popID,const vmesh::GlobalID block) const;
      Real get_velocity_block_vx_max(const uint popID,const vmesh::GlobalID block) const;
      Real get_velocity_block_vy_min(const uint popID,const vmesh::GlobalID block) const;
      Real get_velocity_block_vy_max(const uint popID,const vmesh::GlobalID block) const;
      Real get_velocity_block_vz_min(const uint popID,const vmesh::GlobalID block) const;
      Real get_velocity_block_vz_max(const uint popID,const vmesh::GlobalID block) const;

      static unsigned int invalid_block_index();
      static vmesh::GlobalID invalid_global_id();
      static vmesh::LocalID invalid_local_id();

      size_t count(const vmesh::GlobalID& block,const uint popID) const;

      void add_values(const vmesh::GlobalID& targetGID,
		      std::unordered_map<vmesh::GlobalID,Realf[(WID+2)*(WID+2)*(WID+2)]>& sourceData,
                      const uint popID);

      void printMeshSizes();
      static bool setCommunicatedSpecies(const uint popID);

      // Following functions adjust velocity blocks stored on the cell //
      bool add_velocity_block(const vmesh::GlobalID& block,const uint popID);
      void adjustSingleCellVelocityBlocks(const uint popID, bool doDeleteEmpty=false);
      void adjust_velocity_blocks(const std::vector<SpatialCell*>& spatial_neighbors,
                                  const uint popID,
                                  bool doDeleteEmptyBlocks=true);
      // Templated function for storing a v-space read from a file or generated elsewhere
      template <typename fileReal> void add_velocity_blocks(const uint popID,const std::vector<vmesh::GlobalID>& blocks,fileReal* avgBuffer);

      void update_velocity_block_content_lists(const uint popID);
      bool checkMesh(const uint popID);
      void clear(const uint popID, bool shrink=false);
      void setNewSizeClear(const uint popID, const vmesh::LocalID& newSize);
      void setNewSizeClear(const uint popID);

      uint64_t get_cell_memory_capacity();
      uint64_t get_cell_memory_size();
      void merge_values(const uint popID);
      void prepare_to_receive_blocks(const uint popID);
      bool shrink_to_fit();
      size_t size(const uint popID) const;
      void remove_velocity_block(const vmesh::GlobalID& block,const uint popID);
      vmesh::VelocityMesh* get_velocity_mesh(const size_t& popID);
      vmesh::VelocityBlockContainer* get_velocity_blocks(const size_t& popID);
      const vmesh::VelocityBlockContainer* get_velocity_blocks(const size_t& popID) const;

      void set_max_r_dt(const uint popID,const Real& value);
      void set_max_v_dt(const uint popID,const Real& value);

      // Following functions are related to MPI //
      std::tuple<void*, int, MPI_Datatype> get_mpi_datatype(const CellID cellID,const int sender_rank,const int receiver_rank,
                                                            const bool receiving,const int neighborhood);
      static uint64_t get_mpi_transfer_type(void);
      static void set_mpi_transfer_type(const uint64_t type,bool atSysBoundaries=false);
      static void set_mpi_transfer_direction(const int dimension);
      void set_mpi_transfer_enabled(bool transferEnabled);
      void updateSparseMinValue(const uint popID);
      Real getVelocityBlockMinValue(const uint popID) const;

      // Member variables //
      std::array<Real, vderivatives::N_V_DERIVATIVES> derivativesV;           /**< Derivatives of V for vorticity AMR.*/
      std::array<Real, bvolderivatives::N_BVOL_DERIVATIVES> derivativesBVOL;  /**< Derivatives of BVOL needed by the acceleration.
                                                                               * Separate array because it does not need to be communicated.*/
      std::array<Real, CellParams::N_SPATIAL_CELL_PARAMS> parameters;         /**< Bulk variables in this spatial cell.*/
      std::array<Realf, WID3> null_block_data;

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

      //SpatialCell& operator=(const SpatialCell& other);
    private:
      //SpatialCell& operator=(const SpatialCell&);

      bool compute_block_has_content(const vmesh::GlobalID& block,const uint popID) const;

      static int activePopID;
      bool initialized;
      bool mpiTransferEnabled;

      std::vector<spatial_cell::Population> populations;                        /**< Particle population variables.*/
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
   inline void SpatialCell::debug_population_check(const uint popID, const vmesh::LocalID blockLID) const {
      debug_population_check(popID);
      #ifdef DEBUG_SPATIAL_CELL
      if (blockLID >= populations[popID].blockContainer->size()) {
         std::cerr << "ERROR, block LID out of bounds, blockContainer->size() " << populations[popID].blockContainer->size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
         exit(1);
      }
      #endif
   }

   inline vmesh::GlobalID SpatialCell::find_velocity_block(vmesh::GlobalID cellIndices[3],const uint popID) {
      debug_population_check(popID);
      return populations[popID].vmesh->findBlock(cellIndices);
   }

   inline Realf* SpatialCell::get_data(const uint popID) {
      debug_population_check(popID);
      return populations[popID].blockContainer->getData();
   }

   inline const Realf* SpatialCell::get_data(const uint popID) const {
      debug_population_check(popID);
      return populations[popID].blockContainer->getData();
   }

   inline Realf* SpatialCell::get_data(const vmesh::LocalID& blockLID,const uint popID) {
      debug_population_check(popID,blockLID);
      if (blockLID == vmesh::VelocityMesh::invalidLocalID()) {
         return null_block_data.data();
      }
      return populations[popID].blockContainer->getData(blockLID);
   }

   inline const Realf* SpatialCell::get_data(const vmesh::LocalID& blockLID,const uint popID) const {
      debug_population_check(popID,blockLID);
      if (blockLID == vmesh::VelocityMesh::invalidLocalID()) {
         return null_block_data.data();
      }
      return populations[popID].blockContainer->getData(blockLID);
   }

   inline Real* SpatialCell::get_block_parameters(const uint popID) {
      debug_population_check(popID);
      return populations[popID].blockContainer->getParameters();
   }

   inline const Real* SpatialCell::get_block_parameters(const uint popID) const {
      debug_population_check(popID);
      return populations[popID].blockContainer->getParameters();
   }

   inline Real* SpatialCell::get_block_parameters(const vmesh::LocalID& blockLID,const uint popID) {
      debug_population_check(popID,blockLID);
      return populations[popID].blockContainer->getParameters(blockLID);
   }

   inline const Real* SpatialCell::get_block_parameters(const vmesh::LocalID& blockLID,const uint popID) const {
      debug_population_check(popID,blockLID);
      return populations[popID].blockContainer->getParameters(blockLID);
   }

   inline Real* SpatialCell::get_cell_parameters() {
      return parameters.data();
   }

   inline const Real* SpatialCell::get_cell_parameters() const {
      return parameters.data();
   }

   inline vmesh::LocalID SpatialCell::get_number_of_velocity_blocks(const uint popID) const {
      debug_population_check(popID);
      return populations[popID].blockContainer->size();
   }

    /** Get the total number of velocity blocks in this cell, summed over
     * all existing particle populations.
     * @return Total number of velocity blocks in the cell.*/
    inline vmesh::LocalID SpatialCell::get_number_of_all_velocity_blocks() const {
        vmesh::LocalID N_blocks = 0;
        for (size_t p=0; p<populations.size(); ++p)
            N_blocks += populations[p].blockContainer->size();
        return N_blocks;
    }

   inline int SpatialCell::get_number_of_populations() const {
      return populations.size();
   }

   inline Population & SpatialCell::get_population(const uint popID) {
      return populations[popID];
   }

   inline const Population & SpatialCell::get_population(const uint popID) const {
      return populations[popID];
   }

   inline void SpatialCell::set_population(const Population& pop, cuint popID) {
      this->populations[popID] = pop;
   }
   inline void SpatialCell::scale_population(creal factor, cuint popID) {
      (this->populations[popID]).Scale(factor);
   }
   inline void SpatialCell::increment_population(const Population& pop, creal factor, cuint popID) {
      (this->populations[popID]).Increment(pop, factor);
   }
   inline std::vector<Population>& SpatialCell::get_populations() {
      return populations;
   }
   inline const std::vector<Population>& SpatialCell::get_populations() const {
      return populations;
   }

   inline const std::array<uint32_t, 3>& SpatialCell::get_velocity_grid_length(const uint popID) {
      return populations[popID].vmesh->getGridLength();
   }

   inline void SpatialCell::get_velocity_block_coordinates(const uint popID,const vmesh::GlobalID& globalID,Real* coords) {
      populations[popID].vmesh->getBlockCoordinates(globalID,coords);
   }

   /*!
    Returns the indices of given velocity block
    */
   inline velocity_block_indices_t SpatialCell::get_velocity_block_indices(const uint popID,const vmesh::GlobalID block) {
      velocity_block_indices_t indices;
      populations[popID].vmesh->getIndices(block,indices[0],indices[1],indices[2]);
      return indices;
   }

   /*!
    Returns the velocity block at given indices or error_velocity_block
    */
   inline vmesh::GlobalID SpatialCell::get_velocity_block(const uint popID,const velocity_block_indices_t indices) const {
      return populations[popID].vmesh->getGlobalID(indices[0],indices[1],indices[2]);
   }

   inline vmesh::GlobalID SpatialCell::get_velocity_block(const uint popID,vmesh::GlobalID blockIndices[3]) const {
      return populations[popID].vmesh->getGlobalID(blockIndices[0],blockIndices[1],blockIndices[2]);
   }

   /*!
    Returns the velocity block at given location or
    error_velocity_block if outside of the velocity grid
    */
   inline vmesh::GlobalID SpatialCell::get_velocity_block(const uint popID,const Real vx,const Real vy,const Real vz) const {
      Real coords[3] = {vx,vy,vz};
      return populations[popID].vmesh->getGlobalID(coords);
   }

   inline vmesh::GlobalID SpatialCell::get_velocity_block(const uint popID,const Real* coords) const {
      return populations[popID].vmesh->getGlobalID(coords);
   }

   inline vmesh::GlobalID SpatialCell::get_velocity_block_global_id(const vmesh::LocalID& blockLID,const uint popID) const {
      debug_population_check(popID);
      return populations[popID].vmesh->getGlobalID(blockLID);
   }

   inline vmesh::LocalID SpatialCell::get_velocity_block_local_id(const vmesh::GlobalID& blockGID,const uint popID) const {
      debug_population_check(popID);
      return populations[popID].vmesh->getLocalID(blockGID);
   }

   inline void SpatialCell::get_velocity_block_size(const uint popID,const vmesh::GlobalID block,Real blockSize[3]) {
      populations[popID].vmesh->getBlockSize(block,blockSize);
   }

   /*!
    Returns the edge where given velocity block starts.
    */
   inline Real SpatialCell::get_velocity_block_vx_min(const uint popID,const vmesh::GlobalID block) const {
      Real coords[3];
      populations[popID].vmesh->getBlockCoordinates(block,coords);
      return coords[0];
   }

   /*!
    Returns the edge where given velocity block ends.
    */
   inline Real SpatialCell::get_velocity_block_vx_max(const uint popID,const vmesh::GlobalID block) const {
      Real coords[3];
      populations[popID].vmesh->getBlockCoordinates(block,coords);

      Real size[3];
      populations[popID].vmesh->getBlockSize(block,size);
      return coords[0]+size[0];
   }

   /*!
    Returns the edge where given velocity block starts.
    */
   inline Real SpatialCell::get_velocity_block_vy_min(const uint popID,const vmesh::GlobalID block) const {
      Real coords[3];
      populations[popID].vmesh->getBlockCoordinates(block,coords);
      return coords[1];
   }

   /*!
    Returns the edge where given velocity block ends.
    */
   inline Real SpatialCell::get_velocity_block_vy_max(const uint popID,const vmesh::GlobalID block) const {
      Real coords[3];
      populations[popID].vmesh->getBlockCoordinates(block,coords);

      Real size[3];
      populations[popID].vmesh->getBlockSize(block,size);
      return coords[1]+size[1];
   }

   /*!
    Returns the edge where given velocity block starts.
    */
   inline Real SpatialCell::get_velocity_block_vz_min(const uint popID,const vmesh::GlobalID block) const {
      Real coords[3];
      populations[popID].vmesh->getBlockCoordinates(block,coords);
      return coords[2];
   }

   /*!
    Returns the edge where given velocity block ends.
    */
   inline Real SpatialCell::get_velocity_block_vz_max(const uint popID,const vmesh::GlobalID block) const {
      Real coords[3];
      populations[popID].vmesh->getBlockCoordinates(block,coords);

      Real size[3];
      populations[popID].vmesh->getBlockSize(block,size);
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
   inline size_t SpatialCell::count(const vmesh::GlobalID& block,const uint popID) const {
      debug_population_check(popID);
      return populations[popID].vmesh->count(block);
   }

   /*!
    Returns the number of existing velocity blocks.
    */
   inline size_t SpatialCell::size(const uint popID) const {
      debug_population_check(popID);
      return populations[popID].vmesh->size();
   }

   inline vmesh::VelocityMesh* SpatialCell::get_velocity_mesh(const size_t& popID) {
      debug_population_check(popID);
      return populations[popID].vmesh;
   }

   inline vmesh::VelocityBlockContainer* SpatialCell::get_velocity_blocks(const size_t& popID) {
      debug_population_check(popID);
      return populations[popID].blockContainer;
   }
   inline const vmesh::VelocityBlockContainer* SpatialCell::get_velocity_blocks(const size_t& popID) const {
      debug_population_check(popID);
      return populations[popID].blockContainer;
   }

   inline bool SpatialCell::checkMesh(const uint popID) {
      debug_population_check(popID);
      return populations[popID].vmesh->check();
   }

   /*!
    * Removes all velocity blocks from this spatial cell and frees memory in the cell
   */
   inline void SpatialCell::clear(const uint popID, bool shrink) {
      debug_population_check(popID);
      populations[popID].vmesh->clear(shrink);
      populations[popID].blockContainer->clear(shrink);
    }

   /*!
     Ensures the selected population vmesh localToGlobalMap has sufficient capacity and is of correct size. Does not alter
     its contents. Ensures the selected population vmesh globalToLocalMap has sufficient capacity and is empty.
     Ensures the selected population VBC has sufficient capacity and is of correct size.
   */
   inline void SpatialCell::setNewSizeClear(const uint popID, const vmesh::LocalID& newSize) {
      populations[popID].ResizeClear(newSize);
   }
   inline void SpatialCell::setNewSizeClear(const uint popID) {
      populations[popID].ResizeClear(populations[popID].N_blocks);
   }

   /*!
    Return the memory consumption in bytes as reported using the size()
    functions of the containers in spatial cell
    */
   inline uint64_t SpatialCell::get_cell_memory_size() {
      uint64_t size = 0;
      size += 2 * WID3 * sizeof(Realf);
      size += velocity_block_with_content_list->size() * sizeof(vmesh::GlobalID);
      size += velocity_block_with_no_content_list->size() * sizeof(vmesh::GlobalID);
      size += CellParams::N_SPATIAL_CELL_PARAMS * sizeof(Real);
      size += bvolderivatives::N_BVOL_DERIVATIVES * sizeof(Real);

      for (size_t p=0; p<populations.size(); ++p) {
          size += populations[p].vmesh->sizeInBytes();
          size += populations[p].blockContainer->sizeInBytes();
      }

      return size;
   }

   /*!
    Return the memory consumption in bytes as reported using
    the size() functions of the containers in spatial cell
    */
   inline uint64_t SpatialCell::get_cell_memory_capacity() {
      uint64_t capacity = 0;

      capacity += 2 * WID3 * sizeof(Realf);
      capacity += velocity_block_with_content_list->capacity()  * sizeof(vmesh::GlobalID);
      capacity += velocity_block_with_no_content_list->capacity()  * sizeof(vmesh::GlobalID);
      capacity += CellParams::N_SPATIAL_CELL_PARAMS * sizeof(Real);
      capacity += bvolderivatives::N_BVOL_DERIVATIVES * sizeof(Real);

      for (size_t p=0; p<populations.size(); ++p) {
        capacity += populations[p].vmesh->capacityInBytes();
        capacity += populations[p].blockContainer->capacityInBytes();
      }

      return capacity;
   }

   /*!
    Adds an empty velocity block into this spatial cell.
    Returns true if given block was added or already exists.
    Returns false if given block is invalid or would be outside
    of the velocity grid.
    */
   inline bool SpatialCell::add_velocity_block(const vmesh::GlobalID& block,const uint popID) {
      debug_population_check(popID);
      // Block insert will fail, if the block already exists, or if
      // there are too many blocks in the spatial cell
      bool success = true;
      if (populations[popID].vmesh->push_back(block) == false) {
         return false;
      }

      const vmesh::LocalID VBC_LID = populations[popID].blockContainer->push_back_and_zero();

      // Set block parameters:
      Real* parameters = get_block_parameters(VBC_LID,popID);
      populations[popID].vmesh->getBlockInfo(block, parameters);

      // The following call 'should' be the fastest, but is actually
      // much slower that the parameter setting above
      //vmesh::VelocityMesh::getBlockInfo(block,get_block_parameters( blockContainer->push_back() ));
      return success;
   }

   /** Adds a vector of velocity blocks to the population, sets the parameters, and fills the data
       with phase-space densities from the provided buffer (which was read from a file).
   */
   template <typename fileReal> void SpatialCell::add_velocity_blocks(const uint popID,const std::vector<vmesh::GlobalID>& blocks,fileReal* avgBuffer) {
      debug_population_check(popID);
      // Return if no blocks to add
      const uint nBlocks = blocks.size();
      if (nBlocks==0) {
         return;
      }
      populations[popID].vmesh->setNewCapacity(nBlocks);
      populations[popID].blockContainer->setNewCapacity(nBlocks);

      // Add blocks to mesh
      const uint8_t adds = populations[popID].vmesh->push_back(blocks);
      if (adds == 0) {
         std::cerr << "Failed to add blocks" << std::endl;
         return;
      }

      // Add blocks to block container
      vmesh::LocalID startLID = populations[popID].blockContainer->push_back(nBlocks);
      Real* parameters = populations[popID].blockContainer->getParameters(startLID);

      #ifdef DEBUG_SPATIAL_CELL
      if (populations[popID].vmesh->size() != populations[popID].blockContainer->size()) {
         std::cerr << "size mismatch in " << __FILE__ << ' ' << __LINE__ << std::endl; exit(1);
      }
      #endif

      // Set block parameters
      for (size_t b=0; b<nBlocks; ++b) {
         vmesh::GlobalID VBC_GID = blocks.at(b);
         populations[popID].vmesh->getBlockInfo(VBC_GID, parameters);
         parameters += BlockParams::N_VELOCITY_BLOCK_PARAMS;
      }

      //copy avgs data, here a conversion may happen between float and double
      Realf *cellBlockData = populations[popID].blockContainer->getData(startLID);
      for(uint64_t i = 0; i< WID3 * nBlocks ; i++){
         cellBlockData[i] = avgBuffer[i];
      }
   }

   /*!
    Removes given block from the velocity grid.
    Does nothing if given block doesn't exist.
    */
   inline void SpatialCell::remove_velocity_block(const vmesh::GlobalID& block,const uint popID) {
      debug_population_check(popID);
      if (block == invalid_global_id()) {
         //std::cerr << "not removing, block " << block << " is invalid" << std::endl;
         return;
      }

      const vmesh::LocalID removedLID = populations[popID].vmesh->getLocalID(block);
      if (removedLID == invalid_local_id()) {
         //std::cerr << "not removing since block " << block << " does not exist" << std::endl;
         return;
      }

      // Get local ID of the last block:
      const vmesh::LocalID lastLID = populations[popID].vmesh->size()-1;
      // If block to remove is already last:
      if (lastLID == removedLID) {
         // Just remove the block
         populations[popID].vmesh->pop();
         populations[popID].blockContainer->pop();
      } else {
         // Move the last block to the removed position
         populations[popID].vmesh->move(lastLID,removedLID);
         populations[popID].blockContainer->move(lastLID,removedLID);
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

} // namespaces

#endif
