#include "cpu_1d_ppm_nonuniform.hpp"
//#include "cpu_1d_ppm_nonuniform_conserving.hpp"
#include "vec.h"
#include "../grid.h"
#include "../object_wrapper.h"
#include "../memoryallocation.h"
#include "cpu_trans_map_amr.hpp"
#include "cpu_trans_pencils.hpp"

using namespace std;
using namespace spatial_cell;

// indices in padded source block, which is of type Vec with VECL
// elements in each vector.

#define i_trans_ps_blockv_pencil(planeVectorIndex, planeIndex, blockIndex, lengthOfPencil) ( (blockIndex)  +  ( (planeVectorIndex) + (planeIndex) * VEC_PER_PLANE ) * ( lengthOfPencil) )

inline bool check_skip_remapping(Vec* values) {
   for (int index=-VLASOV_STENCIL_WIDTH; index<VLASOV_STENCIL_WIDTH+1; ++index) {
      if (horizontal_or(values[index] > Vec(0))) return false;
   }
   return true;
}

/* Propagate a given velocity block in all spatial cells of a pencil by a time step dt using a PPM reconstruction.
 *
 * @param dz Width of spatial cells in the direction of the pencil, vector datatype
 * @param values Density values of the block, vector datatype
 * @param dimension Satial dimension
 * @param blockGID Global ID of the velocity block.
 * @param dt Time step
 * @param vmesh Velocity mesh object
 * @param lengthOfPencil Number of cells in the pencil
 */
void propagatePencil(
   Realf* dz,
   Vec* values, // Vec-ordered block data values for pencils
   const uint dimension,
   const uint blockGID,
   const Realv dt,
   const vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID> &vmesh,
   const uint lengthOfPencil,
   const Realv threshold,
   Realf** blockDataPointer, // Spacing is for sources, but will be written into
   Realf* targetRatios, // Vector holding target ratios
   const unsigned int* const cellid_transpose
) {
   // Get velocity data from vmesh that we need later to calculate the translation
   velocity_block_indices_t block_indices;
   uint8_t refLevel;
   vmesh.getIndices(blockGID,refLevel, block_indices[0], block_indices[1], block_indices[2]);
   Realv dvz = vmesh.getCellSize(refLevel)[dimension];
   Realv vz_min = vmesh.getMeshMinLimits()[dimension];

   // Assuming 1 neighbor in the target array because of the CFL condition
   // In fact propagating to > 1 neighbor will give an error
   // Also defined in the calling function for the allocation of targetValues

   // Go over length of propagated cells
   for (int i = VLASOV_STENCIL_WIDTH; i < (int)lengthOfPencil-VLASOV_STENCIL_WIDTH; i++){
      // Get pointers to block data used for output.
      // CUDATODO: use blockGID to get pointers here
      Realf* block_data_m1 = blockDataPointer[i - 1];
      Realf* block_data    = blockDataPointer[i];
      Realf* block_data_p1 = blockDataPointer[i + 1];
      // Cells which shouldn't be written to (e.g. sysboundary cells) have a targetRatio of 0
      // Also need to check if pointer is valid, because a cell can be missing an elsewhere propagated block
      Realf areaRatio_m1 = targetRatios[i - 1];
      Realf areaRatio    = targetRatios[i];
      Realf areaRatio_p1 = targetRatios[i + 1];

      Realf vector[VECL];
      // Loop over planes
      for (uint k = 0; k < WID; ++k) {
         const Realv cell_vz = (block_indices[dimension] * WID + k + 0.5) * dvz + vz_min; //cell centered velocity
         const Vec z_translation = cell_vz * dt / dz[i]; // how much it moved in time dt (reduced units)

         // Determine direction of translation
         // part of density goes here (cell index change along spatial direcion)
         Vecb positiveTranslationDirection = (z_translation > Vec(0.0));

         // Calculate normalized coordinates in current cell.
         // The coordinates (scaled units from 0 to 1) between which we will
         // integrate to put mass in the target  neighboring cell.
         // Normalize the coordinates to the origin cell. Then we scale with the difference
         // in volume between target and origin later when adding the integrated value.
         Vec z_1,z_2;
         z_1 = select(positiveTranslationDirection, 1.0 - z_translation, 0.0);
         z_2 = select(positiveTranslationDirection, 1.0, - z_translation);

         // if( horizontal_or(abs(z_1) > Vec(1.0)) || horizontal_or(abs(z_2) > Vec(1.0)) ) {
         //    std::cout << "Error, CFL condition violated\n";
         //    std::cout << "Exiting\n";
         //    std::exit(1);
         // }

         // Loop over Vec's in current plance
         for (uint planeVector = 0; planeVector < VEC_PER_PLANE; planeVector++) {
            // Check if all values are 0:
            if (check_skip_remapping(values + i_trans_ps_blockv_pencil(planeVector, k, i, lengthOfPencil))) continue;

            // Compute polynomial coefficients
            Vec a[3];
            // Silly indexing into coefficient calculation necessary due to built-in assumptions of unsigned indexing.
            compute_ppm_coeff_nonuniform(dz + i - VLASOV_STENCIL_WIDTH,
                                         values + i_trans_ps_blockv_pencil(planeVector, k, i, lengthOfPencil) - VLASOV_STENCIL_WIDTH,
                                         h4, VLASOV_STENCIL_WIDTH, a, threshold);

            // Compute integral
            const Vec ngbr_target_density =
               z_2 * ( a[0] + z_2 * ( a[1] + z_2 * a[2] ) ) -
               z_1 * ( a[0] + z_1 * ( a[1] + z_1 * a[2] ) );

            // Store mapped density in two target cells
            // in the current original cells we will put the rest of the original density
            if (areaRatio && block_data) {
               const Vec selfContribution = (values[i_trans_ps_blockv_pencil(planeVector, k, i, lengthOfPencil)] - ngbr_target_density) * areaRatio;
               selfContribution.store(vector);
               // Loop over 3rd (vectorized) vspace dimension
               #pragma omp simd
               for (uint iv = 0; iv < VECL; iv++) {
                  block_data[cellid_transpose[iv + planeVector * VECL + k * WID2]] += vector[iv];
               }
            }
            if (areaRatio_p1 && block_data_p1) {
               const Vec p1Contribution = select(positiveTranslationDirection, ngbr_target_density
                                                 * dz[i] / dz[i + 1], Vec(0.0)) * areaRatio_p1;
               p1Contribution.store(vector);
               // Loop over 3rd (vectorized) vspace dimension
               #pragma omp simd
               for (uint iv = 0; iv < VECL; iv++) {
                  block_data_p1[cellid_transpose[iv + planeVector * VECL + k * WID2]] += vector[iv];
               }
            }
            if (areaRatio_m1 && block_data_m1) {
               const Vec m1Contribution = select(!positiveTranslationDirection, ngbr_target_density
                                                 * dz[i] / dz[i - 1], Vec(0.0)) * areaRatio_m1;
               m1Contribution.store(vector);
               // Loop over 3rd (vectorized) vspace dimension
               #pragma omp simd
               for (uint iv = 0; iv < VECL; iv++) {
                  block_data_m1[cellid_transpose[iv + planeVector * VECL + k * WID2]] += vector[iv];
               }
            }
         }
      }
   }
}

/* Copy the pencil source data to the temporary values array, so that the
 * dimensions are correctly swapped.
 *
 * This function must be thread-safe.
 *
 * @param blockDataPointer Vector of pre-prepared pointers to input (cell) block data
 * @param start Index from blockDataPointer to start at
 * @param int lengthOfPencil Number of spatial cells in pencil (not inclusive 2*VLASOV_STENCIL_WIDTH
 * @param values Vector into which the data should be loaded
 * @param cellid_transpose
 * @param popID ID of the particle species.
 */
bool copy_trans_block_data_amr(
   Realf** pencilBlockData,
   const int lengthOfPencil,
   Vec* values,
   const unsigned int* const cellid_transpose,
   const uint popID) {

   //  Copy volume averages of this block from all spatial cells:
   for (int b = 0; b < lengthOfPencil; b++) {
      if(pencilBlockData[b] != NULL) {
         Realf blockValues[WID3];
         Realf* block_data = pencilBlockData[b];
         // Copy data to a temporary array and transpose values so that mapping is along k direction.
         #pragma omp simd
         for (uint i=0; i<WID3; ++i) {
            blockValues[i] = block_data[cellid_transpose[i]];
         }
         // now load values into the actual values table..
         uint offset =0;
         for (uint k=0; k<WID; k++) {
            for(uint planeVector = 0; planeVector < VEC_PER_PLANE; planeVector++){
               // store data, when reading data from data we swap dimensions
               // using precomputed plane_index_to_id and cell_indices_to_id
               values[i_trans_ps_blockv_pencil(planeVector, k, b, lengthOfPencil)].load(blockValues + offset);
               offset += VECL;
            }
         }
      } else {
         for (uint k=0; k<WID; ++k) {
            for(uint planeVector = 0; planeVector < VEC_PER_PLANE; planeVector++) {
               values[i_trans_ps_blockv_pencil(planeVector, k, b, lengthOfPencil)] = Vec(0);
            }
         }
      }
   }
   return true;
}

/* Map velocity blocks in all local cells forward by one time step in one spatial dimension.
 * This function uses 1-cell wide pencils to update cells in-place to avoid allocating large
 * temporary buffers.
 *
 * @param [in] mpiGrid DCCRG grid object
 * @param [in] localPropagatedCells List of local cells that get propagated
 * ie. not boundary or DO_NOT_COMPUTE
 * @param [in] remoteTargetCells List of non-local target cells
 * @param [in] dimension Spatial dimension
 * @param [in] dt Time step
 * @param [in] popId Particle population ID
 */
bool trans_map_1d_amr(const dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                      const vector<CellID>& localPropagatedCells,
                      const vector<CellID>& remoteTargetCells,
                      std::vector<uint>& nPencils,
                      const uint dimension,
                      const Realv dt,
                      const uint popID) {
   /***********************/
   phiprof::Timer setupTimer {"trans-amr-setup"};
   /***********************/
   // return if there's no cells to propagate
   if(localPropagatedCells.size() == 0) {
      return false;
   }

   uint cell_indices_to_id[3]; /*< used when computing id of target cell in block*/
   unsigned int cellid_transpose[WID3]; /*< defines the transpose for the solver internal (transposed) id: i + j*WID + k*WID2 to actual one*/
   // Fiddle indices x,y,z in VELOCITY SPACE
   switch (dimension) {
   case 0:
      // set values in array that is used to convert block indices
      // to global ID using a dot product.
      cell_indices_to_id[0]=WID2;
      cell_indices_to_id[1]=WID;
      cell_indices_to_id[2]=1;
      break;
   case 1:
      // set values in array that is used to convert block indices
      // to global ID using a dot product
      cell_indices_to_id[0]=1;
      cell_indices_to_id[1]=WID2;
      cell_indices_to_id[2]=WID;
      break;
   case 2:
      // set values in array that is used to convert block indices
      // to global id using a dot product.
      cell_indices_to_id[0]=1;
      cell_indices_to_id[1]=WID;
      cell_indices_to_id[2]=WID2;
      break;
   default:
      cerr << __FILE__ << ":"<< __LINE__ << " Wrong dimension, abort"<<endl;
      abort();
      break;
   }

   // Vector with all cell ids
   vector<CellID> allCells(localPropagatedCells);
   allCells.insert(allCells.end(), remoteTargetCells.begin(), remoteTargetCells.end());

   // init cellid_transpose (moved here to take advantage of the omp parallel region)
   #pragma omp parallel for collapse(2) schedule(static)
   for (uint k=0; k<WID; ++k) {
      for (uint j=0; j<WID; ++j) {
         for (uint i=0; i<WID; ++i) {
            const uint cell =
               i * cell_indices_to_id[0] +
               j * cell_indices_to_id[1] +
               k * cell_indices_to_id[2];
            cellid_transpose[ i + j * WID + k * WID2] = cell;
         }
      }
   }

   // Only needed if pencil counts are used as weight multiplier in load balance
   if (Parameters::prepareForRebalance == true) {
      for (uint i=0; i<localPropagatedCells.size(); i++) {
         for (uint ip=0; ip<DimensionPencils[dimension].N; ip++) {
            // Read only central IDs for each pencil
            std::vector<CellID> centerIds = DimensionPencils[dimension].getIds(ip);
            cuint myPencilCount = std::count(centerIds.begin(), centerIds.end(), localPropagatedCells[i]);
            nPencils[i] += myPencilCount;
            nPencils[nPencils.size()-1] += myPencilCount;
         }
      }
   }

   // Get a pointer to the velocity mesh of the first spatial cell
   const vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& vmesh = mpiGrid[localPropagatedCells[0]]->get_velocity_mesh(popID);

   phiprof::Timer buildBlockListTimer {"trans-amr-buildBlockList"};
   // Get a unique sorted list of blockids that are in any of the
   // local target cells.
   // Note: Now only considers local propagated cells, as there's no point in
   // Propagating data which does not have an existing target block.
   std::vector<vmesh::GlobalID> unionOfBlocks;
   std::unordered_set<vmesh::GlobalID> unionOfBlocksSet;
#pragma omp parallel
   {
      std::unordered_set<vmesh::GlobalID> thread_unionOfBlocksSet;

#pragma omp for
      for(unsigned int i=0; i<allCells.size(); i++) {
         CellID cellid = allCells[i];
         const vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& cvmesh = mpiGrid[cellid]->get_velocity_mesh(popID);
         for (vmesh::LocalID block_i=0; block_i< cvmesh.size(); ++block_i) {
            thread_unionOfBlocksSet.insert(cvmesh.getGlobalID(block_i));
         }
      }
#pragma omp critical
      {
         unionOfBlocksSet.insert(thread_unionOfBlocksSet.begin(), thread_unionOfBlocksSet.end());
      } // pragma omp critical
   } // pragma omp parallel
   unionOfBlocks.insert(unionOfBlocks.end(), unionOfBlocksSet.begin(), unionOfBlocksSet.end());
   buildBlockListTimer.stop();
   /***********************/
   setupTimer.stop();
   /***********************/

   int mappingTimerId = phiprof::initializeTimer("trans-amr-mapping");
   int loadTimerId = phiprof::initializeTimer("trans-amr-load source data");
   int memsetTimerId = phiprof::initializeTimer("trans-amr-MemSet");
   int propagateTimerId = phiprof::initializeTimer("trans-amr-propagatePencil");

#pragma omp parallel
   {
      // Vector of pointers to cell block data, used for both reading and writing
      std::vector<Vec> blockDataBuffer(DimensionPencils[dimension].sumOfLengths*WID3/VECL);
      std::vector<Realf*> cellBlockData(DimensionPencils[dimension].sumOfLengths);
      std::vector<uint> pencilBlocksCount(DimensionPencils[dimension].N);

      // Loop over velocity space blocks (threaded).
#pragma omp for schedule(dynamic,1)
      for(uint blocki = 0; blocki < unionOfBlocks.size(); blocki++) {
         // Get global id of the velocity block
         vmesh::GlobalID blockGID = unionOfBlocks[blocki];

         phiprof::Timer mappingTimer {mappingTimerId}; // mapping (top-level)

         // Load data for pencils.
         phiprof::Timer loadTimer {loadTimerId};
         for (uint pencili = 0; pencili < DimensionPencils[dimension].N; ++pencili){
            int nonEmptyBlocks = 0;
            int L = DimensionPencils[dimension].lengthOfPencils[pencili];
            int start = DimensionPencils[dimension].idsStart[pencili];
            // Loop over cells in pencil
            for (int b = 0; b < L; b++) {
               // Get cell pointer and local block id
               SpatialCell* srcCell = mpiGrid[DimensionPencils[dimension].ids[start + b]];
               const vmesh::LocalID blockLID = srcCell->get_velocity_block_local_id(blockGID,popID);
               // Store block data pointer for both loading of data and writing back to the cell
               if (blockLID != srcCell->invalid_local_id()) {
                  // Get data pointer
                  cellBlockData[start + b] = srcCell->get_data(blockLID,popID);
                  nonEmptyBlocks++;
               } else {
                  cellBlockData[start + b] = NULL;
               }
            }
            if(nonEmptyBlocks == 0) {
               continue;
            }
            pencilBlocksCount.at(pencili) = nonEmptyBlocks;
            // Transpose and copy block data from cells to source buffer
            Vec* blockDataSource = blockDataBuffer.data() + start*WID3/VECL;
            Realf** pencilBlockData = cellBlockData.data() + start;
            copy_trans_block_data_amr(pencilBlockData, L, blockDataSource,
                                      cellid_transpose, popID);
         }
         loadTimer.stop();

         phiprof::Timer memsetTimer {memsetTimerId};
         // reset blocks in all non-sysboundary neighbor spatial cells for this block id
         for (CellID target_cell_id: DimensionTargetCells[dimension]) {
            SpatialCell* target_cell = mpiGrid[target_cell_id];
            if (target_cell) {
               // Get local velocity block id
               const vmesh::LocalID blockLID = target_cell->get_velocity_block_local_id(blockGID, popID);
               // Check for invalid block id
               if (blockLID != vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>::invalidLocalID()) {
                  // Get a pointer to the block data
                  Realf* blockData = target_cell->get_data(blockLID, popID);
                  memset(blockData, 0, WID3*sizeof(Realf));
               }
            }
         }
         memsetTimer.stop();

         phiprof::Timer propagateTimer {propagateTimerId};
         for(uint pencili = 0; pencili < DimensionPencils[dimension].N; ++pencili){
            // Skip pencils without blocks
            if (pencilBlocksCount.at(pencili) == 0) {
               continue;
            }
            // sourceVecData => targetBlockData[this pencil])
            int L = DimensionPencils[dimension].lengthOfPencils[pencili];
            int start = DimensionPencils[dimension].idsStart[pencili];
            // Dz and sourceVecData are both padded by VLASOV_STENCIL_WIDTH
            // Dz has 1 value/cell, sourceVecData has WID3 values/cell
            // vmesh is required just for general indexes and accessors
            Realv scalingthreshold = mpiGrid[DimensionPencils[dimension].ids[start + VLASOV_STENCIL_WIDTH]]->getVelocityBlockMinValue(popID);
            Realf* pencilDZ = DimensionPencils[dimension].sourceDZ.data() + start;
            Realf* pencilRatios = DimensionPencils[dimension].targetRatios.data() + start;
            Realf** pencilBlockData = cellBlockData.data() + start;
            Vec* blockDataSource = blockDataBuffer.data() +start*WID3/VECL;
            propagatePencil(pencilDZ,
                            blockDataSource,
                            dimension,
                            blockGID,
                            dt,
                            vmesh,
                            L,
                            scalingthreshold,
                            pencilBlockData,
                            pencilRatios,
                            cellid_transpose
               );
         }
         propagateTimer.stop();

         mappingTimer.stop(); // mapping (top-level)
      } // Closes loop over blocks
   } // closes pragma omp parallel

   return true;
}


/* Get an index that identifies which cell in the list of sibling cells this cell is.
 *
 * @param mpiGrid DCCRG grid object
 * @param cellid DCCRG id of this cell
 */
int get_sibling_index(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid, const CellID& cellid) {

   const int NO_SIBLINGS = 0;
   if(mpiGrid.get_refinement_level(cellid) == 0) {
      return NO_SIBLINGS;
   }

   //CellID parent = mpiGrid.mapping.get_parent(cellid);
   CellID parent = mpiGrid.get_parent(cellid);

   if (parent == INVALID_CELLID) {
      std::cerr<<"Invalid parent id"<<std::endl;
      abort();
   }

   // get_all_children returns an array instead of a vector now, need to map it to a vector for find and distance
   // std::array<uint64_t, 8> siblingarr = mpiGrid.mapping.get_all_children(parent);
   // vector<CellID> siblings(siblingarr.begin(), siblingarr.end());
   vector<CellID> siblings = mpiGrid.get_all_children(parent);
   auto location = std::find(siblings.begin(),siblings.end(),cellid);
   auto index = std::distance(siblings.begin(), location);
   if (index>7) {
      std::cerr<<"Invalid parent id"<<std::endl;
      abort();
   }
   return index;

}

/* This function communicates the mapping on process boundaries, and then updates the data to their correct values.
 * When sending data between neighbors of different refinement levels, special care has to be taken to ensure that
 * The sending and receiving ranks allocate the correct size arrays for neighbor_block_data.
 * This is partially due to DCCRG defining neighborhood size relative to the host cell. For details, see
 * https://github.com/fmihpc/dccrg/issues/12
 *
 * This function is not used if local ghost translation is active.
 *
 * @param mpiGrid DCCRG grid object
 * @param dimension Spatial dimension
 * @param direction Direction of communication (+ or -)
 * @param popId Particle population ID
 */
void update_remote_mapping_contribution_amr(
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   const uint dimension,
   int direction,
   const uint popID) {

   const vector<CellID> local_cells = mpiGrid.get_local_cells_on_process_boundary(VLASOV_SOLVER_NEIGHBORHOOD_ID);
   const vector<CellID> remote_cells = mpiGrid.get_remote_cells_on_process_boundary(VLASOV_SOLVER_NEIGHBORHOOD_ID);
   vector<CellID> receive_cells;
   set<CellID> send_cells;

   vector<CellID> receive_origin_cells;
   vector<uint> receive_origin_index;

   int neighborhood = 0;

   //normalize and set neighborhoods
   if(direction > 0) {
      direction = 1;
      switch (dimension) {
      case 0:
         neighborhood = SHIFT_P_X_NEIGHBORHOOD_ID;
         break;
      case 1:
         neighborhood = SHIFT_P_Y_NEIGHBORHOOD_ID;
         break;
      case 2:
         neighborhood = SHIFT_P_Z_NEIGHBORHOOD_ID;
         break;
      default:
         cerr << __FILE__ << ":"<< __LINE__ << " Wrong dimension, abort"<<endl;
         abort();
      }
   }
   if(direction < 0) {
      direction = -1;
      switch (dimension) {
      case 0:
         neighborhood = SHIFT_M_X_NEIGHBORHOOD_ID;
         break;
      case 1:
         neighborhood = SHIFT_M_Y_NEIGHBORHOOD_ID;
         break;
      case 2:
         neighborhood = SHIFT_M_Z_NEIGHBORHOOD_ID;
         break;
      default:
         cerr << __FILE__ << ":"<< __LINE__ << " Wrong dimension, abort"<<endl;
         abort();
      }
   }

   // MPI_Barrier(MPI_COMM_WORLD);
   // cout << "begin update_remote_mapping_contribution_amr, dimension = " << dimension << ", direction = " << direction << endl;
   // MPI_Barrier(MPI_COMM_WORLD);

   // Initialize remote cells
   for (auto rc : remote_cells) {
      SpatialCell *ccell = mpiGrid[rc];
      // Initialize number of blocks to 0 and block data to a default value.
      // We need the default for 1 to 1 communications
      if(ccell) {
         for (uint i = 0; i < MAX_NEIGHBORS_PER_DIM; ++i) {
            ccell->neighbor_block_data[i] = ccell->get_data(popID);
            ccell->neighbor_number_of_blocks[i] = 0;
         }
      }
   }

   // Initialize local cells
   for (auto lc : local_cells) {
      SpatialCell *ccell = mpiGrid[lc];
      if(ccell) {
         // Initialize number of blocks to 0 and neighbor block data pointer to the local block data pointer
         for (uint i = 0; i < MAX_NEIGHBORS_PER_DIM; ++i) {
            ccell->neighbor_block_data[i] = ccell->get_data(popID);
            ccell->neighbor_number_of_blocks[i] = 0;
         }
      }
   }

   vector<Realf*> receiveBuffers;
   vector<Realf*> sendBuffers;

   for (auto c : local_cells) {

      SpatialCell *ccell = mpiGrid[c];

      if (!ccell) continue;

      vector<CellID> p_nbrs;
      vector<CellID> n_nbrs;

      for (const auto& [neighbor, dir] : mpiGrid.get_face_neighbors_of(c)) {
         if(dir == ((int)dimension + 1) * direction) {
            p_nbrs.push_back(neighbor);
         }

         if(dir == -1 * ((int)dimension + 1) * direction) {
            n_nbrs.push_back(neighbor);
         }
      }

      uint sendIndex = 0;
      uint recvIndex = 0;

      int mySiblingIndex = get_sibling_index(mpiGrid,c);

      // Set up sends if any neighbor cells in p_nbrs are non-local.
      if (!all_of(p_nbrs.begin(), p_nbrs.end(), [&mpiGrid](CellID i){return mpiGrid.is_local(i);})) {

         // ccell adds a neighbor_block_data block for each neighbor in the positive direction to its local data
         for (const auto& nbr : p_nbrs) {

            //Send data in nbr target array that we just mapped to, if
            // 1) it is a valid target,
            // 2) the source cell in center was translated,
            // 3) Cell is remote.
            if(nbr != INVALID_CELLID && do_translate_cell(ccell) && !mpiGrid.is_local(nbr)) {

               /*
                 Select the index to the neighbor_block_data and neighbor_number_of_blocks arrays
                 1) Ref_c == Ref_nbr == 0, index = 0
                 2) Ref_c == Ref_nbr != 0, index = c sibling index
                 3) Ref_c >  Ref_nbr     , index = c sibling index
                 4) Ref_c <  Ref_nbr     , index = nbr sibling index
                */

               if(mpiGrid.get_refinement_level(c) >= mpiGrid.get_refinement_level(nbr)) {
                  sendIndex = mySiblingIndex;
               } else {
                  sendIndex = get_sibling_index(mpiGrid,nbr);
               }

               SpatialCell *pcell = mpiGrid[nbr];

               // 4) it exists and is not a boundary cell,
               if(pcell && pcell->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {

                  ccell->neighbor_number_of_blocks.at(sendIndex) = pcell->get_number_of_velocity_blocks(popID);

                  if(send_cells.find(nbr) == send_cells.end()) {
                     // 5 We have not already sent data from this rank to this cell.

                     ccell->neighbor_block_data.at(sendIndex) = pcell->get_data(popID);
                     send_cells.insert(nbr);

                  } else {

                     // The receiving cell can't know which cell is sending the data from this rank.
                     // Therefore, we have to send 0's from other cells in the case where multiple cells
                     // from one rank are sending to the same remote cell so that all sent cells can be
                     // summed for the correct result.

                     ccell->neighbor_block_data.at(sendIndex) =
                        (Realf*) aligned_malloc(ccell->neighbor_number_of_blocks.at(sendIndex) * WID3 * sizeof(Realf), WID3);
                     sendBuffers.push_back(ccell->neighbor_block_data.at(sendIndex));
                     for (uint j = 0; j < ccell->neighbor_number_of_blocks.at(sendIndex) * WID3; ++j) {
                        ccell->neighbor_block_data.at(sendIndex)[j] = 0.0;

                     } // closes for(uint j = 0; j < ccell->neighbor_number_of_blocks.at(sendIndex) * WID3; ++j)

                  } // closes if(send_cells.find(nbr) == send_cells.end())

               } // closes if(pcell && pcell->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY)

            } // closes if(nbr != INVALID_CELLID && do_translate_cell(ccell) && !mpiGrid.is_local(nbr))

         } // closes for(uint i_nbr = 0; i_nbr < nbrs_to.size(); ++i_nbr)

      } // closes if(!all_of(nbrs_to.begin(), nbrs_to.end(),[&mpiGrid](CellID i){return mpiGrid.is_local(i);}))

      // Set up receives if any neighbor cells in n_nbrs are non-local.
      if (!all_of(n_nbrs.begin(), n_nbrs.end(), [&mpiGrid](CellID i){return mpiGrid.is_local(i);})) {

         // ccell adds a neighbor_block_data block for each neighbor in the positive direction to its local data
         for (const auto& nbr : n_nbrs) {

            if (nbr != INVALID_CELLID && !mpiGrid.is_local(nbr) &&
                ccell->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
               //Receive data that ncell mapped to this local cell data array,
               //if 1) ncell is a valid source cell, 2) center cell is to be updated (normal cell) 3) ncell is remote

               SpatialCell *ncell = mpiGrid[nbr];

               // Check for null pointer
               if(!ncell) {
                  continue;
               }

               /*
                 Select the index to the neighbor_block_data and neighbor_number_of_blocks arrays
                 1) Ref_nbr == Ref_c == 0, index = 0
                 2) Ref_nbr == Ref_c != 0, index = nbr sibling index
                 3) Ref_nbr >  Ref_c     , index = nbr sibling index
                 4) Ref_nbr <  Ref_c     , index = c   sibling index
                */

               if(mpiGrid.get_refinement_level(nbr) >= mpiGrid.get_refinement_level(c)) {

                  // Allocate memory for one sibling at recvIndex.

                  recvIndex = get_sibling_index(mpiGrid,nbr);

                  ncell->neighbor_number_of_blocks.at(recvIndex) = ccell->get_number_of_velocity_blocks(popID);
                  ncell->neighbor_block_data.at(recvIndex) =
                     (Realf*) aligned_malloc(ncell->neighbor_number_of_blocks.at(recvIndex) * WID3 * sizeof(Realf), WID3);
                  receiveBuffers.push_back(ncell->neighbor_block_data.at(recvIndex));

               } else {

                  recvIndex = mySiblingIndex;

                  // std::array<uint64_t, 8> siblingarr = mpiGrid.mapping.get_all_children(mpiGrid.mapping.get_parent(c));
                  // vector<CellID> mySiblings(siblingarr.begin(), siblingarr.end());
                  auto mySiblings = mpiGrid.get_all_children(mpiGrid.get_parent(c));
                  auto myIndices = mpiGrid.mapping.get_indices(c);

                  // Allocate memory for each sibling to receive all the data sent by coarser ncell.
                  // only allocate blocks for face neighbors.
                  for (uint i_sib = 0; i_sib < MAX_NEIGHBORS_PER_DIM; ++i_sib) {

                     auto sibling = mySiblings.at(i_sib);
                     auto sibIndices = mpiGrid.mapping.get_indices(sibling);
                     auto* scell = mpiGrid[sibling];


                     // Only allocate siblings that are remote face neighbors to ncell
                     // Also take care to have these consistent with the sending process neighbor checks!
                     if(sibling != INVALID_CELLID
                        && scell
                        && mpiGrid.get_process(sibling) != mpiGrid.get_process(nbr)
                        && myIndices.at(dimension) == sibIndices.at(dimension)
                        && ncell->neighbor_number_of_blocks.at(i_sib) != scell->get_number_of_velocity_blocks(popID)
                        && scell->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {


                        ncell->neighbor_number_of_blocks.at(i_sib) = scell->get_number_of_velocity_blocks(popID);
                        ncell->neighbor_block_data.at(i_sib) =
                           (Realf*) aligned_malloc(ncell->neighbor_number_of_blocks.at(i_sib) * WID3 * sizeof(Realf), WID3);
                        receiveBuffers.push_back(ncell->neighbor_block_data.at(i_sib));
                     }
                  }
               }

               receive_cells.push_back(c);
               receive_origin_cells.push_back(nbr);
               receive_origin_index.push_back(recvIndex);

            } // closes (nbr != INVALID_CELLID && !mpiGrid.is_local(nbr) && ...)

         } // closes for(uint i_nbr = 0; i_nbr < nbrs_of.size(); ++i_nbr)

      } // closes if(!all_of(nbrs_of.begin(), nbrs_of.end(),[&mpiGrid](CellID i){return mpiGrid.is_local(i);}))

   } // closes for (auto c : local_cells) {

   MPI_Barrier(MPI_COMM_WORLD);

   // Do communication
   SpatialCell::setCommunicatedSpecies(popID);
   SpatialCell::set_mpi_transfer_type(Transfer::NEIGHBOR_VEL_BLOCK_DATA);
   mpiGrid.update_copies_of_remote_neighbors(neighborhood);

   MPI_Barrier(MPI_COMM_WORLD);

   // Reduce data: sum received data in the data array to
   // the target grid in the temporary block container
   //#pragma omp parallel
   {
      for (size_t c = 0; c < receive_cells.size(); ++c) {
         SpatialCell* receive_cell = mpiGrid[receive_cells[c]];
         SpatialCell* origin_cell = mpiGrid[receive_origin_cells[c]];

         if(!receive_cell || !origin_cell) {
            continue;
         }

         Realf *blockData = receive_cell->get_data(popID);
         Realf *neighborData = origin_cell->neighbor_block_data[receive_origin_index[c]];

         //#pragma omp for
         for(uint vCell = 0; vCell < WID3 * receive_cell->get_number_of_velocity_blocks(popID); ++vCell) {
            blockData[vCell] += neighborData[vCell];
         }
      }

      // send cell data is set to zero. This is to avoid double copy if
      // one cell is the neighbor on bot + and - side to the same process
      for (auto c : send_cells) {
         SpatialCell* spatial_cell = mpiGrid[c];
         Realf * blockData = spatial_cell->get_data(popID);
         //#pragma omp for nowait
         for(unsigned int vCell = 0; vCell < WID3 * spatial_cell->get_number_of_velocity_blocks(popID); ++vCell) {
            // copy received target data to temporary array where target data is stored.
            blockData[vCell] = 0;
         }
      }
   }

   for (auto p : receiveBuffers) {
      aligned_free(p);
   }
   for (auto p : sendBuffers) {
      aligned_free(p);
   }

   // MPI_Barrier(MPI_COMM_WORLD);
   // cout << "end update_remote_mapping_contribution_amr, dimension = " << dimension << ", direction = " << direction << endl;
   // MPI_Barrier(MPI_COMM_WORLD);

}
