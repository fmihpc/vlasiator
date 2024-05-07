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

/** GPU kernel for incrementing one set of blocks with another
    used for accruing remote neighbor contributions
*/
__global__ static void remote_increment_kernel (
   Realf* blockData,
   Realf* neighborData,
   vmesh::LocalID nBlocks
   ) {
   const int gpuBlocks = gridDim.x;
   const int blocki = blockIdx.x; // ==LID
   const int i = threadIdx.x;
   const int j = threadIdx.y;
   const int k = threadIdx.z;
   const uint ti = k*WID2 + j*WID + i;
   // Increment value
   atomicAdd(&blockData[blocki * WID3 + ti],neighborData[blocki * WID3 + ti]);
}

/* This function communicates the mapping on process boundaries, and then updates the data to their correct values.
 * When sending data between neighbors of different refinement levels, special care has to be taken to ensure that
 * The sending and receiving ranks allocate the correct size arrays for neighbor_block_data.
 * This is partially due to DCCRG defining neighborhood size relative to the host cell. For details, see
 * https://github.com/fmihpc/dccrg/issues/12
 *
 * @param mpiGrid DCCRG grid object
 * @param dimension Spatial dimension
 * @param direction Direction of communication (+ or -)
 * @param popId Particle population ID
 */
void gpu_update_remote_mapping_contribution_amr(
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   const uint dimension,
   int direction,
   const uint popID) {

   // GPUTODO: First attempts at using managed memory for remote neighbours
   // Should move to re-using managed memory buffers and ensuring size is suitable?
   // If that path is taken, it should also check for any local cells *not* on process
   // boundary and free the buffers from those cells.
   int device = gpu_getDevice();

   int neighborhood = 0;
   int both_neighborhood = 0;

   //normalize and set neighborhoods
   if(direction > 0) {
      direction = 1;
      switch (dimension) {
      case 0:
         neighborhood = SHIFT_P_X_NEIGHBORHOOD_ID;
         both_neighborhood = VLASOV_SOLVER_TARGET_X_NEIGHBORHOOD_ID;
         break;
      case 1:
         neighborhood = SHIFT_P_Y_NEIGHBORHOOD_ID;
         both_neighborhood = VLASOV_SOLVER_TARGET_Y_NEIGHBORHOOD_ID;
         break;
      case 2:
         neighborhood = SHIFT_P_Z_NEIGHBORHOOD_ID;
         both_neighborhood = VLASOV_SOLVER_TARGET_Z_NEIGHBORHOOD_ID;
         break;
      }
   }
   if(direction < 0) {
      direction = -1;
      switch (dimension) {
      case 0:
         neighborhood = SHIFT_M_X_NEIGHBORHOOD_ID;
         both_neighborhood = VLASOV_SOLVER_TARGET_X_NEIGHBORHOOD_ID;
         break;
      case 1:
         neighborhood = SHIFT_M_Y_NEIGHBORHOOD_ID;
         both_neighborhood = VLASOV_SOLVER_TARGET_Y_NEIGHBORHOOD_ID;
         break;
      case 2:
         neighborhood = SHIFT_M_Z_NEIGHBORHOOD_ID;
         both_neighborhood = VLASOV_SOLVER_TARGET_Z_NEIGHBORHOOD_ID;
         break;
      }
   }

   const vector<CellID>& local_cells_pb = mpiGrid.get_local_cells_on_process_boundary(both_neighborhood);
   const vector<CellID> remote_cells_pb = mpiGrid.get_remote_cells_on_process_boundary(both_neighborhood);

   // Fast return if no cells to process
   if (local_cells_pb.size() == remote_cells_pb.size() == 0) {
      return;
   }

   const vector<CellID>& local_cells = getLocalCells();
   const vector<CellID> remote_cells = mpiGrid.get_remote_cells_on_process_boundary(VLASOV_SOLVER_NEIGHBORHOOD_ID);

   vector<CellID> receive_cells;
   set<CellID> send_cells;
   vector<CellID> receive_origin_cells;
   vector<uint> receive_origin_index;

   phiprof::Timer updateRemoteTimerPre {"trans-amr-remotes-setup-getcells"};
   // Initialize remote cells
   #pragma omp parallel for
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
   #pragma omp parallel for
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
   updateRemoteTimerPre.stop();

   vector<Realf*> receiveBuffers;
   vector<Realf*> sendBuffers;

   phiprof::Timer updateRemoteTimer0 {"trans-amr-remotes-setup-localcells"};
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
         phiprof::Timer updateRemoteTimer1 {"trans-amr-remotes-setup-sends"};
         // ccell adds a neighbor_block_data block for each neighbor in the positive direction to its local data
         for (const auto nbr : p_nbrs) {
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

                     if (ccell->neighbor_number_of_blocks.at(sendIndex) == 0) {
                        ccell->neighbor_block_data.at(sendIndex) = 0;
                        sendBuffers.push_back(0);
                        continue;
                     }
                     // GPUTODO: This is now unified memory. With GPU-aware MPI it could be on-device.
                     CHK_ERR( gpuMallocManaged((void**)&ccell->neighbor_block_data.at(sendIndex), ccell->neighbor_number_of_blocks.at(sendIndex) * WID3 * sizeof(Realf)) );
                     // CHK_ERR( gpuMemPrefetchAsync(ccell->neighbor_block_data.at(sendIndex),ccell->neighbor_number_of_blocks.at(sendIndex) * WID3 * sizeof(Realf),device,0) );
                     CHK_ERR( gpuMemset(ccell->neighbor_block_data.at(sendIndex), 0, ccell->neighbor_number_of_blocks.at(sendIndex) * WID3 * sizeof(Realf)) );
                     sendBuffers.push_back(ccell->neighbor_block_data.at(sendIndex));
                  } // closes if(send_cells.find(nbr) == send_cells.end())
               } // closes if(pcell && pcell->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY)
            } // closes if(nbr != INVALID_CELLID && do_translate_cell(ccell) && !mpiGrid.is_local(nbr))
         } // closes for(uint i_nbr = 0; i_nbr < nbrs_to.size(); ++i_nbr)
      } // closes if(!all_of(nbrs_to.begin(), nbrs_to.end(),[&mpiGrid](CellID i){return mpiGrid.is_local(i);}))

      // Set up receives if any neighbor cells in n_nbrs are non-local.
      if (!all_of(n_nbrs.begin(), n_nbrs.end(), [&mpiGrid](CellID i){return mpiGrid.is_local(i);})) {
         phiprof::Timer updateRemoteTimer2 {"trans-amr-remotes-setup-receives"};
         // ccell adds a neighbor_block_data block for each neighbor in the positive direction to its local data
         for (const auto nbr : n_nbrs) {
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
                  if (ncell->neighbor_number_of_blocks.at(recvIndex) == 0) {
                     receiveBuffers.push_back(0);
                  } else {
                     // GPUTODO: This is now unified memory. With GPU-aware MPI it could be on-device.
                     CHK_ERR( gpuMallocManaged((void**)&ncell->neighbor_block_data.at(recvIndex), ncell->neighbor_number_of_blocks.at(recvIndex) * WID3 * sizeof(Realf)) );
                     CHK_ERR( gpuMemPrefetchAsync(ncell->neighbor_block_data.at(recvIndex), ncell->neighbor_number_of_blocks.at(recvIndex) * WID3 * sizeof(Realf), device,0) );
                     receiveBuffers.push_back(ncell->neighbor_block_data.at(recvIndex));
                  }
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
                        if (ncell->neighbor_number_of_blocks.at(i_sib) == 0) {
                           receiveBuffers.push_back(0);
                        } else {
                           // GPUTODO: This is now unified memory. With GPU-aware MPI it could be on-device.
                           CHK_ERR( gpuMallocManaged((void**)&ncell->neighbor_block_data.at(i_sib), ncell->neighbor_number_of_blocks.at(i_sib) * WID3 * sizeof(Realf)) );
                           CHK_ERR( gpuMemPrefetchAsync(ncell->neighbor_block_data.at(i_sib), ncell->neighbor_number_of_blocks.at(i_sib) * WID3 * sizeof(Realf), device,0) );
                           receiveBuffers.push_back(ncell->neighbor_block_data.at(i_sib));
                        }
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
   updateRemoteTimer0.stop();

   MPI_Barrier(MPI_COMM_WORLD);
   phiprof::Timer updateRemoteTimer3 {"trans-amr-remotes-MPI"};

   // Do communication
   SpatialCell::setCommunicatedSpecies(popID);
   SpatialCell::set_mpi_transfer_type(Transfer::NEIGHBOR_VEL_BLOCK_DATA);
   mpiGrid.update_copies_of_remote_neighbors(neighborhood);
   updateRemoteTimer3.stop();

   MPI_Barrier(MPI_COMM_WORLD);

   // Reduce data: sum received data in the data array to
   // the target grid in the temporary block container
   // #pragma omp parallel
   if (receive_cells.size() != 0) {
      phiprof::Timer updateRemoteTimerIncrement {"trans-amr-remotes-increment"};
      for (size_t c = 0; c < receive_cells.size(); ++c) {
         SpatialCell* receive_cell = mpiGrid[receive_cells[c]];
         SpatialCell* origin_cell = mpiGrid[receive_origin_cells[c]];
         if(!receive_cell || !origin_cell) {
            continue;
         }

         Realf *blockData = receive_cell->get_data(popID);
         Realf *neighborData = origin_cell->neighbor_block_data[receive_origin_index[c]];
         vmesh::LocalID nBlocks = receive_cell->get_number_of_velocity_blocks(popID);
         const uint maxThreads = gpu_getMaxThreads();
         // Increment needs to be parallel-safe, so use modulo of cellid as stream number
         gpuStream_t cellStream = gpuStreamList[receive_cells[c] % maxThreads];
         if (nBlocks>0) {
            dim3 block(WID,WID,WID);
            remote_increment_kernel<<<nBlocks, block, 0, cellStream>>> (
               blockData,
               neighborData,
               nBlocks
               );
            CHK_ERR( gpuPeekAtLastError() );
            //CHK_ERR( gpuStreamSynchronize(cellStream) );
         }
      }
      // Since all increment kernel streams were launched from outside an openmp region, use device sync here.
      CHK_ERR( gpuDeviceSynchronize() );

      // send cell data is set to zero. This is to avoid double copy if
      // one cell is the neighbor on bot + and - side to the same process
      vector<CellID> send_cells_vector(send_cells.begin(), send_cells.end());
      #pragma omp parallel for
      for (uint c = 0; c < send_cells_vector.size(); c++) {
         SpatialCell* send_cell = mpiGrid[send_cells_vector[c]];
         gpuStream_t stream = gpu_getStream();
         Realf* blockData = send_cell->get_data(popID);
         CHK_ERR( gpuMemsetAsync(blockData, 0, WID3*send_cell->get_number_of_velocity_blocks(popID)*sizeof(Realf),stream) );
      }
   }

   phiprof::Timer updateRemoteTimerFree {"trans-amr-remotes-free"};
   for (auto p : receiveBuffers) {
      CHK_ERR( gpuFree(p) );
   }
   for (auto p : sendBuffers) {
      CHK_ERR( gpuFree(p) );
   }
   updateRemoteTimerFree.stop();
}
