/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
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

#include "common.h"
#include <array>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <tuple>
#include <stdint.h>
#include <bitset>

#ifdef _OPENMP
   #include <omp.h>
#endif

#include <zoltan.h>
#include <dccrg.hpp>
#include <phiprof.hpp>

#include "vlasovmover.h"
#include "../grid.h"
#include "../definitions.h"
#include "../object_wrapper.h"
#include "../mpiconversion.h"

#include "arch_moments.h"

#include "cpu_trans_pencils.hpp"
#include "cpu_acc_transform.hpp" // for updateAccelerationMaxdt
#ifdef USE_GPU
#include "gpu_moments.h"
#include "gpu_acc_semilag.hpp"
#include "gpu_trans_map_amr.hpp"
#else
#include "cpu_acc_semilag.hpp"
#include "cpu_trans_pencils.hpp"
#include "cpu_trans_map_amr.hpp"
#include "cpu_acc_transform.hpp"

#include "spline.h"
#endif

using namespace spatial_cell;

/** Propagates the distribution function in spatial space.

    Based on SLICE-3D algorithm: Zerroukat, M., and T. Allen. "A
    three-dimensional monotone and conservative semi-Lagrangian scheme
    (SLICE-3D) for transport problems." Quarterly Journal of the Royal
    Meteorological Society 138.667 (2012): 1640-1651.

 */
void calculateSpatialTranslation(
        dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
        const vector<CellID>& local_propagated_cells,
        const vector<CellID>& remoteTargetCellsx,
        const vector<CellID>& remoteTargetCellsy,
        const vector<CellID>& remoteTargetCellsz,
        vector<uint>& nPencils,
        const Realf dt,
        const uint popID,
        Real &time
) {

    double t1;

    int myRank;
    MPI_Comm_rank(MPI_COMM_WORLD,&myRank);

   phiprof::Timer btzTimer {"barrier-trans-pre-z", {"Barriers","MPI"}};
   MPI_Barrier(MPI_COMM_WORLD);
   btzTimer.stop();

   // ------------- SLICE - map dist function in Z --------------- //
   if(P::zcells_ini > 1){

      phiprof::Timer transTimer {"transfer-stencil-data-z", {"MPI"}};
      //updateRemoteVelocityBlockLists(mpiGrid,popID,VLASOV_SOLVER_Z);
      SpatialCell::set_mpi_transfer_type(Transfer::VEL_BLOCK_DATA,false);
      mpiGrid.update_copies_of_remote_neighbors(Neighborhoods::VLASOV_SOLVER_Z);
      transTimer.stop();

      // bt=phiprof::initializeTimer("barrier-trans-pre-trans_map_1d-z","Barriers","MPI");
      // phiprof::start(bt);
      // MPI_Barrier(MPI_COMM_WORLD);
      // phiprof::stop(bt);

      t1 = MPI_Wtime();
      phiprof::Timer computeTimer {"compute-mapping-z"};
      trans_map_1d_amr(mpiGrid, DimensionPencils, local_propagated_cells, remoteTargetCellsz, nPencils, 2, dt, 0, popID); // map along z//
      computeTimer.stop();
      time += MPI_Wtime() - t1;

      phiprof::Timer btTimer {"barrier-trans-pre-update_remote-z", {"Barriers","MPI"}};
      MPI_Barrier(MPI_COMM_WORLD);
      btTimer.stop();

      phiprof::Timer updateRemoteTimer {"update_remote-z", {"MPI"}};
      update_remote_mapping_contribution_amr(mpiGrid, 2,+1,popID);
      update_remote_mapping_contribution_amr(mpiGrid, 2,-1,popID);
      updateRemoteTimer.stop();

   }

   phiprof::Timer btxTimer {"barrier-trans-pre-x", {"Barriers","MPI"}};
   MPI_Barrier(MPI_COMM_WORLD);
   btxTimer.stop();

   // ------------- SLICE - map dist function in X --------------- //
   if(P::xcells_ini > 1){

      phiprof::Timer transTimer {"transfer-stencil-data-x", {"MPI"}};
      //updateRemoteVelocityBlockLists(mpiGrid,popID,VLASOV_SOLVER_X);
      SpatialCell::set_mpi_transfer_type(Transfer::VEL_BLOCK_DATA,false);
      mpiGrid.update_copies_of_remote_neighbors(Neighborhoods::VLASOV_SOLVER_X);
      transTimer.stop();

      // bt=phiprof::initializeTimer("barrier-trans-pre-trans_map_1d-x","Barriers","MPI");
      // phiprof::start(bt);
      // MPI_Barrier(MPI_COMM_WORLD);
      // phiprof::stop(bt);

      t1 = MPI_Wtime();
      phiprof::Timer computeTimer {"compute-mapping-x"};
      trans_map_1d_amr(mpiGrid, DimensionPencils, local_propagated_cells, remoteTargetCellsx, nPencils, 0, dt, 0, popID); // map along x//
      computeTimer.stop();
      time += MPI_Wtime() - t1;

      phiprof::Timer btTimer {"barrier-trans-pre-update_remote-x", {"Barriers","MPI"}};
      MPI_Barrier(MPI_COMM_WORLD);
      btTimer.stop();

      phiprof::Timer updateRemoteTimer {"update_remote-x", {"MPI"}};
      update_remote_mapping_contribution_amr(mpiGrid, 0,+1,popID);
      update_remote_mapping_contribution_amr(mpiGrid, 0,-1,popID);
      updateRemoteTimer.stop();

   }

   phiprof::Timer btyTimer {"barrier-trans-pre-y", {"Barriers","MPI"}};
   MPI_Barrier(MPI_COMM_WORLD);
   btyTimer.stop();

   // ------------- SLICE - map dist function in Y --------------- //
   if(P::ycells_ini > 1) {

      phiprof::Timer transTimer {"transfer-stencil-data-y", {"MPI"}};
      //updateRemoteVelocityBlockLists(mpiGrid,popID,VLASOV_SOLVER_Y);
      SpatialCell::set_mpi_transfer_type(Transfer::VEL_BLOCK_DATA,false);
      mpiGrid.update_copies_of_remote_neighbors(Neighborhoods::VLASOV_SOLVER_Y);
      transTimer.stop();

      // bt=phiprof::initializeTimer("barrier-trans-pre-trans_map_1d-y","Barriers","MPI");
      // phiprof::start(bt);
      // MPI_Barrier(MPI_COMM_WORLD);
      // phiprof::stop(bt);

      t1 = MPI_Wtime();
      phiprof::Timer computeTimer {"compute-mapping-y"};
      trans_map_1d_amr(mpiGrid, DimensionPencils, local_propagated_cells, remoteTargetCellsy, nPencils, 1, dt, 0, popID); // map along y//
      computeTimer.stop();
      time += MPI_Wtime() - t1;

      phiprof::Timer btTimer {"barrier-trans-pre-update_remote-y", {"Barriers","MPI"}};
      MPI_Barrier(MPI_COMM_WORLD);
      btTimer.stop();

      phiprof::Timer updateRemoteTimer {"update_remote-y", {"MPI"}};
      update_remote_mapping_contribution_amr(mpiGrid, 1,+1,popID);
      update_remote_mapping_contribution_amr(mpiGrid, 1,-1,popID);
      updateRemoteTimer.stop();

   }

   for(CellID c : local_propagated_cells)
   {
      // if (c == 16) std::cout << c << " at TIME_R " << mpiGrid[c]->parameters[CellParams::TIME_R] << " + " << dt <<"\n";
      if(mpiGrid[c]->get_timeclass_turn_r())
         mpiGrid[c]->parameters[CellParams::TIME_R] += dt;
   }
   phiprof::Timer btpostimer {"barrier-trans-post-trans",{"Barriers","MPI"}};
   MPI_Barrier(MPI_COMM_WORLD);
   btpostimer.stop();

   // MPI_Barrier(MPI_COMM_WORLD);
   // bailout(true, "", __FILE__, __LINE__);
}

/** Communicate the distribution function from necessary remote neighbors
 *  before ghost translation.
    Now does all required calculations on ghost cells,
    coalescing all interim MPI communication into one call..
 */

void communicatePreSpatialGhostTranslation(
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   const vector<CellID>& local_propagated_cells,
   vector<uint>& nPencils,
   const creal dt,
   const uint popID,
   Real &time,
   const uint tictoc,
   int tc
   ) {

   // Ghost translation, need all cell information, not just for a single direction.
   // No need for remote target cells; pass a dummy list.
   const vector<CellID> dummy_cells;
   uint neighborhood;

   // Select the pencil set and the required neighborhood based on the tictoc value
   std::array<setOfPencils,3>* pencilSet;
   if (P::currentMaxTimeclass == 0) {
      neighborhood = Neighborhoods::VLASOV_SOLVER_GHOST;
   }
   else if (tictoc == Timeclasses::TIC) {
      neighborhood = Neighborhoods::VLASOV_SOLVER_TIMEGHOST_OUTER_HALO;
   }
   else if (tictoc == Timeclasses::TOC) {
      neighborhood = Neighborhoods::VLASOV_SOLVER_GHOST;
   }
   else {
      std::cerr << __FILE__<<":"<<__LINE__<< " Unknown state: Current maxtimeclass=" << P::currentMaxTimeclass << ", tictoc = " << tictoc << std::endl;
      abort();
   }

   updateRemoteVelocityBlockLists(mpiGrid,popID,neighborhood, tc);
   // Need to re-do in case block lists of boundary cells change after
   // the block adjustment just after ACC.

   phiprof::Timer prepreBarrierTimer {"MPI barrier-pre-trans-comm"};
   MPI_Barrier(MPI_COMM_WORLD);
   prepreBarrierTimer.stop();

   phiprof::Timer transferTimer {"transfer-stencil-data-all",{"MPI"}};
   SpatialCell::set_mpi_transfer_type(Transfer::VEL_BLOCK_DATA,false);
   mpiGrid.update_copies_of_remote_neighbors(neighborhood);
   transferTimer.stop();

   phiprof::Timer preBarrierTimer {"MPI barrier-post-comm-pre-trans"};
   MPI_Barrier(MPI_COMM_WORLD);
   preBarrierTimer.stop();
}


/** Communicate the distribution function from necessary remote neighbors
 *  before ghost translation.
    Now does all required calculations on ghost cells,
    coalescing all interim MPI communication into one call..
 */

void communicatePreSpatialGhostTranslationCoalesced(
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   std::vector<std::pair<uint, int>>& population_indexes,
   const uint tictoc,
   int tc
   ) {

   // Ghost translation, need all cell information, not just for a single direction.
   // No need for remote target cells; pass a dummy list.
   const vector<CellID> dummy_cells;
   uint neighborhood;

   // Select the pencil set and the required neighborhood based on the tictoc value
   std::array<setOfPencils,3>* pencilSet;
   if (P::currentMaxTimeclass == 0) {
      neighborhood = Neighborhoods::VLASOV_SOLVER_GHOST;
   }
   else if (tictoc == Timeclasses::TIC) {
      neighborhood = Neighborhoods::VLASOV_SOLVER_TIMEGHOST_OUTER_HALO;
   }
   else if (tictoc == Timeclasses::TOC) {
      neighborhood = Neighborhoods::VLASOV_SOLVER_TIMEGHOST_EXACT_HALO;
   }
   else {
      std::cerr << __FILE__<<":"<<__LINE__<< " Unknown state: Current maxtimeclass=" << P::currentMaxTimeclass << ", tictoc = " << tictoc << std::endl;
      abort();
   }

   // updateRemoteVelocityBlockListsCoalesced(mpiGrid,population_indexes,neighborhood);
   P::coalesceGhostComms = false;
   for (auto& index : population_indexes) {
      updateRemoteVelocityBlockLists(mpiGrid, index.first, neighborhood, index.second);
   }
   P::coalesceGhostComms = true;
   // std::cerr << __FILE__<<":"<<__LINE__<< " Done updateing remote velocity block lists" << std::endl;
   // Need to re-do in case block lists of boundary cells change after
   // the block adjustment just after ACC.

   phiprof::Timer prepreBarrierTimer {"MPI barrier-pre-trans-comm"};
   MPI_Barrier(MPI_COMM_WORLD);
   prepreBarrierTimer.stop();
   SpatialCell::setCommunicatedSpecies(population_indexes);

   phiprof::Timer transferTimer {"transfer-stencil-data-all",{"MPI"}};
   SpatialCell::set_mpi_transfer_type(Transfer::VEL_BLOCK_DATA,false);
   mpiGrid.update_copies_of_remote_neighbors(neighborhood);
   transferTimer.stop();

   //phiprof::Timer preBarrierTimer {"MPI barrier-post-comm-pre-trans"};
   //MPI_Barrier(MPI_COMM_WORLD);
   //preBarrierTimer.stop();
}



/** Propagates the distribution function in spatial space.
    Now does all required calculations on ghost cells,
    coalescing all interim MPI communication into one call..

    Based on SLICE-3D algorithm: Zerroukat, M., and T. Allen. "A
    three-dimensional monotone and conservative semi-Lagrangian scheme
    (SLICE-3D) for transport problems." Quarterly Journal of the Royal
    Meteorological Society 138.667 (2012): 1640-1651.
 */

void calculateSpatialGhostTranslation(
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   const vector<CellID>& local_propagated_cells,
   vector<uint>& nPencils,
   const creal dt,
   const uint popID,
   Real &time,
   const uint tictoc,
   int tc
   ) {

   // Ghost translation, need all cell information, not just for a single direction.
   // No need for remote target cells; pass a dummy list.
   const vector<CellID> dummy_cells;
   
   //uint neighborhood;

   // Select the pencil set and the required neighborhood based on the tictoc value
   std::array<setOfPencils,3>* pencilSet;
   if (P::currentMaxTimeclass == 0) {
     // neighborhood = Neighborhoods::VLASOV_SOLVER_GHOST;
      pencilSet = &DimensionPencils;
      // std::cerr << __FILE__<<":"<<__LINE__<< " tictoc = " << tictoc<< ", pencilSet->Nx = " << (*pencilSet)[0].N << ", pencilSet->Ny = " << (*pencilSet)[1].N << ", pencilSet->Nz = " << (*pencilSet)[2].N << std::endl;
   }
   else if (tictoc == Timeclasses::TIC) {
      //neighborhood = Neighborhoods::VLASOV_SOLVER_TIMEGHOST_OUTER_HALO;
      pencilSet = &DimensionPencils;
      // std::cerr << __FILE__<<":"<<__LINE__<< "tc " << tc << (tictoc==Timeclasses::TIC?" tic":" toc")<< ", pencilSet->Nx = " << (*pencilSet)[0].N << ", pencilSet->Ny = " << (*pencilSet)[1].N << ", pencilSet->Nz = " << (*pencilSet)[2].N << std::endl;
   }
   else if (tictoc == Timeclasses::TOC) {
      //neighborhood = Neighborhoods::VLASOV_SOLVER_TIMEGHOST_EXACT_HALO;
      pencilSet = &DimensionPencils_toc;
      // std::cerr << __FILE__<<":"<<__LINE__<< "tc " << tc<< (tictoc==Timeclasses::TIC?" tic":" toc")<< ", pencilSet->Nx = " << (*pencilSet)[0].N << ", pencilSet->Ny = " << (*pencilSet)[1].N << ", pencilSet->Nz = " << (*pencilSet)[2].N << std::endl;
   }
   else {
      std::cerr << __FILE__<<":"<<__LINE__<< " Unknown state: Current maxtimeclass=" << P::currentMaxTimeclass << ", tictoc = " << tictoc << std::endl;
      abort();
   }
/* --- this is now done ahead of this function
   updateRemoteVelocityBlockLists(mpiGrid,popID,neighborhood, tc);
   // Need to re-do in case block lists of boundary cells change after
   // the block adjustment just after ACC.

   phiprof::Timer prepreBarrierTimer {"MPI barrier-pre-trans-comm"};
   MPI_Barrier(MPI_COMM_WORLD);
   prepreBarrierTimer.stop();

   phiprof::Timer transferTimer {"transfer-stencil-data-all",{"MPI"}};
   SpatialCell::set_mpi_transfer_type(Transfer::VEL_BLOCK_DATA,false);
   mpiGrid.update_copies_of_remote_neighbors(neighborhood);
   transferTimer.stop();
*/
//   phiprof::Timer preBarrierTimer {"MPI barrier-pre-trans"}; // TODO is this even needed?
//   MPI_Barrier(MPI_COMM_WORLD);
//   preBarrierTimer.stop();

   // Create temporary copies of copy-timeghosts so they are not overwritten and invalidated by the mapping
   // on the fringes. This placement tries to minimize the temporary memory footprint.

   std::cerr << "GT, tictoc = " << tictoc << ", pop " << popID << ", tc " << tc << "\n";
   std::map<CellID, Population> copy_ghosts;
   if (P::currentMaxTimeclass > 0 && tictoc == Timeclasses::TIC){
      std::cerr << "Tic copy\n";
      phiprof::Timer copyghostOuter {"copy-ghost temporaries"};
      phiprof::Timer copyghostInner {"gather copies"};
      for (auto c : getLocalCells()) {
         // std::cerr << " check copy for cellid " << c << ": " << mpiGrid[c]->requested_timeclass_copy_ghosts.count(tc) << "\n"; 
         if(mpiGrid[c]->requested_timeclass_copy_ghosts.count(tc) > 0){
            copy_ghosts[c] = mpiGrid[c]->get_population(popID, tc);
            // std::cerr << &copy_ghosts[c] << " vs " << &mpiGrid[c]->get_population(popID, tc) << "\n";
         }
         mpiGrid[c]->get_population(popID, tc).state = mpiGrid[c]->get_population(popID, tc).state & ~(1u << LastPropagationState::TRANSLATED_X);
         mpiGrid[c]->get_population(popID, tc).state = mpiGrid[c]->get_population(popID, tc).state & ~(1u << LastPropagationState::TRANSLATED_Y);
         mpiGrid[c]->get_population(popID, tc).state = mpiGrid[c]->get_population(popID, tc).state & ~(1u << LastPropagationState::TRANSLATED_Z);

      }
   }
   setOfPencils* pencils;
   std::set<CellID> target_cells = {}; // diagnostics and T_R updating (T_R should update the cells in active_list)
   std::set<CellID> source_cells = {}; 

   
   //#warning TODO: Implement also 2D / non-AMR ghost translation?
   // ------------- SLICE - map dist function in Z --------------- //
   phiprof::Timer mappingZTimer {"compute-mapping-z"};
   trans_map_1d_amr(mpiGrid, *pencilSet, local_propagated_cells, dummy_cells, nPencils, 2, dt, tc, popID); // map along z//
   
   pencils = &pencilSet->at(2);
   uint pencili = 0;
   uint curend = pencils->idsStart[pencili]+pencils->lengthOfPencils[pencili];
   for (size_t i = 0; i < pencils->sumOfLengths; ++i){
      if (i >=curend){ // in the next pencil already
         pencili++;
         curend = pencils->idsStart[pencili]+pencils->lengthOfPencils[pencili];
      }
      int tc = pencils->timeclasses[pencili];
      if (pencils->targetRatios[i] > 0.0){
         mpiGrid[pencils->ids[i]]->get_population(popID, tc).state = mpiGrid[pencils->ids[i]]->get_population(popID, tc).state | (1u << LastPropagationState::TRANSLATED_Z);
         target_cells.insert(pencils->ids[i]);
      }
      source_cells.insert(pencils->ids[i]);
   }   
   
   mappingZTimer.stop();

   // ------------- SLICE - map dist function in X --------------- //
   phiprof::Timer mappingXTimer {"compute-mapping-x"};
   trans_map_1d_amr(mpiGrid, *pencilSet, local_propagated_cells, dummy_cells, nPencils, 0,dt, tc, popID); // map along x//
   pencils = &pencilSet->at(0);
   pencili = 0;
   curend = pencils->idsStart[pencili]+pencils->lengthOfPencils[pencili];
   for (size_t i = 0; i < pencils->sumOfLengths; ++i){
      if (i >=curend){ // in the next pencil already
         pencili++;
         curend = pencils->idsStart[pencili]+pencils->lengthOfPencils[pencili];
      }
      int tc = pencils->timeclasses[pencili];
      if (pencils->targetRatios[i] > 0.0){
         mpiGrid[pencils->ids[i]]->get_population(popID, tc).state = mpiGrid[pencils->ids[i]]->get_population(popID, tc).state | (1u << LastPropagationState::TRANSLATED_X);
         target_cells.insert(pencils->ids[i]);
      }
      source_cells.insert(pencils->ids[i]);
   }
   mappingXTimer.stop();

   // ------------- SLICE - map dist function in Y --------------- //
   phiprof::Timer mappingYTimer {"compute-mapping-y"};
   trans_map_1d_amr(mpiGrid, *pencilSet, local_propagated_cells, dummy_cells, nPencils, 1,dt, tc, popID); // map along y//
   pencils = &pencilSet->at(1);
   pencili = 0;
   curend = pencils->idsStart[pencili]+pencils->lengthOfPencils[pencili];
   for (size_t i = 0; i < pencils->sumOfLengths; ++i){
      if (i >=curend){ // in the next pencil already
         pencili++;
         curend = pencils->idsStart[pencili]+pencils->lengthOfPencils[pencili];
      }
      int tc = pencils->timeclasses[pencili];
      if (pencils->targetRatios[i] > 0.0){
         mpiGrid[pencils->ids[i]]->get_population(popID, tc).state = mpiGrid[pencils->ids[i]]->get_population(popID, tc).state | (1u << LastPropagationState::TRANSLATED_Y);
         target_cells.insert(pencils->ids[i]);
      }
      source_cells.insert(pencils->ids[i]);
   }
   mappingYTimer.stop();

//   phiprof::Timer postBarrierTimer {"MPI barrier-post-trans"};
//   MPI_Barrier(MPI_COMM_WORLD);
//   postBarrierTimer.stop();

   Real t0 = 1e18;
   for(CellID c : source_cells){
      if(tictoc == Timeclasses::TIC){
         if (t0 == 1e18 && do_translate_cell(mpiGrid[c]) && mpiGrid[c]->has_source_timeclass_tic(tc)){
            t0 = mpiGrid[c]->get_population(popID, tc).T_V;
         }
      }else{
         if (t0 == 1e18 && do_translate_cell(mpiGrid[c]) && mpiGrid[c]->has_source_timeclass_toc(tc)){
            t0 = mpiGrid[c]->get_population(popID, tc).T_V;
         }
      }
   }
   for(CellID c : source_cells){
      if (!do_translate_cell(mpiGrid[c])) continue;
      if (tictoc == Timeclasses::TIC && !mpiGrid[c]->has_source_timeclass_tic(tc)) continue;
      if (tictoc == Timeclasses::TOC && !mpiGrid[c]->has_source_timeclass_toc(tc)) continue;
      Population& pop = mpiGrid[c]->get_population(popID, tc);
      if(t0 != pop.T_V){
         std::cerr << *source_cells.begin() << " and " << c << " have differing source times: " << t0 << " and " << pop.T_V << ", respectively. Delta " << pop.T_V - t0 << std::endl;
      }
   }
   for(CellID c : target_cells)
   {
      Population& pop = mpiGrid[c]->get_population(popID, tc);
      if(pop.state & (1u << LastPropagationState::TRANSLATED_X)){
         pop.T_R += dt/3.0;
      }
      if(pop.state & (1u << LastPropagationState::TRANSLATED_Y)){
         pop.T_R += dt/3.0;
      }
      if(pop.state & (1u << LastPropagationState::TRANSLATED_Z)){
         pop.T_R += dt/3.0;
      }
      if((pop.state & (1u << LastPropagationState::TRANSLATED_X)) && 
         (pop.state & (1u << LastPropagationState::TRANSLATED_Y)) &&
         (pop.state & (1u << LastPropagationState::TRANSLATED_Z))){
         
         pop.state = pop.state | (1u << LastPropagationState::TRANSLATED);
         pop.state = pop.state & ~(1u << LastPropagationState::TRANSLATED_X);
         pop.state = pop.state & ~(1u << LastPropagationState::TRANSLATED_Y);
         pop.state = pop.state & ~(1u << LastPropagationState::TRANSLATED_Z);
         }
      else{
         string tictocstr = tictoc == Timeclasses::TIC ? "tic" : "toc";
         if(do_translate_cell(mpiGrid[c]) && (
            (tictoc == Timeclasses::TIC && mpiGrid[c]->has_active_timeclass_tic(tc)) ||
            (tictoc == Timeclasses::TOC && mpiGrid[c]->has_active_timeclass_toc(tc))
            )
         ){
            std::cerr << __FILE__ << ":" << __LINE__ << ": incomplete translation for " << mpiGrid[c]->parameters[CellParams::CELLID] << "t" << mpiGrid[c]->parameters[CellParams::TIMECLASS] << " t = " << P::t << " at timeclass " << tc << tictocstr << ", state = " << std::bitset<8>(pop.state) << ", tG flags " << mpiGrid[c]->requested_timeclass_ghosts.count(tc) << ", " << mpiGrid[c]->requested_timeclass_copy_ghosts.count(tc) << ", " << mpiGrid[c]->requested_timeclass_outer_ghosts.count(tc)<< std::endl;
         }
      }

   }

   if(P::currentMaxTimeclass > 0 && tictoc == Timeclasses::TIC ){
      phiprof::Timer copyghostOuter {"copy-ghost temporaries"};
      phiprof::Timer copyghostInner {"restore copies"};
      for (auto c : copy_ghosts) {
         mpiGrid[c.first]->get_population(popID, tc) = c.second;
         mpiGrid[c.first]->get_population(popID, tc).T_R += dt;
      }
   }
   copy_ghosts.clear();

   return;
}

/*!

  Propagates the distribution function in spatial space.

  Based on SLICE-3D algorithm: Zerroukat, M., and T. Allen. "A
  three‐dimensional monotone and conservative semi‐Lagrangian scheme
  (SLICE‐3D) for transport problems." Quarterly Journal of the Royal
  Meteorological Society 138.667 (2012): 1640-1651.

*/
void calculateSpatialTranslation(
        dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
        creal dt,
        const bool initializationOrLB) {
   typedef Parameters P;
   // std::cerr << std::scientific << "calculateSpatialTranslation at t="<<P::t << "\n";//", for dtfactor="<<dt<<"\n";
   phiprof::Timer semilagTimer {"semilag-trans"};

   //double t1 = MPI_Wtime();

   const vector<CellID>& localCells = getLocalCells();
   vector<CellID> remoteTargetCellsx;
   vector<CellID> remoteTargetCellsy;
   vector<CellID> remoteTargetCellsz;
   vector<CellID> local_propagated_cells;
   vector<vector<CellID>> tc_propagated_cells = vector<vector<CellID>>();
   vector<CellID> local_target_cells;
   vector<vector<CellID>> tc_target_cells = vector<vector<CellID>>();
   vector<set<CellID>> tc_propagated_cell_sets = vector<set<CellID>>();
   vector<set<CellID>> tc_target_cell_sets = vector<set<CellID>>();

   vector<vector<uint>> nPencils(P::currentMaxTimeclass+1); // one for each timeclass
   Real time=0.0;

   // If dt=0 we are either initializing or distribution functions are not translated.
   if (dt == 0.0 && initializationOrLB == true) {
      // if dt=0.0 and initialization == true, we are either initializing or load balancing.
      // hence we can just calculate the moments and return.
      calculateMoments_R(mpiGrid,localCells,true,true);
      return;
   }

   // TC propagation lists, TODO move out of here somewhere sensible and less often called
   if (P::currentMaxTimeclass > 0) {
      for (int tc = 0; tc <= P::currentMaxTimeclass; tc++)
      {
         tc_propagated_cell_sets.push_back(set<CellID>());
         for (auto c : localCells){
            if(mpiGrid[c]->parameters[CellParams::TIMECLASS] == tc){
               tc_propagated_cell_sets[tc].insert(c);
            }
         }
         set<CellID> exact_halo_cells;
         for (auto c : tc_propagated_cell_sets[tc]){
            auto neighbors = mpiGrid.get_neighbors_of(c, Neighborhoods::VLASOV_SOLVER_TIMEGHOST_EXACT_HALO);
            for (auto n : *neighbors){
               exact_halo_cells.insert(n.first); // NB insert_range in C++23
            }
         }
         for(auto c : exact_halo_cells) {
         //   tc_propagated_cell_sets[tc].insert(c);
         }

         tc_propagated_cells.push_back(vector<CellID>(tc_propagated_cell_sets[tc].begin(),tc_propagated_cell_sets[tc].end()));
         // std::cout << "initing up tc " << tc << " vectors \n";
         // tc_propagated_cells.push_back(vector<CellID>(ghostTranslate_source[tc].begin(),ghostTranslate_source[tc].end()));
         // tc_target_cells.push_back(vector<CellID>(ghostTranslate_active[tc].begin(),ghostTranslate_active[tc].end()));
         // tc_propagated_cell_sets.push_back(set<CellID>(tc_propagated_cells_set));
         // tc_target_cell_sets.push_back(set<CellID>(ghostTranslate_active[tc].begin(),ghostTranslate_active[tc].end()));
      }
   }

   phiprof::Timer computeTimer {"compute_cell_lists"};
   if (!P::vlasovSolverGhostTranslate) {
      remoteTargetCellsx = mpiGrid.get_remote_cells_on_process_boundary(Neighborhoods::VLASOV_SOLVER_TARGET_X);
      remoteTargetCellsy = mpiGrid.get_remote_cells_on_process_boundary(Neighborhoods::VLASOV_SOLVER_TARGET_Y);
      remoteTargetCellsz = mpiGrid.get_remote_cells_on_process_boundary(Neighborhoods::VLASOV_SOLVER_TARGET_Z);
   }

   // Figure out which spatial cells are translated,
   // result independent of particle species.
   // If performing ghost translation, this is used for LB.
   for (size_t c=0; c<localCells.size(); ++c) {
      if (do_translate_cell(mpiGrid[localCells[c]])) {
         local_propagated_cells.push_back(localCells[c]);
      }
   }


   // TC propagation lists, TODO move out of here somewhere sensible and less often called
   computeTimer.stop();
   if (P::currentMaxTimeclass > 0) {
      for (int tc = 0; tc <= P::currentMaxTimeclass; tc++)
      {
         // std::cout << "initing up tc " << tc << " vectors \n";
         tc_propagated_cells[tc] = vector<CellID>(tc_propagated_cell_sets[tc].begin(),tc_propagated_cell_sets[tc].end());
         // tc_target_cells[tc] = vector<CellID>(tc_target_cell_sets[tc].begin(),tc_target_cell_sets[tc].end());
      }
   } else {
      // If we are not using timeclasses, we can just use the local_propagated_cells vector
      // for all the calculations.
      
      tc_propagated_cells.push_back(local_propagated_cells);
   }
   computeTimer.stop();

   //
   if (dt == 0.0 && initializationOrLB == false) {
      // If we are not initializing, we are in the main loop, and we are not propagating.
      // Here we have to update _R_PREV moments, update _R moments on a timeclass basis, and return.

      for (int tc=0; tc <= P::currentMaxTimeclass; tc++) {
         int mod = 1 << (P::currentMaxTimeclass - tc);
         if ((P::fractionalTimestep % mod) == 0) {
            for (size_t c=0; c<tc_propagated_cells.at(tc).size(); ++c) {
               const CellID cellID = tc_propagated_cells.at(tc).at(c);
               SpatialCell* SC = mpiGrid[cellID];

               SC->parameters[CellParams::RHOM_R_PREV_PREV] = SC->parameters[CellParams::RHOM_R_PREV];
               SC->parameters[CellParams::VX_R_PREV_PREV] = SC->parameters[CellParams::VX_R_PREV];
               SC->parameters[CellParams::VY_R_PREV_PREV] = SC->parameters[CellParams::VY_R_PREV];
               SC->parameters[CellParams::VZ_R_PREV_PREV] = SC->parameters[CellParams::VZ_R_PREV];
               SC->parameters[CellParams::RHOQ_R_PREV_PREV] = SC->parameters[CellParams::RHOQ_R_PREV];
               SC->parameters[CellParams::P_11_R_PREV_PREV] = SC->parameters[CellParams::P_11_R_PREV];
               SC->parameters[CellParams::P_22_R_PREV_PREV] = SC->parameters[CellParams::P_22_R_PREV];
               SC->parameters[CellParams::P_33_R_PREV_PREV] = SC->parameters[CellParams::P_33_R_PREV];
               SC->parameters[CellParams::P_23_R_PREV_PREV] = SC->parameters[CellParams::P_23_R_PREV];
               SC->parameters[CellParams::P_13_R_PREV_PREV] = SC->parameters[CellParams::P_13_R_PREV];
               SC->parameters[CellParams::P_12_R_PREV_PREV] = SC->parameters[CellParams::P_12_R_PREV];

               SC->parameters[CellParams::RHOM_R_PREV] = SC->parameters[CellParams::RHOM_R];
               SC->parameters[CellParams::VX_R_PREV] = SC->parameters[CellParams::VX_R];
               SC->parameters[CellParams::VY_R_PREV] = SC->parameters[CellParams::VY_R];
               SC->parameters[CellParams::VZ_R_PREV] = SC->parameters[CellParams::VZ_R];
               SC->parameters[CellParams::RHOQ_R_PREV] = SC->parameters[CellParams::RHOQ_R];
               SC->parameters[CellParams::P_11_R_PREV] = SC->parameters[CellParams::P_11_R];
               SC->parameters[CellParams::P_22_R_PREV] = SC->parameters[CellParams::P_22_R];
               SC->parameters[CellParams::P_33_R_PREV] = SC->parameters[CellParams::P_33_R];
               SC->parameters[CellParams::P_23_R_PREV] = SC->parameters[CellParams::P_23_R];
               SC->parameters[CellParams::P_13_R_PREV] = SC->parameters[CellParams::P_13_R];
               SC->parameters[CellParams::P_12_R_PREV] = SC->parameters[CellParams::P_12_R];
            }
         }
      }

      //calculating moments based on timeclass
      for (int tc=0; tc <= P::currentMaxTimeclass; tc++) {
         int mod = 1 << (P::currentMaxTimeclass - tc);
         if ((P::fractionalTimestep % mod) == 0) {
            calculateMoments_R(mpiGrid,tc_propagated_cells.at(tc),true);
         }
      }
      return;
   }
   
   if (P::prepareForRebalance == true) {
      // One more element to count the sums
      for(int tc = 0; tc <= P::currentMaxTimeclass; tc++) {
         nPencils.at(tc).resize(tc_propagated_cells.at(tc).size()+1, 0);
      }
   }

   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
   P::coalesceGhostComms = true;
   if (P::vlasovSolverGhostTranslate){
      phiprof::Timer pretrans {"Pre-ghostTranslate comms"};
      std::vector<std::pair<uint, int>> population_indexes;
      for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
      // if ghost translating, do the comms ahead of translations, which are now local ops
         for(int tc = 0; tc <= P::currentMaxTimeclass; tc++){
            SpatialCell::setCommunicatedSpecies(popID,tc);
            int mod = 1 << (P::currentMaxTimeclass - tc);
            int mod2 = 2 << (P::currentMaxTimeclass - tc);
            if((P::fractionalTimestep % mod) == 0){
               // std::cout << "rank " << myRank << ": " << tc_propagated_cells[tc].size() << " cells: calculateSpatialTranslation tc " << tc << " by dt " << P::timeclassDt[tc] <<"\n";
               // Local translation without interim communication
               uint tictoc = Timeclasses::TIC;
               if (tc == 0) { // we are the coarsest timeclass, so we always have fine timeclass data to fetch and sync
                  tictoc = Timeclasses::TOC;
               }
               else { // Other timeclasses need to get time ghost halo layers from the coarser timeclasses
                     // to have the intermediate translation step (TIC)
                  if ((P::fractionalTimestep % mod2) == 0){
                     // Initial timeghost halo layer needs to have a larger source region
                     tictoc = Timeclasses::TIC;
                  }
                  else{
                     // Subsequent timeghost halo layers can do with a smaller layer, using the translation
                     // targets from the TIC
                     tictoc = Timeclasses::TOC;
                  }
               }
               string tictocstr = (tictoc == Timeclasses::TIC ? "tic" : "toc");
               //string profName = "pre-translate comm "+getObjectWrapper().particleSpecies[popID].name+" tc "+std::to_string(tc) + tictocstr;
               //phiprof::Timer timer {profName};
               population_indexes.push_back(std::make_pair(popID, tc));
/*               if (!P::coalesceGhostComms) {
                  communicatePreSpatialGhostTranslation(
                     mpiGrid,
                     tc_propagated_cells[tc], // Used for LB
                     nPencils,
                     P::timeclassDt[tc],
                     popID,
                     time,
                     tictoc,
                     tc
                     );
               }*/
            }
         }
      }
      pretrans.stop();
      phiprof::Timer ccomms {"Coalesced comms"};
      if (P::coalesceGhostComms) {
         communicatePreSpatialGhostTranslationCoalesced(
            mpiGrid,
            population_indexes,
            Timeclasses::TIC,
            0
            );
      }
      P::coalesceGhostComms = false;
      ccomms.stop();
   }


   // Translate all particle species
   for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
      string profName = "translate "+getObjectWrapper().particleSpecies[popID].name;
      phiprof::Timer timer {profName};
      //SpatialCell::setCommunicatedSpecies(popID);
      for(int tc = 0; tc <= P::currentMaxTimeclass; tc++){
         SpatialCell::setCommunicatedSpecies(popID,tc);
         int mod = 1 << (P::currentMaxTimeclass - tc);
         int mod2 = 2 << (P::currentMaxTimeclass - tc);
         if((P::fractionalTimestep % mod) == 0){
            // std::cout << "rank " << myRank << ": " << tc_propagated_cells[tc].size() << " cells: calculateSpatialTranslation tc " << tc << " by dt " << P::timeclassDt[tc] <<"\n";
            if (P::vlasovSolverGhostTranslate) {
               // Local translation without interim communication
               uint tictoc = Timeclasses::TIC;
               if(P::timeclassTictoc){
                  if (tc == 0) { // we are the coarsest timeclass, so we always have fine timeclass data to fetch and sync
                     tictoc = Timeclasses::TOC;
                  }
                  else { // Other timeclasses need to get time ghost halo layers from the coarser timeclasses
                        // to have the intermediate translation step (TIC)
                     if ((P::fractionalTimestep % mod2) == 0){
                        // Initial timeghost halo layer needs to have a larger source region
                        tictoc = Timeclasses::TIC;
                     }
                     else{
                        // Subsequent timeghost halo layers can do with a smaller layer, using the translation
                        // targets from the TIC
                        tictoc = Timeclasses::TOC;
                     }
                  }
               }
               string tictocstr = (tictoc == Timeclasses::TIC ? " tic" : " toc");
               string profNamet = profName + tictocstr;
               phiprof::Timer timer {profNamet}; 
               // std::cerr << "rank " << myRank << ", time " << P::t << " " << tc_propagated_cells[tc].size() << " cells: calculateSpatialTranslation tc " << tc << tictocstr << " by dt " << P::timeclassDt[tc] <<"\n";
               calculateSpatialGhostTranslation(
                  mpiGrid,
                  tc_propagated_cells[tc], // Used for LB
                  nPencils.at(tc),
                  P::timeclassDt[tc],
                  popID,
                  time,
                  tictoc,
                  tc
                  );
            } else {
               calculateSpatialTranslation(
                  mpiGrid,
                  tc_propagated_cells[tc], //local_propagated_cells,
                  remoteTargetCellsx,
                  remoteTargetCellsy,
                  remoteTargetCellsz,
                  nPencils.at(tc),
                  P::timeclassDt[tc],//dt,
                  popID,
                  time
               );
            }
         }
      }
   }

   if (Parameters::prepareForRebalance == true) {
      // clear weight on all local cells
      for (size_t c=0; c<localCells.size(); ++c) {
         SpatialCell* SC = mpiGrid[localCells[c]];
         SC->parameters[CellParams::LBWEIGHTCOUNTER] = 0;
      }
      for (size_t c=0; c<local_propagated_cells.size(); ++c) {
         // Gather total blocks in cell
         SpatialCell* SC = mpiGrid[local_propagated_cells[c]];
         Real counter = 0;
         for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
            for(int timeclass = 0; timeclass <= P::currentMaxTimeclass; ++timeclass){
              counter += pow(P::timeclassLBmantissa,timeclass)*SC->get_number_of_velocity_blocks(popID, timeclass);
            }
         }

         // int accelerationsteps = 0; // Account for time spent in acceleration as well
         // if (mpiGrid[local_propagated_cells[c]]->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
         //    accelerationsteps = 3;
         // }
         if (SC->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY) {
            // Set sysb cells to a small weight
            SC->parameters[CellParams::LBWEIGHTCOUNTER] = counter * 0.5;
         } else {
            if (P::amrMaxSpatialRefLevel == 0) {
               SC->parameters[CellParams::LBWEIGHTCOUNTER] = 3 * counter;
            } else {
               // SC->parameters[CellParams::LBWEIGHTCOUNTER] = nPencils[c] * counter;
               SC->parameters[CellParams::LBWEIGHTCOUNTER] = 3 * counter;
               // SC->parameters[CellParams::LBWEIGHTCOUNTER] += (nPencils[c]+accelerationsteps) * counter;
               // SC->parameters[CellParams::LBWEIGHTCOUNTER] += time / localCells.size();
            }
         }
      }
   }

   // This loop saves the _R-moments before updating into a previous buffer so they can be used for interpolating
   // for timeclasses.
   for (int tc=0; tc <= P::currentMaxTimeclass; tc++) {
      int mod = 1 << (P::currentMaxTimeclass - tc);
      if ((P::fractionalTimestep % mod) == 0) {
         for (size_t c=0; c<tc_propagated_cells.at(tc).size(); ++c) {
            const CellID cellID = tc_propagated_cells.at(tc).at(c);
            SpatialCell* SC = mpiGrid[cellID];

            SC->parameters[CellParams::RHOM_R_PREV_PREV] = SC->parameters[CellParams::RHOM_R_PREV];
            SC->parameters[CellParams::VX_R_PREV_PREV] = SC->parameters[CellParams::VX_R_PREV];
            SC->parameters[CellParams::VY_R_PREV_PREV] = SC->parameters[CellParams::VY_R_PREV];
            SC->parameters[CellParams::VZ_R_PREV_PREV] = SC->parameters[CellParams::VZ_R_PREV];
            SC->parameters[CellParams::RHOQ_R_PREV_PREV] = SC->parameters[CellParams::RHOQ_R_PREV];
            SC->parameters[CellParams::P_11_R_PREV_PREV] = SC->parameters[CellParams::P_11_R_PREV];
            SC->parameters[CellParams::P_22_R_PREV_PREV] = SC->parameters[CellParams::P_22_R_PREV];
            SC->parameters[CellParams::P_33_R_PREV_PREV] = SC->parameters[CellParams::P_33_R_PREV];
            SC->parameters[CellParams::P_23_R_PREV_PREV] = SC->parameters[CellParams::P_23_R_PREV];
            SC->parameters[CellParams::P_13_R_PREV_PREV] = SC->parameters[CellParams::P_13_R_PREV];
            SC->parameters[CellParams::P_12_R_PREV_PREV] = SC->parameters[CellParams::P_12_R_PREV];

            SC->parameters[CellParams::RHOM_R_PREV] = SC->parameters[CellParams::RHOM_R];
            SC->parameters[CellParams::VX_R_PREV] = SC->parameters[CellParams::VX_R];
            SC->parameters[CellParams::VY_R_PREV] = SC->parameters[CellParams::VY_R];
            SC->parameters[CellParams::VZ_R_PREV] = SC->parameters[CellParams::VZ_R];
            SC->parameters[CellParams::RHOQ_R_PREV] = SC->parameters[CellParams::RHOQ_R];
            SC->parameters[CellParams::P_11_R_PREV] = SC->parameters[CellParams::P_11_R];
            SC->parameters[CellParams::P_22_R_PREV] = SC->parameters[CellParams::P_22_R];
            SC->parameters[CellParams::P_33_R_PREV] = SC->parameters[CellParams::P_33_R];
            SC->parameters[CellParams::P_23_R_PREV] = SC->parameters[CellParams::P_23_R];
            SC->parameters[CellParams::P_13_R_PREV] = SC->parameters[CellParams::P_13_R];
            SC->parameters[CellParams::P_12_R_PREV] = SC->parameters[CellParams::P_12_R];
         }
      }
   }
   // Mapping complete, update moments and maximum dt limits //

   //calculating moments based on timeclass
   for (int tc=0; tc <= P::currentMaxTimeclass; tc++) {
      int mod = 1 << (P::currentMaxTimeclass - tc);
      if ((P::fractionalTimestep % mod) == 0) {
         calculateMoments_R(mpiGrid,tc_propagated_cells.at(tc),true);
      }
   }
   if (P::vlasovSolverGhostTranslate){
      phiprof::Timer postBarrierTimer {"MPI barrier-post-trans"};
      MPI_Barrier(MPI_COMM_WORLD);
      postBarrierTimer.stop();
   }
}

/*
  --------------------------------------------------
  Acceleration (velocity space propagation)
  --------------------------------------------------
*/

/**
 * Compute the number of subcycles needed for the acceleration of the particle
 * species in the spatial cell. Note that one should first prepare to
 * accelerate the cell with updateAccelerationMaxdt.
 *
 * @param spatial_cell Spatial cell containing the accelerated population.
 * @param popID ID of the accelerated particle species.
*/

int getAccelerationSubcycles(const SpatialCell* spatial_cell, const Real dt, const uint popID)
{
   //return max( convert<uint>(ceil(dt*spatial_cell->CellParams[CELLPARAMS::TIMECLASSDT] / spatial_cell->get_max_v_dt(popID))), 1u);
   return max( convert<int>(ceil(dt / spatial_cell->get_max_v_dt(popID))), 1);
}
/*!
  Compute the number of subcycles needed from maxVdt and target dt.

 * @param spatial_cell Spatial cell containing the accelerated population.
 * @param popID ID of the accelerated particle species.
*/

int getAccelerationSubcycles(Real maxVdt, Real dt)
{
   //return max( convert<uint>(ceil(dt*spatial_cell->CellParams[CELLPARAMS::TIMECLASSDT] / spatial_cell->get_max_v_dt(popID))), 1u);
   return max( convert<int>(ceil(dt / maxVdt)), 1);
}

/** Accelerate the given population to new time t+dt.
 * @param popID Particle population ID.
 * @param globalMaxSubcycles Number of times acceleration is subcycled.
 * @param step The current subcycle step.
 * @param mpiGrid Parallel grid library.
 * @param acceleratedCells List of cells in which the population is accelerated.
 * @param dt Timestep factor.
 * @param tc Timeclass */
 
void calculateAcceleration(const uint popID,const int globalMaxSubcycles,const int step,
                           dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                           // const std::vector<CellID>& propagatedCells,
                           const std::vector<AccelerationPayload> acceleratedCells,
                           const Real& dt, const int tc = -1) {
   // Set active population
   SpatialCell::setCommunicatedSpecies(popID,tc);

   // Calculate velocity moments, these are needed to
   // calculate the transforms used in the accelerations.
   // Calculated moments are stored in the "_V" variables.
   std::set<CellID> propagatedSet = {}; // Find the unique cells from the payload vector
   std::map<int,set<CellID>> tcPropagatedSets;
   for (int timeclass = 0; timeclass <= P::currentMaxTimeclass; ++timeclass){
      tcPropagatedSets[timeclass] = {};
   }
   for(auto& payload:acceleratedCells){
      propagatedSet.insert(payload.cellptr->get_cellid());
      tcPropagatedSets[payload.timeclass].insert(payload.cellptr->get_cellid());
   }

   calculateMoments_V(mpiGrid, std::vector<CellID>(propagatedSet.begin(),propagatedSet.end()), false);

   // set seed, initialise generator and get value. The order is the same
   // for all cells, but varies with timestep.
   std::default_random_engine rndState;
   //rndState.seed(P::tstep);
   rndState.seed(P::tstep * (int)(pow(2.0, P::currentMaxTimeclass)) + P::fractionalTimestep);
   //rndState.seed(P::tstep + P::fractionalTimestep); // WARNING this formulation actually has some correlations (P::tstep + P::fractionalTimestep can do aliasing...)
   #ifndef DEBUG_TIMECLASSES
      uint map_order=std::uniform_int_distribution<>(0,2)(rndState);
   #else
      uint map_order=1;
   #endif
   // Semi-Lagrangian acceleration for those cells which are subcycled
   // std::cout << "Propagating (ACC) cells ";
   // #pragma omp parallel for schedule(dynamic,1)
   // for (size_t c=0; c<propagatedCells.size(); ++c) {
   //    const CellID cellID = propagatedCells[c];
   //    std::cout << cellID << " ";
   // }
   // std::cout << "\n";
   // for (size_t c=0; c<propagatedCells.size(); ++c) {
      // std::cout << "Accelerating payload at fract " << P::fractionalTimestep << ", cells(timeclass,dt,): \n\t"; }
   // for(auto& payload:acceleratedCells){
      
   //    std::cout << payload.cellptr->get_cellid() <<"("<< payload.timeclass << ", " << payload.dt <<") ";
   // }
   // std::cout<<"\n";
   for(auto& payload:acceleratedCells){
      
      
      const CellID cellID = payload.cellptr->get_cellid();
      const Real maxVdt = payload.cellptr->get_max_v_dt(popID);//mpiGrid[cellID]->get_max_v_dt(popID);
      Real celldt = payload.dt;//dt*mpiGrid[cellID]->get_tc_dt();

      
      //compute subcycle dt. The length is maxVdt on all steps
      //except the last one. This is to keep the neighboring
      //spatial cells in sync, so that two neighboring cells with
      //different number of subcycles have similar timesteps,
      //except that one takes an additional short step. This keeps
      //spatial block neighbors as much in sync as possible for
      //adjust blocks.
      Real thisSubcycleDt;
      if( (step + 1) * maxVdt > fabs(celldt)) {
	      thisSubcycleDt = max(fabs(celldt) - step * maxVdt, 0.0);
      } else{
         thisSubcycleDt = maxVdt;
      }
      if (dt<0) {
         thisSubcycleDt = -thisSubcycleDt;
      }
      spatial_cell::Population& pop = mpiGrid[cellID]->get_population(popID, payload.timeclass);
      pop.subcycleDt = thisSubcycleDt;
   }

      //generate pseudo-random order which is always the same irrespective of parallelization, restarts, etc.
      //std::default_random_engine rndState;
      // set seed, initialise generator and get value. The order is the same
      // for all cells, but varies with timestep.
      //rndState.seed(P::tstep + P::fractionalTimestep); // WARNING this formulation actually has some correlations (P::tstep + P::fractionalTimestep can do aliasing...)


   phiprof::Timer semilagAccTimer {"cell-semilag-acc"};
#ifdef USE_GPU
   gpu_accelerate_cells(mpiGrid[cellID],popID,map_order,subcycleDt);
#else
   cpu_accelerate_cells(mpiGrid, acceleratedCells, popID, map_order);
#endif

   for (auto& payload:acceleratedCells){
      spatial_cell::Population& pop = payload.cellptr->get_population(popID, payload.timeclass);
      pop.T_V += pop.subcycleDt;
   }
   semilagAccTimer.stop();

   //global adjust after each subcycle to keep number of blocks managable. Even the ones not
   //accelerating anyore participate. It is important to keep
   //the spatial dimension to make sure that we do not loose
   //stuff streaming in from other cells, perhaps not connected
   //to the existing distribution function in the cell.
   //- All cells update and communicate their lists of content blocks
   //- Only cells which were accerelated on this step need to be adjusted (blocks removed or added).
   //- Not done here on last step (done after loop)

   if(step < (globalMaxSubcycles - 1)){
      for (int timeclass = 0; timeclass <= P::currentMaxTimeclass; timeclass++){ // Filter to necessary tcs - adjustVelocityBlocks comms need COMM-WORLD yet     
         adjustVelocityBlocks(mpiGrid, std::vector<CellID>(tcPropagatedSets[timeclass].begin(),tcPropagatedSets[timeclass].end()), false, popID, timeclass);
      }
   }
}

/* Set timeghosts from the native VDF state before translation on the corresponding timeclass
*/
void setAccelerationTimeGhosts( SpatialCell* spatial_cell, const uint popID){

   bool addPayload = false;
   for(auto i : spatial_cell->get_all_ghosts()) {
         // Example: On tc-0 cell, tc-1 ghosts requested ghosts of tc-0
         /*               |0--1/4-2/4-3/4-4/4-5/4-6/4--|
         tc-1 after-acc   |----x-------x-------x-------|
         tc-1 after-trans |x-------x-------x-------x---|
         tc-0 after-acc   |----o---x---o-------o---x---|
         tc-0 after-trans |x---------------x-----------|
         legend:
         - x: true state 
         - o: ghost state
         */
      int tc_d = i-spatial_cell->get_tc();
      if(tc_d > 0) {
         // Ghosts for finer cadences
         if (true /*spatial_cell->requested_timeclass_ghosts.count(i) > 0*/){ // folding copy-ghosts here, they are copied over at translation if needed
            
            if (!P::tc_leapfrog_init) {
               // Leapfrogging not yet initialized - initialize
               spatial_cell->set_ghost_population(spatial_cell->get_population(popID,spatial_cell->get_tc()),popID,i);
            }
            else if (spatial_cell->get_timeclass_turn_v()) {
               // This cell has been natively propagated. We update the ghost population to the current state, and then acc by half the timeclass dt to sync with the finer timeclass leapfrog state.
               spatial_cell->set_ghost_population(spatial_cell->get_population(popID,spatial_cell->get_tc()),popID,i);
            }
            else if (spatial_cell->get_timeclass_turn_v(i)) {
               // The finer timeclass is to be propagated, and we do not have a native state to update the VDF from. Just propagate.
            }
            else{
               //do nothing
            }
         }

         else{
             // No ghosts requested from at this timeclass, do nothing
         }
         
      }
      else if (tc_d < 0){
         // Ghosts for coarser cadences.
         // Example: on tc-1 cell, tc-0 ghosts requested from us
         /*                    |0--1/4-2/4-3/4-4/4-5/4-6/4--|
         tc-1      after-acc   |----x-------x---.---x-------|
         tc-1      after-trans |x-------x-------x-------x---|
         tc-1ghost after-acc   |--------o---------------o---|
         tc-1ghost after-trans |o---------------o-----------|
         tc-0      after-acc   |--------x---------------x---|
         tc-0      after-trans |x---------------x-----------|
         legend:
         - x: true state 
         - o: ghost state
         */
         if (true /*spatial_cell->requested_timeclass_ghosts.count(i) > 0*/){

            if (!P::tc_leapfrog_init){
                              
               // if it is slower-tc's turn, we are synced at after-trans state
               // -> Copy state, but needs to acc by half-tc-0-dt ("always init")
               spatial_cell->set_ghost_population(spatial_cell->get_population(popID,spatial_cell->get_tc()),popID,i);

            }
            else if (spatial_cell->get_timeclass_turn_v(i) || (P::tstep == 0 && P::fractionalTimestep == 0)) {
               
               // if it is slower-tc's turn, we are synced at after-trans state
               // -> Copy state, but needs to acc by half-tc-0-dt ("always init")
               spatial_cell->set_ghost_population(spatial_cell->get_population(popID,spatial_cell->get_tc()),popID,i);
            }
            else{ //do nothing
            }
         }

         else{
             // No ghosts requested from at this timeclass, do nothing
         }
      }
   } // for over ghost requests
}


/* Get structs (timeclass, SpatialCell*, dt, subcycle step) for acceleration
*/
void setAccelerationPayloads(vector<AccelerationPayload>& outvec, SpatialCell* spatial_cell, const uint popID, const Real& dt){

   // if(spatial_cell->get_timeclass_turn_v() || (accelerateSpecificCells && std::find(cellsToAccelerate.begin(), cellsToAccelerate.end(), spatial_cell->get_cellid()) != cellsToAccelerate.end())){
   //    AccelerationPayload payload = {spatial_cell->get_tc(),spatial_cell,dt,0};
   //    outvec.push_back(payload);
   // }

   if(spatial_cell->get_timeclass_turn_v()){
      outvec.push_back(AccelerationPayload(spatial_cell->get_tc(),spatial_cell,dt));
   }

   // std::cout << "outvec size: " << outvec.size() << "\n";

   bool addPayload = false;
   for(auto i : spatial_cell->get_all_ghosts()) {
         // Example: On tc-0 cell, tc-1 ghosts requested ghosts of tc-0
         /*               |0--1/4-2/4-3/4-4/4-5/4-6/4--|
         tc-1 after-acc   |----x-------x-------x-------|
         tc-1 after-trans |x-------x-------x-------x---|
         tc-0 after-acc   |----o---x---o-------o---x---|
         tc-0 after-trans |x---------------x-----------|
         legend:
         - x: true state 
         - o: ghost state
         */
      int tc_d = i-spatial_cell->get_tc();
      if(tc_d > 0) {
         // Ghosts for finer cadences
         if (true /*spatial_cell->requested_timeclass_ghosts.count(i) > 0*/){ // folding copy-ghosts here, they are copied over at translation if needed
            
            if (!P::tc_leapfrog_init) {
               // Leapfrogging not yet initialized - initialize
               // std::cerr << "leapfrog init for " << spatial_cell->parameters[CellParams::CELLID] << " for tc " << i << std::endl;
               // spatial_cell->set_ghost_population(spatial_cell->get_population(popID,spatial_cell->get_tc()),popID,i);
               
               outvec.push_back(AccelerationPayload(i, spatial_cell, dt/pow(2,tc_d)));
            }
            else if (spatial_cell->get_timeclass_turn_v()) {
               // This cell has been natively propagated. We update the ghost population to the current state, and then acc by half the timeclass dt to sync with the finer timeclass leapfrog state.
               // spatial_cell->set_ghost_population(spatial_cell->get_population(popID,spatial_cell->get_tc()),popID,i);

               outvec.push_back(AccelerationPayload(i, spatial_cell, dt/pow(2,tc_d)*1/2));
            }
            else if (spatial_cell->get_timeclass_turn_v(i)) {
               // The finer timeclass is to be propagated, and we do not have a native state to update the VDF from. Just propagate.
               outvec.push_back(AccelerationPayload(i, spatial_cell, dt/pow(2,tc_d)));
            }
            else{
               //do nothing
            }
         }
         // else if(spatial_cell->requested_timeclass_copy_ghosts.count(i) > 0 || spatial_cell->requested_timeclass_outer_ghosts.count(i) > 0){
         //    // Copy ghosts needed, and regular ghosts are not available, so initialize and accelerate those as needed.
         //    if (!P::tc_leapfrog_init) {
         //       // Leapfrogging not yet initialized - initialize
         //       spatial_cell->set_ghost_population(spatial_cell->get_population(popID,spatial_cell->get_tc()),popID,i);

         //       outvec.push_back(AccelerationPayload(i, spatial_cell, dt/pow(2,tc_d)));
         //    }
         //    else if (spatial_cell->get_timeclass_turn_v(i)){
         //       spatial_cell->set_ghost_population(spatial_cell->get_population(popID,spatial_cell->get_tc()),popID,i);

         //       outvec.push_back(AccelerationPayload(i, spatial_cell, dt/pow(2,tc_d)*3/2));
         //    }
         //    else{
         //       spatial_cell->set_ghost_population(spatial_cell->get_population(popID,spatial_cell->get_tc()),popID,i);
         //    }
         // }
         else{
             // No ghosts requested from at this timeclass, do nothing
         }
         
      }
      else if (tc_d < 0){
         // Ghosts for coarser cadences.
         // Example: on tc-1 cell, tc-0 ghosts requested from us
         /*                    |0--1/4-2/4-3/4-4/4-5/4-6/4--|
         tc-1      after-acc   |----x-------x---.---x-------|
         tc-1      after-trans |x-------x-------x-------x---|
         tc-1ghost after-acc   |--------o---------------o---|
         tc-1ghost after-trans |o---------------o-----------|
         tc-0      after-acc   |--------x---------------x---|
         tc-0      after-trans |x---------------x-----------|
         legend:
         - x: true state 
         - o: ghost state
         */
         if (true /*spatial_cell->requested_timeclass_ghosts.count(i) > 0*/){

            if (!P::tc_leapfrog_init){
                              
               // if it is slower-tc's turn, we are synced at after-trans state
               // -> Copy state, but needs to acc by half-tc-0-dt ("always init")
               // spatial_cell->set_ghost_population(spatial_cell->get_population(popID,spatial_cell->get_tc()),popID,i);

               outvec.push_back(AccelerationPayload(i, spatial_cell,  dt/pow(2,tc_d)));
            }
            else if (spatial_cell->get_timeclass_turn_v(i) || (P::tstep == 0 && P::fractionalTimestep == 0)) {
               
               // if it is slower-tc's turn, we are synced at after-trans state
               // -> Copy state, but needs to acc by half-tc-0-dt ("always init")
               // spatial_cell->set_ghost_population(spatial_cell->get_population(popID,spatial_cell->get_tc()),popID,i);

               outvec.push_back(AccelerationPayload(i, spatial_cell,  5./4.*dt/pow(2,tc_d)));
            }
            else{ //do nothing
            }
         }
         // else if(spatial_cell->requested_timeclass_copy_ghosts.count(i) > 0 || spatial_cell->requested_timeclass_outer_ghosts.count(i) > 0){
         //    // Copy ghosts needed, and regular ghosts are not available, so initialize and accelerate those as needed.
         //    if (!P::tc_leapfrog_init){
                              
         //       // if it is slower-tc's turn, we are synced at after-trans state
         //       // -> Copy state, but needs to acc by half-tc-0-dt ("always init")
         //       spatial_cell->set_ghost_population(spatial_cell->get_population(popID,spatial_cell->get_tc()),popID,i);

         //       outvec.push_back(AccelerationPayload(i, spatial_cell,  dt/pow(2,tc_d)));
         //    }
         //    else if (spatial_cell->get_timeclass_turn_v(i)){
         //       spatial_cell->set_ghost_population(spatial_cell->get_population(popID,spatial_cell->get_tc()),popID,i);

         //       outvec.push_back(AccelerationPayload(i, spatial_cell, dt/pow(2,tc_d)*5/4));
         //    }
         //    else{
         //       spatial_cell->set_ghost_population(spatial_cell->get_population(popID,spatial_cell->get_tc()),popID,i);
         //    }
         // }
         else{
             // No ghosts requested from at this timeclass, do nothing
         }
      }
   } // for over ghost requests
}

/** Accelerate all particle populations to new time t+dt.
 * This function is AMR safe. NB: acceleration by dt = 0 is not a no-op.
 * @param mpiGrid Parallel grid library.
 * @param dt Time step factor: cells will propagated by dt*CellParams[CellParams::CELLDT] if needed.*/
void calculateAcceleration(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                           const Real dt,
                           bool accelerateSpecificCells,
                           const vector<CellID>& cellsToAccelerate
                          ) {
   typedef Parameters P;
   const vector<CellID>& cells = getLocalCells();
   set<CellID> cellsToPropagateSet;
   vector<CellID> cellsToPropagateVector;
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD,&myRank);

   if (dt==0.0) {
      for (auto c : cells) {
         cellsToPropagateSet.insert(c);
      }
   }

   // std::cout << "-----------------calculateAcceleration at t="<<P::t << ", for dtfactor="<<dt<<"\n";
   if (dt == 0.0 && (P::tstep > 0 || P::fractionalTimestep > 0)) {

      // Even if acceleration is turned off we need to adjust velocity blocks
      // because the boundary conditions may have altered the velocity space,
      // and to update changes in no-content blocks during translation.
      for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
// std::cerr << __FILE__<<":"<<__LINE__<< " calling adjustVelocityBlocks at t = "
//          << P::t << ", preparing to receive; len cells = " << cells.size() <<
//          "\n";
         for (int tc = 0; tc <= P::currentMaxTimeclass; tc++){ // Filter to necessary tcs
            adjustVelocityBlocks(mpiGrid, cells, true, popID, tc);
         }
      }
   //    // Needed to calculate _V moments afterwards for the local cells
   //    // TODO also for ghost cells?
   //    for (auto c : cells){
   //       cellsToPropagateSet.insert(c);
   //    }
   // } else if (dt == 0.0 && (P::tstep == 0 && P::fractionalTimestep == 0)) { // zero dt step to initialize moments at the start of the sim, collect all cells
   //    for (auto c : cells){
   //       cellsToPropagateSet.insert(c);
   //    }
   } else {
      // Fairly ugly but no goto
      phiprof::Timer accTimer {"semilag-acc"};

      std::set<int> timeclasses_handled;
      // Accelerate all particle species
      for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
         int maxSubcycles=0;
         int globalMaxSubcycles;

         // Set active population
         SpatialCell::setCommunicatedSpecies(popID); // FIXME: should this have a timeclass argument?

         vector<AccelerationPayload> propagatePayloads; // <timeclass, SpatialCell*, dt, subcycle step>

         // Iterate through all local cells and collect cells to propagate.
         // Ghost cells (spatial cells at the boundary of the simulation
         // volume) do not need to be propagated:
         // vector<CellID> propagatedCells;


         if (!accelerateSpecificCells) {
            //std::cout << "Accelerating all cells for popID " << popID << "\n";
            for (size_t c=0; c<cells.size(); ++c) {
               SpatialCell* SC = mpiGrid[cells[c]];
               Real dt_cell;
               if(dt < 0.0) { // Revert to previous real time, stored in cell
                     //dt_cell = SC->parameters[CellParams::TIME_R] - P::t;
                     dt_cell = dt*SC->get_tc_dt();

                  }
               else { // dt is a factor of 0.5 or 1.0; so dt_cell is local timestep * dt factor
                     dt_cell = dt*SC->get_tc_dt();
                  }
               // if ( (SC->parameters[CellParams::CELLID] == 9 || SC->parameters[CellParams::CELLID] == 11 || SC->parameters[CellParams::CELLID] == 12))
                  // std::cout << "vdt on tc  " << SC->parameters[CellParams::TIMECLASS] << " on ftstep " << P::fractionalTimestep << ", dt " << dt_cell <<"\n";

               const vmesh::VelocityMesh* vmesh = SC->get_velocity_mesh(popID);
               // disregard boundary cells, in preparation for acceleration
               if (  (SC->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) ||
                     // Include inflow-Maxwellian
                     (P::vlasovAccelerateMaxwellianBoundaries && (SC->sysBoundaryFlag == sysboundarytype::MAXWELLIAN)) ) {
                        if (vmesh->size() != 0){   //do not propagate spatial cells with no blocks

                              setAccelerationPayloads(propagatePayloads, mpiGrid[cells[c]], 0, dt_cell);
                              // propagatePayloads.insert(propagatePayloads.end(), outvec.begin(),outvec.end());
                  
                        }
                        //prepare for acceleration, updates max dt for each cell, it
                        //needs to be set to something sensible for _all_ cells, even if
                        //they are not propagated
                        updateAccelerationMaxdt(SC, popID);

                     //update max subcycles for all cells in this process - population-specific
                     // maxSubcycles = max((int)getAccelerationSubcycles(SC->get_max_v_dt(popID), dt_cell), maxSubcycles);
                     // spatial_cell::Population& pop = SC->get_population(popID);
                     // pop.ACCSUBCYCLES = getAccelerationSubcycles(SC->get_max_v_dt(popID), dt_cell);
                     #ifdef USE_GPU
                     #pragma omp critical
                     {
                        if (blockCount > gpuMaxBlockCount) {
                           gpuMaxBlockCount = blockCount;
                        }
                     }
                     #endif                  
               } // if propagate loop
            } // for loop over cells
         } // if !accelerateSpecificCells
         else {
            //std::cout << "Accelerating specific cells for popID " << popID << "\n";
            for (size_t c=0; c<cellsToAccelerate.size(); ++c) {
               SpatialCell* SC = mpiGrid[cellsToAccelerate[c]];
               Real dt_cell;
               if(dt < 0.0) { // Revert to previous real time, stored in cell
                     // dt_cell = SC->parameters[CellParams::TIME_R] - P::t;
                     dt_cell = dt*SC->get_tc_dt();

                  }
               else { // dt is a factor of 0.5 or 1.0; so dt_cell is local timestep * dt factor
                     dt_cell = dt*SC->get_tc_dt();
                  }
                  std::cout << "dt_cell: " << dt_cell << "\n";

               // if ( (SC->parameters[CellParams::CELLID] == 9 || SC->parameters[CellParams::CELLID] == 11 || SC->parameters[CellParams::CELLID] == 12))
                  // std::cout << "vdt on tc  " << SC->parameters[CellParams::TIMECLASS] << " on ftstep " << P::fractionalTimestep << ", dt " << dt_cell <<"\n";

               const vmesh::VelocityMesh* vmesh = SC->get_velocity_mesh(popID);
               // disregard boundary cells, in preparation for acceleration
               if (  (SC->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) ||
                     // Include inflow-Maxwellian
                     (P::vlasovAccelerateMaxwellianBoundaries && (SC->sysBoundaryFlag == sysboundarytype::MAXWELLIAN)) ) {
                        if (vmesh->size() != 0){   //do not propagate spatial cells with no blocks
                              setAccelerationPayloads(propagatePayloads, mpiGrid[cellsToAccelerate[c]], 0, dt_cell);
                        }
                        //prepare for acceleration, updates max dt for each cell, it
                        //needs to be set to something sensible for _all_ cells, even if
                        //they are not propagated
                        updateAccelerationMaxdt(SC, popID);

                     #ifdef USE_GPU
                     #pragma omp critical
                     {
                        if (blockCount > gpuMaxBlockCount) {
                           gpuMaxBlockCount = blockCount;
                        }
                     }
                     #endif                  
               } // if propagate loop
            } // for loop over cells
         } // else accelerateSpecificCells
         // This can possibly have funky dts, lets check those 
         for(auto& payload : propagatePayloads){
            //
            //update max subcycles for all cells in this process - population-specific
            // std::cerr << __FILE__ <<":"<<__LINE__<<" "<<myRank <<" : cell " << payload.cellptr->parameters[CellParams::CELLID] << " acc dt " << payload.cellptr->get_max_v_dt(popID) << "; payload.dt " << payload.dt << " maxSubcycles " << maxSubcycles << "\n";
            maxSubcycles = max(getAccelerationSubcycles(payload.cellptr->get_max_v_dt(popID), payload.dt), maxSubcycles);
            spatial_cell::Population& pop = payload.cellptr->get_population(popID,payload.timeclass);
            pop.ACCSUBCYCLES = getAccelerationSubcycles(payload.cellptr->get_max_v_dt(popID), payload.dt);
            timeclasses_handled.insert(payload.cellptr->get_tc());
         }

#ifdef USE_GPU
         // Ensure accelerator has enough temporary memory allocated
         phiprof::Timer verificationTimer {"gpu allocation verifications"};
         gpu_vlasov_allocate(gpuMaxBlockCount);
         gpu_acc_allocate(gpuMaxBlockCount);
         verificationTimer.stop();
#endif

         // Compute global maximum for number of subcycles
         // std::cerr << __FILE__ <<":"<<__LINE__<<" " <<myRank <<": maxSubcycles = " << maxSubcycles <<"\n";
         MPI_Allreduce(&maxSubcycles, &globalMaxSubcycles, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
         // std::cerr << __FILE__ <<":"<<__LINE__<<" " <<myRank <<": globalMaxSubcycles = " << globalMaxSubcycles <<"\n";
         // propagatedCells = std::vector<CellID>(cellsToPropagateSet.begin(),cellsToPropagateSet.end());
         std::vector<CellID> propagatedCells = std::vector<CellID>(); // so that it compiles, REMOVE

         // propagatedCells = std::vector<AccelerationPayload>(propagatePayloads.begin(),propagatePayloads.end());
         // substep global max times
         for(int step=0; step<globalMaxSubcycles; ++step) {
            if(step > 0) {
               // prune list of cells to propagate to only contained those which are now subcycled
               std::vector<AccelerationPayload> temp;
               for (const auto& payload: propagatePayloads) {
                  // if (step < getAccelerationSubcycles(mpiGrid[cell], dt*mpiGrid[cell]->get_tc_dt(), popID) ) {
                  if (step < getAccelerationSubcycles(payload.cellptr->get_max_v_dt(popID), payload.dt) ) {
                     temp.push_back(payload);
                  }
               }

               propagatePayloads.swap(temp);
            }
            // Accelerate population over one subcycle step
            // std::cout << "--------------------- an ACC subcycle " << step << " ------------------\n";
            calculateAcceleration(popID,globalMaxSubcycles,step,mpiGrid,propagatePayloads,dt);
            if (step == 0){ // we're not getting new cells after the first subcycle step
               for(auto& payload:propagatePayloads){
                  cellsToPropagateSet.insert(payload.cellptr->get_cellid());
               }
            }
         } // for-loop over acceleration substeps

         // final adjust for all cells, also fixing remote cells.
// std::cerr << __FILE__<<":"<<__LINE__<< " calling adjustVelocityBlocks at t = "
//          << P::t << ", preparing to receive; len cells = " << cells.size() <<
//          "\n";

         for (int timeclass : timeclasses_handled){
            for(auto& payload:propagatePayloads){
               spatial_cell::Population& pop = payload.cellptr->get_population(popID, payload.timeclass);
               if (pop.state & (1u << LastPropagationState::TRANSLATED)||(pop.state == LastPropagationState::UNINITIALIZED)){
                  pop.state = LastPropagationState::ACCELERATED;
               }
               else{
                  if(do_translate_cell(payload.cellptr) && dt != 0.0 && pop.state != LastPropagationState::ACCELERATED){
                     std::cerr << __FILE__ << ":" << __LINE__ << ": Accelerated a badly-translated cell " << payload.cellptr->parameters[CellParams::CELLID] << "t" << payload.cellptr->parameters[CellParams::TIMECLASS] << " t = " << P::t << " at timeclass " << timeclass << ", state = " << std::bitset<8>(pop.state) << " tG flags " << payload.cellptr->requested_timeclass_ghosts.count(timeclass) << ", " << payload.cellptr->requested_timeclass_copy_ghosts.count(timeclass) << ", " << payload.cellptr->requested_timeclass_outer_ghosts.count(timeclass) << std::endl;
                  }
                  pop.state = LastPropagationState::ACCELERATED;
               }
            }
         }
         for (int tc = 0; tc <= P::currentMaxTimeclass; tc++){ // Filter to necessary tcs - adjustVelocityBlocks comms need COMM-WORLD yet
            adjustVelocityBlocks(mpiGrid, cells, true, popID, tc);
         }
      } // for-loop over particle species

      //now cellsToPropagateSet contains all cells which have been propagated and whose moments need updating

   } //else

   // std::cout << "---------------------------- ACC finished --------------------\n";


   //converting cellsToPropagateSet to vector
   for (auto cell : cellsToPropagateSet) {
      cellsToPropagateVector.push_back(cell);
   }

   // Recalculate "_V" velocity moments
   calculateMoments_V(mpiGrid,cellsToPropagateVector,true,(dt==0.0));

   // std::cerr << "calculated V moments...\n";

   // Set CellParams::MAXVDT to be the minimum dt of all per-species values
   #pragma omp parallel for
   for (size_t c=0; c<cells.size(); ++c) {
      SpatialCell* cell = mpiGrid[cells[c]];
      cell->parameters[CellParams::MAXVDT] = numeric_limits<Real>::max();
      for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
         cell->parameters[CellParams::MAXVDT]
            = min(cell->get_max_v_dt(popID), cell->parameters[CellParams::MAXVDT]);
      }
      if(abs(cell->parameters[CellParams::XCRD]+cell->parameters[CellParams::DX]/2) < 6*P::dx_ini) {
            //cell->parameters[CellParams::MAXVDT] /= 2;
            // cout << "maxvdt \n";
         }
      if(abs(cell->parameters[CellParams::XCRD]+cell->parameters[CellParams::DX]/2) < 2*P::dx_ini) {
            //cell->parameters[CellParams::MAXVDT] /= 2;
            // cout << "maxvdt \n";
         }
   }
   // std::cerr << "reaches end of accel top level\n";

}



/*--------------------------------------------------
  Functions for computing moments
  --------------------------------------------------*/

void calculateInterpolatedVelocityMoments(
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   const int cp_rhom,
   const int cp_vx,
   const int cp_vy,
   const int cp_vz,
   const int cp_rhoq,
   const int cp_p11,
   const int cp_p22,
   const int cp_p33,
   const int cp_p23,
   const int cp_p13,
   const int cp_p12
) {
   const vector<CellID>& cells = getLocalCells();

   //Iterate through all local cells
    #pragma omp parallel for
   for (size_t c=0; c<cells.size(); ++c) {
      const CellID cellID = cells[c];
      SpatialCell* SC = mpiGrid[cellID];
      const Real tr = SC->parameters[CellParams::TIME_R];
      const Real tv = SC->parameters[CellParams::TIME_V];
      SC->parameters[cp_rhom  ] = 0.5* ( SC->parameters[CellParams::RHOM_R] + SC->parameters[CellParams::RHOM_V] );
      SC->parameters[cp_vx]   = 0.5* ( SC->parameters[CellParams::VX_R] + SC->parameters[CellParams::VX_V] );
      SC->parameters[cp_vy] = 0.5* ( SC->parameters[CellParams::VY_R] + SC->parameters[CellParams::VY_V] );
      SC->parameters[cp_vz] = 0.5* ( SC->parameters[CellParams::VZ_R] + SC->parameters[CellParams::VZ_V] );
      SC->parameters[cp_rhoq  ] = 0.5* ( SC->parameters[CellParams::RHOQ_R] + SC->parameters[CellParams::RHOQ_V] );
      SC->parameters[cp_p11]   = 0.5* ( SC->parameters[CellParams::P_11_R] + SC->parameters[CellParams::P_11_V] );
      SC->parameters[cp_p22]   = 0.5* ( SC->parameters[CellParams::P_22_R] + SC->parameters[CellParams::P_22_V] );
      SC->parameters[cp_p33]   = 0.5* ( SC->parameters[CellParams::P_33_R] + SC->parameters[CellParams::P_33_V] );
      SC->parameters[cp_p23]   = 0.5* ( SC->parameters[CellParams::P_23_R] + SC->parameters[CellParams::P_23_V] );
      SC->parameters[cp_p13]   = 0.5* ( SC->parameters[CellParams::P_13_R] + SC->parameters[CellParams::P_13_V] );
      SC->parameters[cp_p12]   = 0.5* ( SC->parameters[CellParams::P_12_R] + SC->parameters[CellParams::P_12_V] );

      // for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
      //    spatial_cell::Population& pop = SC->get_population(popID);
      //    pop.RHO = 0.5 * ( pop.RHO_R + pop.RHO_V );
      //    for(int i=0; i<3; i++) {
      //       pop.V[i] = 0.5 * ( pop.V_R[i] + pop.V_V[i] );
      //       pop.P[i]    = 0.5 * ( pop.P_R[i] + pop.P_V[i] );
      //   }
      //}
   }
}

Real linearInterpolation(Real x0, Real y0, Real x1, Real y1, Real x) {
   // https://en.wikipedia.org/wiki/Linear_interpolation
   // this is used in the function below.
   if (!(x >= x0 && x<=x1)) {
      std::cerr << "linearInterpolation: x = " << x << " is out of bounds [" << x0 << ", " << x1 << "]\n";
      assert(false);
   }

   return (y0 * (x1 - x) + y1 * (x - x0))/(x1 - x0);
}

Real lagrangeInterpolation2order(Real x0, Real y0, Real x1, Real y1, Real x2, Real y2, Real x) {
   // https://mathworld.wolfram.com/LagrangeInterpolatingPolynomial.html
   // Lagrange polynomial for interpolation between three points.
   
   if (!(x >= x0 && x<=x2)) {
      std::cerr << "lagrangeInterpolation2order: x = " << x << " is out of bounds [" << x0 << ", " << x2 << "]\n";
      assert(false);
   }

   return (y0 * (x - x1) * (x - x2) / ((x0 - x1) * (x0 - x2)) +
         y1 * (x - x0) * (x - x2) / ((x1 - x0) * (x1 - x2)) +
         y2 * (x - x0) * (x - x1) / ((x2 - x0) * (x2 - x1)));
}

Real lagrangeInterpolation3order(Real x0, Real y0, Real x1, Real y1, Real x2, Real y2, Real x3, Real y3, Real x) {
   // https://mathworld.wolfram.com/LagrangeInterpolatingPolynomial.html
   // Lagrange polynomial for interpolation between four points.

   if (!(x >= x0 && x<=x3)) {
      std::cerr << "lagrangeInterpolation3order: x = " << x << " is out of bounds [" << x0 << ", " << x3 << "]\n";
      assert(false);
   }

   return (y0 * (x - x1) * (x - x2) * (x - x3) / ((x0 - x1) * (x0 - x2) * (x0 - x3)) +
         y1 * (x - x0) * (x - x2) * (x - x3) / ((x1 - x0) * (x1 - x2) * (x1 - x3)) +
         y2 * (x - x0) * (x - x1) * (x - x3) / ((x2 - x0) * (x2 - x1) * (x2 - x3)) +
         y3 * (x - x0) * (x - x1) * (x - x2) / ((x3 - x0) * (x3 - x1) * (x3 - x2)));
}

Real cubicHermiteSplineInterpolation(Real x0, Real y0, Real x1, Real y1, Real x2, Real y2, Real x3, Real y3, Real x) {
   // https://kluge.in-chemnitz.de/opensource/spline/
   // Cubic Hermite spline interpolation between four points.

   if (!(x >= x0 && x<=x3)) {
      std::cerr << "cubicHermiteSplineInterpolation: x = " << x << " is out of bounds [" << x0 << ", " << x3 << "]\n";
      assert(false);
   }

   vector<Real> xvals = {x0, x1, x2, x3};
   vector<Real> yvals = {y0, y1, y2, y3};
   tk::spline s(xvals,yvals,tk::spline::cspline_hermite);

   return s(x);
}


void interpolateMomentsForTimeclasses(
  dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   const int cp_rhom,
   const int cp_rhoq,
   const int cp_p11,
   const int cp_p22,
   const int cp_p33,
   const int cp_p23,
   const int cp_p13,
   const int cp_p12,
   const int cp_vx,
   const int cp_vy,
   const int cp_vz,
   const bool dt2, // true if second moment / dt2
   const bool lastUpdateWasTranslation
) {

   const vector<CellID>& cells = getLocalCells();

   #pragma omp parallel for
   for (size_t c=0; c<cells.size(); ++c) {
      const CellID cellID = cells[c];
      SpatialCell* SC = mpiGrid[cellID];
      const int timeclass = SC->parameters[CellParams::TIMECLASS];

      // this function is called in the main loop after translation, but before acceleration.
      // therefore, if SC->get_timeclass_turn_v() is true, it means that the
      // spatial cell has been translated and the spatial distribution is temporally 
      // ahead of the velocity distribution. If it is false, the velocity distribution
      // is ahead of the spatial distribution (we are on a fractionaltimestep where the cell
      // is not translated / accelerated).

      // ^^ this determines which true moments are used for interpolation.

      if (P::currentMaxTimeclass == 0 || timeclass == P::currentMaxTimeclass) {
         // calculateInterpolatedVelocityMoments functionality here, if timeclass is the max one.
         if (dt2) {

            SC->parameters[cp_rhom  ] = 0.5* ( SC->parameters[CellParams::RHOM_R] + SC->parameters[CellParams::RHOM_V] );
            SC->parameters[cp_rhoq  ] = 0.5* ( SC->parameters[CellParams::RHOQ_R] + SC->parameters[CellParams::RHOQ_V] );
            SC->parameters[cp_p11]   = 0.5* ( SC->parameters[CellParams::P_11_R] + SC->parameters[CellParams::P_11_V] );
            SC->parameters[cp_p22]   = 0.5* ( SC->parameters[CellParams::P_22_R] + SC->parameters[CellParams::P_22_V] );
            SC->parameters[cp_p33]   = 0.5* ( SC->parameters[CellParams::P_33_R] + SC->parameters[CellParams::P_33_V] );
            SC->parameters[cp_p23]   = 0.5* ( SC->parameters[CellParams::P_23_R] + SC->parameters[CellParams::P_23_V] );
            SC->parameters[cp_p13]   = 0.5* ( SC->parameters[CellParams::P_13_R] + SC->parameters[CellParams::P_13_V] );
            SC->parameters[cp_p12]   = 0.5* ( SC->parameters[CellParams::P_12_R] + SC->parameters[CellParams::P_12_V] );
            SC->parameters[cp_vx]   = 0.5* ( SC->parameters[CellParams::VX_R] + SC->parameters[CellParams::VX_V] );
            SC->parameters[cp_vy] = 0.5* ( SC->parameters[CellParams::VY_R] + SC->parameters[CellParams::VY_V] );
            SC->parameters[cp_vz] = 0.5* ( SC->parameters[CellParams::VZ_R] + SC->parameters[CellParams::VZ_V] );

         } else {

            SC->parameters[cp_rhom] = 0.5* ( SC->parameters[CellParams::RHOM_R_PREV] + SC->parameters[CellParams::RHOM_V] );
            SC->parameters[cp_rhoq] = 0.5* ( SC->parameters[CellParams::RHOQ_R_PREV] + SC->parameters[CellParams::RHOQ_V] );
            SC->parameters[cp_p11]   = 0.5* ( SC->parameters[CellParams::P_11_R_PREV] + SC->parameters[CellParams::P_11_V] );
            SC->parameters[cp_p22]   = 0.5* ( SC->parameters[CellParams::P_22_R_PREV] + SC->parameters[CellParams::P_22_V] );
            SC->parameters[cp_p33]   = 0.5* ( SC->parameters[CellParams::P_33_R_PREV] + SC->parameters[CellParams::P_33_V] );
            SC->parameters[cp_p23]   = 0.5* ( SC->parameters[CellParams::P_23_R_PREV] + SC->parameters[CellParams::P_23_V] );
            SC->parameters[cp_p13]   = 0.5* ( SC->parameters[CellParams::P_13_R_PREV] + SC->parameters[CellParams::P_13_V] );
            SC->parameters[cp_p12]   = 0.5* ( SC->parameters[CellParams::P_12_R_PREV] + SC->parameters[CellParams::P_12_V] );
            SC->parameters[cp_vx]   = 0.5* ( SC->parameters[CellParams::VX_R_PREV] + SC->parameters[CellParams::VX_V] );
            SC->parameters[cp_vy] = 0.5* ( SC->parameters[CellParams::VY_R_PREV] + SC->parameters[CellParams::VY_V] );
            SC->parameters[cp_vz] = 0.5* ( SC->parameters[CellParams::VZ_R_PREV] + SC->parameters[CellParams::VZ_V] );

            }

     } else { // this block if timeclass != maxTC

         Real RTCpow = pow(2, P::currentMaxTimeclass - timeclass);
         Real modul = P::fractionalTimestep % (int)RTCpow;
         Real normModul = modul/RTCpow;
         if (dt2) {
            normModul += 0.5/RTCpow; // for dt2, we need to shift the interpolation by 0.5
         }

         // in this function the moments are calculated to an interval [0, 1] that matches a single timestep in real space, or [t, t + dt].
         // for example, the moment "at 0.5" means the moment at t+dt/2.
         // normModul is the point in the interval [0, 1] where the moment is calculated.

         // !! if translation and acceleration are changed to not update both on fractimestep 0, this will break

         if ((P::tcMomentInterpolationType != -1 && P::tcMomentInterpolationType != 1 &&
             P::tcMomentInterpolationType != 2 && P::tcMomentInterpolationType != 3)) {
            std::cerr << "ERROR: Invalid value for P::tcMomentInterpolationType: " << P::tcMomentInterpolationType << "\n";
            exit(1);
         }

         // std::cout << "Interpolating moments for timeclass " << timeclass << " with normModul = " << normModul << " , cellID " << cellID << "\n";
         // std::cout << "Interpolation type: " << P::tcMomentInterpolationType << "\n";
         // std::cout << "Timeclass turn v: " << SC->get_timeclass_turn_v() << "\n";
         // std::cout << "Timeclass: " << timeclass << ", max timeclass: " << P::currentMaxTimeclass << "\n";

         // temporary arrays for true moments.
         const int nMomentsToInterp = 11;

         std::array<Real, nMomentsToInterp> avgMoments1;
         std::array<Real, nMomentsToInterp> avgMoments2;
         std::array<Real, nMomentsToInterp> avgMoments3;
         std::array<Real, nMomentsToInterp> avgMoments4;
         std::array<Real, nMomentsToInterp> avgMoments5;

         // for type of interpolation (P::tcMomentInterpolationType) -1 is cubic C^1 Hermite spline, 1 is linear, 2 is lagrange 2nd order, 3 is lagrange 3rd order.

         if (lastUpdateWasTranslation && SC->get_timeclass_turn_v()) { // aka if translation moments are ahead of acceleration moments
            // here, temporal order from newest to oldest is _R, _V, _R_PREV, _V_PREV, _R_PREV_PREV, _V_PREV_PREV
            for (int i=0; i<nMomentsToInterp; i++) {
               avgMoments1[i] = 0.5*( SC->parameters[CellParams::RHOM_R+i] + SC->parameters[CellParams::RHOM_V+i] ); // at 0.5
               avgMoments2[i] = 0.5*( SC->parameters[CellParams::RHOM_V+i] + SC->parameters[CellParams::RHOM_R_PREV+i] ); // at 0.0
               avgMoments3[i] = 0.5*( SC->parameters[CellParams::RHOM_R_PREV+i] + SC->parameters[CellParams::RHOM_V_PREV+i] );  // at -0.5
               avgMoments4[i] = 0.5*( SC->parameters[CellParams::RHOM_V_PREV+i] + SC->parameters[CellParams::RHOM_R_PREV_PREV+i] ); // at -1.0
               avgMoments5[i] = 0.5*( SC->parameters[CellParams::RHOM_R_PREV_PREV+i] + SC->parameters[CellParams::RHOM_V_PREV_PREV+i] ); // at -1.5
            }

            for (int i=0; i<nMomentsToInterp; i++) {

               switch(P::tcMomentInterpolationType) {
                  case 1:
                     SC->parameters[cp_rhom+i] = linearInterpolation(0.0, avgMoments2[i], 0.5, avgMoments1[i], normModul);
                     break;
                  case 2:
                     SC->parameters[cp_rhom+i] = lagrangeInterpolation2order(-0.5, avgMoments3[i], 0.0, avgMoments2[i], 0.5, avgMoments1[i], normModul);
                     break;
                  case 3:
                     SC->parameters[cp_rhom+i] = lagrangeInterpolation3order(-1.0, avgMoments4[i], -0.5, avgMoments3[i], 0.0, avgMoments2[i], 0.5, avgMoments1[i], normModul);
                     break;
                  case -1:
                     SC->parameters[cp_rhom+i] = cubicHermiteSplineInterpolation(-1.0, avgMoments4[i], -0.5, avgMoments3[i], 0.0, avgMoments2[i], 0.5, avgMoments1[i], normModul);
                     break;
               }
            }
         } else { // if acceleration moments are ahead of translation moments
            // here, temporal order from newest to oldest is _V, _R, _V_PREV, _R_PREV, _V_PREV_PREV, _R_PREV_PREV
            for (int i=0; i<nMomentsToInterp; i++) {
               avgMoments1[i] = 0.5*( SC->parameters[CellParams::RHOM_V+i] + SC->parameters[CellParams::RHOM_R+i] ); // at 1.0
               avgMoments2[i] = 0.5*( SC->parameters[CellParams::RHOM_R+i] + SC->parameters[CellParams::RHOM_V_PREV+i] ); // at 0.5
               avgMoments3[i] = 0.5*( SC->parameters[CellParams::RHOM_V_PREV+i] + SC->parameters[CellParams::RHOM_R_PREV+i] ); // at 0.0
               avgMoments4[i] = 0.5*( SC->parameters[CellParams::RHOM_R_PREV+i] + SC->parameters[CellParams::RHOM_V_PREV_PREV+i] ); // at -0.5
               avgMoments5[i] = 0.5*( SC->parameters[CellParams::RHOM_V_PREV_PREV+i] + SC->parameters[CellParams::RHOM_R_PREV_PREV+i] ); // at -1.0
            }

            for (int i=0; i<nMomentsToInterp; i++) {
               switch(P::tcMomentInterpolationType) {
                  case 1:
                     if (normModul < 0.5) {
                        SC->parameters[cp_rhom+i] = linearInterpolation(0.0, avgMoments3[i], 0.5, avgMoments2[i], normModul);
                     } else {
                        SC->parameters[cp_rhom+i] = linearInterpolation(0.5, avgMoments2[i], 1.0, avgMoments1[i], normModul);
                     }
                     break;
                  case 2:
                     SC->parameters[cp_rhom+i] = lagrangeInterpolation2order(0.0, avgMoments3[i], 0.5, avgMoments2[i], 1.0, avgMoments1[i], normModul);
                     break;
                  case 3:
                     SC->parameters[cp_rhom+i] = lagrangeInterpolation3order(-0.5, avgMoments4[i], 0.0, avgMoments3[i], 0.5, avgMoments2[i], 1.0, avgMoments1[i], normModul);
                     break;
                  case -1:
                     SC->parameters[cp_rhom+i] = cubicHermiteSplineInterpolation(-0.5, avgMoments4[i], 0.0, avgMoments3[i], 0.5, avgMoments2[i], 1.0, avgMoments1[i], normModul);
                     break;
               }
            }
         }
      } // if timeclass = maxTC
   } // for loop over cells
}

void updateParticlePopulations(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid, const bool lastUpdateWasTranslation) {

   const vector<CellID>& cells = getLocalCells();
   #pragma omp parallel for
   for (size_t c=0; c<cells.size(); ++c) {

      const CellID cellID = cells[c];
      SpatialCell* SC = mpiGrid[cellID];
      const int timeclass = SC->parameters[CellParams::TIMECLASS];

      if (P::currentMaxTimeclass == 0 || timeclass == P::currentMaxTimeclass) {

         for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
            spatial_cell::Population& pop = SC->get_population(popID);
            pop.RHO = 0.5 * ( pop.RHO_R + pop.RHO_V );
            for(int i=0; i<3; i++) {
               pop.V[i] = 0.5 * ( pop.V_R[i] + pop.V_V[i] );
            }
            for(int i=0; i<6; i++) {
               pop.P[i] = 0.5 * ( pop.P_R[i] + pop.P_V[i] );
            }
         }

      } else {

         //interpolate output data
         Real RTCpow = pow(2, P::currentMaxTimeclass - timeclass);
         Real modul = P::fractionalTimestep % (int)RTCpow;
         Real normModul = modul/RTCpow; // point at which to interpolate to, between [0, 1]

         for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
            spatial_cell::Population& pop = SC->get_population(popID);

            std::array<Real, 5> avgMomentsRho; // logic of these lists is flipped in comparison to the other interpolation function
            std::array<std::array<Real, 3>, 5> avgMomentsV; 
            std::array<std::array<Real, 6>, 5> avgMomentsP;

            if (lastUpdateWasTranslation && SC->get_timeclass_turn_v()) { // aka if translation moments are ahead of acceleration moments
               // here, temporal order from newest to oldest is _R, _V, _R_PREV, _V_PREV, _R_PREV_PREV, _V_PREV_PREV

               avgMomentsRho.at(0) = 0.5 * ( pop.RHO_R + pop.RHO_V ); // 0.5
               avgMomentsRho.at(1) = 0.5 * ( pop.RHO_V + pop.RHO_R_PREV ); // 0.0
               avgMomentsRho.at(2) = 0.5 * ( pop.RHO_R_PREV + pop.RHO_V_PREV ); // -0.5
               avgMomentsRho.at(3) = 0.5 * ( pop.RHO_V_PREV + pop.RHO_R_PREV_PREV ); // -1.0
               avgMomentsRho.at(4) = 0.5 * ( pop.RHO_R_PREV_PREV + pop.RHO_V_PREV_PREV ); // -1.5

               for(int i=0; i<3; i++) {
                  avgMomentsV.at(0).at(i) = 0.5 * (pop.V_R[i] + pop.V_V[i]);
                  avgMomentsV.at(1).at(i) = 0.5 * (pop.V_V[i] + pop.V_R_PREV[i]);
                  avgMomentsV.at(2).at(i) = 0.5 * (pop.V_R_PREV[i] + pop.V_V_PREV[i]);
                  avgMomentsV.at(3).at(i) = 0.5 * (pop.V_V_PREV[i] + pop.V_R_PREV_PREV[i]);
                  avgMomentsV.at(4).at(i) = 0.5 * (pop.V_R_PREV_PREV[i] + pop.V_V_PREV_PREV[i]);
               }
               for(int i=0; i<6; i++) {
                  avgMomentsP.at(0).at(i) = 0.5 * (pop.P_R[i] + pop.P_V[i]);
                  avgMomentsP.at(1).at(i) = 0.5 * (pop.P_V[i] + pop.P_R_PREV[i]);
                  avgMomentsP.at(2).at(i) = 0.5 * (pop.P_R_PREV[i] + pop.P_V_PREV[i]);
                  avgMomentsP.at(3).at(i) = 0.5 * (pop.P_V_PREV[i] + pop.P_R_PREV_PREV[i]);
                  avgMomentsP.at(4).at(i) = 0.5 * (pop.P_R_PREV_PREV[i] + pop.P_V_PREV_PREV[i]);               
               }

               switch(P::tcMomentInterpolationType) {
                  case 1:
                     pop.RHO = linearInterpolation(0.0, avgMomentsRho.at(1), 0.5, avgMomentsRho.at(0), normModul);
                     for(int i=0; i<3; i++) {
                        pop.V[i] = linearInterpolation(0.0, avgMomentsV.at(1).at(i), 0.5, avgMomentsV.at(0).at(i), normModul);
                     }
                     for(int i=0; i<6; i++) {
                        pop.P[i] = linearInterpolation(0.0, avgMomentsP.at(1).at(i), 0.5, avgMomentsP.at(0).at(i), normModul);
                     }
                     break;
                  case 2:
                     pop.RHO = lagrangeInterpolation2order(-0.5, avgMomentsRho.at(2), 0.0, avgMomentsRho.at(1), 0.5, avgMomentsRho.at(0), normModul);
                     for(int i=0; i<3; i++) {
                        pop.V[i] = lagrangeInterpolation2order(-0.5, avgMomentsV.at(2).at(i), 0.0, avgMomentsV.at(1).at(i), 0.5, avgMomentsV.at(0).at(i), normModul);
                     }
                     for(int i=0; i<6; i++) {
                        pop.P[i] = lagrangeInterpolation2order(-0.5, avgMomentsP.at(2).at(i), 0.0, avgMomentsP.at(1).at(i), 0.5, avgMomentsP.at(0).at(i), normModul);
                     }
                     break;
                  case 3:
                     pop.RHO = lagrangeInterpolation3order(-1.0, avgMomentsRho.at(3), -0.5, avgMomentsRho.at(2), 0.0, avgMomentsRho.at(1), 0.5, avgMomentsRho.at(0), normModul);
                     for(int i=0; i<3; i++) {
                        pop.V[i] = lagrangeInterpolation3order(-1.0, avgMomentsV.at(3).at(i), -0.5, avgMomentsV.at(2).at(i), 0.0, avgMomentsV.at(1).at(i), 0.5, avgMomentsV.at(0).at(i), normModul);
                     }
                     for(int i=0; i<6; i++) {
                        pop.P[i] = lagrangeInterpolation3order(-1.0, avgMomentsP.at(3).at(i), -0.5, avgMomentsP.at(2).at(i), 0.0, avgMomentsP.at(1).at(i), 0.5, avgMomentsP.at(0).at(i), normModul);
                     }
                     break;
                  case -1:
                     pop.RHO = cubicHermiteSplineInterpolation(-1.0, avgMomentsRho.at(3), -0.5, avgMomentsRho.at(2), 0.0, avgMomentsRho.at(1), 0.5, avgMomentsRho.at(0), normModul);
                     for(int i=0; i<3; i++) {
                        pop.V[i] = cubicHermiteSplineInterpolation(-1.0, avgMomentsV.at(3).at(i), -0.5, avgMomentsV.at(2).at(i), 0.0, avgMomentsV.at(1).at(i), 0.5, avgMomentsV.at(0).at(i), normModul);
                     }
                     for(int i=0; i<6; i++) {
                        pop.P[i] = cubicHermiteSplineInterpolation(-1.0, avgMomentsP.at(3).at(i), -0.5, avgMomentsP.at(2).at(i), 0.0, avgMomentsP.at(1).at(i), 0.5, avgMomentsP.at(0).at(i), normModul);
                     }  
                     break;
               } 
   

            } else { // aka if acceleration moments are ahead of translation moments
               // here, temporal order from newest to oldest is _V, _R, _V_PREV, _R_PREV, _V_PREV_PREV, _R_PREV_PREV
               
               avgMomentsRho.at(0) = 0.5 * ( pop.RHO_V + pop.RHO_R ); // 1.0
               avgMomentsRho.at(1) = 0.5 * ( pop.RHO_R + pop.RHO_V_PREV ); // 0.5
               avgMomentsRho.at(2) = 0.5 * ( pop.RHO_V_PREV + pop.RHO_R_PREV ); // 0.0
               avgMomentsRho.at(3) = 0.5 * ( pop.RHO_R_PREV + pop.RHO_V_PREV_PREV ); // -0.5
               avgMomentsRho.at(4) = 0.5 * ( pop.RHO_V_PREV_PREV + pop.RHO_R_PREV_PREV ); // -1.0

               for(int i=0; i<3; i++) {
                  avgMomentsV.at(0).at(i) = 0.5 * (pop.V_V[i] + pop.V_R[i]);
                  avgMomentsV.at(1).at(i) = 0.5 * (pop.V_R[i] + pop.V_V_PREV[i]);
                  avgMomentsV.at(2).at(i) = 0.5 * (pop.V_V_PREV[i] + pop.V_R_PREV[i]);
                  avgMomentsV.at(3).at(i) = 0.5 * (pop.V_R_PREV[i] + pop.V_V_PREV_PREV[i]);
                  avgMomentsV.at(4).at(i) = 0.5 * (pop.V_V_PREV_PREV[i] + pop.V_R_PREV_PREV[i]);
               }
               for(int i=0; i<6; i++) {
                  avgMomentsP.at(0).at(i) = 0.5 * (pop.P_V[i] + pop.P_R[i]);
                  avgMomentsP.at(1).at(i) = 0.5 * (pop.P_R[i] + pop.P_V_PREV[i]);
                  avgMomentsP.at(2).at(i) = 0.5 * (pop.P_V_PREV[i] + pop.P_R_PREV[i]);
                  avgMomentsP.at(3).at(i) = 0.5 * (pop.P_R_PREV[i] + pop.P_V_PREV_PREV[i]);
                  avgMomentsP.at(4).at(i) = 0.5 * (pop.P_V_PREV_PREV[i] + pop.P_R_PREV_PREV[i]);
               }

               switch(P::tcMomentInterpolationType) {
                  case 1:
                     if (normModul < 0.5) {
                        pop.RHO = linearInterpolation(0.0, avgMomentsRho.at(2), 0.5, avgMomentsRho.at(1), normModul);
                     } else {
                        pop.RHO = linearInterpolation(0.5, avgMomentsRho.at(1), 1.0, avgMomentsRho.at(0), normModul);
                     }
                     for(int i=0; i<3; i++) {
                        if (normModul < 0.5) {
                           pop.V[i] = linearInterpolation(0.0, avgMomentsV.at(2).at(i), 0.5, avgMomentsV.at(1).at(i), normModul);
                        } else {
                           pop.V[i] = linearInterpolation(0.5, avgMomentsV.at(1).at(i), 1.0, avgMomentsV.at(0).at(i), normModul);
                        }
                     }
                     for(int i=0; i<6; i++) {
                        if (normModul < 0.5) {
                           pop.P[i] = linearInterpolation(0.0, avgMomentsP.at(2).at(i), 0.5, avgMomentsP.at(1).at(i), normModul);
                        } else {
                           pop.P[i] = linearInterpolation(0.5, avgMomentsP.at(1).at(i), 1.0, avgMomentsP.at(0).at(i), normModul);
                        }
                     }
                     break;
                  case 2:
                     pop.RHO = lagrangeInterpolation2order(0.0, avgMomentsRho.at(2), 0.5, avgMomentsRho.at(1), 1.0, avgMomentsRho.at(0), normModul);
                     for(int i=0; i<3; i++) {
                        pop.V[i] = lagrangeInterpolation2order(0.0, avgMomentsV.at(2).at(i), 0.5, avgMomentsV.at(1).at(i), 1.0, avgMomentsV.at(0).at(i), normModul);
                     }
                     for(int i=0; i<6; i++) {
                        pop.P[i] = lagrangeInterpolation2order(0.0, avgMomentsP.at(2).at(i), 0.5, avgMomentsP.at(1).at(i), 1.0, avgMomentsP.at(0).at(i), normModul);
                     }
                     break;
                  case 3:
                     pop.RHO = lagrangeInterpolation3order(-0.5, avgMomentsRho.at(3), 0.0, avgMomentsRho.at(2), 0.5, avgMomentsRho.at(1), 1.0, avgMomentsRho.at(0), normModul);
                     for(int i=0; i<3; i++) {
                        pop.V[i] = lagrangeInterpolation3order(-0.5, avgMomentsV.at(3).at(i), 0.0, avgMomentsV.at(2).at(i), 0.5, avgMomentsV.at(1).at(i), 1.0, avgMomentsV.at(0).at(i), normModul);
                     }
                     for(int i=0; i<6; i++) {
                        pop.P[i] = lagrangeInterpolation3order(-0.5, avgMomentsP.at(3).at(i), 0.0, avgMomentsP.at(2).at(i), 0.5, avgMomentsP.at(1).at(i), 1.0, avgMomentsP.at(0).at(i), normModul);
                     }
                     break;
                  case -1:
                     pop.RHO = cubicHermiteSplineInterpolation(-0.5, avgMomentsRho.at(3), 0.0, avgMomentsRho.at(2), 0.5, avgMomentsRho.at(1), 1.0, avgMomentsRho.at(0), normModul);
                     for(int i=0; i<3; i++) {
                        pop.V[i] = cubicHermiteSplineInterpolation(-0.5, avgMomentsV.at(3).at(i), 0.0, avgMomentsV.at(2).at(i), 0.5, avgMomentsV.at(1).at(i), 1.0, avgMomentsV.at(0).at(i), normModul);
                     }
                     for(int i=0; i<6; i++) {
                        pop.P[i] = cubicHermiteSplineInterpolation(-0.5, avgMomentsP.at(3).at(i), 0.0, avgMomentsP.at(2).at(i), 0.5, avgMomentsP.at(1).at(i), 1.0, avgMomentsP.at(0).at(i), normModul);
                     }
                     break;
               }
            } // else for acceleration moments ahead of translation moments
         } // for loop over populations
      } // if using timeclasses
   } // for loop over cells 
}


void calculateInitialVelocityMoments(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid) {
   const vector<CellID>& cells = getLocalCells();
   phiprof::Timer timer {"Calculate moments"};

   // Iterate through all local cells (incl. system boundary cells):
   // Setting the GPU device inside the moment call itself, because
   // it's being called from so many different projects etc
   #pragma omp parallel for
   for (size_t c=0; c<cells.size(); ++c) {
      const CellID cellID = cells[c];
      SpatialCell* SC = mpiGrid[cellID];
      calculateCellMoments(SC,true,false);
      // WARNING the following is sane as this function is only called by initializeGrid.
      // We need initialized _DT2 values for the dt=0 field propagation done in the beginning.
      // Later these will be set properly.
      SC->parameters[CellParams::RHOM_DT2] = SC->parameters[CellParams::RHOM];
      SC->parameters[CellParams::VX_DT2] = SC->parameters[CellParams::VX];
      SC->parameters[CellParams::VY_DT2] = SC->parameters[CellParams::VY];
      SC->parameters[CellParams::VZ_DT2] = SC->parameters[CellParams::VZ];
      SC->parameters[CellParams::RHOQ_DT2] = SC->parameters[CellParams::RHOQ];
      SC->parameters[CellParams::P_11_DT2] = SC->parameters[CellParams::P_11];
      SC->parameters[CellParams::P_22_DT2] = SC->parameters[CellParams::P_22];
      SC->parameters[CellParams::P_33_DT2] = SC->parameters[CellParams::P_33];
      SC->parameters[CellParams::P_23_DT2] = SC->parameters[CellParams::P_23];
      SC->parameters[CellParams::P_13_DT2] = SC->parameters[CellParams::P_13];
      SC->parameters[CellParams::P_12_DT2] = SC->parameters[CellParams::P_12];
   } // for-loop over spatial cells
}

void updatePreviousVMoments(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid, bool isInitialization) {

   const vector<CellID>& cells = getLocalCells();
   for (size_t c=0; c<cells.size(); c++) {
      const CellID cellID = cells[c];
      SpatialCell* SC = mpiGrid[cellID];
      Real tdiff = P::timeclassDt[SC->parameters[CellParams::TIMECLASS]];

      if (isInitialization == true) {
         // initializiing all _PREV and _PREV_PREV moments to same values as _V moments
         // this block only gets executed once, when the simulation is started
         SC->parameters[CellParams::RHOM_V_PREV_PREV] = SC->parameters[CellParams::RHOM_V];
         SC->parameters[CellParams::RHOQ_V_PREV_PREV] = SC->parameters[CellParams::RHOQ_V];
         SC->parameters[CellParams::P_11_V_PREV_PREV] = SC->parameters[CellParams::P_11_V];
         SC->parameters[CellParams::P_22_V_PREV_PREV] = SC->parameters[CellParams::P_22_V];
         SC->parameters[CellParams::P_33_V_PREV_PREV] = SC->parameters[CellParams::P_33_V];
         SC->parameters[CellParams::P_23_V_PREV_PREV] = SC->parameters[CellParams::P_23_V];
         SC->parameters[CellParams::P_13_V_PREV_PREV] = SC->parameters[CellParams::P_13_V];
         SC->parameters[CellParams::P_12_V_PREV_PREV] = SC->parameters[CellParams::P_12_V];
         

         SC->parameters[CellParams::RHOM_V_PREV] = SC->parameters[CellParams::RHOM_V];
         SC->parameters[CellParams::RHOQ_V_PREV] = SC->parameters[CellParams::RHOQ_V];
         SC->parameters[CellParams::P_11_V_PREV] = SC->parameters[CellParams::P_11_V];
         SC->parameters[CellParams::P_22_V_PREV] = SC->parameters[CellParams::P_22_V];
         SC->parameters[CellParams::P_33_V_PREV] = SC->parameters[CellParams::P_33_V];
         SC->parameters[CellParams::P_23_V_PREV] = SC->parameters[CellParams::P_23_V];
         SC->parameters[CellParams::P_13_V_PREV] = SC->parameters[CellParams::P_13_V];
         SC->parameters[CellParams::P_12_V_PREV] = SC->parameters[CellParams::P_12_V];

         SC->parameters[CellParams::RHOM_R_PREV_PREV] = SC->parameters[CellParams::RHOM_R];
         SC->parameters[CellParams::RHOQ_R_PREV_PREV] = SC->parameters[CellParams::RHOQ_R];
         SC->parameters[CellParams::P_11_R_PREV_PREV] = SC->parameters[CellParams::P_11_R];
         SC->parameters[CellParams::P_22_R_PREV_PREV] = SC->parameters[CellParams::P_22_R];
         SC->parameters[CellParams::P_33_R_PREV_PREV] = SC->parameters[CellParams::P_33_R];
         SC->parameters[CellParams::P_23_R_PREV_PREV] = SC->parameters[CellParams::P_23_R];
         SC->parameters[CellParams::P_13_R_PREV_PREV] = SC->parameters[CellParams::P_13_R];
         SC->parameters[CellParams::P_12_R_PREV_PREV] = SC->parameters[CellParams::P_12_R];

         SC->parameters[CellParams::RHOM_R_PREV] = SC->parameters[CellParams::RHOM_R];
         SC->parameters[CellParams::RHOQ_R_PREV] = SC->parameters[CellParams::RHOQ_R];
         SC->parameters[CellParams::P_11_R_PREV] = SC->parameters[CellParams::P_11_R];
         SC->parameters[CellParams::P_22_R_PREV] = SC->parameters[CellParams::P_22_R];
         SC->parameters[CellParams::P_33_R_PREV] = SC->parameters[CellParams::P_33_R];
         SC->parameters[CellParams::P_23_R_PREV] = SC->parameters[CellParams::P_23_R];
         SC->parameters[CellParams::P_13_R_PREV] = SC->parameters[CellParams::P_13_R];
         SC->parameters[CellParams::P_12_R_PREV] = SC->parameters[CellParams::P_12_R];

         for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
            spatial_cell::Population& pop = SC->get_population(popID);
            pop.RHO_R_PREV = pop.RHO_R;
            pop.RHO_V_PREV = pop.RHO_V;
            pop.RHO_R_PREV_PREV = pop.RHO_R;
            pop.RHO_V_PREV_PREV = pop.RHO_V;
            for(int i=0; i<3; i++) {
               pop.V_R_PREV[i] = pop.V_R[i];
               pop.V_V_PREV[i] = pop.V_V[i];
               pop.V_R_PREV_PREV[i] = pop.V_R[i];
               pop.V_V_PREV_PREV[i] = pop.V_V[i];
            }
            for(int i=0; i<6; i++) {
               pop.P_R_PREV[i] = pop.P_R[i];
               pop.P_V_PREV[i] = pop.P_V[i];
               pop.P_R_PREV_PREV[i] = pop.P_R[i];
               pop.P_V_PREV_PREV[i] = pop.P_V[i];
            }
         }

         if (false) {
                  
            Eigen::Transform<Real,3,Eigen::Affine> vUpdateMatrixRP = compute_acceleration_transformation(SC, 0, tdiff);
            Eigen::Transform<Real,3,Eigen::Affine> vUpdateMatrixRPP = compute_acceleration_transformation(SC, 0, tdiff*2);
            Eigen::Transform<Real,3,Eigen::Affine> vUpdateMatrixVP = compute_acceleration_transformation(SC, 0, tdiff);
            Eigen::Transform<Real,3,Eigen::Affine> vUpdateMatrixVPP = compute_acceleration_transformation(SC, 0, tdiff*2);

            auto bwd_RP = vUpdateMatrixRP.inverse();
            auto bwd_RPP = vUpdateMatrixRPP.inverse();
            auto bwd_VP = vUpdateMatrixVP.inverse();
            auto bwd_VPP = vUpdateMatrixVPP.inverse();

            Eigen::Matrix<Real,3,1> V_R_P(SC->parameters[CellParams::VX_R_PREV], SC->parameters[CellParams::VY_R_PREV], SC->parameters[CellParams::VZ_R_PREV]);
            Eigen::Matrix<Real,3,1> V_R_PP(SC->parameters[CellParams::VX_R_PREV_PREV], SC->parameters[CellParams::VY_R_PREV_PREV], SC->parameters[CellParams::VZ_R_PREV_PREV]);
            Eigen::Matrix<Real,3,1> V_V_P(SC->parameters[CellParams::VX_V_PREV], SC->parameters[CellParams::VY_V_PREV], SC->parameters[CellParams::VZ_V_PREV]);
            Eigen::Matrix<Real,3,1> V_V_PP(SC->parameters[CellParams::VX_V_PREV_PREV], SC->parameters[CellParams::VY_V_PREV_PREV], SC->parameters[CellParams::VZ_V_PREV_PREV]);

            Eigen::Matrix<Real,3,1> V_R_P_updated = bwd_RP * V_R_P;
            Eigen::Matrix<Real,3,1> V_R_PP_updated = bwd_RPP * V_R_PP;
            Eigen::Matrix<Real,3,1> V_V_P_updated = bwd_VP * V_V_P;
            Eigen::Matrix<Real,3,1> V_V_PP_updated = bwd_VPP * V_V_PP;

            SC->parameters[CellParams::VX_R_PREV] = V_R_P_updated(0);
            SC->parameters[CellParams::VY_R_PREV] = V_R_P_updated(1);
            SC->parameters[CellParams::VZ_R_PREV] = V_R_P_updated(2);

            SC->parameters[CellParams::VX_R_PREV_PREV] = V_R_PP_updated(0);
            SC->parameters[CellParams::VY_R_PREV_PREV] = V_R_PP_updated(1);
            SC->parameters[CellParams::VZ_R_PREV_PREV] = V_R_PP_updated(2);

            SC->parameters[CellParams::VX_V_PREV] = V_V_P_updated(0);
            SC->parameters[CellParams::VY_V_PREV] = V_V_P_updated(1);
            SC->parameters[CellParams::VZ_V_PREV] = V_V_P_updated(2);

            SC->parameters[CellParams::VX_V_PREV_PREV] = V_V_PP_updated(0);
            SC->parameters[CellParams::VY_V_PREV_PREV] = V_V_PP_updated(1);
            SC->parameters[CellParams::VZ_V_PREV_PREV] = V_V_PP_updated(2);

         } else {

            SC->parameters[CellParams::VX_V_PREV_PREV] = SC->parameters[CellParams::VX_V];
            SC->parameters[CellParams::VY_V_PREV_PREV] = SC->parameters[CellParams::VY_V];
            SC->parameters[CellParams::VZ_V_PREV_PREV] = SC->parameters[CellParams::VZ_V];
            
            SC->parameters[CellParams::VX_V_PREV] = SC->parameters[CellParams::VX_V];
            SC->parameters[CellParams::VY_V_PREV] = SC->parameters[CellParams::VY_V];
            SC->parameters[CellParams::VZ_V_PREV] = SC->parameters[CellParams::VZ_V];

            SC->parameters[CellParams::VX_R_PREV_PREV] = SC->parameters[CellParams::VX_R];
            SC->parameters[CellParams::VY_R_PREV_PREV] = SC->parameters[CellParams::VY_R];
            SC->parameters[CellParams::VZ_R_PREV_PREV] = SC->parameters[CellParams::VZ_R];

            SC->parameters[CellParams::VX_R_PREV] = SC->parameters[CellParams::VX_R];
            SC->parameters[CellParams::VY_R_PREV] = SC->parameters[CellParams::VY_R];
            SC->parameters[CellParams::VZ_R_PREV] = SC->parameters[CellParams::VZ_R];
         }

      } else {
         if (SC->get_timeclass_turn_v() == true) {
            // updating _PREV_PREV moments
            SC->parameters[CellParams::RHOM_V_PREV_PREV] = SC->parameters[CellParams::RHOM_V_PREV];
            SC->parameters[CellParams::VX_V_PREV_PREV] = SC->parameters[CellParams::VX_V_PREV];
            SC->parameters[CellParams::VY_V_PREV_PREV] = SC->parameters[CellParams::VY_V_PREV];
            SC->parameters[CellParams::VZ_V_PREV_PREV] = SC->parameters[CellParams::VZ_V_PREV];
            SC->parameters[CellParams::RHOQ_V_PREV_PREV] = SC->parameters[CellParams::RHOQ_V_PREV];
            SC->parameters[CellParams::P_11_V_PREV_PREV] = SC->parameters[CellParams::P_11_V_PREV];
            SC->parameters[CellParams::P_22_V_PREV_PREV] = SC->parameters[CellParams::P_22_V_PREV];
            SC->parameters[CellParams::P_33_V_PREV_PREV] = SC->parameters[CellParams::P_33_V_PREV];
            SC->parameters[CellParams::P_23_V_PREV_PREV] = SC->parameters[CellParams::P_23_V_PREV];
            SC->parameters[CellParams::P_13_V_PREV_PREV] = SC->parameters[CellParams::P_13_V_PREV];
            SC->parameters[CellParams::P_12_V_PREV_PREV] = SC->parameters[CellParams::P_12_V_PREV];

            // updating _PREV moments
            SC->parameters[CellParams::RHOM_V_PREV] = SC->parameters[CellParams::RHOM_V];
            SC->parameters[CellParams::VX_V_PREV] = SC->parameters[CellParams::VX_V];
            SC->parameters[CellParams::VY_V_PREV] = SC->parameters[CellParams::VY_V];
            SC->parameters[CellParams::VZ_V_PREV] = SC->parameters[CellParams::VZ_V];
            SC->parameters[CellParams::RHOQ_V_PREV] = SC->parameters[CellParams::RHOQ_V];
            SC->parameters[CellParams::P_11_V_PREV] = SC->parameters[CellParams::P_11_V];
            SC->parameters[CellParams::P_22_V_PREV] = SC->parameters[CellParams::P_22_V];
            SC->parameters[CellParams::P_33_V_PREV] = SC->parameters[CellParams::P_33_V];
            SC->parameters[CellParams::P_23_V_PREV] = SC->parameters[CellParams::P_23_V];
            SC->parameters[CellParams::P_13_V_PREV] = SC->parameters[CellParams::P_13_V];
            SC->parameters[CellParams::P_12_V_PREV] = SC->parameters[CellParams::P_12_V];
         }
      }
   }  
}
