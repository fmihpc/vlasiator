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

#include <cstdlib>
#include <iostream>
#include <vector>
#include <stdint.h>

#ifdef _OPENMP
   #include <omp.h>
#endif
#include <zoltan.h>
#include <dccrg.hpp>
#include <phiprof.hpp>

#include "../spatial_cell_wrapper.hpp"
#include "../vlasovmover.h"
#include "../grid.h"
#include "../definitions.h"
#include "../object_wrapper.h"
#include "../mpiconversion.h"

#include "cpu_moments.h"
#include "cpu_acc_semilag.hpp"
#include "cpu_trans_map.hpp"
#include "cpu_trans_pencils.hpp"
#include "cpu_trans_map_amr.hpp"
#include "cpu_acc_transform.hpp"

#include "spline.h"

using namespace std;
using namespace spatial_cell;

creal ZERO    = 0.0;
creal HALF    = 0.5;
creal FOURTH  = 1.0/4.0;
creal SIXTH   = 1.0/6.0;
creal ONE     = 1.0;
creal TWO     = 2.0;
creal EPSILON = 1.0e-25;

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
        creal dt,
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
      //updateRemoteVelocityBlockLists(mpiGrid,popID,VLASOV_SOLVER_Z_NEIGHBORHOOD_ID);
      SpatialCell::set_mpi_transfer_type(Transfer::VEL_BLOCK_DATA,false);
      mpiGrid.update_copies_of_remote_neighbors(VLASOV_SOLVER_Z_NEIGHBORHOOD_ID);
      transTimer.stop();

      // bt=phiprof::initializeTimer("barrier-trans-pre-trans_map_1d-z","Barriers","MPI");
      // phiprof::start(bt);
      // MPI_Barrier(MPI_COMM_WORLD);
      // phiprof::stop(bt);

      t1 = MPI_Wtime();
      phiprof::Timer computeTimer {"compute-mapping-z"};
      if(P::amrMaxSpatialRefLevel == 0) {
         trans_map_1d(mpiGrid,local_propagated_cells, remoteTargetCellsz, 2, dt,popID); // map along z//
      } else {
         trans_map_1d_amr(mpiGrid,local_propagated_cells, remoteTargetCellsz, nPencils, 2, dt, 0, popID); // map along z//
      }
      computeTimer.stop();
      time += MPI_Wtime() - t1;

      phiprof::Timer btTimer {"barrier-trans-pre-update_remote-z", {"Barriers","MPI"}};
      MPI_Barrier(MPI_COMM_WORLD);
      btTimer.stop();

      phiprof::Timer updateRemoteTimer {"update_remote-z", {"MPI"}};
      if(P::amrMaxSpatialRefLevel == 0) {
         update_remote_mapping_contribution(mpiGrid, 2,+1,popID);
         update_remote_mapping_contribution(mpiGrid, 2,-1,popID);
      } else {
         update_remote_mapping_contribution_amr(mpiGrid, 2,+1,popID);
         update_remote_mapping_contribution_amr(mpiGrid, 2,-1,popID);
      }
      updateRemoteTimer.stop();

   }

   phiprof::Timer btxTimer {"barrier-trans-pre-x", {"Barriers","MPI"}};
   MPI_Barrier(MPI_COMM_WORLD);
   btxTimer.stop();
   
   // ------------- SLICE - map dist function in X --------------- //
   if(P::xcells_ini > 1){
      
      phiprof::Timer transTimer {"transfer-stencil-data-x", {"MPI"}};
      //updateRemoteVelocityBlockLists(mpiGrid,popID,VLASOV_SOLVER_X_NEIGHBORHOOD_ID);
      SpatialCell::set_mpi_transfer_type(Transfer::VEL_BLOCK_DATA,false);
      mpiGrid.update_copies_of_remote_neighbors(VLASOV_SOLVER_X_NEIGHBORHOOD_ID);
      transTimer.stop();
      
      // bt=phiprof::initializeTimer("barrier-trans-pre-trans_map_1d-x","Barriers","MPI");
      // phiprof::start(bt);
      // MPI_Barrier(MPI_COMM_WORLD);
      // phiprof::stop(bt);

      t1 = MPI_Wtime();
      phiprof::Timer computeTimer {"compute-mapping-x"};
      if(P::amrMaxSpatialRefLevel == 0) {
         trans_map_1d(mpiGrid,local_propagated_cells, remoteTargetCellsx, 0,dt,popID); // map along x//
      } else {
         trans_map_1d_amr(mpiGrid,local_propagated_cells, remoteTargetCellsx, nPencils, 0, dt, 0, popID); // map along x//
      }
      computeTimer.stop();
      time += MPI_Wtime() - t1;

      phiprof::Timer btTimer {"barrier-trans-pre-update_remote-x", {"Barriers","MPI"}};
      MPI_Barrier(MPI_COMM_WORLD);
      btTimer.stop();

      phiprof::Timer updateRemoteTimer {"update_remote-x", {"MPI"}};
      if(P::amrMaxSpatialRefLevel == 0) {
         update_remote_mapping_contribution(mpiGrid, 0,+1,popID);
         update_remote_mapping_contribution(mpiGrid, 0,-1,popID);
      } else {
         update_remote_mapping_contribution_amr(mpiGrid, 0,+1,popID);
         update_remote_mapping_contribution_amr(mpiGrid, 0,-1,popID);
      }
      updateRemoteTimer.stop();
   }

   phiprof::Timer btyTimer {"barrier-trans-pre-y", {"Barriers","MPI"}};
   MPI_Barrier(MPI_COMM_WORLD);
   btyTimer.stop();

   // ------------- SLICE - map dist function in Y --------------- //
   if(P::ycells_ini > 1) {
      
      phiprof::Timer transTimer {"transfer-stencil-data-y", {"MPI"}};
      //updateRemoteVelocityBlockLists(mpiGrid,popID,VLASOV_SOLVER_Y_NEIGHBORHOOD_ID);
      SpatialCell::set_mpi_transfer_type(Transfer::VEL_BLOCK_DATA,false);
      mpiGrid.update_copies_of_remote_neighbors(VLASOV_SOLVER_Y_NEIGHBORHOOD_ID);
      transTimer.stop();
      
      // bt=phiprof::initializeTimer("barrier-trans-pre-trans_map_1d-y","Barriers","MPI");
      // phiprof::start(bt);
      // MPI_Barrier(MPI_COMM_WORLD);
      // phiprof::stop(bt);

      t1 = MPI_Wtime();
      phiprof::Timer computeTimer {"compute-mapping-y"};
      if(P::amrMaxSpatialRefLevel == 0) {
         trans_map_1d(mpiGrid,local_propagated_cells, remoteTargetCellsy, 1,dt,popID); // map along y//
      } else {
         trans_map_1d_amr(mpiGrid,local_propagated_cells, remoteTargetCellsy, nPencils, 1, dt, 0, popID); // map along y//
      }
      computeTimer.stop();
      time += MPI_Wtime() - t1;
      
      phiprof::Timer btTimer {"barrier-trans-pre-update_remote-y", {"Barriers","MPI"}};
      MPI_Barrier(MPI_COMM_WORLD);
      btTimer.stop();

      phiprof::Timer updateRemoteTimer {"update_remote-y", {"MPI"}};
      if(P::amrMaxSpatialRefLevel == 0) {
         update_remote_mapping_contribution(mpiGrid, 1,+1,popID);
         update_remote_mapping_contribution(mpiGrid, 1,-1,popID);
      } else {
         update_remote_mapping_contribution_amr(mpiGrid, 1,+1,popID);
         update_remote_mapping_contribution_amr(mpiGrid, 1,-1,popID);
      }
      updateRemoteTimer.stop();
   }

   for(CellID c : local_propagated_cells)
   {
      // if (c == 16) std::cout << c << " at TIME_R " << mpiGrid[c]->parameters[CellParams::TIME_R] << " + " << dt <<"\n";
      mpiGrid[c]->parameters[CellParams::TIME_R] += dt;
   }
   phiprof::Timer btpostimer {"barrier-trans-post-trans",{"Barriers","MPI"}};
   MPI_Barrier(MPI_COMM_WORLD);
   btpostimer.stop();

   // MPI_Barrier(MPI_COMM_WORLD);
   // bailout(true, "", __FILE__, __LINE__);
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
   creal dt,
   const uint popID,
   Real &time,
   int tc
   ) {

   // Ghost translation, need all cell information, not just for a single direction.
   // No need for remote target cells; pass a dummy list.
   const vector<CellID> dummy_cells;

   updateRemoteVelocityBlockLists(mpiGrid,popID,VLASOV_SOLVER_GHOST_NEIGHBORHOOD_ID);
   // Need to re-do in case block lists of boundary cells change after
   // the block adjustment just after ACC.

   phiprof::Timer prepreBarrierTimer {"MPI barrier-pre-trans-comm"};
   MPI_Barrier(MPI_COMM_WORLD);
   prepreBarrierTimer.stop();

   phiprof::Timer transferTimer {"transfer-stencil-data-all",{"MPI"}};
   SpatialCell::set_mpi_transfer_type(Transfer::VEL_BLOCK_DATA,false);
   mpiGrid.update_copies_of_remote_neighbors(VLASOV_SOLVER_GHOST_NEIGHBORHOOD_ID);
   transferTimer.stop();

   phiprof::Timer preBarrierTimer {"MPI barrier-pre-trans"};
   MPI_Barrier(MPI_COMM_WORLD);
   preBarrierTimer.stop();

   //#warning TODO: Implement also 2D / non-AMR ghost translation?
   // ------------- SLICE - map dist function in Z --------------- //
   phiprof::Timer mappingZTimer {"compute-mapping-z"};
   trans_map_1d_amr(mpiGrid,local_propagated_cells, dummy_cells, nPencils, 2, dt, tc, popID); // map along z//
   mappingZTimer.stop();

   // ------------- SLICE - map dist function in X --------------- //
   phiprof::Timer mappingXTimer {"compute-mapping-x"};
   trans_map_1d_amr(mpiGrid,local_propagated_cells, dummy_cells, nPencils, 0,dt, tc, popID); // map along x//
   mappingXTimer.stop();

   // ------------- SLICE - map dist function in Y --------------- //
   phiprof::Timer mappingYTimer {"compute-mapping-y"};
   trans_map_1d_amr(mpiGrid,local_propagated_cells, dummy_cells, nPencils, 1,dt, tc, popID); // map along y//
   mappingYTimer.stop();

   phiprof::Timer postBarrierTimer {"MPI barrier-post-trans"};
   MPI_Barrier(MPI_COMM_WORLD);
   postBarrierTimer.stop();
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
   std::cerr << std::scientific << "calculateSpatialTranslation at t="<<P::t << "\n";//", for dtfactor="<<dt<<"\n";
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

   vector<uint> nPencils;
   Real time=0.0;

   // If dt=0 we are either initializing or distribution functions are not translated.
   if (dt == 0.0 && initializationOrLB == true) {
      // if dt=0.0 and initialization == true, we are either initializing or load balancing.
      // hence we can just calculate the moments and return.
      calculateMoments_R(mpiGrid,localCells,true);
      return;
   }

   // TC propagation lists, TODO move out of here somewhere sensible and less often called
   for (int tc = 0; tc <= P::maxTimeclass; tc++)
   {
      // std::cout << "initing up tc " << tc << " vectors \n";
      tc_propagated_cells.push_back(vector<CellID>());
      tc_target_cells.push_back(vector<CellID>());
      tc_propagated_cell_sets.push_back(set<CellID>());
      tc_target_cell_sets.push_back(set<CellID>());
   }
   
   phiprof::Timer computeTimer {"compute_cell_lists"};
   if (!P::vlasovSolverGhostTranslate) {
      remoteTargetCellsx = mpiGrid.get_remote_cells_on_process_boundary(VLASOV_SOLVER_TARGET_X_NEIGHBORHOOD_ID);
      remoteTargetCellsy = mpiGrid.get_remote_cells_on_process_boundary(VLASOV_SOLVER_TARGET_Y_NEIGHBORHOOD_ID);
      remoteTargetCellsz = mpiGrid.get_remote_cells_on_process_boundary(VLASOV_SOLVER_TARGET_Z_NEIGHBORHOOD_ID);
   }

   // std::cout << "setting up tc cell vectors \n";
   // Figure out which spatial cells are translated,
   // result independent of particle species.
   // If performing ghost translation, this is used for LB.
   for (size_t c=0; c<localCells.size(); ++c) {
      int cellTC = (int)mpiGrid[localCells[c]]->parameters[CellParams::TIMECLASS];
      
      if (do_translate_cell(mpiGrid[localCells[c]],-1)) {
         tc_propagated_cell_sets[cellTC].insert(localCells[c]);
      }
      for (auto i: mpiGrid[localCells[c]]->requested_timeclass_ghosts){
         tc_propagated_cell_sets[i].insert(localCells[c]);
      }
      if (do_translate_cell(mpiGrid[localCells[c]],-1) &&
          mpiGrid[localCells[c]]->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
            tc_target_cell_sets[cellTC].insert(localCells[c]);
      }
      for (auto i: mpiGrid[localCells[c]]->requested_timeclass_ghosts){
         if(mpiGrid[localCells[c]]->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY)
            tc_propagated_cell_sets[i].insert(localCells[c]);
      }
   }


   // TC propagation lists, TODO move out of here somewhere sensible and less often called
   for (int tc = 0; tc <= P::maxTimeclass; tc++)
   {
      // std::cout << "initing up tc " << tc << " vectors \n";
      tc_propagated_cells[tc] = vector<CellID>(tc_propagated_cell_sets[tc].begin(),tc_propagated_cell_sets[tc].end());
      tc_target_cells[tc] = vector<CellID>(tc_target_cell_sets[tc].begin(),tc_target_cell_sets[tc].end());
   }

   // for (int tc = 0; tc <= P::maxTimeclass; tc++)
   // {
   //    std::cout << "" << tc << " cellvector size "<< tc_propagated_cells[tc].size()<<" \n";
      
   // }
   // Figure out which spatial cells are translated,
   // result independent of particle species.
   // If performing ghost translation, this is used for LB.
   for (size_t c=0; c<localCells.size(); ++c) {
      if (do_translate_cell(mpiGrid[localCells[c]],-1)) {
         local_propagated_cells.push_back(localCells[c]);
      }
   }

   if (P::prepareForRebalance == true && P::amrMaxSpatialRefLevel != 0) {
      // One more element to count the sums
      nPencils.resize(local_propagated_cells.size()+1, 0);
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

               SC->parameters[CellParams::RHOM_R_PREV] = SC->parameters[CellParams::RHOM_R];
               SC->parameters[CellParams::VX_R_PREV] = SC->parameters[CellParams::VX_R];
               SC->parameters[CellParams::VY_R_PREV] = SC->parameters[CellParams::VY_R];
               SC->parameters[CellParams::VZ_R_PREV] = SC->parameters[CellParams::VZ_R];
               SC->parameters[CellParams::RHOQ_R_PREV] = SC->parameters[CellParams::RHOQ_R];
               SC->parameters[CellParams::P_11_R_PREV] = SC->parameters[CellParams::P_11_R];
               SC->parameters[CellParams::P_22_R_PREV] = SC->parameters[CellParams::P_22_R];
               SC->parameters[CellParams::P_33_R_PREV] = SC->parameters[CellParams::P_33_R];
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

   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
   // Translate all particle species
   for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
      string profName = "translate "+getObjectWrapper().particleSpecies[popID].name;
      phiprof::Timer timer {profName};
      SpatialCell::setCommunicatedSpecies(popID);
      for(int tc = 0; tc <= P::currentMaxTimeclass; tc++){
         
         int mod = 1 << (P::currentMaxTimeclass - tc);
         if((P::fractionalTimestep % mod) == 0){
            
            std::cout << "rank " << myRank << ": " << tc_propagated_cells[tc].size() << " cells: calculateSpatialTranslation tc " << tc << " by dt " << P::timeclassDt[tc] <<"\n";
              if (P::vlasovSolverGhostTranslate && (P::amrMaxSpatialRefLevel > 0) ) {
               // Local translation without interim communication
               // Not yet implemented for non-AMR solver
               calculateSpatialGhostTranslation(
                  mpiGrid,
                  tc_propagated_cells[tc], // Used for LB
                  nPencils,
                  P::timeclassDt[tc],
                  popID,
                  time,
                  tc
                  );
            } else {
               calculateSpatialTranslation(
                  mpiGrid,
                  tc_propagated_cells[tc], //local_propagated_cells,
                  remoteTargetCellsx,
                  remoteTargetCellsy,
                  remoteTargetCellsz,
                  nPencils,
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
            counter += SC->get_number_of_velocity_blocks(popID);
         }

         // int accelerationsteps = 0; // Account for time spent in acceleration as well
         // if (mpiGrid[local_propagated_cells[c]]->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) accelerationsteps = 3;
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

            SC->parameters[CellParams::RHOM_R_PREV] = SC->parameters[CellParams::RHOM_R];
            SC->parameters[CellParams::VX_R_PREV] = SC->parameters[CellParams::VX_R];
            SC->parameters[CellParams::VY_R_PREV] = SC->parameters[CellParams::VY_R];
            SC->parameters[CellParams::VZ_R_PREV] = SC->parameters[CellParams::VZ_R];
            SC->parameters[CellParams::RHOQ_R_PREV] = SC->parameters[CellParams::RHOQ_R];
            SC->parameters[CellParams::P_11_R_PREV] = SC->parameters[CellParams::P_11_R];
            SC->parameters[CellParams::P_22_R_PREV] = SC->parameters[CellParams::P_22_R];
            SC->parameters[CellParams::P_33_R_PREV] = SC->parameters[CellParams::P_33_R];
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
}

/*
  --------------------------------------------------
  Acceleration (velocity space propagation)
  --------------------------------------------------
*/

/** Accelerate the given population to new time t+dt.
 * This function is AMR safe.
 * @param popID Particle population ID.
 * @param globalMaxSubcycles Number of times acceleration is subcycled.
 * @param step The current subcycle step.
 * @param mpiGrid Parallel grid library.
 * @param propagatedCells List of cells in which the population is accelerated.
 * @param dt Timestep factor.*/
void calculateAcceleration(const uint popID,const uint globalMaxSubcycles,const uint step,
                           dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                           const std::vector<CellID>& propagatedCells,
                           const Real& dt) {
   // Set active population
   SpatialCell::setCommunicatedSpecies(popID);

   // Calculate velocity moments, these are needed to
   // calculate the transforms used in the accelerations.
   // Calculated moments are stored in the "_V" variables.
   calculateMoments_V(mpiGrid, propagatedCells, false);

   int timerId {phiprof::initializeTimer("cell-semilag-acc")};

   // Semi-Lagrangian acceleration for those cells which are subcycled
   std::cout << "Propagating (ACC) cells ";
   #pragma omp parallel for schedule(dynamic,1)
   for (size_t c=0; c<propagatedCells.size(); ++c) {
      const CellID cellID = propagatedCells[c];
      std::cout << cellID << " ";
   }
   std::cout << "\n";
   for (size_t c=0; c<propagatedCells.size(); ++c) {
      const CellID cellID = propagatedCells[c];
      const Real maxVdt = mpiGrid[cellID]->get_max_v_dt(popID);//*mpiGrid[cellID]->parameters[CellParams::TIMECLASSDT];
      Real celldt = dt*mpiGrid[cellID]->get_tc_dt();

      
      //compute subcycle dt. The length is maxVdt on all steps
      //except the last one. This is to keep the neighboring
      //spatial cells in sync, so that two neighboring cells with
      //different number of subcycles have similar timesteps,
      //except that one takes an additional short step. This keeps
      //spatial block neighbors as much in sync as possible for
      //adjust blocks.
      Real subcycleDt;
      if( (step + 1) * maxVdt > fabs(celldt)) {
	      subcycleDt = max(fabs(celldt) - step * maxVdt, 0.0);
      } else{
         subcycleDt = maxVdt;
      }
      if (dt<0) subcycleDt = -subcycleDt;

      //generate pseudo-random order which is always the same irrespective of parallelization, restarts, etc.
      std::default_random_engine rndState;
      // set seed, initialise generator and get value. The order is the same
      // for all cells, but varies with timestep.
      rndState.seed(P::tstep + P::fractionalTimestep); // WARNING this formulation actually has some correlations (P::tstep + P::fractionalTimestep can do aliasing...)

      #ifndef DEBUG_TIMECLASSES
         uint map_order=std::uniform_int_distribution<>(0,2)(rndState);
      #else
         uint map_order=1;
      #endif
      map_order=1;
      phiprof::Timer semilagAccTimer {"cell-semilag-acc"};
      cpu_accelerate_cell(mpiGrid[cellID],popID,map_order,subcycleDt, 0);
      // if (cellID == 16)
      //    #pragma omp critical(output)
      //    {
      //       std::cout<< "Cellid " << cellID << " at t=" << mpiGrid[cellID]->parameters[CellParams::TIME_V] <<": subcycledt " << subcycleDt << " maxvdt "<< maxVdt << " step " << step << " globalmax " <<  globalMaxSubcycles << " step: " << step << "\n";
      //    }

      mpiGrid[cellID]->parameters[CellParams::TIME_V] += subcycleDt;

      semilagAccTimer.stop();
   }
   //global adjust after each subcycle to keep number of blocks managable. Even the ones not
   //accelerating anyore participate. It is important to keep
   //the spatial dimension to make sure that we do not loose
   //stuff streaming in from other cells, perhaps not connected
   //to the existing distribution function in the cell.
   //- All cells update and communicate their lists of content blocks
   //- Only cells which were accerelated on this step need to be adjusted (blocks removed or added).
   //- Not done here on last step (done after loop)
   if(step < (globalMaxSubcycles - 1))
   {
      adjustVelocityBlocks(mpiGrid, propagatedCells, false, popID);
   }
}

/** Accelerate all particle populations to new time t+dt.
 * This function is AMR safe. NB: acceleration by dt = 0 is not an idempotent operation.
 * @param mpiGrid Parallel grid library.
 * @param dt Time step factor: cells will propagated by dt*CellParams[CellParams::CELLDT] if needed.*/
void calculateAcceleration(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                           Real dt
                          ) {
   typedef Parameters P;
   const vector<CellID>& cells = getLocalCells();   
   set<CellID> cellsToPropagateSet;
   vector<CellID> cellsToPropagateVector;
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
   std::cout << "-----------------calculateAcceleration at t="<<P::t << ", for dtfactor="<<dt<<"\n";
   if (dt == 0.0 && (P::tstep > 0 || P::fractionalTimestep > 0)) {

      // Even if acceleration is turned off we need to adjust velocity blocks
      // because the boundary conditions may have altered the velocity space,
      // and to update changes in no-content blocks during translation.
      for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
// std::cerr << __FILE__<<":"<<__LINE__<< " calling adjustVelocityBlocks at t = " 
//          << P::t << ", preparing to receive; len cells = " << cells.size() <<
//          "\n";
         adjustVelocityBlocks(mpiGrid, cells, true, popID);
      }
   } else {
      // Fairly ugly but no goto
      phiprof::Timer timer {"semilag-acc"};

      
      // Accelerate all particle species
      for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
         int maxSubcycles=0;
         int globalMaxSubcycles;

         // Set active population
         SpatialCell::setCommunicatedSpecies(popID);
         
         // Iterate through all local cells and collect cells to propagate.
         // Ghost cells (spatial cells at the boundary of the simulation 
         // volume) do not need to be propagated:
         vector<CellID> propagatedCells;
         for (size_t c=0; c<cells.size(); ++c) {
            SpatialCell* SC = mpiGrid[cells[c]];
            Real dt_cell;
            if(dt < 0.0) { // Revert to previous real time, stored in cell
                  dt_cell = SC->parameters[CellParams::TIME_R] - P::t;
               }
               else { // dt is a factor of 0.5 or 1.0; so dt_cell is local timestep * dt factor
                  dt_cell = dt*SC->get_tc_dt();
               }
            // if ( (SC->parameters[CellParams::CELLID] == 9 || SC->parameters[CellParams::CELLID] == 11 || SC->parameters[CellParams::CELLID] == 12))
               // std::cout << "vdt on tc  " << SC->parameters[CellParams::TIMECLASS] << " on ftstep " << P::fractionalTimestep << ", dt " << dt_cell <<"\n";

            const vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& vmesh = SC->get_velocity_mesh(popID);
            // disregard boundary cells, in preparation for acceleration
            if (  (SC->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) ||
                  // Include inflow-Maxwellian
                  (P::vlasovAccelerateMaxwellianBoundaries && (SC->sysBoundaryFlag == sysboundarytype::MAXWELLIAN)) ) {
                     if (vmesh.size() != 0){   //do not propagate spatial cells with no blocks
                           if(SC->requested_timeclass_ghosts.size() > 0) {
                              std::cerr << myRank <<": CellID " << cells[c] << " requested timeghosts: ";
                              for (auto i : SC->requested_timeclass_ghosts){
                                 std::cerr << i << " ";
                                 if (SC->get_timeclass_turn_v(i)){
                                    propagatedCells.push_back(cells[c]);
                                    cellsToPropagateSet.insert(cells[c]);
                                 }
                              }
                              std::cerr << "\n";
                           }
                           
                           if ( SC->get_timeclass_turn_v() == true){ // propagate only if it is the cell's turn)
                              propagatedCells.push_back(cells[c]);
                              cellsToPropagateSet.insert(cells[c]);
                           }
                     }
                     //prepare for acceleration, updates max dt for each cell, it
                     //needs to be set to something sensible for _all_ cells, even if
                     //they are not propagated
                     prepareAccelerateCell(SC, popID);

               //update max subcycles for all cells in this process
               maxSubcycles = max((int)getAccelerationSubcycles(SC, dt_cell, popID), maxSubcycles);
               spatial_cell::Population& pop = SC->get_population(popID);
               pop.ACCSUBCYCLES = getAccelerationSubcycles(SC, dt_cell, popID);
            }
            
         } // for loop over cells
         // Compute global maximum for number of subcycles
         MPI_Allreduce(&maxSubcycles, &globalMaxSubcycles, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
         propagatedCells = std::vector<CellID>(cellsToPropagateSet.begin(),cellsToPropagateSet.end());
         // substep global max times
         for(uint step=0; step<(uint)globalMaxSubcycles; ++step) {
            if(step > 0) {
               // prune list of cells to propagate to only contained those which are now subcycled
               vector<CellID> temp;
               for (const auto& cell: cellsToPropagateSet) {
                  if (step < getAccelerationSubcycles(mpiGrid[cell], dt*mpiGrid[cell]->get_tc_dt(), popID) ) {
                     temp.push_back(cell);
                  }
               }
               
               propagatedCells.swap(temp);
            }
            // Accelerate population over one subcycle step
            std::cout << "--------------------- an ACC subcycle ------------------\n";
            calculateAcceleration(popID,(uint)globalMaxSubcycles,step,mpiGrid,propagatedCells,dt);
         } // for-loop over acceleration substeps
         
         // final adjust for all cells, also fixing remote cells.
// std::cerr << __FILE__<<":"<<__LINE__<< " calling adjustVelocityBlocks at t = " 
//          << P::t << ", preparing to receive; len cells = " << cells.size() <<
//          "\n";        
         adjustVelocityBlocks(mpiGrid, cells, true, popID);
      } // for-loop over particle species
      timer.stop();

      //now cellsToPropagateSet contains all cells which have been propagated and whose moments need updating

   } //else

   std::cout << "---------------------------- ACC finished --------------------\n";

   //converting cellsToPropagateSet to vector
   for (auto cell : cellsToPropagateSet) {
      cellsToPropagateVector.push_back(cell);
   }

   // Recalculate "_V" velocity moments
   calculateMoments_V(mpiGrid,cellsToPropagateVector,true);

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
   const int cp_p33
) {
   const vector<CellID>& cells = getLocalCells();

   //Iterate through all local cells
    #pragma omp parallel for
   for (size_t c=0; c<cells.size(); ++c) {
      const CellID cellID = cells[c];
      SpatialCell* SC = mpiGrid[cellID];
      const double tr = SC->parameters[CellParams::TIME_R];
      const double tv = SC->parameters[CellParams::TIME_V];
      SC->parameters[cp_rhom  ] = 0.5* ( SC->parameters[CellParams::RHOM_R] + SC->parameters[CellParams::RHOM_V] );
      SC->parameters[cp_vx]   = 0.5* ( SC->parameters[CellParams::VX_R] + SC->parameters[CellParams::VX_V] );
      SC->parameters[cp_vy] = 0.5* ( SC->parameters[CellParams::VY_R] + SC->parameters[CellParams::VY_V] );
      SC->parameters[cp_vz] = 0.5* ( SC->parameters[CellParams::VZ_R] + SC->parameters[CellParams::VZ_V] );
      SC->parameters[cp_rhoq  ] = 0.5* ( SC->parameters[CellParams::RHOQ_R] + SC->parameters[CellParams::RHOQ_V] );
      SC->parameters[cp_p11]   = 0.5* ( SC->parameters[CellParams::P_11_R] + SC->parameters[CellParams::P_11_V] );
      SC->parameters[cp_p22]   = 0.5* ( SC->parameters[CellParams::P_22_R] + SC->parameters[CellParams::P_22_V] );
      SC->parameters[cp_p33]   = 0.5* ( SC->parameters[CellParams::P_33_R] + SC->parameters[CellParams::P_33_V] );

      for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
         spatial_cell::Population& pop = SC->get_population(popID);
         pop.RHO = 0.5 * ( pop.RHO_R + pop.RHO_V );
         for(int i=0; i<3; i++) {
            pop.V[i] = 0.5 * ( pop.V_R[i] + pop.V_V[i] );
            pop.P[i]    = 0.5 * ( pop.P_R[i] + pop.P_V[i] );
         }
      }
   }
}

double linearInterpolation(double x0, double y0, double x1, double y1, double x) {
    // https://en.wikipedia.org/wiki/Linear_interpolation
    // this is used in the function below.
    return (y0 * (x1 - x) + y1 * (x - x0))/(x1 - x0);
}

double lagrangeInterpolation2order(double x0, double y0, double x1, double y1, double x2, double y2, double x) {
    // https://mathworld.wolfram.com/LagrangeInterpolatingPolynomial.html
    // Lagrange polynomial for interpolation between three points.

    return (y0 * (x - x1) * (x - x2) / ((x0 - x1) * (x0 - x2)) +
            y1 * (x - x0) * (x - x2) / ((x1 - x0) * (x1 - x2)) +
            y2 * (x - x0) * (x - x1) / ((x2 - x0) * (x2 - x1)));
}

double lagrangeInterpolation3order(double x0, double y0, double x1, double y1, double x2, double y2, double x3, double y3, double x) {
    // https://mathworld.wolfram.com/LagrangeInterpolatingPolynomial.html
    // Lagrange polynomial for interpolation between four points.

    return (y0 * (x - x1) * (x - x2) * (x - x3) / ((x0 - x1) * (x0 - x2) * (x0 - x3)) +
            y1 * (x - x0) * (x - x2) * (x - x3) / ((x1 - x0) * (x1 - x2) * (x1 - x3)) +
            y2 * (x - x0) * (x - x1) * (x - x3) / ((x2 - x0) * (x2 - x1) * (x2 - x3)) +
            y3 * (x - x0) * (x - x1) * (x - x2) / ((x3 - x0) * (x3 - x1) * (x3 - x2)));
}

double cubicHermiteSplineInterpolation(double x0, double y0, double x1, double y1, double x2, double y2, double x3, double y3, double x) {
   // https://kluge.in-chemnitz.de/opensource/spline/
   // Cubic Hermite spline interpolation between four points.

   vector<double> xvals = {x0, x1, x2, x3};
   vector<double> yvals = {y0, y1, y2, y3};
   tk::spline s(xvals,yvals,tk::spline::cspline_hermite);

   return s(x);
}


void interpolateMomentsForTimeclasses(
  dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   const int cp_rhom,
   const int cp_vx,
   const int cp_vy,
   const int cp_vz,
   const int cp_rhoq,
   const int cp_p11,
   const int cp_p22,
   const int cp_p33,
   const int fracTimeStep, // goes from 0 to 2^maxtimeclass-1
   const int maxTC,
   const bool dt2 // true if second moment / dt2
) {

   const vector<CellID>& cells = getLocalCells();

   #pragma omp parallel for
   for (size_t c=0; c<cells.size(); ++c) {
      const CellID cellID = cells[c];
      SpatialCell* SC = mpiGrid[cellID];
      const int timeclass = SC->parameters[CellParams::TIMECLASS];
      // const double tr = SC->parameters[CellParams::TIME_R];
      // const double tv = SC->parameters[CellParams::TIME_V];

      if (timeclass == maxTC) {
         // calculateInterpolatedVelocityMoments functionality here, if timeclass is the max one.
         if (!dt2) {
            SC->parameters[cp_rhom] = 0.5* ( SC->parameters[CellParams::RHOM_R_PREV] + SC->parameters[CellParams::RHOM_V] );
            SC->parameters[cp_vx]   = 0.5* ( SC->parameters[CellParams::VX_R_PREV] + SC->parameters[CellParams::VX_V] );
            SC->parameters[cp_vy] = 0.5* ( SC->parameters[CellParams::VY_R_PREV] + SC->parameters[CellParams::VY_V] );
            SC->parameters[cp_vz] = 0.5* ( SC->parameters[CellParams::VZ_R_PREV] + SC->parameters[CellParams::VZ_V] );
            SC->parameters[cp_rhoq] = 0.5* ( SC->parameters[CellParams::RHOQ_R_PREV] + SC->parameters[CellParams::RHOQ_V] );
            SC->parameters[cp_p11]   = 0.5* ( SC->parameters[CellParams::P_11_R_PREV] + SC->parameters[CellParams::P_11_V] );
            SC->parameters[cp_p22]   = 0.5* ( SC->parameters[CellParams::P_22_R_PREV] + SC->parameters[CellParams::P_22_V] );
            SC->parameters[cp_p33]   = 0.5* ( SC->parameters[CellParams::P_33_R_PREV] + SC->parameters[CellParams::P_33_V] );

            // for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
            //    spatial_cell::Population& pop = SC->get_population(popID);
            //    pop.RHO = 0.5 * ( pop.RHO_R_PREV + pop.RHO_V );
            //    for(int i=0; i<3; i++) {
            //       pop.V[i] = 0.5 * ( pop.V_R_PREV[i] + pop.V_V[i] );
            //       pop.P[i] = 0.5 * ( pop.P_R_PREV[i] + pop.P_V[i] );
            //    }
            // }

         } else {
            SC->parameters[cp_rhom  ] = 0.5* ( SC->parameters[CellParams::RHOM_R] + SC->parameters[CellParams::RHOM_V] );
            SC->parameters[cp_vx]   = 0.5* ( SC->parameters[CellParams::VX_R] + SC->parameters[CellParams::VX_V] );
            SC->parameters[cp_vy] = 0.5* ( SC->parameters[CellParams::VY_R] + SC->parameters[CellParams::VY_V] );
            SC->parameters[cp_vz] = 0.5* ( SC->parameters[CellParams::VZ_R] + SC->parameters[CellParams::VZ_V] );
            SC->parameters[cp_rhoq  ] = 0.5* ( SC->parameters[CellParams::RHOQ_R] + SC->parameters[CellParams::RHOQ_V] );
            SC->parameters[cp_p11]   = 0.5* ( SC->parameters[CellParams::P_11_R] + SC->parameters[CellParams::P_11_V] );
            SC->parameters[cp_p22]   = 0.5* ( SC->parameters[CellParams::P_22_R] + SC->parameters[CellParams::P_22_V] );
            SC->parameters[cp_p33]   = 0.5* ( SC->parameters[CellParams::P_33_R] + SC->parameters[CellParams::P_33_V] );

         //    for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
         //       spatial_cell::Population& pop = SC->get_population(popID);
         //       pop.RHO = 0.5 * ( pop.RHO_R + pop.RHO_V );
         //       for(int i=0; i<3; i++) {
         //          pop.V[i] = 0.5 * ( pop.V_R[i] + pop.V_V[i] );
         //          pop.P[i] = 0.5 * ( pop.P_R[i] + pop.P_V[i] );
         //       }
         //    }
            }

     } else { // this block if timeclass != maxTC

         int reverseTC = maxTC - timeclass;
         double RTCpow = pow(2, reverseTC);
         double modul = fracTimeStep % (int)RTCpow;
         double normModul = modul/RTCpow;
         if (!dt2) {
            normModul += 0.25/RTCpow;
         } else {
            normModul += 0.25/RTCpow+0.5/RTCpow;
         }

         // !! if translation and acceleration are changed to not update both on fractimestep 0, this will break

         if (P::tcVMomentPropagation == true) {
            if (SC->get_timeclass_turn_v()) { // clamping down values as they get updated
               SC->parameters[cp_vx] = lagrangeInterpolation2order(-0.25, 0.5*(SC->parameters[CellParams::VX_V_PREV]+SC->parameters[CellParams::VX_R_PREV]), 0.25, 0.5*(SC->parameters[CellParams::VX_V]+SC->parameters[CellParams::VX_R_PREV]), 0.75, 0.5*(SC->parameters[CellParams::VX_V]+SC->parameters[CellParams::VX_R]), normModul);
               SC->parameters[cp_vy] = lagrangeInterpolation2order(-0.25, 0.5*(SC->parameters[CellParams::VY_V_PREV]+SC->parameters[CellParams::VY_R_PREV]), 0.25, 0.5*(SC->parameters[CellParams::VY_V]+SC->parameters[CellParams::VY_R_PREV]), 0.75, 0.5*(SC->parameters[CellParams::VY_V]+SC->parameters[CellParams::VY_R]), normModul);
               SC->parameters[cp_vz] = lagrangeInterpolation2order(-0.25, 0.5*(SC->parameters[CellParams::VZ_V_PREV]+SC->parameters[CellParams::VZ_R_PREV]), 0.25, 0.5*(SC->parameters[CellParams::VZ_V]+SC->parameters[CellParams::VZ_R_PREV]), 0.75, 0.5*(SC->parameters[CellParams::VZ_V]+SC->parameters[CellParams::VZ_R]), normModul);
            } else {

            Eigen::Transform<Real,3,Eigen::Affine> vUpdateMatrix = compute_acceleration_transformation(SC, 0, SC->parameters[CellParams::TIMECLASSDT]*(1.0/RTCpow));
            const Eigen::Matrix<Real,3,1> V(SC->parameters[cp_vx], SC->parameters[cp_vy], SC->parameters[cp_vz]);
            Eigen::Matrix<Real,3,1> V_updated = vUpdateMatrix * V;

            SC->parameters[cp_vx] = V_updated(0);
            SC->parameters[cp_vy] = V_updated(1);
            SC->parameters[cp_vz] = V_updated(2);
            }
         }

         // temporary arrays for true moments.
         int nMomentsToInterp = (P::tcVMomentPropagation) ? 5 : 8;

         double avgMoments1[nMomentsToInterp];
         double avgMoments2[nMomentsToInterp];
         double avgMoments3[nMomentsToInterp];
         double avgMoments4[nMomentsToInterp];
         double avgMoments5[nMomentsToInterp];

         // for type of interpolation (P::tcMomentInterpolationType) -1 is cubic C^1 Hermite spline, 1 is linear, 2 is lagrange 2nd order, 3 is lagrange 3rd order.

         if (SC->get_timeclass_turn_v()) { // aka if translation moments are ahead of acceleration moments
            // here, temporal order from newest to oldest is _R, _V, _R_PREV, _V_PREV, _R_PREV_PREV, _V_PREV_PREV
            for (int i=0; i<nMomentsToInterp; i++) {
               avgMoments1[i] = 0.5*( SC->parameters[CellParams::RHOM_R+i] + SC->parameters[CellParams::RHOM_V+i] ); // at 0.75
               avgMoments2[i] = 0.5*( SC->parameters[CellParams::RHOM_V+i] + SC->parameters[CellParams::RHOM_R_PREV+i] ); // at 0.25
               avgMoments3[i] = 0.5*( SC->parameters[CellParams::RHOM_R_PREV+i] + SC->parameters[CellParams::RHOM_V_PREV+i] );  // at -0.25
               avgMoments4[i] = 0.5*( SC->parameters[CellParams::RHOM_V_PREV+i] + SC->parameters[CellParams::RHOM_R_PREV_PREV+i] ); // at -0.75
               avgMoments5[i] = 0.5*( SC->parameters[CellParams::RHOM_R_PREV_PREV+i] + SC->parameters[CellParams::RHOM_V_PREV_PREV+i] ); // at -1.25
            }

            for (int i=0; i<nMomentsToInterp; i++) {

               switch(P::tcMomentInterpolationType) {
                  case 1:
                     if (normModul < 0.25) {
                        SC->parameters[cp_rhom+i] = linearInterpolation(-0.25, avgMoments3[i], 0.25, avgMoments2[i], normModul);
                     } else if (normModul < 0.75) {
                        SC->parameters[cp_rhom+i] = linearInterpolation(0.25, avgMoments2[i], 0.75, avgMoments1[i], normModul);
                     } else {
                     // this is never reached
                     }
                     break;
                  case 2:
                     SC->parameters[cp_rhom+i] = lagrangeInterpolation2order(-0.25, avgMoments3[i], 0.25, avgMoments2[i], 0.75, avgMoments1[i], normModul);
                     break;
                  case 3:
                     SC->parameters[cp_rhom+i] = lagrangeInterpolation3order(-0.75, avgMoments4[i], -0.25, avgMoments3[i], 0.25, avgMoments2[i], 0.75, avgMoments1[i], normModul);
                     break;
                  case -1:
                     SC->parameters[cp_rhom+i] = cubicHermiteSplineInterpolation(-0.75, avgMoments4[i], -0.25, avgMoments3[i], 0.25, avgMoments2[i], 0.75, avgMoments1[i], normModul);
                     break;
               }
            }
         } else { // if acceleration moments are ahead of translation moments
            // here, temporal order from newest to oldest is _V, _R, _V_PREV, _R_PREV, _V_PREV_PREV, _R_PREV_PREV
            for (int i=0; i<nMomentsToInterp; i++) {
               avgMoments1[i] = 0.5*( SC->parameters[CellParams::RHOM_V+i] + SC->parameters[CellParams::RHOM_R+i] ); // at 1.25
               avgMoments2[i] = 0.5*( SC->parameters[CellParams::RHOM_R+i] + SC->parameters[CellParams::RHOM_V_PREV+i] ); // at 0.75
               avgMoments3[i] = 0.5*( SC->parameters[CellParams::RHOM_V_PREV+i] + SC->parameters[CellParams::RHOM_R_PREV+i] ); // at 0.25
               avgMoments4[i] = 0.5*( SC->parameters[CellParams::RHOM_R_PREV+i] + SC->parameters[CellParams::RHOM_V_PREV_PREV+i] ); // at -0.25
               avgMoments5[i] = 0.5*( SC->parameters[CellParams::RHOM_V_PREV_PREV+i] + SC->parameters[CellParams::RHOM_R_PREV_PREV+i] ); // at -0.75
            }

            for (int i=0; i<nMomentsToInterp; i++) {
               switch(P::tcMomentInterpolationType) {
                  case 1:
                     if (normModul < 0.25) {
                        SC->parameters[cp_rhom+i] = linearInterpolation(-0.25, avgMoments4[i], 0.25, avgMoments3[i], normModul);
                     } else if (normModul < 0.75) {
                        SC->parameters[cp_rhom+i] = linearInterpolation(0.25, avgMoments3[i], 0.75, avgMoments2[i], normModul);
                     } else {
                        SC->parameters[cp_rhom+i] = linearInterpolation(0.75, avgMoments2[i], 1.25, avgMoments1[i], normModul);
                     }
                     break;
                  case 2:
                     SC->parameters[cp_rhom+i] = lagrangeInterpolation2order(0.25, avgMoments3[i], 0.75, avgMoments2[i], 1.25, avgMoments1[i], normModul);
                     break;
                  case 3:
                     SC->parameters[cp_rhom+i] = lagrangeInterpolation3order(-0.25, avgMoments4[i], 0.25, avgMoments3[i], 0.75, avgMoments2[i], 1.25, avgMoments1[i], normModul);
                     break;
                  case -1:
                     SC->parameters[cp_rhom+i] = cubicHermiteSplineInterpolation(-0.25, avgMoments4[i], 0.25, avgMoments3[i], 0.75, avgMoments2[i], 1.25, avgMoments1[i], normModul);
                     break;
               }
            }
         }
      } // if timeclass = maxTC
   } // for loop over cells
}

void updateParticlePopulations(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid) {

   const vector<CellID>& cells = getLocalCells();
   #pragma omp parallel for
   for (size_t c=0; c<cells.size(); ++c) {

      const CellID cellID = cells[c];
      SpatialCell* SC = mpiGrid[cellID];

      if (SC->get_timeclass_turn_v() == true) {
         for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
            spatial_cell::Population& pop = SC->get_population(popID);
            pop.RHO = 0.5 * ( pop.RHO_R + pop.RHO_V );
            for(int i=0; i<3; i++) {
               pop.V[i] = 0.5 * ( pop.V_R[i] + pop.V_V[i] );
               pop.P[i] = 0.5 * ( pop.P_R[i] + pop.P_V[i] );
            }
         }
      }
   }
}


void calculateInitialVelocityMoments(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid) {
   const vector<CellID>& cells = getLocalCells();
   phiprof::Timer timer {"Calculate moments"};

   // Iterate through all local cells (incl. system boundary cells):
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
   } // for-loop over spatial cells
}

void updatePreviousVMoments(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid, bool isInitialization) {

   const vector<CellID>& cells = getLocalCells();
   for (size_t c=0; c<cells.size(); c++) {
      const CellID cellID = cells[c];
      SpatialCell* SC = mpiGrid[cellID];
      double tdiff = SC->parameters[CellParams::TIMECLASSDT];

      if (isInitialization == true) {
         // initializiing all _PREV and _PREV_PREV moments to same values as _V moments
         // this block only gets executed once, when the simulation is started
         SC->parameters[CellParams::RHOM_V_PREV_PREV] = SC->parameters[CellParams::RHOM_V];
         SC->parameters[CellParams::RHOQ_V_PREV_PREV] = SC->parameters[CellParams::RHOQ_V];
         SC->parameters[CellParams::P_11_V_PREV_PREV] = SC->parameters[CellParams::P_11_V];
         SC->parameters[CellParams::P_22_V_PREV_PREV] = SC->parameters[CellParams::P_22_V];
         SC->parameters[CellParams::P_33_V_PREV_PREV] = SC->parameters[CellParams::P_33_V];

         SC->parameters[CellParams::RHOM_V_PREV] = SC->parameters[CellParams::RHOM_V];
         SC->parameters[CellParams::RHOQ_V_PREV] = SC->parameters[CellParams::RHOQ_V];
         SC->parameters[CellParams::P_11_V_PREV] = SC->parameters[CellParams::P_11_V];
         SC->parameters[CellParams::P_22_V_PREV] = SC->parameters[CellParams::P_22_V];
         SC->parameters[CellParams::P_33_V_PREV] = SC->parameters[CellParams::P_33_V];

         SC->parameters[CellParams::RHOM_R_PREV_PREV] = SC->parameters[CellParams::RHOM_R];
         SC->parameters[CellParams::RHOQ_R_PREV_PREV] = SC->parameters[CellParams::RHOQ_R];
         SC->parameters[CellParams::P_11_R_PREV_PREV] = SC->parameters[CellParams::P_11_R];
         SC->parameters[CellParams::P_22_R_PREV_PREV] = SC->parameters[CellParams::P_22_R];
         SC->parameters[CellParams::P_33_R_PREV_PREV] = SC->parameters[CellParams::P_33_R];

         SC->parameters[CellParams::RHOM_R_PREV] = SC->parameters[CellParams::RHOM_R];
         SC->parameters[CellParams::RHOQ_R_PREV] = SC->parameters[CellParams::RHOQ_R];
         SC->parameters[CellParams::P_11_R_PREV] = SC->parameters[CellParams::P_11_R];
         SC->parameters[CellParams::P_22_R_PREV] = SC->parameters[CellParams::P_22_R];
         SC->parameters[CellParams::P_33_R_PREV] = SC->parameters[CellParams::P_33_R];

         if (true) {
                  
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

            // updating _PREV moments
            SC->parameters[CellParams::RHOM_V_PREV] = SC->parameters[CellParams::RHOM_V];
            SC->parameters[CellParams::VX_V_PREV] = SC->parameters[CellParams::VX_V];
            SC->parameters[CellParams::VY_V_PREV] = SC->parameters[CellParams::VY_V];
            SC->parameters[CellParams::VZ_V_PREV] = SC->parameters[CellParams::VZ_V];
            SC->parameters[CellParams::RHOQ_V_PREV] = SC->parameters[CellParams::RHOQ_V];
            SC->parameters[CellParams::P_11_V_PREV] = SC->parameters[CellParams::P_11_V];
            SC->parameters[CellParams::P_22_V_PREV] = SC->parameters[CellParams::P_22_V];
            SC->parameters[CellParams::P_33_V_PREV] = SC->parameters[CellParams::P_33_V];   
         }
      }
   }  
}