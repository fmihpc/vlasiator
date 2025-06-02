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

// #include <dccrg.hpp>
// #include <dccrg_cartesian_geometry.hpp>
#include <phiprof.hpp>
#include "../definitions.h"

#include "cpu_acc_semilag.hpp"
#include "cpu_acc_intersections.hpp"
#include "cpu_acc_map.hpp"
#include "cpu_acc_transform.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

/*!
  Calls semi-lagrangian acceleration routines for the provided list of cells

 * @param mpiGrid DCCRG container of spatial cells
 * @param acceleratedCells vector of cells for which to perform acceleration
 * @param popID ID of the accelerated particle species.
 * @param map_order Order in which vx,vy,vz mappings are performed.
*/

void cpu_accelerate_cells(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                          const std::vector<CellID>& acceleratedCells,
                          const uint popID,
                          const uint map_order,
                          const int tc
   ) {
   int timerId {phiprof::initializeTimer("cell-semilag-acc")};
   int intersections_id {phiprof::initializeTimer("cell-compute-intersections")};

   #pragma omp parallel // Launch workshare region
   {
      // Calculate intersections (should be constant cost per cell)
      #pragma omp for schedule(static,1)
      for (size_t c=0; c<acceleratedCells.size(); ++c) {
         const CellID cellID = acceleratedCells[c];
         SpatialCell* SC = mpiGrid[cellID];
         Population& pop = SC->get_population(popID);
         compute_cell_intersections(SC, popID, map_order, pop.subcycleDt, intersections_id);
      }
      #pragma omp barrier
      // Semi-Lagrangian acceleration for all cells active in this subcycle,
      // dimension-by-dimension. Dynamic cost due to varying block counts.
      #pragma omp for schedule(dynamic,1)
      for (size_t c=0; c<acceleratedCells.size(); ++c) {
         const CellID cellID = acceleratedCells[c];
         SpatialCell* SC = mpiGrid[cellID];
         phiprof::Timer semilagAccTimer {timerId};
         cpu_accelerate_cell(SC,popID,map_order,SC->get_tc_dt(),tc);
         semilagAccTimer.stop();
      }
   }
}


/*!
  Compute the number of subcycles needed for the acceleration of the particle
  species in the spatial cell. Note that one should first prepare to
  accelerate the cell with prepareAccelerateCell.

 * @param spatial_cell Spatial cell containing the accelerated population.
 * @param popID ID of the accelerated particle species.
*/

uint getAccelerationSubcycles(SpatialCell* spatial_cell, Real dt, const uint popID)
{
   //return max( convert<uint>(ceil(dt*spatial_cell->CellParams[CELLPARAMS::TIMECLASSDT] / spatial_cell->get_max_v_dt(popID))), 1u);
   return max( convert<uint>(ceil(dt / spatial_cell->get_max_v_dt(popID))), 1u);
}
/*!
  Compute the number of subcycles needed from maxVdt and target dt.

 * @param spatial_cell Spatial cell containing the accelerated population.
 * @param popID ID of the accelerated particle species.
*/

uint getAccelerationSubcycles(Real maxVdt, Real dt)
{
   //return max( convert<uint>(ceil(dt*spatial_cell->CellParams[CELLPARAMS::TIMECLASSDT] / spatial_cell->get_max_v_dt(popID))), 1u);
   return max( convert<uint>(ceil(dt / maxVdt)), 1u);
}

/*!
  Propagates the distribution function in velocity space of given real
  space cell.

  Based on SLICE-3D algorithm: Zerroukat, M., and T. Allen. "A
  three‐dimensional monotone and conservative semi‐Lagrangian scheme
  (SLICE‐3D) for transport problems." Quarterly Journal of the Royal
  Meteorological Society 138.667 (2012): 1640-1651.

 * @param spatial_cell Spatial cell containing the accelerated population.
 * @param popID ID of the accelerated particle species.
 * @param map_order Order in which vx,vy,vz mappings are performed.
*/

void cpu_accelerate_cell(SpatialCell* spatial_cell,
                         const uint popID,     
                         const uint map_order,
                         const Real& dt,
                         int timeclass) {
   //double t1 = MPI_Wtime();
   vmesh::VelocityMesh* vmeshPtr = NULL; 
   vmesh::VelocityBlockContainer* blockContainerPtr = NULL;
   // Main branch: check if time-ghost data required, handle that and recurse as needed
   // .... of course this is now substepping unnecessarily...
   // std::cout << spatial_cell->parameters[CellParams::CELLID] << "c Accelerate, initial refs" << " vmesh " << &vmesh << " blockContainer " << &blockContainer << "\n";

   int tc_delta = -999;

   if(false && tc_delta == 0) // Handles the default case (called from vlasovmover)
   {
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
            
            if (!P::tc_leapfrog_init) {
               spatial_cell->set_velocity_mesh_ghost(popID, i);
               spatial_cell->set_velocity_blocks_ghost(popID, i); 
               std::cerr << "c"<< spatial_cell->parameters[CellParams::CELLID]<<" with " << spatial_cell->get_velocity_blocks(popID)->size() <<" blocks\n";

               cpu_accelerate_cell(spatial_cell, popID, map_order, dt/pow(2,tc_d), tc_d);
               if (spatial_cell->parameters[CellParams::CELLID]  == 16){
                  std::cout << "16c Initial nudge" << "\n";
               }
               double sum = 0;
               for (auto bd: spatial_cell->get_velocity_blocks(popID)->getDataVector_raw()){
                  sum+=bd;
               }
               std::cerr << "c"<< spatial_cell->parameters[CellParams::CELLID]<< " post-acc sum " << sum << " with " << spatial_cell->get_velocity_blocks(popID)->size() <<" blocks\n";
            }
            else if (spatial_cell->get_timeclass_turn_v()) {
               spatial_cell->set_velocity_mesh_ghost(popID, i);
               spatial_cell->set_velocity_blocks_ghost(popID, i); 
               std::cerr << "c"<< spatial_cell->parameters[CellParams::CELLID]<<" with " << spatial_cell->get_velocity_blocks(popID)->size() <<" blocks\n";

               cpu_accelerate_cell(spatial_cell, popID, map_order, dt/pow(2,tc_d)*1/2, tc_d);
               if (spatial_cell->parameters[CellParams::CELLID]  == 16){
                  std::cout << "16c tc-0 copy and nudge" << "\n";
               }
               double sum = 0;
               for (auto bd: spatial_cell->get_velocity_blocks(popID)->getDataVector_raw()){
                  sum+=bd;
               }
               std::cerr << "c"<< spatial_cell->parameters[CellParams::CELLID]<< " post-acc sum " << sum << " with " << spatial_cell->get_velocity_blocks(popID)->size() <<" blocks\n";
            }
            else if (spatial_cell->get_timeclass_turn_v(i)) {
               cpu_accelerate_cell(spatial_cell, popID, map_order, dt/pow(2,tc_d), tc_d);
               if (spatial_cell->parameters[CellParams::CELLID]  == 16){
                  std::cout << "16c tc-1 nudge" << "\n";
               }
            }
            else if (spatial_cell->requested_timeclass_copy_ghosts.count(i) && spatial_cell->get_timeclass_turn_v(i)){
               spatial_cell->set_velocity_mesh_ghost(popID, i);
               spatial_cell->set_velocity_blocks_ghost(popID, i); 
               std::cerr << "c"<< spatial_cell->parameters[CellParams::CELLID]<<" with " << spatial_cell->get_velocity_blocks(popID)->size() <<" blocks\n";
               cpu_accelerate_cell(spatial_cell, popID, map_order, dt/pow(2,tc_d)*3/2, tc_d);
               if (spatial_cell->parameters[CellParams::CELLID]  == 15){
                  std::cout << "15c tc-1 copy and nudge again" << "\n";
               }
               double sum = 0;

               for (auto bd: spatial_cell->get_velocity_blocks(popID)->getDataVector_raw()){
                  sum+=bd;
               }
               std::cerr << "c " << spatial_cell->parameters[CellParams::CELLID]<< " post-acc sum " << sum << " with " << spatial_cell->get_velocity_blocks(popID)->size() <<" blocks\n";
            }
            else{
               //do nothing
            }
            
         }
         else if (tc_d < 0){

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
            if (!P::tc_leapfrog_init){
                              
               // if it is slower-tc's turn, we are synced at after-trans state
               // -> Copy state, but needs to acc by half-tc-0-dt ("always init")
               std::cerr << "c"<< spatial_cell->parameters[CellParams::CELLID]<<"\n";
               spatial_cell->set_velocity_mesh_ghost(popID, i);
               spatial_cell->set_velocity_blocks_ghost(popID, i);
               // cpu_accelerate_cell(spatial_cell, popID, map_order, P::timeclassDt[i], tc_d);
               double dtt = dt/pow(2,tc_d);
               if (spatial_cell->parameters[CellParams::CELLID]  == 17){
                  std::cout << "17c Copy-and-init-nudge from "<< dt << " div by " << pow(2,tc_d) << "\n";
                  std::cout << __FILE__<<":"<<__LINE__<< "\t17c Copying and propagating ghost at tc " << spatial_cell->get_tc() + tc_delta << " by dt = " << dtt << " being run at cell " << "\n";
               }
               else{
                  std::cout << __FILE__<<":"<<__LINE__<< "\t Copying and propagating ghost at tc " << spatial_cell->get_tc() + tc_delta << " by dt = " << dtt << " being run at cell " << "\n";
               }
               
               cpu_accelerate_cell(spatial_cell, popID, map_order, dtt, tc_d);
            
            }
            else if (spatial_cell->get_timeclass_turn_v(i) || (P::tstep == 0 && P::fractionalTimestep == 0)) {
               
               // if it is slower-tc's turn, we are synced at after-trans state
               // -> Copy state, but needs to acc by half-tc-0-dt ("always init")
               std::cerr << "c"<< spatial_cell->parameters[CellParams::CELLID]<<"\n";
               spatial_cell->set_velocity_mesh_ghost(popID, i);
               spatial_cell->set_velocity_blocks_ghost(popID, i);
               // cpu_accelerate_cell(spatial_cell, popID, map_order, P::timeclassDt[i], tc_d);
               double dtt = 5./4.*dt/pow(2,tc_d);
               if (spatial_cell->parameters[CellParams::CELLID]  == 17){
                  std::cout << "17c Copy-and-init-nudge from "<< dt << " div by " << pow(2,tc_d) << "\n";
                  std::cout << __FILE__<<":"<<__LINE__<< "\t17c Copying and propagating ghost at tc " << spatial_cell->get_tc() + tc_delta << " by dt = " << dtt << " being run at cell " << "\n";
               }
               else{
                  std::cout << __FILE__<<":"<<__LINE__<< "\t Copying and propagating ghost at tc " << spatial_cell->get_tc() + tc_delta << " by dt = " << dtt << " being run at cell " << "\n";
               }
               

               cpu_accelerate_cell(spatial_cell, popID, map_order, dtt, tc_d);
            }
            else{
               if (spatial_cell->parameters[CellParams::CELLID]  == 17){
                  std::cout << "17c Do-nothing" << "\n";
               }
               // ghost can just stay put!
               std::cout << __FILE__<<":"<<__LINE__<< "\tLeaving ghost as is at tc " << spatial_cell->get_tc() + tc_delta << " by dt = " << dt << " being run at cell " << "\n";
            }

            // cpu_accelerate_cell(spatial_cell, popID, map_order, dt/pow(2,tc_d), tc_d);

         }
      } // for over ghost requests
      if (!spatial_cell->get_timeclass_turn_v())
      { // Not our turn and ghosts have been handled -> skidaddle
         return;
      }

      vmeshPtr = spatial_cell->get_velocity_mesh(popID);
      blockContainerPtr = spatial_cell->get_velocity_blocks(popID);


      std::cout << __FILE__<<":"<<__LINE__ << "\t"<<spatial_cell->parameters[CellParams::CELLID] << "c Accelerate tc " << spatial_cell->get_tc() << " with tcdelta " 
      << tc_delta << " at dt = " << dt << "\n";//"; vmesh " << &vmesh << " blockContainer " << &blockContainer << "\n";

   }
   // Ghost branch: select the ghost vmesh defined by tc_delta instead for acc
   // dt we already have adjusted as needed
   else{

      // std::cout << __FILE__<<":"<<__LINE__<<"\t"<< spatial_cell->parameters[CellParams::CELLID] << "c\tRequest-accelerate at tc " << spatial_cell->get_tc() + tc_delta << " at dt = " << dt << "\n";
      vmeshPtr = spatial_cell->get_velocity_mesh(popID,timeclass);
      blockContainerPtr = spatial_cell->get_velocity_blocks(popID,timeclass);
    }
    
   vmesh::VelocityMesh& vmesh = *vmeshPtr;
   vmesh::VelocityBlockContainer& blockContainer = *blockContainerPtr;


   // vmesh::VelocityBlockContainer<vmesh::LocalID>& blockContainer = spatial_cell->get_velocity_blocks(popID);

   // compute transform, forward in time and backward in time
   //compute the transform performed in this acceleration
   Eigen::Transform<Real,3,Eigen::Affine> fwd_transform= compute_acceleration_transformation(spatial_cell,popID,dt);
   Eigen::Transform<Real,3,Eigen::Affine> bwd_transform= fwd_transform.inverse();

   // if (spatial_cell->parameters[CellParams::CELLID] == 16)
   // #pragma omp critical(output)
   // {
   //    std::cout<< "Cellid " << spatial_cell->parameters[CellParams::CELLID] << 
   //    " at t=" << spatial_cell->parameters[CellParams::TIME_V] <<": dt " << 
   //    dt <<  "\n";
   //    std::stringstream ss;
   //    ss.precision(64);
   //    ss << std::scientific;
   //    ss << fwd_transform(0,0) << " " << fwd_transform(0,1) << " " << fwd_transform(0,2) << "\n" <<
   //          fwd_transform(1,0) << " " << fwd_transform(1,1) << " " << fwd_transform(1,2) << "\n" <<
   //          fwd_transform(2,0) << " " << fwd_transform(2,1) << " " << fwd_transform(2,2) << "\n";
   //    std::cout << ss.str() << "\n";
   // }


   const uint8_t refLevel = 0;
   Real intersection_z,intersection_z_di,intersection_z_dj,intersection_z_dk;
   Real intersection_x,intersection_x_di,intersection_x_dj,intersection_x_dk;
   Real intersection_y,intersection_y_di,intersection_y_dj,intersection_y_dk;
   int intersections_id {phiprof::initializeTimer("compute-intersections")};
   int mapping_id {phiprof::initializeTimer("compute-mapping")};
//    switch(map_order){
//       case 0: {
//          //Map order XYZ
//          phiprof::Timer intersectionsTimer {"compute-intersections"};
//          compute_intersections_1st(&vmesh,bwd_transform, fwd_transform, 0,
//                                    intersection_x,intersection_x_di,intersection_x_dj,intersection_x_dk);
//          compute_intersections_2nd(&vmesh,bwd_transform, fwd_transform, 1,
//                                    intersection_y,intersection_y_di,intersection_y_dj,intersection_y_dk);
//          compute_intersections_3rd(&vmesh,bwd_transform, fwd_transform, 2,
//                                    intersection_z,intersection_z_di,intersection_z_dj,intersection_z_dk);
//          intersectionsTimer.stop();
//          phiprof::Timer mappingTimer {mapping_id};
//          map_1d(spatial_cell, popID, intersection_x,intersection_x_di,intersection_x_dj,intersection_x_dk,0); // map along x
//          map_1d(spatial_cell, popID, intersection_y,intersection_y_di,intersection_y_dj,intersection_y_dk,1); // map along y
//          map_1d(spatial_cell, popID, intersection_z,intersection_z_di,intersection_z_dj,intersection_z_dk,2); // map along z
//          mappingTimer.stop();
// =======

   Population& pop = spatial_cell->get_population(popID);
   switch(map_order){
      case 0: {
         //Map order XYZ
         map_1d(spatial_cell, popID, pop.intersection_x,
                pop.intersection_x_di,pop.intersection_x_dj,pop.intersection_x_dk,0); // map along x
         map_1d(spatial_cell, popID, pop.intersection_y,
                pop.intersection_y_di,pop.intersection_y_dj,pop.intersection_y_dk,1); // map along y
         map_1d(spatial_cell, popID, pop.intersection_z,
                pop.intersection_z_di,pop.intersection_z_dj,pop.intersection_z_dk,2); // map along z
         break;
      }
      case 1: {
         //Map order YZX
         map_1d(spatial_cell, popID, pop.intersection_y,
                pop.intersection_y_di,pop.intersection_y_dj,pop.intersection_y_dk,1); // map along y
         map_1d(spatial_cell, popID, pop.intersection_z,
                pop.intersection_z_di,pop.intersection_z_dj,pop.intersection_z_dk,2); // map along z
         map_1d(spatial_cell, popID, pop.intersection_x,
                pop.intersection_x_di,pop.intersection_x_dj,pop.intersection_x_dk,0); // map along x
         break;
      }
      case 2: {
         //Map order Z X Y
         map_1d(spatial_cell, popID, pop.intersection_z,
                pop.intersection_z_di,pop.intersection_z_dj,pop.intersection_z_dk,2); // map along z
         map_1d(spatial_cell, popID, pop.intersection_x,
                pop.intersection_x_di,pop.intersection_x_dj,pop.intersection_x_dk,0); // map along x
         map_1d(spatial_cell, popID, pop.intersection_y,
                pop.intersection_y_di,pop.intersection_y_dj,pop.intersection_y_dk,1); // map along y
         break;
      }
   }
}
