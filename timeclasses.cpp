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

#include "timeclasses.hpp"
#include "grid.h"
#include "object_wrapper.h"

bool isDtTooLarge(Real dt, Real rdt, Real vdt, Real fsdt){
   return (dt > rdt * P::vlasovSolverMaxCFL ||
           dt > vdt * P::vlasovSolverMaxCFL * P::maxSlAccelerationSubcycles ||
           dt > fsdt * P::fieldSolverMaxCFL * P::maxFieldSolverSubcycles);
}

bool isDtTooSmall(Real dt, Real rdt, Real vdt, Real fsdt){
   return (dt < rdt * P::vlasovSolverMinCFL &&
           dt < vdt * P::vlasovSolverMinCFL * P::maxSlAccelerationSubcycles &&
           dt < fsdt * P::fieldSolverMinCFL * P::maxFieldSolverSubcycles);
}


// returns empty vector if all timeclasses are fine (= their timestep fits their timeclass)
// if not, returns those cells which need a bigger timeclass
// recalculates all cellwise tc limits

// skips over cells that are certain boundaries as defined above, as those should not affect timestep length
std::vector<CellID> checkCellTimeclasses(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid) {

   std::vector<CellID> retVec = {};
   const vector<CellID>& cells = getLocalCells();

   for (vector<CellID>::const_iterator cell_id=cells.begin(); cell_id!=cells.end(); ++cell_id) {
   
      if (mpiGrid[*cell_id]->cellIsTimeclassRelevant()) {
         if (!(mpiGrid[*cell_id]->cellTimeclassIsCorrect())) {
            retVec.push_back(*cell_id);
         }
      }
   }

   return retVec;
}

void updateTimeclassDts(Real fsdt, const bool applyModifier) {

   // reduce fsdt by buffer amount

   //fsdt /= pow(2.0, P::timeclassBuffer);

   std::vector<Real> newTimeclassDts(P::currentMaxTimeclass+1);
   //logFile << std::endl;
   //logFile << "(TC) timeclassDts set to " << std::endl;
   for(int i = 0; i <= P::currentMaxTimeclass; ++i){
      if (applyModifier) {
         newTimeclassDts[i] = fsdt*pow(2,P::currentMaxTimeclass - i)*P::dtSettingModifier;
      } else {
         newTimeclassDts[i] = fsdt*pow(2,P::currentMaxTimeclass - i);
      }
      //logFile << newTimeclassDts[i] << "s, ";
   }
   //logFile << std::endl;
   //logFile << std::endl;
   P::timeclassDt = newTimeclassDts;

}

void increaseTimeclass(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                              const std::vector<CellID>& cellsToIncreaseTimeclass,
                              bool& additionalTimeclassCreated) {
   phiprof::Timer increaseTimeclassTimer {"increase-timeclass"};

   additionalTimeclassCreated = false;

   // Increase timeclass for given cells

   if (P::fractionalTimestep == 0) {
      // first we step them back 
      //calculateAcceleration(mpiGrid, -0.5, true, cellsToIncreaseTimeclass);


      for (size_t c=0; c<cellsToIncreaseTimeclass.size(); ++c) {
         const CellID cell = cellsToIncreaseTimeclass[c];
         SpatialCell* spatialCell = mpiGrid[cell];
         for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {



            // before we increase timeclass, we copy the cell's ghost population of tc+1 into its actual population
            // then we put its current population into a coarser ghost
            // basically swapping main population and one tc level finer ghost population
            // this assumes that higher level ghost exists
            // TODO add error handling and/or an alternate way to increase timeclass later 

            auto newCoarserPop = spatialCell->get_population(popID);
            auto newFinerPop = spatialCell->get_population(popID, spatialCell->parameters[CellParams::TIMECLASS]+1);

            if (spatialCell->parameters[CellParams::TIMECLASS] != P::currentMaxTimeclass) {
               // If the cell is not at the maximum timeclass, we can increase it
               //std::cerr << "Increasing timeclass for cell " << cell << " with tc " << spatialCell->parameters[CellParams::TIMECLASS] << " by one"<< "\n";
               //std::cerr << "current max timeclass is " << P::currentMaxTimeclass << "\n";
               spatialCell->parameters[CellParams::TIMECLASS] += 1;
               spatialCell->parameters[CellParams::TIMECLASSDT] = spatialCell->get_tc_dt();
            } else {

               // If the cell is already at the maximum timeclass, we must create a new timeclass one higher
               std::cerr << "Cell " << cell << " is already at the maximum timeclass, creating a new one" << "\n";
               std::cerr << "current max timeclass is " << P::currentMaxTimeclass << "\n";

               std::cerr << "this is not supported yet, aborting" << "\n";
               abort();

               additionalTimeclassCreated = true;
               P::currentMaxTimeclass += 1;
               spatialCell->parameters[CellParams::TIMECLASS] = P::currentMaxTimeclass;
            
               P::timeclassDt.resize(P::currentMaxTimeclass + 1);
               P::timeclassDt.end()[-1] = P::timeclassDt.end()[-2]/2.0;

               spatialCell->parameters[CellParams::TIMECLASSDT] = spatialCell->get_tc_dt();
            }

            spatialCell->set_population(newFinerPop, popID);
            spatialCell->set_ghost_population(newCoarserPop, popID, spatialCell->parameters[CellParams::TIMECLASS]-1);
            spatialCell->requested_timeclass_ghosts.insert(spatialCell->parameters[CellParams::TIMECLASS]-1);         
            spatialCell->requested_timeclass_copy_ghosts.insert(spatialCell->parameters[CellParams::TIMECLASS]-1);
            // change cell time
            spatialCell->parameters[CellParams::TIME_V] -= P::timeclassDt[spatialCell->parameters[CellParams::TIMECLASS]]*0.5;         
         }
      }

      prepareAMRLists(mpiGrid);

      //calculateAcceleration(mpiGrid, 0.5, true, cellsToIncreaseTimeclass);

      //std::cerr << "calling prepareAMRLists after increasing timeclass\n";
      //std::cerr << "current max timeclass is " << P::currentMaxTimeclass << "\n";
      // this might be overkill, but for initial testing
      // prepareAMRLists(mpiGrid);
      // calculateAcceleration(mpiGrid, 0.0);
      // calculateSpatialTranslation(mpiGrid, 0.0, false);

      //remove extra ghosts from accelerated cells

      for (size_t c=0; c<cellsToIncreaseTimeclass.size(); ++c) {
         const CellID cell = cellsToIncreaseTimeclass[c];
         SpatialCell* spatialCell = mpiGrid[cell];
         for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {

            spatialCell->requested_timeclass_ghosts.erase(spatialCell->parameters[CellParams::TIMECLASS]);
            spatialCell->requested_timeclass_copy_ghosts.erase(spatialCell->parameters[CellParams::TIMECLASS]);
            spatialCell->remove_ghost_population(popID, spatialCell->parameters[CellParams::TIMECLASS]);
         }
      }

   } else {
      std::cout << "not implemented yet, aborting...\n";
      abort();
   }

}


//calculates currentmaxtimeclass
void calculateGlobalTcVariables(Real fsdt, Real globalMaxDt) {

   //setting fsdt smaller by the buffer amount
   //fsdt = fsdt / pow(2, P::timeclassBuffer);

   // This is the full range of timeclasses that could be used based on the physical environment
   int timeclassRange = max(int(log2(globalMaxDt/fsdt)),0);
   if (timeclassRange < P::initialMaxTimeclass) {
      // TODO figure this out if needed
      //std::cerr << "timeclassrange (" << (timeclassRange) << ") bigger than initialmaxtimeclass (" << P::initialMaxTimeclass << "), aborting" << std::endl;
      //abort();
   }

   if (P::tc_test_type == 1) {
      P::currentMaxTimeclass = P::initialMaxTimeclass;
      return;
   }

   if(P::tcOverrideTimeclass > -1 && P::tc_test_type != 0 && P::tc_test_type != 6){
      //P::currentMaxTimeclass = min(P::initialMaxTimeclass,P::tcOverrideTimeclass);
      // if we want a special test, just set the current timeclass to the initial one, and trust the programmer knows what they are doing.
      P::currentMaxTimeclass = P::initialMaxTimeclass;
      return;
   }

   // ... and we need to clamp that with the parameter for number of MaxTimeclasses
   P::currentMaxTimeclass = min(P::initialMaxTimeclass, timeclassRange);

}


void initiateAllCellTimeclasses(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid) {

   if (P::tc_test_type == 0) {
      // normal case, assign timeclasses based on CFL
      const vector<CellID>& cells = getLocalCells();
   
      for (vector<CellID>::const_iterator cell_id=cells.begin(); cell_id!=cells.end(); ++cell_id) {

         SpatialCell* cell = mpiGrid[*cell_id];
         // if (cell->parameters[CellParams::MAXVDT] != 0.0) {
         //    cellMaxDt = min(cell->parameters[CellParams::MAXRDT], cell->parameters[CellParams::MAXVDT] * P::maxSlAccelerationSubcycles);
         // } else {
         //    cellMaxDt = cell->parameters[CellParams::MAXRDT];
         // }
         //std::cerr << "cellMaxDt for cell " << *cell_id << " is " << cellMaxDt << std::endl;
         cell->assignCellTimeclass();
      }
   } else if(P::tc_test_type == 1){

      // if (P::dynamicTimestep) {
      //    std::cerr << "using dynamic timestep and special test not supported, aborting...\n";
      //    abort();
      // }

      if (P::initialMaxTimeclass != 1) {
         std::cerr << "not supported, aborting...\n";
         abort();
      }

      if (P::dynamicTimestep) {
         updateTimeclassDts(P::timeclassDt[1]*0.5); // halve the timestep lenghts as we want the longest dt to still be viable
         P::dt = P::timeclassDt[P::currentMaxTimeclass];
      }
      
      // set cell TCs such that one half is tc0 and one half is tc1.
      auto cells = getLocalCells();
      for (vector<CellID>::const_iterator cell_id=cells.begin(); cell_id!=cells.end(); ++cell_id) {

         SpatialCell* cell = mpiGrid[*cell_id];
         if (cell->parameters[CellParams::XCRD] <= -100.0) {
            cell->parameters[CellParams::TIMECLASS] = 1;
            cell->parameters[CellParams::TIMECLASSDT] = P::timeclassDt[1];
         } else {
            cell->parameters[CellParams::TIMECLASS] = 0;
            cell->parameters[CellParams::TIMECLASSDT] = P::timeclassDt[0];
         }
      }

   }
   else if(P::tc_test_type == 2) { 

   //constant timeclass in whole simulation domain

      if (P::tcOverrideTimeclass == -1 || P::dynamicTimestep) {
         std::cerr << "please set timeclass for overriding and use static timestep...\n";
         abort();
      }

      auto cells = getLocalCells();
      for (vector<CellID>::const_iterator cell_id=cells.begin(); cell_id!=cells.end(); ++cell_id) {
         SpatialCell* cell = mpiGrid[*cell_id];
         
         cell->parameters[CellParams::TIMECLASS] = P::tcOverrideTimeclass;
         cell->parameters[CellParams::TIMECLASSDT] = P::timeclassDt[P::tcOverrideTimeclass];
      
      }

   } else if(P::tc_test_type == 3) { 

      // static TC sphere areas up to some R_E
      // hardcoded up to 4 different levels

      assert(P::tcStaticSphereRadiusLvl1 > 0.0 && P::tcStaticSphereRadiusLvl2 > 0.0 && P::tcStaticSphereRadiusLvl3 > 0.0);
      assert(P::tcStaticSphereRadiusLvl1 < P::tcStaticSphereRadiusLvl2 && P::tcStaticSphereRadiusLvl2 < P::tcStaticSphereRadiusLvl3);

      int myRank;
      MPI_Comm_rank(MPI_COMM_WORLD,&myRank);

      if (myRank == 0) {
         std::cout << "Using static pre-programmed TC sphere areas, three level values: " << P::tcStaticSphereRadiusLvl1 << ", " << P::tcStaticSphereRadiusLvl2 << ", " << P::tcStaticSphereRadiusLvl3 << std::endl;
      }

      auto cells = getLocalCells();
      for (vector<CellID>::const_iterator cell_id=cells.begin(); cell_id!=cells.end(); ++cell_id) {
         SpatialCell* cell = mpiGrid[*cell_id];

         // calculate position of cell center
         const Real cellRadius = sqrt(pow(cell->parameters[CellParams::XCRD]+0.5*cell->parameters[CellParams::DX],2) + pow(cell->parameters[CellParams::YCRD]+0.5*cell->parameters[CellParams::DY],2) + pow(cell->parameters[CellParams::ZCRD]+0.5*cell->parameters[CellParams::DZ],2));
         

         if (P::currentMaxTimeclass==0) {
            cell->parameters[CellParams::TIMECLASS] = 0;
            cell->parameters[CellParams::TIMECLASSDT] = P::timeclassDt[0];
         } else if (P::currentMaxTimeclass==1) {
            if (cellRadius < P::tcStaticSphereRadiusLvl1) {
               cell->parameters[CellParams::TIMECLASS] = 1;
               cell->parameters[CellParams::TIMECLASSDT] = P::timeclassDt[1];
            } else {
               cell->parameters[CellParams::TIMECLASS] = 0;
               cell->parameters[CellParams::TIMECLASSDT] = P::timeclassDt[0];
            }
         } else if (P::currentMaxTimeclass==2) {
            if (cellRadius < P::tcStaticSphereRadiusLvl1) {
               cell->parameters[CellParams::TIMECLASS] = 2;
               cell->parameters[CellParams::TIMECLASSDT] = P::timeclassDt[2];
            } else if (cellRadius < P::tcStaticSphereRadiusLvl2) {
               cell->parameters[CellParams::TIMECLASS] = 1;
               cell->parameters[CellParams::TIMECLASSDT] = P::timeclassDt[1];
            } else {
               cell->parameters[CellParams::TIMECLASS] = 0;
               cell->parameters[CellParams::TIMECLASSDT] = P::timeclassDt[0];
            }
         } else if (P::currentMaxTimeclass==3) {
            if (cellRadius < P::tcStaticSphereRadiusLvl1) {
               cell->parameters[CellParams::TIMECLASS] = 3;
               cell->parameters[CellParams::TIMECLASSDT] = P::timeclassDt[3];
            } else if (cellRadius < P::tcStaticSphereRadiusLvl2) {
               cell->parameters[CellParams::TIMECLASS] = 2;
               cell->parameters[CellParams::TIMECLASSDT] = P::timeclassDt[2];
            } else if (cellRadius < P::tcStaticSphereRadiusLvl3) {
               cell->parameters[CellParams::TIMECLASS] = 1;
               cell->parameters[CellParams::TIMECLASSDT] = P::timeclassDt[1];
            } else {
               cell->parameters[CellParams::TIMECLASS] = 0;
               cell->parameters[CellParams::TIMECLASSDT] = P::timeclassDt[0];
            }
         } else {
            std::cerr << "not supported, aborting...\n";
            abort();
         }
      }
   } else {
      
      std::cerr << "not supported tc test, aborting...\n";
      abort();
   }
   P::timeclassesInitialized = true;
}

//check that timeclass settings are sensible
void timeclassDebugAssertions(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid) {

   assert(P::currentMaxTimeclass >= 0 && P::initialMaxTimeclass >= 0 && "Current and initial max timeclass must be non-negative");
   // fair assumption that no more than 20 timeclass need to exist (in reality more like 10)
   assert(P::currentMaxTimeclass < 20 && P::initialMaxTimeclass < 20 && "Do you really need more than 20 timeclasses?");
   for (int i=0; i<=P::currentMaxTimeclass; ++i) {
      assert(P::timeclassDt[i] >= 0.0 && "Timeclass dt must be non-negative");
   }

   for (const auto& cell_id : getLocalCells()) {
      SpatialCell* cell = mpiGrid[cell_id];
      assert(cell->parameters[CellParams::TIMECLASS] >= 0 && cell->parameters[CellParams::TIMECLASS] <= P::currentMaxTimeclass && "Cell timeclass must be within valid range");
      // assert(cell->parameters[CellParams::TIMECLASSDT] == P::timeclassDt[cell->parameters[CellParams::TIMECLASS]] && "Cell timeclass dt must match global timeclass dt");
   }
}
