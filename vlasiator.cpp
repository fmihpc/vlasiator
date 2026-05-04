/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 * Copyright 2024 CSC - IT Center for Science
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
#include "spatial_cells/spatial_cell_cpu.hpp"
#include "vlasovsolver/cpu_acc_transform.hpp"
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <limits>
#include <type_traits>
#include <vector>
#include <sstream>
#include <ctime>

#ifdef _OPENMP
   #include <omp.h>
#endif

#ifdef USE_GPU
#include "arch/gpu_base.hpp"
#endif

#include <fsgrid.hpp>

#include "vlasovsolver/vlasovmover.h"
#include "vlasovsolver/vec.h"
#include "definitions.h"
#include "mpiconversion.h"
#include "logger.h"
#include "parameters.h"
#include "readparameters.h"
#include "spatial_cells/spatial_cell_wrapper.hpp"
#include "datareduction/datareducer.h"

#include "sysboundary/sysboundary.h"
#include "vlasovsolver/common_pitch_angle_diffusion.hpp"

#include "fieldtracing/fieldtracing.h"

#include "fieldsolver/fs_common.h"
#include "projects/project.h"
#include "grid.h"
#include "iowrite.h"
#include "ioread.h"
#include "memory_report.h"

#include "object_wrapper.h"
#include "velocity_mesh_parameters.h"
#include "fieldsolver/gridGlue.hpp"
#include "fieldsolver/derivatives.hpp"

#include <signal.h>

#ifdef CATCH_FPE
#include <fenv.h>
/*! Function used to abort the program upon detecting a floating point exception. Which exceptions are caught is defined using the function feenableexcept.
 */
void fpehandler(int sig_num)
{
   signal(SIGFPE, fpehandler);
   printf("SIGFPE: floating point exception occured, exiting.\n");
   abort();
}
#endif

#include "phiprof.hpp"

Logger logFile, diagnostic;

using namespace std;

int globalflags::bailingOut = 0;
bool globalflags::writeRestart = false;
bool globalflags::writeRecover = false;
bool globalflags::balanceLoad = false;
bool globalflags::doRefine = false;
bool globalflags::ionosphereJustSolved = false;

#ifdef CATCH_SIGTERM
// The normal behaviour on SIGTERM is to simply abort the simulation in place.
// This implementation instead attempts to write a restart file and then quit,
// to work nicely with slurm's job preemption mechanism.
void termhandler(int sig_num) {
   logFile << "Caught SIGTERM. Writing recover and initiating bailout." << endl << flush;
   globalflags::bailingOut = 1;
   globalflags::writeRecover = 1;
}
#endif

ObjectWrapper objectWrapper;

ObjectWrapper& getObjectWrapper() {
   return objectWrapper;
}

/** Get local cell IDs. This function creates a cached copy of the 
 * cell ID lists to significantly improve performance. The cell ID 
 * cache is recalculated every time the mesh partitioning changes.
 * @return Local cell IDs.*/
const std::vector<CellID>& getLocalCells() {
   return Parameters::localCells;
}

void addTimedBarrier(string name){
#ifndef DEBUG_VLASIATOR
//let's not do a barrier
   return;
#endif
   phiprof::Timer btimer {name, {"Barriers", "MPI"}};
   MPI_Barrier(MPI_COMM_WORLD);
}

inline bool isDtTooLarge(Real dt, Real rdt, Real vdt, Real fsdt){
   return (dt > rdt * P::vlasovSolverMaxCFL ||
           dt > vdt * P::vlasovSolverMaxCFL * P::maxSlAccelerationSubcycles ||
           dt > fsdt * P::fieldSolverMaxCFL * P::maxFieldSolverSubcycles);
}

inline bool isDtTooSmall(Real dt, Real rdt, Real vdt, Real fsdt){
   return (dt < rdt * P::vlasovSolverMinCFL &&
           dt < vdt * P::vlasovSolverMinCFL * P::maxSlAccelerationSubcycles &&
           dt < fsdt * P::fieldSolverMinCFL * P::maxFieldSolverSubcycles);
}

// assume that maxrdt and maxvdt are updated before calling
bool cellTimeclassIsCorrect(SpatialCell* cell) {

   Real cellDt;
   if (cell->parameters[CellParams::MAXVDT] != 0.0) {
      cellDt = min(cell->parameters[CellParams::MAXRDT], cell->parameters[CellParams::MAXVDT] * P::maxSlAccelerationSubcycles);
   } else {
      cellDt = cell->parameters[CellParams::MAXRDT];
   }

   //std::cerr << "comparing celldt " << cellDt << " and timeclassdt " << P::timeclassDt[cell->parameters[CellParams::TIMECLASS]] << " for cell " << cell->get_cellid() << std::endl;
   //std::cerr << "their ratio: " << cellDt / P::timeclassDt[cell->parameters[CellParams::TIMECLASS]] << std::endl;
   if (cellDt > P::timeclassDt[cell->parameters[CellParams::TIMECLASS]]) {
      //std::cerr << "cell timeclass is correct" << std::endl;
      return true;
   } else {
      //std::cerr << "bad cell found!" << std::endl;
      //std::cerr << "cells sysboundaryflag and sysboundarylayer: " << cell->sysBoundaryFlag << ", " << cell->sysBoundaryLayer << std::endl;
      return false;
   }
}

// checks if cells' boundarytype is such that it should be taken into consideration for 
// timestep limiting and such
// checks taken from reduce_vlasov_dt

// should be changed into a member function
bool cellIsTimeclassRelevant(SpatialCell* cell) {

   // relevancy for acceleration
   if (!(cell->parameters[CellParams::MAXVDT] != 0 &&
      (cell->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY ||
      (P::vlasovAccelerateMaxwellianBoundaries && cell->sysBoundaryFlag == sysboundarytype::MAXWELLIAN)))) {
         return false;
      }

   // relevancy for translation
   if (!(cell->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY ||
      (cell->sysBoundaryLayer == 1 && cell->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY))) {
         return false;
      }

   return true;
}

// returns empty vector if all timeclasses are fine (= their timestep fits their timeclass)
// if not, returns those cells which need a bigger timeclass
// recalculates all cellwise tc limits

// skips over cells that are certain boundaries as defined above, as those should not affect timestep length
std::vector<CellID> checkCellTimeclasses(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid) {

   std::vector<CellID> retVec = {};
   const vector<CellID>& cells = getLocalCells();

   for (vector<CellID>::const_iterator cell_id=cells.begin(); cell_id!=cells.end(); ++cell_id) {
   
      if ((cellIsTimeclassRelevant(mpiGrid[*cell_id]))) {
         if (!(cellTimeclassIsCorrect(mpiGrid[*cell_id]))) {
            retVec.push_back(*cell_id);
         }
      }
   }

   return retVec;
}

// should be a member function of SC class
// sets cell parameters
void assignCellTimeclass(SpatialCell* cell, const double cellDt) {

   double baseTcDt = P::timeclassDt[P::currentMaxTimeclass - P::timeclassBuffer];

   if (P::tcOverrideTimeclass > -1) {
      cell->parameters[CellParams::TIMECLASS] = P::tcOverrideTimeclass;
      cell->parameters[CellParams::TIMECLASSDT] = cell->get_tc_dt();
      return;
   }

   if (cell->sysBoundaryFlag == sysboundarytype::COPYSPHERE || 
      cell->sysBoundaryFlag == sysboundarytype::IONOSPHERE ||
      cell->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) { // Copysphere and ionosphere cells always use the maximum timeclass
      cell->parameters[CellParams::TIMECLASS] = P::currentMaxTimeclass - P::timeclassBuffer;
      cell->parameters[CellParams::TIMECLASSDT] = cell->get_tc_dt();
      return;
   }

   // should this be a ceiling instead of floor??
   double dtdiff = int(log2((cellDt * P::timeclassDtModifier)/baseTcDt));
   int cellTimeClass = max(0.0,(P::currentMaxTimeclass - P::timeclassBuffer) - max(0.0, dtdiff));

   //std::cout << "assigning tc " << cellTimeClass << " for cell " << cell->get_cellid() << " with tcdt " << P::timeclassDt[cellTimeClass] << std::endl;

   cell->parameters[CellParams::TIMECLASS] = cellTimeClass;
   cell->parameters[CellParams::TIMECLASSDT] = cell->get_tc_dt();

}

void updateCellDtLimits() {


}

// goes through all cells, and sets their timeclasses according to some baseDt. also sets all timeclass--related cell parameters
void assingCellTimeclasses(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid) {

   const vector<CellID>& cells = getLocalCells();
   
   double cellMaxDt;
   for (vector<CellID>::const_iterator cell_id=cells.begin(); cell_id!=cells.end(); ++cell_id) {

      SpatialCell* cell = mpiGrid[*cell_id];
      if (cell->parameters[CellParams::MAXVDT] != 0.0) {
         cellMaxDt = min(cell->parameters[CellParams::MAXRDT], cell->parameters[CellParams::MAXVDT] * P::maxSlAccelerationSubcycles);
      } else {
         cellMaxDt = cell->parameters[CellParams::MAXRDT];
      }
      //std::cerr << "cellMaxDt for cell " << *cell_id << " is " << cellMaxDt << std::endl;
      assignCellTimeclass(cell, cellMaxDt);
   }
}

void updateTimeclassDts(Real fsdt) {

   // reduce fsdt by buffer amount

   fsdt /= pow(2.0, P::timeclassBuffer);

   std::vector<Real> newTimeclassDts(P::currentMaxTimeclass+1);
   logFile << std::endl;
   logFile << "(TC) timeclassDts set to " << std::endl;
   for(int i = 0; i <= P::currentMaxTimeclass; ++i){
      newTimeclassDts[i] = fsdt*pow(2,P::currentMaxTimeclass - i)*P::timeclassDtModifier;
      logFile << newTimeclassDts[i] << "s, ";
   }
   logFile << std::endl;
   logFile << std::endl;
   P::timeclassDt = newTimeclassDts;

}

void increaseTimeclass(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                              const std::vector<CellID>& cellsToIncreaseTimeclass,
                              bool& additionalTimeclassCreated) {
   phiprof::Timer increaseTimeclassTimer {"increase-timeclass"};

   additionalTimeclassCreated = false;

   // Increase timeclass for given cells

   for (size_t c=0; c<cellsToIncreaseTimeclass.size(); ++c) {
      const CellID cell = cellsToIncreaseTimeclass[c];
      SpatialCell* spatialCell = mpiGrid[cell];
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
   }

   //std::cerr << "calling prepareAMRLists after increasing timeclass\n";
   //std::cerr << "current max timeclass is " << P::currentMaxTimeclass << "\n";
   // this might be overkill, but for initial testing
   prepareAMRLists(mpiGrid);
   calculateAcceleration(mpiGrid, 0.0);


   auto lCells = getLocalCells();
  // std::cerr << "timeclass increase done, now cells have following ghost distributions:\n";
   for (size_t c=0; c<lCells.size(); ++c) {
      const CellID cell = lCells[c];
      SpatialCell* spatialCell = mpiGrid[cell];
      //std::cerr << "cell " << cell << " has timeclass " << spatialCell->parameters[CellParams::TIMECLASS] << " and ghost distribution: ";
      for (auto& ghost : spatialCell->requested_timeclass_ghosts) {
         //std::cerr << ghost << " ";
      }
      //std::cerr << "\n";

   }

}

// returns vector of timestep values

std::vector<Real> computeNewTimeStep(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
			FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
         std::vector<Real>& dtMaxLocal, std::vector<Real>& dtMaxGlobal, 
         std::vector<Real>& dtMinMaxLocal, std::vector<Real>& dtMinMaxGlobal) {

   phiprof::Timer computeTimestepTimer {"compute-timestep"};
   // Compute maximum time step. This cannot be done at the first step as the solvers compute the limits for each cell.

   const vector<CellID>& cells = getLocalCells();
   // newTimeclassDts = std::vector<Real>(P::maxTimeclass+1);

   // Compute max dt for Vlasov solver
   reduce_vlasov_dt(mpiGrid, cells, dtMaxLocal, dtMinMaxLocal);

   // compute max dt for fieldsolver
   const std::array<FsGridTools::FsIndex_t, 3> gridDims(technicalGrid.getLocalSize());
   for (FsGridTools::FsIndex_t k = 0; k < gridDims[2]; k++) {
      for (FsGridTools::FsIndex_t j = 0; j < gridDims[1]; j++) {
         for (FsGridTools::FsIndex_t i = 0; i < gridDims[0]; i++) {
            fsgrids::technical* cell = technicalGrid.get(i, j, k);
            if (cell->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY ||
               (cell->sysBoundaryLayer == 1 && cell->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY)) {
               dtMaxLocal[2] = min(dtMaxLocal[2], cell->maxFsDt);
               dtMinMaxLocal[2] = max(dtMinMaxLocal[2], cell->maxFsDt);
            }
         }
      }
   }
   int myRank;MPI_Comm_rank(MPI_COMM_WORLD,&myRank);

   MPI_Allreduce(dtMaxLocal.data(), dtMaxGlobal.data(), 3, MPI_Type<Real>(), MPI_MIN, MPI_COMM_WORLD);
   MPI_Allreduce(dtMinMaxLocal.data(), dtMinMaxGlobal.data(), 3, MPI_Type<Real>(), MPI_MAX, MPI_COMM_WORLD);

   // If any of the solvers are disabled there should be no limits in timespace from it
   if (!P::propagateVlasovTranslation)
      dtMaxGlobal[0] = numeric_limits<Real>::max();
   if (!P::propagateVlasovAcceleration)
      dtMaxGlobal[1] = numeric_limits<Real>::max();
   if (!P::propagateField)
      dtMaxGlobal[2] = numeric_limits<Real>::max();

   creal meanVlasovCFL = 0.5 * (P::vlasovSolverMaxCFL + P::vlasovSolverMinCFL);
   creal meanFieldsCFL = 0.5 * (P::fieldSolverMaxCFL + P::fieldSolverMinCFL);

   Real localDt, baseDt, fsdt;

   // localDt: max dt in the local MPI domain
   localDt = meanVlasovCFL * dtMaxLocal[0];
   localDt = min(localDt,meanVlasovCFL * dtMaxLocal[1] * P::maxSlAccelerationSubcycles);
   localDt = min(localDt,meanFieldsCFL * dtMaxLocal[2] * P::maxFieldSolverSubcycles);
   if (myRank == MASTER_RANK) cout << "localDt " << localDt <<"\n";
   // newDt: max dt globally, at the highest timeclass
   fsdt = meanVlasovCFL * dtMaxGlobal[0];
   fsdt = min(fsdt,meanVlasovCFL * dtMaxGlobal[1] * P::maxSlAccelerationSubcycles);
   fsdt = min(fsdt,meanFieldsCFL * dtMaxGlobal[2] * P::maxFieldSolverSubcycles);
   if (myRank == MASTER_RANK) cout << "fsdt " << fsdt <<"\n";
   // baseDt: longest max dt of any rank
   baseDt = meanVlasovCFL * dtMinMaxGlobal[0];
   baseDt = min(baseDt,meanVlasovCFL * dtMinMaxGlobal[1] * P::maxSlAccelerationSubcycles);
   baseDt = min(baseDt,meanFieldsCFL * dtMinMaxGlobal[2] * P::maxFieldSolverSubcycles);   
   if (myRank == MASTER_RANK) cout << "baseDt " << baseDt <<"\n";

   std::vector<Real> retVec = {localDt, fsdt, baseDt};

   return retVec;
}


// // called when a cell changes its timeclass. checks if that cells' neighbors possess all required timeghosts. If not, returns vector of pairs, each containing a cellid and timeclass of timeghost needed to be created.
// std::vector<std::pair<CellID, int>> checkNeighborTimeghosts(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid, SpatialCell* cell) {

//    std::vector<std::pair<CellID, int>> ret;

//    auto neighbors = mpiGrid.get_neighbors_of(cell->get_cellid(), Neighborhoods::VLASOV_SOLVER_TIMEGHOST_EXACT_HALO);

//    for (auto n : *neighbors){
//       CellID cid = n.first;

      

//    }


// }

//creates timeghosts in cells where they were requested
void createNewTimeghosts(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid, CellID cell) {

   getGhostNeighborsforTC(mpiGrid, {cell}, std::set<CellID> &active_cells, int timeclass)

}

// check goodness of current used fsdt, if it isnt good, changes newDt to good one and sets isChanged to true. Also sets subcycling number.
void handleChangingofDt(const std::vector<Real>& dtMaxGlobal, bool& isChanged, Real& newDt) {

   int myRank;MPI_Comm_rank(MPI_COMM_WORLD,&myRank);

   creal meanVlasovCFL = 0.5 * (P::vlasovSolverMaxCFL + P::vlasovSolverMinCFL);
   creal meanFieldsCFL = 0.5 * (P::fieldSolverMaxCFL + P::fieldSolverMinCFL);

   Real subcycleDt;

      // reduce/increase dt if it is too high for any of the three propagators or too low for all propagators
   if (isDtTooLarge(P::timeclassDt[P::currentMaxTimeclass - P::timeclassBuffer], dtMaxGlobal[0],dtMaxGlobal[1],dtMaxGlobal[2]) ||
          isDtTooSmall(P::timeclassDt[P::currentMaxTimeclass - P::timeclassBuffer], dtMaxGlobal[0],dtMaxGlobal[1],dtMaxGlobal[2])) {

      if (P::staticTimeclasses) {
         std::cerr << "aborting since timestep wants to change with static timeclasses" << std::endl;
         std::cerr << "if this happens at the beginning of the run, then try to tweak P::timeclassDtModifier" << std::endl;
         abort();
      }

      // new dt computed
      isChanged = true;

      // set new timestep to the lowest one of all interval-midpoints
      newDt = meanVlasovCFL * dtMaxGlobal[0];
      newDt = min(newDt, meanVlasovCFL * dtMaxGlobal[1] * P::maxSlAccelerationSubcycles);
      newDt = min(newDt, meanFieldsCFL * dtMaxGlobal[2] * P::maxFieldSolverSubcycles);

      logFile << "(TIMESTEP) New dt = " << newDt << " computed on step " << P::tstep << " at " << P::t
              << "s   Maximum possible dt (not including  vlasovsolver CFL " << P::vlasovSolverMinCFL << "-"
              << P::vlasovSolverMaxCFL << " or fieldsolver CFL " << P::fieldSolverMinCFL << "-" << P::fieldSolverMaxCFL
              << ") in {r, v, BE} was " << dtMaxGlobal[0] << " " << dtMaxGlobal[1] << " " << dtMaxGlobal[2] << " "
              << " Including subcycling { v, BE}  was " << dtMaxGlobal[1] * P::maxSlAccelerationSubcycles << " "
              << dtMaxGlobal[2] * P::maxFieldSolverSubcycles << " " << endl
              << writeVerbose;

      if (P::dynamicTimestep) {
         // Check if the calculated value was and continues to be above the ceiling
         if (P::dt_ceil > 0.0 && newDt >= P::dt_ceil && P::dt == P::dt_ceil) {
            isChanged = false;
            newDt = P::dt_ceil;
            return;
         }
         // Check if we at this time exceeded the ceiling
         if (P::dt_ceil > 0.0 && newDt > P::dt_ceil) {
            newDt = P::dt_ceil;
            logFile << "(TIMESTEP) However, ceiling timestep in config overrides larger dynamic and dt = " << P::dt_ceil << endl << writeVerbose;
         } 
         subcycleDt = newDt;
      } else {
         logFile << "(TIMESTEP) However, fixed timestep in config overrides and dt = " << P::dt << endl << writeVerbose;
         subcycleDt = P::dt;
      }
   } else {
      subcycleDt = P::dt;
   }

   // Subcycle if field solver dt < global dt (including CFL) (new or old dt hence the hassle with subcycleDt
   if (meanFieldsCFL * dtMaxGlobal[2] < subcycleDt && P::propagateField) {
      P::fieldSolverSubcycles =
          min(convert<uint>(ceil(subcycleDt / (meanFieldsCFL * dtMaxGlobal[2]))), P::maxFieldSolverSubcycles);
   } else {
      P::fieldSolverSubcycles = 1;
   }
}

//calculates currentmaxtimeclass
void calculateGlobalTcVariables(Real fsdt, Real globalMaxDt) {

   //setting fsdt smaller by the buffer amount
   //fsdt = fsdt / pow(2, P::timeclassBuffer);

   // This is the full range of timeclasses that could be used based on the physical environment
   int timeclassRange = int(log2(globalMaxDt/fsdt));

   if (timeclassRange < P::initialMaxTimeclass) {
      // TODO figure this out if needed
      //std::cerr << "timeclassrange (" << (timeclassRange) << ") bigger than initialmaxtimeclass (" << P::initialMaxTimeclass << "), aborting" << std::endl;
      //abort();
   }

   // ... and we need to clamp that with the parameter for number of MaxTimeclasses
   P::currentMaxTimeclass = min(P::initialMaxTimeclass, timeclassRange);

   if(P::tcOverrideTimeclass > -1){
      std::cerr << "Setting all tc to " << P::tcOverrideTimeclass << "\n";
      P::currentMaxTimeclass = min(P::initialMaxTimeclass,P::tcOverrideTimeclass);
   }
}

// sets all simulation cells to some timeclass, calculates currentMaxTimeClass, calculates timeclassDts, sets dt
void initiateAllCellTimeclasses(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid) {

   // if(P::tcDebugBox){
   //    P::currentMaxTimeclass = 1;
   // }

   // if (myRank == MASTER_RANK) cout << "dtrange: " << timeclassRange << ", newDt = " << newDt <<
   //  ", baseDt = " << baseDt << ", fsdt " << fsdt << ", current max tc "<< P::currentMaxTimeclass<< std::endl;
   // baseDt = fsdt*pow(2, P::currentMaxTimeclass);
   // if (myRank == MASTER_RANK) cout << "for new baseDt = " << baseDt << std::endl;
   
   // TODO handle changing P::currentMaxTimeclass!

   // We need dts and timeclasses relative to the shortest viable maxDt:
   // int dtdiff = int(log2(localDt/fsdt));
   // int localTimeClass = max(0,P::currentMaxTimeclass - max(0, dtdiff)); // this shouldn't actually matter anymore?
   
   // if(P::tcDebugBox){
   //    localTimeClass = 0;
   //    for (vector<CellID>::const_iterator cell_id=cells.begin(); cell_id!=cells.end(); ++cell_id) {
   //       SpatialCell* cell = mpiGrid[*cell_id];
   //       const Real x = cell->parameters[CellParams::XCRD];
   //       const Real y = cell->parameters[CellParams::YCRD];
   //       const Real z = cell->parameters[CellParams::ZCRD];
   //       if (abs(x-P::tcBoxCenterX) < P::tcBoxHalfWidthX &&
   //           abs(y-P::tcBoxCenterY) < P::tcBoxHalfWidthY &&
   //           abs(z-P::tcBoxCenterZ) < P::tcBoxHalfWidthZ ){
   //          cell->parameters[CellParams::TIMECLASS] = 1;
   //          localTimeClass = 1;
   //       }
   //       else{
   //          cell->parameters[CellParams::TIMECLASS] = 0;
   //          localTimeClass = max(0, localTimeClass);
   //       }
   //    }
      
   //    for (vector<CellID>::const_iterator cell_id=cells.begin(); cell_id!=cells.end(); ++cell_id) {
   //       SpatialCell* cell = mpiGrid[*cell_id];
   //       if(P::tcRankwise){
   //          cell->parameters[CellParams::TIMECLASS] = localTimeClass;
   //       }
   //       cell->parameters[CellParams::TIMECLASSDT] = cell->get_tc_dt();
   //    }
   
   // }
   if(P::tc_test_type == 1){
      // if(P::tcOverrideTimeclass > -1){
      //    localTimeClass = P::tcOverrideTimeclass;
      // }
      // else {
      //    localTimeClass = 0;   // For the MPI-rank based timeclasses. Implement to CellParams if cell-based.
   // Move to params.
      // }
      // P::currentMaxTimeclass = P::initialMaxTimeclass;

      // // for(int i = 0; i <= P::initialMaxTimeclass; ++i){
      // //    newTimeclassDts[i] = fsdt*pow(2,P::currentMaxTimeclass - min(i,P::currentMaxTimeclass));
      // // }
      // // P::timeclassDt = newTimeclassDts;
      
      // for (vector<CellID>::const_iterator cell_id=cells.begin(); cell_id!=cells.end(); ++cell_id) {
      //    SpatialCell* cell = mpiGrid[*cell_id];

      //    cell->parameters[CellParams::TIMECLASS] = min(localTimeClass, P::initialMaxTimeclass);
      //    cell->parameters[CellParams::TIMECLASSDT] = cell->get_tc_dt();
      // }
   }
   else if(P::tc_test_type == 2 || P::tc_test_type == 3){ 
      // std::cerr << "TC test 2\n";
      // if(P::initialMaxTimeclass > 2){
      //    std::cerr << "This test works best with timeclass 1 or 2\n";
      //    // abort();
      // }
      // if(P::initialMaxTimeclass > 0) {
      //    P::currentMaxTimeclass = P::initialMaxTimeclass;
      // }
      // else{
      //    P::currentMaxTimeclass = 0;
      // }
      // // for(int i = 0; i <= P::initialMaxTimeclass; ++i){
      // //    newTimeclassDts[i] = fsdt*pow(2,P::currentMaxTimeclass - min(i,P::currentMaxTimeclass));
      // // }
      // // P::timeclassDt = newTimeclassDts;
      // // if(P::tcOverrideTimeclass > -1){
      // //    localTimeClass = P::tcOverrideTimeclass;
      // // }
      // // else {
      // //    localTimeClass = myRank % 2;
      // // }
      
      // for (vector<CellID>::const_iterator cell_id=cells.begin(); cell_id!=cells.end(); ++cell_id) {
      //    SpatialCell* cell = mpiGrid[*cell_id];
      //    // cell->parameters[CellParams::TIMECLASS] = min(localTimeClass, P::maxTimeclass);

      //    // first block: one half timeclass0 and other timeclassmax
      //    // second block: three parts: timeclassmax-2, timeclassmax-1, timeclassmax
      //    if (P::initialMaxTimeclass == 1) {
      //       cell->parameters[CellParams::TIMECLASS] = min(int(cell->parameters[CellParams::XCRD] > -100/*epsilon*/)*P::initialMaxTimeclass, P::initialMaxTimeclass);
      //    } else if (P::initialMaxTimeclass == 2) {
      //       if (cell->parameters[CellParams::XCRD] < -15*(cell->parameters[CellParams::DX])) {
      //          cell->parameters[CellParams::TIMECLASS] = 0;
      //       } else if (cell->parameters[CellParams::XCRD] > 15*(cell->parameters[CellParams::DX])) {
      //          cell->parameters[CellParams::TIMECLASS] = 2;
      //       } else {
      //          cell->parameters[CellParams::TIMECLASS] = 1;
      //       }
      //    } else if (P::initialMaxTimeclass == 3) {
      //       cell->parameters[CellParams::TIMECLASS] = min(int(cell->parameters[CellParams::XCRD] > -100/*epsilon*/)*P::initialMaxTimeclass, P::initialMaxTimeclass);
      //    }
      //    cell->parameters[CellParams::TIMECLASSDT] = cell->get_tc_dt();
      // }
      
   } else if (P::tc_test_type ==4) {
      // //for 2d testing with tc box in the middle

      // int sidelen = P::xcells_ini;

      // if(P::initialMaxTimeclass > 0) {
      //    P::currentMaxTimeclass = P::initialMaxTimeclass;
      // }
      // else{
      //    P::currentMaxTimeclass = 0;
      // }
      // for(int i = 0; i <= P::initialMaxTimeclass; ++i){
      //    newTimeclassDts[i] = fsdt*pow(2,P::currentMaxTimeclass - min(i,P::currentMaxTimeclass));
      // }
      // P::timeclassDt = newTimeclassDts;
      // // if(P::tcOverrideTimeclass > -1){
      // //    localTimeClass = P::tcOverrideTimeclass;
      // // }
      // // else {
      // //    localTimeClass = myRank % 2;
      // // }
      
      // for (vector<CellID>::const_iterator cell_id=cells.begin(); cell_id!=cells.end(); ++cell_id) {
      //    SpatialCell* cell = mpiGrid[*cell_id];

      //    if (cell->parameters[CellParams::XCRD] > (sidelen/4.0)*cell->parameters[CellParams::DX] && cell->parameters[CellParams::XCRD] < (3.0*sidelen/4.0)*cell->parameters[CellParams::DX] &&
      //    cell->parameters[CellParams::YCRD] > (sidelen/4.0)*cell->parameters[CellParams::DX] && cell->parameters[CellParams::YCRD] < (3.0*sidelen/4.0)*cell->parameters[CellParams::DX]) {
      //       cell->parameters[CellParams::TIMECLASS] = 1;   
      //    } else {
      //       cell->parameters[CellParams::TIMECLASS] = 0;
      //    }
      //    cell->parameters[CellParams::TIMECLASSDT] = cell->get_tc_dt();
      // }


   } else if (P::tc_test_type == 5) {
      // //for 2d testing with tc box in the middle, only with one side longer than the other

      // int sidelenX = P::xcells_ini;
      // int sidelenY = P::ycells_ini;

      // if(P::initialMaxTimeclass > 0) {
      //    P::currentMaxTimeclass = P::initialMaxTimeclass;
      // }
      // else{
      //    P::currentMaxTimeclass = 0;
      // }
      // for(int i = 0; i <= P::initialMaxTimeclass; ++i){
      //    newTimeclassDts[i] = fsdt*pow(2,P::currentMaxTimeclass - min(i,P::currentMaxTimeclass));
      // }
      // P::timeclassDt = newTimeclassDts;
      // // if(P::tcOverrideTimeclass > -1){
      // //    localTimeClass = P::tcOverrideTimeclass;
      // // }
      // // else {
      // //    localTimeClass = myRank % 2;
      // // }
      
      // for (vector<CellID>::const_iterator cell_id=cells.begin(); cell_id!=cells.end(); ++cell_id) {
      //    SpatialCell* cell = mpiGrid[*cell_id];

      //    if (cell->parameters[CellParams::XCRD] > (sidelenX/4.0)*cell->parameters[CellParams::DX] && cell->parameters[CellParams::XCRD] < (3.0*sidelenX/4.0)*cell->parameters[CellParams::DX] &&
      //    cell->parameters[CellParams::YCRD] > (sidelenY/4.0)*cell->parameters[CellParams::DX] && cell->parameters[CellParams::YCRD] < (3.0*sidelenY/4.0)*cell->parameters[CellParams::DX]) {
      //       cell->parameters[CellParams::TIMECLASS] = 1;   
      //    } else {
      //       cell->parameters[CellParams::TIMECLASS] = 0;
      //    }
      //    cell->parameters[CellParams::TIMECLASSDT] = cell->get_tc_dt();
      // }

   
   } else if (P::tc_test_type == 8) {

      // normal operation except cells are assigned with a modifier to cellwise dt limit
      //std::cout << "dt modifier: " << P::timeclassDtModifier << std::endl;
      assingCellTimeclasses(mpiGrid);         

   } else {

      // For normal operation, set the timeclass and timeclassdt for each cell according to their timestep lenght calculated by reducers.

      assingCellTimeclasses(mpiGrid);         

   }
}

// void getGhostNeighborsforTC(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
//                               const std::vector<CellID>& cellsToCheckNeighbors) {
//    /*
//    1st version
//    every timestep, go through every cell c, and get its ghost neighbours. 
//    Then, for every ghost neighbour, send c's timeclass to its requested_timeclass_ghosts
//    */
//    /*
//    2nd version TODO:
//    every timestep, check if computeNewTimestep changes any cells' timeclass. Then go through v1 functionality.
//    */

//    for (size_t c=0; c<cellsToCheckNeighbors.size(); ++c) {
//       const CellID cell = cellsToCheckNeighbors[c];
//       auto neighbors = mpiGrid.get_neighbors_of(cell, VLASOV_SOLVER_GHOST_NEIGHBORHOOD_ID);
//       auto& neighborsRef = *neighbors;
//       auto neighborsRemote = mpiGrid.get_remote_neighbors_of(cell, VLASOV_SOLVER_GHOST_NEIGHBORHOOD_ID);

//       // get_neighbours_of returns a pointer to a vector of pairs, and each pairs' first element is the CellID
//       // get_remote_neighbors_of returns a vector of CellIDs

//       for (size_t i=0; i<neighborsRef.size(); ++i) {
//          if (mpiGrid[(neighborsRef)[i].first]->parameters[CellParams::TIMECLASS] != mpiGrid[cell]->parameters[CellParams::TIMECLASS]) {
//             mpiGrid[(neighborsRef)[i].first]->requested_timeclass_ghosts.insert(mpiGrid[cell]->parameters[CellParams::TIMECLASS]);
//          }
//       }
//       for (size_t i=0; i<neighborsRemote.size(); ++i) {
//          if (mpiGrid[(neighborsRemote)[i]]->parameters[CellParams::TIMECLASS] != mpiGrid[cell]->parameters[CellParams::TIMECLASS]) {
//             mpiGrid[neighborsRemote[i]]->requested_timeclass_ghosts.insert(mpiGrid[cell]->parameters[CellParams::TIMECLASS]);
//          }
//       }
//    }
// }

int simulate(int argn,char* args[]) {
   int myRank, doBailout=0;
   const creal DT_EPSILON=1e-12;
   typedef Parameters P;
   Real newDt;
   bool dtIsChanged {false};
   bool additionalTimeclassCreated {false};
   bool aCellHadTimeclassChanged {false};

   /* Arrays for storing local (per process) and global max dt
   0th position stores ordinary space propagation dt
   1st position stores velocity space propagation dt
   2nd position stores field propagation dt
   */
   std::vector<Real> dtMaxLocal(3);
   std::vector<Real> dtMaxGlobal(3);
   std::vector<Real> dtMinMaxLocal(3);
   std::vector<Real> dtMinMaxGlobal(3);
   
   dtMaxLocal[0] = numeric_limits<Real>::max();
   dtMaxLocal[1] = numeric_limits<Real>::max();
   dtMaxLocal[2] = numeric_limits<Real>::max();

   dtMaxGlobal[0] = numeric_limits<Real>::max();
   dtMaxGlobal[1] = numeric_limits<Real>::max();
   dtMaxGlobal[2] = numeric_limits<Real>::max();

   dtMinMaxLocal[0] = numeric_limits<Real>::min();
   dtMinMaxLocal[1] = numeric_limits<Real>::min();
   dtMinMaxLocal[2] = numeric_limits<Real>::min();

   dtMinMaxGlobal[0] = numeric_limits<Real>::min();
   dtMinMaxGlobal[1] = numeric_limits<Real>::min();
   dtMinMaxGlobal[2] = numeric_limits<Real>::min();

   MPI_Comm_rank(MPI_COMM_WORLD,&myRank);

   phiprof::initialize();

   double initialWtime =  MPI_Wtime();
   SysBoundary& sysBoundaryContainer = getObjectWrapper().sysBoundaryContainer;
   
   #ifdef CATCH_FPE
   // WARNING FE_INEXACT is too sensitive to be used. See man fenv.
   //feenableexcept(FE_DIVBYZERO|FE_INVALID|FE_OVERFLOW|FE_UNDERFLOW);
   feenableexcept(FE_DIVBYZERO|FE_INVALID|FE_OVERFLOW);
   //feenableexcept(FE_DIVBYZERO|FE_INVALID);
   signal(SIGFPE, fpehandler);
   #endif

   #ifdef CATCH_SIGTERM
   signal(SIGTERM, termhandler);
   #endif

   // Initialize memory allocator configuration.
   memory_configurator();

   phiprof::Timer mainTimer {"main"};
   phiprof::Timer initTimer {"Initialization"};

   phiprof::Timer readParamsTimer {"Read parameters"};
   // Allocate host-side velocity mesh wrapper
   vmesh::allocateMeshWrapper();
   // init parameter file reader
   Readparameters readparameters(argn,args);

   P::addParameters();

   // Add parameters for number of populations
   getObjectWrapper().addParameters();
   readparameters.parse();
   P::getParameters();

   getObjectWrapper().addPopulationParameters();
   sysBoundaryContainer.addParameters();
   projects::Project::addParameters();

   Project* project = projects::createProject();
   getObjectWrapper().project = project;
   readparameters.parse(true, false); // 2nd parsing for specific population parameters
   readparameters.helpMessage(); // Call after last parse, exits after printing help if help requested
   getObjectWrapper().getPopulationParameters();
   sysBoundaryContainer.getParameters();
   project->getParameters();

   #ifdef USE_GPU
   // Activate device, create streams
   gpu_init_device();
   #endif
   // Fill in rest of velocity meshes data, upload GPU version
   vmesh::getMeshWrapper()->initVelocityMeshes(getObjectWrapper().particleSpecies.size());
   readParamsTimer.stop();

   // Check for correct application of vectorclass values:
   if ( (VECL<WID) ||
        (VECL*VEC_PER_PLANE != WID2) ||
        (VECL*VEC_PER_BLOCK != WID3) ||
        //(VPREC > VECL) ||
        (VECL != (int)VECL) ||
        (VPREC != (int)VPREC) ||
        (VEC_PER_PLANE != (int)VEC_PER_PLANE) ||
        (VEC_PER_BLOCK != (int)VEC_PER_BLOCK) ) {
      if (myRank == MASTER_RANK) {
         cerr << "(MAIN) ERROR: Vectorclass definition mismatch!" << endl;
         cerr << "VECL " << VECL <<" VEC_PER_PLANE " << VEC_PER_PLANE <<" WID " << WID <<" VEC_PER_BLOCK " << VEC_PER_BLOCK << " VPREC "<< VPREC<<endl;
      }
      exit(1);
   }

   if (P::initialMaxTimeclass > 0 && P::vlasovSolverGhostTranslate == false) {
      // if we are using timeclasses, we need to use ghost translation
      cerr << "(MAIN) Warning: Using timeclasses requires ghost translation, please turn GT on. exiting..." << endl;
      logFile << "(MAIN) Warning: Using timeclasses requires ghost translation, please turn GT on. exiting..." << endl;
      exit(1);
   }

   // Verify correct handling of floating point exceptions
   // see https://github.com/fmihpc/vlasiator/pull/845
   {
      double qnan = std::numeric_limits<double>::quiet_NaN();
      double pinf = std::numeric_limits<double>::infinity();
      double ninf = -std::numeric_limits<double>::infinity();
      bool isnan1 = std::isnan(qnan);
      bool isinf2 = std::isinf(pinf);
      bool isinf3 = std::isinf(ninf);
      bool isfinite1 = std::isfinite(qnan);
      bool isfinite2 = std::isfinite(pinf);
      bool isfinite3 = std::isfinite(ninf);
      if (!isnan1||!isinf2||!isinf3||isfinite1||isfinite2||isfinite3) {
      if (myRank == MASTER_RANK) {
         cerr << "(MAIN) ERROR: Floating point exceptions not being caught!" << endl;
      }
      exit(1);
      }
   }
   //Get version and config info here
   std::string version;
   std::string config;
   //Only master needs the info
   if (myRank==MASTER_RANK){
      version=readparameters.versionInfo();
      config=readparameters.configInfo();
   }

   // Init parallel logger:

   phiprof::Timer openLoggerTimer {"open logFile & diagnostic"};
   //if restarting we will append to logfiles
   if(!P::writeFullBGB) {
      if (logFile.open(MPI_COMM_WORLD,MASTER_RANK,"logfile.txt",P::isRestart) == false) {
         if(myRank == MASTER_RANK) cerr << "(MAIN) ERROR: Logger failed to open logfile!" << endl;
         exit(1);
      }
   } else {
      // If we are out to write the full background field and derivatives, we don't want to overwrite the existing run's logfile.
      if (logFile.open(MPI_COMM_WORLD,MASTER_RANK,"logfile_fullbgbio.txt",false) == false) {
         if(myRank == MASTER_RANK) cerr << "(MAIN) ERROR: Logger failed to open logfile_fullbgbio!" << endl;
         exit(1);
      }
   }
   if (P::diagnosticInterval != 0) {
      if (diagnostic.open(MPI_COMM_WORLD,MASTER_RANK,"diagnostic.txt",P::isRestart) == false) {
         if(myRank == MASTER_RANK) cerr << "(MAIN) ERROR: Logger failed to open diagnostic file!" << endl;
         exit(1);
      }
   }

   int mpiProcs;
   MPI_Comm_size(MPI_COMM_WORLD,&mpiProcs);

   char nodename[MPI_MAX_PROCESSOR_NAME]; 
   int namelength, nodehash;
   int nodeRank, interRank;
   int nNodes;

   hash<string> hasher; 
   MPI_Comm nodeComm;
   MPI_Comm interComm;
  
   //get name of this node
   MPI_Get_processor_name(nodename,&namelength);   
   nodehash=(int)(hasher(string(nodename)) % std::numeric_limits<int>::max());
   
   //intra-node communicator
   MPI_Comm_split(MPI_COMM_WORLD, nodehash, myRank, &nodeComm);
   MPI_Comm_rank(nodeComm,&nodeRank);
   //create communicator for inter-node communication
   MPI_Comm_split(MPI_COMM_WORLD, nodeRank, myRank, &interComm);
   MPI_Comm_rank(interComm, &interRank);
   MPI_Comm_size(interComm, &nNodes);

   MPI_Comm_free(&interComm);
   MPI_Comm_free(&nodeComm);

   logFile << "(MAIN) Starting simulation with " << mpiProcs << " MPI processes ";
   #ifdef _OPENMP
      logFile << "and " << omp_get_max_threads();
   #else
      logFile << "and 0";
   #endif
   logFile << " OpenMP threads per process on " << nNodes << " nodes" << endl << writeVerbose;      
   openLoggerTimer.stop();
   
   // Init project
   phiprof::Timer initProjectimer {"Init project"};
   if (project->initialize() == false) {
      if(myRank == MASTER_RANK) cerr << "(MAIN): Project did not initialize correctly!" << endl;
      exit(1);
   }
   if (project->initialized() == false) {
      if (myRank == MASTER_RANK) {
         cerr << "(MAIN): Project base class was not initialized!" << endl;
         cerr << "\t Call Project::initialize() in your project's initialize()-function." << endl;
         exit(1);
      }
   }
   initProjectimer.stop();

   // Initialize simplified Fieldsolver grids.
   // Needs to be done here already ad the background field will be set right away, before going to initializeGrid even
   phiprof::Timer initFsTimer {"Init fieldsolver grids"};

   std::array<FsGridTools::FsSize_t, 3> fsGridDimensions = {convert<FsGridTools::FsSize_t>(P::xcells_ini * pow(2,P::amrMaxSpatialRefLevel)),
							    convert<FsGridTools::FsSize_t>(P::ycells_ini * pow(2,P::amrMaxSpatialRefLevel)),
							    convert<FsGridTools::FsSize_t>(P::zcells_ini * pow(2,P::amrMaxSpatialRefLevel))};

   std::array<bool,3> periodicity{sysBoundaryContainer.isPeriodic(0),
                                  sysBoundaryContainer.isPeriodic(1),
                                  sysBoundaryContainer.isPeriodic(2)};

   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> perBGrid(fsGridDimensions, MPI_COMM_WORLD, periodicity,P::manualFsGridDecomposition);
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> perBDt2Grid(fsGridDimensions, MPI_COMM_WORLD, periodicity,P::manualFsGridDecomposition);
   FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> EGrid(fsGridDimensions, MPI_COMM_WORLD, periodicity,P::manualFsGridDecomposition);
   FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> EDt2Grid(fsGridDimensions, MPI_COMM_WORLD, periodicity,P::manualFsGridDecomposition);
   FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> EHallGrid(fsGridDimensions, MPI_COMM_WORLD, periodicity,P::manualFsGridDecomposition);
   FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> EGradPeGrid(fsGridDimensions, MPI_COMM_WORLD, periodicity,P::manualFsGridDecomposition);
   FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> EGradPeDt2Grid(fsGridDimensions, MPI_COMM_WORLD, periodicity,P::manualFsGridDecomposition);
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> momentsGrid(fsGridDimensions, MPI_COMM_WORLD, periodicity,P::manualFsGridDecomposition);
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> momentsDt2Grid(fsGridDimensions, MPI_COMM_WORLD, periodicity,P::manualFsGridDecomposition);
   FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> dPerBGrid(fsGridDimensions, MPI_COMM_WORLD, periodicity,P::manualFsGridDecomposition);
   FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> dMomentsGrid(fsGridDimensions, MPI_COMM_WORLD, periodicity,P::manualFsGridDecomposition);
   FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> dMomentsDt2Grid(fsGridDimensions, MPI_COMM_WORLD, periodicity,P::manualFsGridDecomposition);
   FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> BgBGrid(fsGridDimensions, MPI_COMM_WORLD, periodicity,P::manualFsGridDecomposition);
   FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> volGrid(fsGridDimensions, MPI_COMM_WORLD, periodicity,P::manualFsGridDecomposition);
   FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> technicalGrid(fsGridDimensions, MPI_COMM_WORLD, periodicity,P::manualFsGridDecomposition);

   // Set DX, DY and DZ
   // TODO: This is currently just taking the values from cell 1, and assuming them to be
   // constant throughout the simulation.
   perBGrid.DX = perBDt2Grid.DX = EGrid.DX = EDt2Grid.DX = EHallGrid.DX = EGradPeGrid.DX = EGradPeDt2Grid.DX = momentsGrid.DX
      = momentsDt2Grid.DX = dPerBGrid.DX = dMomentsGrid.DX = dMomentsDt2Grid.DX = BgBGrid.DX = volGrid.DX = technicalGrid.DX
      = P::dx_ini / pow(2, P::amrMaxSpatialRefLevel);
   perBGrid.DY = perBDt2Grid.DY = EGrid.DY = EDt2Grid.DY = EHallGrid.DY = EGradPeGrid.DY = EGradPeDt2Grid.DY = momentsGrid.DY
      = momentsDt2Grid.DY = dPerBGrid.DY = dMomentsGrid.DY = dMomentsDt2Grid.DY = BgBGrid.DY = volGrid.DY = technicalGrid.DY
      = P::dy_ini / pow(2, P::amrMaxSpatialRefLevel);
   perBGrid.DZ = perBDt2Grid.DZ = EGrid.DZ = EDt2Grid.DZ = EHallGrid.DZ = EGradPeGrid.DZ = EGradPeDt2Grid.DZ = momentsGrid.DZ
      = momentsDt2Grid.DZ = dPerBGrid.DZ = dMomentsGrid.DZ = dMomentsDt2Grid.DZ = BgBGrid.DZ = volGrid.DZ = technicalGrid.DZ
      = P::dz_ini / pow(2, P::amrMaxSpatialRefLevel);
   // Set the physical start (lower left corner) X, Y, Z
   perBGrid.physicalGlobalStart = perBDt2Grid.physicalGlobalStart = EGrid.physicalGlobalStart = EDt2Grid.physicalGlobalStart
      = EHallGrid.physicalGlobalStart = EGradPeGrid.physicalGlobalStart = EGradPeDt2Grid.physicalGlobalStart = momentsGrid.physicalGlobalStart
      = momentsDt2Grid.physicalGlobalStart = dPerBGrid.physicalGlobalStart = dMomentsGrid.physicalGlobalStart = dMomentsDt2Grid.physicalGlobalStart
      = BgBGrid.physicalGlobalStart = volGrid.physicalGlobalStart = technicalGrid.physicalGlobalStart
      = {P::xmin, P::ymin, P::zmin};

   // Checking that spatial cells are cubic, otherwise field solver is incorrect (cf. derivatives in E, Hall term)
   constexpr Real uniformTolerance=1e-3;
   if ((abs((technicalGrid.DX - technicalGrid.DY) / technicalGrid.DX) >uniformTolerance) ||
       (abs((technicalGrid.DX - technicalGrid.DZ) / technicalGrid.DX) >uniformTolerance) ||
       (abs((technicalGrid.DY - technicalGrid.DZ) / technicalGrid.DY) >uniformTolerance)) {
      if (myRank == MASTER_RANK) {
         std::cerr << "WARNING: Your spatial cells seem not to be cubic. The simulation will now abort!" << std::endl;
      }
      //just abort sending SIGTERM to all tasks
      MPI_Abort(MPI_COMM_WORLD, -1);
   }
   initFsTimer.stop();

   // Initialize grid.  After initializeGrid local cells have dist
   // functions, and B fields set. Cells have also been classified for
   // the various sys boundary conditions.  All remote cells have been
   // created. All spatial date computed this far is up to date for
   // FULL_NEIGHBORHOOD. Block lists up to date for
   // VLASOV_SOLVER_NEIGHBORHOOD (but dist function has not been communicated)
   phiprof::Timer initGridsTimer {"Init grids"};
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry> mpiGrid;

   initializeGrids(
      argn,
      args,
      mpiGrid,
      perBGrid,
      BgBGrid,
      momentsGrid,
      momentsDt2Grid,
      dMomentsGrid,
      EGrid,
      EGradPeGrid,
      volGrid,
      technicalGrid,
      sysBoundaryContainer,
      *project
   );
   const std::vector<CellID>& cells = getLocalCells();

   phiprof::Timer reportMemoryTimer {"report-memory-consumption"};
   if (myRank == MASTER_RANK){
      cout << "(MAIN): Completed grid initialization." << endl;
      logFile << "(MAIN): Completed grid initialization." << endl << writeVerbose;
   }
   report_memory_consumption(mpiGrid);
   reportMemoryTimer.stop();
   
   // There are projects that have non-uniform and non-zero perturbed B, e.g. Magnetosphere with dipole type 4.
   // For inflow cells (e.g. maxwellian), we cannot take a FSgrid perturbed B value from the templateCell,
   // because we need a copy of the value from initialization in both perBGrid and perBDt2Grid and it isn't
   // touched as we are in boundary cells for components that aren't solved. We do a straight full copy instead
   // of looping and detecting boundary types here.
   perBDt2Grid.copyData(perBGrid);

   initGridsTimer.stop();
   
   // Initialize data reduction operators. This should be done elsewhere in order to initialize 
   // user-defined operators:
   phiprof::Timer initDROsTimer {"Init DROs"};
   DataReducer outputReducer, diagnosticReducer;

   if(P::writeFullBGB) {
      // We need the following variables for this, let's just erase and replace the entries in the list
      P::outputVariableList.clear();
      P::outputVariableList= {"fg_b_background", "fg_b_background_vol", "fg_derivs_b_background"};
   }

   initializeDataReducers(&outputReducer, &diagnosticReducer);
   initDROsTimer.stop();

   // Free up memory:
   readparameters.~Readparameters();

   if(P::writeFullBGB) {
      logFile << "Writing out full BGB components and derivatives and exiting." << endl << writeVerbose;

      // initialize the communicators so we can write out ionosphere grid metadata.
      SBC::ionosphereGrid.updateIonosphereCommunicator(mpiGrid, technicalGrid);

      P::systemWriteDistributionWriteStride.push_back(0);
      P::systemWriteName.push_back("bgb");
      P::systemWriteDistributionWriteXlineStride.push_back(0);
      P::systemWriteDistributionWriteYlineStride.push_back(0);
      P::systemWriteDistributionWriteZlineStride.push_back(0);
      P::systemWritePath.push_back("./");
      P::systemWriteFsGrid.push_back(true);

      for(uint si=0; si<P::systemWriteName.size(); si++) {
         P::systemWrites.push_back(0);
      }

      const bool writeGhosts = true;
      if( writeGrid(mpiGrid,
            perBGrid,
            EGrid,
            EHallGrid,
            EGradPeGrid,
            momentsGrid,
            dPerBGrid,
            dMomentsGrid,
            BgBGrid,
            volGrid,
            technicalGrid,
            version,
            config,
            &outputReducer,
            P::systemWriteName.size()-1,
            P::restartStripeFactor,
            writeGhosts
         ) == false
      ) {
         cerr << "FAILED TO WRITE GRID AT " << __FILE__ << " " << __LINE__ << endl;
      }
      initTimer.stop();
      mainTimer.stop();
      
      phiprof::print(MPI_COMM_WORLD,"phiprof");
      
      if (myRank == MASTER_RANK) logFile << "(MAIN): Exiting." << endl << writeVerbose;
      logFile.close();
      if (P::diagnosticInterval != 0) diagnostic.close();
      
      perBGrid.finalize();
      perBDt2Grid.finalize();
      EGrid.finalize();
      EDt2Grid.finalize();
      EHallGrid.finalize();
      EGradPeGrid.finalize();
      EGradPeDt2Grid.finalize();
      momentsGrid.finalize();
      momentsDt2Grid.finalize();
      dPerBGrid.finalize();
      dMomentsGrid.finalize();
      dMomentsDt2Grid.finalize();
      BgBGrid.finalize();
      volGrid.finalize();
      technicalGrid.finalize();

      MPI_Finalize();
      return 0;
   }

   // Run the field solver once with zero dt. This will initialize
   // Fieldsolver dt limits, and also calculate volumetric B-fields.
   // At restart, all we need at this stage has been read from the restart, the rest will be recomputed in due time.
   if(P::isRestart == false) {
      propagateFields(
         perBGrid,
         perBDt2Grid,
         EGrid,
         EDt2Grid,
         EHallGrid,
         EGradPeGrid,
         EGradPeDt2Grid,
         momentsGrid,
         momentsDt2Grid,
         dPerBGrid,
         dMomentsGrid,
         dMomentsDt2Grid,
         BgBGrid,
         volGrid,
         technicalGrid,
         sysBoundaryContainer, 0.0, 1.0
      );
   }

   phiprof::Timer getFieldsTimer {"getFieldsFromFsGrid"};
   volGrid.updateGhostCells();
   getFieldsFromFsGrid(volGrid, BgBGrid, EGradPeGrid, dMomentsGrid, technicalGrid, mpiGrid, cells);
   getFieldsTimer.stop();

   // Build communicator for ionosphere solving
   SBC::ionosphereGrid.updateIonosphereCommunicator(mpiGrid, technicalGrid);
   // If not a restart, perBGrid and dPerBGrid are up to date after propagateFields just above. Otherwise, we should compute them.
   if(P::isRestart) {
      calculateDerivativesSimple(
         perBGrid,
         perBDt2Grid,
         momentsGrid,
         momentsDt2Grid,
         dPerBGrid,
         dMomentsGrid,
         dMomentsDt2Grid,
         technicalGrid,
         sysBoundaryContainer,
         RK_ORDER1, // Update and compute on non-dt2 grids.
         false // Don't communicate moments, they are not needed here.
      );
      dPerBGrid.updateGhostCells();
   }
   FieldTracing::calculateIonosphereFsgridCoupling(technicalGrid, perBGrid, dPerBGrid, SBC::ionosphereGrid.nodes, SBC::Ionosphere::radius);
   SBC::ionosphereGrid.initSolver(!P::isRestart); // If it is a restart we do not want to zero out everything
   if(SBC::Ionosphere::couplingInterval > 0 && P::isRestart) {
      SBC::Ionosphere::solveCount = floor(P::t / SBC::Ionosphere::couplingInterval)+1;
   } else {
      SBC::Ionosphere::solveCount = 1;
   }
   
   if(P::isRestart) {
      // If it is a restart, we want to regenerate proper ig_inplanecurrent as well in case there's IO before the next solver step.
      SBC::ionosphereGrid.calculateConductivityTensor(SBC::Ionosphere::F10_7, SBC::Ionosphere::recombAlpha, SBC::Ionosphere::backgroundIonisation, true);
   }

   phiprof::Timer dttimer {"compute-dt"};
   // Run Vlasov solver once with zero dt to initialize
   // per-cell dt limits. Also compute initial _R and _V moments at restart.
   calculateSpatialTranslation(mpiGrid,0.0,true);
   calculateAcceleration(mpiGrid,0.0);

   sysBoundaryContainer.setupL2OutflowAtRestart(mpiGrid);

   dttimer.stop();

   // Save restart data
   if (P::writeInitialState) {
      // Calculate these so refinement parameters can be tuned based on the vlsv
      calculateScaledDeltasSimple(mpiGrid);

      FieldTracing::reduceData(technicalGrid, perBGrid, dPerBGrid, mpiGrid, SBC::ionosphereGrid.nodes); /*!< Call the reductions (e.g. field tracing) */
      
      phiprof::Timer timer {"write-initial-state"};
      
      if (myRank == MASTER_RANK)
         logFile << "(IO): Writing initial state to disk, tstep = "  << endl << writeVerbose;
      P::systemWriteDistributionWriteStride.push_back(1);
      P::systemWriteName.push_back("initial-grid");
      P::systemWriteDistributionWriteXlineStride.push_back(0);
      P::systemWriteDistributionWriteYlineStride.push_back(0);
      P::systemWriteDistributionWriteZlineStride.push_back(0);
      P::systemWritePath.push_back("./");
      P::systemWriteFsGrid.push_back(true);

      for(uint si=0; si<P::systemWriteName.size(); si++) {
         P::systemWrites.push_back(0);
      }

      const bool writeGhosts = true;
      if( writeGrid(mpiGrid,
            perBGrid, // TODO: Merge all the fsgrids passed here into one meta-object
            EGrid,
            EHallGrid,
            EGradPeGrid,
            momentsGrid,
            dPerBGrid,
            dMomentsGrid,
            BgBGrid,
            volGrid,
            technicalGrid,
            version,
            config,
            &outputReducer,
            P::systemWriteName.size()-1,
            P::restartStripeFactor,
            writeGhosts
         ) == false
      ) {
         cerr << "FAILED TO WRITE GRID AT " << __FILE__ << " " << __LINE__ << endl;
      }

      P::systemWriteDistributionWriteStride.pop_back();
      P::systemWriteName.pop_back();
      P::systemWriteDistributionWriteXlineStride.pop_back();
      P::systemWriteDistributionWriteYlineStride.pop_back();
      P::systemWriteDistributionWriteZlineStride.pop_back();
      P::systemWritePath.pop_back();
      P::systemWriteFsGrid.pop_back();
   }

   //std::cerr << "calculating timesteps" << std::endl;

   P::tc_leapfrog_init = false; // not used
   if (P::isRestart == false) {
      //compute new dt
      phiprof::Timer computeDtimer {"compute-dt"};

      auto timeStepVector = computeNewTimeStep(mpiGrid, technicalGrid, dtMaxLocal, dtMaxGlobal, dtMinMaxLocal, dtMinMaxGlobal);


      calculateGlobalTcVariables(timeStepVector.at(1), timeStepVector.at(2));

      // this is called, because the next function checks against the smallest tcdt
      updateTimeclassDts(timeStepVector.at(1));

      // checks if smallest tcdt is good
      handleChangingofDt(dtMaxGlobal, dtIsChanged, newDt);

      if (P::dynamicTimestep == true && dtIsChanged) {
         // Only actually update the timestep if dynamicTimestep is on
         updateTimeclassDts(newDt);
         P::dt=P::timeclassDt[P::currentMaxTimeclass];
      } else if (P::dynamicTimestep == true && !dtIsChanged) {
         //updateTimeclassDts(timeStepVector.at(1));
         P::dt=P::timeclassDt[P::currentMaxTimeclass];
      } else {
         dtIsChanged = false;
         updateTimeclassDts(P::dt);
      }

      initiateAllCellTimeclasses(mpiGrid);

      if(myRank == MASTER_RANK){
         //std::cout << "timeclass dts = ";
         for(int i = 0; i <= P::currentMaxTimeclass; ++i){
            //std::cout << i <<": "<<P::timeclassDt[i] << "s, ";
         }
         //std::cout << endl;
      }

      for (vector<CellID>::const_iterator cell_id=cells.begin(); cell_id!=cells.end(); ++cell_id) {

         SpatialCell* cell = mpiGrid[*cell_id];
         // std::cerr << "timeclass and tcdt of cell " << cell->get_cellid() << ": " << cell->parameters[CellParams::TIMECLASS] << ", " << cell->parameters[CellParams::TIMECLASSDT] << std::endl;
      }
      //std::cerr << __FILE__ << " " << __LINE__ << std::endl;
      computeDtimer.stop();

      balanceLoad(mpiGrid, sysBoundaryContainer, technicalGrid);
      
      //std::cerr << __FILE__ << " " << __LINE__ << std::endl;

      auto lCells = getLocalCells();
      //std::cerr << "checking cell ghost distributions:\n";
      for (size_t c=0; c<lCells.size(); ++c) {
         const CellID cell = lCells[c];
         SpatialCell* spatialCell = mpiGrid[cell];
         //std::cerr << "cell " << cell << " has timeclass " << spatialCell->parameters[CellParams::TIMECLASS] << " and ghost distribution: ";
         for (auto& ghost : spatialCell->requested_timeclass_ghosts) {
            //std::cerr << ghost << " ";
         }
         //std::cerr << "\n";

      }

      //go forward by dt/2 in V, initializes leapfrog split. In restarts the
      //the distribution function is already propagated forward in time by dt/2
      phiprof::Timer propagateHalfTimer {"propagate-velocity-space-dt/2"};
      if (P::propagateVlasovAcceleration) {
         calculateAcceleration(mpiGrid, 0.5);
      } else {
         //zero step to set up moments _v
         calculateAcceleration(mpiGrid, 0.0);
      }
      P::tc_leapfrog_init = true;

      //std::cerr << __FILE__ << " " << __LINE__ << std::endl;

      propagateHalfTimer.stop();

      updatePreviousVMoments(mpiGrid, true);
      // std::cerr <<__FILE__<<":"<<__LINE__<<" ("<<myRank <<") Calling balanceLoad\n";

      // Apply boundary conditions
      if (P::propagateVlasovTranslation || P::propagateVlasovAcceleration ) {
         phiprof::Timer updateBoundariesTimer {("update system boundaries (Vlasov post-acceleration)")};
         sysBoundaryContainer.applySysBoundaryVlasovConditions(mpiGrid, 0.5*P::dt, true);
         updateBoundariesTimer.stop();
         addTimedBarrier("barrier-boundary-conditions");
      }
      // std::cerr <<__FILE__<<":"<<__LINE__<<" ("<<myRank <<")\n";
      // Also update all moments. They won't be transmitted to FSgrid until the field solver is called, though.
      phiprof::Timer computeMomentsTimer {"Compute interp moments"};
      //std::cout << "for initial interpolated moments\n";
      calculateInterpolatedVelocityMoments(
         mpiGrid,
         CellParams::RHOM,
         CellParams::VX,
         CellParams::VY,
         CellParams::VZ,
         CellParams::RHOQ,
         CellParams::P_11,
         CellParams::P_22,
         CellParams::P_33,
         CellParams::P_23,
         CellParams::P_13,
         CellParams::P_12
      );
      
      updateParticlePopulations(mpiGrid);

      computeMomentsTimer.stop();
   }
// std::cerr <<__FILE__<<":"<<__LINE__<<" ("<<myRank <<")\n";
   initTimer.stop();

   // ***********************************
   // ***** INITIALIZATION COMPLETE *****
   // ***********************************

   // Main simulation loop:
   if (myRank == MASTER_RANK){
      cout << "(MAIN): Starting main simulation loop." << endl;
      logFile << "(MAIN): Starting main simulation loop." << endl << writeVerbose;
      //report filtering if we are in an AMR run
      if (P::amrMaxSpatialRefLevel>0){
         logFile<<"Filtering Report: "<<endl;
         for (int refLevel=0 ; refLevel<= P::amrMaxSpatialRefLevel; refLevel++){
            logFile<<"\tRefinement Level " <<refLevel<<"==> Passes "<<P::numPasses.at(refLevel)<<endl;
         }
            logFile<<endl;
      }
   }

   phiprof::Timer reportMemTimer {"report-memory-consumption"};
   report_memory_consumption(mpiGrid);
   reportMemTimer.stop();
   
   uint64_t computedCells=0;
   //Compute here based on time what the file intervals are
   P::systemWrites.clear();
   for(uint i=0;i< P::systemWriteTimeInterval.size();i++){
      int index=(int)(P::t_min/P::systemWriteTimeInterval[i]);
      //if we are already over 1% further than the time interval time that
      //is requested for writing, then jump to next writing index. This is to
      //make sure that at restart we do not write in the middle of
      //the interval.
      if(P::t_min>(index+0.01)*P::systemWriteTimeInterval[i]) {
         index++;
         // Special case for large timesteps
         int index2=(int)((P::t_min+P::dt)/P::systemWriteTimeInterval[i]);
         if (index2>index) index=index2;
      }
      P::systemWrites.push_back(index);
   }

   // Invalidate cached cell lists just to be sure (might not be needed)
   P::meshRepartitioned = true;

   uint wallTimeRestartCounter=1;
   uint recoverCounter=0;

   int doNow[donow::N_DONOW] = {0}; // 0: writeRestartNow, 1: writeRecoverNow, 2: balanceLoadNow, 3: refineNow ; declared outside main loop
   bool overrideRebalanceNow = false; // declared outside main loop
   bool refineNow = false; // declared outside main loop

   addTimedBarrier("barrier-end-initialization");
   
   phiprof::Timer simulationTimer {"Simulation"};
   double startTime=  MPI_Wtime();
   double beforeTime = MPI_Wtime();
   double beforeSimulationTime=P::t_min;
   double beforeStep=P::tstep_min;

   while(P::tstep <= P::tstep_max  &&
         P::t-P::dt <= P::t_max+DT_EPSILON &&
         wallTimeRestartCounter <= P::exitAfterRestarts) {

      //std::cout << "start of main simulation loop, below dt, timeclassDts, currentmaxtimeclass" << std::endl;
      //std::cout << P::dt << std::endl;
      for (auto i: P::timeclassDt) {
         //std::cout << i << " ";
      }
      //std::cout << endl;
      //std::cout << P::currentMaxTimeclass << std::endl;
      
      addTimedBarrier("barrier-loop-start");
      
      phiprof::Timer ioTimer {"IO"};

      phiprof::Timer externalsTimer {"checkExternalCommands"};
      if(myRank ==  MASTER_RANK) {
         // check whether STOP or KILL or SAVE has been passed, should be done by MASTER_RANK only as it can reset P::bailout_write_restart
         checkExternalCommands();
      }
      externalsTimer.stop();

      //write out phiprof profiles and logs with a lower interval than normal
      //diagnostic (every 10 diagnostic intervals).
      phiprof::Timer loggingTimer {"logfile-io"};
      logFile << "---------- tstep = " << P::tstep << " (" <<P::fractionalTimestep<<"/"<<(2 << (P::currentMaxTimeclass-1)) <<") t = " << P::t <<" dt = " << P::dt << " FS cycles = " << P::fieldSolverSubcycles << " ----------" << endl;
      if (P::diagnosticInterval != 0 &&
          P::tstep % (P::diagnosticInterval*10) == 0 &&
          P::tstep-P::tstep_min >0) {

         phiprof::print(MPI_COMM_WORLD,"phiprof");

         double currentTime=MPI_Wtime();
         double timePerStep=double(currentTime  - beforeTime) / (P::tstep-beforeStep);
         double timePerSecond=double(currentTime  - beforeTime) / (P::t-beforeSimulationTime + DT_EPSILON);
         double remainingTime=min(timePerStep*(P::tstep_max-P::tstep),timePerSecond*(P::t_max-P::t));
         time_t finalWallTime=time(NULL)+(time_t)remainingTime; //assume time_t is in seconds, as it is almost always
         struct tm *finalWallTimeInfo=localtime(&finalWallTime);
         logFile << "(TIME) current " << nNodes*(currentTime - startTime)/3600 << " node-hours" << endl;
         #if _OPENMP
            logFile << "(TIME) current " << omp_get_max_threads()*mpiProcs*(currentTime - startTime)/3600 << " thread-hours" << endl;
         #endif
         logFile << "(TIME) current walltime/step " << timePerStep<< " s" <<endl;
         logFile << "(TIME) current walltime/simusecond " << timePerSecond<<" s" <<endl;
         logFile << "(TIME) Estimated completion time is " <<asctime(finalWallTimeInfo)<<endl;
         //reset before values, we want to report speed since last report of speed.
         beforeTime = MPI_Wtime();
         beforeSimulationTime=P::t;
         beforeStep=P::tstep;

      }
      logFile << writeVerbose;
      loggingTimer.stop();

      // Check whether diagnostic output has to be produced
      if (P::diagnosticInterval != 0 && (P::tstep % P::diagnosticInterval == 0) && P::fractionalTimestep == 0) {
         phiprof::Timer memTimer {"memory-report"};
         memTimer.start();
         report_memory_consumption(mpiGrid);
         memTimer.stop();
         phiprof::Timer cellTimer {"cell-count-report"};
         cellTimer.start();
         report_cell_and_block_counts(mpiGrid);
         cellTimer.stop();

         phiprof::Timer diagnosticTimer {"diagnostic-io"};
         if (writeDiagnostic(mpiGrid, diagnosticReducer) == false) {
            if(myRank == MASTER_RANK)  cerr << "ERROR with diagnostic computation" << endl;

         }
      }

      // write system, loop through write classes
      for (uint i = 0; i < P::systemWriteTimeInterval.size(); i++) {
         // if (true || (P::systemWriteTimeInterval[i] >= 0.0 &&
         //     P::t >= P::systemWrites[i] * P::systemWriteTimeInterval[i] - DT_EPSILON)) {
            if (P::systemWriteTimeInterval[i] >= 0.0 &&
                P::t >= P::systemWrites[i] * P::systemWriteTimeInterval[i] - DT_EPSILON) {
            // If we have only just restarted, the bulk file should already exist from the previous slot.
            if ((P::tstep == P::tstep_min) && (P::tstep>0)) {
               P::systemWrites[i]++;
               // Special case for large timesteps
               int index2=(int)((P::t+P::dt)/P::systemWriteTimeInterval[i]);
               if (index2>P::systemWrites[i]) P::systemWrites[i]=index2;
               continue;
            }

            // Calculate these so refinement parameters can be tuned based on the vlsv
            calculateScaledDeltasSimple(mpiGrid);

            FieldTracing::reduceData(technicalGrid, perBGrid, dPerBGrid, mpiGrid, SBC::ionosphereGrid.nodes); /*!< Call the reductions (e.g. field tracing) */
            
            phiprof::Timer writeSysTimer {"write-system"};
            logFile << "(IO): Writing spatial cell and reduced system data to disk, tstep = " << P::tstep << " t = " << P::t << endl << writeVerbose;
            const bool writeGhosts = true;
            if(writeGrid(mpiGrid,
               perBGrid, // TODO: Merge all the fsgrids passed here into one meta-object
               EGrid,
               EHallGrid,
               EGradPeGrid,
               momentsGrid,
               dPerBGrid,
               dMomentsGrid,
               BgBGrid,
               volGrid,
               technicalGrid,
               version,
               config,
               &outputReducer,
               i,
               P::systemStripeFactor,
               writeGhosts
               ) == false
            ) {
               cerr << "FAILED TO WRITE GRID AT" << __FILE__ << " " << __LINE__ << endl;
            }
            P::systemWrites[i]++;
            // Special case for large timesteps
            int index2=(int)((P::t+P::dt)/P::systemWriteTimeInterval[i]);
            if (index2>P::systemWrites[i]) P::systemWrites[i]=index2;
            logFile << "(IO): .... done!" << endl << writeVerbose;
         }
      }

      // Reduce globalflags::bailingOut from all processes
      phiprof::Timer bailoutReduceTimer {"Bailout-allreduce"};
      MPI_Allreduce(&(globalflags::bailingOut), &(doBailout), 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      bailoutReduceTimer.stop();

      // Write restart data if needed
      // Combined with checking of additional load balancing to have only one collective call.
      phiprof::Timer restartCheckTimer {"compute-is-restart-written-and-extra-LB"};
      if (myRank == MASTER_RANK) {
         doNow[donow::SAVE] = 0;
         doNow[donow::DORC] = 0;
         if (  (P::saveRestartWalltimeInterval >= 0.0
            && (P::saveRestartWalltimeInterval*wallTimeRestartCounter <=  MPI_Wtime()-initialWtime
               || (P::tstep == P::tstep_max && P::fractionalTimestep == 0)
               || (P::t >= P::t_max && P::fractionalTimestep == 0)))
            || (doBailout > 0 && P::bailout_write_restart)
            || globalflags::writeRestart
         ) {
            doNow[donow::SAVE] = 1;
            if (globalflags::writeRestart == true) {
               doNow[donow::SAVE] = 2; // Setting to 2 so as to not increment the restart count below.
               globalflags::writeRestart = false; // This flag is only used by MASTER_RANK here and it needs to be reset after a restart write has been issued.
            }
         }
         if (  (P::saveRecoverTstepInterval > 0
            && P::tstep % P::saveRecoverTstepInterval == 0
            && P::tstep != P::tstep_min)
            || globalflags::writeRecover
         ) {
            doNow[donow::DORC] = 1;
            if (globalflags::writeRecover == true) {
               globalflags::writeRecover = false; // This flag is only used by MASTER_RANK here and it needs to be reset after a recover write has been issued.
            }
         }
         if (globalflags::balanceLoad || globalflags::doRefine) {
            doNow[donow::DOLB] = 1;
            globalflags::balanceLoad = false;
            if (globalflags::doRefine) {
               doNow[donow::DOMR] = 1;
               globalflags::doRefine = false;
            }
         }
      }
      MPI_Bcast( &doNow, 4 , MPI_INT , MASTER_RANK ,MPI_COMM_WORLD);
      if (doNow[donow::DOLB] == 1) {
         P::prepareForRebalance = true;
      }
      if (doNow[donow::DOMR] == 1) {
         refineNow = true;
      }
      restartCheckTimer.stop();

      if (doNow[donow::SAVE] >= 1){ // write restart
         phiprof::Timer timer {"write-restart"};
         if (doNow[donow::SAVE] == 1) { // write restart
            wallTimeRestartCounter++;
         }

         // Refinement params for restart refinement
         calculateScaledDeltasSimple(mpiGrid);
         
         if (myRank == MASTER_RANK)
            logFile << "(IO): Writing restart data to disk, tstep = " << P::tstep << " t = " << P::t << endl << writeVerbose;
         //Write the restart:
         if( writeRestart(mpiGrid,
                  perBGrid, // TODO: Merge all the fsgrids passed here into one meta-object
                  EGrid,
                  EHallGrid,
                  EGradPeGrid,
                  momentsGrid,
                  dPerBGrid,
                  dMomentsGrid,
                  BgBGrid,
                  volGrid,
                  technicalGrid,
                  version,
                  config,
                  outputReducer,
                  "restart",
                  (uint)P::t,
                  true, // add the date of the file to the name
                  P::restartStripeFactor) == false ) {
            logFile << "(IO): ERROR Failed to write restart!" << endl << writeVerbose;
            cerr << "FAILED TO WRITE RESTART" << endl;
         }
         if (myRank == MASTER_RANK) {
            logFile << "(IO): .... done!"<< endl << writeVerbose;
         }
         timer.stop();
      }
      if (doNow[donow::DORC] == 1){ // write recover
         phiprof::Timer timer {"write-recover"};

         // Refinement params for restart refinement
         calculateScaledDeltasSimple(mpiGrid);

         if (myRank == MASTER_RANK)
            logFile << "(IO): Writing recover data to disk, index = " << recoverCounter % P::recoverMaxFiles << ", tstep = " << P::tstep << " t = " << P::t << endl << writeVerbose;
         //Write the recover:
         if( writeRestart(mpiGrid,
                  perBGrid, // TODO: Merge all the fsgrids passed here into one meta-object
                  EGrid,
                  EHallGrid,
                  EGradPeGrid,
                  momentsGrid,
                  dPerBGrid,
                  dMomentsGrid,
                  BgBGrid,
                  volGrid,
                  technicalGrid,
                  version,
                  config,
                  outputReducer,
                  "recover",
                  recoverCounter % P::recoverMaxFiles,
                  false, // overwrite so do not put date in file name
                  P::restartStripeFactor) == false ) {
            logFile << "(IO): ERROR Failed to write recover!" << endl << writeVerbose;
            cerr << "FAILED TO WRITE RECOVER" << endl;
         }
         recoverCounter++;
         if (myRank == MASTER_RANK) {
            logFile << "(IO): .... done!"<< endl << writeVerbose;
         }
         timer.stop();
      }

      ioTimer.stop();
      addTimedBarrier("barrier-end-io");

      // reset these for next time around
      doNow[donow::SAVE] = doNow[donow::DORC] = doNow[donow::DOLB] = doNow[donow::DOMR] = 0;

      //no need to propagate if we are on the final step, we just
      //wanted to make sure all IO is done even for final step
      if(P::tstep == P::tstep_max ||
         P::t >= P::t_max ||
         doBailout > 0) {
         break;
      }

      // std::cout << "main loop at" << __FILE__ << " " << __LINE__ << " " << P::tstep << " " << P::fractionalTimestep << std::endl;

      //Re-loadbalance if needed
      //TODO - add LB measure and do LB if it exceeds threshold
      if(((P::tstep % P::rebalanceInterval == 0 && P::tstep > P::tstep_min && P::fractionalTimestep == 0) || overrideRebalanceNow)) {
         logFile << "(LB): Start load balance, tstep = " << P::tstep << " t = " << P::t << endl << writeVerbose;

         phiprof::Timer shrinkTimer {"Shrink_to_fit"};
         // * shrink to fit before LB * //
         shrink_to_fit_grid_data(mpiGrid);
         shrinkTimer.stop();

         if (refineNow || (!dtIsChanged && P::adaptRefinement && P::tstep % (P::rebalanceInterval * P::refineCadence) == 0 && P::t > P::refineAfter)) { 
            logFile << "(AMR): Adapting refinement!"  << endl << writeVerbose;
            refineNow = false;
            if (!adaptRefinement(mpiGrid, technicalGrid, sysBoundaryContainer, *project)) {
               // OOM, rebalance and try again
               logFile << "(LB) AMR rebalancing with heavier refinement weights." << endl;
               globalflags::bailingOut = false; // Reset this
               for (auto id : mpiGrid.get_local_cells_to_refine()) {
                  mpiGrid[id]->parameters[CellParams::LBWEIGHTCOUNTER] *= 8.0;
               }
               balanceLoad(mpiGrid, sysBoundaryContainer, technicalGrid);
               // We can /= 8.0 now as cells have potentially migrated. Go back to block-based count for now.
               for (auto id : mpiGrid.get_local_cells_to_refine()) {
                  mpiGrid[id]->parameters[CellParams::LBWEIGHTCOUNTER] = 0;
                  for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
                     mpiGrid[id]->parameters[CellParams::LBWEIGHTCOUNTER] += mpiGrid[id]->get_number_of_velocity_blocks(popID);
                  }
               }

               mpiGrid.cancel_refining();
               if (!adaptRefinement(mpiGrid, technicalGrid, sysBoundaryContainer, *project)) {
                  for (auto id : mpiGrid.get_local_cells_to_refine()) {
                     mpiGrid[id]->parameters[CellParams::LBWEIGHTCOUNTER] *= 8.0;
                  }
                  continue;   // Refinement failed and we're bailing out
               } else {
                  globalflags::bailingOut = false; // Reset this
               }
            }

            // Calculate new dt limits since we might break CFL when refining
            phiprof::Timer computeDtimer {"compute-dt-amr"};
            calculateSpatialTranslation(mpiGrid,0.0,true);
            calculateAcceleration(mpiGrid,0.0);
         }
         // This now uses the block-based count just copied between the two refinement calls above.
         balanceLoad(mpiGrid, sysBoundaryContainer, technicalGrid);
         addTimedBarrier("barrier-end-load-balance");
         logFile << "(LB): ... done!"  << endl << writeVerbose;
         P::prepareForRebalance = false;

         overrideRebalanceNow = false;

         // Make sure the ionosphere communicator is up-to-date, in case inner boundary cells
         // moved.
         SBC::ionosphereGrid.updateIonosphereCommunicator(mpiGrid, technicalGrid);
      }
   
      //get local cells
      const vector<CellID>& cells = getLocalCells();

      //compute how many spatial cells we solve for this step
      computedCells=0;
      for(size_t i=0; i<cells.size(); i++) {
         for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID)
            computedCells += (uint64_t)mpiGrid[cells[i]]->get_number_of_velocity_blocks(popID)*WID3;
      }

      //Check if dt needs to be changed, and propagate V back a half-step to change dt and set up new situation
      //do not compute new dt on first step (in restarts dt comes from file, otherwise it was initialized before we entered
      //simulation loop
      // FIXME what if dt changes at a restart?? ---- dunno about this, but...
      //
      // TODO timeclasses: We need to check for a new timestep at each smallest timestep - might very well
      // change, 4 levels -> smallest does 16 steps to one long step
      // This redefines TCs actually everywhere
      // WHAT IF THIS MEANS SLOWER REGION CHANGES TC, instead of having just an adjusted timestep?
      // -make sure TC is immutable during slowest step - this is a no-go, since the fastest steps can explode in between
      // -calculate the back-step+new-half-step correction only when about to propagate slower cells
      //   -> no need to redo acceleration each time fastest class changes dt?
      //      ... but updated V moments needed?
      //   -> do the acc shuffle for all cells to begin with

      // std::vector<CellID> cellsToTestIncrease = {25};

      // if (P::tstep == 20 && P::fractionalTimestep == 0) {
      //    std::cout << "TESTING INCREASE TIMECLASS" << std::endl;
      //    // we must roll back all cells that get dropped a timeclass.

      //    std::cout << "rolling back all cells that get raised a timeclass" << std::endl;
      //    calculateAcceleration(mpiGrid, -0.5, true, cellsToTestIncrease); // This sets cells back to previous TIME_R

      //    increaseTimeclass(mpiGrid, cellsToTestIncrease, additionalTimeclassCreated);

      //    calculateAcceleration(mpiGrid, 0.5, true, cellsToTestIncrease); // This propagates by 0.5

      // }

      // for (auto c: cells) {

         // std::cout << "\n";
         // std::cout << "Cell " << c << " has timeclass " << mpiGrid[c]->parameters[CellParams::TIMECLASS] << std::endl;
         // std::cout << "sysboundaryflag: " << mpiGrid[c]->sysBoundaryFlag << std::endl;
         // std::cout << "regular vspace size: " << mpiGrid[c]->get_velocity_mesh(0)->size() << std::endl;
         // std::cout << "its requested timeclass ghosts: " << std::endl;
         // for (int req: mpiGrid[c]->requested_timeclass_ghosts) {
         //    std::cout << req << " " << std::endl;
         //    std::cout << "vspace size of that requested: " << mpiGrid[c]->get_velocity_mesh(0, req)->size() << std::endl;
         // }

         // std::cout << "its requested timeclass copy ghosts: " << std::endl;
         // for (int req: mpiGrid[c]->requested_timeclass_copy_ghosts) {
         //    std::cout << req << " " << std::endl;
         // }

      //    std::cout << "maxdt in v and r: " << mpiGrid[c]->parameters[CellParams::MAXVDT] << ", " << mpiGrid[c]->parameters[CellParams::MAXRDT] << std::endl;
      //    std::cout << "is not sysboundary cell: " << (mpiGrid[c]->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) << std::endl;
      //    std::cout << "timeclassExactHaloExtent: " << P::timeclassExactHaloExtent << std::endl;
      //    std::cout << "size of timeghost exact neighborhood:" << mpiGrid.get_neighbors_of(c, VLASOV_SOLVER_TIMEGHOST_EXACT_HALO_NEIGHBORHOOD_ID)->size() << std::endl;


      //    std::cout << "\n";
      // }

      // std::vector<Real> newTimeclassDts = std::vector<Real>(P::initialMaxTimeclass+1);

      // if (P::initialMaxTimeclass >= 4 && P::fractionalTimestep == 0 && P::tstep != 0 && P::dynamicTimestep) {
      //    logFile << "(DT): Checking dynamic timestep on fractimestep0, tstep = " << P::tstep << " t = " << P::t << endl << writeVerbose;
      //    logFile << "fractimestep is 0, recomputing timeclasses and dts" << endl << writeVerbose;
      //    //computeNewTimeStep(mpiGrid, technicalGrid, dtIsChanged);
         
      //    if (P::propagateVlasovAcceleration) {
      //       // Back half dt to real time, forward by new half dt
      //       calculateAcceleration(mpiGrid,-0.5); //This sets cells back to previous TIME_R
      //       calculateAcceleration(mpiGrid, 0.5); //This propagates by 0.5 
      //    } else {
      //       //zero step to set up moments _v
      //       calculateAcceleration(mpiGrid, 0.0);
      //    }

      // P::dt=newDt;
      // P::timeclassDt = newTimeclassDts;
      
      // logFile <<" dt changed to "<<P::dt <<"s, distribution function was half-stepped to real-time and back"<<endl<<writeVerbose;
      // goto endOfDynDTCheck; // All dt's and timeclasses have been recomputed, so we can continue with the next step.

      // }

      /*
      steps:
      1. Check if any cell is on the wrong timeclass, i.e. compare its cell--wise dt with its timeclass timestep
      2. If any cells are on any timeclass, do one of two things
         2.a. if fractionaltimestep == 0, recalculate and assign all cells freshly
         2.b. if not, increase the timeclass of all bad cells by one
            2.b.I. if individual cell timeclasses are increased, we must then on the next fractimestep freshly assign all cells
                   !! it might be worth considering doing this fresh reassignemt to all cells in production runs.
      */


      if (P::dynamicTimestep) {

         dtIsChanged = false;
         // this calls tooLarge and tooSmall both for the smallest tc timestep
         handleChangingofDt(dtMaxGlobal, dtIsChanged, newDt);
         // checks if dt is good

         //std::cout << "doing dynamic timestep checks" << std::endl;
         auto timestepvector = computeNewTimeStep(mpiGrid, technicalGrid, dtMaxLocal, dtMaxGlobal, dtMinMaxLocal, dtMinMaxGlobal);
         //std::cout << "global minimun timestep is " << timestepvector.at(1) << std::endl;
         //std::cout << "ratio of global minimum and smallest tcdt: " << timestepvector.at(1) / P::timeclassDt[P::currentMaxTimeclass - P::timeclassBuffer] << std::endl;

         std::vector<Real> placeholder1(3), placeholder2(3);

         // update maxrdt
         reduce_vlasov_dt(mpiGrid, cells, placeholder1, placeholder2);

         // update maxvdt
         // this is done when calculateAcceleration is called and is redundant here, but for testing
         for (CellID c: cells) {
            for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
               updateAccelerationMaxdt(mpiGrid[c], popID);
               mpiGrid[c]->parameters[CellParams::MAXVDT] = min(mpiGrid[c]->parameters[CellParams::MAXVDT], mpiGrid[c]->get_max_v_dt(popID));
            }
         }

         std::vector<CellID> badTcCells = checkCellTimeclasses(mpiGrid);

         if (badTcCells.size() != 0) {

            //std::cout << badTcCells.size() << " bad cells found:" << std::endl;
            logFile << "\nBAD CELLS FOUND, FRACTIMESTEP = " << P::fractionalTimestep << "\n" << std::endl;

            for (auto c: badTcCells) {
               //std::cout << mpiGrid[c]->get_cellid() << std::endl;
            }

            if (P::staticTimeclasses) {
               std::cout << "aborting due to timeclass changing with static timeclasses" << std::endl;
               abort();
            }

            if (P::fractionalTimestep == 0) {

               //std::cout << "recomputing all timeclasses as fractimestep is 0" << std::endl;
               logFile << "\nrecomputing all timeclasses as fractimestep is 0\n" << std::endl;

               // step back all cells in acceleration space
               calculateAcceleration(mpiGrid, -0.5);

               // recompute all timeclasses again, and step them forward in space
               
               calculateGlobalTcVariables(timestepvector.at(1), timestepvector.at(2));

               // this is called, because the next function checks against the smallest tcdt
               updateTimeclassDts(timestepvector.at(1));

               // checks if smallest tcdt is good
               handleChangingofDt(dtMaxGlobal, dtIsChanged, newDt);

               if (dtIsChanged) {
                  // the above might change dt
                  updateTimeclassDts(newDt);
                  P::dt=P::timeclassDt[P::currentMaxTimeclass];
               } else {
                  updateTimeclassDts(timestepvector.at(1));
                  P::dt=P::timeclassDt[P::currentMaxTimeclass];
               }

               assingCellTimeclasses(mpiGrid);         

               if(myRank == MASTER_RANK){
                  //std::cout << "timeclass dts = ";
                  for(int i = 0; i <= P::currentMaxTimeclass; ++i){
                     //std::cout << i <<": "<<P::timeclassDt[i] << "s, ";
                  }
                  //std::cout << endl;
               }

               // restep cells
               calculateAcceleration(mpiGrid, 0.5);
               logFile << "\ncompleted timeclass reassignemtn\n" << std::endl;

            } else {
               aCellHadTimeclassChanged = true;
               //std::cout << "increasing timeclasses of bad cells" << std::endl;
               logFile << "\nincreasing timeclasses of bad cells\n" << std::endl;

               //rolling back
               calculateAcceleration(mpiGrid, -0.5, true, badTcCells);

               increaseTimeclass(mpiGrid, badTcCells, additionalTimeclassCreated);

               // rolling forward

               calculateAcceleration(mpiGrid, 0.5, true, badTcCells);
               logFile << "\ncompleted timeclass increase\n" << std::endl;

            }

         } else {
            //std::cout << "all cells pass check" << std::endl;
         }

      }

      if (aCellHadTimeclassChanged == true && P::fractionalTimestep == 0) {


         auto timestepvector = computeNewTimeStep(mpiGrid, technicalGrid, dtMaxLocal, dtMaxGlobal, dtMinMaxLocal, dtMinMaxGlobal);

         //if we upgraded a cell timeclass, we do some stuff on the next 0th fractional timestep

         aCellHadTimeclassChanged = false;
         //std::cout << "recomputing all timeclasses as fractimestep is 0, and we increased a cell timeclass on the previous full timestep." << std::endl;

         // step back all cells in acceleration space
         calculateAcceleration(mpiGrid, -0.5);

         // recompute all timeclasses again, and step them forward in space
         
         calculateGlobalTcVariables(timestepvector.at(1), timestepvector.at(2));

         // this is called, because the next function checks against the smallest tcdt
         updateTimeclassDts(timestepvector.at(1));

         // checks if smallest tcdt is good
         handleChangingofDt(dtMaxGlobal, dtIsChanged, newDt);

         if (dtIsChanged) {
            // the above might change dt
            updateTimeclassDts(newDt);
            P::dt=P::timeclassDt[P::currentMaxTimeclass];
         } else {
            updateTimeclassDts(timestepvector.at(1));
            P::dt=P::timeclassDt[P::currentMaxTimeclass];
         }

         assingCellTimeclasses(mpiGrid);         

         if(myRank == MASTER_RANK){
            //std::cout << "timeclass dts = ";
            for(int i = 0; i <= P::currentMaxTimeclass; ++i){
               //std::cout << i <<": "<<P::timeclassDt[i] << "s, ";
            }
            //std::cout << endl;
         }

         // restep cells
         calculateAcceleration(mpiGrid, 0.5);

         logFile << "\ncompleted timeclass reassignment after cell tc increase\n" << std::endl;

      }      

      

      //    phiprof::Timer dyndt_check {"check-and-update-dynamic-dt"};
         
      //    //check if any cell is has maximum dt lower than its timeclass dt

      //    Real tempArr1[3];
      //    Real tempArr2[3];

      //    //updating dt_r maximums
      //    //reduce_vlasov_dt(mpiGrid, cells, tempArr1, tempArr2);

      //    for (auto cID: cells) {
      //       SpatialCell* cell = mpiGrid[cID];

      //       //only check cells whose timeclass turn it is
      //       if (!(cell->get_timeclass_turn_v())) {
      //          continue;
      //       }

      //       for (ui
      //       //computeNewTimeStep(mpiGrid, technicalGrid, newDt, dtIsChanged, newTimeclassDts);
      //       balanceLoad(mpiGrid, sysBoundaryContainer, technicalGrid);

      //       P::timeclassDt.resize(P::currentMaxTimeclass+1);
      //       P::timeclassDt = newTimeclassDts;

      //       P::dt=newDt;
      //       P::timeclassDt = newTimeclassDts;

      //       additionalTimeclassCreated = false;
      //    }

      // }nt popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
      //          //updating dt_v maximums
      //          updateAccelerationMaxdt(cell, popID);

      //          // If the cell's max r_dt is lower than the timeclass dt, we need to adjust the timeclass
      //          if (cell->parameters[CellParams::MAXRDT] < cell->get_tc_dt() || 
      //              cell->parameters[CellParams::MAXVDT]* P::maxSlAccelerationSubcycles < cell->get_tc_dt()) {
      //             logFile << "(DT): Cell " << cID << " has max r_dt/v_dt lower than timeclass dt!" << endl << writeVerbose;
      //             logFile << "(DT): Cell " << cID << " has max r_dt: " << cell->parameters[CellParams::MAXRDT] << ", max v_dt times maxsubcycles: " << cell->parameters[CellParams::MAXVDT] << " * " << P::maxSlAccelerationSubcycles << " = " << cell->parameters[CellParams::MAXVDT]* P::maxSlAccelerationSubcycles << ", timeclass dt: " << cell->get_tc_dt() << endl << writeVerbose;
      //             if (myRank == MASTER_RANK) {
      //                std::cout << "Cell " << cID << " has max r_dt/v_dt lower than timeclass dt!" << std::endl;
      //                std::cout << "fractimestep: " << P::fractionalTimestep << ", timeclass of cell: " << cell->parameters[CellParams::TIMECLASS] << " get_timeclass_turn: " << cell->get_timeclass_turn_v() <<  std::endl;
      //             }

      //             if (P::fractionalTimestep == 0) {
      //                // compute new timesteps and assign all cells to timeclasses freshly
      //                logFile << "fractimestep is 0, recomputing timeclasses and dts" << endl << writeVerbose;
      //                //computeNewTimeStep(mpiGrid, technicalGrid, newDt, dtIsChanged, newTimeclassDts);
                     
      //                if (P::propagateVlasovAcceleration) {
      //                   // Back half dt to real time, forward by new half dt
      //                   calculateAcceleration(mpiGrid,-0.5); //This sets cells back to previous TIME_R
      //                   calculateAcceleration(mpiGrid, 0.5); //This propagates by 0.5 
      //                } else {
      //                   //zero step to set up moments _v
      //                   calculateAcceleration(mpiGrid, 0.0);
      //                }

      //             P::dt=newDt;
      //             P::timeclassDt = newTimeclassDts;
                  
      //             logFile <<" dt changed to "<<P::dt <<"s, distribution function was half-stepped to real-time and back"<<endl<<writeVerbose;
      //             goto endOfDynDTCheck; // All dt's and timeclasses have been recomputed, so we can continue with the next step.

      //             } else {
      //                logFile << "fractimestep is not 0, rolling back cell " << cID << " and increasing timeclass" << endl << writeVerbose;
      //                // if in the middle of a timestep, we need to roll back the offending cells and increase their timeclass
      //                calculateAcceleration(mpiGrid, -0.5, true, {cID}); // This sets cells back to previous TIME_R
      //                increaseTimeclass(mpiGrid, {cID}, additionalTimeclassCreated);
      //                calculateAcceleration(mpiGrid, 0.5, true, {cID}); // This propagates by 0.5

      //             }
      //          }
      //       }
      //    }

      //    if (additionalTimeclassCreated && P::fractionalTimestep == 0) {
      //       // If we created a new timeclass this step, we need to recompute the timeclasses and their dt
      //       // we also decrease the amount of timeclasses by one, since we just created a new one

      //       P::currentMaxTimeclass -= 1;
            
      //       //computeNewTimeStep(mpiGrid, technicalGrid, newDt, dtIsChanged, newTimeclassDts);
      //       balanceLoad(mpiGrid, sysBoundaryContainer, technicalGrid);

      //       P::timeclassDt.resize(P::currentMaxTimeclass+1);
      //       P::timeclassDt = newTimeclassDts;

      //       P::dt=newDt;
      //       P::timeclassDt = newTimeclassDts;

      //       additionalTimeclassCreated = false;
      //    }

      // }

      endOfDynDTCheck:

      // if(P::dynamicTimestep  && P::tstep > P::tstep_min && P::fractionalTimestep == 0) {
      //    std::cout << "Computing new dts\n";
      //    computeNewTimeStep(mpiGrid, technicalGrid, newDt, dtIsChanged, newTimeclassDts);
      //    // if (P::vlasovSolverGhostTranslate) {
      //    //    getGhostNeighborsforTC(mpiGrid, cells);
      //    // }
      //    if(myRank == MASTER_RANK){
      //       std::cout << "timeclass dts = ";
      //       for(int i = 0; i <= P::maxTimeclass; ++i){
      //          std::cout << i <<": " << newTimeclassDts[i] << "s, ";
      //       }
      //       std::cout << endl;
      //    }
      //    addTimedBarrier("barrier-check-dt");
      //    if(dtIsChanged) {
      //       phiprof::Timer updateDtimer {"update-dt"};
      //       //propagate velocity space back to real-time
      //       if( P::propagateVlasovAcceleration) {
      //          // Back half dt to real time, forward by new half dt
      //          calculateAcceleration(mpiGrid,-0.5); //This sets cells back to previous TIME_R
      //          calculateAcceleration(mpiGrid, 0.5); //This propagates by 0.5 
      //       }
      //       else {
      //          //zero step to set up moments _v
      //          calculateAcceleration(mpiGrid, 0.0);
      //       }

      //       P::dt=newDt;
      //       P::timeclassDt = newTimeclassDts;
            
      //       logFile <<" dt changed to "<<P::dt <<"s, distribution function was half-stepped to real-time and back"<<endl<<writeVerbose;
      //       updateDtimer.stop();
      //       continue; //
      //       //addTimedBarrier("barrier-new-dt-set");
      //    }
      //    balanceLoad(mpiGrid, sysBoundaryContainer, technicalGrid);

      // }
      
      if (((P::tstep % P::rebalanceInterval == P::rebalanceInterval-1) && (P::fractionalTimestep == 0)) || P::prepareForRebalance == true) {
         if(P::prepareForRebalance == true) {
            overrideRebalanceNow = true;
         } else {
            P::prepareForRebalance = true;
         }
         #pragma omp parallel for
         for (size_t c=0; c<cells.size(); ++c) {
            mpiGrid[cells[c]]->get_cell_parameters()[CellParams::LBWEIGHTCOUNTER] = 0;
         }
      }
      
      phiprof::Timer propagateTimer {"Propagate"};
      //Propagate the state of simulation forward in time by dt:
      
      // Update boundary condition states (time-varying)
      if (P::propagateVlasovTranslation || P::propagateVlasovAcceleration) {
         phiprof::Timer timer {"Update system boundaries (Vlasov pre-translation)"};
         sysBoundaryContainer.updateState(mpiGrid, technicalGrid, perBGrid, BgBGrid, P::t + 0.5 * P::dt);

         // updateState leaves mpiGrid and fsgrid in mismatching states, interpolated moments need to be recalculated
         // TODO: Check whether updated state is the same as previously so synchronization can be skipped when not needed?
         calculateInterpolatedVelocityMoments(
            mpiGrid,
            CellParams::RHOM,
            CellParams::VX,
            CellParams::VY,
            CellParams::VZ,
            CellParams::RHOQ,
            CellParams::P_11,
            CellParams::P_22,
            CellParams::P_33,
            CellParams::P_23,
            CellParams::P_13,
            CellParams::P_12
         );
         timer.stop();
         addTimedBarrier("barrier-boundary-conditions");
      }

      phiprof::Timer spatialSpaceTimer {"Spatial-space"};
      if( P::propagateVlasovTranslation) {
         calculateSpatialTranslation(mpiGrid,1.0,false);
      } else {
         calculateSpatialTranslation(mpiGrid,0.0,false);
      }
      spatialSpaceTimer.stop(computedCells, "Cells");
      
      // Apply boundary conditions
      if (P::propagateVlasovTranslation || P::propagateVlasovAcceleration ) {
         phiprof::Timer timer {"Update system boundaries (Vlasov post-translation)"};
         sysBoundaryContainer.applySysBoundaryVlasovConditions(mpiGrid, P::t+0.5*P::dt, false);
         timer.stop();
         addTimedBarrier("barrier-boundary-conditions");
      }
      
      phiprof::Timer momentsTimer {"Compute interp moments"};
      //std::cout << "for dt2 in main loop at t="<<P::t<<"\n";

      interpolateMomentsForTimeclasses(
         mpiGrid,
         CellParams::RHOM,
         CellParams::RHOQ,
         CellParams::P_11,
         CellParams::P_22,
         CellParams::P_33,
         CellParams::P_23,
         CellParams::P_13,
         CellParams::P_12,
         CellParams::VX,
         CellParams::VY,
         CellParams::VZ,
         false
      );

      interpolateMomentsForTimeclasses(
         mpiGrid,
         CellParams::RHOM_DT2,
         CellParams::RHOQ_DT2,
         CellParams::P_11_DT2,
         CellParams::P_22_DT2,
         CellParams::P_33_DT2,
         CellParams::P_23_DT2,
         CellParams::P_13_DT2,
         CellParams::P_12_DT2,
         CellParams::VX_DT2,
         CellParams::VY_DT2,
         CellParams::VZ_DT2,
         true
      );

      updateParticlePopulations(mpiGrid);

      // auto cell1 = mpiGrid[cells[5]];
      // auto cell2 = mpiGrid[cells[20]];

      // std::cout << "maxtc: " << P::maxTimeclass << std::endl;

      // std::cout << "cell 1 tc "<< cell1->parameters[CellParams::TIMECLASS] << " VX " << cell1->parameters[CellParams::VX] << " VY " << cell1->parameters[CellParams::VY] << " VZ " << cell1->parameters[CellParams::VZ] << std::endl;
      // std::cout << "cell 2 tc "<< cell2->parameters[CellParams::TIMECLASS] << " VX " << cell2->parameters[CellParams::VX] << " VY " << cell2->parameters[CellParams::VY] << " VZ " << cell2->parameters[CellParams::VZ] << std::endl;

      // std::cout << "cell1 vx_v " << cell1->parameters[CellParams::VX_V] << " vx_r " << cell1->parameters[CellParams::VX_R] << std::endl;
      // std::cout << "cell2 vx_v " << cell2->parameters[CellParams::VX_V] << " vx_r " << cell2->parameters[CellParams::VX_R] << std::endl;

      momentsTimer.stop();
      
      // Propagate fields forward in time by dt. This needs to be done before the
      // moments for t + dt are computed (field uses t and t+0.5dt)
      if (P::propagateField) {
         phiprof::Timer propagateTimer {"Propagate Fields"};

         phiprof::Timer couplingInTimer {"fsgrid-coupling-in"};
         // Copy moments over into the fsgrid.
         //setupTechnicalFsGrid(mpiGrid, cells, technicalGrid);
         feedMomentsIntoFsGrid(mpiGrid, cells, momentsGrid, technicalGrid, false);
         feedMomentsIntoFsGrid(mpiGrid, cells, momentsDt2Grid, technicalGrid, true);
         couplingInTimer.stop();
         
         propagateFields(
            perBGrid,
            perBDt2Grid,
            EGrid,
            EDt2Grid,
            EHallGrid,
            EGradPeGrid,
            EGradPeDt2Grid,
            momentsGrid,
            momentsDt2Grid,
            dPerBGrid,
            dMomentsGrid,
            dMomentsDt2Grid,
            BgBGrid,
            volGrid,
            technicalGrid,
            sysBoundaryContainer,
            P::dt,
            P::fieldSolverSubcycles
         );

         phiprof::Timer getFieldsTimer {"getFieldsFromFsGrid"};
         // Copy results back from fsgrid.
         volGrid.updateGhostCells();
         technicalGrid.updateGhostCells();
         getFieldsFromFsGrid(volGrid, BgBGrid, EGradPeGrid, dMomentsGrid, technicalGrid, mpiGrid, cells);
         getFieldsTimer.stop();
         propagateTimer.stop(cells.size(),"SpatialCells");
         addTimedBarrier("barrier-after-field-solver");
      }

      if(FieldTracing::fieldTracingParameters.useCache) {
         FieldTracing::resetReconstructionCoefficientsCache();
      }

      // Map current data down into the ionosphere
      // momentsGrid was ghost-updated in the field solver above, volGrid just after a few lines above.
      // perBGrid was ghost-updated before derivatives were computed in the field solver.
      // dPerBGrid was updated before the electric fields.
      if(SBC::ionosphereGrid.nodes.size() > 0 && ((P::t > SBC::Ionosphere::solveCount * SBC::Ionosphere::couplingInterval && SBC::Ionosphere::couplingInterval > 0) || SBC::Ionosphere::couplingInterval == 0)) {
         FieldTracing::calculateIonosphereFsgridCoupling(technicalGrid, perBGrid, dPerBGrid, SBC::ionosphereGrid.nodes, SBC::Ionosphere::radius);
         SBC::ionosphereGrid.mapDownBoundaryData(perBGrid, dPerBGrid, momentsGrid, volGrid, technicalGrid);
         SBC::ionosphereGrid.calculateConductivityTensor(SBC::Ionosphere::F10_7, SBC::Ionosphere::recombAlpha, SBC::Ionosphere::backgroundIonisation);

         // Solve ionosphere
         int nIterations, nRestarts;
         Real residual, minPotentialN, maxPotentialN, minPotentialS, maxPotentialS;
         SBC::ionosphereGrid.solve(nIterations, nRestarts, residual, minPotentialN, maxPotentialN, minPotentialS, maxPotentialS);
         logFile << "tstep = " << P::tstep << "("<<P::fractionalTimestep<<"/"<< (1u << (P::currentMaxTimeclass)) << ")"
         << " t = " << P::t
         << " ionosphere iterations = " << nIterations
         << " restarts = " << nRestarts
         << " residual = " << std::scientific << residual << std::defaultfloat
         << " N potential min " << minPotentialN
         << " max " << maxPotentialN
         << " difference " << maxPotentialN - minPotentialN
         << " S potential min " << minPotentialS
         << " max " << maxPotentialS
         << " difference " << maxPotentialS - minPotentialS
         << endl;
         SBC::Ionosphere::solveCount++;
         globalflags::ionosphereJustSolved = true;
      }

      // updating _V_PREV moments here, before _V moments are updated 
      updatePreviousVMoments(mpiGrid, false);
      
      phiprof::Timer vspaceTimer {"Velocity-space"};
      if ( P::propagateVlasovAcceleration ) {
      // calculateAcceleration(mpiGrid,P::dt);
         calculateAcceleration(mpiGrid,1.0);
         addTimedBarrier("barrier-after-ad just-blocks");
      } else {
         //zero step to set up moments _v
         calculateAcceleration(mpiGrid, 0.0);
      }
      vspaceTimer.stop(computedCells, "Cells");
      addTimedBarrier("barrier-after-acceleration");

      if (P::artificialPADiff){
         phiprof::Timer diffusionTimer {"Pitch-angle diffusion"};
         for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
	      pitchAngleDiffusion(mpiGrid,popID);
         }
         diffusionTimer.stop(computedCells, "Cells");
      }

      if (P::propagateVlasovTranslation || P::propagateVlasovAcceleration ) {
         phiprof::Timer timer {"Update system boundaries (Vlasov post-acceleration)"};
         sysBoundaryContainer.applySysBoundaryVlasovConditions(mpiGrid, P::t + 0.5 * P::dt, true);
         timer.stop();
         addTimedBarrier("barrier-boundary-conditions");
      }
      
      // momentsTimer.start();
      // *here we compute rho and rho_v for timestep t + dt, so next
      // timestep * //
      // calculateInterpolatedVelocityMoments(
      //    mpiGrid,
      //    CellParams::RHOM,
      //    CellParams::VX,
      //    CellParams::VY,
      //    CellParams::VZ,
      //    CellParams::RHOQ,
      //    CellParams::P_11,
      //    CellParams::P_22,
      //    CellParams::P_33,
      //    CellParams::P_23,
      //    CellParams::P_13,
      //    CellParams::P_12
      // );
      
      // updateParticlePopulations(mpiGrid);

      // momentsTimer.stop();

      propagateTimer.stop(computedCells,"Cells");
      
      phiprof::Timer endStepTimer {"Project endTimeStep"};
      project->hook(hook::END_OF_TIME_STEP, mpiGrid, perBGrid);
      endStepTimer.stop();

      // Check timestep
      if (P::dt < P::bailout_min_dt) {
         stringstream s;
         s << "The timestep dt=" << P::dt << " went below bailout.min_dt (" << to_string(P::bailout_min_dt) << ")." << endl;
         bailout(true, s.str(), __FILE__, __LINE__);
      }
      //Move forward in time
      P::meshRepartitioned = false;
      globalflags::ionosphereJustSolved = false;
      ++P::fractionalTimestep;
      if(P::fractionalTimestep % (1u << (P::currentMaxTimeclass)) == 0){
         ++P::tstep;
         P::fractionalTimestep = 0;
      }
      P::t += P::dt;
   } // End main loop ----------------------------------------------------------

   double after = MPI_Wtime();

   simulationTimer.stop();
   phiprof::Timer finalizationTimer {"Finalization"};
   if (myRank == MASTER_RANK) {
      if (doBailout > 0) {
         logFile << "(BAILOUT): Bailing out, see error log for details." << endl;
      }

      double timePerStep;

      if (P::tstep == P::tstep_min) {
         timePerStep=0.0;
      } else {
         timePerStep=double(after  - startTime) / (P::tstep-P::tstep_min);
      }
      double timePerSecond=double(after  - startTime) / (P::t-P::t_min+DT_EPSILON);
      logFile << "(MAIN): All timesteps calculated." << endl;
      logFile << "\t (TIME) total run time " << after - startTime << " s, total simulated time " << P::t -P::t_min<< " s" << endl;
      logFile << "\t (TIME) total " << nNodes*(after - startTime)/3600 << " node-hours" << endl;
      #if _OPENMP
         logFile << "\t (TIME) total " << omp_get_max_threads()*mpiProcs*(after - startTime)/3600 << " thread-hours" << endl;
      #endif

      if(P::t != 0.0) {
         logFile << "\t (TIME) seconds per timestep " << timePerStep  <<
         ", seconds per simulated second " <<  timePerSecond << endl;
      }
      logFile << writeVerbose;
   }
   
   finalizationTimer.stop();
   mainTimer.stop();

   #ifdef USE_GPU
   // Deallocate buffers, clear device
   vmesh::deallocateMeshWrapper();
   gpu_clear_device();
   getObjectWrapper().sysBoundaryContainer.clear();
   #endif
   
   phiprof::print(MPI_COMM_WORLD,"phiprof");
   
   if (myRank == MASTER_RANK) {
      logFile << "(MAIN): Completed requested simulation. Exiting." << endl << writeVerbose;
      cout << "(MAIN): Completed requested simulation. Exiting." << endl;
   }
   logFile.close();
   if (P::diagnosticInterval != 0) {
      diagnostic.close();
   }

   return 0;
}

int main(int argn, char* args[]) {
   // Before MPI_Init we hardwire some settings, if we are in OpenMPI
   int myRank;
   int required=MPI_THREAD_FUNNELED;
   int provided, resultlen;
   char mpiversion[MPI_MAX_LIBRARY_VERSION_STRING];
   bool overrideMCAompio = false;

   MPI_Get_library_version(mpiversion, &resultlen);
   string versionstr = string(mpiversion);
   stringstream mpiioMessage;

   if(versionstr.find("Open MPI") != std::string::npos) {
      #ifdef VLASIATOR_ALLOW_MCA_OMPIO
         mpiioMessage << "We detected OpenMPI but the compilation flag VLASIATOR_ALLOW_MCA_OMPIO was set so we do not override the default MCA io flag." << endl;
      #else
         overrideMCAompio = true;
         int index, count;
         char io_value[64];
         MPI_T_cvar_handle io_handle;
         
         MPI_T_init_thread(required, &provided);
         MPI_T_cvar_get_index("io", &index);
         MPI_T_cvar_handle_alloc(index, NULL, &io_handle, &count);
         MPI_T_cvar_write(io_handle, "^ompio");
         MPI_T_cvar_read(io_handle, io_value);
         MPI_T_cvar_handle_free(&io_handle);
         mpiioMessage << "We detected OpenMPI so we set the cvars value to disable ompio, MCA io: " << io_value << endl;
      #endif
   }
   
   // After the MPI_T settings we can init MPI all right.
   MPI_Init_thread(&argn,&args,required,&provided);
   MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
   if (required > provided){
      if(myRank==MASTER_RANK) {
         cerr << "(MAIN): MPI_Init_thread failed! Got " << provided << ", need "<<required <<endl;
      }
      exit(1);
   }
   if (myRank == MASTER_RANK) {
      const char* mpiioenv = std::getenv("OMPI_MCA_io");
      if(mpiioenv != nullptr) {
         std::string mpiioenvstr(mpiioenv);
         if(mpiioenvstr.find("^ompio") == std::string::npos) {
            cout << mpiioMessage.str();
         }
      }
   }

   int ret {simulate(argn, args)};

   if(overrideMCAompio) {
      MPI_T_finalize();
   }
   MPI_Finalize();

   return ret;
}
