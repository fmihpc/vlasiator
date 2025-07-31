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
#include <cstdlib>
#include <iostream>
#include <cmath>
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

void computeNewTimeStep(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
			FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid, Real &newDt, bool &isChanged) {

   phiprof::Timer computeTimestepTimer {"compute-timestep"};
   // Compute maximum time step. This cannot be done at the first step as the solvers compute the limits for each cell.

   isChanged = false;

   const vector<CellID>& cells = getLocalCells();
   /* Arrays for storing local (per process) and global max dt
      0th position stores ordinary space propagation dt
      1st position stores velocity space propagation dt
      2nd position stores field propagation dt
   */
   Real dtMaxLocal[3];
   Real dtMaxGlobal[3];

   dtMaxLocal[0] = numeric_limits<Real>::max();
   dtMaxLocal[1] = numeric_limits<Real>::max();
   dtMaxLocal[2] = numeric_limits<Real>::max();

   // Compute max dt for Vlasov solver
   reduce_vlasov_dt(mpiGrid, cells, dtMaxLocal);

   // compute max dt for fieldsolver
   const std::array<FsGridTools::FsIndex_t, 3> gridDims(technicalGrid.getLocalSize());
   for (FsGridTools::FsIndex_t k = 0; k < gridDims[2]; k++) {
      for (FsGridTools::FsIndex_t j = 0; j < gridDims[1]; j++) {
         for (FsGridTools::FsIndex_t i = 0; i < gridDims[0]; i++) {
            fsgrids::technical* cell = technicalGrid.get(i, j, k);
            if (cell->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY ||
               (cell->sysBoundaryLayer == 1 && cell->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY)) {
               dtMaxLocal[2] = min(dtMaxLocal[2], cell->maxFsDt);
            }
         }
      }
   }

   MPI_Allreduce(&(dtMaxLocal[0]), &(dtMaxGlobal[0]), 3, MPI_Type<Real>(), MPI_MIN, MPI_COMM_WORLD);

   // If any of the solvers are disabled there should be no limits in timespace from it
   if (!P::propagateVlasovTranslation)
      dtMaxGlobal[0] = numeric_limits<Real>::max();
   if (!P::propagateVlasovAcceleration)
      dtMaxGlobal[1] = numeric_limits<Real>::max();
   if (!P::propagateField)
      dtMaxGlobal[2] = numeric_limits<Real>::max();

   creal meanVlasovCFL = 0.5 * (P::vlasovSolverMaxCFL + P::vlasovSolverMinCFL);
   creal meanFieldsCFL = 0.5 * (P::fieldSolverMaxCFL + P::fieldSolverMinCFL);
   Real subcycleDt;

   // reduce/increase dt if it is too high for any of the three propagators or too low for all propagators
   if ((P::dt > dtMaxGlobal[0] * P::vlasovSolverMaxCFL ||
        P::dt > dtMaxGlobal[1] * P::vlasovSolverMaxCFL * P::maxSlAccelerationSubcycles ||
        P::dt > dtMaxGlobal[2] * P::fieldSolverMaxCFL * P::maxFieldSolverSubcycles) ||
       (P::dt < dtMaxGlobal[0] * P::vlasovSolverMinCFL &&
        P::dt < dtMaxGlobal[1] * P::vlasovSolverMinCFL * P::maxSlAccelerationSubcycles &&
        P::dt < dtMaxGlobal[2] * P::fieldSolverMinCFL * P::maxFieldSolverSubcycles)) {

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

int simulate(int argn,char* args[]) {
   int myRank, doBailout=0;
   const creal DT_EPSILON=1e-12;
   typedef Parameters P;
   Real newDt;
   bool dtIsChanged {false};
   
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
   calculateSpatialTranslation(mpiGrid,0.0);
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

   if (P::isRestart == false) {
      //compute new dt
      phiprof::Timer computeDtTimer {"compute-dt"};
      computeNewTimeStep(mpiGrid, technicalGrid, newDt, dtIsChanged);
      if (P::dynamicTimestep == true && dtIsChanged == true) {
         // Only actually update the timestep if dynamicTimestep is on
         P::dt=newDt;
      } else {
         dtIsChanged = false;
      }
      computeDtTimer.stop();
      
      //go forward by dt/2 in V, initializes leapfrog split. In restarts the
      //the distribution function is already propagated forward in time by dt/2
      phiprof::Timer propagateHalfTimer {"propagate-velocity-space-dt/2"};
      if (P::propagateVlasovAcceleration) {
         calculateAcceleration(mpiGrid, 0.5*P::dt);
      } else {
         //zero step to set up moments _v
         calculateAcceleration(mpiGrid, 0.0);
      }
      propagateHalfTimer.stop();

      // Apply boundary conditions
      if (P::propagateVlasovTranslation || P::propagateVlasovAcceleration ) {
         phiprof::Timer updateBoundariesTimer {("update system boundaries (Vlasov post-acceleration)")};
         sysBoundaryContainer.applySysBoundaryVlasovConditions(mpiGrid, 0.5*P::dt, true);
         updateBoundariesTimer.stop();
         addTimedBarrier("barrier-boundary-conditions");
      }
      // Also update all moments. They won't be transmitted to FSgrid until the field solver is called, though.
      phiprof::Timer computeMomentsTimer {"Compute interp moments"};
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
      computeMomentsTimer.stop();
   }

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
      logFile << "---------- tstep = " << P::tstep << " t = " << P::t <<" dt = " << P::dt << " FS cycles = " << P::fieldSolverSubcycles << " ----------" << endl;
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
      if (P::diagnosticInterval != 0 && P::tstep % P::diagnosticInterval == 0) {
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
               || P::tstep == P::tstep_max
               || P::t >= P::t_max))
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

      //Re-loadbalance if needed
      //TODO - add LB measure and do LB if it exceeds threshold
      if(((P::tstep % P::rebalanceInterval == 0 && P::tstep > P::tstep_min) || overrideRebalanceNow)) {
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
            calculateSpatialTranslation(mpiGrid,0.0);
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
      // FIXME what if dt changes at a restart??
      if(P::dynamicTimestep  && P::tstep > P::tstep_min) {
         computeNewTimeStep(mpiGrid, technicalGrid, newDt, dtIsChanged);
         addTimedBarrier("barrier-check-dt");
         if(dtIsChanged) {
            phiprof::Timer updateDtimer {"update-dt"};
            //propagate velocity space back to real-time
            if( P::propagateVlasovAcceleration ) {
               // Back half dt to real time, forward by new half dt
               calculateAcceleration(mpiGrid,-0.5*P::dt + 0.5*newDt);
            }
            else {
               //zero step to set up moments _v
               calculateAcceleration(mpiGrid, 0.0);
            }

            P::dt=newDt;

            logFile <<" dt changed to "<<P::dt <<"s, distribution function was half-stepped to real-time and back"<<endl<<writeVerbose;
            updateDtimer.stop();
            continue; //
            //addTimedBarrier("barrier-new-dt-set");
         }
      }

      if (P::tstep % P::rebalanceInterval == P::rebalanceInterval-1 || P::prepareForRebalance == true) {
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
         sysBoundaryContainer.updateState(mpiGrid, perBGrid, BgBGrid, P::t + 0.5 * P::dt);
         timer.stop();
         addTimedBarrier("barrier-boundary-conditions");
      }

      phiprof::Timer spatialSpaceTimer {"Spatial-space"};
      if( P::propagateVlasovTranslation) {
         calculateSpatialTranslation(mpiGrid,P::dt);
      } else {
         calculateSpatialTranslation(mpiGrid,0.0);
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
      calculateInterpolatedVelocityMoments(
         mpiGrid,
         CellParams::RHOM_DT2,
         CellParams::VX_DT2,
         CellParams::VY_DT2,
         CellParams::VZ_DT2,
         CellParams::RHOQ_DT2,
         CellParams::P_11_DT2,
         CellParams::P_22_DT2,
         CellParams::P_33_DT2,
         CellParams::P_23_DT2,
         CellParams::P_13_DT2,
         CellParams::P_12_DT2
      );
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
         logFile << "tstep = " << P::tstep
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
      
      phiprof::Timer vspaceTimer {"Velocity-space"};
      if ( P::propagateVlasovAcceleration ) {
         calculateAcceleration(mpiGrid,P::dt);
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
      
      momentsTimer.start();
      // *here we compute rho and rho_v for timestep t + dt, so next
      // timestep * //
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
      momentsTimer.stop();

      propagateTimer.stop(computedCells,"Cells");
      
      phiprof::Timer endStepTimer {"Project endTimeStep"};
      project->hook(hook::END_OF_TIME_STEP, mpiGrid, perBGrid);
      endStepTimer.stop();

      // Check timestep
      if (P::dt < P::bailout_min_dt) {
         stringstream s;
         s << "The timestep dt=" << P::dt << " went below bailout.bailout_min_dt (" << to_string(P::bailout_min_dt) << ")." << endl;
         bailout(true, s.str(), __FILE__, __LINE__);
      }
      //Move forward in time
      P::meshRepartitioned = false;
      globalflags::ionosphereJustSolved = false;
      ++P::tstep;
      P::t += P::dt;

   }

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
