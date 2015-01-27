/*
This file is part of Vlasiator.

Copyright 2010, 2011, 2012, 2013, 2014 Finnish Meteorological Institute

*/

#include <cstdlib>
#include <iostream>
#include <cmath>
#include <vector>
#include <sstream>
#include <ctime>
#include <omp.h>

#include "vlasovmover.h"
#include "definitions.h"
#include "mpiconversion.h"
#include "logger.h"
#include "parameters.h"
#include "readparameters.h"
#include "spatial_cell.hpp"
#include "datareduction/datareducer.h"
#include "sysboundary/sysboundary.h"

#include "fieldsolver/fs_common.h"
#include "poisson_solver/poisson_solver.h"
#include "projects/project.h"
#include "grid.h"
#include "iowrite.h"
#include "ioread.h"

#include "object_wrapper.h"

#ifdef CATCH_FPE
#include <fenv.h>
#include <signal.h>
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
static dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry> mpiGrid;

using namespace std;
using namespace phiprof;

int globalflags::bailingOut = 0;

ObjectWrapper objectWrapper;

void addTimedBarrier(string name){
#ifdef NDEBUG
//let's not do  a barrier
   return; 
#endif
   int bt=phiprof::initializeTimer(name,"Barriers","MPI");
   phiprof::start(bt);
   MPI_Barrier(MPI_COMM_WORLD);
   phiprof::stop(bt);
}

bool computeNewTimeStep(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,Real &newDt, bool &isChanged) {

   phiprof::start("compute-timestep");
   //compute maximum time-step, this cannot be done at the first
   //step as the solvers compute the limits for each cell

   isChanged=false;

   vector<uint64_t> cells = mpiGrid.get_cells();
   /* Arrays for storing local (per process) and global max dt
      0th position stores ordinary space propagation dt
      1st position stores velocity space propagation dt
      2nd position stores field propagation dt
   */
   Real dtMaxLocal[3];
   Real dtMaxGlobal[3];
   
   dtMaxLocal[0]=std::numeric_limits<Real>::max();
   dtMaxLocal[1]=std::numeric_limits<Real>::max();
   dtMaxLocal[2]=std::numeric_limits<Real>::max();

   for (std::vector<uint64_t>::const_iterator cell_id = cells.begin(); cell_id != cells.end(); ++cell_id) {
      SpatialCell* cell = mpiGrid[*cell_id];
      if ( cell->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY ||
           (cell->sysBoundaryLayer == 1 && cell->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY )) {
         //spatial fluxes computed also for boundary cells
         dtMaxLocal[0]=min(dtMaxLocal[0], cell->parameters[CellParams::MAXRDT]);
         dtMaxLocal[2]=min(dtMaxLocal[2], cell->parameters[CellParams::MAXFDT]);
      }
      
      if (cell->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
         //Acceleration only done on non sysboundary cells
         dtMaxLocal[1]=min(dtMaxLocal[1], cell->parameters[CellParams::MAXVDT]);
      }
   }
   MPI_Allreduce(&(dtMaxLocal[0]), &(dtMaxGlobal[0]), 3, MPI_Type<Real>(), MPI_MIN, MPI_COMM_WORLD);
   
   //If any of the solvers are disabled there should be no limits in timespace from it
   if (P::propagateVlasovTranslation == false)
      dtMaxGlobal[0]=std::numeric_limits<Real>::max();
   if (P::propagateVlasovAcceleration == false)
      dtMaxGlobal[1]=std::numeric_limits<Real>::max();
   if (P::propagateField == false)
      dtMaxGlobal[2]=std::numeric_limits<Real>::max();
   
   
   creal meanVlasovCFL = 0.5*(P::vlasovSolverMaxCFL+ P::vlasovSolverMinCFL);
   creal meanFieldsCFL = 0.5*(P::fieldSolverMaxCFL+ P::fieldSolverMinCFL);
   Real subcycleDt;
   
   //reduce dt if it is too high for any of the three propagators, or too low for all propagators
   if(( P::dt > dtMaxGlobal[0] * P::vlasovSolverMaxCFL ||
        P::dt > dtMaxGlobal[1] * P::vlasovSolverMaxCFL * P::maxSlAccelerationSubcycles ||
        P::dt > dtMaxGlobal[2] * P::fieldSolverMaxCFL * P::maxFieldSolverSubcycles ) ||
      ( P::dt < dtMaxGlobal[0] * P::vlasovSolverMinCFL && 
        P::dt < dtMaxGlobal[1] * P::vlasovSolverMinCFL * P::maxSlAccelerationSubcycles &&
        P::dt < dtMaxGlobal[2] * P::fieldSolverMinCFL * P::maxFieldSolverSubcycles )
      ) {
      //new dt computed
      isChanged=true;

      //set new timestep to the lowest one of all interval-midpoints
      const Real half = 0.5;
      newDt = meanVlasovCFL * dtMaxGlobal[0];
      newDt = min(newDt,meanVlasovCFL * dtMaxGlobal[1] * P::maxSlAccelerationSubcycles);
      newDt = min(newDt,meanFieldsCFL * dtMaxGlobal[2] * P::maxFieldSolverSubcycles);
   
      logFile <<"(TIMESTEP) New dt = " << newDt << " computed on step "<<  P::tstep <<" at " <<P::t << 
         "s   Maximum possible dt (not including  vlasovsolver CFL "<< 
         P::vlasovSolverMinCFL <<"-"<<P::vlasovSolverMaxCFL<<
         " or fieldsolver CFL "<< 
         P::fieldSolverMinCFL <<"-"<<P::fieldSolverMaxCFL<<
         ") in {r, v, BE} was " <<
         dtMaxGlobal[0] << " " <<
         dtMaxGlobal[1] << " " <<
         dtMaxGlobal[2] << " " <<
         " Including subcycling { v, BE}  was " <<
         dtMaxGlobal[1] * P::maxSlAccelerationSubcycles << " " <<
         dtMaxGlobal[2] * P::maxFieldSolverSubcycles<< " " <<
         endl << writeVerbose;
      subcycleDt = newDt;
   } else {
      subcycleDt = P::dt;
   }
   
   // Subcycle if field solver dt < global dt (including CFL) (new or old dt hence the hassle with subcycleDt
   if (meanFieldsCFL*dtMaxGlobal[2] < subcycleDt && P::propagateField) {
      P::fieldSolverSubcycles = min(convert<int>(ceil(subcycleDt / (meanFieldsCFL*dtMaxGlobal[2]))), P::maxFieldSolverSubcycles);
   } else {
      P::fieldSolverSubcycles = 1;
   }
   
   phiprof::stop("compute-timestep");
   return true;
}

ObjectWrapper& getObjectWrapper() {
   return objectWrapper;
}

/** Get local cell IDs. This function creates a cached copy of the 
 * cell ID lists to significantly improve performance. The cell ID 
 * cache is recalculated every time the mesh partitioning changes.
 * @return Local cell IDs.*/
const std::vector<CellID>& getLocalCells() {
   if (Parameters::meshRepartitioned == true) {
      if (Parameters::localCellsCalculated != Parameters::tstep) {
         {
            vector<CellID> dummy;
            dummy.swap(Parameters::localCells);
         }
         Parameters::localCells = mpiGrid.get_cells();
      }
   }
   return Parameters::localCells;
}

int main(int argn,char* args[]) {
   bool success = true;
   int myRank, doBailout;
   const creal DT_EPSILON=1e-12;
   typedef Parameters P;
   Real newDt;
   bool dtIsChanged;
   
// Init MPI:
   int required=MPI_THREAD_FUNNELED;
   int provided;
   MPI_Init_thread(&argn,&args,required,&provided);
   if (required > provided){
      MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
      if(myRank==MASTER_RANK)
         cerr << "(MAIN): MPI_Init_thread failed! Got " << provided << ", need "<<required <<endl;
      exit(1);
   }
   
   double initialWtime =  MPI_Wtime();
   
   MPI_Comm comm = MPI_COMM_WORLD;
   MPI_Comm_rank(comm,&myRank);
   SysBoundary sysBoundaries;
   bool isSysBoundaryCondDynamic;

   #ifdef CATCH_FPE
   // WARNING FE_INEXACT is too sensitive to be used. See man fenv.
   //feenableexcept(FE_DIVBYZERO|FE_INVALID|FE_OVERFLOW|FE_UNDERFLOW);
   feenableexcept(FE_DIVBYZERO|FE_INVALID|FE_OVERFLOW);
   //feenableexcept(FE_DIVBYZERO|FE_INVALID);
   signal(SIGFPE, fpehandler);
   #endif

   phiprof::start("main");
   phiprof::start("Initialization");
   phiprof::start("Read parameters");
   //init parameter file reader
   Readparameters readparameters(argn,args,MPI_COMM_WORLD);
   P::addParameters();
   projects::Project::addParameters();
   sysBoundaries.addParameters();
   readparameters.parse();
   P::getParameters();
   
   Project* project = projects::createProject();
   
   project->getParameters();
   sysBoundaries.getParameters();
   phiprof::stop("Read parameters");
   
   phiprof::start("Init project");
   if (project->initialize() == false) {
      if(myRank == MASTER_RANK) cerr << "(MAIN): Project did not initialize correctly!" << endl;
      exit(1);
   }
   phiprof::stop("Init project");

   // Init parallel logger:
   phiprof::start("open logFile & diagnostic");
   //if restarting we will append to logfiles
   if (logFile.open(MPI_COMM_WORLD,MASTER_RANK,"logfile.txt",P::isRestart) == false) {
      if(myRank == MASTER_RANK) cerr << "(MAIN) ERROR: Logger failed to open logfile!" << endl;
      exit(1);
   } else {
      int mpiProcesses;
      int threads = omp_get_max_threads();
      MPI_Comm_size(comm,&mpiProcesses);      
      logFile << "(MAIN) Starting simulation with " << mpiProcesses << " MPI processes and ";
      logFile << threads << " OpenMP threads per process" << endl << writeVerbose;
   }
   if (P::diagnosticInterval != 0) {
      if (diagnostic.open(MPI_COMM_WORLD,MASTER_RANK,"diagnostic.txt",P::isRestart) == false) {
         if(myRank == MASTER_RANK) cerr << "(MAIN) ERROR: Logger failed to open diagnostic file!" << endl;
         exit(1);
      }
   }
   phiprof::stop("open logFile & diagnostic");
   
   // Add AMR refinement criterias:
   amr_ref_criteria::addRefinementCriteria();
   
   phiprof::start("Init grid");
   /* Initialize grid.  After initializeGrid local cells have dist
      functions, and B fields set. Cells have also been classified for
      the various sys boundary conditions.  All remote cells have been
      created. All spatial date computed this far is up to date for
      FULL_NEIGHBORHOOD. Block lists up to date for
      VLASOV_SOLVER_NEIGHBORHOOD (but dist function has not been communicated)
   */
   initializeGrid(argn,args,mpiGrid,sysBoundaries,*project);
   isSysBoundaryCondDynamic = sysBoundaries.isDynamic();
   phiprof::stop("Init grid");
   phiprof::start("Init DROs");
   // Initialize data reduction operators. This should be done elsewhere in order to initialize 
   // user-defined operators:
   DataReducer outputReducer, diagnosticReducer;
   initializeDataReducers(&outputReducer, &diagnosticReducer);
   phiprof::stop("Init DROs");

   // Initialize field propagator:
   if (P::propagateField ) { 
      phiprof::start("Init field propagator");
      if (initializeFieldPropagator(mpiGrid, sysBoundaries) == false) {
         logFile << "(MAIN): Field propagator did not initialize correctly!" << endl << writeVerbose;
         exit(1);
      }
      phiprof::stop("Init field propagator");
   }
   
   // Initialize Poisson solver (if used)
   if (P::propagatePotential == true) {
      phiprof::start("Init Poisson solver");
      if (poisson::initialize(mpiGrid) == false) {
         logFile << "(MAIN): Poisson solver did not initialize correctly!" << endl << writeVerbose;
         exit(1);
      }
      phiprof::stop("Init Poisson solver");
   }

   // Free up memory:
   readparameters.finalize();
   
   // Save restart data
   if (P::writeInitialState) {
      phiprof::start("write-initial-state");
      if (myRank == MASTER_RANK)
         logFile << "(IO): Writing initial state to disk, tstep = "  << endl << writeVerbose;
//    P::systemWriteDistributionWriteStride[i], P::systemWriteName[i], P::systemWrites[i]
      P::systemWriteDistributionWriteStride.push_back(1);
      P::systemWriteName.push_back("initial-grid");
      P::systemWriteDistributionWriteXlineStride.push_back(0);
      P::systemWriteDistributionWriteYlineStride.push_back(0);
      P::systemWriteDistributionWriteZlineStride.push_back(0);
      
      for(uint si=0; si<P::systemWriteName.size(); si++) {
         P::systemWrites.push_back(0);
      }
      
      const bool writeGhosts = true;
      if( writeGrid(mpiGrid,outputReducer,P::systemWriteName.size()-1, writeGhosts) == false ) {
         cerr << "FAILED TO WRITE GRID AT" << __FILE__ << " " << __LINE__ << endl;
      }
      
      P::systemWriteDistributionWriteStride.pop_back();
      P::systemWriteName.pop_back();
      P::systemWriteDistributionWriteXlineStride.pop_back();
      P::systemWriteDistributionWriteYlineStride.pop_back();
      P::systemWriteDistributionWriteZlineStride.pop_back();
      
      phiprof::stop("write-initial-state");
   }

   if (P::dynamicTimestep && !P::isRestart) {
      //compute vlasovsolver once with zero dt, this is to initialize
      //per-cell dt limits. In restarts, we read in dt from file
      phiprof::start("compute-dt");
      calculateSpatialTranslation(mpiGrid,0.0);
      calculateAcceleration(mpiGrid,0.0);

      
      if(P::propagateField) {
         propagateFields(mpiGrid, sysBoundaries, 0.0, 1.0);
      }
      //compute new dt
      computeNewTimeStep(mpiGrid,newDt,dtIsChanged);
      if(dtIsChanged)
         P::dt=newDt;
      phiprof::stop("compute-dt");
   }

   if (!P::isRestart) {
      //go forward by dt/2 in V, initializes leapfrog split. In restarts the
      //the distribution function is already propagated forward in time by dt/2
      phiprof::start("propagate-velocity-space-dt/2");
      if (P::propagateVlasovAcceleration) {
         calculateAcceleration(mpiGrid, 0.5*P::dt);
      } else {
         //zero step to set up moments _v
         calculateAcceleration(mpiGrid, 0.0);
      }
      phiprof::stop("propagate-velocity-space-dt/2");

   }
   phiprof::stop("Initialization");


   // ***********************************
   // ***** INITIALIZATION COMPLETE *****
   // ***********************************
   
   // Main simulation loop:
   if (myRank == MASTER_RANK) logFile << "(MAIN): Starting main simulation loop." << endl << writeVerbose;

   unsigned int computedCells=0;
   unsigned int computedTotalCells=0;
  //Compute here based on time what the file intervals are
   P::systemWrites.clear();
   for(uint i=0;i< P::systemWriteTimeInterval.size();i++){
      int index=(int)(P::t_min/P::systemWriteTimeInterval[i]);
      //if we are already over 1% further than the time interval time that
      //is requested for writing, then jump to next writing index. This is to
      //make sure that at restart we do not write in the middle of
      //the interval.
      if(P::t_min>(index+0.01)*P::systemWriteTimeInterval[i])
         index++;
      P::systemWrites.push_back(index);
   }

   // Invalidate cached cell lists just to be sure (might not be needed)
   P::meshRepartitioned = true;

   unsigned int wallTimeRestartCounter=1;
   
   addTimedBarrier("barrier-end-initialization");
   
   phiprof::start("Simulation");
   double startTime=  MPI_Wtime();
   double beforeTime = MPI_Wtime();
   double beforeSimulationTime=P::t_min;
   double beforeStep=P::tstep_min;
   
   while(P::tstep <=P::tstep_max  &&
         P::t-P::dt <= P::t_max+DT_EPSILON &&
         wallTimeRestartCounter <= P::exitAfterRestarts) {
      
      addTimedBarrier("barrier-loop-start");
      
      phiprof::start("IO");
      
      if(myRank ==  MASTER_RANK) {
         // check whether STOP or KILL has been passed, should be done by MASTER_RANK only as it can reset P::bailout_write_restart
         checkExternalCommands();
      }
      
      //write out phiprof profiles and logs with a lower interval than normal
      //diagnostic (every 10 diagnostic intervals).
      logFile << "---------- tstep = " << P::tstep << " t = " << P::t <<" dt = " << P::dt << " FS cycles = " << P::fieldSolverSubcycles << " ----------" << endl;
      if (P::diagnosticInterval != 0 &&
          P::tstep % (P::diagnosticInterval*10) == 0 &&
          P::tstep-P::tstep_min >0) {
         MPI_Barrier(MPI_COMM_WORLD);
         phiprof::print(MPI_COMM_WORLD,"phiprof_reduced",0.01);
         // MPI_Barrier(MPI_COMM_WORLD);
         //phiprof::print(MPI_COMM_WORLD,"phiprof_full");
         // phiprof::printLogProfile(MPI_COMM_WORLD,P::tstep,"phiprof_log"," ",7);
         
         double currentTime=MPI_Wtime();
         double timePerStep=double(currentTime  - beforeTime) / (P::tstep-beforeStep);
         double timePerSecond=double(currentTime  - beforeTime) / (P::t-beforeSimulationTime + DT_EPSILON);
         double remainingTime=min(timePerStep*(P::tstep_max-P::tstep),timePerSecond*(P::t_max-P::t));
         time_t finalWallTime=time(NULL)+(time_t)remainingTime; //assume time_t is in seconds, as it is almost always
         struct tm *finalWallTimeInfo=localtime(&finalWallTime);
         logFile << "(TIME) current walltime/step " << timePerStep<< " s" <<endl;
         logFile << "(TIME) current walltime/simusecond " << timePerSecond<<" s" <<endl;
         logFile << "(TIME) Estimated completion time is " <<asctime(finalWallTimeInfo)<<endl;
         //reset before values, we want to report speed since last report of speed.
         beforeTime = MPI_Wtime();
         beforeSimulationTime=P::t;
         beforeStep=P::tstep;
         //report_grid_memory_consumption(mpiGrid);
         report_process_memory_consumption();
      }
      logFile << writeVerbose;      
      
      // Check whether diagnostic output has to be produced
      if (P::diagnosticInterval != 0 && P::tstep % P::diagnosticInterval == 0) {
         phiprof::start("Diagnostic");
         if (writeDiagnostic(mpiGrid, diagnosticReducer) == false) {
            if(myRank == MASTER_RANK)  cerr << "ERROR with diagnostic computation" << endl;
            
         }
         phiprof::stop("Diagnostic");
      }
      // write system, loop through write classes
      for (uint i = 0; i < P::systemWriteTimeInterval.size(); i++) {
         if (P::systemWriteTimeInterval[i] >= 0.0 &&
                 P::t >= P::systemWrites[i] * P::systemWriteTimeInterval[i] - DT_EPSILON) {
            
            phiprof::start("write-system");
            logFile << "(IO): Writing spatial cell and reduced system data to disk, tstep = " << P::tstep << " t = " << P::t << endl << writeVerbose;
            const bool writeGhosts = true;
            if( writeGrid(mpiGrid,outputReducer, i, writeGhosts) == false ) {
               cerr << "FAILED TO WRITE GRID AT" << __FILE__ << " " << __LINE__ << endl;
            }
            P::systemWrites[i]++;
            logFile << "(IO): .... done!" << endl << writeVerbose;
            phiprof::stop("write-system");
         }
      }
      
      // Reduce globalflags::bailingOut from all processes
      MPI_Allreduce(&(globalflags::bailingOut), &(doBailout), 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      
      // Write restart data if needed
      int writeRestartNow;
      
      if (myRank == MASTER_RANK) {
         if (P::saveRestartWalltimeInterval >=0.0 && (
               P::saveRestartWalltimeInterval*wallTimeRestartCounter <=  MPI_Wtime()-initialWtime ||
               P::tstep ==P::tstep_max ||
               P::t >= P::t_max ||
               (doBailout > 0 && P::bailout_write_restart)
            )
         ) {
            writeRestartNow = 1;
         }
         else {
            writeRestartNow = 0;
         }
      }
      MPI_Bcast( &writeRestartNow, 1 , MPI_INT , MASTER_RANK ,MPI_COMM_WORLD);
      
      if (writeRestartNow == 1){
         phiprof::start("write-restart");
         wallTimeRestartCounter++;
         
         if (myRank == MASTER_RANK)
            logFile << "(IO): Writing restart data to disk, tstep = " << P::tstep << " t = " << P::t << endl << writeVerbose;
         //Write the restart:
         if( writeRestart(mpiGrid,outputReducer,"restart",(uint)P::t, P::restartStripeFactor) == false ) {
            logFile << "(IO): ERROR Failed to write restart!" << endl << writeVerbose;
            cerr << "FAILED TO WRITE RESTART" << endl;
         }
         if (myRank == MASTER_RANK)
            logFile << "(IO): .... done!"<< endl << writeVerbose;
         phiprof::stop("write-restart");
      }
      
      phiprof::stop("IO");
      addTimedBarrier("barrier-end-io");
      
      //no need to propagate if we are on the final step, we just
      //wanted to make sure all IO is done even for final step
      if(P::tstep ==P::tstep_max ||
         P::t >= P::t_max ||
         doBailout > 0) {
         break;
      }
      
      //Re-loadbalance if needed
      //TODO - add LB measure and do LB if it exceeds threshold
      if( P::tstep%P::rebalanceInterval == 0 && P::tstep> P::tstep_min) {
         logFile << "(LB): Start load balance, tstep = " << P::tstep << " t = " << P::t << endl << writeVerbose;
         balanceLoad(mpiGrid);
         addTimedBarrier("barrier-end-load-balance");
         phiprof::start("Shrink_to_fit");
         // * shrink to fit after LB * //
         shrink_to_fit_grid_data(mpiGrid);
         phiprof::stop("Shrink_to_fit");
         logFile << "(LB): ... done!"  << endl << writeVerbose;
      }
      
      //get local cells
      const vector<CellID>& cells = getLocalCells();

      //compute how many spatial cells we solve for this step
      computedCells=0;
      for(uint i=0;i<cells.size();i++)  computedCells+=mpiGrid[cells[i]]->get_number_of_velocity_blocks()*WID3;
      computedTotalCells+=computedCells;
      
      //Check if dt needs to be changed, and propagate V back a half-step to change dt and set up new situation
      //do not compute new dt on first step (in restarts dt comes from file, otherwise it was initialized before we entered
      //simulation loop
      // FIXME what if dt changes at a restart??
      if(P::dynamicTimestep  && P::tstep> P::tstep_min) {
         computeNewTimeStep(mpiGrid,newDt,dtIsChanged);
         addTimedBarrier("barrier-check-dt");
         if(dtIsChanged) {
            phiprof::start("update-dt");
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
            phiprof::stop("update-dt");
            continue; //
            addTimedBarrier("barrier-new-dt-set");
         }
      }
      
      phiprof::start("Propagate");
      //Propagate the state of simulation forward in time by dt:
      
      phiprof::start("Spatial-space");
      if( P::propagateVlasovTranslation) {
         calculateSpatialTranslation(mpiGrid,P::dt);
      } else {
         calculateSpatialTranslation(mpiGrid,0.0);
      }
      phiprof::stop("Spatial-space",computedCells,"Cells");

      phiprof::start("Compute interp moments");
      calculateInterpolatedVelocityMoments(
         mpiGrid,
         CellParams::RHO_DT2,
         CellParams::RHOVX_DT2,
         CellParams::RHOVY_DT2,
         CellParams::RHOVZ_DT2,
         CellParams::P_11_DT2,
         CellParams::P_22_DT2,
         CellParams::P_33_DT2
      );
      phiprof::stop("Compute interp moments");
      
      if (P::propagateVlasovTranslation || P::propagateVlasovAcceleration ) {
         phiprof::start("Update system boundaries (Vlasov)");
         sysBoundaries.applySysBoundaryVlasovConditions(mpiGrid, P::t+0.5*P::dt); 
         phiprof::stop("Update system boundaries (Vlasov)");
         addTimedBarrier("barrier-boundary-conditions");
      }
      
      // Propagate fields forward in time by dt. This needs to be done before the
      // moments for t + dt are computed (field uses t and t+0.5dt)
      if (P::propagateField) {
         phiprof::start("Propagate Fields");
         propagateFields(mpiGrid, sysBoundaries, P::dt, P::fieldSolverSubcycles);
         phiprof::stop("Propagate Fields",cells.size(),"SpatialCells");
         addTimedBarrier("barrier-after-field-solver");
      }
      
      if (P::propagatePotential == true) {
         poisson::solve(mpiGrid);
      }
      
      phiprof::start("Velocity-space");
      if ( P::propagateVlasovAcceleration ) {
         calculateAcceleration(mpiGrid,P::dt);
         addTimedBarrier("barrier-after-ad just-blocks");
      } else {
         //zero step to set up moments _v
         calculateAcceleration(mpiGrid, 0.0);
      }

      phiprof::stop("Velocity-space",computedCells,"Cells");
      addTimedBarrier("barrier-after-acceleration");

      // *here we compute rho and rho_v for timestep t + dt, so next
      // timestep * //
      calculateInterpolatedVelocityMoments(
         mpiGrid,
         CellParams::RHO,
         CellParams::RHOVX,
         CellParams::RHOVY,
         CellParams::RHOVZ,
         CellParams::P_11,
         CellParams::P_22,
         CellParams::P_33
      );

      phiprof::stop("Propagate",computedCells,"Cells");
      
      // Check timestep
      if(P::dt < P::bailout_min_dt) {
         string message = "The timestep went below bailout.bailout_min_dt (" + to_string(P::bailout_min_dt) + ").";
         bailout(true, message, __FILE__, __LINE__);
      }
      //Move forward in time
      P::meshRepartitioned = false;
      ++P::tstep;
      P::t += P::dt;
   }

   double after = MPI_Wtime();
   
   phiprof::stop("Simulation");
   phiprof::start("Finalization");
   if (P::propagateField ) { 
      finalizeFieldPropagator(mpiGrid);
   }
   if (P::propagatePotential == true) {
      poisson::finalize();
   }
   if (myRank == MASTER_RANK) {
      if(doBailout > 0) {
         logFile << "(BAILOUT): Bailing out, see error log for details." << endl;
      }
      
      double timePerStep;
      if(P::tstep == P::tstep_min) {
         timePerStep=0.0;
      } else {
         timePerStep=double(after  - startTime) / (P::tstep-P::tstep_min);
      }
      double timePerSecond=double(after  - startTime) / (P::t-P::t_min+DT_EPSILON);
      logFile << "(MAIN): All timesteps calculated." << endl;
      logFile << "\t (TIME) total run time " << after - startTime << " s, total simulated time " << P::t -P::t_min<< " s" << endl;
      if(P::t != 0.0) {
         logFile << "\t (TIME) seconds per timestep " << timePerStep  <<
         ", seconds per simulated second " <<  timePerSecond << endl;
      }
      logFile << writeVerbose;
   }
   
   phiprof::stop("Finalization");
   phiprof::stop("main");
   
   phiprof::print(MPI_COMM_WORLD,"phiprof_full");
   phiprof::print(MPI_COMM_WORLD,"phiprof_reduced",0.01);
   phiprof::printLogProfile(MPI_COMM_WORLD,P::tstep,"phiprof_log"," ",7);
   
   if (myRank == MASTER_RANK) logFile << "(MAIN): Exiting." << endl << writeVerbose;
   logFile.close();
   if (P::diagnosticInterval != 0) diagnostic.close();

   MPI_Finalize();
   return 0;
}
