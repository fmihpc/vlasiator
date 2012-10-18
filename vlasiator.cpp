/*
This file is part of Vlasiator.

Copyright 2010, 2011, 2012 Finnish Meteorological Institute

Vlasiator is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3
as published by the Free Software Foundation.

Vlasiator is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <cstdlib>
#include <iostream>
#include <cmath>
#include <vector>
#include <sstream>
#include <ctime>

#include "vlasovmover.h"
#include "definitions.h"
#include "mpiconversion.h"
#include "logger.h"
#include "parameters.h"
#include "readparameters.h"
#include "spatial_cell.hpp"
#include "datareduction/datareducer.h"
#include "sysboundary/sysboundary.h"
#include "transferstencil.h"

#include "vlsvwriter2.h" 
#include "fieldsolver.h"
#include "project.h"
#include "grid.h"
#include "fileio.h"

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

using namespace std;
using namespace phiprof;


bool changeTimeStep(dccrg::Dccrg<SpatialCell>& mpiGrid,Real &newDt, bool &isChanged) {

   phiprof::start("compute-timestep");      
   //compute maximum time-step, this cannot be done at the first
   //step as the solvers compute the limits for each cell

   isChanged=false;

   vector<uint64_t> cells = mpiGrid.get_cells();         
   Real dtMaxLocal[3];
   Real dtMaxGlobal[3];
   dtMaxLocal[0]=std::numeric_limits<Real>::max();
   dtMaxLocal[1]=std::numeric_limits<Real>::max();
   dtMaxLocal[2]=std::numeric_limits<Real>::max();
   for (std::vector<uint64_t>::const_iterator cell_id = cells.begin(); cell_id != cells.end(); ++cell_id) {
      SpatialCell* cell = mpiGrid[*cell_id];
      if (cell->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY) continue;
      dtMaxLocal[0]=min(dtMaxLocal[0], cell->parameters[CellParams::MAXRDT]);
      dtMaxLocal[1]=min(dtMaxLocal[1], cell->parameters[CellParams::MAXVDT]);
      dtMaxLocal[2]=min(dtMaxLocal[2], cell->parameters[CellParams::MAXFDT]);
   }
   MPI_Allreduce(&(dtMaxLocal[0]), &(dtMaxGlobal[0]), 3, MPI_Type<Real>(), MPI_MIN, MPI_COMM_WORLD);
   


      
   Real maxVDtNoSubstepping=dtMaxGlobal[1];
   
   //Increase timestep limit in velocity space based on
   //maximum number of substeps we are allowed to take. As
   //the length of each substep is unknown beforehand the
   //limit is not hard, it may be exceeded in some cases.
   // A value of 0 means that there is no limit on substepping
   if(P::maxAccelerationSubsteps==0)
      dtMaxGlobal[1]=std::numeric_limits<Real>::max();
   else
      dtMaxGlobal[1]*=P::maxAccelerationSubsteps;
   
   //If fieldsolver is off there should be no limit on time-step from it
   if (P::propagateField == false) {
      dtMaxGlobal[2]=std::numeric_limits<Real>::max();
   }
   
   //If vlasovsolver is off there should be no limit on time-step from it
   if (P::propagateVlasov == false) {
      dtMaxGlobal[0]=std::numeric_limits<Real>::max();
      dtMaxGlobal[1]=std::numeric_limits<Real>::max();
      maxVDtNoSubstepping=std::numeric_limits<Real>::max();
   }
   
   Real dtMax=std::numeric_limits<Real>::max();
   dtMax=min(dtMax, dtMaxGlobal[0]);
   dtMax=min(dtMax, dtMaxGlobal[1]); 
   dtMax=min(dtMax, dtMaxGlobal[2]);
   
   
   if(P::dt < dtMax*P::CFL_min || P::dt > dtMax*P::CFL_max ) {
      newDt = 0.5*(P::CFL_min + P::CFL_max)*dtMax;
      isChanged=true;
      logFile <<"(TIMESTEP) New dt computed: "<< newDt <<
         " Max dt (not including CFL "<< P::CFL_min <<"-"<<P::CFL_max<<" ) in {r, v+subs, v, BE} was " <<
         dtMaxGlobal[0] << " " <<
         dtMaxGlobal[1] << " " <<
         maxVDtNoSubstepping << " " <<
         dtMaxGlobal[2] << endl << writeVerbose;
   }

   phiprof::stop("compute-timestep");

   return true;
}


int main(int argn,char* args[]) {
   bool success = true;
   int myRank;
   const creal DT_EPSILON=1e-12;
   typedef Parameters P;
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
   MPI_Comm comm = MPI_COMM_WORLD;
   MPI_Comm_rank(comm,&myRank);
   dccrg::Dccrg<SpatialCell> mpiGrid;
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
   addProjectParameters();
   sysBoundaries.addParameters();
   readparameters.parse();
   P::getParameters();
   getProjectParameters();
   sysBoundaries.getParameters();
   phiprof::stop("Read parameters");
   
   phiprof::start("Init project");
   if (initializeProject() == false) {
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
   }
   if (P::diagnosticInterval != 0) {
      if (diagnostic.open(MPI_COMM_WORLD,MASTER_RANK,"diagnostic.txt",P::isRestart) == false) {
         if(myRank == MASTER_RANK) cerr << "(MAIN) ERROR: Logger failed to open diagnostic file!" << endl;
         exit(1);
      }
   }
   phiprof::stop("open logFile & diagnostic");
   
   phiprof::start("Init grid");
   initializeGrid(argn,args,mpiGrid,sysBoundaries);
   isSysBoundaryCondDynamic = sysBoundaries.isDynamic();
   phiprof::stop("Init grid");
   
   phiprof::start("Init DROs");
   // Initialize data reduction operators. This should be done elsewhere in order to initialize 
   // user-defined operators:
   DataReducer outputReducer, diagnosticReducer;
   initializeDataReducers(&outputReducer, &diagnosticReducer);
   phiprof::stop("Init DROs");
   
   //VlsWriter vlsWriter;
   phiprof::start("Init vlasov propagator");
   // Initialize Vlasov propagator:
   if (initializeMover(mpiGrid) == false) {
      logFile << "(MAIN): Vlasov propagator did not initialize correctly!" << endl << writeVerbose;
      exit(1);
   }
   //compute moments, and set them in RHO* and RHO*_DT2
   calculateVelocityMoments(mpiGrid,false);
   phiprof::stop("Init vlasov propagator");
   
   phiprof::start("Init field propagator");
   // Initialize field propagator:
   if (initializeFieldPropagator(mpiGrid) == false) {
       logFile << "(MAIN): Field propagator did not initialize correctly!" << endl << writeVerbose;
       exit(1);
   }
   phiprof::stop("Init field propagator");
   
   

   
   // Free up memory:
   readparameters.finalize();
   
   
   // ***********************************
   // ***** INITIALIZATION COMPLETE *****
   // ***********************************
   

   
   
   // Main simulation loop:
   if (myRank == MASTER_RANK) logFile << "(MAIN): Starting main simulation loop." << endl << writeVerbose;
   
   double before = MPI_Wtime();
   unsigned int totalComputedSpatialCells=0;
   unsigned int computedSpatialCells=0;
   unsigned int totalComputedBlocks=0;
   unsigned int computedBlocks=0;
   unsigned int restartWrites=(int)(P::t_min/P::saveRestartTimeInterval);
   unsigned int systemWrites=(int)(P::t_min/P::saveSystemTimeInterval);
   Real newDt;
   bool dtIsChanged;
         
   


   if(P::dynamicTimestep && !P::isRestart) {
      //compute vlasovsolver once with zero dt, this is  to initialize
      //per-cell dt limits. In restarts, we read in dt from file
      phiprof::start("compute-dt");
      if(P::propagateVlasov) {
         calculateSpatialFluxes(mpiGrid,0.0);
         calculateSpatialPropagation(mpiGrid,true,0.0);
      }
      if(P::propagateField) {
         propagateFields(mpiGrid,0.0);
      }
      //compute new dt
      changeTimeStep(mpiGrid,newDt,dtIsChanged);
      P::dt=newDt;
      phiprof::stop("compute-dt");
      
   }

   
   if(P::propagateVlasov && !P::isRestart) {
      //go forward by dt/2 in x, initializes leapfrog split. In restarts the
      //the distribution function is already propagated forward in time by dt/2
      phiprof::start("propagate-spatial-space-dt/2");
      calculateSpatialFluxes(mpiGrid,0.5*P::dt);
      calculateSpatialPropagation(mpiGrid,false,0.0);
      phiprof::stop("propagate-spatial-space-dt/2");
   }
   
   phiprof::stop("Initialization");
   MPI_Barrier(MPI_COMM_WORLD);   

   phiprof::start("Simulation");

   while(P::tstep <=P::tstep_max  &&
         P::t-P::dt <= P::t_max+DT_EPSILON) {

      phiprof::start("IO");
      
      //write out phiprof profiles and logs with a lower interval than normal
      //diagnostic (every 10 diagnostic intervals). Also, do not print
      //first step or until we have rebalanced.
      logFile << "------------------ tstep = " << P::tstep << " t = " << P::t <<" dt = " << P::dt << " ------------------" << endl;
      if (P::diagnosticInterval != 0 &&
          P::tstep % (P::diagnosticInterval*10) == 0 &&
          P::tstep-P::tstep_min >0  &&
          P::tstep-P::tstep_min > P::rebalanceInterval) {
         phiprof::print(MPI_COMM_WORLD,"phiprof_reduced",0.01);
         phiprof::printLogProfile(MPI_COMM_WORLD,P::tstep,"phiprof_log"," ",7);
         
         double currentTime=MPI_Wtime();
         double timePerStep=double(currentTime  - before) / (P::tstep-P::tstep_min);
         double timePerSecond=double(currentTime  - before) / (P::t-P::t_min + DT_EPSILON);
         double remainingTime=min(timePerStep*(P::tstep_max-P::tstep),timePerSecond*(P::t_max-P::t));
         time_t finalWallTime=time(NULL)+(time_t)remainingTime; //assume time_t is in seconds, as it is almost always
         struct tm *finalWallTimeInfo=localtime(&finalWallTime);
         logFile << "(TIME) walltime/step " << timePerStep<< " s" <<endl;
         logFile << "(TIME) walltime/simusecond (s)" << timePerSecond<<" s" <<endl;
         logFile << "(TIME) Estimated completion time is " <<asctime(finalWallTimeInfo)<<endl;
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
      

      
      
      // Save reduced data
      if ( P::saveSystemTimeInterval >=0.0 && 
           P::t >= systemWrites*P::saveSystemTimeInterval-DT_EPSILON ){
         
         phiprof::start("write-system");
         if (myRank == MASTER_RANK)
            logFile << "(IO): Writing spatial cell and reduced system data to disk, tstep = " << P::tstep << " t = " << P::t << endl << writeVerbose;
         writeGrid(mpiGrid,outputReducer,"grid",systemWrites,false);
         systemWrites++;
         phiprof::stop("write-system");
      }

      // Save restart data        
      if (P::saveRestartTimeInterval >=0.0 &&
          P::t >= restartWrites*P::saveRestartTimeInterval-DT_EPSILON){
         phiprof::start("write-restart");
         if (myRank == MASTER_RANK)
            logFile << "(IO): Writing spatial cell and restart data to disk, tstep = " << P::tstep << " t = " << P::t << endl << writeVerbose;
         writeGrid(mpiGrid,outputReducer,"restart",restartWrites,true);
         restartWrites++;
         phiprof::stop("write-restart");
      }
      
      
      phiprof::stop("IO");
      
      
      //no need to propagate if we are on the final step, we just
      //wanted to make sure all IO is done even for final step
      if(P::tstep ==P::tstep_max ||
         P::t >= P::t_max) {
         break;
      }

            
      //Re-loadbalance if needed, not done on first step
      if(P::tstep%P::rebalanceInterval == 0 && P::tstep> P::tstep_min)
         balanceLoad(mpiGrid);
      
      //get local cells       
      vector<uint64_t> cells = mpiGrid.get_cells();      
      //compute how many spatial cells we solve for this step
      computedSpatialCells=cells.size();
      computedBlocks=0;
      for(uint i=0;i<cells.size();i++)  computedBlocks+=mpiGrid[cells[i]]->number_of_blocks;
      totalComputedSpatialCells+=computedSpatialCells;
      totalComputedBlocks+=computedBlocks;
      


            
      //Check if dt needs to be changed, and propagate half-steps properly to change dt and set up new situation
      //do not compute new dt on first step (in restarts dt comes from file, otherwise it was initialized before we entered
      //simulation loop
      if(P::dynamicTimestep  && P::tstep> P::tstep_min) {
         changeTimeStep(mpiGrid,newDt,dtIsChanged);         
         if(dtIsChanged) {
            phiprof::start("update-dt");
            //propagate velocity space to real-time, do not do it if dt is zero
            if(P::dt >0)
               calculateAcceleration(mpiGrid,0.5*P::dt);
            //re-compute rho for real time for fieldsolver, and compute _DT2 values as an average of old and new values
            calculateVelocityMoments(mpiGrid,true);
            
            // Propagate fields forward in time by 0.5*dt
            if (P::propagateField == true) {
               phiprof::start("Propagate Fields");
               propagateFields(mpiGrid,0.5*P::dt);
               phiprof::stop("Propagate Fields",computedSpatialCells,"SpatialCells");
            } else {
               // TODO Whatever field updating/volume averaging/etc. needed in test particle and other test cases have to be put here.
               // In doing this be sure the needed components have been updated.
            }
            ++P::tstep;
            P::t += P::dt*0.5;
            P::dt=newDt;
            
            //go forward by dt/2 again in x
            calculateSpatialFluxes(mpiGrid,0.5*P::dt);
            calculateSpatialPropagation(mpiGrid,false,0.0);
            logFile <<" dt changed to "<<P::dt <<"s, distribution function was half-stepped to real-time and back"<<endl<<writeVerbose;
            phiprof::stop("update-dt");
            continue; //
         }
      }
      
      
      bool updateVelocityBlocksAfterAcceleration=false;
#ifdef SEMILAG
      updateVelocityBlocksAfterAcceleration=true;
#endif
      if(P::maxAccelerationSubsteps!=1)
         updateVelocityBlocksAfterAcceleration=true;
      
      phiprof::start("Propagate");
      //Propagate the state of simulation forward in time by dt:      
      if (P::propagateVlasov == true) {
         phiprof::start("Propagate Vlasov");

         phiprof::start("Velocity-space");
         calculateAcceleration(mpiGrid,P::dt);
         phiprof::stop("Velocity-space",computedBlocks,"Blocks");
         if(updateVelocityBlocksAfterAcceleration){
            //need to do a update of block lists as all cells have made local changes
            updateRemoteVelocityBlockLists(mpiGrid);
            adjustVelocityBlocks(mpiGrid);
         }         
         calculateInterpolatedVelocityMoments(mpiGrid,CellParams::RHO_DT2,CellParams::RHOVX_DT2,CellParams::RHOVY_DT2,CellParams::RHOVZ_DT2);

         phiprof::start("Spatial-space");
         calculateSpatialFluxes(mpiGrid,P::dt);
         calculateSpatialPropagation(mpiGrid,false,0.0);
         phiprof::stop("Spatial-space",computedBlocks,"Blocks");

         calculateInterpolatedVelocityMoments(mpiGrid,CellParams::RHO,CellParams::RHOVX,CellParams::RHOVY,CellParams::RHOVZ);
         
         if(!updateVelocityBlocksAfterAcceleration){
            //if no semi-lagrangean or substepping in leveque   
            //acceleration then we do the adjustment here. In that
            //case no local block adjustments have been done, so
            //remote velocty blocks do not need to be updated
            if(P::tstep%P::blockAdjustmentInterval == 0){
               adjustVelocityBlocks(mpiGrid);
            }
         }
         phiprof::stop("Propagate Vlasov",computedBlocks,"Blocks");
      }
      
      // Propagate fields forward in time by dt.
      if (P::propagateField == true) {
         phiprof::start("Propagate Fields");
         propagateFields(mpiGrid,P::dt);
         phiprof::stop("Propagate Fields",computedSpatialCells,"SpatialCells");
      } else {
         // TODO Whatever field updating/volume averaging/etc. needed in test particle and other test cases have to be put here.
         // In doing this be sure the needed components have been updated.
      }
      phiprof::stop("Propagate",computedBlocks,"Blocks");
      
      //Move forward in time      
      ++P::tstep;
      P::t += P::dt;
   }

   double after = MPI_Wtime();
   
   phiprof::stop("Simulation",totalComputedBlocks,"Blocks");
   phiprof::start("Finalization");   
   finalizeMover();
   finalizeFieldPropagator(mpiGrid);
   
   if (myRank == MASTER_RANK) {
      double timePerStep;
      if(P::tstep == P::tstep_min) {
         timePerStep=0.0;
      } else {
         timePerStep=double(after  - before) / (P::tstep-P::tstep_min);	
      }
      double timePerSecond=double(after  - before) / (P::t-P::t_min+DT_EPSILON);
      logFile << "(MAIN): All timesteps calculated." << endl;
      logFile << "\t (TIME) total run time " << after - before << " s, total simulated time " << P::t -P::t_min<< " s" << endl;
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
   
   if (myRank == MASTER_RANK) logFile << "(MAIN): Exiting." << endl << writeVerbose;
   logFile.close();
   if (P::diagnosticInterval != 0) diagnostic.close();

   MPI_Finalize();
   return 0;
}
