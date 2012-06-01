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

#include "boost/mpi.hpp"
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
#include "datareducer.h"
#include "datareductionoperator.h"
#include "transferstencil.h"

#include "vlsvwriter2.h" 
#include "fieldsolver.h"
#include "project.h"
#include "grid.h"

#ifdef CATCH_FPE
#include <fenv.h>
#include <signal.h>
void fpehandler(int sig_num)
{
   signal(SIGFPE, fpehandler);
   printf("SIGFPE: floating point exception occured, exiting.\n");
   abort();
}
#endif

#include "phiprof.hpp"


Logger logfile, diagnostic;

using namespace std;
using namespace phiprof;


int main(int argn,char* args[]) {
   bool success = true;
   const int MASTER_RANK = 0;
   int myrank;
   typedef Parameters P;
   // Init MPI: 
#ifdef _OPENMP
   //init threaded MPI when compiled using openmp
   int required=MPI_THREAD_FUNNELED;
   int provided;
   MPI_Init_thread(&argn,&args,required,&provided);
   if ( required >provided){
      MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
      if(myrank==MASTER_RANK)
         cerr << "(MAIN): MPI_Init_thread failed!" << endl;
      exit(1);
   }    
#endif
   //Init boost-mpi
   boost::mpi::environment env(argn,args);
   boost::mpi::communicator comm;
   
   MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
   
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
   readparameters.parse();
   P::getParameters();
   getProjectParameters();
   phiprof::stop("Read parameters");

   
   phiprof::start("Init project");
   if (initializeProject() == false) {
      logfile << "(MAIN): Project did not initialize correctly!" << endl << writeVerbose;
      exit(1);
   }
   phiprof::stop("Init project");
   

// Init parallel logger:
   phiprof::start("open logfile & diagnostic");
   if (logfile.open(MPI_COMM_WORLD,MASTER_RANK,"logfile.txt") == false) {
      cerr << "(MAIN) ERROR: Logger failed to open logfile!" << endl;
      exit(1);
   }
   if (P::diagnosticInterval != 0) {
      if (diagnostic.open(MPI_COMM_WORLD,MASTER_RANK,"diagnostic.txt") == false) {
         cerr << "(MAIN) ERROR: Logger failed to open diagnostic file!" << endl;
         exit(1);
      }
   }
   phiprof::stop("open logfile & diagnostic");
   
   phiprof::start("Init grid");
   dccrg::Dccrg<SpatialCell> mpiGrid;
   if (initializeGrid(argn,args,mpiGrid) == false) {
      logfile << "(MAIN): Grid did not initialize correctly!" << endl << writeVerbose;
      exit(1);
   }  
   phiprof::stop("Init grid");

   // Initialize data reduction operators. This should be done elsewhere in order to initialize 
   // user-defined operators:
   DataReducer outputReducer, diagnosticReducer;
   initializeDataReducers(&outputReducer, &diagnosticReducer);
   
   //VlsWriter vlsWriter;
   phiprof::start("Init vlasov propagator");
   // Initialize Vlasov propagator:
   if (initializeMover(mpiGrid) == false) {
      logfile << "(MAIN): Vlasov propagator did not initialize correctly!" << endl << writeVerbose;
      exit(1);
   }   
   calculateVelocityMoments(mpiGrid);
   phiprof::stop("Init vlasov propagator");
   
   phiprof::start("Init field propagator");
   // Initialize field propagator:
   if (initializeFieldPropagator(mpiGrid,P::propagateField) == false) {
       logfile << "(MAIN): Field propagator did not initialize correctly!" << endl << writeVerbose;
       exit(1);
   }
   phiprof::stop("Init field propagator");
   




   
   // ***********************************
   // ***** INITIALIZATION COMPLETE *****
   // ***********************************
   
   // Free up memory:
   readparameters.finalize();

   phiprof::start("Save initial state");
   // Write initial state:
   if (P::save_spatial_grid) {
      if (myrank == MASTER_RANK) {
	 logfile << "(MAIN): Saving initial state of variables to disk." << endl << writeVerbose;
      }

   if (writeGrid(mpiGrid,outputReducer,true) == false) {
	 logfile << "(MAIN): ERROR occurred while writing spatial cell and restart data!" << endl << writeVerbose;
      }
   }
   phiprof::stop("Save initial state");
   
   if (P::diagnosticInterval != 0)
   {
      phiprof::start("Diagnostic");
      if (computeDiagnostic(mpiGrid, diagnosticReducer, P::tstep_min) == false) {
	 cerr << "ERROR with diagnostic computation" << endl;
      }
      phiprof::stop("Diagnostic");
   }
   
   phiprof::stop("Initialization");
   MPI_Barrier(MPI_COMM_WORLD);


   // Main simulation loop:
   if (myrank == MASTER_RANK) logfile << "(MAIN): Starting main simulation loop." << endl << writeVerbose;

   double before = MPI_Wtime();
   unsigned int totalComputedSpatialCells=0;
   unsigned int computedSpatialCells=0;
   unsigned int totalComputedBlocks=0;
   unsigned int computedBlocks=0;

   phiprof::start("Simulation");
   for (luint tstep=P::tstep_min; tstep < P::tsteps; ++tstep) {
      //compute how many spatial cells we solve for this step
      vector<uint64_t> cells = mpiGrid.get_cells();
      computedSpatialCells=cells.size();
      computedBlocks=0;
      for(uint i=0;i<cells.size();i++)
         computedBlocks+=mpiGrid[cells[i]]->number_of_blocks;
      
         
      totalComputedSpatialCells+=computedSpatialCells;
      totalComputedBlocks+=computedBlocks;
      
      phiprof::start("Propagate");
      // Recalculate (maybe) spatial cell parameters
      calculateSimParameters(mpiGrid, P::t, P::dt);
      
      // use globally minimum timestep
      P::dt = all_reduce(comm, P::dt, boost::mpi::minimum<Real>());
      // Propagate the state of simulation forward in time by dt:      
      if (P::propagateVlasov == true) {
         phiprof::start("Propagate Vlasov");
         switch (P::splitMethod){
              case 0:
                 phiprof::start("propagation");
                 calculateSpatialFluxes(mpiGrid,P::dt);
#ifdef SEMILAG
                 adjustVelocityBlocks(mpiGrid);
#endif
                 calculateSpatialPropagation(mpiGrid,true,P::dt);
                 phiprof::stop("propagation",computedBlocks,"Blocks");
                 break;
              case 1:
                 phiprof::start("First propagation");
                 calculateSpatialFluxes(mpiGrid,0.5*P::dt);
                 calculateSpatialPropagation(mpiGrid,true,P::dt);
                 phiprof::stop("First propagation",computedBlocks,"Blocks");
#ifdef SEMILAG
                 adjustVelocityBlocks(mpiGrid);
#endif
                 phiprof::start("Second propagation");
                 calculateSpatialFluxes(mpiGrid,0.5*P::dt);
                 calculateSpatialPropagation(mpiGrid,false,0);
                 phiprof::stop("Second propagation",computedBlocks,"Blocks");
                 break;
              case 2:
                 phiprof::start("Acceleration");
                 calculateAcceleration(mpiGrid,0.5*P::dt);
                 phiprof::stop("Acceleration",computedBlocks,"Blocks");
#ifdef SEMILAG
                 adjustVelocityBlocks(mpiGrid);
#endif
                 phiprof::start("Trans + acc");
                 calculateSpatialFluxes(mpiGrid,P::dt);
                 calculateSpatialPropagation(mpiGrid,true,0.5*P::dt);
                 phiprof::stop("Trans + acc",computedBlocks,"Blocks");
#ifdef SEMILAG
                 //if we are careful this might not be needed, if the
                 //next timestep acceleration is the next step, then
                 //it does not matter what other cells have done
                 adjustVelocityBlocks(mpiGrid);
#endif
                 break;

         }

#ifndef SEMILAG
         //if no semi-lagrangean then we do the adjustment here.
         if(P::tstep%P::blockAdjustmentInterval == 0){
            adjustVelocityBlocks(mpiGrid);
         }
#endif        
	  
         phiprof::stop("Propagate Vlasov",computedBlocks,"Blocks");
      }
      
      // Propagate fields forward in time by dt. If field is not 
      // propagated self-consistently (test-Vlasov simulation), then 
      // re-calculate face-averaged E,B fields. This requires that 
      // edge-E and face-B have been shared with remote neighbours 
      // (not done by calculateFaceAveragedFields).

      if (P::propagateField == true) {
          phiprof::start("Propagate Fields");
          propagateFields(mpiGrid,P::dt);
          phiprof::stop("Propagate Fields",computedSpatialCells,"SpatialCells");
      } else {
	 calculateFaceAveragedFields(mpiGrid);
      }
      phiprof::stop("Propagate",computedBlocks,"Blocks");
      
      ++P::tstep;
      P::t += P::dt;
      
      // Check if data needs to be written to disk:
      if (P::tstep % P::saveRestartInterval == 0 || P::tstep % P::saveInterval == 0) {
         phiprof::start("IO");
	 bool writeRestartData = false;
	 if (P::tstep % P::saveRestartInterval == 0) {
	   writeRestartData = true;
	   if (myrank == MASTER_RANK)
	   logfile << "(MAIN): Writing spatial cell and restart data to disk, tstep = " << P::tstep << " t = " << P::t << endl << writeVerbose;
	 } else
	   if (myrank == MASTER_RANK)
	     logfile << "(MAIN): Writing spatial cell data to disk, tstep = " << P::tstep << " t = " << P::t << endl << writeVerbose;
	 
	   if (writeGrid(mpiGrid,outputReducer,writeRestartData) == false) {
	    if (myrank == MASTER_RANK)
	      logfile << "(MAIN): ERROR occurred while writing spatial cell and restart data!" << endl << writeVerbose;
	 }
         phiprof::stop("IO");
      }
      
      // Check whether diagnostic output has to be produced
      if (P::diagnosticInterval != 0 && P::tstep % P::diagnosticInterval == 0) {
	 phiprof::start("Diagnostic");
	 if (computeDiagnostic(mpiGrid, diagnosticReducer, tstep+1) == false) {
	    cerr << "ERROR with diagnostic computation" << endl;
	 }
	 phiprof::stop("Diagnostic");
         //also print out phiprof log
         phiprof::printLogProfile(MPI_COMM_WORLD,P::tstep,"phiprof_log"," ",7);
      }

      
      //compute maximum time-step
      phiprof::start("compute-timestep");
      Real dtmax_local[3];
      Real dtmax_global[3];
      dtmax_local[0]=std::numeric_limits<Real>::max();
      dtmax_local[1]=std::numeric_limits<Real>::max();
      dtmax_local[2]=std::numeric_limits<Real>::max();
      for (std::vector<uint64_t>::const_iterator
              cell_id = cells.begin();
           cell_id != cells.end();
           ++cell_id
           ) {
         SpatialCell* cell = mpiGrid[*cell_id];
         dtmax_local[0]=min(dtmax_local[0],cell->parameters[CellParams::MAXRDT]);
         dtmax_local[1]=min(dtmax_local[1],cell->parameters[CellParams::MAXVDT]);
         dtmax_local[2]=min(dtmax_local[2],cell->parameters[CellParams::MAXFDT]);
      }
      MPI_Allreduce(&(dtmax_local[0]),&(dtmax_global[0]),3,MPI_Type<Real>(), MPI_MIN, MPI_COMM_WORLD);
      if (myrank == MASTER_RANK)
         logfile << "(MAIN) tstep = " << P::tstep << " dt = " << P::dt  <<
            " max timestep in (r,v,BE) is "<< dtmax_global[0] <<" " <<dtmax_global[1] <<" " << dtmax_global[2]  << endl << writeVerbose;
      
      if(P::dynamicTimestep){
         Real dtmax=std::numeric_limits<Real>::max();
         switch (P::splitMethod){
             case 0:
                dtmax=min(dtmax,dtmax_global[0]);
                dtmax=min(dtmax,dtmax_global[1]);
                dtmax=min(dtmax,dtmax_global[2]);
                P::dt=P::CFL*dtmax;
                break;
             case 1:
                dtmax=min(dtmax,2*dtmax_global[0]); //half-steps in ordinary space
                dtmax=min(dtmax,dtmax_global[1]);
                dtmax=min(dtmax,dtmax_global[2]);
                P::dt=P::CFL*dtmax;
                break;
             case 2:
                dtmax=min(dtmax,dtmax_global[0]);
                dtmax=min(dtmax,2*dtmax_global[1]);  //half-steps in velocity space
                dtmax=min(dtmax,dtmax_global[2]);
                P::dt=P::CFL*dtmax;
                break;
         }
      }
      
      phiprof::stop("compute-timestep");
      
      //Last-but-not-least, re-loadbalance 
      if(P::tstep%P::rebalanceInterval == 0 && P::tstep> P::tstep_min)
         balanceLoad(mpiGrid);
         

   }
   double after = MPI_Wtime();

   phiprof::stop("Simulation",totalComputedBlocks,"Blocks");
   phiprof::start("Finalization");   
   finalizeMover();
   finalizeFieldPropagator(mpiGrid);
   
   if (myrank == MASTER_RANK) {
      logfile << "(MAIN): All timesteps calculated." << endl;
      logfile << "\t (TIME) total run time " << after - before << " s, total simulated time " << P::t << " s" << endl;
      if(P::t != 0.0) {
	 logfile << "\t (TIME) seconds per timestep " << double(after - before) / P::tsteps <<
	 ", seconds per simulated second " << double(after - before) / P::t << endl;
      }
      logfile << writeVerbose;
   }
   
   phiprof::stop("Finalization");   
   phiprof::stop("main");
   
   phiprof::print(MPI_COMM_WORLD,"phiprof_full");
   phiprof::print(MPI_COMM_WORLD,"phiprof_reduced",0.01);
   
   if (myrank == MASTER_RANK) logfile << "(MAIN): Exiting." << endl << writeVerbose;
   logfile.close();
   if (P::diagnosticInterval != 0) diagnostic.close();
   return 0;
}

