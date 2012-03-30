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
#include "mpilogger.h"
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

#include "profile.hpp"


MPILogger mpilogger;

using namespace std;
using namespace profile;

   



int main(int argn,char* args[]) {
   bool success = true;
   const int MASTER_RANK = 0;
   int myrank;
   typedef Parameters P;
   // Init MPI: 
#ifdef _OPENMP
   //init threaded MPI when comppiled using openmp
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
   
   profile::start("main");
   profile::start("Initialization");
   profile::start("Read parameters");
   //init parameter file reader
   Readparameters readparameters(argn,args,MPI_COMM_WORLD);
   P::addParameters();
   addProjectParameters();
   readparameters.parse();
   P::getParameters();
   getProjectParameters();
   profile::stop("Read parameters");


// Init parallel logger:
   profile::start("open mpilogger");
   if (mpilogger.open(MPI_COMM_WORLD,MASTER_RANK,"logfile.txt") == false) {
      cerr << "(MAIN) ERROR: MPILogger failed to open output file!" << endl;
      exit(1);
   }
   profile::stop("open mpilogger");
   
   profile::start("Init grid");
   dccrg::Dccrg<SpatialCell> mpiGrid;
   if (initializeGrid(argn,args,mpiGrid) == false) {
      mpilogger << "(MAIN): Grid did not initialize correctly!" << endl << write;
      exit(1);
   }  
   profile::stop("Init grid");

   // Initialize data reduction operators. This should be done elsewhere in order to initialize 
   // user-defined operators:
   DataReducer reducer;
   reducer.addOperator(new DRO::VariableB);
   reducer.addOperator(new DRO::VariableE);
   reducer.addOperator(new DRO::VariableRho);
   reducer.addOperator(new DRO::VariableRhoV);
   reducer.addOperator(new DRO::MPIrank);
   reducer.addOperator(new DRO::Blocks);
   reducer.addOperator(new DRO::VariableVolE);
   reducer.addOperator(new DRO::VariableVolB);
   reducer.addOperator(new DRO::VariablePressure);
      
   //VlsWriter vlsWriter;
   profile::start("Init vlasov propagator");
   // Initialize Vlasov propagator:
   if (initializeMover(mpiGrid) == false) {
      mpilogger << "(MAIN): Vlasov propagator did not initialize correctly!" << endl << write;
      exit(1);
   }   
   calculateVelocityMoments(mpiGrid);
   profile::stop("Init vlasov propagator");
   
   profile::start("Init field propagator");
   // Initialize field propagator:
   if (initializeFieldPropagator(mpiGrid,P::propagateField) == false) {
       mpilogger << "(MAIN): Field propagator did not initialize correctly!" << endl << write;
       exit(1);
   }
   profile::stop("Init field propagator");
   


   profile::start("Init project");
   if (initializeProject() == false) {
      mpilogger << "(MAIN): Project did not initialize correctly!" << endl << write;
      exit(1);
   }
   profile::stop("Init project");


   
   // ***********************************
   // ***** INITIALIZATION COMPLETE *****
   // ***********************************
   
   // Free up memory:
   readparameters.finalize();

   profile::start("Save initial state");
   // Write initial state:
   if (P::save_spatial_grid) {
      if (myrank == MASTER_RANK) {
	 mpilogger << "(MAIN): Saving initial state of variables to disk." << endl << write;
      }

      if (writeGrid(mpiGrid,reducer,true) == false) {
	 mpilogger << "(MAIN): ERROR occurred while writing spatial cell and restart data!" << endl << write;
      }
   }
   profile::stop("Save initial state");
   profile::stop("Initialization");
   MPI_Barrier(MPI_COMM_WORLD);


   // Main simulation loop:
   if (myrank == MASTER_RANK) mpilogger << "(MAIN): Starting main simulation loop." << endl << write;

   double before = MPI_Wtime();
   unsigned int totalComputedSpatialCells=0;
   unsigned int computedSpatialCells=0;
   unsigned int totalComputedBlocks=0;
   unsigned int computedBlocks=0;

   profile::start("Simulation");
   for (luint tstep=P::tstep_min; tstep < P::tsteps; ++tstep) {

      if (myrank == MASTER_RANK)  cout << "On step " << tstep << endl;
      //compute how many spatial cells we solve for this step
      vector<uint64_t> cells = mpiGrid.get_cells();
      computedSpatialCells=cells.size();
      computedBlocks=0;
      for(uint i=0;i<cells.size();i++)
         computedBlocks+=mpiGrid[cells[i]]->number_of_blocks;
      
         
      totalComputedSpatialCells+=computedSpatialCells;
      totalComputedBlocks+=computedBlocks;
      
      profile::start("Propagate");
      // Recalculate (maybe) spatial cell parameters
      calculateSimParameters(mpiGrid, P::t, P::dt);
      
      // use globally minimum timestep
      P::dt = all_reduce(comm, P::dt, boost::mpi::minimum<Real>());
      // Propagate the state of simulation forward in time by dt:      
      if (P::propagateVlasov == true) {
         profile::start("Propagate Vlasov");
         switch (P::splitMethod){
              case 0:
                 profile::start("propagation");
                 calculateSpatialFluxes(mpiGrid,P::dt);
#ifdef SEMILAG
                 adjustVelocityBlocks(mpiGrid);
#endif
                 calculateSpatialPropagation(mpiGrid,true,P::dt);
                 profile::stop("propagation",computedBlocks,"Blocks");
                 break;
              case 1:
                 profile::start("First propagation");
                 calculateSpatialFluxes(mpiGrid,0.5*P::dt);
                 calculateSpatialPropagation(mpiGrid,true,P::dt);
                 profile::stop("First propagation",computedBlocks,"Blocks");
#ifdef SEMILAG
                 adjustVelocityBlocks(mpiGrid);
#endif
                 profile::start("Second propagation");
                 calculateSpatialFluxes(mpiGrid,0.5*P::dt);
                 calculateSpatialPropagation(mpiGrid,false,0);
                 profile::stop("Second propagation",computedBlocks,"Blocks");
                 break;
              case 2:
                 profile::start("Acceleration");
                 calculateAcceleration(mpiGrid,0.5*P::dt);
                 profile::stop("Acceleration",computedBlocks,"Blocks");
#ifdef SEMILAG
                 adjustVelocityBlocks(mpiGrid);
#endif
                 profile::start("Trans + acc");
                 calculateSpatialFluxes(mpiGrid,P::dt);
                 calculateSpatialPropagation(mpiGrid,true,0.5*P::dt);
                 profile::stop("Trans + acc",computedBlocks,"Blocks");
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
	  
          if(P::tstep%P::rebalanceInterval == 0 && P::tstep> P::tstep_min)
             balanceLoad(mpiGrid);
          
          profile::stop("Propagate Vlasov",computedBlocks,"Blocks");
      }
 
      // Propagate fields forward in time by dt. If field is not 
      // propagated self-consistently (test-Vlasov simulation), then 
      // re-calculate face-averaged E,B fields. This requires that 
      // edge-E and face-B have been shared with remote neighbours 
      // (not done by calculateFaceAveragedFields).

      if (P::propagateField == true) {
          profile::start("Propagate Fields");
          propagateFields(mpiGrid,P::dt);
          profile::stop("Propagate Fields",computedSpatialCells,"SpatialCells");
      } else {
	 calculateFaceAveragedFields(mpiGrid);
      }
      profile::stop("Propagate",computedBlocks,"Blocks");

      ++P::tstep;
      P::t += P::dt;
      
      // Check if data needs to be written to disk:
      if (P::tstep % P::saveRestartInterval == 0 || P::tstep % P::diagnInterval == 0) {
         profile::start("IO");
	 bool writeRestartData = false;
	 if (P::tstep % P::saveRestartInterval == 0) {
	   writeRestartData = true;
	   if (myrank == MASTER_RANK)
	   mpilogger << "(MAIN): Writing spatial cell and restart data to disk, tstep = " << P::tstep << " t = " << P::t << endl << write;
	 } else
	   if (myrank == MASTER_RANK)
	     mpilogger << "(MAIN): Writing spatial cell data to disk, tstep = " << P::tstep << " t = " << P::t << endl << write;
	 
	 if (writeGrid(mpiGrid,reducer,writeRestartData) == false) {
	    if (myrank == MASTER_RANK)
	      mpilogger << "(MAIN): ERROR occurred while writing spatial cell and restart data!" << endl << write;
	 }
         profile::stop("IO");
      }
   }
   double after = MPI_Wtime();

   profile::stop("Simulation",totalComputedBlocks,"Blocks");
   profile::start("Finalization");   
   finalizeMover();
   finalizeFieldPropagator(mpiGrid);
   
   if (myrank == MASTER_RANK) {
      mpilogger << "(MAIN): All timesteps calculated." << endl;
      mpilogger << "\t (TIME) total run time " << after - before << " s, total simulated time " << P::t << " s" << endl;
      if(P::t != 0.0) {
	 mpilogger << "\t (TIME) seconds per timestep " << double(after - before) / P::tsteps <<
	 ", seconds per simulated second " << double(after - before) / P::t << endl;
      }
      mpilogger << write;
   }
   
   profile::stop("Finalization");   
   profile::stop("main");
   profile::print(MPI_COMM_WORLD);
   profile::print(MPI_COMM_WORLD,0.01);

   
   

   if (myrank == MASTER_RANK) mpilogger << "(MAIN): Exiting." << endl << write;
   mpilogger.close();
   return 0;
}

