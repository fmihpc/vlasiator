#include <cstdlib>
#include <iostream>
#include <cmath>
#include <vector>
#include <sstream>

#include "definitions.h"
#include "logger.h"
#include "parameters.h"
#include "grid.h"
#include "writevars.h"

#ifndef NOCUDA
  #include "devicegrid.h"
  #include "cudafuncs.h"
#endif

using namespace std;

#ifndef NOCUDA
  extern void buildBaseGrid(Grid& grid,DeviceGrid& deviceGrid);
  extern bool acceleration(Grid& grid,DeviceGrid& deviceGrid,creal& DT);
  extern bool translation_step1(Grid& grid,DeviceGrid& deviceGrid,creal& DT);
  extern bool translation_step2(Grid& grid,DeviceGrid& deviceGrid,creal& DT);
  extern bool translation_step3(Grid& grid,DeviceGrid& deviceGrid,creal& DT);
#else
  extern void buildBaseGrid(Grid& grid);
#endif

extern bool cpu_acceleration(cuint& cellIndex,SpatialCell& cell,Grid& grid,creal& DT);
extern bool cpu_translation1(cuint& cellIndex,SpatialCell& cell,Grid& grid,creal& DT);
extern bool cpu_translation2(cuint& cellIndex,SpatialCell& cell,Grid& grid,creal& DT);
extern bool cpu_translation3(cuint& cellIndex,SpatialCell& cell,Grid& grid,creal& DT);

Logger logger("logfile.txt");

int main(int argn,char* args[]) {
   typedef Parameters P;
   
   logger << "(MAIN): Starting up." << endl;
   Parameters parameters;

   // Allocate memory for grids 
   Grid grid;
   #ifndef NOCUDA
   DeviceGrid deviceGrid;
   #endif

   // Check that requested memory was allocated:
   bool success = true;
   if (grid.isInitialized() == false) success = false;
   #ifndef NOCUDA
   if (deviceGrid.isInitialized() == false) success = false;
   #endif
   if (success == false) {
      cerr << "An error has occurred, aborting. See logfile for details." << endl;
      logger << "Aborting" << endl;
      return 1;
   }
   
   // Create initial state to CPU memory:
   #ifndef NOCUDA
   buildBaseGrid(grid,deviceGrid);
   #else
   buildBaseGrid(grid);
   #endif
   
   #ifndef NOCUDA
     // Copy volume averages, spatial cell parameters, velocity grid block parameters, 
     // and neighbour lists to device:
     for (uint i=0; i<grid.size(); ++i) {
	if (grid[i]->devSync(Cell::Blocks,Cell::CpuToDev) == false) cerr << "Error while copying blockArray" << endl;
	if (grid[i]->devSync(Cell::BlockParams,Cell::CpuToDev) == false) cerr << "Error while copying blockParams" << endl;
	if (grid[i]->devSync(Cell::CellParams,Cell::CpuToDev) == false) cerr << "Error while copying cellParams" << endl;
	if (grid[i]->devSync(Cell::NbrsVel,Cell::CpuToDev) == false) cerr << "Error while copying nbrsVel" << endl;
     }
   #endif
   
   // Write initial state:
   writeSpatialCells(grid,0);   

   vector<uint> spaIndices;
   spaIndices.push_back(0);
   //spaIndices.push_back(9);
   //spaIndices.push_back(10);   
   writeBlocks(grid,0,spaIndices);

     {
	double dens = 0.0;
	for (uint cell=0; cell<grid.size(); ++cell) {
	   for (uint b=0; b<grid[cell]->N_blocks; ++b) {
	      for (uint k=0; k<WID3; ++k) {
		 dens += grid[cell]->cpu_avgs[b*SIZE_VELBLOCK+k];
	      }
	   }
	}  
	cout << "Total density = " << dens << endl;
     }
   
   // Main simulation loop:
   logger << "(MAIN): Starting main simulation loop." << endl;
   for (uint tstep=0; tstep < P::tsteps; ++tstep) {
      #ifndef NOCUDA
        acceleration(grid,deviceGrid,0.025);
        translation_step1(grid,deviceGrid,0.025);
        translation_step2(grid,deviceGrid,0.025);
        translation_step3(grid,deviceGrid,0.025);
      #else
        for (uint i=0; i<grid.size(); ++i) cpu_acceleration(i,*(grid[i]),grid,0.025);
        for (uint i=0; i<grid.size(); ++i) cpu_translation1(i,*(grid[i]),grid,0.025);
        for (uint i=0; i<grid.size(); ++i) cpu_translation2(i,*(grid[i]),grid,0.025);
        for (uint i=0; i<grid.size(); ++i) cpu_translation3(i,*(grid[i]),grid,0.025);
      #endif
      
      // Check if the full simulation state should be written to disk
      if (tstep % P::saveInterval == 0) {
	 logger << "(MAIN): Saving full state to disk at tstep = " << tstep << endl;
      }

      // Check if variables and derived quantities should be written to disk
      if (tstep % P::diagnInterval == 0) {
	 logger << "(MAIN): Saving variables to disk at tstep = " << tstep << endl;

	 #ifndef NOCUDA
	   // Load data from device before plotting:
	   for (uint cell=0; cell<grid.size(); ++cell) {
	      grid[cell]->devSync(Cell::Blocks,Cell::DevToCpu);
	   }
	 #endif
	 writeBlocks(grid,tstep+1,spaIndices);
	 writeSpatialCells(grid,tstep+1);
      }
   }
   
   // Write final state:
   #ifndef NOCUDA
     for (uint cell=0; cell<grid.size(); ++cell) {
	grid[cell]->devSync(Cell::Blocks,Cell::DevToCpu);
     }
   #endif
     {
	double dens = 0.0;
	for (uint cell=0; cell<grid.size(); ++cell) {
	   for (uint b=0; b<grid[cell]->N_blocks; ++b) {
	      for (uint k=0; k<WID3; ++k) {
		 dens += grid[cell]->cpu_avgs[b*SIZE_VELBLOCK+k];
	      }
	   }
	}
	cout << "Total density = " << dens << endl;
     }
   writeSpatialCells(grid,P::tsteps);
   writeBlocks(grid,P::tsteps,spaIndices);

   logger << "(MAIN): Exiting." << endl;
   return 0;
}




