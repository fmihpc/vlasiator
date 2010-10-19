#include <cstdlib>
#include <iostream>
#include <cmath>
#include <vector>
#include <sstream>

#include "definitions.h"
#include "logger.h"
#include "parameters.h"
#include "grid.h"
#include "devicegrid.h"
#include "cudafuncs.h"
#include "writevars.h"

using namespace std;

extern void buildBaseGrid(Grid& grid,DeviceGrid& deviceGrid);
extern bool acceleration(Grid& grid,DeviceGrid& deviceGrid,creal& DT);
extern bool translation_step1(Grid& grid,DeviceGrid& deviceGrid,creal& DT);
extern bool translation_step2(Grid& grid,DeviceGrid& deviceGrid,creal& DT);
extern bool translation_step3(Grid& grid,DeviceGrid& deviceGrid,creal& DT);

Logger logger("logfile.txt");

int main(int argn,char* args[]) {
   typedef Parameters P;
   
   logger << "(MAIN): Starting up." << endl;
   Parameters parameters;

   // Allocate memory for grids 
   Grid grid;
   DeviceGrid deviceGrid;
   
   // Create initial state
   buildBaseGrid(grid,deviceGrid);

   // Copy volume averages, spatial cell parameters, velocity grid block parameters, 
   // and neighbour lists to device:
   for (uint i=0; i<grid.size(); ++i) {
      if (grid[i]->devSync(Cell::Blocks,Cell::CpuToDev) == false) cerr << "Error while copying blockArray" << endl;
      if (grid[i]->devSync(Cell::BlockParams,Cell::CpuToDev) == false) cerr << "Error while copying blockParams" << endl;
      if (grid[i]->devSync(Cell::CellParams,Cell::CpuToDev) == false) cerr << "Error while copying cellParams" << endl;
      if (grid[i]->devSync(Cell::NbrsVel,Cell::CpuToDev) == false) cerr << "Error while copying nbrsVel" << endl;
   }
   
   // Write initial state:
   writeSpatialCells(grid,deviceGrid,0);   
   vector<uint> spaIndices;
   spaIndices.push_back(84);
   spaIndices.push_back(118);
   spaIndices.push_back(152);   
   writeBlocks(grid,deviceGrid,0,spaIndices);
   
   // Main simulation loop:
   logger << "(MAIN): Starting main simulation loop." << endl;
   for (uint tstep=0; tstep < P::tsteps; ++tstep) {
      
      acceleration(grid,deviceGrid,0.025);

      translation_step1(grid,deviceGrid,0.025);
      translation_step2(grid,deviceGrid,0.025);
      translation_step3(grid,deviceGrid,0.025);

      // Check if the full simulation state should be written to disk
      if (tstep % P::saveInterval == 0) {
	 logger << "(MAIN): Saving full state to disk at tstep = " << tstep << endl;
      }

      // Check if variables and derived quantities should be written to disk
      if (tstep % P::diagnInterval == 0) {
	 logger << "(MAIN): Saving variables to disk at tstep = " << tstep << endl;

	 writeBlocks(grid,deviceGrid,tstep+1,spaIndices);
	 writeSpatialCells(grid,deviceGrid,tstep+1);
      }      
   }

   // Write final state:
   writeSpatialCells(grid,deviceGrid,P::tsteps);
   writeBlocks(grid,deviceGrid,P::tsteps,spaIndices);
   
   logger << "(MAIN): Exiting." << endl;
   return 0;
}




