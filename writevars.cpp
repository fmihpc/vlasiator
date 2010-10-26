#include <cstdlib>
#include <iostream>
#include <sstream>

#include "logger.h"
#include "writevars.h"
#include "silowriter.h"
#ifndef NOCUDA
  #include "cudafuncs.h"
#endif

using namespace std;

extern Logger logger;

void openOutputFile(const std::string& fileName) {
   openOutputFile(fileName.c_str(),"test");
}

void closeOutputFile(const std::string& fileName) {
   closeOutputFile();
}

void writeVelocityBlockGrid(const string& gridName,cuint& BLOCKS,real* blockParams) {
   //cout << "Writing " << BLOCKS << " blocks to silo" << endl;
   writeVelocityBlockGrid3D(gridName,BLOCKS,blockParams);
}

void writeVelocityBlockScalar(const string& varName,const string& gridName,cuint& BLOCKS,real* array) {
   //cout << "Writing data from address " << array << endl;
   writeVelocityBlockGridScalar3D(varName,gridName,BLOCKS,array);
}

void writeSpatialCells(Grid& grid,cuint& STEP) {
   stringstream fname;
   fname.width(5);
   fname.fill('0');
   fname << "cells" << STEP << ".silo";
   openOutputFile(fname.str().c_str(),"test");

   double n = 0.0;
   reserveSpatialCells(grid.size());
   for (uint cell=0; cell<grid.size(); ++cell) {
      // Integrate over all velocity blocks:
      n = 0.0;
      for (uint b=0; b<grid[cell]->N_blocks; ++b) {
	 for (uint i=0; i<WID3; ++i) n += (grid[cell]->cpu_avgs)[b*SIZE_VELBLOCK + i];
      }
      addSpatialCell(grid[cell]->cpu_cellParams,n);
   }
   writeSpatialCells("spatgrid","n");
   
   closeOutputFile();
   freeCells();
}

void writeBlocks(Grid& grid,cuint& STEP,const vector<uint>& spaIndices) {
   stringstream fname;
   fname.width(5);
   fname.fill('0');
   fname << "blocks" << STEP << ".silo";
   openOutputFile(fname.str().c_str(),"test");
   
   writeVelocityBlockGrid3D("velgrid",grid[spaIndices[0]]->N_blocks,grid[spaIndices[0]]->cpu_blockParams);
   
   for (size_t cell=0; cell<spaIndices.size(); ++cell) {
	{
	   stringstream cmp;
	   cmp.width(4);
	   cmp.fill('0');
	   cmp << "Avg" << spaIndices[cell];
	   writeVelocityBlockGridScalar3D(cmp.str().c_str(),"velgrid",grid[0]->N_blocks,grid[0]->cpu_avgs);
	}
	{
	   stringstream cmp;
	   cmp.width(4);
	   cmp.fill('0');
	   cmp << "D1x" << spaIndices[cell];
	   writeVelocityBlockGridScalar3D(cmp.str().c_str(),"velgrid",grid[0]->N_blocks,grid[0]->cpu_d1x);
	}
	{
	   stringstream cmp;
	   cmp.width(4);
	   cmp.fill('0');
	   cmp << "D1y" << spaIndices[cell];
	   writeVelocityBlockGridScalar3D(cmp.str().c_str(),"velgrid",grid[0]->N_blocks,grid[0]->cpu_d1y);
	}
	{
	   stringstream cmp;
	   cmp.width(4);
	   cmp.fill('0');
	   cmp << "D1z" << spaIndices[cell];
	   writeVelocityBlockGridScalar3D(cmp.str().c_str(),"velgrid",grid[0]->N_blocks,grid[0]->cpu_d1z);
	}
	{
	   stringstream cmp;
	   cmp.width(4);
	   cmp.fill('0');
	   cmp << "Fx" << spaIndices[cell];
	   writeVelocityBlockGridScalar3D(cmp.str().c_str(),"velgrid",grid[0]->N_blocks,grid[0]->cpu_fx);
	}
   }
   
   
   /*
   gpuCopyArray("Avg",MAX_VEL_BLOCKS*SIZE_VELBLOCK*sizeof(real),grid.getBlockArray(),deviceGrid.getBlockArray(),false);
   for (uint i=0; i<spaIndices.size(); ++i) {
      stringstream cmp; cmp.width(3); cmp.fill('0'); cmp << "Avg" << spaIndices[i];
      writeVelocityBlockGridScalar3D(cmp.str().c_str(),"velgrid",grid[spaIndices[i]]->N_blocks,grid[spaIndices[i]]->cpu_avgs);
   }
   gpuCopyArray("Fz",MAX_VEL_BLOCKS*SIZE_VELBLOCK*sizeof(real),grid.getBlockArray(),deviceGrid.getFz(),false);
   for (uint i=0; i<spaIndices.size(); ++i) {
      stringstream cmp; cmp.width(3); cmp.fill('0'); cmp << "Fz" << spaIndices[i];
      writeVelocityBlockGridScalar3D(cmp.str().c_str(),"velgrid",grid[spaIndices[i]]->N_blocks,grid[spaIndices[i]]->cpu_avgs);
   }
   gpuCopyArray("Fy",MAX_VEL_BLOCKS*SIZE_VELBLOCK*sizeof(real),grid.getBlockArray(),deviceGrid.getFy(),false);
   for (uint i=0; i<spaIndices.size(); ++i) {
      stringstream cmp; cmp.width(3); cmp.fill('0'); cmp << "Fy" << spaIndices[i];
      writeVelocityBlockGridScalar3D(cmp.str().c_str(),"velgrid",grid[spaIndices[i]]->N_blocks,grid[spaIndices[i]]->cpu_avgs);
   }
   gpuCopyArray("Fx",MAX_VEL_BLOCKS*SIZE_VELBLOCK*sizeof(real),grid.getBlockArray(),deviceGrid.getFx(),false);
   for (uint i=0; i<spaIndices.size(); ++i) {
      stringstream cmp; cmp.width(3); cmp.fill('0'); cmp << "Fx" << spaIndices[i];
      writeVelocityBlockGridScalar3D(cmp.str().c_str(),"velgrid",grid[spaIndices[i]]->N_blocks,grid[spaIndices[i]]->cpu_avgs);
   }
    */
   /*
   for (uint i=0; i<grid.size(); ++i) {
      for (uint b=0; b<grid[0]->N_blocks; ++b) {
	 for (uint k=0; k<WID3; ++k) {
	    if (grid[0]->cpu_avgs[b*SIZE_VELBLOCK+k] != grid[0]->cpu_avgs[b*SIZE_VELBLOCK+k]) {
	       cout << "NAN detected cell = " << i << " block = " << b << " cell = " << k << endl;
	    }
	 }
      }
   }
   
   for (uint i=1; i<grid.size(); ++i) {
      for (uint b=0; b<grid[0]->N_blocks; ++b) {
	 for (uint k=0; k<WID3; ++k) {
	    grid[0]->cpu_avgs[b*SIZE_VELBLOCK+k] += grid[i]->cpu_avgs[b*SIZE_VELBLOCK+k];
	 }
      }
   }
   writeVelocityBlockGrid3D("velgrid",grid[0]->N_blocks,grid[0]->cpu_blockParams);
     {
	stringstream cmp; cmp.width(3); cmp.fill('0'); cmp << "Avg" << 0;
	writeVelocityBlockGridScalar3D(cmp.str().c_str(),"velgrid",grid[0]->N_blocks,grid[0]->cpu_avgs);
     }
   
   for (uint i=1; i<grid.size(); ++i) {
      for (uint b=0; b<grid[0]->N_blocks; ++b) {
	 for (uint k=0; k<WID3; ++k) grid[0]->cpu_avgs[b*SIZE_VELBLOCK+k] += grid[i]->cpu_avgs[b*SIZE_VELBLOCK+k];
      }
   }
     {
	stringstream cmp; cmp.width(3); cmp.fill('0'); cmp << "Fx" << 0;
	writeVelocityBlockGridScalar3D(cmp.str().c_str(),"velgrid",grid[0]->N_blocks,grid[0]->cpu_avgs);
     }
   */
   closeOutputFile();
}






