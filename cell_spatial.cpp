/*
This file is part of Vlasiator.

Copyright 2010, 2011 Finnish Meteorological Institute

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
#include <limits>

#include "mpilogger.h"
#include "mpiconversion.h"
#include "common.h"
#include "cell_spatial.h"
#include "parameters.h"
#include "memalloc.h"

using namespace std;
using namespace fieldsolver;

extern Grid grid;
extern MPILogger mpilogger;

SpatialCell::SpatialCell() {
   //cout << "Spatial cell default constructor called" << endl;
   N_blocks = Parameters::vzblocks_ini * Parameters::vyblocks_ini * Parameters::vxblocks_ini;
   cpuIndex = numeric_limits<uint>::max();
   isGhostCell = false;
   
   allocateArray(&cpu_cellParams,CellParams::SIZE_CELLPARAMS);
   try {
      cpu_derivatives = new Real[SIZE_DERIVATIVES];
   }
   catch (exception& e) {
      cerr << __FILE__ << ":" << __LINE__
         << "Couldn't allocate memory for cpu_derivatives: " << e.what()
         << endl;
      abort();
   }
   cpu_nbrsSpa     = NULL;
   cpu_nbrsVel     = NULL;
   cpu_blockParams = NULL;
   cpu_avgs        = NULL;
   #ifndef CUDA
   cpu_fx          = NULL;
   cpu_fy          = NULL;
   cpu_fz          = NULL;
   #ifdef SOLVER_KT
   cpu_d1x         = NULL;
   cpu_d1y         = NULL;
   cpu_d1z         = NULL;
   cpu_d2x         = NULL;
   cpu_d2y         = NULL;
   cpu_d2z         = NULL;
   #endif
   #endif
}

SpatialCell::SpatialCell(const SpatialCell& s) {
   //cout << "Spatial cell copy constructor called" << endl;

   try {
      cpu_derivatives = new Real[SIZE_DERIVATIVES];
   }
   catch (exception& e) {
      cerr << __FILE__ << ":" << __LINE__
         << "Couldn't allocate memory for cpu_derivatives: " << e.what()
         << endl;
      abort();
   }
   
   // Copy variables related to the spatial cell:
   N_blocks       = s.N_blocks;
   isGhostCell    = s.isGhostCell;
   allocateArray(&cpu_cellParams,CellParams::SIZE_CELLPARAMS);
   for (uint i=0; i<CellParams::SIZE_CELLPARAMS; ++i) cpu_cellParams[i] = s.cpu_cellParams[i];
   for (uint i=0; i<SIZE_DERIVATIVES; ++i) cpu_derivatives[i] = s.cpu_derivatives[i];
   
   cpuIndex = s.cpuIndex;
   // If the SpatialCell to copy has allocated memory, increase the reference count 
   // and set pointers:
   if (cpuIndex == numeric_limits<uint>::max()) return;
   if (grid.addReference(cpuIndex) == false) {
      mpilogger << "SpatialCell: reference increase failed, aborting." << endl << write;
      exit(1);
   }
   cpu_nbrsSpa     = s.cpu_nbrsSpa;
   cpu_nbrsVel     = s.cpu_nbrsVel;
   cpu_blockParams = s.cpu_blockParams;
   cpu_avgs        = s.cpu_avgs;
   #ifndef CUDA
   cpu_fx          = s.cpu_fx;
   cpu_fy          = s.cpu_fy;
   cpu_fz          = s.cpu_fz;
   #ifdef SOLVER_KT
   cpu_d1x         = s.cpu_d1x;
   cpu_d1y         = s.cpu_d1y;
   cpu_d1z         = s.cpu_d1z;
   cpu_d2x         = s.cpu_d2x;
   cpu_d2y         = s.cpu_d2y;
   cpu_d2z         = s.cpu_d2z;
   #endif
   #endif
}

SpatialCell& SpatialCell::operator=(const SpatialCell& s) {
   //cout << "Spatial cell assignment operator called" << endl;

   // Clear previous memory:
   if (cpuIndex != numeric_limits<uint>::max()) {
      if (grid.removeReference(cpuIndex) == false) {
	 mpilogger << "SpatialCell operator=: Failed to remove reference." << endl << write;
      }
   }
   if (cpu_cellParams != NULL) freeArray(cpu_cellParams);
   allocateArray(&cpu_cellParams,CellParams::SIZE_CELLPARAMS);
   
   // Copy variables related to the spatial cell:
   N_blocks = s.N_blocks;
   isGhostCell = s.isGhostCell;
   for (uint i=0; i<CellParams::SIZE_CELLPARAMS; ++i) cpu_cellParams[i] = s.cpu_cellParams[i];
   for (uint i=0; i<SIZE_DERIVATIVES; ++i) cpu_derivatives[i] = s.cpu_derivatives[i];
   
   // Copy variables related to the velocity grid:
   cpuIndex = s.cpuIndex;
   if (cpuIndex == numeric_limits<uint>::max()) return *this;
   
   if (grid.addReference(cpuIndex) == false) {
      mpilogger << "SpatialCell: reference increase failed, aborting." << endl;
      mpilogger << "\t operator= got cpuIndex " << s.cpuIndex << " from copied SpatialCell." << endl << write;
      exit(1);
   }
   cpu_nbrsSpa     = s.cpu_nbrsSpa;
   cpu_nbrsVel     = s.cpu_nbrsVel;
   cpu_blockParams = s.cpu_blockParams;
   cpu_avgs        = s.cpu_avgs;
   #ifndef CUDA
   cpu_fx          = s.cpu_fx;
   cpu_fy          = s.cpu_fy;
   cpu_fz          = s.cpu_fz;
   #ifdef SOLVER_KT
   cpu_d1x         = s.cpu_d1x;
   cpu_d1y         = s.cpu_d1y;
   cpu_d1z         = s.cpu_d1z;
   cpu_d2x         = s.cpu_d2x;
   cpu_d2y         = s.cpu_d2y;
   cpu_d2z         = s.cpu_d2z;
   #endif
   #endif
   return *this;
}

SpatialCell::~SpatialCell() {
   //cout << "Spatial cell destructor called" << endl;
   // Free CPU memory:
   finalize();
}

bool SpatialCell::initialize(cuint& N_blocks) {
   bool success = true;
   // Deallocate previous memory (if required):
   if (cpuIndex != numeric_limits<uint>::max()) {
      if (grid.removeReference(cpuIndex) == false) success = false;
      cpu_nbrsSpa     = NULL;
      cpu_nbrsVel     = NULL;
      cpu_blockParams = NULL;
      cpu_avgs        = NULL;
      #ifndef CUDA
      cpu_fx          = NULL;
      cpu_fy          = NULL;
      cpu_fz          = NULL;
      #ifdef SOLVER_KT
      cpu_d1x         = NULL;
      cpu_d1x         = NULL;
      cpu_d1y         = NULL;
      cpu_d2z         = NULL;
      cpu_d2y         = NULL;
      cpu_d2z         = NULL;
      #endif
      #endif
   }
   
   // Attempt to get pointers to available memory from Grid:
   cpuIndex = grid.getFreeMemory(N_blocks);
   if (cpuIndex == numeric_limits<uint>::max()) {
      // FIXME: following doesn't produce output if vlasiator segfaults right afterwards
      mpilogger << "SpatialCell " << cpuIndex << " : initialization failed." << endl;
      return false;
   }
   this->N_blocks  = N_blocks;
   isGhostCell     = false;
   cpu_nbrsSpa     = grid.getNbrsSpa()     + cpuIndex*SIZE_NBRS_SPA;
   cpu_nbrsVel     = grid.getNbrsVel()     + cpuIndex*SIZE_NBRS_VEL;
   cpu_blockParams = grid.getBlockParams() + cpuIndex*SIZE_BLOCKPARAMS;
   cpu_avgs        = grid.getAvgs()        + cpuIndex*SIZE_VELBLOCK;
   #ifndef CUDA
   cpu_fx          = grid.getFx()          + cpuIndex*SIZE_FLUXS;
   cpu_fy          = grid.getFy()          + cpuIndex*SIZE_FLUXS;
   cpu_fz          = grid.getFz()          + cpuIndex*SIZE_FLUXS;
   #ifdef SOLVER_KT
   cpu_d1x         = grid.getD1x()         + cpuIndex*SIZE_DERIV;
   cpu_d1y         = grid.getD1y()         + cpuIndex*SIZE_DERIV;
   cpu_d1z         = grid.getD1z()         + cpuIndex*SIZE_DERIV;
   cpu_d2x         = grid.getD2x()         + cpuIndex*SIZE_DERIV;
   cpu_d2y         = grid.getD2y()         + cpuIndex*SIZE_DERIV;
   cpu_d2z         = grid.getD2z()         + cpuIndex*SIZE_DERIV;
   #endif
   #endif
   return success;
}

bool SpatialCell::finalize() {
   if (cpuIndex != numeric_limits<uint>::max()) {
      if (grid.removeReference(cpuIndex) == false) {
	 mpilogger << "SpatialCell ERROR: Reference removal failed" << endl << write;
      }
      cpuIndex = numeric_limits<uint>::max();
   }
   freeArray(cpu_cellParams);
   delete cpu_derivatives;
   cpu_derivatives = NULL;
   return true;
}

void SpatialCell::getMemInfo() {
   mpilogger << "cpuIndex = " << cpuIndex << endl << write;
}

bool SpatialCell::clone(const SpatialCell& s) {
   //cout << "Spatial cell cloned" << endl;
   N_blocks = s.N_blocks;
   isGhostCell = s.isGhostCell;
   // Copy cell contents to new memory locations:
   for (uint i=0; i<CellParams::SIZE_CELLPARAMS; ++i) cpu_cellParams[i] = s.cpu_cellParams[i];
   for (uint i=0; i<SIZE_DERIVATIVES; ++i) cpu_derivatives[i] = s.cpu_derivatives[i];
   for (uint i=0; i<SIZE_NBRS_SPA; ++i  ) cpu_nbrsSpa[i]    = s.cpu_nbrsSpa[i];
   
   for (uint i=0; i<N_blocks*SIZE_NBRS_VEL; ++i)    cpu_nbrsVel[i]     = s.cpu_nbrsVel[i];
   for (uint i=0; i<N_blocks*SIZE_BLOCKPARAMS; ++i) cpu_blockParams[i] = s.cpu_blockParams[i];
   for (uint i=0; i<N_blocks*SIZE_VELBLOCK; ++i)    cpu_avgs[i]        = s.cpu_avgs[i];
   #ifndef CUDA
   for (uint i=0; i<N_blocks*SIZE_FLUXS; ++i)       cpu_fx[i]          = s.cpu_fx[i];
   for (uint i=0; i<N_blocks*SIZE_FLUXS; ++i)       cpu_fy[i]          = s.cpu_fy[i];
   for (uint i=0; i<N_blocks*SIZE_FLUXS; ++i)       cpu_fz[i]          = s.cpu_fz[i];
   #ifdef SOLVER_KT
   for (uint i=0; i<N_blocks*SIZE_DERIV; ++i)       cpu_d1x[i]         = s.cpu_d1x[i];
   for (uint i=0; i<N_blocks*SIZE_DERIV; ++i)       cpu_d1y[i]         = s.cpu_d1y[i];
   for (uint i=0; i<N_blocks*SIZE_DERIV; ++i)       cpu_d1z[i]         = s.cpu_d1z[i];
   for (uint i=0; i<N_blocks*SIZE_DERIV; ++i)       cpu_d2x[i]         = s.cpu_d2x[i];
   for (uint i=0; i<N_blocks*SIZE_DERIV; ++i)       cpu_d2y[i]         = s.cpu_d2y[i];
   for (uint i=0; i<N_blocks*SIZE_DERIV; ++i)       cpu_d2z[i]         = s.cpu_d2z[i];
   #endif
   #endif
   return  true;
}
#ifndef PARGRID
uint SpatialCell::base_address_identifier = 0;

size_t SpatialCell::size(void) {
	return 1;
}

MPI_Datatype SpatialCell::mpi_datatype(void)
{
   MPI_Datatype data_type;

   SpatialCell::getMPIdatatype(SpatialCell::base_address_identifier, data_type);
   return data_type;
}

void* SpatialCell::at(void) {
   return this->getBaseAddress(SpatialCell::base_address_identifier);
}
#endif

void* SpatialCell::getBaseAddress(cuint identifier) {
   switch (identifier) {
    case 0:
      return this->cpu_avgs;
      break;
    #ifndef CUDA
    #ifdef SOLVER_KT    
    case 1:
      return this->cpu_d1x;
      break;
    #endif  
    case 2:
      return this->cpu_fx;
      break;
    case 3:
      return this->cpu_cellParams + CellParams::BX;
      break;
    case 4:
      return this->cpu_derivatives;
      break;
    case 5:
      return &(this->N_blocks);
    case 6:
      return this->cpu_cellParams + CellParams::EX;
      break;
    #endif
    case 7:
      return &(isGhostCell);
      break;
    default:
      cerr << "Unsupported base address" << endl;
      abort();
      break;
   }
   return cpu_avgs;
}

void SpatialCell::getMPIdatatype(cuint identifier,MPI_Datatype& dataType) {
   
   MPI_Datatype dataTypes[7];
   int blockLengths[7];
   MPI_Aint displacements[7];
   
   switch (identifier) {
    case 0: // Transfer averages:
      dataTypes[0] = MPI_Type<Real>();
      blockLengths[0] = this->N_blocks*SIZE_VELBLOCK;
      displacements[0] = 0;                   // Base address is cpu_avgs
      if (MPI_Type_create_struct(1,blockLengths,displacements,dataTypes,&dataType) != MPI_SUCCESS) {
	 #ifndef NDEBUG
	    mpilogger << "SpatialCell::getMPIdatatype ERROR failed to create MPI_Datatype!" << endl << write;
         #endif
      }
      break;
    case 1: // Transfer 1st derivatives:
      for (int i=0; i<3; ++i) dataTypes[i] =  MPI_Type<Real>();;
      for (int i=0; i<3; ++i) blockLengths[i] = this->N_blocks*SIZE_VELBLOCK;
      displacements[0] = 0;                       // Base address is cpu_d1x
      displacements[1] = 1*MAX_VEL_BLOCKS*SIZE_VELBLOCK*sizeof(Real); // d1y
      displacements[2] = 2*MAX_VEL_BLOCKS*SIZE_VELBLOCK*sizeof(Real); // d1z
      if (MPI_Type_create_struct(3,blockLengths,displacements,dataTypes,&dataType) != MPI_SUCCESS) {
         #ifndef NDEBUG
	    mpilogger << "SpatialCell::getMPIdatatype ERROR failed to create MPI_Datatype!" << endl << write;
         #endif
      }      
      break;
    case 2: // Transfer fluxes:
      for (int i=0; i<3; ++i) dataTypes[i] =  MPI_Type<Real>();
      for (int i=0; i<3; ++i) blockLengths[i] = this->N_blocks*SIZE_VELBLOCK;
      displacements[0] = 0;                   // Base address is cpu_fx
      displacements[1] = 1*MAX_VEL_BLOCKS*SIZE_VELBLOCK*sizeof(Real);
      displacements[2] = 2*MAX_VEL_BLOCKS*SIZE_VELBLOCK*sizeof(Real);
      if (MPI_Type_create_struct(3,blockLengths,displacements,dataTypes,&dataType) != MPI_SUCCESS) {
	 #ifndef NDEBUG
	    mpilogger << "SpatialCell::getMPIdatatype ERROR failed to create MPI_Datatype!" << endl << write;
	 #endif
      }
      break;
    case 3: // transfer cell BX,BY,BZ,RHO,RHOVX,RHOVY,RHOVZ
      dataTypes[0] = MPI_BYTE;
      blockLengths[0] = sizeof(Real) * 7;
      displacements[0] = 0;	// Base address starts at CellParams::BX
      if (MPI_Type_create_struct(1,blockLengths,displacements,dataTypes,&dataType) != MPI_SUCCESS) {
	 #ifndef NDEBUG
	    mpilogger << "SpatialCell::getMPIdatatype ERROR failed to create MPI_Datatype!" << endl << write;
	 #endif
      }
      break;
   case 4: // transfer cell derivatives
      dataTypes[0] = MPI_BYTE;
      blockLengths[0] = sizeof(Real) * SIZE_DERIVATIVES;
      displacements[0] = 0;	// Base address is cpu_derivatives
      if (MPI_Type_create_struct(1,blockLengths,displacements,dataTypes,&dataType) != MPI_SUCCESS) {
	 #ifndef NDEBUG
	    mpilogger << "SpatialCell::getMPIdatatype ERROR failed to create MPI_Datatype!" << endl << write;
	 #endif
      }
      break;
   case 5: // transfer N_blocks
      dataTypes[0] = MPI_BYTE;
      blockLengths[0] = sizeof(uint);
      displacements[0] = 0;	// Base address is N_blocks
      if (MPI_Type_create_struct(1,blockLengths,displacements,dataTypes,&dataType) != MPI_SUCCESS) {
	 #ifndef NDEBUG
	    mpilogger << "SpatialCell::getMPIdatatype ERROR failed to create MPI_Datatype!" << endl << write;
	 #endif
      }
      break;
    case 6: // transfer cell EX,EY,EZ
      dataTypes[0] = MPI_BYTE;
      blockLengths[0] = sizeof(Real) * 3;
      displacements[0] = 0;	// Base address starts at CellParams::EX
      if (MPI_Type_create_struct(1,blockLengths,displacements,dataTypes,&dataType) != MPI_SUCCESS) {
	 #ifndef NDEBUG
	    mpilogger << "SpatialCell::getMPIdatatype ERROR failed to create MPI_Datatype!" << endl << write;
	 #endif
      }
      break;
    case 7: // Transfer isGhostCell
      dataTypes[0]     = MPI_BYTE;
      blockLengths[0]  = sizeof(bool);
      displacements[0] = 0;
      if (MPI_Type_create_struct(1,blockLengths,displacements,dataTypes,&dataType) != MPI_SUCCESS) {
	 #ifndef NDEBUG
	    mpilogger << "SpatialCell::getMPIdatatype ERROR failed to create MPI_Datatype!" << endl << write;
         #endif
      }
      break;
   }
}

