#include <cstdlib>
#include <iostream>
#include <cmath>
#include <limits>

#include "mpilogger.h"
#include "common.h"
#include "cell_spatial.h"
#include "parameters.h"
#include "memalloc.h"

using namespace std;

extern Grid grid;
extern MPILogger mpilogger;

SpatialCell::SpatialCell() {

   //cout << "Spatial cell default constructor called" << endl;

   //cpu_cellParams = new Real[SIZE_CELLPARAMS];
   allocateArray(&cpu_cellParams,SIZE_CELLPARAMS);
   //cpu_nbrsSpa    = new uint[SIZE_NBRS_SPA];

   N_blocks = Parameters::vxblocks_ini * Parameters::vyblocks_ini * Parameters::vzblocks_ini;
   cpuIndex = grid.getFreeMemory(N_blocks);
   if (cpuIndex == numeric_limits<uint>::max()) {
      cerr << "Couldn't reserve memory for spatial cell" << endl;
      exit(1);
   }

   cpu_nbrsSpa     = grid.getNbrsSpa()     + cpuIndex*SIZE_NBRS_SPA;
   cpu_nbrsVel     = grid.getNbrsVel()     + cpuIndex*SIZE_NBRS_VEL;
   cpu_blockParams = grid.getBlockParams() + cpuIndex*SIZE_BLOCKPARAMS;
   cpu_avgs        = grid.getAvgs()        + cpuIndex*SIZE_VELBLOCK;
   #ifndef CUDA
   cpu_fx          = grid.getFx()          + cpuIndex*SIZE_FLUXS;
   cpu_fy          = grid.getFy()          + cpuIndex*SIZE_FLUXS;
   cpu_fz          = grid.getFz()          + cpuIndex*SIZE_FLUXS;
   cpu_d1x         = grid.getD1x()         + cpuIndex*SIZE_DERIV;
   cpu_d1y         = grid.getD1y()         + cpuIndex*SIZE_DERIV;
   cpu_d1z         = grid.getD1z()         + cpuIndex*SIZE_DERIV;
   cpu_d2x         = grid.getD2x()         + cpuIndex*SIZE_DERIV;
   cpu_d2y         = grid.getD2y()         + cpuIndex*SIZE_DERIV;
   cpu_d2z         = grid.getD2z()         + cpuIndex*SIZE_DERIV;
   #endif
}

SpatialCell::SpatialCell(const SpatialCell& s) {

   //cout << "Spatial cell copy constructor called" << endl;

   // Copy variables related to the spatial cell:
   N_blocks       = s.N_blocks;
   //cpu_cellParams = new Real[SIZE_CELLPARAMS];
   allocateArray(&cpu_cellParams,SIZE_CELLPARAMS);
   //cpu_nbrsSpa    = new uint[SIZE_NBRS_SPA];
   for (uint i=0; i<SIZE_CELLPARAMS; ++i) cpu_cellParams[i] = s.cpu_cellParams[i];
   //for (uint i=0; i<SIZE_NBRS_SPA; ++i  ) cpu_nbrsSpa[i]    = s.cpu_nbrsSpa[i];
   
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
   cpu_d1x         = s.cpu_d1x;
   cpu_d1y         = s.cpu_d1y;
   cpu_d1z         = s.cpu_d1z;
   cpu_d2x         = s.cpu_d2x;
   cpu_d2y         = s.cpu_d2y;
   cpu_d2z         = s.cpu_d2z;
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
   // Copy variables related to the spatial cell:
   N_blocks = s.N_blocks;
   for (uint i=0; i<SIZE_CELLPARAMS; ++i) cpu_cellParams[i] = s.cpu_cellParams[i];
   //for (uint i=0; i<SIZE_NBRS_SPA; ++i  ) cpu_nbrsSpa[i]    = s.cpu_nbrsSpa[i];
   
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
   cpu_d1x         = s.cpu_d1x;
   cpu_d1y         = s.cpu_d1y;
   cpu_d1z         = s.cpu_d1z;
   cpu_d2x         = s.cpu_d2x;
   cpu_d2y         = s.cpu_d2y;
   cpu_d2z         = s.cpu_d2z;
   #endif
   return *this;
}

SpatialCell::~SpatialCell() {
   //cout << "Spatial cell destructor called" << endl;
   // Free CPU memory:
   finalize();
   //freeMemory();
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
      cpu_d1x         = NULL;
      cpu_d1x         = NULL;
      cpu_d1y         = NULL;
      cpu_d2z         = NULL;
      cpu_d2y         = NULL;
      cpu_d2z         = NULL;
      #endif
   }
   
   // Attempt to get pointers to available memory from Grid:
   cpuIndex = grid.getFreeMemory(N_blocks);
   if (cpuIndex == numeric_limits<uint>::max()) return false;
   this->N_blocks = N_blocks;
   cpu_nbrsSpa     = grid.getNbrsSpa()     + cpuIndex*SIZE_NBRS_SPA;
   cpu_nbrsVel     = grid.getNbrsVel()     + cpuIndex*SIZE_NBRS_VEL;
   cpu_blockParams = grid.getBlockParams() + cpuIndex*SIZE_BLOCKPARAMS;
   cpu_avgs        = grid.getAvgs()        + cpuIndex*SIZE_VELBLOCK;
   #ifndef CUDA
   cpu_fx          = grid.getFx()          + cpuIndex*SIZE_FLUXS;
   cpu_fy          = grid.getFy()          + cpuIndex*SIZE_FLUXS;
   cpu_fz          = grid.getFz()          + cpuIndex*SIZE_FLUXS;
   cpu_d1x         = grid.getD1x()         + cpuIndex*SIZE_DERIV;
   cpu_d1y         = grid.getD1y()         + cpuIndex*SIZE_DERIV;
   cpu_d1z         = grid.getD1z()         + cpuIndex*SIZE_DERIV;
   cpu_d2x         = grid.getD2x()         + cpuIndex*SIZE_DERIV;
   cpu_d2y         = grid.getD2y()         + cpuIndex*SIZE_DERIV;
   cpu_d2z         = grid.getD2z()         + cpuIndex*SIZE_DERIV;
   #endif
   return success;
}

bool SpatialCell::finalize() {
//bool SpatialCell::freeMemory() {
   if (cpuIndex != numeric_limits<uint>::max()) {
      if (grid.removeReference(cpuIndex) == false) {
	 mpilogger << "SpatialCell ERROR: Reference removal failed" << endl << write;
      }
      cpuIndex = numeric_limits<uint>::max();
   }
   //delete cpu_cellParams;
   freeArray(cpu_cellParams);
   //delete cpu_nbrsSpa;
   return true;
}

void SpatialCell::getMemInfo() {
   mpilogger << "cpuIndex = " << cpuIndex << endl << write;
}

bool SpatialCell::clone(const SpatialCell& s) {
   //cout << "Spatial cell cloned" << endl;
   N_blocks = s.N_blocks;
   // Copy cell contents to new memory locations:
   for (uint i=0; i<SIZE_CELLPARAMS; ++i) cpu_cellParams[i] = s.cpu_cellParams[i];
   for (uint i=0; i<SIZE_NBRS_SPA; ++i  ) cpu_nbrsSpa[i]    = s.cpu_nbrsSpa[i];
   
   for (uint i=0; i<N_blocks*SIZE_NBRS_VEL; ++i)    cpu_nbrsVel[i]     = s.cpu_nbrsVel[i];
   for (uint i=0; i<N_blocks*SIZE_BLOCKPARAMS; ++i) cpu_blockParams[i] = s.cpu_blockParams[i];
   for (uint i=0; i<N_blocks*SIZE_VELBLOCK; ++i)    cpu_avgs[i]        = s.cpu_avgs[i];
   #ifndef CUDA
   for (uint i=0; i<N_blocks*SIZE_FLUXS; ++i)       cpu_fx[i]          = s.cpu_fx[i];
   for (uint i=0; i<N_blocks*SIZE_FLUXS; ++i)       cpu_fy[i]          = s.cpu_fy[i];
   for (uint i=0; i<N_blocks*SIZE_FLUXS; ++i)       cpu_fz[i]          = s.cpu_fz[i];
   for (uint i=0; i<N_blocks*SIZE_DERIV; ++i)       cpu_d1x[i]         = s.cpu_d1x[i];
   for (uint i=0; i<N_blocks*SIZE_DERIV; ++i)       cpu_d1y[i]         = s.cpu_d1y[i];
   for (uint i=0; i<N_blocks*SIZE_DERIV; ++i)       cpu_d1z[i]         = s.cpu_d1z[i];
   for (uint i=0; i<N_blocks*SIZE_DERIV; ++i)       cpu_d2x[i]         = s.cpu_d2x[i];
   for (uint i=0; i<N_blocks*SIZE_DERIV; ++i)       cpu_d2y[i]         = s.cpu_d2y[i];
   for (uint i=0; i<N_blocks*SIZE_DERIV; ++i)       cpu_d2z[i]         = s.cpu_d2z[i];
   #endif
   return  true;
}
#ifndef PARGRID
uint SpatialCell::base_address_identifier = 0;

size_t SpatialCell::size(void) {
   cuint N_BLOCKS = Parameters::vxblocks_ini * Parameters::vyblocks_ini * Parameters::vzblocks_ini;
   switch (SpatialCell::base_address_identifier) {
   case 0:
      return sizeof(Real) * N_BLOCKS * SIZE_VELBLOCK;
   case 1:
      return sizeof(Real) * 7 * MAX_VEL_BLOCKS * SIZE_VELBLOCK;
   case 2:
      return sizeof(Real) * 3 * MAX_VEL_BLOCKS * SIZE_VELBLOCK;
      break;
   default:
      return sizeof(Real) * N_BLOCKS * SIZE_VELBLOCK;
      break;
   }
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

/*
#ifdef PARGRID
#include <mpi.h>
void SpatialCell::allocate() {
   
}

#endif
*/

void* SpatialCell::getBaseAddress(cuint identifier) {
   // Hack: make sure that the pointers in SpatialCell are correct:
   switch (identifier) {
    case 0:
      return cpu_avgs;
      break;
    #ifndef CUDA
    case 1:
      return cpu_d1x;
      break;
    case 2:
      return cpu_fx;
      break;
    #endif
   }
   return cpu_avgs;
}

void SpatialCell::getMPIdatatype(cuint identifier,MPI_Datatype& dataType) {
   // Another hack: this is a static member function, so SpatialCell::N_blocks 
   // cannot be used below:
   typedef Parameters P;
   cuint N_BLOCKS = P::vxblocks_ini*P::vyblocks_ini*P::vzblocks_ini;
   
   MPI_Datatype dataTypes[7];
   int blockLengths[7];
   MPI_Aint displacements[7];
   
   switch (identifier) {
    case 0: // Transfer averages:
      dataTypes[0] = MPI_FLOAT;
      blockLengths[0] = N_BLOCKS*SIZE_VELBLOCK;
      displacements[0] = 0;                   // Base address is cpu_avgs
      if (MPI_Type_create_struct(1,blockLengths,displacements,dataTypes,&dataType) != MPI_SUCCESS) {
	 #ifndef NDEBUG
	    mpilogger << "SpatialCell::getMPIdatatype ERROR failed to create MPI_Datatype!" << endl << write;
         #endif
      }
      break;
    case 1: // Transfer 1st derivatives:
      for (int i=0; i<3; ++i) dataTypes[i] = MPI_FLOAT;
      for (int i=0; i<3; ++i) blockLengths[i] = N_BLOCKS*SIZE_VELBLOCK;
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
      for (int i=0; i<3; ++i) dataTypes[i] = MPI_FLOAT;
      for (int i=0; i<3; ++i) blockLengths[i] = N_BLOCKS*SIZE_VELBLOCK;
      displacements[0] = 0;                   // Base address is cpu_fx
      displacements[1] = 1*MAX_VEL_BLOCKS*SIZE_VELBLOCK*sizeof(Real);
      displacements[2] = 2*MAX_VEL_BLOCKS*SIZE_VELBLOCK*sizeof(Real);
      if (MPI_Type_create_struct(3,blockLengths,displacements,dataTypes,&dataType) != MPI_SUCCESS) {
	 #ifndef NDEBUG
	    mpilogger << "SpatialCell::getMPIdatatype ERROR failed to create MPI_Datatype!" << endl << write;
	 #endif
      }
      break;
    case 3: // transfer everything
      cerr << "getMPIdatatype for identifier 3 not implemented yet" << endl;
      abort();
      break;
   }
}

