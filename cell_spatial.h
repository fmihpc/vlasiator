#ifndef CELL_SPATIAL
#define CELL_SPATIAL

#include <utility>
#include <vector>

#include "definitions.h"
#include "common.h"
#include "parameters.h"

namespace Cell {
   enum CellType {INNER,BOUNDARY,GHOST,UNINITIALIZED};
   enum Array {Blocks,BlockParams,NbrsVel,CellParams,Fx,Fy,Fz,D1x,D1y,D1z};
   enum Dir {CpuToDev,DevToCpu};
};

struct SpatialCell {
   /** A function for loading (and saving) spatial cell data from archive.
    * Boost::mpi uses archives to transmit data between processes, and 
    * this function really just determines the data that is exchanged between
    * processes.
    * The transmitted data is determined from the value of Parameters::transmit.
    * Set its value with bitwise or operations of the constants stored in 
    * namespace Transmit. For example, if Parameters::transmit == 
    * Transmit::AVGS | Transmit::FLUXES, then volume averages of distrib. function 
    * and the fluxes are transmitted.
    * @param ar Archive where the cell data is loaded from or saved to.
    * @param version The version number of the Archive, may be useful for
    * backward compability.
    */
   template<typename Archiver> void serialize(Archiver& ar,cuint& version) {
      typedef Parameters P;
      
      if ((P::transmit & Transmit::CELL_PARAMS)  > 0) for (uint i=0; i<SIZE_CELLPARAMS; ++i) ar & cpu_cellParams[i];
      if ((P::transmit & Transmit::BLOCK_PARAMS) > 0) for (uint i=0; i<SIZE_BLOCKPARAMS; ++i) ar & cpu_blockParams[i];
      if ((P::transmit & Transmit::AVGS)         > 0) for (uint i=0; i<N_blocks*SIZE_VELBLOCK; ++i) ar & cpu_avgs[i];
      if ((P::transmit & Transmit::FLUXES)       > 0) {
	 for (uint i=0; i<N_blocks*SIZE_FLUXS; ++i) ar & cpu_fx[i];
	 for (uint i=0; i<N_blocks*SIZE_FLUXS; ++i) ar & cpu_fy[i];
	 for (uint i=0; i<N_blocks*SIZE_FLUXS; ++i) ar & cpu_fz[i];
      }
      if ((P::transmit & Transmit::DERIV1)       > 0) {
	 for (uint i=0; i<N_blocks*SIZE_DERIV; ++i) ar & cpu_d1x[i];
	 for (uint i=0; i<N_blocks*SIZE_DERIV; ++i) ar & cpu_d1y[i];
	 for (uint i=0; i<N_blocks*SIZE_DERIV; ++i) ar & cpu_d1z[i];
      }
      if ((P::transmit & Transmit::DERIV2)       > 0) {
	 for (uint i=0; i<N_blocks*SIZE_DERIV; ++i) ar & cpu_d2x[i];
	 for (uint i=0; i<N_blocks*SIZE_DERIV; ++i) ar & cpu_d2y[i];
	 for (uint i=0; i<N_blocks*SIZE_DERIV; ++i) ar & cpu_d2z[i];
      }
      if ((P::transmit & Transmit::NBRSVEL)      > 0) for (uint i=0; i<N_blocks*SIZE_NBRS_VEL; ++i) ar & cpu_nbrsVel[i];
      /*
      for (uint i=0; i<SIZE_CELLPARAMS; ++i) ar & cpu_cellParams[i];
      for (uint i=0; i<SIZE_BLOCKPARAMS; ++i) ar & cpu_blockParams[i];
      for (uint i=0; i<N_blocks*SIZE_VELBLOCK; ++i) ar & cpu_avgs[i];
      for (uint i=0; i<N_blocks*SIZE_FLUXS; ++i) ar & cpu_fx[i];
      for (uint i=0; i<N_blocks*SIZE_FLUXS; ++i) ar & cpu_fy[i];
      for (uint i=0; i<N_blocks*SIZE_FLUXS; ++i) ar & cpu_fz[i];
      for (uint i=0; i<N_blocks*SIZE_DERIV; ++i) ar & cpu_d1x[i];
      for (uint i=0; i<N_blocks*SIZE_DERIV; ++i) ar & cpu_d1y[i];
      for (uint i=0; i<N_blocks*SIZE_DERIV; ++i) ar & cpu_d1z[i];
      for (uint i=0; i<N_blocks*SIZE_DERIV; ++i) ar & cpu_d2x[i];
      for (uint i=0; i<N_blocks*SIZE_DERIV; ++i) ar & cpu_d2y[i];
      for (uint i=0; i<N_blocks*SIZE_DERIV; ++i) ar & cpu_d2z[i];
      for (uint i=0; i<N_blocks*SIZE_NBRS_VEL; ++i) ar & cpu_nbrsVel[i];
      */
   }
   
   bool pageLocked; /**< If true, the cell's memory is page-locked.*/
   uint cpuIndex;
   uint gpuIndex;
   uint N_blocks;

   Cell::CellType cellType; 
   uint cellIndex; 
   uint velBlockIndex; 
   uint i_ind; /**< The index of the spatial cell in x-coordinate.*/
   uint j_ind; /**< The index of the spatial cell in y-coordinate.*/
   uint k_ind; /**< The index of the spatial cell in z-coordinate.*/
   
   // Pointers to arrays containing spatial cell parameters in CPU memory
   Real* cpu_cellParams;  /**< Pointer to physical cell parameters in cpu memory.*/
   uint* cpu_nbrsSpa;     /**< Pointer to spatial  neighbout list in CPU memory.*/
   
   // Pointers to arrays containing velocity block parameters in CPU memory
   uint* cpu_nbrsVel;     /**< Pointer to velocity neighbour list in CPU memory.*/
   Real* cpu_avgs;        /**< Pointer to velocity block array in CPU memory.*/
   Real* cpu_blockParams; /**< Pointer to velocity block parameter array in CPU memory.*/
   Real* cpu_fx;
   Real* cpu_fy;
   Real* cpu_fz;
   Real* cpu_d1x;
   Real* cpu_d2x;
   Real* cpu_d1y;
   Real* cpu_d2y;
   Real* cpu_d1z;
   Real* cpu_d2z;
   
   // Pointers to arrays containing velocity block parameters in device memory
   Real* gpu_cellParams;  /**< Pointer to spatial cell parameter array in GPU memory.*/
   Real* gpu_avgs;        /**< Pointer to velocity block array in GPU memory.*/
   Real* gpu_blockParams; /**< Pointer to velocity block parameter array in GPU memory.*/
   uint* gpu_nbrsVel;     /**< Pointer to velocity neighbout list in GPU memory.*/
   Real* gpu_fx;
   Real* gpu_fy;
   Real* gpu_fz;
   Real* gpu_d1x;
   Real* gpu_d1y;
   Real* gpu_d1z;
   
   SpatialCell();
   SpatialCell(const SpatialCell& s);
   ~SpatialCell();
   SpatialCell& operator=(const SpatialCell& s);

   #ifndef NOCUDA
     bool devSync(const Cell::Array& array,const Cell::Dir direction);
   #endif
   
 private:
   void allocateMemory();
   void freeMemory();
   void swapMemoryType();
   void copyMemory(const SpatialCell& s);
   void copyVariables(const SpatialCell& s);
};

#endif

