#ifndef CELL_SPATIAL
#define CELL_SPATIAL

#include <utility>
#include <vector>

#include "definitions.h"

namespace Cell {
   enum CellType {INNER,BOUNDARY,GHOST,UNINITIALIZED};
   enum Array {Blocks,BlockParams,NbrsVel,CellParams,Fx,Fy,Fz,D1x,D1y,D1z};
   enum Dir {CpuToDev,DevToCpu};
};

struct SpatialCell {
   uint gpuIndex;
   uint N_blocks;

   Cell::CellType cellType; 
   uint cellIndex; 
   uint velBlockIndex; 
   uint i_ind; /**< The index of the spatial cell in x-coordinate.*/
   uint j_ind; /**< The index of the spatial cell in y-coordinate.*/
   uint k_ind; /**< The index of the spatial cell in z-coordinate.*/
   
   // Pointers to arrays containing spatial cell parameters in CPU memory
   real* cpu_cellParams;  /**< Pointer to physical cell parameters in cpu memory.*/
   uint* cpu_nbrsSpa;     /**< Pointer to spatial  neighbout list in CPU memory.*/
   
   // Pointers to arrays containing velocity block parameters in CPU memory
   uint* cpu_nbrsVel;     /**< Pointer to velocity neighbour list in CPU memory.*/
   real* cpu_avgs;        /**< Pointer to velocity block array in CPU memory.*/
   real* cpu_blockParams; /**< Pointer to velocity block parameter array in CPU memory.*/
   real* cpu_fx;
   real* cpu_fy;
   real* cpu_fz;
   real* cpu_d1x;
   real* cpu_d2x;
   real* cpu_d1y;
   real* cpu_d2y;
   real* cpu_d1z;
   real* cpu_d2z;
   
   // Pointers to arrays containing velocity block parameters in device memory
   real* gpu_cellParams;  /**< Pointer to spatial cell parameter array in GPU memory.*/
   real* gpu_avgs;        /**< Pointer to velocity block array in GPU memory.*/
   real* gpu_blockParams; /**< Pointer to velocity block parameter array in GPU memory.*/
   uint* gpu_nbrsVel;     /**< Pointer to velocity neighbout list in GPU memory.*/
   real* gpu_fx;
   real* gpu_fy;
   real* gpu_fz;
   real* gpu_d1x;
   real* gpu_d1y;
   real* gpu_d1z;
   
   SpatialCell();
   SpatialCell(const SpatialCell& s);
   ~SpatialCell();

   bool devSync(const Cell::Array& array,const Cell::Dir direction);
};

#endif

