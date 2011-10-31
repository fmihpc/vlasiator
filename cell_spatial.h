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

#ifndef CELL_SPATIAL
#define CELL_SPATIAL

#include <cstring>
#include <utility>
#include <vector>
#include "definitions.h"
#include "common.h"
#include "parameters.h"
#include "grid.h"

extern Grid grid;

namespace Cell {
   enum CellType {INNER,BOUNDARY,GHOST,UNINITIALIZED};
   enum Array {Blocks,BlockParams,NbrsVel,CellParams,Fx,Fy,Fz,D1x,D1y,D1z};
   enum Dir {CpuToDev,DevToCpu};
}

struct SpatialCell {
   size_t size(void);
   static uint base_address_identifier;
   void* at(void);
   MPI_Datatype mpi_datatype(void);
   void* getBaseAddress(cuint identifier);
   void getMPIdatatype(cuint identifier,MPI_Datatype& dataType);
   
   uint cpuIndex;         /**< An index to Grid which is used to calculate the data array pointers (cpu_avgs etc.).*/
   bool isGhostCell;      /**< If true, this cell is a ghost cell, i.e. Vlasov boundary condition functions 
			   * need to be called for this cell.*/
   uint N_blocks;         /**< Number of velocity blocks in this cell.*/

   // Pointers to arrays containing spatial cell parameters in CPU memory
   Real* cpu_cellParams;  /**< Pointer to physical cell parameters in cpu memory.*/
   uint* cpu_nbrsSpa;     /**< Pointer to spatial  neighbout list in CPU memory.*/
   Real* cpu_derivatives;
   
   // Pointers to arrays containing velocity block parameters in CPU memory
   uint* cpu_nbrsVel;     /**< Pointer to velocity neighbour list in CPU memory.*/
   Real* cpu_avgs;        /**< Pointer to velocity block array in CPU memory.*/
   Real* cpu_blockParams; /**< Pointer to velocity block parameter array in CPU memory.*/
   #ifndef CUDA
   Real* cpu_fx;          /**< Pointer to x-flux array in CPU memory.*/
   Real* cpu_fy;          /**< Pointer to y-flux array in CPU memory.*/
   Real* cpu_fz;          /**< Pointer to z-flux array in CPU memory.*/
   #ifdef SOLVER_KT
   Real* cpu_d1x;         /**< Pointer to array in CPU memory that contains 1st derivatives to x-direction.*/
   Real* cpu_d2x;         /**< Pointer to array in CPU memory that contains 2nd derivatives to x-direction.*/
   Real* cpu_d1y;         /**< Pointer to array in CPU memory that contains 1st derivatives to y-direction.*/
   Real* cpu_d2y;         /**< Pointer to array in CPU memory that contains 2nd derivatives to y-direction.*/
   Real* cpu_d1z;         /**< Pointer to array in CPU memory that contains 1st derivatives to z-direction.*/
   Real* cpu_d2z;         /**< Pointer to array in CPU memory that contains 2nd derivatives to z-direction.*/
   #endif
   #endif 
   SpatialCell();
   SpatialCell(const SpatialCell& s);
   ~SpatialCell();
   SpatialCell& operator=(const SpatialCell& s);

   #ifndef NOCUDA
     bool devSync(const Cell::Array& array,const Cell::Dir direction);
   #endif

   bool clone(const SpatialCell& s);
   void getMemInfo();
   
   bool initialize(cuint& N_blocks);
   bool finalize();
};

#endif

