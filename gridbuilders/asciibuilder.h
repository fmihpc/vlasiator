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

#ifndef ASCIIREADER_H
#define ASCIIREADER_H

#include <vector>
#include "mpibuilder.h"

class AsciiBuilder: public MPIBuilder {
 public:
   AsciiBuilder();
   ~AsciiBuilder();

   bool calculatesAnalyticInitialState();
   bool finalize();
   bool getCellBlockData(const VirtualCell::ID& cellID,cuint& N_blocks,Real* blocks,Real* blockParams,uint* nbrsVel);
   bool getCellIDs(std::vector<VirtualCell::ID>& cellIDs,std::vector<uchar>& N_nbrs);
   bool getCellNbrData(const VirtualCell::ID& N_cells,VirtualCell::ID* cellIDs,Real* coords,VirtualCell::ID* spatNbrIDs,uchar* nbrTypes);
   bool getCellNumberOfBlocks(const VirtualCell::ID& N_cells,const VirtualCell::ID* const cellIDs,uint* N_blocks);
   bool getCellParams(const VirtualCell::ID& N_cells,const VirtualCell::ID* const cellIDs,Real* cellParams);
   bool getParameter(const std::string& parameterName,std::string& value);
   bool getTotalNumberOfCells(VirtualCell::ID& N_cells);
   bool initialize(MPI_Comm comm,const int& MASTER_RANK);
   
 protected:
   bool initialized;
   uint N_cells;
   uint N_dimensions;

   uint vx_blocks;    /**< Initial number of velocity blocks in vx-direction.*/
   uint vy_blocks;    /**< Initial number of velocity blocks in vy-direction.*/
   uint vz_blocks;    /**< Initial number of velocity blocks in vz-direction.*/
   uint xsize;        /**< Initial number of cells in x-direction.*/
   uint ysize;        /**< Initial number of cells in y-direction.*/
   uint zsize;        /**< Initial number of cells in z-direction.*/

   Real dx;           /**< Size of spatial cells in x-direction.*/
   Real dy;           /**< Size of spatial cells in y-direction.*/
   Real dz;           /**< Size of spatial cells in z-direction.*/
   Real xmin;         /**< Minimum value of spatial grid x-coordinate.*/
   Real xmax;         /**< Maximum value of spatial grid x-coordinate.*/
   Real ymin;         /**< Minimum value of spatial grid y-coordinate.*/
   Real ymax;         /**< Maximum value of spatial grid y-coordinate.*/
   Real zmin;         /**< Minimum value of spatial grid z-coordinate.*/
   Real zmax;         /**< Maximum value of spatial grid z-coordinate.*/
   Real vx_min;       /**< Minimum value of velocity block vx coordinate.*/
   Real vx_max;       /**< Maximum value of velocity block vx coordinate.*/
   Real vy_min;       /**< Minimum value of velocity block vy coordinate.*/
   Real vy_max;       /**< Maximum value of velocity block vy coordinate.*/
   Real vz_min;       /**< Minimum value of velocity block vz coordinate.*/
   Real vz_max;       /**< Maximum value of velocity block vz coordinate.*/
   
   uint velBlockIndex(cuint& iv,cuint& jv,cuint& kv);
};

#endif
