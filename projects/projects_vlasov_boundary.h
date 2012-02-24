/*
This file is part of Vlasiator.

Copyright 2011 Finnish Meteorological Institute

Vlasiator is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3
as published by the Free Software Foundation.

Vlasiator is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Vlasiator. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef PROJECTS_VLASOVBOUNDARY_H
#define PROJECTS_VLASOVBOUNDARY_H

#include "projects_common.h"

// *********************************
// ***** TEMPLATE DECLARATIONS *****
// *********************************


template<typename CELLID,typename UINT>
void vlasovBoundaryCopyFromExistingFaceNbr(const CELLID& cellID,const UINT& boundaryFlags,const dccrg::Dccrg<SpatialCell>& mpiGrid);


// ********************************
// ***** TEMPLATE DEFINITIONS *****
// ********************************

//if the spatialcells are neighbors
inline void copyCellData(SpatialCell *from, SpatialCell *to){
   for(unsigned int block_i=0; block_i<to->number_of_blocks;block_i++){
      unsigned int block = to->velocity_block_list[block_i];         
      Velocity_Block* to_block_ptr = to->at(block);
      //if block does not exist in from cell, then the empty null_block will be returned        
      Velocity_Block* from_block_ptr = from->at(block); 
      for (uint i=0; i<SIZE_VELBLOCK; i++) to_block_ptr->data[i]=from_block_ptr->data[i];
   }
}

inline void zeroCellData( SpatialCell *to){
   for(unsigned int block_i=0; block_i<to->number_of_blocks;block_i++){
      unsigned int block = to->velocity_block_list[block_i];         
      Velocity_Block* to_block_ptr = to->at(block);
      //if block does not exist in from cell, then the empty null_block will be returned        
      for (uint i=0; i<SIZE_VELBLOCK; i++) to_block_ptr->data[i]=0.0;
   }
}


template<typename CELLID,typename UINT>
    void vlasovBoundaryCopyFromExistingFaceNbr(const CELLID& cellID,const UINT& existingCells,
                                               const UINT& nonExistingCells,const dccrg::Dccrg<SpatialCell>& mpiGrid) {
   const UINT missingPosX = (1 << 14);
   const UINT missingNegX = (1 << 12);
   const UINT missingPosY = (1 << 16);
   const UINT missingNegY = (1 << 10);
   const UINT missingPosZ = (1 << 22);
   const UINT missingNegZ = (1 << 4);
   const UINT mask = ((1 << 4) | (1 << 10) | (1 << 12) | (1 << 14) | (1 << 16) | (1 << 22));

   SpatialCell* spatialCell=mpiGrid[cellID];
   switch ((nonExistingCells & mask)) {
       case missingNegZ:
          // If +z neighbour exists, copy values from there:
          if (((existingCells >> 22) & 1) == 1){
             copyCellData(mpiGrid[getNeighbour(mpiGrid,cellID,0,0,+1)],mpiGrid[cellID]);
          }
          else{
             zeroCellData(mpiGrid[cellID]);
          }
          break;
       case missingNegY:
      // If +y neighbour exists, copy values from there:
          if (((existingCells >> 16) & 1) == 1){
             copyCellData(mpiGrid[getNeighbour(mpiGrid,cellID,0,1,0)],mpiGrid[cellID]);
          }
          else{
             zeroCellData(mpiGrid[cellID]);
          }
          break;
       case missingNegX:
          // If +x  neighbour exists, copy values from there:
          if (((existingCells >> 14) & 1) == 1) {
             copyCellData(mpiGrid[getNeighbour(mpiGrid,cellID,1,0,0)],mpiGrid[cellID]);
          }
          else{
             zeroCellData(mpiGrid[cellID]);
          }
          
          break;
       case missingPosX:
          //If  -x neighbour exists, copy values from there:
          if (((existingCells >> 12) & 1) == 1) {
             copyCellData(mpiGrid[getNeighbour(mpiGrid,cellID,-1,0,0)],mpiGrid[cellID]);
          }
          else{
             zeroCellData(mpiGrid[cellID]);
          }  
          
          break;
       case missingPosY:
          // If -y neighbour exists, copy values from there:
          if (((existingCells >> 10) & 1) == 1) {
             copyCellData(mpiGrid[getNeighbour(mpiGrid,cellID,0,-1,0)],mpiGrid[cellID]);
          }
          else{
             zeroCellData(mpiGrid[cellID]);
          }
          break;
       case missingPosZ:
          // If -z neighbour exists, copy values from there:
          if (((existingCells >>  4) & 1) == 1) {
             copyCellData(mpiGrid[getNeighbour(mpiGrid,cellID,0,0,-1)],mpiGrid[cellID]);
          }
          else{
             zeroCellData(mpiGrid[cellID]);
          }
          break;
       default:
          // Set boundary cells to zero value:
          zeroCellData(mpiGrid[cellID]);
          break;
   }
}

#endif
