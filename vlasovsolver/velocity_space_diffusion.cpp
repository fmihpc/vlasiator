/*
 * This file is part of Vlasiator.
 * Copyright 2010-2020 University of Helsinki
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://www.physics.helsinki.fi/vlasiator/
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#include <zoltan.h>
#include <dccrg.hpp>
#include "../common.h"
#include "../spatial_cell.hpp"
#include <dccrg_cartesian_geometry.hpp>
#include <vector3d.h>
#include "../parameters.h"
#include <math.h>


using namespace spatial_cell;

static bool checkExistingNeighbour(SpatialCell* cell, Realf VX, Realf VY, Realf VZ, const uint popID) {

      const vmesh::GlobalID blockGID = cell->get_velocity_block(popID,VX, VY, VZ);
      vmesh::LocalID blockLID = cell->get_population(popID).vmesh.getLocalID(blockGID);
      return blockLID != vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>::invalidLocalID();
}

void velocitySpaceDiffusion(
        dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,const uint popID){

    const auto LocalCells=getLocalCells(); 
    for (auto & CellID: LocalCells) { //Iterate through spatial cell

        SpatialCell& cell = *mpiGrid[CellID];

        std::vector<std::array<Realf,3>> arraydf(cell.get_number_of_velocity_blocks(popID)*WID3); // Array of vspace size

        for (int coord = 0; coord < 3; coord++) { // First derivative loop

	   Vec3d B(cell.parameters[CellParams::PERBXVOL] +  cell.parameters[CellParams::BGBXVOL],
                   cell.parameters[CellParams::PERBYVOL] +  cell.parameters[CellParams::BGBYVOL],
	           cell.parameters[CellParams::PERBZVOL] +  cell.parameters[CellParams::BGBZVOL]);
           Vec3d b = normalize_vector(B);


	   const Real* parameters  = cell.get_block_parameters(popID);

           const vmesh::LocalID* nBlocks = cell.get_velocity_grid_length(popID);

            for (vmesh::LocalID n=0; n<cell.get_number_of_velocity_blocks(popID); n++) { //Iterate through velocity blocks
               for (uint k = 0; k < WID; ++k) for (uint j = 0; j < WID; ++j) for (uint i = 0; i < WID; ++i) {

                  //Get velocity space coordinates                    
	          const Real VX 
                     =          parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VXCRD] 
                     + (i + 0.5)*parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVX];
                  const Real VY 
                     =          parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VYCRD] 
                     + (j + 0.5)*parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVY];
                  const Real VZ 
                     =          parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VZCRD]
                     + (k + 0.5)*parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVZ];
                  
                  Vec3d bulkV(cell.parameters[CellParams::VX], cell.parameters[CellParams::VY], cell.parameters[CellParams::VZ]);

                  Vec3d V(VX - bulkV[0] , VY - bulkV[1] , VZ - bulkV[2]); //Velocity in the cell, in the plasma frame
                  
                  Vec3d NeighbourVcoord(VX,VY,VZ); //Coordinates of Neighbour Cell
                  
                  const Real DV 
                     = parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVX];
	
                  Real DVarray[3] = {0,0,0};
                  DVarray[coord] = DV;
                  Vec3d VDV;
                  VDV.load(DVarray);

                  NeighbourVcoord += VDV; 

                  Realf CellValue      = cell.get_value(VX,VY,VZ,popID);
                  Realf NeighbourValue = cell.get_value(NeighbourVcoord[0],NeighbourVcoord[1],NeighbourVcoord[2],popID);

                  Realf dfdcoord = (NeighbourValue - CellValue)/DV;
                  arraydf[WID3*n+i+WID*j+WID*WID*k][coord] = dfdcoord;
               }

            } 

        }


	for (int coord = 0; coord < 3; coord++) { // Second derivative loop

           SpatialCell& cell = *mpiGrid[CellID];

	   Vec3d B(cell.parameters[CellParams::PERBXVOL] +  cell.parameters[CellParams::BGBXVOL],
                   cell.parameters[CellParams::PERBYVOL] +  cell.parameters[CellParams::BGBYVOL],
	           cell.parameters[CellParams::PERBZVOL] +  cell.parameters[CellParams::BGBZVOL]);
           Vec3d b = normalize_vector(B);


	   const Real* parameters  = cell.get_block_parameters(popID);

           const vmesh::LocalID* nBlocks = cell.get_velocity_grid_length(popID);

            for (vmesh::LocalID n=0; n<cell.get_number_of_velocity_blocks(popID); n++) { //Iterate through velocity blocks
               for (uint k = 0; k < WID; ++k) for (uint j = 0; j < WID; ++j) for (uint i = 0; i < WID; ++i) {

                  //Get velocity space coordinates                    
	          const Real VX 
                     =          parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VXCRD] 
                     + (i + 0.5)*parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVX];
                  const Real VY 
                     =          parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VYCRD] 
                     + (j + 0.5)*parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVY];
                  const Real VZ 
                     =          parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VZCRD]
                     + (k + 0.5)*parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVZ];
                  
                  Vec3d bulkV(cell.parameters[CellParams::VX], cell.parameters[CellParams::VY], cell.parameters[CellParams::VZ]);

                  Vec3d V(VX - bulkV[0] , VY - bulkV[1] , VZ - bulkV[2]); //Velocity in the cell, in the plasma frame
                  Realf normV = sqrt(dot_product(V,V));
                                    
                  Vec3d NeighbourVcoord(VX,VY,VZ); //Coordinates of Neighbour Cell
                  
                  const Real DV 
                     = parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVX];
	
                  Real DVarray[3] = {0,0,0};
                  DVarray[coord] = DV;
                  Vec3d VDV;
                  VDV.load(DVarray);

                  NeighbourVcoord -= VDV; 
                  auto LeftNGID = cell.get_velocity_block(popID,NeighbourVcoord[0],NeighbourVcoord[1],NeighbourVcoord[2],0);
                  auto LeftNLID = cell.get_velocity_block_local_id(LeftNGID,popID);

                  int Lefti = (i - (coord == 0)? 1:0)%WID; //Looking for the left cell neighbour which could be in another block
                  int Leftj = (j - (coord == 1)? 1:0)%WID;
                  int Leftk = (k - (coord == 2)? 1:0)%WID;

                  Vec3d leftV  = V-1.0/2.0*VDV;
                  Vec3d rightV = V+1.0/2.0*VDV;
                  
                  Realf normVright = sqrt(dot_product(rightV,rightV));
                  Realf normVleft  = sqrt(dot_product(leftV,leftV));

                  Realf ddcoordleft     = 0.0;
                  Realf ddcoordright    = 0.0;
                  Realf precoeffyzleft  = 0.0;
                  Realf precoeffyzright = 0.0;
                  Realf precoeffx       = 0.0;

                  if (coord == 0){
                      precoeffx    = (V[1]*V[1] + V[2]*V[2]);
                      if (LeftNGID == vmesh::INVALID_GLOBALID) {ddcoordleft = 0.0;}
                      else {ddcoordleft  = precoeffx/normVleft * arraydf[WID3*LeftNLID+Lefti+WID*Leftj+WID*WID*Leftk][coord];}
                      ddcoordright = precoeffx/normVright * arraydf[WID3*n+i+WID*j+WID*WID*k][coord];
                  } else{    
                      precoeffyzleft  = sqrt(leftV[1]*leftV[1] + leftV[2]*leftV[2])/normVleft;
                      precoeffyzright = sqrt(rightV[1]*rightV[1] + rightV[2]*rightV[2])/normVright;
                      if (LeftNGID == vmesh::INVALID_GLOBALID) {ddcoordleft = 0.0;}
                      else{ddcoordleft     = precoeffyzleft * arraydf[WID3*LeftNLID+Lefti+WID*Leftj+WID*WID*Leftk][coord];}
                      ddcoordright    = precoeffyzright * arraydf[WID3*n+i+WID*j+WID*WID*k][coord];
                  }

                  Realf Dvv = Parameters::PADcoefficient;
                  Realf dt  = Parameters::dt;

                  Realf ddv = Dvv * (ddcoordright - ddcoordleft)/DV; // Second derivative (left and right so it is centered on cell)

                  Realf termdcoord = 0.0;
                  if (coord == 0){
                      termdcoord = normV * ddv;
                  } else{
                      termdcoord = normV * V[0]*V[0] / (2.0 * sqrt(V[1]*V[1] + V[2]*V[2])) * ddv; 
                  
                  }
                  
                  Realf CellValue = cell.get_value(VX,VY,VZ,popID);
                  CellValue = CellValue + termdcoord * dt ;
                  if (CellValue <= 0.0) { CellValue = 0.0;}

                  cell.set_value(VX,VY,VZ,CellValue,popID);

           } 

	 }
       }
    }
}


