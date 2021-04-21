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

        std::vector<std::array<Realf,3>> arrayDFright(cell.get_number_of_velocity_blocks(popID)*WID3); // Array of vspace size for storing derivatives +DV
        std::vector<std::array<Realf,3>> arrayDFleft(cell.get_number_of_velocity_blocks(popID)*WID3); // Array of vspace size for storing derivatives -DV


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
                  
                  Vec3d V(VX,VY,VZ); // Velocity in the cell, in the simulation frame
                                  
                  const Real DV 
                     = parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVX];
	
                  Real DVarray[3] = {0,0,0};
                  DVarray[coord] = DV;
                  Vec3d VDV;
                  VDV.load(DVarray);

                  Vec3d NeighbourVright = V + VDV; // Cell coordinates in +DV direction
                  Vec3d NeighbourVleft  = V - VDV; // Cell coordinates in -DV direction
 
                  // f values for center, +DV and -DV (= 0 if cell doesnt exist)
                  Realf CellValue      = cell.get_value(VX,VY,VZ,popID);
                  Realf CellValueRight = cell.get_value(NeighbourVright[0],NeighbourVright[1],NeighbourVright[2],popID);
                  Realf CellValueLeft  = cell.get_value(NeighbourVleft[0],NeighbourVleft[1],NeighbourVleft[2],popID);

                  // First derivatives on right and left faces of center cell
                  Realf dfdcoordRight = (CellValueRight - CellValue)/DV;
                  Realf dfdcoordLeft  = (CellValue - CellValueLeft)/DV;

                  arrayDFright[WID3*n+i+WID*j+WID*WID*k][coord] = dfdcoordRight;
                  arrayDFleft[WID3*n+i+WID*j+WID*WID*k][coord]  = dfdcoordLeft;
                     
               }

            } 

        }


	for (int coord = 0; coord < 3; coord++) { // Second derivative loop

           SpatialCell& cell = *mpiGrid[CellID];

	   Vec3d B(cell.parameters[CellParams::PERBXVOL] +  cell.parameters[CellParams::BGBXVOL],
                   cell.parameters[CellParams::PERBYVOL] +  cell.parameters[CellParams::BGBYVOL],
	           cell.parameters[CellParams::PERBZVOL] +  cell.parameters[CellParams::BGBZVOL]);
           Vec3d b = normalize_vector(B);
           Vec3d ey(0.0,1.0,0.0);
           Vec3d evec;
           if (dot_product(b,ey) > 0.1) {Vec3d evec(0.0,1.0,0.0);}
           else {Vec3d evec(0.0,0.0,1.0);} 
           Vec3d c = normalize_vector(cross_product(b,evec));
           Vec3d d = normalize_vector(cross_product(b,c));

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
                  Vec3d V(VX,VY,VZ); // Velocity in the cell, in the simulation frame
                  Vec3d Vplasma(VX - bulkV[0] , VY - bulkV[1] , VZ - bulkV[2]); //Velocity in the cell, in the plasma frame
                  
                  Realf normV = sqrt(dot_product(Vplasma,Vplasma));

                  const Real DV 
                     = parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVX];
	
                  Real DVarray[3] = {0,0,0};
                  DVarray[coord] = DV;
                  Vec3d VDV;
                  VDV.load(DVarray);

                  Realf Dvv = Parameters::PADcoefficient; // Diffusion coefficient taken from cfg file

                  // Calculation of theta at center of the cell
                  Vec3d r     = normalize_vector(Vplasma);
                  Realf rc    = dot_product(r,c);
                  Realf rd    = dot_product(r,d);
                  Realf theta = atan2(rc,rd);  

                  // Calculation of terms inside the second derivative according to Eq. (18) of the PDF
                      // Right terms

                  Vec3d rightVplasma = Vplasma + 1.0/2.0*VDV; // Velocity at the right face of the cell, in the plasma frame
                  Realf normVright   = sqrt(dot_product(rightVplasma,rightVplasma));

                  Realf rightTermDVX = sqrt(rightVplasma[1]*rightVplasma[1] + rightVplasma[2]*rightVplasma[2]) * arrayDFright[WID3*n+i+WID*j+WID*WID*k][0];
                  Realf rightTermDVY = rightVplasma[0] * sin(theta) * arrayDFright[WID3*n+i+WID*j+WID*WID*k][1];
                  Realf rightTermDVZ = rightVplasma[0] * cos(theta) * arrayDFright[WID3*n+i+WID*j+WID*WID*k][2];
            
                  Realf rightTerm = sqrt(rightVplasma[1]*rightVplasma[1] + rightVplasma[2]*rightVplasma[2])/normVright * Dvv * (rightTermDVY + rightTermDVZ - rightTermDVX);

                      // Left terms
                  
                  Vec3d leftVplasma = Vplasma - 1.0/2.0*VDV; // Velocity at the right face of the cell, in the plasma frame
                  Realf normVleft   = sqrt(dot_product(leftVplasma,leftVplasma));

                  Realf leftTermDVX = sqrt(leftVplasma[1]*leftVplasma[1] + leftVplasma[2]*leftVplasma[2]) * arrayDFleft[WID3*n+i+WID*j+WID*WID*k][0];
                  Realf leftTermDVY = leftVplasma[0] * sin(theta) * arrayDFleft[WID3*n+i+WID*j+WID*WID*k][1];
                  Realf leftTermDVZ = leftVplasma[0] * cos(theta) * arrayDFleft[WID3*n+i+WID*j+WID*WID*k][2];
            
                  Realf leftTerm = sqrt(leftVplasma[1]*leftVplasma[1] + leftVplasma[2]*leftVplasma[2])/normVleft * Dvv * (leftTermDVY + leftTermDVZ - leftTermDVX);
                  
                  // Second derivative (centered difference of left and right sides)

                  Realf precoeff = 0.0;
                  if (coord == 0) { precoeff = - normV;}
                  else if (coord == 1) {precoeff = normV * Vplasma[0] * sin(theta) / sqrt(Vplasma[1]*Vplasma[1] + Vplasma[2]*Vplasma[2]);}
                  else if (coord == 2) {precoeff = normV * Vplasma[0] * cos(theta) / sqrt(Vplasma[1]*Vplasma[1] + Vplasma[2]*Vplasma[2]);} 

                  Realf dfdtCoord = precoeff * (rightTerm - leftTerm)/DV; 

                  // Update cell

                  Realf dt = Parameters::dt; // Simulation time step

                  Realf CellValue = cell.get_value(VX,VY,VZ,popID);
                  CellValue = CellValue + dfdtCoord * dt ;
                  if (CellValue <= 0.0) { CellValue = 0.0;}

                  cell.set_value(VX,VY,VZ,CellValue,popID);

           } 

	 }
       }
    }
}


