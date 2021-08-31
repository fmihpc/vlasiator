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
#include <iostream>
#include <fstream>
#include <iomanip>

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
        //std::vector<Realf> theta(cell.get_number_of_velocity_blocks(popID)*WID3); // Array of vspace size for storing theta

        std::vector<Realf> dfdt(cell.get_number_of_velocity_blocks(popID)*WID3); // Array of vspace size to store dfdt

        for (int coord = 0; coord < 3; coord++) { // First derivative loop

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
                  
                  std::vector<Realf> V = {VX,VY,VZ}; // Velocity in the cell, in the simulation frame
                                  
                  const Real DV 
                     = parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVX];
	
                  std::vector<Realf> VDV = {0.0,0.0,0.0};
                  VDV.at(coord) = DV;

                  std::vector<Realf> NeighbourVright; // Cell coordinates in +DV direction
                  for (int indx = 0; indx < V.size(); indx++) {
                      NeighbourVright.push_back(V.at(indx) + VDV.at(indx));
                  }

                  std::vector<Realf> NeighbourVleft; // Cell coordinates in -DV direction
                  for (int indx = 0; indx < V.size(); indx++) {
                      NeighbourVleft.push_back(V.at(indx) - VDV.at(indx));
                  }

                  // f values for center, +DV and -DV (= 0 if cell doesnt exist)
                  // Take log10(f) to minimize oscillations
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

        } // End first derivative 

	for (int coord = 0; coord < 3; coord++) { // Second derivative loop

            SpatialCell& cell = *mpiGrid[CellID];

	    std::vector<Realf> B = {cell.parameters[CellParams::PERBXVOL] +  cell.parameters[CellParams::BGBXVOL],
                                    cell.parameters[CellParams::PERBYVOL] +  cell.parameters[CellParams::BGBYVOL],
	                            cell.parameters[CellParams::PERBZVOL] +  cell.parameters[CellParams::BGBZVOL]};
            Realf normB = sqrt(B.at(0)*B.at(0) + B.at(1)*B.at(1) + B.at(2)*B.at(2));     
            std::vector<Realf> b = {B.at(0)/normB, B.at(1)/normB, B.at(2)/normB};

            std::array<Realf,3> evec = {0.0,1.0,0.0};
            std::array<Realf,3> evecarray = {0.0,0.0,0.0};
            Realf dotProd = b[0]*evec[0] + b[1]*evec[1] + b[2]*evec[2];
            if (abs(dotProd) < 0.1) { evecarray.at(1) = 1.0;}
            else  {evecarray.at(2) = 1.0;}

            std::array<Realf,3> cvec = {b.at(1) * evecarray.at(2) - b.at(2) * evecarray.at(1),
                                       b.at(2) * evecarray.at(0) - b.at(0) * evecarray.at(2),
                                       b.at(0) * evecarray.at(1) - b.at(1) * evecarray.at(0)};             // cvec = b.evecarray
            Realf cvecNorm = sqrt(cvec.at(0)*cvec.at(0) + cvec.at(1)*cvec.at(1) + cvec.at(2)*cvec.at(2));  
            std::array<Realf,3> c = {cvec.at(0)/cvecNorm, cvec.at(1)/cvecNorm, cvec.at(2)/cvecNorm};

            std::array<Realf,3> dvec = {b.at(1) * c.at(2) - b.at(2) * c.at(1),
                                       b.at(2) * c.at(0) - b.at(0) * c.at(2),
                                       b.at(0) * c.at(1) - b.at(1) * c.at(0)};                             // dvec = b.c
            Realf dvecNorm = sqrt(dvec.at(0)*dvec.at(0) + dvec.at(1)*dvec.at(1) + dvec.at(2)*dvec.at(2));
            std::array<Realf,3> d = {dvec.at(0)/dvecNorm, dvec.at(1)/dvecNorm, dvec.at(2)/dvecNorm};

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
                   
                   std::vector<Realf> bulkV = {cell.parameters[CellParams::VX], cell.parameters[CellParams::VY], cell.parameters[CellParams::VZ]};
                   std::vector<Realf> V     = {VX,VY,VZ}; // Velocity in the cell, in the simulation frame
                   std::vector<Realf> Vplasma;            //Velocity in the cell, in the plasma frame
                   for (int indx = 0; indx < V.size(); indx++) {
                       Vplasma.push_back(V.at(indx) - bulkV.at(indx));
                   }

                   //VPCoords[WID3*n+i+WID*j+WID*WID*k][0] = V[0];
                   //VPCoords[WID3*n+i+WID*j+WID*WID*k][1] = V[1];
                   //VPCoords[WID3*n+i+WID*j+WID*WID*k][2] = V[2];

                   Realf normV = sqrt(Vplasma.at(0)*Vplasma.at(0) + Vplasma.at(1)*Vplasma.at(1) + Vplasma.at(2)*Vplasma.at(2));

                   const Real DV 
                      = parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVX];
	
                   std::vector<Realf> VDV = {0.0,0.0,0.0};
                   VDV.at(coord) = DV;

                   Realf Dvv = Parameters::PADcoefficient; // Diffusion coefficient taken from cfg file

                   // Calculation of theta at center of the cell
                   std::array<Realf,3> r  = {Vplasma.at(0)/normV, Vplasma.at(1)/normV, Vplasma.at(2)/normV};
                   Realf rc = r.at(0)*c.at(0) + r.at(1)*c.at(1) + r.at(2)*c.at(2);
                   Realf rd = r.at(0)*d.at(0) + r.at(1)*d.at(1) + r.at(2)*d.at(2);
                   //theta[WID3*n+i+WID*j+WID*WID*k] = atan2(rc,rd);  
                   Realf theta = atan2(rc,rd);

                   // Calculation of terms inside the second derivative according to Eq. (18) of the PDF
                       // Right terms

                   std::vector<Realf> rightVplasma; // Velocity at the right face of the cell, in the plasma frame
                   for (int indx = 0; indx < V.size(); indx++) {
                       rightVplasma.push_back(V.at(indx) + 1.0/2.0*VDV.at(indx));
                   }

                   Realf normVright   = sqrt(rightVplasma.at(0)*rightVplasma.at(0) + rightVplasma.at(1)*rightVplasma.at(1) + rightVplasma.at(2)*rightVplasma.at(2));

                   Realf rightTermDVX = sqrt(rightVplasma[1]*rightVplasma[1] + rightVplasma[2]*rightVplasma[2]) * arrayDFright[WID3*n+i+WID*j+WID*WID*k][0];
                   Realf rightTermDVY = rightVplasma[0] * sin(theta) * arrayDFright[WID3*n+i+WID*j+WID*WID*k][1];
                   Realf rightTermDVZ = rightVplasma[0] * cos(theta) * arrayDFright[WID3*n+i+WID*j+WID*WID*k][2];
            
                   Realf rightTerm = sqrt(rightVplasma[1]*rightVplasma[1] + rightVplasma[2]*rightVplasma[2])/normVright * Dvv * (rightTermDVY + rightTermDVZ - rightTermDVX);

                       // Left terms
                   
                   std::vector<Realf> leftVplasma; // Velocity at the left face of the cell, in the plasma frame
                   for (int indx = 0; indx < V.size(); indx++) {
                       leftVplasma.push_back(V.at(indx) - 1.0/2.0*VDV.at(indx));
                   }
                   
                   Realf normVleft   = sqrt(leftVplasma.at(0)*leftVplasma.at(0) + leftVplasma.at(1)*leftVplasma.at(1) + leftVplasma.at(2)*leftVplasma.at(2));

                   Realf leftTermDVX = sqrt(leftVplasma[1]*leftVplasma[1] + leftVplasma[2]*leftVplasma[2]) * arrayDFleft[WID3*n+i+WID*j+WID*WID*k][0];
                   Realf leftTermDVY = leftVplasma[0] * sin(theta) * arrayDFleft[WID3*n+i+WID*j+WID*WID*k][1];
                   Realf leftTermDVZ = leftVplasma[0] * cos(theta) * arrayDFleft[WID3*n+i+WID*j+WID*WID*k][2];
            
                   Realf leftTerm = sqrt(leftVplasma[1]*leftVplasma[1] + leftVplasma[2]*leftVplasma[2])/normVleft * Dvv * (leftTermDVY + leftTermDVZ - leftTermDVX);
                   
                   // Second derivative (centered difference of left and right sides)

                   Realf precoeff = 0.0;
                   if (coord == 0) { precoeff = - normV;}
                   else if (coord == 1) {precoeff = normV * Vplasma[0] * sin(theta) / sqrt(Vplasma[1]*Vplasma[1] + Vplasma[2]*Vplasma[2]);}
                   else if (coord == 2) {precoeff = normV * Vplasma[0] * cos(theta) / sqrt(Vplasma[1]*Vplasma[1] + Vplasma[2]*Vplasma[2]);} 

                   //dfdtCoord[WID3*n+i+WID*j+WID*WID*k][coord] = precoeff * (rightTerm - leftTerm)/DV; 
                   Realf dfdtCoord = precoeff * (rightTerm - leftTerm)/DV;

                   //Sum dfdtCoord
                   dfdt[WID3*n+i+WID*j+WID*WID*k] = dfdt[WID3*n+i+WID*j+WID*WID*k] + dfdtCoord;
           
                   // Update cell
                   Realf dt = Parameters::dt; // Simulation time step

                   Realf CellValue = cell.get_value(VX,VY,VZ,popID);
                   //CellValue = CellValue + dfdtCoord[WID3*n+i+WID*j+WID*WID*k][coord] * dt ;
                   Realf NewCellValue = CellValue + dfdtCoord * dt ;
                   if (NewCellValue <= 0.0) { NewCellValue = 0.0;}

                   cell.set_value(VX,VY,VZ,NewCellValue,popID);
              

                   }
            }
 
        }

        //Loop to check CFL and update cell
        for (vmesh::LocalID n=0; n<cell.get_number_of_velocity_blocks(popID); n++) { //Iterate through velocity blocks
            for (uint k = 0; k < WID; ++k) for (uint j = 0; j < WID; ++j) for (uint i = 0; i < WID; ++i) {

                //Get velocity space coordinates                   
                const Real VX = parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VXCRD] + (i + 0.5)*parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVX];
                const Real VY = parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VYCRD] + (j + 0.5)*parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVY];
                const Real VZ = parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VZCRD] + (k + 0.5)*parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVZ];

                //Check CFL
                Realf CellValue = cell.get_value(VX,VY,VZ,popID);
                checkCFL = dfdt[WID3*n+i+WID*j+WID*WID*k] / CellValue;

             }
         }


    }
}


