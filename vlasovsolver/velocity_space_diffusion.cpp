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
	const Real* parameters  = cell.get_block_parameters(popID);
        const vmesh::LocalID* nBlocks = cell.get_velocity_grid_length(popID);
        
        Realf Sparsity    = cell.getVelocityBlockMinValue(popID);

        Realf dtTotalDiff = 0.0; // Diffusion time elapsed

        while (dtTotalDiff < Parameters::dt) {

            Realf RemainT = Parameters::dt - dtTotalDiff; //Remaining time before reaching simulation time step

            std::vector<std::array<Realf,3>> arrayDFdcoordFirst(cell.get_number_of_velocity_blocks(popID)*WID3); // Array of vspace size for storing first derivatives at center of cell
            std::vector<std::array<Realf,3>> arrayDFdcoordSecond(cell.get_number_of_velocity_blocks(popID)*WID3); // Array of vspace size for storing second derivatives at center of cell
            std::vector<Realf> arrayDFdVXdVY(cell.get_number_of_velocity_blocks(popID)*WID3); // Array of vspace size to store mix derivative dvxdvy
            std::vector<Realf> arrayDFdVXdVZ(cell.get_number_of_velocity_blocks(popID)*WID3); // Array of vspace size to store mix derivative dvxdvz
            std::vector<Realf> arrayDFdVYdVZ(cell.get_number_of_velocity_blocks(popID)*WID3); // Array of vspace size to store mix derivative dvydvz
            std::vector<Realf> dfdt(cell.get_number_of_velocity_blocks(popID)*WID3); // Array of vspace size to store dfdt
            std::vector<Realf> checkCFL(cell.get_number_of_velocity_blocks(popID)*WID3); // Array of vspace size to store checkCFl

            for (int coord = 0; coord < 3; coord++) { // Loop for first and second derivatives per coordinates
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
                      Realf CellValue      = cell.get_value(VX,VY,VZ,popID);
                      Realf CellValueRight = cell.get_value(NeighbourVright[0],NeighbourVright[1],NeighbourVright[2],popID);
                      Realf CellValueLeft  = cell.get_value(NeighbourVleft[0],NeighbourVleft[1],NeighbourVleft[2],popID);


                      // First derivatives on right and left faces of cell
                      Realf dfdcoordRight = (CellValueRight - CellValue)/DV;
                      Realf dfdcoordLeft  = (CellValue - CellValueLeft)/DV;
     
                      // First derivatives at the center of the cell
                      arrayDFdcoordFirst[WID3*n+i+WID*j+WID*WID*k][coord] = (dfdcoordRight + dfdcoordLeft)/2.0; 
                    
                      // Second derivatives at the center of the cell
                      arrayDFdcoordSecond[WID3*n+i+WID*j+WID*WID*k][coord] = (CellValueRight - 2*CellValue + CellValueLeft)/(DV*DV);     

                   } // End cell loop

                } // End block loop 

            } // End coord loop and derivatives

            for (vmesh::LocalID n=0; n<cell.get_number_of_velocity_blocks(popID); n++) { // Loop for mixed derivatives
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
	    
                    // 19 cell values
                    Realf CellValue     = cell.get_value(VX,VY,VZ,popID);
                    Realf CellValuePX   = cell.get_value(VX+DV,VY,   VZ,   popID);
                    Realf CellValueMX   = cell.get_value(VX-DV,VY,   VZ,   popID);
                    Realf CellValuePY   = cell.get_value(VX,   VY+DV,VZ,   popID);
                    Realf CellValueMY   = cell.get_value(VX,   VY-DV,VZ,   popID);
                    Realf CellValuePZ   = cell.get_value(VX,   VY,   VZ+DV,popID);
                    Realf CellValueMZ   = cell.get_value(VX,   VY,   VZ-DV,popID);
                    Realf CellValuePXPY = cell.get_value(VX+DV,VY+DV,VZ,   popID);
                    Realf CellValuePXMY = cell.get_value(VX+DV,VY-DV,VZ,   popID);                 
                    Realf CellValueMXPY = cell.get_value(VX-DV,VY+DV,VZ,   popID);
                    Realf CellValueMXMY = cell.get_value(VX-DV,VY-DV,VZ,   popID);
                    Realf CellValuePXPZ = cell.get_value(VX+DV,VY,   VZ+DV,popID);
                    Realf CellValuePXMZ = cell.get_value(VX+DV,VY,   VZ-DV,popID); 
                    Realf CellValueMXPZ = cell.get_value(VX-DV,VY,   VZ+DV,popID);
                    Realf CellValueMXMZ = cell.get_value(VX-DV,VY,   VZ-DV,popID); 
                    Realf CellValuePYPZ = cell.get_value(VX,   VY+DV,VZ+DV,popID);
                    Realf CellValuePYMZ = cell.get_value(VX,   VY+DV,VZ-DV,popID);
                    Realf CellValueMYPZ = cell.get_value(VX,   VY-DV,VZ+DV,popID);
                    Realf CellValueMYMZ = cell.get_value(VX,   VY-DV,VZ-DV,popID);
                                
                    
                    //First derivatives
                        //dvx        
                    Realf dfdxPX     = (CellValuePX - CellValue)/DV;
                    Realf dfdxMX     = (CellValue - CellValueMX)/DV;
                    Realf dfdxCenter = (dfdxPX + dfdxMX)/2.0;
                        //dvy
                    Realf dfdyPY     = (CellValuePY - CellValue)/DV;
                    Realf dfdyMY     = (CellValue - CellValueMY)/DV;
                    Realf dfdyCenter = (dfdyPY + dfdyMY)/2.0;
                        //dvz
                    Realf dfdzPZ     = (CellValuePZ - CellValue)/DV;
                    Realf dfdzMZ     = (CellValue - CellValueMZ)/DV;
                    Realf dfdzCenter = (dfdzPZ + dfdzMZ)/2.0;


                    // d2f/dvxdvy
                        //dvx first
                    Realf dfdxPXPY   = (CellValuePXPY - CellValuePY)/DV;
                    Realf dfdxMXPY   = (CellValuePY - CellValueMXPY)/DV;
                    Realf dfdxPY     = (dfdxPXPY + dfdxMXPY)/2.0;
 
                    Realf dfdxPXMY   = (CellValuePXMY - CellValueMY)/DV;
                    Realf dfdxMXMY   = (CellValueMY - CellValueMXMY)/DV;
                    Realf dfdxMY     = (dfdxPXMY + dfdxMXMY)/2.0;

                    Realf dfdxdyPY = (dfdxPY - dfdxCenter)/DV;
                    Realf dfdxdyMY = (dfdxCenter - dfdxMY)/DV;
                    Realf dfdxdy   = (dfdxdyPY + dfdxdyMY)/2.0;
                        //dvy first
                    Realf dfdyPXPY   = (CellValuePXPY - CellValuePX)/DV;
                    Realf dfdyPXMY   = (CellValuePX - CellValuePXMY)/DV;
                    Realf dfdyPX     = (dfdyPXPY + dfdyPXMY)/2.0;
 
                    Realf dfdyMXPY   = (CellValueMXPY - CellValueMX)/DV;
                    Realf dfdyMXMY   = (CellValueMX - CellValueMXMY)/DV;
                    Realf dfdyMX     = (dfdyMXPY + dfdyMXMY)/2.0;

                    Realf dfdydxPX = (dfdyPX - dfdyCenter)/DV;
                    Realf dfdydxMX = (dfdyCenter - dfdyMX)/DV;
                    Realf dfdydx   = (dfdydxPX + dfdydxMX)/2.0;
                        //average
                    arrayDFdVXdVY[WID3*n+i+WID*j+WID*WID*k] = (dfdxdy + dfdydx)/2.0;
   

                    // d2f/dvxdvz
                        //dvx first
                    Realf dfdxPXPZ   = (CellValuePXPZ - CellValuePZ)/DV;
                    Realf dfdxMXPZ   = (CellValuePZ - CellValueMXPZ)/DV;
                    Realf dfdxPZ     = (dfdxPXPZ + dfdxMXPZ)/2.0;
 
                    Realf dfdxPXMZ   = (CellValuePXMZ - CellValueMZ)/DV;
                    Realf dfdxMXMZ   = (CellValueMZ - CellValueMXMZ)/DV;
                    Realf dfdxMZ     = (dfdxPXMZ + dfdxMXMZ)/2.0;

                    Realf dfdxdzPZ = (dfdxPZ - dfdxCenter)/DV;
                    Realf dfdxdzMZ = (dfdxCenter - dfdxMZ)/DV;
                    Realf dfdxdz   = (dfdxdzPZ + dfdxdzMZ)/2.0;
                        //dvz first
                    Realf dfdzPXPZ   = (CellValuePXPZ - CellValuePX)/DV;
                    Realf dfdzPXMZ   = (CellValuePX - CellValuePXMZ)/DV;
                    Realf dfdzPX     = (dfdzPXPZ + dfdzPXMZ)/2.0;
 
                    Realf dfdzMXPZ   = (CellValueMXPZ - CellValueMX)/DV;
                    Realf dfdzMXMZ   = (CellValueMX - CellValueMXMZ)/DV;
                    Realf dfdzMX     = (dfdzMXPZ + dfdzMXMZ)/2.0;

                    Realf dfdzdxPX = (dfdzPX - dfdzCenter)/DV;
                    Realf dfdzdxMX = (dfdzCenter - dfdzMX)/DV;
                    Realf dfdzdx   = (dfdzdxPX + dfdzdxMX)/2.0;
                        //average
	            arrayDFdVXdVZ[WID3*n+i+WID*j+WID*WID*k] = (dfdxdz + dfdzdx)/2.0;


                    // d2f/dvydvz
                        //dvy first
                    Realf dfdyPYPZ   = (CellValuePYPZ - CellValuePZ)/DV;
                    Realf dfdyMYPZ   = (CellValuePZ - CellValueMYPZ)/DV;
                    Realf dfdyPZ     = (dfdyPYPZ + dfdyMYPZ)/2.0;
 
                    Realf dfdyPYMZ   = (CellValuePYMZ - CellValueMZ)/DV;
                    Realf dfdyMYMZ   = (CellValueMZ - CellValueMYMZ)/DV;
                    Realf dfdyMZ     = (dfdyPYMZ + dfdyMYMZ)/2.0;

                    Realf dfdydzPZ = (dfdyPZ - dfdyCenter)/DV;
                    Realf dfdydzMZ = (dfdyCenter - dfdyMZ)/DV;
                    Realf dfdydz   = (dfdydzPZ + dfdydzMZ)/2.0;
                        //dvz first
                    Realf dfdzPYPZ   = (CellValuePYPZ - CellValuePY)/DV;
                    Realf dfdzPYMZ   = (CellValuePY - CellValuePYMZ)/DV;
                    Realf dfdzPY     = (dfdzPYPZ + dfdzPYMZ)/2.0;
 
                    Realf dfdzMYPZ   = (CellValueMYPZ - CellValueMY)/DV;
                    Realf dfdzMYMZ   = (CellValueMY - CellValueMYMZ)/DV;
                    Realf dfdzMY     = (dfdzMYPZ + dfdzMYMZ)/2.0;

                    Realf dfdzdyPY = (dfdzPY - dfdzCenter)/DV;
                    Realf dfdzdyMY = (dfdzCenter - dfdzMY)/DV;
                    Realf dfdzdy   = (dfdzdyPY + dfdzdyMY)/2.0;
                        //average
	            arrayDFdVYdVZ[WID3*n+i+WID*j+WID*WID*k] = (dfdydz + dfdzdy)/2.0;

                } // End cell loop

            } // End block loop and mixed derivatives

	    // dfdt loop
            SpatialCell& cell = *mpiGrid[CellID];

	    std::vector<Realf> B = {cell.parameters[CellParams::PERBXVOL] +  cell.parameters[CellParams::BGBXVOL],
                                    cell.parameters[CellParams::PERBYVOL] +  cell.parameters[CellParams::BGBYVOL],
	                            cell.parameters[CellParams::PERBZVOL] +  cell.parameters[CellParams::BGBZVOL]};
            Realf normB = sqrt(B.at(0)*B.at(0) + B.at(1)*B.at(1) + B.at(2)*B.at(2));     
            std::vector<Realf> b = {B.at(0)/normB, B.at(1)/normB, B.at(2)/normB};
            
            // Building orthonormal basis for phi calculation
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

                   Realf normV = sqrt(Vplasma.at(0)*Vplasma.at(0) + Vplasma.at(1)*Vplasma.at(1) + Vplasma.at(2)*Vplasma.at(2));

                   const Real DV 
                      = parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVX];
	    
                   // Calculation of phi at center of the cell
                   std::array<Realf,3> r  = {Vplasma.at(0)/normV, Vplasma.at(1)/normV, Vplasma.at(2)/normV};
                   Realf rc = r.at(0)*c.at(0) + r.at(1)*c.at(1) + r.at(2)*c.at(2);
                   Realf rd = r.at(0)*d.at(0) + r.at(1)*d.at(1) + r.at(2)*d.at(2);
                   //phi[WID3*n+i+WID*j+WID*WID*k] = atan2(rc,rd);  
                   Realf phi = atan2(rc,rd);

                   // dfdt according to Eq. 28 of the PDF
                   dfdt[WID3*n+i+WID*j+WID*WID*k] = Parameters::PADcoefficient * 
                   (-2.0 * arrayDFdcoordFirst[WID3*n+i+WID*j+WID*WID*k][0] + (2.0 * Vplasma[0]*Vplasma[0] - normV*normV) / sqrt(Vplasma[1]*Vplasma[1] + Vplasma[2]*Vplasma[2]) * (sin(phi) * arrayDFdcoordFirst[WID3*n+i+WID*j+WID*WID*k][1] + cos(phi) * arrayDFdcoordFirst[WID3*n+i+WID*j+WID*WID*k][2])
                     - 2.0 * Vplasma[0] * sqrt(Vplasma[1]*Vplasma[1] + Vplasma[2]*Vplasma[2]) * (sin(phi) * arrayDFdVXdVY[WID3*n+i+WID*j+WID*WID*k] + cos(phi) * arrayDFdVXdVZ[WID3*n+i+WID*j+WID*WID*k]) + 2.0 * Vplasma[0]*Vplasma[0] * Vplasma[2] * Vplasma[1] / (Vplasma[2]*Vplasma[2] + Vplasma[1]*Vplasma[1]) * arrayDFdVYdVZ[WID3*n+i+WID*j+WID*WID*k] 
                     + (Vplasma[1]*Vplasma[1] + Vplasma[2]*Vplasma[2]) * arrayDFdcoordSecond[WID3*n+i+WID*j+WID*WID*k][0] + Vplasma[0]*Vplasma[0] * Vplasma[2]*Vplasma[2] / (Vplasma[2]*Vplasma[2] + Vplasma[1]*Vplasma[1]) * arrayDFdcoordSecond[WID3*n+i+WID*j+WID*WID*k][2] + Vplasma[0]*Vplasma[0] * Vplasma[1]*Vplasma[1] / (Vplasma[2]*Vplasma[2] + Vplasma[1]*Vplasma[1]) * arrayDFdcoordSecond[WID3*n+i+WID*j+WID*WID*k][1]);
                   // CFL condition
                   Realf CellValue = cell.get_value(VX,VY,VZ,popID);
                   if (CellValue == 0.0) {CellValue = Sparsity;} //Set CellValue to sparsity Threshold for empty cells otherwise div by 0
                   checkCFL[WID3*n+i+WID*j+WID*WID*k] = CellValue * Parameters::PADCFL * (1.0 / abs(dfdt[WID3*n+i+WID*j+WID*WID*k]));

                } // End cell loop

            } // End block and dfdt loop

            //Calculate Diffusion time step based on min of CFL condition  
            Realf mincheckCFL = *min_element(checkCFL.begin(),checkCFL.end());
            Realf Ddt = mincheckCFL; // Diffusion time step
            if (Ddt > RemainT) { Ddt = RemainT; }
            std::cout << "Diffusion dt = " << Ddt << std::endl;
            dtTotalDiff = dtTotalDiff + Ddt;

            //Loop to check CFL and update cell
            for (vmesh::LocalID n=0; n<cell.get_number_of_velocity_blocks(popID); n++) { //Iterate through velocity blocks
                for (uint k = 0; k < WID; ++k) for (uint j = 0; j < WID; ++j) for (uint i = 0; i < WID; ++i) {
                    const Real* parameters  = cell.get_block_parameters(popID);

                    //Get velocity space coordinates                   
                    const Real VX = parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VXCRD] + (i + 0.5)*parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVX];
                    const Real VY = parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VYCRD] + (j + 0.5)*parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVY];
                    const Real VZ = parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VZCRD] + (k + 0.5)*parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVZ];

                    //Check CFL
                    Realf CellValue = cell.get_value(VX,VY,VZ,popID);

                    //Update cell
                    Realf NewCellValue = CellValue + dfdt[WID3*n+i+WID*j+WID*WID*k] * Ddt;
                    if (NewCellValue <= 0.0) { NewCellValue = 0.0;}

                    cell.set_value(VX,VY,VZ,NewCellValue,popID);

                    
                } // End cell loop

            } // End block loop

        } // End Time loop

    } // End spatial cell loop

} // End function


