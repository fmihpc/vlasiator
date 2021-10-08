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
#include "../object_wrapper.h"
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
        const vmesh::MeshParameters& vMesh = getObjectWrapper().velocityMeshes[0];      

        Realf Sparsity    = cell.getVelocityBlockMinValue(popID);

        Realf dtTotalDiff = 0.0; // Diffusion time elapsed

        int nbins_v  = Parameters::PADvbins;
        int nbins_mu = Parameters::PADmubins;

        Realf mumin   = -1.0;
        Realf mumax   = +1.0;
        Realf dmubins = (mumax - mumin)/nbins_mu;

        Realf Vmin = vMesh.meshLimits[0];
        Realf Vmax = vMesh.meshLimits[1]; 
        std::cout << "Vmin, Vmax  = " << Vmin << ", " << Vmax << std::endl;
        Realf dVbins = (Vmax - Vmin)/(4.0*(nbins_v));  
        
        while (dtTotalDiff < Parameters::dt) {

            Realf RemainT = Parameters::dt - dtTotalDiff; //Remaining time before reaching simulation time step

            std::vector<Realf> dfdt(cell.get_number_of_velocity_blocks(popID)*WID3); // Array of vspace size to store dfdt
            std::vector<Realf> checkCFL(cell.get_number_of_velocity_blocks(popID)*WID3); // Array of vspace size to store checkCFl

            int fcount[nbins_v][nbins_mu];    // Array to count number of f stored
            Realf fmu[nbins_v][nbins_mu];     // Array to store f(v,mu)
            Realf dfdmu[nbins_v][nbins_mu];   // Array to store dfdmu
            Realf dfdt_mu[nbins_v][nbins_mu]; // Array to store dfdt_mu
            
            // Build 2d array of f(v,mu)
            for (vmesh::LocalID n=0; n<cell.get_number_of_velocity_blocks(popID); n++) { // Iterate through velocity blocks
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
  
                   std::vector<Realf> bulkV = {cell.parameters[CellParams::VX], cell.parameters[CellParams::VY], cell.parameters[CellParams::VZ]};
                   std::vector<Realf> Vplasma;            // Velocity in the cell, in the plasma frame
                   for (int indx = 0; indx < V.size(); indx++) {
                       Vplasma.push_back(V.at(indx) - bulkV.at(indx));
                   }


                   Realf normV = sqrt(Vplasma.at(0)*Vplasma.at(0) + Vplasma.at(1)*Vplasma.at(1) + Vplasma.at(2)*Vplasma.at(2));

                   Realf Vpara = Vplasma.at(0);
                   Realf Vperp = sqrt(Vplasma.at(1)*Vplasma.at(1) + Vplasma.at(2)*Vplasma.at(2));

                   Realf theta = atan2(Vperp,Vpara);
                   Realf mu    = cos(theta);
 
                   int Vcount = static_cast<int>(floor(normV / dVbins));                      

                   int mucount = static_cast<int>(floor(mu+1.0 / dmubins));                      

                   fcount[Vcount][mucount] += 1;

                   Realf CellValue      = cell.get_value(VX,VY,VZ,popID);
                   fmu[Vcount][mucount] += CellValue;

                }
            } // End blocks

            for (int indv = 0; indv < nbins_v; indv++) { // Divide f by count 
                for(int indmu = 0; indmu < nbins_mu; indmu++) {
                    fmu[indv][indmu] = fmu[indv][indmu]/fcount[indv][indmu]; 
                }
            }

            // Compute dfdmu
            for (int indv = 0; indv < nbins_v; indv++) { 
                for(int indmu = 0; indmu < nbins_mu; indmu++) {
                    if (indmu == 0) { dfdmu[indv][indmu] = (fmu[indv][indmu+1] - fmu[indv][indmu])/dmubins; }
                    else if (indmu == nbins_mu-1) { dfdmu[indv][indmu] = (fmu[indv][indmu] - fmu[indv][indmu-1])/dmubins; }
                    else { dfdmu[indv][indmu] = (fmu[indv][indmu+1] - fmu[indv][indmu-1])/(2.0*dmubins); } 
                }
            } 

            // Compute dfdt_mu
            for (int indv = 0; indv < nbins_v; indv++) { 
                for(int indmu = 0; indmu < nbins_mu; indmu++) {
                    if (indmu == 0) { dfdt_mu[indv][indmu] = Parameters::PADcoefficient * ( (1.0 - (3.0/2.0*dmubins - 1.0)*(3.0/2.0*dmubins - 1.0)) * dfdmu[indv][indmu+1] - (1.0 - (dmubins/2.0 - 1.0)*(dmubins/2.0 - 1.0)) * dfdmu[indv][indmu])/dmubins; }
                    else if (indmu == nbins_mu-1) { dfdt_mu[indv][indmu] = Parameters::PADcoefficient * ( (1.0 - (dmubins*(indmu+1.0/2.0)-1.0)*(dmubins*(indmu+1.0/2.0)-1.0)) * dfdmu[indv][indmu] - (1.0  - (dmubins*(indmu-1.0/2.0)-1.0)*(dmubins*(indmu-1.0/2.0)-1.0) * dfdmu[indv][indmu-1]) )/dmubins; }
                    else { dfdt_mu[indv][indmu] = Parameters::PADcoefficient * ( (1.0 - (dmubins*(indmu+3.0/2.0)-1.0)*(dmubins*(indmu+3.0/2.0)-1.0)) * dfdmu[indv][indmu+1] - (1.0 - (dmubins*(indmu-1.0/2.0)-1.0)*(dmubins*(indmu-1.0/2.0)-1.0)) * dfdmu[indv][indmu-1] )/(2.0*dmubins); } 
                }
            }

            // Compute dfdt
            for (vmesh::LocalID n=0; n<cell.get_number_of_velocity_blocks(popID); n++) { // Iterate through velocity blocks             
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
  
                   std::vector<Realf> bulkV = {cell.parameters[CellParams::VX], cell.parameters[CellParams::VY], cell.parameters[CellParams::VZ]};
                   std::vector<Realf> Vplasma;            // Velocity in the cell, in the plasma frame
                   for (int indx = 0; indx < V.size(); indx++) {
                       Vplasma.push_back(V.at(indx) - bulkV.at(indx));
                   }   
                
                   Realf normV = sqrt(Vplasma.at(0)*Vplasma.at(0) + Vplasma.at(1)*Vplasma.at(1) + Vplasma.at(2)*Vplasma.at(2));
                
                   Realf Vpara = Vplasma.at(0);
                   Realf Vperp = sqrt(Vplasma.at(1)*Vplasma.at(1) + Vplasma.at(2)*Vplasma.at(2));
  
                   Realf theta = atan2(Vperp,Vpara);
                   Realf mu    = cos(theta);

                   int Vcount = static_cast<int>(floor(normV / dVbins));    

                   int mucount = static_cast<int>(floor(mu+1.0 / dmubins));

                   Realf CellValue = cell.get_value(VX,VY,VZ,popID);
                   dfdt[WID3*n+i+WID*j+WID*WID*k] = dfdt_mu[Vcount][mucount] * CellValue / fmu[Vcount][mucount];

                   if (CellValue < Sparsity) {CellValue = Sparsity;} //Set CellValue to sparsity Threshold for empty cells otherwise div by 0
                   if (abs(dfdt[WID3*n+i+WID*j+WID*WID*k]) > 0.0) {
                   checkCFL[WID3*n+i+WID*j+WID*WID*k] = CellValue * Parameters::PADCFL * (1.0 / abs(dfdt[WID3*n+i+WID*j+WID*WID*k]));} else {
                   checkCFL[WID3*n+i+WID*j+WID*WID*k] = std::numeric_limits<Realf>::max();}
                }
            }

            //Calculate Diffusion time step based on min of CFL condition
            std::vector<Realf>::iterator mincheckCFL;
            mincheckCFL = std::min_element(checkCFL.begin(),checkCFL.end());
            if (mincheckCFL == checkCFL.end()) {break;}
            Realf Ddt = *mincheckCFL; // Diffusion time step
            if (Ddt > RemainT) { Ddt = RemainT; }
            //std::cout << "Diffusion dt = " << Ddt << std::endl;
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
               }
           }

        } // End Time loop

    } // End spatial cell loop

} // End function
