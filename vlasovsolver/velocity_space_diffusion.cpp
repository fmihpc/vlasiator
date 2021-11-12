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
#include <iterator>

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

        Realf Sparsity    = 0.01 * cell.getVelocityBlockMinValue(popID);

        Realf dtTotalDiff = 0.0; // Diffusion time elapsed

        int nbins_v  = Parameters::PADvbins;
        int nbins_mu = Parameters::PADmubins;

        Realf mumin   = -1.0;
        Realf mumax   = +1.0;
        Realf dmubins = (mumax - mumin)/nbins_mu;

        Realf Vmin = 0.0;
        Realf Vmax = 2*sqrt(3)*vMesh.meshLimits[1];
        Realf dVbins = (Vmax - Vmin)/nbins_v;  
        
        int subCount = 0; // Counter for substeps, used to print out. To be removed.

        while (dtTotalDiff < Parameters::dt) {

            Realf RemainT = Parameters::dt - dtTotalDiff; //Remaining time before reaching simulation time step

            std::vector<Realf> dfdt(cell.get_number_of_velocity_blocks(popID)*WID3); // Array of vspace size to store dfdt
            std::vector<Realf> checkCFL(cell.get_number_of_velocity_blocks(popID)*WID3, std::numeric_limits<Realf>::max()); // Array of vspace size to store checkCFl
            int fcount[nbins_v][nbins_mu]    = {0};   // Array to count number of f stored
            Realf fmu[nbins_v][nbins_mu]     = {0.0}; // Array to store f(v,mu)
            Realf dfdmu[nbins_v][nbins_mu]   = {0.0}; // Array to store dfdmu
            Realf dfdmu2[nbins_v][nbins_mu]  = {0.0}; // Array to store dfdmumu
            Realf dfdt_mu[nbins_v][nbins_mu] = {0.0}; // Array to store dfdt_mu
            std::vector<Realf> ratio(cell.get_number_of_velocity_blocks(popID)*WID3); // Array to store CellValue / fmu

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
 
                   int Vcount;
                   if (normV < Vmin) { continue; }
                   else { Vcount = static_cast<int>(floor((normV-Vmin) / dVbins)); }                      

                   int mucount = static_cast<int>(floor( (mu+1.0) / dmubins));                      

                   fcount[Vcount][mucount] += 1;

                   Realf CellValue      = cell.get_value(VX,VY,VZ,popID);
                   fmu[Vcount][mucount] += CellValue;

                }
            } // End blocks

            for (int indv = 0; indv < nbins_v; indv++) { // Divide f by count 
                for(int indmu = 0; indmu < nbins_mu; indmu++) {
                    if (fcount[indv][indmu] == 0) { fmu[indv][indmu] = 0.0;}
                    else {fmu[indv][indmu] = fmu[indv][indmu] / fcount[indv][indmu];} 
                }
            }
            
            if (subCount == 0) {
                // Save muspace to text
                std::string path_save = "/wrk/users/dubart/300_test/proc_test/mu_files/";
                std::ostringstream tmp;
                tmp << std::setw(7) << std::setfill('0') << P::tstep;
                std::string tstepString = tmp.str();
                std::ofstream muv_array(path_save + "muv_array_" + tstepString + ".txt");
                for (int indv = 0; indv < nbins_v; indv++) {
                    for(int indmu = 0; indmu < nbins_mu; indmu++) {
                        muv_array << fmu[indv][indmu] << ' ';
                    }
                    muv_array << std::endl;
                }
            }

            int cRight;
            int cLeft;

            // Compute dfdmu and dfdmu2 (take first non-zero neighbours)
            for (int indv = 0; indv < nbins_v; indv++) { 
                for(int indmu = 0; indmu < nbins_mu; indmu++) {
                    if (indmu == 0) {
                        cLeft  = 0;
                        cRight = 1;
                        while( (fcount[indv][indmu + cRight] == 0) && (indmu + cRight < nbins_mu-1) ) { cRight += 1; }
                        if( (fcount[indv][indmu + cRight] == 0) && (indmu + cRight == nbins_mu-1) ) { cRight = 0;}
                    } else if (indmu == nbins_mu-1) {
                        cLeft  = 1;
                        cRight = 0;
                        while( (fcount[indv][indmu - cLeft] == 0) && (indmu - cLeft > 0) ) { cLeft += 1; }
                        if( (fcount[indv][indmu - cLeft] == 0) && (indmu - cLeft == 0) ) { cLeft = 0;}
                    } else {
                        cLeft  = 1;
                        cRight = 1;
                        while( (fcount[indv][indmu + cRight] == 0) && (indmu + cRight < nbins_mu-1) ) { cRight += 1; }
                        if( (fcount[indv][indmu + cRight] == 0) && (indmu + cRight == nbins_mu-1) ) { cRight = 0;}
                        while( (fcount[indv][indmu - cLeft] == 0) && (indmu - cLeft > 0) ) { cLeft += 1; }
                        if( (fcount[indv][indmu - cLeft] == 0) && (indmu - cLeft == 0) ) { cLeft = 0;} 
                    } 
                    if( (cRight == 0) && (cLeft != 0) ) { 
                        dfdmu[indv][indmu]  = (fmu[indv][indmu + cRight] - fmu[indv][indmu-cLeft])/((cRight + cLeft)*dmubins) ;
                        dfdmu2[indv][indmu] = 0.0;
                    } else if( (cLeft == 0) && (cRight != 0) ) { 
                        dfdmu[indv][indmu]  = (fmu[indv][indmu + cRight] - fmu[indv][indmu-cLeft])/((cRight + cLeft)*dmubins) ;
                        dfdmu2[indv][indmu] = 0.0;
                    } else if( (cLeft == 0) && (cRight == 0) ) {
                        dfdmu[indv][indmu]  = 0.0;
                        dfdmu2[indv][indmu] = 0.0;
                    } else {
                        dfdmu[indv][indmu]  = (fmu[indv][indmu + cRight] - fmu[indv][indmu-cLeft])/((cRight + cLeft)*dmubins) ;
                        dfdmu2[indv][indmu] = ( (fmu[indv][indmu + cRight] - fmu[indv][indmu])/(cRight*dmubins) - (fmu[indv][indmu] - fmu[indv][indmu-cLeft])/(cLeft*dmubins) ) / (0.5 * dmubins * (cRight + cLeft)); 
                    }
                }
            } 

            if (subCount == 0) {
                // Save dfdmu to text
                std::string path_save = "/wrk/users/dubart/300_test/proc_test/mu_files/";
                std::ostringstream tmp;
                tmp << std::setw(7) << std::setfill('0') << P::tstep;
                std::string tstepString = tmp.str();
                std::ofstream dfdmu_array(path_save + "dfdmu_array_" + tstepString + ".txt");
                for (int indv = 0; indv < nbins_v; indv++) {
                    for(int indmu = 0; indmu < nbins_mu; indmu++) {
                        dfdmu_array << dfdmu[indv][indmu] << ' ';
                    }
                    dfdmu_array << std::endl;
                }
            }
            
            // Compute dfdt_mu
            for (int indv = 0; indv < nbins_v; indv++) { 
                for(int indmu = 0; indmu < nbins_mu; indmu++) {
                    dfdt_mu[indv][indmu] = Parameters::PADcoefficient * (
                                           - 2.0 * (dmubins * (indmu+0.5) - 1.0) * dfdmu[indv][indmu]
                                           + (1.0 - (dmubins * (indmu+0.5) - 1.0)*(dmubins * (indmu+0.5) - 1.0)) * dfdmu2[indv][indmu] );
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

                   int Vcount;
                   Realf VminApply = 2.0*DV;
                   if (normV < VminApply) { continue; }
                   else { Vcount = static_cast<int>(floor((normV-Vmin) / dVbins)); }    

                   int mucount = static_cast<int>(floor((mu+1.0) / dmubins));

                   Realf CellValue = cell.get_value(VX,VY,VZ,popID);
                   if (CellValue < Sparsity) {CellValue = Sparsity;} //Set CellValue to sparsity Threshold for empty cells otherwise div by 0
                   //Realf CellCalc  = 1.0; // f is always < 1.0
                   //if ((CellValue == 0.0) && (dfdmu[Vcount][mucount] != 0.0)) {
                   //    Realf CellValuePDX = cell.get_value(VX+DV,VY,VZ,popID); 
                   //    Realf CellValueMDX = cell.get_value(VX-DV,VY,VZ,popID); 
                   //    Realf CellValuePDY = cell.get_value(VX,VY+DV,VZ,popID); 
                   //    Realf CellValueMDY = cell.get_value(VX,VY-DV,VZ,popID);  
                   //    Realf CellValuePDZ = cell.get_value(VX,VY,VZ+DV,popID); 
                   //    Realf CellValueMDZ = cell.get_value(VX,VY,VZ-DV,popID);
                   //    std::array<Realf,6> Compare = {CellValuePDX,CellValueMDX,CellValuePDY,CellValueMDY,CellValuePDZ,CellValueMDZ};
                   //    for (int indx = 0; indx < Compare.size(); indx++) { 
                   //        if ((Compare[indx] < CellCalc) && (Compare[indx] != 0.0)) {CellCalc = Compare[indx];}
                   //        else { continue;}
                   //    }
                   //    if (CellCalc == 1.0) {continue;} // means all CellValues = 0.0
                   //} else if ((CellValue == 0.0) && (dfdmu[Vcount][mucount] == 0.0)) {continue;}
                   //else {CellCalc = CellValue;}

                   if (fmu[Vcount][mucount] == 0.0) { ratio[WID3*n+i+WID*j+WID*WID*k] = 1.0; }
                   else { ratio[WID3*n+i+WID*j+WID*WID*k] = CellValue / fmu[Vcount][mucount]; }                  

                   dfdt[WID3*n+i+WID*j+WID*WID*k] = dfdt_mu[Vcount][mucount] * ratio[WID3*n+i+WID*j+WID*WID*k];
                   
                   if (abs(dfdt[WID3*n+i+WID*j+WID*WID*k]) > 0.0) {
                   checkCFL[WID3*n+i+WID*j+WID*WID*k] = CellValue * Parameters::PADCFL * (1.0 / abs(dfdt[WID3*n+i+WID*j+WID*WID*k]));}
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

            if (subCount == 0) {
                // Save dfdt to text
                std::string path_save = "/wrk/users/dubart/300_test/proc_test/mu_files/";
                std::ostringstream tmp;
                tmp << std::setw(7) << std::setfill('0') << P::tstep;
                std::string tstepString = tmp.str();
                std::ofstream dfdt_array(path_save + "dfdt_array_" + tstepString + ".txt");
                std::ofstream checks(path_save + "checks_" + tstepString + ".txt");
                for (vmesh::LocalID n=0; n<cell.get_number_of_velocity_blocks(popID); n++) {
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

                        dfdt_array << VX << " " << VY << " " << VZ << " " << dfdt[WID3*n+i+WID*j+WID*WID*k] << std::endl;
                        checks << VX << " " << VY << " " << VZ << " " << ratio[WID3*n+i+WID*j+WID*WID*k] << " " << Ddt << std::endl;
                    }
                }
                std::ofstream DifftStep(path_save + "Ddt.txt", std::ios_base::app);
                DifftStep << P::tstep << " " << Ddt << std::endl;           
            }

            //Loop to update cell
            for (vmesh::LocalID n=0; n<cell.get_number_of_velocity_blocks(popID); n++) { //Iterate through velocity blocks
                for (uint k = 0; k < WID; ++k) for (uint j = 0; j < WID; ++j) for (uint i = 0; i < WID; ++i) {
                    const Real* parameters  = cell.get_block_parameters(popID);

                    //Get velocity space coordinates                   
                    const Real VX = parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VXCRD] + (i + 0.5)*parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVX];
                    const Real VY = parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VYCRD] + (j + 0.5)*parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVY];
                    const Real VZ = parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VZCRD] + (k + 0.5)*parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVZ];

                    Realf CellValue = cell.get_value(VX,VY,VZ,popID);

                    //Update cell
                    Realf NewCellValue = CellValue + dfdt[WID3*n+i+WID*j+WID*WID*k] * Ddt;
                    if (NewCellValue <= 0.0) { NewCellValue = 0.0;}

                    cell.set_value(VX,VY,VZ,NewCellValue,popID);
               }
           }

        subCount += 1;
        } // End Time loop

    } // End spatial cell loop

} // End function
