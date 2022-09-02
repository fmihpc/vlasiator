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
#include "vectorclass.h"

// TODO: FIX if WID is not 4

using namespace spatial_cell;

static bool checkExistingNeighbour(SpatialCell* cell, Realf VX, Realf VY, Realf VZ, const uint popID) {

      const vmesh::GlobalID blockGID = cell->get_velocity_block(popID,VX, VY, VZ);
      vmesh::LocalID blockLID        = cell->get_population(popID).vmesh.getLocalID(blockGID);
      return blockLID               != vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>::invalidLocalID();
}

void velocitySpaceDiffusion(
        dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,const uint popID){

    const auto LocalCells=getLocalCells();
    #pragma omp parallel for
    for (int CellIdx = 0; CellIdx < LocalCells.size(); CellIdx++) { //Iterate through spatial cell

        auto CellID                        = LocalCells[CellIdx];
        SpatialCell& cell                  = *mpiGrid[CellID];
	const Real* parameters             = cell.get_block_parameters(popID);
        const vmesh::LocalID* nBlocks      = cell.get_velocity_grid_length(popID);
        const vmesh::MeshParameters& vMesh = getObjectWrapper().velocityMeshes[0];      

        Realf Sparsity    = 0.01 * cell.getVelocityBlockMinValue(popID);

        Realf dtTotalDiff = 0.0; // Diffusion time elapsed

        int nbins_v  = Parameters::PADvbins;
        int nbins_mu = Parameters::PADmubins;

        Realf mumin   = -1.0;
        Realf mumax   = +1.0;
        Realf dmubins = (mumax - mumin)/nbins_mu;

        Realf Vmin   = 0.0; // In case we need to avoid center cells
        Realf Vmax   = 2*sqrt(3)*vMesh.meshLimits[1];
        Realf dVbins = (Vmax - Vmin)/nbins_v;  
        
        std::vector<Realf> dfdt       (cell.get_number_of_velocity_blocks(popID)*WID3); // Array of vspace size to store dfdt
        std::vector<Realf> checkCFL   (cell.get_number_of_velocity_blocks(popID)*WID3); // Array of vspace size to store checkCFl
        std::vector<int> Vcount_array (cell.get_number_of_velocity_blocks(popID)*WID3); // Array to store vcount per cell
        std::vector<int> mucount_array(cell.get_number_of_velocity_blocks(popID)*WID3); // Array to store mucount per cell
        std::vector<std::vector<int>>   fcount (nbins_v,std::vector<int>(nbins_mu));    // Array to count number of f stored
        std::vector<std::vector<Realf>> fmu    (nbins_v,std::vector<Realf>(nbins_mu));  // Array to store f(v,mu)
        std::vector<std::vector<Realf>> dfdmu  (nbins_v,std::vector<Realf>(nbins_mu));  // Array to store dfdmu
        std::vector<std::vector<Realf>> dfdmu2 (nbins_v,std::vector<Realf>(nbins_mu));  // Array to store dfdmumu
        std::vector<std::vector<Realf>> dfdt_mu(nbins_v,std::vector<Realf>(nbins_mu));  // Array to store dfdt_mu

        std::array<Realf,3> bulkV = {cell.parameters[CellParams::VX], cell.parameters[CellParams::VY], cell.parameters[CellParams::VZ]};

        while (dtTotalDiff < Parameters::dt) { // Substep loop

            Realf RemainT = Parameters::dt - dtTotalDiff; //Remaining time before reaching simulation time step

            dfdt.assign(cell.get_number_of_velocity_blocks(popID)*WID3,0.0);                                    // Array of vspace size to store dfdt
            checkCFL.assign(cell.get_number_of_velocity_blocks(popID)*WID3, std::numeric_limits<Realf>::max()); // Array of vspace size to store checkCFl initialized with max value to not mess up checkCFL for empty cells
            Vcount_array.assign (cell.get_number_of_velocity_blocks(popID)*WID3,0); // Array to store vcount per cell
            mucount_array.assign(cell.get_number_of_velocity_blocks(popID)*WID3,0); // Array to store mucount per cell

            // Initialized back to 0 at each substep
            for (int i=0; i<nbins_v; i++) {
                for (int j=0; j<nbins_mu; j++) {
                    fcount .at(i).at(j) = 0;
                    fmu    .at(i).at(j) = 0.0;
                    dfdmu  .at(i).at(j) = 0.0;
                    dfdmu2 .at(i).at(j) = 0.0;
                    dfdt_mu.at(i).at(j) = 0.0;
                }
            }

            phiprof::start("fmu building");
            // Build 2d array of f(v,mu)
            for (vmesh::LocalID n=0; n<cell.get_number_of_velocity_blocks(popID); n++) { // Iterate through velocity blocks
                for (uint k = 0; k < WID; ++k) for (uint j = 0; j < WID; ++j) { // Iterate through coordinates (z,y)

                   //Get velocity space coordinates                    
	           const Vec4d VX(parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VXCRD] 
                                  + (0 + 0.5)*parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVX],
                                  parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VXCRD] 
                                  + (1 + 0.5)*parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVX],
                                  parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VXCRD] 
                                  + (2 + 0.5)*parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVX],
                                  parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VXCRD] 
                                  + (3 + 0.5)*parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVX]);

                   const Vec4d VY(parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VYCRD] 
                                  + (j + 0.5)*parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVY]);

                   const Vec4d VZ(parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VZCRD]
                                  + (k + 0.5)*parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVZ]);
                   
                   std::array<Vec4d,3> V = {VX,VY,VZ}; // Velocity in the cell, in the simulation frame
                                   
                   std::array<Vec4d,3> Vplasma; // Velocity in the cell, in the plasma frame

                   for (int indx = 0; indx < 3; indx++) { Vplasma.at(indx) = (V.at(indx) - Vec4d(bulkV.at(indx))); }
              
                   Vec4d normV = sqrt(Vplasma.at(0)*Vplasma.at(0) + Vplasma.at(1)*Vplasma.at(1) + Vplasma.at(2)*Vplasma.at(2));

                   Vec4d Vpara = Vplasma.at(0);

                   Vec4d mu = Vpara/(normV+std::numeric_limits<Realf>::min()); // + min value to avoid division by 0
 
                   Vec4i Vcount;
                   Vcount = round_to_int(floor((normV-Vmin) / dVbins));

                   // Dont know how to handble that with vectors
                   // if (normV < Vmin) { continue; } // To avoid center cells if needed
                   // else { Vcount = round_to_int(floor((normV-Vmin) / dVbins)); }                      

                   Vec4i mucount;
                   mucount = round_to_int(floor((mu+1.0) / dmubins));                      

                   Vec4d Vmu = dVbins * (to_double(Vcount)+0.5); // Take value at the center of the mu cell

                   #ifdef DPF
                   Vec4d CellValue;
                   #else
                   Vec4f CellValue;
                   #endif
                   CellValue.load(&cell.get_data(n,popID)[WID*j+WID*WID*k]);

                   Vcount .store(&Vcount_array .at(WID*j+WID*WID*k));
                   mucount.store(&mucount_array.at(WID*j+WID*WID*k));

                   for (uint i = 0; i<WID; i++) {
                       fmu   .at(Vcount[i]).at(mucount[i]) += 2.0 * M_PI * Vmu[i]*Vmu[i] * CellValue[i];
                       fcount.at(Vcount[i]).at(mucount[i]) += 1;
                   }
                   
                } // End coordinates
            } // End blocks
            phiprof::stop("fmu building");

            for (int indv = 0; indv < nbins_v; indv++) { // Divide f by count 
                for(int indmu = 0; indmu < nbins_mu; indmu++) {
                    if (fcount.at(indv).at(indmu) == 0) { fmu.at(indv).at(indmu) = 0.0;}
                    else {fmu.at(indv).at(indmu) = fmu.at(indv).at(indmu) / fcount.at(indv).at(indmu);} 
                }
            }
            
            int cRight;
            int cLeft;

            phiprof::start("spatial derivatives");
            // Compute dfdmu and dfdmu2 (take first non-zero neighbours)
            for (int indv = 0; indv < nbins_v; indv++) { 
                for(int indmu = 0; indmu < nbins_mu; indmu++) {
                    if (indmu == 0) {
                        cLeft  = 0;
                        cRight = 1;
                        while( (fcount.at(indv).at(indmu + cRight) == 0) && (indmu + cRight < nbins_mu-1) ) { cRight += 1; }
                        if( (fcount.at(indv).at(indmu + cRight) == 0) && (indmu + cRight == nbins_mu-1) ) { cRight = 0;}
                    } else if (indmu == nbins_mu-1) {
                        cLeft  = 1;
                        cRight = 0;
                        while( (fcount.at(indv).at(indmu - cLeft) == 0) && (indmu - cLeft > 0) ) { cLeft += 1; }
                        if( (fcount.at(indv).at(indmu - cLeft) == 0) && (indmu - cLeft == 0) ) { cLeft = 0;}
                    } else {
                        cLeft  = 1;
                        cRight = 1;
                        while( (fcount.at(indv).at(indmu + cRight) == 0) && (indmu + cRight < nbins_mu-1) ) { cRight += 1; }
                        if( (fcount.at(indv).at(indmu + cRight) == 0) && (indmu + cRight == nbins_mu-1) ) { cRight = 0;}
                        while( (fcount.at(indv).at(indmu - cLeft) == 0) && (indmu - cLeft > 0) ) { cLeft += 1; }
                        if( (fcount.at(indv).at(indmu - cLeft) == 0) && (indmu - cLeft == 0) ) { cLeft = 0;} 
                    } 
                    if( (cRight == 0) && (cLeft != 0) ) { 
                        dfdmu .at(indv).at(indmu) = (fmu.at(indv).at(indmu + cRight) - fmu.at(indv).at(indmu-cLeft))/((cRight + cLeft)*dmubins) ;
                        dfdmu2.at(indv).at(indmu) = 0.0;
                    } else if( (cLeft == 0) && (cRight != 0) ) { 
                        dfdmu .at(indv).at(indmu) = (fmu.at(indv).at(indmu + cRight) - fmu.at(indv).at(indmu-cLeft))/((cRight + cLeft)*dmubins) ;
                        dfdmu2.at(indv).at(indmu) = 0.0;
                    } else if( (cLeft == 0) && (cRight == 0) ) {
                        dfdmu .at(indv).at(indmu) = 0.0;
                        dfdmu2.at(indv).at(indmu) = 0.0;
                    } else {
                        dfdmu .at(indv).at(indmu) = (fmu.at(indv).at(indmu + cRight) - fmu.at(indv).at(indmu-cLeft))/((cRight + cLeft)*dmubins) ;
                        dfdmu2.at(indv).at(indmu) = ( (fmu.at(indv).at(indmu + cRight) - fmu.at(indv).at(indmu))/(cRight*dmubins) - (fmu.at(indv).at(indmu) - fmu.at(indv).at(indmu-cLeft))/(cLeft*dmubins) ) / (0.5 * dmubins * (cRight + cLeft)); 
                    }
                } // End mu loop
            } // End v loop
            phiprof::stop("spatial derivatives");

            phiprof::start("mu time derivatives");
            // Compute dfdt_mu
            for (int indv = 0; indv < nbins_v; indv++) { 
                for(int indmu = 0; indmu < nbins_mu; indmu++) {
                    dfdt_mu.at(indv).at(indmu) = Parameters::PADcoefficient * (
                                           - 2.0 * (dmubins * (indmu+0.5) - 1.0) * dfdmu.at(indv).at(indmu)
                                           + (1.0 - (dmubins * (indmu+0.5) - 1.0)*(dmubins * (indmu+0.5) - 1.0)) * dfdmu2.at(indv).at(indmu) );
                }
            } 
            phiprof::stop("mu time derivatives");

            phiprof::start("diffusion time derivative");
            // Compute dfdt
            for (vmesh::LocalID n=0; n<cell.get_number_of_velocity_blocks(popID); n++) { // Iterate through velocity blocks             
                for (uint k = 0; k < WID; ++k) for (uint j = 0; j < WID; ++j) for (uint i = 0; i < WID; ++i) {
                
                   Realf CellValue = cell.get_data(n,popID)[i+WID*j+WID*WID*k];

                   int Vcount  = Vcount_array .at(WID3*n+i+WID*j+WID*WID*k);
                   int mucount = mucount_array.at(WID3*n+i+WID*j+WID*WID*k);

                   Realf Vmu = dVbins * (Vcount+0.5);

                   dfdt.at(WID3*n+i+WID*j+WID*WID*k) = dfdt_mu.at(Vcount).at(mucount) / (2.0 * M_PI * Vmu*Vmu);
                   
                   if (CellValue < Sparsity) {CellValue = Sparsity;} // Set CellValue to sparsity Threshold for empty cells otherwise div by 0
                   if (abs(dfdt.at(WID3*n+i+WID*j+WID*WID*k)) > 0.0) { checkCFL.at(WID3*n+i+WID*j+WID*WID*k) = CellValue * Parameters::PADCFL * (1.0 / abs(dfdt.at(WID3*n+i+WID*j+WID*WID*k))); }
                } // End coordinates 
            } // End Blocks
            phiprof::stop("diffusion time derivative");

            phiprof::start("calculate CFL");
            //Calculate Diffusion time step based on min of CFL condition
            std::vector<Realf>::iterator mincheckCFL;
            mincheckCFL = std::min_element(checkCFL.begin(),checkCFL.end());
            if (mincheckCFL == checkCFL.end()) {break;}
            Realf Ddt = *mincheckCFL; // Diffusion time step
            if (Ddt > RemainT) { Ddt = RemainT; }
            dtTotalDiff = dtTotalDiff + Ddt;
            phiprof::stop("calculate CFL");

            phiprof::start("update cell");
            //Loop to update cell
            for (vmesh::LocalID n=0; n<cell.get_number_of_velocity_blocks(popID); n++) { //Iterate through velocity blocks
                for (uint k = 0; k < WID; ++k) for (uint j = 0; j < WID; ++j) for (uint i = 0; i < WID; ++i) {
                    const Real* parameters  = cell.get_block_parameters(popID);

                    Realf CellValue = cell.get_data(n,popID)[i+WID*j+WID*WID*k];
                    
                    //Update cell
                    Realf NewCellValue;
                    NewCellValue = CellValue + dfdt.at(WID3*n+i+WID*j+WID*WID*k) * Ddt;
                    if (NewCellValue <= 0.0) { NewCellValue = 0.0;}

                    cell.get_data(n,popID)[i+WID*j+WID*WID*k] = NewCellValue;
               } // End coordinates
           } // End block
           phiprof::stop("update cell");
        
        } // End Time loop

    } // End spatial cell loop

} // End function
