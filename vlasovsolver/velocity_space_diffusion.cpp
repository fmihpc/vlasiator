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

        phiprof::start("Initialisation");
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

        std::vector<Realf> dfdt       (cell.get_number_of_velocity_blocks(popID) * WID3); // Array of vspace size to store dfdt
        std::vector<int> Vcount_array (cell.get_number_of_velocity_blocks(popID) * WID3); // Array to store vcount per cell
        std::vector<int> mucount_array(cell.get_number_of_velocity_blocks(popID) * WID3); // Array to store mucount per cell
        int fcount   [nbins_v][nbins_mu];                                                 // Array to count number of f stored
        Realf fmu    [nbins_v][nbins_mu];                                                 // Array to store f(v,mu)
        Realf dfdmu  [nbins_v][nbins_mu];                                                 // Array to store dfdmu
        Realf dfdmu2 [nbins_v][nbins_mu];                                                 // Array to store dfdmumu
        Realf dfdt_mu[nbins_v][nbins_mu];                                                 // Array to store dfdt_mu
        
        std::vector<Realf> checkCFL (cell.get_number_of_velocity_blocks(popID)*WID3, std::numeric_limits<Realf>::max()); // Array of vspace size to store checkCFL (This one hates me)
        
        std::array<Realf,3> bulkV = {cell.parameters[CellParams::VX], cell.parameters[CellParams::VY], cell.parameters[CellParams::VZ]};
        phiprof::stop("Initialisation");

        while (dtTotalDiff < Parameters::dt) { // Substep loop

            phiprof::start("Zeroing");
            Realf RemainT = Parameters::dt - dtTotalDiff; //Remaining time before reaching simulation time step
 
            // Initialised back to zero at each substep
            memset(fmu          , 0.0, sizeof(fmu));
            memset(fcount       , 0.0, sizeof(fcount));

            phiprof::stop("Zeroing");

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

                   for (int indx = 0; indx < 3; indx++) { Vplasma[indx] = (V[indx] - Vec4d(bulkV[indx])); }
              
                   Vec4d normV = sqrt(Vplasma[0]*Vplasma[0] + Vplasma[1]*Vplasma[1] + Vplasma[2]*Vplasma[2]);

                   Vec4d Vpara = Vplasma[0];

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

                   Vcount .store(&Vcount_array [WID3*n+WID*j+WID*WID*k]);
                   mucount.store(&mucount_array[WID3*n+WID*j+WID*WID*k]);

                   for (uint i = 0; i<WID; i++) {
                       fmu   [Vcount[i]][mucount[i]] += 2.0 * M_PI * Vmu[i]*Vmu[i] * CellValue[i];
                       fcount[Vcount[i]][mucount[i]] += 1;
                   }
                   
                } // End coordinates
            } // End blocks
            phiprof::stop("fmu building");

            for (int indv = 0; indv < nbins_v; indv++) { // Divide f by count 
                for(int indmu = 0; indmu < nbins_mu; indmu++) {
                    if (fcount[indv][indmu] == 0) { fmu[indv][indmu] = 0.0;}
                    else {fmu[indv][indmu] = fmu[indv][indmu] / fcount[indv][indmu];} 
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
                        while( (fcount[indv][indmu + cRight] == 0) && (indmu + cRight < nbins_mu-1) )  { cRight += 1; }
                        if(    (fcount[indv][indmu + cRight] == 0) && (indmu + cRight == nbins_mu-1) ) { cRight  = 0; }
                    } else if (indmu == nbins_mu-1) {
                        cLeft  = 1;
                        cRight = 0;
                        while( (fcount[indv][indmu - cLeft] == 0) && (indmu - cLeft > 0) )  { cLeft += 1; }
                        if(    (fcount[indv][indmu - cLeft] == 0) && (indmu - cLeft == 0) ) { cLeft  = 0; }
                    } else {
                        cLeft  = 1;
                        cRight = 1;
                        while( (fcount[indv][indmu + cRight] == 0) && (indmu + cRight < nbins_mu-1) )  { cRight += 1; }
                        if(    (fcount[indv][indmu + cRight] == 0) && (indmu + cRight == nbins_mu-1) ) { cRight  = 0; }
                        while( (fcount[indv][indmu - cLeft ] == 0) && (indmu - cLeft  > 0) )           { cLeft  += 1; }
                        if(    (fcount[indv][indmu - cLeft ] == 0) && (indmu - cLeft  == 0) )          { cLeft   = 0; } 
                    } 
                    if( (cRight == 0) && (cLeft != 0) ) { 
                        dfdmu [indv][indmu] = (fmu[indv][indmu + cRight] - fmu[indv][indmu-cLeft])/((cRight + cLeft)*dmubins) ;
                        dfdmu2[indv][indmu] = 0.0;
                    } else if( (cLeft == 0) && (cRight != 0) ) { 
                        dfdmu [indv][indmu] = (fmu[indv][indmu + cRight] - fmu[indv][indmu-cLeft])/((cRight + cLeft)*dmubins) ;
                        dfdmu2[indv][indmu] = 0.0;
                    } else if( (cLeft == 0) && (cRight == 0) ) {
                        dfdmu [indv][indmu] = 0.0;
                        dfdmu2[indv][indmu] = 0.0;
                    } else {
                        dfdmu [indv][indmu] = (  fmu[indv][indmu + cRight] - fmu[indv][indmu-cLeft])/((cRight + cLeft)*dmubins) ;
                        dfdmu2[indv][indmu] = ( (fmu[indv][indmu + cRight] - fmu[indv][indmu])/(cRight*dmubins) - (fmu[indv][indmu] - fmu[indv][indmu-cLeft])/(cLeft*dmubins) ) / (0.5 * dmubins * (cRight + cLeft)); 
                    }
                } // End mu loop
            } // End v loop
            phiprof::stop("spatial derivatives");

            phiprof::start("mu time derivatives");
            // Compute dfdt_mu
            for (int indv = 0; indv < nbins_v; indv++) { 
                for(int indmu = 0; indmu < nbins_mu; indmu++) {
                    dfdt_mu[indv][indmu] = Parameters::PADcoefficient * (
                                           - 2.0 * (dmubins * (indmu+0.5) - 1.0) * dfdmu[indv][indmu]
                                           + (1.0 - (dmubins * (indmu+0.5) - 1.0)*(dmubins * (indmu+0.5) - 1.0)) * dfdmu2[indv][indmu] );
                }
            } 
            phiprof::stop("mu time derivatives");

            phiprof::start("diffusion time derivative");
            // Compute dfdt
            for (vmesh::LocalID n=0; n<cell.get_number_of_velocity_blocks(popID); n++) { // Iterate through velocity blocks             
                for (uint k = 0; k < WID; ++k) for (uint j = 0; j < WID; ++j) {
                
                   #ifdef DPF
                   Vec4d CellValue;
                   #else 
                   Vec4f CellValue;
                   #endif
                   CellValue.load(&cell.get_data(n,popID)[WID*j+WID*WID*k]);

                   Vec4db lessSpars = CellValue < Sparsity;
                   CellValue = select(lessSpars, Sparsity, CellValue);

                   Vec4i Vcount;
                   Vec4i mucount;

                   Vcount .load(&Vcount_array [WID3*n+WID*j+WID*WID*k]);
                   mucount.load(&mucount_array[WID3*n+WID*j+WID*WID*k]);                   

                   Vec4d Vmu = dVbins * (to_double(Vcount)+0.5);

                   for (uint i = 0; i < WID; i++) {dfdt[WID3*n+i+WID*j+WID*WID*k] = dfdt_mu[Vcount[i]][mucount[i]] / (2.0 * M_PI * Vmu[i]*Vmu[i]);}
                   
                   Vec4d dfdtCheck;
                   dfdtCheck.load(&dfdt[WID3*n+WID*j+WID*WID*k]);

                   Vec4d checkCFLTemp;                   

                   Vec4db dfdtABS = abs(dfdtCheck) > 0.0;
                   
                   checkCFLTemp = select(dfdtABS, CellValue * Parameters::PADCFL * (1.0 / abs(dfdtCheck)), std::numeric_limits<Realf>::max());
                   checkCFLTemp.store(&checkCFL[WID3*n+WID*j+WID*WID*k]);
                   
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
                for (uint k = 0; k < WID; ++k) for (uint j = 0; j < WID; ++j) { 
                    const Real* parameters  = cell.get_block_parameters(popID);
 
                    #ifdef DPF
                    Vec4d CellValue;
                    Vec4d NewCellValue;
                    #else
                    Vec4f CellValue;
                    Vec4f NewCellValue;
                    #endif
                    CellValue.load(&cell.get_data(n,popID)[WID*j+WID*WID*k]);

                    Vec4d dfdtUpdate;
                    dfdtUpdate.load(&dfdt[WID3*n+WID*j+WID*WID*k]);

                    //Update cell
                    NewCellValue    = CellValue + dfdtUpdate * Ddt;
                    Vec4db lessZero = NewCellValue < 0.0;
                    NewCellValue    = select(lessZero,0.0,NewCellValue);
                    NewCellValue.store(&cell.get_data(n,popID)[WID*j+WID*WID*k]);
               } // End coordinates
           } // End block
           phiprof::stop("update cell");
        
        } // End Time loop

    } // End spatial cell loop

} // End function
