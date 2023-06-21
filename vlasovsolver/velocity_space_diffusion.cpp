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

#ifdef DPF
typedef Vec4d Vec4;
typedef Vec4db Vec4b; 
#define to_real(x) to_double(x)
#else
typedef Vec4f Vec4;
typedef Vec4fb Vec4b;
#define to_real(x) to_float(x)

// TODO: FIX if WID is not 4

using namespace spatial_cell;

static bool checkExistingNeighbour(SpatialCell* cell, Realf VX, Realf VY, Realf VZ, const uint popID) {

      const vmesh::GlobalID blockGID = cell->get_velocity_block(popID,VX, VY, VZ);
      vmesh::LocalID blockLID        = cell->get_population(popID).vmesh.getLocalID(blockGID);
      return blockLID               != vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>::invalidLocalID();
}

void velocitySpaceDiffusion(
        dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,const uint popID){


    int nbins_v  = Parameters::PADvbins;
    int nbins_mu = Parameters::PADmubins;
 
    Realf dmubins = 2.0/nbins_mu;

    int fcount   [nbins_v][nbins_mu]; // Array to count number of f stored
    Realf fmu    [nbins_v][nbins_mu]; // Array to store f(v,mu)
    Realf dfdmu  [nbins_v][nbins_mu]; // Array to store dfdmu
    Realf dfdmu2 [nbins_v][nbins_mu]; // Array to store dfdmumu
    Realf dfdt_mu[nbins_v][nbins_mu]; // Array to store dfdt_mu

    std::string PATHfile = Parameters::PADnu0;

    // Read from nu0Box.txt to get nu0 data depending on betaPara and Taniso
    std::ifstream FILEDmumu;
    FILEDmumu.open(PATHfile);

        // Get betaPara
    std::string lineBeta;
    for (int i = 0; i < 2; i++) { std::getline(FILEDmumu,lineBeta); }
         
    std::istringstream issBeta(lineBeta);
    std::vector<Realf> betaParaArray;

    float numBeta;
    while ((issBeta >> numBeta)) { betaParaArray.push_back(numBeta); }

        // Get Taniso
    std::string lineTaniso;
    for (int i = 0; i < 2; i++) { std::getline(FILEDmumu,lineTaniso); }
         
    std::istringstream issTaniso(lineTaniso);
    std::vector<Realf> TanisoArray;

    float numTaniso;
    while ((issTaniso >> numTaniso)) { TanisoArray.push_back(numTaniso); }
    
        // Discard one line 
    std::string lineDUMP;
    for (int i = 0; i < 1; i++) { std::getline(FILEDmumu,lineDUMP); }

        // Get nu0
    std::string linenu0;
    Realf nu0Array[betaParaArray.size()][TanisoArray.size()]; 
    
    for (int i = 0; i < betaParaArray.size(); i++) {
        std::getline(FILEDmumu,linenu0);
        std::istringstream issnu0(linenu0);
        std::vector<Realf> tempLINE;
        float numTEMP;
        while((issnu0 >> numTEMP)) { tempLINE.push_back(numTEMP); }
        for (int j = 0; j < tempLINE.size(); j++) {nu0Array[i][j] = tempLINE [j];}
    }



    const auto LocalCells=getLocalCells();
    #pragma omp parallel for private(fcount,fmu,dfdmu,dfdmu2,dfdt_mu)
    for (int CellIdx = 0; CellIdx < LocalCells.size(); CellIdx++) { //Iterate through spatial cell

        phiprof::start("Initialisation");
        auto CellID                        = LocalCells[CellIdx];
        SpatialCell& cell                  = *mpiGrid[CellID];
	const Real* parameters             = cell.get_block_parameters(popID);
        const vmesh::LocalID* nBlocks      = cell.get_velocity_grid_length(popID);
        const vmesh::MeshParameters& vMesh = getObjectWrapper().velocityMeshes[0];      

        Realf density_pre_adjust  = 0.0;
        Realf density_post_adjust = 0.0;

        for (size_t i=0; i<cell.get_number_of_velocity_blocks(popID)*WID3; ++i) {
            density_pre_adjust += cell.get_data(popID)[i];
        }

        Realf Sparsity    = 0.01 * cell.getVelocityBlockMinValue(popID);

        Realf dtTotalDiff = 0.0; // Diffusion time elapsed

        Realf Vmin   = 0.0; // In case we need to avoid center cells
        Realf Vmax   = 2*sqrt(3)*vMesh.meshLimits[1];
        Realf dVbins = (Vmax - Vmin)/nbins_v;  
    
        std::array<Realf,4> dfdt;

        std::array<Realf,3> bulkV = {cell.parameters[CellParams::VX], cell.parameters[CellParams::VY], cell.parameters[CellParams::VZ]};
        phiprof::stop("Initialisation");

        std::array<Realf,3> B = {cell.parameters[CellParams::PERBXVOL] +  cell.parameters[CellParams::BGBXVOL],
                                 cell.parameters[CellParams::PERBYVOL] +  cell.parameters[CellParams::BGBYVOL],
	                         cell.parameters[CellParams::PERBZVOL] +  cell.parameters[CellParams::BGBZVOL]};

        Realf Bnorm           = sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);
        std::array<Realf,3> b = {B[0]/Bnorm, B[1]/Bnorm, B[2]/Bnorm};
  

        // There is probably a simpler way to do all of the following but I was too lazy to try to understand how to do it in C++ so I wrote everything myself.       
        // Build Pressure tensor for calculation of anisotropy and beta
        Realf Ptensor[3][3] = {{cell.parameters[CellParams::P_11],cell.parameters[CellParams::P_12],cell.parameters[CellParams::P_13]},
                               {cell.parameters[CellParams::P_12],cell.parameters[CellParams::P_22],cell.parameters[CellParams::P_23]},
                               {cell.parameters[CellParams::P_13],cell.parameters[CellParams::P_23],cell.parameters[CellParams::P_33]}};

        // Build Rotation Matrix
        std::array<Realf,3> eZ       = {0.0,0.0,1.0}; // Rotation around z-axis
        std::array<Realf,3> BeZCross = { B[1]*eZ[2] - B[2]*eZ[1], B[2]*eZ[0] - B[0]*eZ[2], B[0]*eZ[1] - B[1]*eZ[0] };
        Realf vecNorm                = sqrt(BeZCross[0]*BeZCross[0] + BeZCross[1]*BeZCross[1] + BeZCross[2]*BeZCross[2]);
        std::array<Realf,3> uR       = {BeZCross[0]/vecNorm, BeZCross[1]/vecNorm, BeZCross[2]/vecNorm};
        Realf beZDot                 = b[0]*eZ[0] + b[1]*eZ[1] + b[2]*eZ[2];
        Realf angle                  = acos(beZDot);
        
        Realf RMatrix[3][3] = { {cos(angle)+uR[0]*uR[0]*(1.0-cos(angle))      , uR[0]*uR[1]*(1.0-cos(angle))-uR[2]*sin(angle), uR[0]*uR[2]*(1.0-cos(angle))+uR[1]*sin(angle)},
                                {uR[1]*uR[0]*(1.0-cos(angle))+uR[2]*sin(angle), cos(angle)+uR[1]*uR[1]*(1.0-cos(angle))      , uR[1]*uR[2]*(1.0-cos(angle))-uR[0]*sin(angle)},
                                {uR[2]*uR[0]*(1.0-cos(angle))-uR[1]*sin(angle), uR[2]*uR[1]*(1.0-cos(angle))+uR[0]*sin(angle), cos(angle)+uR[2]*uR[2]*(1.0-cos(angle))      } };

        // Rotate Tensor (T' = RTR^{-1})
        Realf RT[3][3] = { {RMatrix[0][0]*Ptensor[0][0] + RMatrix[0][1]*Ptensor[1][0] + RMatrix[0][2]*Ptensor[2][0], RMatrix[0][0]*Ptensor[0][1] + RMatrix[0][1]*Ptensor[1][1] + RMatrix[0][2]*Ptensor[2][1], RMatrix[0][0]*Ptensor[0][2] + RMatrix[0][1]*Ptensor[1][2] + RMatrix[0][2]*Ptensor[2][2]},
                           {RMatrix[1][0]*Ptensor[0][0] + RMatrix[1][1]*Ptensor[1][0] + RMatrix[1][2]*Ptensor[2][0], RMatrix[1][0]*Ptensor[0][1] + RMatrix[1][1]*Ptensor[1][1] + RMatrix[1][2]*Ptensor[2][1], RMatrix[1][0]*Ptensor[0][2] + RMatrix[1][1]*Ptensor[1][2] + RMatrix[1][2]*Ptensor[2][2]},
                           {RMatrix[2][0]*Ptensor[0][0] + RMatrix[2][1]*Ptensor[1][0] + RMatrix[2][2]*Ptensor[2][0], RMatrix[2][0]*Ptensor[0][1] + RMatrix[2][1]*Ptensor[1][1] + RMatrix[2][2]*Ptensor[2][1], RMatrix[2][0]*Ptensor[0][2] + RMatrix[2][1]*Ptensor[1][2] + RMatrix[2][2]*Ptensor[2][2]} };

        Realf Rtranspose[3][3] ={ {cos(angle)+uR[0]*uR[0]*(1.0-cos(angle))      , uR[0]*uR[1]*(1.0-cos(angle))+uR[2]*sin(angle), uR[0]*uR[2]*(1.0-cos(angle))-uR[1]*sin(angle)},
                                  {uR[1]*uR[0]*(1.0-cos(angle))-uR[2]*sin(angle), cos(angle)+uR[1]*uR[1]*(1.0-cos(angle))      , uR[1]*uR[2]*(1.0-cos(angle))+uR[0]*sin(angle)},
                                  {uR[2]*uR[0]*(1.0-cos(angle))+uR[1]*sin(angle), uR[2]*uR[1]*(1.0-cos(angle))-uR[0]*sin(angle), cos(angle)+uR[2]*uR[2]*(1.0-cos(angle))      } };

        Realf PtensorRotated[3][3] = {{RT[0][0]*Rtranspose[0][0] + RT[0][1]*Rtranspose[1][0] + RT[0][2]*Rtranspose[2][0], RT[0][0]*Rtranspose[0][1] + RT[0][1]*Rtranspose[1][1] + RT[0][2]*Rtranspose[2][1], RT[0][0]*Rtranspose[0][2] + RT[0][1]*Rtranspose[1][2] + RT[0][2]*Rtranspose[2][2]},
                                     {RT[1][0]*Rtranspose[0][0] + RT[1][1]*Rtranspose[1][0] + RT[1][2]*Rtranspose[2][0], RT[1][0]*Rtranspose[0][1] + RT[1][1]*Rtranspose[1][1] + RT[1][2]*Rtranspose[2][1], RT[1][0]*Rtranspose[0][2] + RT[1][1]*Rtranspose[1][2] + RT[1][2]*Rtranspose[2][2]},
                                     {RT[2][0]*Rtranspose[0][0] + RT[2][1]*Rtranspose[1][0] + RT[2][2]*Rtranspose[2][0], RT[2][0]*Rtranspose[0][1] + RT[2][1]*Rtranspose[1][1] + RT[2][2]*Rtranspose[2][1], RT[2][0]*Rtranspose[0][2] + RT[2][1]*Rtranspose[1][2] + RT[2][2]*Rtranspose[2][2]} };

        // Anisotropy
        Realf Taniso = 0.5 * (PtensorRotated[0][0] + PtensorRotated[1][1]) / PtensorRotated[2][2];
        // Beta Parallel
        Realf betaParallel = 2.0 * physicalconstants::MU_0 * PtensorRotated[2][2] / (Bnorm*Bnorm);

        // Find indx for first larger value
        int betaIndx   = -1;
        int TanisoIndx = -1;
        for (int i = 0; i < betaParaArray.size(); i++) { if (betaParallel >= betaParaArray[i]) { betaIndx   = i;} }
        for (int i = 0; i < TanisoArray.size()  ; i++) { if (Taniso       >= TanisoArray[i]  ) { TanisoIndx = i;} }

        Realf nu0     = 0.0;
        Realf epsilon = 0.0;

        if ( (betaIndx < 0) || (betaIndx >= betaParaArray.size()-1) || (TanisoIndx < 0) || (TanisoIndx >= TanisoArray.size()-1) ) {continue;}
        else { 
            // bi-linear interpolation with weighted mean to find nu0(betaParallel,Taniso)
            Realf beta1   = betaParaArray[betaIndx];
            Realf beta2   = betaParaArray[betaIndx+1];
            Realf Taniso1 = TanisoArray[TanisoIndx];
            Realf Taniso2 = TanisoArray[TanisoIndx+1];
            Realf nu011   = nu0Array[betaIndx][TanisoIndx];
            Realf nu012   = nu0Array[betaIndx][TanisoIndx+1];
            Realf nu021   = nu0Array[betaIndx+1][TanisoIndx];
            Realf nu022   = nu0Array[betaIndx+1][TanisoIndx+1];
                // Weights
            Realf w11 = (beta2 - betaParallel)*(Taniso2 - Taniso)  / ( (beta2 - beta1)*(Taniso2-Taniso1) ); 
            Realf w12 = (beta2 - betaParallel)*(Taniso  - Taniso1) / ( (beta2 - beta1)*(Taniso2-Taniso1) ); 
            Realf w21 = (betaParallel - beta1)*(Taniso2 - Taniso)  / ( (beta2 - beta1)*(Taniso2-Taniso1) ); 
            Realf w22 = (betaParallel - beta1)*(Taniso  - Taniso1) / ( (beta2 - beta1)*(Taniso2-Taniso1) ); 
                   
            nu0     = w11*nu011 + w12*nu012 + w21*nu021 + w22*nu022;

            if (CellIdx == 0) {
                std::ostringstream tmpText; 
                tmpText << P::tstep << " " << betaParallel << " " << Taniso << " " << beta1 << " " << beta2 << " " << Taniso1 << " " << Taniso2 << " " << nu011 << " " << nu012 << " " << nu021 << " " << nu022 << " " << nu0 << std::endl;
                std::string tmpString = tmpText.str();
                std::cerr << tmpString;
            }
            if (nu0 <= 0.0) {continue;}
        }      
 
        phiprof::start("Subloop");
        while (dtTotalDiff < Parameters::dt) { // Substep loop

            phiprof::start("Zeroing");
            Realf RemainT  = Parameters::dt - dtTotalDiff; //Remaining time before reaching simulation time step
            Realf checkCFL = std::numeric_limits<Realf>::max();

            // Initialised back to zero at each substep
            memset(fmu          , 0.0, sizeof(fmu));
            memset(fcount       , 0.0, sizeof(fcount));

            phiprof::stop("Zeroing");

            phiprof::start("fmu building");
            // Build 2d array of f(v,mu)
            for (vmesh::LocalID n=0; n<cell.get_number_of_velocity_blocks(popID); n++) { // Iterate through velocity blocks
                for (uint k = 0; k < WID; ++k) for (uint j = 0; j < WID; ++j) { // Iterate through coordinates (z,y)

                   //Get velocity space coordinates                    
	           const Vec4 VX(parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VXCRD] 
                                  + (0 + 0.5)*parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVX],
                                  parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VXCRD] 
                                  + (1 + 0.5)*parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVX],
                                  parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VXCRD] 
                                  + (2 + 0.5)*parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVX],
                                  parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VXCRD] 
                                  + (3 + 0.5)*parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVX]);

                   const Vec4 VY(parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VYCRD] 
                                  + (j + 0.5)*parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVY]);

                   const Vec4 VZ(parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VZCRD]
                                  + (k + 0.5)*parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVZ]);

                   std::array<Vec4,3> V = {VX,VY,VZ}; // Velocity in the cell, in the simulation frame
                                   
                   std::array<Vec4,3> Vplasma; // Velocity in the cell, in the plasma frame

                   for (int indx = 0; indx < 3; indx++) { Vplasma[indx] = (V[indx] - Vec4(bulkV[indx])); }
              
                   Vec4 normV = sqrt(Vplasma[0]*Vplasma[0] + Vplasma[1]*Vplasma[1] + Vplasma[2]*Vplasma[2]);

                   Vec4 Vpara = Vplasma[0]*b[0] + Vplasma[1]*b[1] + Vplasma[2]*b[2];

                   Vec4 mu = Vpara/(normV+std::numeric_limits<Realf>::min()); // + min value to avoid division by 0
 
                   Vec4i Vindex;
                   Vindex = round_to_int(floor((normV-Vmin) / dVbins));

                   Vec4i muindex;
                   muindex = round_to_int(floor((mu+1.0) / dmubins));
                   for (uint i = 0; i<WID; i++) {if (muindex[i] == nbins_mu) {muindex[i] == muindex[i] - 1;}} // Safety check to handle edge case where mu = exactly 1.0

                   Vec4 Vmu = dVbins * (to_real(Vindex)+0.5); // Take value at the center of the mu cell

                   Realf* CellValue = &cell.get_data(n,popID)[WID*j+WID*WID*k];

                   for (uint i = 0; i<WID; i++) {
                       fmu   [Vindex[i]][muindex[i]] += 2.0 * M_PI * Vmu[i]*Vmu[i] * CellValue[i];
                       fcount[Vindex[i]][muindex[i]] += 1;
                   }
                   
                } // End coordinates
            } // End blocks
            phiprof::stop("fmu building");

            phiprof::start("space/time derivatives, CFL, Ddt");
            int cRight;
            int cLeft;
            Realf checkCFLtmp = std::numeric_limits<Realf>::max();

            // Compute space/time derivatives (take first non-zero neighbours) & CFL & Ddt
            for (int indv = 0; indv < nbins_v; indv++) { 
    
                // Divide f by count (independent of v but needs to be computed for all mu before derivatives)
                for(int indmu = 0; indmu < nbins_mu; indmu++) { 
                    if (fcount[indv][indmu] == 0 || fmu[indv][indmu] <= 0.0) { fmu[indv][indmu] = std::numeric_limits<Realf>::min();}
                    else {fmu[indv][indmu] = fmu[indv][indmu] / fcount[indv][indmu];} 
                }

                for(int indmu = 0; indmu < nbins_mu; indmu++) {
                    // Compute spatial derivatives
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

                    // Compute time derivative
                    Realf mu    = (indmu+0.5)*dmubins - 1.0;
                    Realf Dmumu = nu0/2.0 * ( abs(mu)/(1.0 + abs(mu)) + epsilon ) * (1.0 - mu*mu);
                    Realf dDmu  = nu0/2.0 * ( (mu/abs(mu)) * ((1.0 - mu*mu)/((1.0 + abs(mu))*(1.0 + abs(mu)))) - 2.0*mu*( abs(mu)/(1.0 + abs(mu)) + epsilon));
                    dfdt_mu[indv][indmu] = dDmu * dfdmu[indv][indmu] + Dmumu * dfdmu2[indv][indmu];

                    // Compute CFL
                    Realf Vmu = dVbins * (float(indv)+0.5);
                    if (fmu[indv][indmu] > Sparsity*(2.0 * M_PI * Vmu*Vmu) && abs(dfdt_mu[indv][indmu]) > 0.0) { checkCFLtmp = fmu[indv][indmu] * Parameters::PADCFL * (1.0/abs(dfdt_mu[indv][indmu])); }
                    if (checkCFLtmp < checkCFL) { checkCFL = checkCFLtmp; }

                } // End mu loop
            } // End v loop

            // Compute Ddt
            Realf Ddt = checkCFL;
            if (Ddt > RemainT) { Ddt = RemainT; }
            dtTotalDiff = dtTotalDiff + Ddt;
            phiprof::stop("space/time derivatives, CFL, Ddt");

            phiprof::start("diffusion time derivative & update cell");
            // Compute dfdt
            for (vmesh::LocalID n=0; n<cell.get_number_of_velocity_blocks(popID); n++) { // Iterate through velocity blocks             
                for (uint k = 0; k < WID; ++k) for (uint j = 0; j < WID; ++j) {
                   //Get velocity space coordinates                    
	           const Vec4 VX(parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VXCRD] 
                                  + (0 + 0.5)*parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVX],
                                  parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VXCRD] 
                                  + (1 + 0.5)*parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVX],
                                  parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VXCRD] 
                                  + (2 + 0.5)*parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVX],
                                  parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VXCRD] 
                                  + (3 + 0.5)*parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVX]);

                   const Vec4 VY(parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VYCRD] 
                                  + (j + 0.5)*parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVY]);

                   const Vec4 VZ(parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VZCRD]
                                  + (k + 0.5)*parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVZ]);
                   
                   std::array<Vec4,3> V = {VX,VY,VZ}; // Velocity in the cell, in the simulation frame
                                   
                   std::array<Vec4,3> Vplasma; // Velocity in the cell, in the plasma frame

                   for (int indx = 0; indx < 3; indx++) { Vplasma[indx] = (V[indx] - Vec4(bulkV[indx])); }
              
                   Vec4 normV = sqrt(Vplasma[0]*Vplasma[0] + Vplasma[1]*Vplasma[1] + Vplasma[2]*Vplasma[2]);

                   Vec4 Vpara = Vplasma[0]*b[0] + Vplasma[1]*b[1] + Vplasma[2]*b[2];

                   Vec4 mu = Vpara/(normV+std::numeric_limits<Realf>::min()); // + min value to avoid division by 0

                   Vec4 CellValue;
                   Vec4 NewCellValue;
                   #endif
                   CellValue.load(&cell.get_data(n,popID)[WID*j+WID*WID*k]);

                   Vec4i Vindex;
                   Vindex = round_to_int(floor((normV-Vmin) / dVbins));
                   
                   Vec4i muindex;
                   muindex = round_to_int(floor((mu+1.0) / dmubins));
                   for (uint i = 0; i<WID; i++) {if (muindex[i] == nbins_mu) {muindex[i] == muindex[i] - 1;}} // Safety check to handle edge case where mu = exactly 1.0

                   Vec4 Vmu = dVbins * (to_real(Vindex)+0.5);

                   for (uint i = 0; i < WID; i++) {
                       dfdt[i] = dfdt_mu[Vindex[i]][muindex[i]] / (2.0 * M_PI * Vmu[i]*Vmu[i]);
                   }

                   //Update cell
                   Vec4 dfdtUpdate;
                   dfdtUpdate.load(&dfdt[0]);
                   NewCellValue    = CellValue + dfdtUpdate * Ddt;
                   Vec4b lessZero = NewCellValue < 0.0;
                   NewCellValue    = select(lessZero,0.0,NewCellValue);
                   NewCellValue.store(&cell.get_data(n,popID)[WID*j+WID*WID*k]);

                } // End coordinates 
            } // End Blocks
            phiprof::stop("diffusion time derivative & update cell");

        } // End Time loop
        phiprof::stop("Subloop");

        for (size_t i=0; i<cell.get_number_of_velocity_blocks(popID)*WID3; ++i) {
            density_post_adjust += cell.get_data(popID)[i];
        }

        if (density_post_adjust != 0.0) {
            for (size_t i=0; i<cell.get_number_of_velocity_blocks(popID)*WID3; ++i) {
                cell.get_data(popID)[i] *= density_pre_adjust/density_post_adjust;
            }
        }

    } // End spatial cell loop

} // End function
