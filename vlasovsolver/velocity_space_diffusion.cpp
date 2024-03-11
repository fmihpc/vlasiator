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
#include "vec.h"

using namespace spatial_cell;

template <typename Lambda> inline static void loop_over_block(Lambda loop_body) {

    for (uint k = 0; k < WID; ++k) {
        for (uint j = 0; j < WID; j +=VECL/WID) { // Iterate through coordinates (z,y)
        
	    // create vectors with the i and j indices in the vector position on the plane.
	    #if VECL == 4 && WID == 4
	    const Veci i_indices = Veci({0, 1, 2, 3});
	    const Veci j_indices = Veci({j, j, j, j});
	    #elif VECL == 4 && WID == 8
	    #error "__FILE__ : __LINE__ : VECL == 4 && WID == 8 cannot work!" 
	    #elif VECL == 8 && WID == 4
	    const Veci i_indices = Veci({0, 1, 2, 3,
		     0, 1, 2, 3});
	    const Veci j_indices = Veci({j, j, j, j,
		     j + 1, j + 1, j + 1, j + 1});
	    #elif VECL == 8 && WID == 8
	    const Veci i_indices = Veci({0, 1, 2, 3, 4, 5, 6, 7});
	    const Veci j_indices = Veci({j, j, j, j, j, j, j, j});
	    #elif VECL == 16 && WID == 4
	    const Veci i_indices = Veci({0, 1, 2, 3,
		     0, 1, 2, 3,
		     0, 1, 2, 3,
		     0, 1, 2, 3});
	    const Veci j_indices = Veci({j, j, j, j,
		     j + 1, j + 1, j + 1, j + 1,
		     j + 2, j + 2, j + 2, j + 2,
		     j + 3, j + 3, j + 3, j + 3});
	    #elif VECL == 16 && WID == 8
	    const Veci i_indices = Veci({0, 1, 2, 3, 4, 5, 6, 7,
		     0, 1, 2, 3, 4, 5, 6, 7});
	    const Veci j_indices = Veci({j,   j,   j,   j,   j,   j,   j,   j,
		     j+1, j+1, j+1, j+1, j+1, j+1, j+1, j+1});
	    #elif VECL == 16 && WID == 16
	    const Veci i_indices = Veci({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15});
	    const Veci j_indices = Veci({j, j, j, j, j, j, j, j, j, j,  j,  j,  j,  j,  j, j});
	    #elif VECL == 32 && WID == 4
	    #error "__FILE__ : __LINE__ : VECL == 4 && WID == 8 cannot work, too long vector for one plane!" 
	    #elif VECL == 32 && WID == 8
	    const Veci i_indices = Veci({0, 1, 2, 3, 4, 5, 6, 7,
		     0, 1, 2, 3, 4, 5, 6, 7,
		     0, 1, 2, 3, 4, 5, 6, 7,
		     0, 1, 2, 3, 4, 5, 6, 7});
	    const Veci j_indices = Veci({j,   j,   j,   j,   j,   j,   j,   j,
		     j+1, j+1, j+1, j+1, j+1, j+1, j+1, j+1,
		     j+2, j+2, j+2, j+2, j+2, j+2, j+2, j+2,
		     j+3, j+3, j+3, j+3, j+3, j+3, j+3, j+3});
	    #elif VECL == 64 && WID == 4
	    #error "__FILE__ : __LINE__ : VECL == 4 && WID == 8 cannot work, too long vector for one plane!" 
	    #elif VECL == 64 && WID == 8
	    const Veci i_indices = Veci({0, 1, 2, 3, 4, 5, 6, 7,
		     0, 1, 2, 3, 4, 5, 6, 7,
		     0, 1, 2, 3, 4, 5, 6, 7,
		     0, 1, 2, 3, 4, 5, 6, 7,
		     0, 1, 2, 3, 4, 5, 6, 7,
		     0, 1, 2, 3, 4, 5, 6, 7,
		     0, 1, 2, 3, 4, 5, 6, 7,
		     0, 1, 2, 3, 4, 5, 6, 7});
	    const Veci j_indices = Veci({j,   j,   j,   j,   j,   j,   j,   j,
		     j+1, j+1, j+1, j+1, j+1, j+1, j+1, j+1,
		     j+2, j+2, j+2, j+2, j+2, j+2, j+2, j+2,
		     j+3, j+3, j+3, j+3, j+3, j+3, j+3, j+3,
		     j+4, j+4, j+4, j+4, j+4, j+4, j+4, j+4,
		     j+5, j+5, j+5, j+5, j+5, j+5, j+5, j+5,
		     j+6, j+6, j+6, j+6, j+6, j+6, j+6, j+6,
		     j+7, j+7, j+7, j+7, j+7, j+7, j+7, j+7});
	    #else
	    #error "This VECL && This WID cannot work!" 
            #define xstr(s) str(s)
            #define str(s) #s
            #pragma message "VECL =" xstr(VECL)
            #pragma message "WID = "xstr(WID) 
	    #endif
            
            loop_body(i_indices,j_indices,k);

        }
    }
}

static bool checkExistingNeighbour(SpatialCell* cell, Real VX, Real VY, Real VZ, const uint popID) {

      const vmesh::GlobalID blockGID = cell->get_velocity_block(popID,VX, VY, VZ);
      vmesh::LocalID blockLID        = cell->get_population(popID).vmesh.getLocalID(blockGID);
      return blockLID               != vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>::invalidLocalID();
}

void velocitySpaceDiffusion(
        dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,const uint popID){


    int nbins_v  = Parameters::PADvbins;
    int nbins_mu = Parameters::PADmubins;
 
    Real dmubins = 2.0/nbins_mu;

    int  fcount [nbins_v*nbins_mu]; // Array to count number of f stored
    int*  fcount_p = reinterpret_cast<int*> (fcount);
    Real fmu    [nbins_v*nbins_mu]; // Array to store f(v,mu)
    Real* fmu_p    = reinterpret_cast<Real*> (fmu);
    Real dfdmu  [nbins_v*nbins_mu]; // Array to store dfdmu
    Real* dfdmu_p    = reinterpret_cast<Real*> (dfdmu);
    Real dfdmu2 [nbins_v*nbins_mu]; // Array to store dfdmumu
    Real* dfdmu2_p    = reinterpret_cast<Real*> (dfdmu2);
    Real dfdt_mu[nbins_v*nbins_mu]; // Array to store dfdt_mu
    Real* dfdt_mu_p    = reinterpret_cast<Real*> (dfdt_mu);

#define MUSPACE(var,v,mu) var##_p[(mu)*nbins_v + (v)]

    std::string PATHfile = Parameters::PADnu0;

    // Read from nu0Box.txt to get nu0 data depending on betaPara and Taniso
    std::ifstream FILEDmumu;
    FILEDmumu.open(PATHfile);

        // Get betaPara
    std::string lineBeta;
    for (int i = 0; i < 2; i++) { std::getline(FILEDmumu,lineBeta); }
         
    std::istringstream issBeta(lineBeta);
    std::vector<Real> betaParaArray;

    float numBeta;
    while ((issBeta >> numBeta)) { betaParaArray.push_back(numBeta); }

        // Get Taniso
    std::string lineTaniso;
    for (int i = 0; i < 2; i++) { std::getline(FILEDmumu,lineTaniso); }
         
    std::istringstream issTaniso(lineTaniso);
    std::vector<Real> TanisoArray;

    float numTaniso;
    while ((issTaniso >> numTaniso)) { TanisoArray.push_back(numTaniso); }
    
        // Discard one line 
    std::string lineDUMP;
    for (int i = 0; i < 1; i++) { std::getline(FILEDmumu,lineDUMP); }

        // Get nu0
    std::string linenu0;
    Real nu0Array[betaParaArray.size()][TanisoArray.size()]; 
    
    for (size_t i = 0; i < betaParaArray.size(); i++) {
        std::getline(FILEDmumu,linenu0);
        std::istringstream issnu0(linenu0);
        std::vector<Real> tempLINE;
        float numTEMP;
        while((issnu0 >> numTEMP)) { tempLINE.push_back(numTEMP); }
        for (size_t j = 0; j < tempLINE.size(); j++) {nu0Array[i][j] = tempLINE [j];}
    }



    const auto LocalCells=getLocalCells();
    #pragma omp parallel for private(fcount,fmu,dfdmu,dfdmu2,dfdt_mu)
    for (size_t CellIdx = 0; CellIdx < LocalCells.size(); CellIdx++) { //Iterate through spatial cell

        phiprof::start("Initialisation");
        auto CellID                        = LocalCells[CellIdx];
        SpatialCell& cell                  = *mpiGrid[CellID];
	const Real* parameters             = cell.get_block_parameters(popID);
        const vmesh::LocalID* nBlocks      = cell.get_velocity_grid_length(popID);
        const vmesh::MeshParameters& vMesh = getObjectWrapper().velocityMeshes[0];      

        Real density_pre_adjust  = 0.0;
        Real density_post_adjust = 0.0;

        for (size_t i=0; i<cell.get_number_of_velocity_blocks(popID)*WID3; ++i) {
            density_pre_adjust += cell.get_data(popID)[i];
        }

        Real Sparsity    = 0.01 * cell.getVelocityBlockMinValue(popID);

        Real dtTotalDiff = 0.0; // Diffusion time elapsed

        Real Vmin   = 0.0; // In case we need to avoid center cells
        Real Vmax   = 2*sqrt(3)*vMesh.meshLimits[1];
        Real dVbins = (Vmax - Vmin)/nbins_v;  
    
        std::array<Realv,VECL> dfdt;

        std::array<Real,3> bulkV = {cell.parameters[CellParams::VX], cell.parameters[CellParams::VY], cell.parameters[CellParams::VZ]};
        phiprof::stop("Initialisation");

        std::array<Real,3> B = {cell.parameters[CellParams::PERBXVOL] +  cell.parameters[CellParams::BGBXVOL],
                                 cell.parameters[CellParams::PERBYVOL] +  cell.parameters[CellParams::BGBYVOL],
	                         cell.parameters[CellParams::PERBZVOL] +  cell.parameters[CellParams::BGBZVOL]};

        Real Bnorm           = sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);
        std::array<Real,3> b = {B[0]/Bnorm, B[1]/Bnorm, B[2]/Bnorm};
  

        // There is probably a simpler way to do all of the following but I was too lazy to try to understand how to do it in C++ so I wrote everything myself.       
        // Build Pressure tensor for calculation of anisotropy and beta
        Real Ptensor[3][3] = {{cell.parameters[CellParams::P_11],cell.parameters[CellParams::P_12],cell.parameters[CellParams::P_13]},
                               {cell.parameters[CellParams::P_12],cell.parameters[CellParams::P_22],cell.parameters[CellParams::P_23]},
                               {cell.parameters[CellParams::P_13],cell.parameters[CellParams::P_23],cell.parameters[CellParams::P_33]}};

        // Build Rotation Matrix
        std::array<Real,3> eZ       = {0.0,0.0,1.0}; // Rotation around z-axis
        std::array<Real,3> BeZCross = { B[1]*eZ[2] - B[2]*eZ[1], B[2]*eZ[0] - B[0]*eZ[2], B[0]*eZ[1] - B[1]*eZ[0] };
        Real vecNorm                = sqrt(BeZCross[0]*BeZCross[0] + BeZCross[1]*BeZCross[1] + BeZCross[2]*BeZCross[2]);
        std::array<Real,3> uR       = {BeZCross[0]/vecNorm, BeZCross[1]/vecNorm, BeZCross[2]/vecNorm};
        Real beZDot                 = b[0]*eZ[0] + b[1]*eZ[1] + b[2]*eZ[2];
        Real angle                  = acos(beZDot);
        
        Real RMatrix[3][3] = { {cos(angle)+uR[0]*uR[0]*(1.0-cos(angle))      , uR[0]*uR[1]*(1.0-cos(angle))-uR[2]*sin(angle), uR[0]*uR[2]*(1.0-cos(angle))+uR[1]*sin(angle)},
                                {uR[1]*uR[0]*(1.0-cos(angle))+uR[2]*sin(angle), cos(angle)+uR[1]*uR[1]*(1.0-cos(angle))      , uR[1]*uR[2]*(1.0-cos(angle))-uR[0]*sin(angle)},
                                {uR[2]*uR[0]*(1.0-cos(angle))-uR[1]*sin(angle), uR[2]*uR[1]*(1.0-cos(angle))+uR[0]*sin(angle), cos(angle)+uR[2]*uR[2]*(1.0-cos(angle))      } };

        // Rotate Tensor (T' = RTR^{-1})
        Real RT[3][3] = { {RMatrix[0][0]*Ptensor[0][0] + RMatrix[0][1]*Ptensor[1][0] + RMatrix[0][2]*Ptensor[2][0], RMatrix[0][0]*Ptensor[0][1] + RMatrix[0][1]*Ptensor[1][1] + RMatrix[0][2]*Ptensor[2][1], RMatrix[0][0]*Ptensor[0][2] + RMatrix[0][1]*Ptensor[1][2] + RMatrix[0][2]*Ptensor[2][2]},
                           {RMatrix[1][0]*Ptensor[0][0] + RMatrix[1][1]*Ptensor[1][0] + RMatrix[1][2]*Ptensor[2][0], RMatrix[1][0]*Ptensor[0][1] + RMatrix[1][1]*Ptensor[1][1] + RMatrix[1][2]*Ptensor[2][1], RMatrix[1][0]*Ptensor[0][2] + RMatrix[1][1]*Ptensor[1][2] + RMatrix[1][2]*Ptensor[2][2]},
                           {RMatrix[2][0]*Ptensor[0][0] + RMatrix[2][1]*Ptensor[1][0] + RMatrix[2][2]*Ptensor[2][0], RMatrix[2][0]*Ptensor[0][1] + RMatrix[2][1]*Ptensor[1][1] + RMatrix[2][2]*Ptensor[2][1], RMatrix[2][0]*Ptensor[0][2] + RMatrix[2][1]*Ptensor[1][2] + RMatrix[2][2]*Ptensor[2][2]} };

        Real Rtranspose[3][3] ={ {cos(angle)+uR[0]*uR[0]*(1.0-cos(angle))      , uR[0]*uR[1]*(1.0-cos(angle))+uR[2]*sin(angle), uR[0]*uR[2]*(1.0-cos(angle))-uR[1]*sin(angle)},
                                  {uR[1]*uR[0]*(1.0-cos(angle))-uR[2]*sin(angle), cos(angle)+uR[1]*uR[1]*(1.0-cos(angle))      , uR[1]*uR[2]*(1.0-cos(angle))+uR[0]*sin(angle)},
                                  {uR[2]*uR[0]*(1.0-cos(angle))+uR[1]*sin(angle), uR[2]*uR[1]*(1.0-cos(angle))-uR[0]*sin(angle), cos(angle)+uR[2]*uR[2]*(1.0-cos(angle))      } };

        Real PtensorRotated[3][3] = {{RT[0][0]*Rtranspose[0][0] + RT[0][1]*Rtranspose[1][0] + RT[0][2]*Rtranspose[2][0], RT[0][0]*Rtranspose[0][1] + RT[0][1]*Rtranspose[1][1] + RT[0][2]*Rtranspose[2][1], RT[0][0]*Rtranspose[0][2] + RT[0][1]*Rtranspose[1][2] + RT[0][2]*Rtranspose[2][2]},
                                     {RT[1][0]*Rtranspose[0][0] + RT[1][1]*Rtranspose[1][0] + RT[1][2]*Rtranspose[2][0], RT[1][0]*Rtranspose[0][1] + RT[1][1]*Rtranspose[1][1] + RT[1][2]*Rtranspose[2][1], RT[1][0]*Rtranspose[0][2] + RT[1][1]*Rtranspose[1][2] + RT[1][2]*Rtranspose[2][2]},
                                     {RT[2][0]*Rtranspose[0][0] + RT[2][1]*Rtranspose[1][0] + RT[2][2]*Rtranspose[2][0], RT[2][0]*Rtranspose[0][1] + RT[2][1]*Rtranspose[1][1] + RT[2][2]*Rtranspose[2][1], RT[2][0]*Rtranspose[0][2] + RT[2][1]*Rtranspose[1][2] + RT[2][2]*Rtranspose[2][2]} };

        // Anisotropy
        Real Taniso = 0.5 * (PtensorRotated[0][0] + PtensorRotated[1][1]) / PtensorRotated[2][2];
        // Beta Parallel
        Real betaParallel = 2.0 * physicalconstants::MU_0 * PtensorRotated[2][2] / (Bnorm*Bnorm);

        // Find indx for first larger value
        int betaIndx   = -1;
        int TanisoIndx = -1;
        for (size_t i = 0; i < betaParaArray.size(); i++) { if (betaParallel >= betaParaArray[i]) { betaIndx   = i;} }
        for (size_t i = 0; i < TanisoArray.size()  ; i++) { if (Taniso       >= TanisoArray[i]  ) { TanisoIndx = i;} }

        Real nu0     = 0.0;
        Real epsilon = 0.0;

        if ( (betaIndx < 0) || (TanisoIndx < 0) ) {continue;}
        else { 
            if (betaIndx >= (int)betaParaArray.size()-1) {betaIndx = (int)betaParaArray.size()-2; betaParallel = betaParaArray[betaIndx+1];}
            if (TanisoIndx >= (int)TanisoArray.size()-1) {TanisoIndx = (int)TanisoArray.size()-2; Taniso = TanisoArray[TanisoIndx+1];}
            // bi-linear interpolation with weighted mean to find nu0(betaParallel,Taniso)
            Real beta1   = betaParaArray[betaIndx];
            Real beta2   = betaParaArray[betaIndx+1];
            Real Taniso1 = TanisoArray[TanisoIndx];
            Real Taniso2 = TanisoArray[TanisoIndx+1];
            Real nu011   = nu0Array[betaIndx][TanisoIndx];
            Real nu012   = nu0Array[betaIndx][TanisoIndx+1];
            Real nu021   = nu0Array[betaIndx+1][TanisoIndx];
            Real nu022   = nu0Array[betaIndx+1][TanisoIndx+1];
                // Weights
            Real w11 = (beta2 - betaParallel)*(Taniso2 - Taniso)  / ( (beta2 - beta1)*(Taniso2-Taniso1) ); 
            Real w12 = (beta2 - betaParallel)*(Taniso  - Taniso1) / ( (beta2 - beta1)*(Taniso2-Taniso1) ); 
            Real w21 = (betaParallel - beta1)*(Taniso2 - Taniso)  / ( (beta2 - beta1)*(Taniso2-Taniso1) ); 
            Real w22 = (betaParallel - beta1)*(Taniso  - Taniso1) / ( (beta2 - beta1)*(Taniso2-Taniso1) ); 
                   
            nu0 = (w11*nu011 + w12*nu012 + w21*nu021 + w22*nu022)/Parameters::PADfudge;
            cell.parameters[CellParams::NU0] = nu0;

            if (nu0 <= 0.0) {continue;}
        }      
 
        phiprof::start("Subloop");
        while (dtTotalDiff < Parameters::dt) { // Substep loop

            phiprof::start("Zeroing");
            Real RemainT  = Parameters::dt - dtTotalDiff; //Remaining time before reaching simulation time step
            Real checkCFL = std::numeric_limits<Real>::max();

            // Initialised back to zero at each substep
            memset(fmu          , 0.0, sizeof(fmu));
            memset(fcount       , 0.0, sizeof(fcount));

            phiprof::stop("Zeroing");

            phiprof::start("fmu building");
            // Build 2d array of f(v,mu)
            for (vmesh::LocalID n=0; n<cell.get_number_of_velocity_blocks(popID); n++) { // Iterate through velocity blocks

               loop_over_block([&](Veci i_indices, Veci j_indices, int k) -> void { // Witchcraft

                   //Get velocity space coordinates                    
                   const Vec VX(parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VXCRD] 
                                + (to_realv(i_indices) + 0.5)*parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVX]);

                   const Vec VY(parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VYCRD] 
                                + (to_realv(j_indices) + 0.5)*parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVY]);

                   const Vec VZ(parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VZCRD]
                                + (k + 0.5)*parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVZ]);

                   std::array<Vec,3> V = {VX,VY,VZ}; // Velocity in the cell, in the simulation frame
                                   
                   std::array<Vec,3> Vplasma; // Velocity in the cell, in the plasma frame

                   for (int indx = 0; indx < 3; indx++) { Vplasma[indx] = (V[indx] - Vec(bulkV[indx])); }
              
                   Vec normV = sqrt(Vplasma[0]*Vplasma[0] + Vplasma[1]*Vplasma[1] + Vplasma[2]*Vplasma[2]);

                   Vec Vpara = Vplasma[0]*b[0] + Vplasma[1]*b[1] + Vplasma[2]*b[2];

                   Vec mu = Vpara/(normV+std::numeric_limits<Real>::min()); // + min value to avoid division by 0
 
                   Veci Vindex;
                   Vindex = roundi(floor((normV-Vmin) / dVbins));

                   Veci muindex;
                   muindex = roundi(floor((mu+1.0) / dmubins));
                   for (uint i = 0; i<VECL; i++) {if (muindex[i] == nbins_mu) {muindex[i] == muindex[i] - 1;}} // Safety check to handle edge case where mu = exactly 1.0

                   Vec Vmu = dVbins * (to_realv(Vindex)+0.5); // Take value at the center of the mu cell

                   for (uint i = 0; i<VECL; i++) {
                       Realf CellValue = cell.get_data(n,popID)[WID*j_indices[i]+WID*WID*k+i_indices[i]];
                       //fmu_p    [muindex[i]*nbins_v + Vindex[i]] += 2.0 * M_PI * Vmu[i]*Vmu[i] * CellValue;
                       MUSPACE(fmu,Vindex[i],muindex[i]) += 2.0 * M_PI * Vmu[i]*Vmu[i] * CellValue;
                       MUSPACE(fcount,Vindex[i],muindex[i]) += 1;
                   }

                }); // End of Witchcraft
            } // End blocks
            phiprof::stop("fmu building");

            phiprof::start("space/time derivatives, CFL, Ddt");
            int cRight;
            int cLeft;
            Real checkCFLtmp = std::numeric_limits<Real>::max();

            // Compute space/time derivatives (take first non-zero neighbours) & CFL & Ddt
            for (int indv = 0; indv < nbins_v; indv++) { 
    
                // Divide f by count (independent of v but needs to be computed for all mu before derivatives)
                for(int indmu = 0; indmu < nbins_mu; indmu++) { 
                    if (MUSPACE(fcount,indv,indmu) == 0 || MUSPACE(fmu,indv,indmu) <= 0.0) { MUSPACE(fmu,indv,indmu) = std::numeric_limits<Real>::min();}
                    else {MUSPACE(fmu,indv,indmu) = MUSPACE(fmu,indv,indmu) / MUSPACE(fcount,indv,indmu);} 
                }

                for(int indmu = 0; indmu < nbins_mu; indmu++) {
                    // Compute spatial derivatives
                    if (indmu == 0) {
                        cLeft  = 0;
                        cRight = 1;
                        while( (MUSPACE(fcount,indv,indmu + cRight) == 0) && (indmu + cRight < nbins_mu-1) )  { cRight += 1; }
                        if(    (MUSPACE(fcount,indv,indmu + cRight) == 0) && (indmu + cRight == nbins_mu-1) ) { cRight  = 0; }
                    } else if (indmu == nbins_mu-1) {
                        cLeft  = 1;
                        cRight = 0;
                        while( (MUSPACE(fcount,indv,indmu - cLeft) == 0) && (indmu - cLeft > 0) )  { cLeft += 1; }
                        if(    (MUSPACE(fcount,indv,indmu - cLeft) == 0) && (indmu - cLeft == 0) ) { cLeft  = 0; }
                    } else {
                        cLeft  = 1;
                        cRight = 1;
                        while( (MUSPACE(fcount,indv,indmu + cRight) == 0) && (indmu + cRight < nbins_mu-1) )  { cRight += 1; }
                        if(    (MUSPACE(fcount,indv,indmu + cRight) == 0) && (indmu + cRight == nbins_mu-1) ) { cRight  = 0; }
                        while( (MUSPACE(fcount,indv,indmu - cLeft ) == 0) && (indmu - cLeft  > 0) )           { cLeft  += 1; }
                        if(    (MUSPACE(fcount,indv,indmu - cLeft ) == 0) && (indmu - cLeft  == 0) )          { cLeft   = 0; } 
                    } 
                    if( (cRight == 0) && (cLeft != 0) ) { 
                        MUSPACE(dfdmu ,indv,indmu) = (MUSPACE(fmu,indv,indmu + cRight) - MUSPACE(fmu,indv,indmu-cLeft))/((cRight + cLeft)*dmubins) ;
                        MUSPACE(dfdmu2,indv,indmu) = 0.0;
                    } else if( (cLeft == 0) && (cRight != 0) ) { 
                        MUSPACE(dfdmu ,indv,indmu) = (MUSPACE(fmu,indv,indmu + cRight) - MUSPACE(fmu,indv,indmu-cLeft))/((cRight + cLeft)*dmubins) ;
                        MUSPACE(dfdmu2,indv,indmu) = 0.0;
                    } else if( (cLeft == 0) && (cRight == 0) ) {
                        MUSPACE(dfdmu ,indv,indmu) = 0.0;
                        MUSPACE(dfdmu2,indv,indmu) = 0.0;
                    } else {
                        MUSPACE(dfdmu ,indv,indmu) = (  MUSPACE(fmu,indv,indmu + cRight) - MUSPACE(fmu,indv,indmu-cLeft))/((cRight + cLeft)*dmubins) ;
                        MUSPACE(dfdmu2,indv,indmu) = ( (MUSPACE(fmu,indv,indmu + cRight) - MUSPACE(fmu,indv,indmu))/(cRight*dmubins) - (MUSPACE(fmu,indv,indmu) - MUSPACE(fmu,indv,indmu-cLeft))/(cLeft*dmubins) ) / (0.5 * dmubins * (cRight + cLeft)); 
                    }

                    // Compute time derivative
                    Real mu    = (indmu+0.5)*dmubins - 1.0;
                    Real Dmumu = nu0/2.0 * ( abs(mu)/(1.0 + abs(mu)) + epsilon ) * (1.0 - mu*mu);
                    Real dDmu  = nu0/2.0 * ( (mu/abs(mu)) * ((1.0 - mu*mu)/((1.0 + abs(mu))*(1.0 + abs(mu)))) - 2.0*mu*( abs(mu)/(1.0 + abs(mu)) + epsilon));
                    MUSPACE(dfdt_mu,indv,indmu) = dDmu * MUSPACE(dfdmu,indv,indmu) + Dmumu * MUSPACE(dfdmu2,indv,indmu);

                    // Compute CFL
                    Real Vmu = dVbins * (float(indv)+0.5);
                    if (MUSPACE(fmu,indv,indmu) > Sparsity*(2.0 * M_PI * Vmu*Vmu) && abs(MUSPACE(dfdt_mu,indv,indmu)) > 0.0) { checkCFLtmp = MUSPACE(fmu,indv,indmu) * Parameters::PADCFL * (1.0/abs(MUSPACE(dfdt_mu,indv,indmu))); }
                    if (checkCFLtmp < checkCFL) { checkCFL = checkCFLtmp; }

                } // End mu loop
            } // End v loop

            // Compute Ddt
            Real Ddt = checkCFL;
            if (Ddt > RemainT) { Ddt = RemainT; }
            dtTotalDiff = dtTotalDiff + Ddt;
            phiprof::stop("space/time derivatives, CFL, Ddt");

            phiprof::start("diffusion time derivative & update cell");
            // Compute dfdt
            for (vmesh::LocalID n=0; n<cell.get_number_of_velocity_blocks(popID); n++) { // Iterate through velocity blocks             

               Real* dfdt_mu_p = reinterpret_cast<Real*> (dfdt_mu);

               loop_over_block([&](Veci i_indices, Veci j_indices, int k) -> void { // Witchcraft

                   //Get velocity space coordinates                    
                   const Vec VX(parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VXCRD] 
                                + (to_realv(i_indices) + 0.5)*parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVX]);

                   const Vec VY(parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VYCRD] 
                                + (to_realv(j_indices) + 0.5)*parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVY]);

                   const Vec VZ(parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VZCRD]
                                + (k + 0.5)*parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVZ]);
                   
                   std::array<Vec,3> V = {VX,VY,VZ}; // Velocity in the cell, in the simulation frame
                                   
                   std::array<Vec,3> Vplasma; // Velocity in the cell, in the plasma frame

                   for (int indx = 0; indx < 3; indx++) { Vplasma[indx] = (V[indx] - Vec(bulkV[indx])); }
              
                   Vec normV = sqrt(Vplasma[0]*Vplasma[0] + Vplasma[1]*Vplasma[1] + Vplasma[2]*Vplasma[2]);

                   Vec Vpara = Vplasma[0]*b[0] + Vplasma[1]*b[1] + Vplasma[2]*b[2];

                   Vec mu = Vpara/(normV+std::numeric_limits<Real>::min()); // + min value to avoid division by 0


                   Veci Vindex;
                   Vindex = roundi(floor((normV-Vmin) / dVbins));
                   
                   Veci muindex;
                   muindex = roundi(floor((mu+1.0) / dmubins));
                   for (uint i = 0; i<VECL; i++) {if (muindex[i] == nbins_mu) {muindex[i] == muindex[i] - 1;}} // Safety check to handle edge case where mu = exactly 1.0

                   Vec Vmu = dVbins * (to_realv(Vindex)+0.5);

                   for (uint i = 0; i < VECL; i++) {
                       dfdt[i] = MUSPACE(dfdt_mu,Vindex[i],muindex[i]) / (2.0 * M_PI * Vmu[i]*Vmu[i]);
                   }

                   //Update cell
                   Vec CellValue;
                   Vec NewCellValue;
                   CellValue.load(&cell.get_data(n,popID)[i_indices[0] + WID*j_indices[0] + WID*WID*k]);
                   Vec dfdtUpdate;
                   dfdtUpdate.load(&dfdt[0]);
                   NewCellValue    = CellValue + dfdtUpdate * Ddt;
                   Vecb   lessZero = NewCellValue < 0.0;
                   NewCellValue    = select(lessZero,0.0,NewCellValue);
                   NewCellValue.store(&cell.get_data(n,popID)[i_indices[0] + WID*j_indices[0] + WID*WID*k]);

                }); // End Witchcraft 
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
