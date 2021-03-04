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

void velocitySpaceDiffusion(
        dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,const uint popID){

    const auto LocalCells=getLocalCells();
    for (auto & CellID: LocalCells) {

        SpatialCell* cell     = mpiGrid[CellID];
        SpatialCell cellCopy  = *mpiGrid[CellID];

	Vec3d B(cellCopy.parameters[CellParams::PERBXVOL] +  cellCopy.parameters[CellParams::BGBXVOL],
                cellCopy.parameters[CellParams::PERBYVOL] +  cellCopy.parameters[CellParams::BGBYVOL],
	        cellCopy.parameters[CellParams::PERBZVOL] +  cellCopy.parameters[CellParams::BGBZVOL]);
        Vec3d b = normalize_vector(B);


	const Real* parameters  = cellCopy.get_block_parameters(popID);

        const vmesh::LocalID* nBlocks = cellCopy.get_velocity_grid_length(popID);

         for (vmesh::LocalID n=0; n<cellCopy.get_number_of_velocity_blocks(popID); n++) {
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
               
               Vec3d bulkV(cellCopy.parameters[CellParams::VX], cellCopy.parameters[CellParams::VY], cellCopy.parameters[CellParams::VZ]);

               Vec3d V(VX - bulkV[0] , VY - bulkV[1] , VZ - bulkV[2]);
               
               const Real DV 
                  = parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVX];
	       
	       //Norm it
               const Real normV = sqrt(V[0]*V[0]+V[1]*V[1]+V[2]*V[2]);

	       if (normV==0){ continue; }

               //Get directions
               Vec3d r     = normalize_vector(V);
               Vec3d phi   = normalize_vector(cross_product(r,b));
               Vec3d theta = normalize_vector(cross_product(r,phi));
	       
               //Get Vmesh
	       std::vector<vmesh::LocalID> Neighbors;
               const Real* MinLimits   = cellCopy.get_population(popID).vmesh.getMeshMinLimits();
               const Real* MaxLimits   = cellCopy.get_population(popID).vmesh.getMeshMaxLimits();
       
               Realf Dmumu = Parameters::PADcoefficient; //from config file
               Realf dt    = Parameters::dt;

               //Get Current cell values
	       Realf CurrentCellValue = cellCopy.get_value(VX,VY,VZ,popID); 
               Realf muSource = dot_product(r,b);


	       //Get Neighbour in +theta direction
               Vec3d VNeighbour1(VX+theta[0]*DV,VY+theta[1]*DV,VZ+theta[2]*DV);
               Vec3d CellCenter1;
               for(int i=0;i<3;i++) {
                   int cellIndex  = (VNeighbour1[i] - MinLimits[i]) / (MaxLimits[i] - MinLimits[i]) * (nBlocks[i] * 4);
                   //CellCenter1[i] = MinLimits[i] + (cellIndex + 0.5)* DV;
                   CellCenter1.insert(i, MinLimits[i] + (cellIndex + 0.5)* DV);
               }
               Realf Neighbour1 = cellCopy.get_value(CellCenter1[0],CellCenter1[1],CellCenter1[2],popID);
               Realf muTarget1  = dot_product(normalize_vector(CellCenter1-bulkV),b);
               Realf dmu1       = abs(muSource - muTarget1);
               //Realf Neighbour1 = cellCopy->get_value(VNeighbour1[0],VNeighbour1[1],VNeighbour1[2],popID);
               //Realf muTarget1  = dot_product(normalize_vector(VNeighbour1-bulkV),b);
        
               Realf diffusionAmount1 = 0.0;
               if (dmu1 != 0.0) {
                       //Diffusion parameter 1
                   Realf erfplus1         = erf((muSource + dmu1/2.0) * 1.0/sqrt(4.0 * Dmumu * dt));
                   Realf erfminus1        = erf((muSource - dmu1/2.0) * 1.0/sqrt(4.0 * Dmumu * dt));
                   diffusionAmount1 = ( 1.0/2.0 * (erfminus1 - erfplus1)) * dmu1;
               } 

	       //Get Neighbour in -theta direction
               Vec3d VNeighbour2(VX-theta[0]*DV,VY-theta[1]*DV,VZ-theta[2]*DV);
               Vec3d CellCenter2;
               for(int i=0;i<3;i++) {
                   int cellIndex  = (VNeighbour2[i] - MinLimits[i]) / (MaxLimits[i] - MinLimits[i]) * (nBlocks[i] * 4);
                   //CellCenter2[i] = MinLimits[i] + (cellIndex + 0.5)* DV;
                   CellCenter2.insert(i, MinLimits[i] + (cellIndex + 0.5)* DV);
               }
               Realf Neighbour2 = cellCopy.get_value(CellCenter2[0],CellCenter2[1],CellCenter2[2],popID);
               Realf muTarget2  = dot_product(normalize_vector(CellCenter2-bulkV),b);
               Realf dmu2       = abs(muSource - muTarget2);
               //Realf Neighbour2 = cellCopy->get_value(VNeighbour2[0],VNeighbour2[1],VNeighbour2[2],popID);
               //Realf muTarget2  = dot_product(normalize_vector(VNeighbour2-bulkV),b);

               Realf diffusionAmount2 = 0.0;
               if (dmu2 != 0.0) {
                       //Diffusion parameter 2
                   Realf erfplus2         = erf((muSource + dmu2/2.0) * 1.0/sqrt(4.0 * Dmumu * dt));
                   Realf erfminus2        = erf((muSource - dmu2/2.0) * 1.0/sqrt(4.0 * Dmumu * dt));
                   diffusionAmount2 = ( 1.0/2.0 * (erfminus2 - erfplus2)) * dmu2;
               }

               //Update current cell value
	       Realf NewCellValue = diffusionAmount1*Neighbour1 + diffusionAmount2*Neighbour2 + (1.0-(diffusionAmount1+diffusionAmount2)) * CurrentCellValue;
	       cell->set_value(VX,VY,VZ,NewCellValue,popID); 

	    }


	 }
    }
}


