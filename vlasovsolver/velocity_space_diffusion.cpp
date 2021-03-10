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

        for (int coord = 0; coord < 3; coord++) {

           SpatialCell* cell     = mpiGrid[CellID];
           SpatialCell cellCopy  = *mpiGrid[CellID];

	   Vec3d B(cellCopy.parameters[CellParams::PERBXVOL] +  cellCopy.parameters[CellParams::BGBXVOL],
                   cellCopy.parameters[CellParams::PERBYVOL] +  cellCopy.parameters[CellParams::BGBYVOL],
	           cellCopy.parameters[CellParams::PERBZVOL] +  cellCopy.parameters[CellParams::BGBZVOL]);
           Vec3d b = normalize_vector(B);


	   const Real* parameters  = cellCopy.get_block_parameters(popID);

           const vmesh::LocalID* nBlocks = cellCopy.get_velocity_grid_length(popID);

            for (vmesh::LocalID n=0; n<cellCopy.get_number_of_velocity_blocks(popID); n++) { //Iterate through velocity blocks
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

                  Vec3d V(VX - bulkV[0] , VY - bulkV[1] , VZ - bulkV[2]); //Velocity in the cell, in the plasma frame
                  Vec3d NeighbourV(VX,VY,VZ);
                 
                  const Real DV 
                     = parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVX];
	
                  Real DVarray[3] = {0,0,0};
                  DVarray[coord] = DV;
                  Vec3d VDV;
                  VDV.load(DVarray);

                  NeighbourV += VDV;              
           
	          //Norm it
                  const Real normV = sqrt(V[0]*V[0]+V[1]*V[1]+V[2]*V[2]);

	          if (normV==0){ continue; }

	          //Get Vmesh
	          std::vector<vmesh::LocalID> Neighbors;
       
                  Realf Dmumu = Parameters::PADcoefficient; //from config file
                  Realf dt    = Parameters::dt; //Overwritten by the simulation dt

                  //Get origin cell values
	          Realf OriginCellValue = cellCopy.get_value(VX,VY,VZ,popID);
                  Vec3d r               = normalize_vector(V);
                  Vec3d phi             = normalize_vector(cross_product(r,b));
                  Vec3d theta           = normalize_vector(cross_product(r,phi));
                  Realf muOrigin        = dot_product(r,b);
                  Realf ONCoeff         = abs(theta[coord]); //How much to diffuse through the interface based on the dot product between

                      //Get cell extents in mu
                          //Along X
                  Vec3d ONPlus(VX,VY,VZ);
                  ONPlus += VDV/2.0;
                  Vec3d ONMinus(VX,VY,VZ);
                  ONMinus -= VDV/2.0;
                  Realf muONMax  = dot_product(normalize_vector(ONPlus-bulkV),b);
                  Realf muONMin  = dot_product(normalize_vector(ONMinus-bulkV),b);
                  Realf ONlength = abs(muONMax - muONMin);
 
                  Realf diffON = 0.0;
                  Realf diffNO = 0.0;
                  Realf NCellValue = 0.0;

                  if (checkExistingNeighbour(&cellCopy,NeighbourV[0],NeighbourV[1],NeighbourV[2],popID)) {
                      //Get +X cell values
                      NCellValue = cellCopy.get_value(NeighbourV[0],NeighbourV[1],NeighbourV[2],popID);
                      Vec3d NVelocity = NeighbourV - bulkV;
                      Vec3d Nr      = normalize_vector(NVelocity);
                      Vec3d Nphi    = normalize_vector(cross_product(Nr,b));
                      Vec3d Ntheta  = normalize_vector(cross_product(Nr,Nphi));
                      Realf Nmu     = dot_product(Nr,b);
                      Realf NOCoeff = abs(Ntheta[coord]); 
                      Realf dmuN    = abs(muOrigin - Nmu);
                          //Get cell extents in mu
                      Vec3d NOPlus(VX,VY,VZ);
                      NOPlus += VDV*3.0/2.0;
                      Vec3d NOMinus(VX,VY,VZ);
                      NOMinus += VDV/2.0;
                      Realf muNOMax = dot_product(normalize_vector(NOPlus-bulkV),b);
                      Realf muNOMin = dot_product(normalize_vector(NOMinus-bulkV),b);
                      Realf Nlength = abs(muNOMax - muNOMin);
                          //Get diffusion from +X to Origin
                      Realf erfplusN  = erf((Nmu + Nlength/2.0)*(1.0/sqrt(4.0*Dmumu*dt)));
                      Realf erfminusN = erf((Nmu - Nlength/2.0)*(1.0/sqrt(4.0*Dmumu*dt)));
                      diffNO          = NOCoeff * 1.0/2.0 * (erfminusN - erfplusN) * dmuN; 
                          //Get diffusion from Origin to +X
                      Realf erfplusON  = erf((muOrigin + ONlength/2.0)*(1.0/sqrt(4.0*Dmumu*dt)));
                      Realf erfminusON = erf((muOrigin - ONlength/2.0)*(1.0/sqrt(4.0*Dmumu*dt)));
                      diffON           = ONCoeff * 1.0/2.0 * (erfminusON - erfplusON) * dmuN; 
                      //Update +X cell value
                      Realf DeltaNValue = - diffNO*NCellValue + diffON*OriginCellValue;
                      cell->increment_value(NeighbourV[0],NeighbourV[1],NeighbourV[2],DeltaNValue,popID);

                      //Update origin cell value
	              Realf DeltaOriginValue = diffNO*NCellValue - diffON*OriginCellValue;      
                      cell->increment_value(VX,VY,VZ,DeltaOriginValue,popID); 

                      std::cerr << "For V(" << V[0] << "," << V[1] << "," << V[2] << ")" << std::endl;
                      std::cerr << "To VNeighbour(" << NeighbourV[0] << "," << NeighbourV[1] << "," << NeighbourV[2] << ")" << std::endl;
                      std::cerr << "theta(" << theta[0] << "," << theta[1] << "," << theta[2] << ")" << std::endl;
                      std::cerr << "Ntheta(" << Ntheta[0] << "," << Ntheta[1] << "," << Ntheta[2] << ")" << std::endl;
                      std::cerr << "muOrigin = " << muOrigin << std::endl;
                      std::cerr << "Nmu = " << Nmu << std::endl;
                      std::cerr << "diffNO = " << diffNO << std::endl;
                      std::cerr << "diffON = " << diffON << std::endl;
                      std::cerr << "OriginCellValue = " << OriginCellValue << std::endl;
                      std::cerr << "NCellValue = " << NCellValue << std::endl;
                      std::cerr << "DeltaNValue = " << DeltaNValue << std::endl;
                      std::cerr << "DeltaOriginValue = " << DeltaOriginValue << std::endl;
                  }
                  
              }
           } 

	 }
    }
}


