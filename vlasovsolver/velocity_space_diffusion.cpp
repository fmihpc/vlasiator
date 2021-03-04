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

               Vec3d V(VX - bulkV[0] , VY - bulkV[1] , VZ - bulkV[2]); //Velocity in the cell, in the plasma frame
               
               const Real DV 
                  = parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVX];
	       
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
               Realf OXCoeff         = abs(theta[0]); //How much to diffuse through the interface based on the dot product between
               Realf OYCoeff         = abs(theta[1]); //the theta vector and the normal to the interface (Same goes for the other cell)
               Realf OZCoeff         = abs(theta[2]); //Number between 0 and 1
                   //Get cell extents in mu
                       //Along X
               Vec3d OXPlus(VX+DV/2.0,VY,VZ);
               Vec3d OXMinus(VX-DV/2.0,VY,VZ);
               Realf muOXMax  = dot_product(normalize_vector(OXPlus-bulkV),b);
               Realf muOXMin  = dot_product(normalize_vector(OXMinus-bulkV),b);
               Realf OXlength = abs(muOXMax - muOXMin);
                       //Along Y
               Vec3d OYPlus(VX,VY+DV/2.0,VZ);
               Vec3d OYMinus(VX,VY-DV/2.0,VZ);
               Realf muOYMax  = dot_product(normalize_vector(OYPlus-bulkV),b);
               Realf muOYMin  = dot_product(normalize_vector(OYMinus-bulkV),b);
               Realf OYlength = abs(muOYMax - muOYMin);
                       //Along Z
               Vec3d OZPlus(VX,VY,VZ+DV/2.0);
               Vec3d OZMinus(VX,VY,VZ-DV/2.0);
               Realf muOZMax  = dot_product(normalize_vector(OZPlus-bulkV),b);
               Realf muOZMin  = dot_product(normalize_vector(OZMinus-bulkV),b);
               Realf OZlength = abs(muOZMax - muOZMin); 
 
               //Get +X cell values
               Realf XCellValue = cellCopy.get_value(VX+DV,VY,VZ,popID);
               Vec3d XVelocity(VX+DV-bulkV[0],VY-bulkV[1],VZ-bulkV[2]);
               Vec3d Xr      = normalize_vector(XVelocity);
               Vec3d Xphi    = normalize_vector(cross_product(Xr,b));
               Vec3d Xtheta  = normalize_vector(cross_product(Xr,Xphi));
               Realf Xmu     = dot_product(Xr,b);
               Realf XOCoeff = abs(Xtheta[0]); 
               Realf dmuX    = abs(muOrigin - Xmu);
                   //Get cell extents in mu
               Vec3d XOPlus(VX+DV*(3.0/2.0),VY,VZ);
               Vec3d XOMinus(VX+DV/2.0,VY,VZ);
               Realf muXOMax = dot_product(normalize_vector(XOPlus-bulkV),b);
               Realf muXOMin = dot_product(normalize_vector(XOMinus-bulkV),b);
               Realf Xlength = abs(muXOMax - muXOMin);
                   //Get diffusion from +X to Origin
               Realf erfplusX  = erf((Xmu + Xlength/2.0)*(1.0/sqrt(4.0*Dmumu*dt)));
               Realf erfminusX = erf((Xmu - Xlength/2.0)*(1.0/sqrt(4.0*Dmumu*dt)));
               Realf diffXO    = XOCoeff * 1.0/2.0 * (erfminusX - erfplusX) * dmuX; 
                   //Get diffusion from Origin to +X
               Realf erfplusOX  = erf((muOrigin + OXlength/2.0)*(1.0/sqrt(4.0*Dmumu*dt)));
               Realf erfminusOX = erf((muOrigin - OXlength/2.0)*(1.0/sqrt(4.0*Dmumu*dt)));
               Realf diffOX     = OXCoeff * 1.0/2.0 * (erfminusOX - erfplusOX) * dmuX; 

               //Get +Y cell values
               Realf YCellValue = cellCopy.get_value(VX,VY+DV,VZ,popID);
               Vec3d YVelocity(VX-bulkV[0],VY+DV-bulkV[1],VZ-bulkV[2]);
	       Vec3d Yr      = normalize_vector(YVelocity);
	       Vec3d Yphi    = normalize_vector(cross_product(Yr,b));
	       Vec3d Ytheta  = normalize_vector(cross_product(Yr,Yphi));
	       Realf Ymu     = dot_product(Yr,b);
               Realf YOCoeff = abs(Ytheta[1]); 
               Realf dmuY    = abs(muOrigin - Ymu);
                   //Get cell extents in mu
               Vec3d YOPlus(VX,VY+DV*(3.0/2.0),VZ);
               Vec3d YOMinus(VX,VY+DV/2.0,VZ);
               Realf muYOMax = dot_product(normalize_vector(YOPlus-bulkV),b);
               Realf muYOMin = dot_product(normalize_vector(YOMinus-bulkV),b);
               Realf Ylength = abs(muYOMax - muYOMin);
                   //Get diffusion from +Y to Origin
               Realf erfplusY  = erf((Ymu + Ylength/2.0)*(1.0/sqrt(4.0*Dmumu*dt)));
               Realf erfminusY = erf((Ymu - Ylength/2.0)*(1.0/sqrt(4.0*Dmumu*dt)));
               Realf diffYO    = YOCoeff * 1.0/2.0 * (erfminusY - erfplusY) * dmuY; 
                   //Get diffusion from Origin to +Y
               Realf erfplusOY  = erf((muOrigin + OYlength/2.0)*(1.0/sqrt(4.0*Dmumu*dt)));
               Realf erfminusOY = erf((muOrigin - OYlength/2.0)*(1.0/sqrt(4.0*Dmumu*dt)));
               Realf diffOY     = OYCoeff * 1.0/2.0 * (erfminusOY - erfplusOY) * dmuY; 

               //Get +Z cell values
               Realf ZCellValue = cellCopy.get_value(VX,VY,VZ+DV,popID);
               Vec3d ZVelocity(VX-bulkV[0],VY-bulkV[1],VZ+DV-bulkV[2]);
	       Vec3d Zr      = normalize_vector(ZVelocity);
	       Vec3d Zphi    = normalize_vector(cross_product(Zr,b));
	       Vec3d Ztheta  = normalize_vector(cross_product(Zr,Zphi));
	       Realf Zmu     = dot_product(Zr,b);
               Realf ZOCoeff = abs(Ztheta[2]); 
               Realf dmuZ    = abs(muOrigin - Zmu);
                   //Get cell extents in mu
               Vec3d ZOPlus(VX,VY,VZ+DV*(3.0/2.0));
               Vec3d ZOMinus(VX,VY,VZ+DV/2.0);
               Realf muZOMax = dot_product(normalize_vector(ZOPlus-bulkV),b);
               Realf muZOMin = dot_product(normalize_vector(ZOMinus-bulkV),b);
               Realf Zlength = abs(muZOMax - muZOMin);
                   //Get diffusion from +Z to Origin
               Realf erfplusZ  = erf((Zmu + Zlength/2.0)*(1.0/sqrt(4.0*Dmumu*dt)));
               Realf erfminusZ = erf((Zmu - Zlength/2.0)*(1.0/sqrt(4.0*Dmumu*dt)));
               Realf diffZO    = ZOCoeff * 1.0/2.0 * (erfminusZ - erfplusZ) * dmuZ; 
                   //Get diffusion from Origin to +Z
               Realf erfplusOZ  = erf((muOrigin + OZlength/2.0)*(1.0/sqrt(4.0*Dmumu*dt)));
               Realf erfminusOZ = erf((muOrigin - OZlength/2.0)*(1.0/sqrt(4.0*Dmumu*dt)));
               Realf diffOZ     = OZCoeff * 1.0/2.0 * (erfminusOZ - erfplusOZ) * dmuZ; 

               
               std::cerr << "diffXO = " << diffXO << std::endl;
               std::cerr << "diffYO = " << diffYO << std::endl;
               std::cerr << "diffZO = " << diffZO << std::endl;
               std::cerr << "diffOX = " << diffOX << std::endl;
               std::cerr << "diffOY = " << diffOY << std::endl;
               std::cerr << "diffOZ = " << diffOZ << std::endl;

               //Update origin cell value
	       Realf DeltaOriginValue = diffXO*XCellValue + diffYO*YCellValue + diffZO*ZCellValue - (diffOX + diffOY + diffOZ)*OriginCellValue;
	       cell->increment_value(VX,VY,VZ,DeltaOriginValue,popID); 

               //Update +X cell value
               Realf DeltaXValue = - diffXO*XCellValue + diffOX*OriginCellValue;
               cell->increment_value(VX+DV,VY,VZ,DeltaXValue,popID);
               //Update +Y cell value
               Realf DeltaYValue = - diffYO*YCellValue + diffOY*OriginCellValue;
               cell->increment_value(VX,VY+DV,VZ,DeltaYValue,popID);
               //Update +Z cell value
               Realf DeltaZValue = - diffZO*ZCellValue + diffOZ*OriginCellValue;
               cell->increment_value(VX,VY,VZ+DV,DeltaZValue,popID);

	    }


	 }
    }
}


