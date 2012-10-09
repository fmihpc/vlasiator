/*
This file is part of Vlasiator.

Copyright 2011, 2012 Finnish Meteorological Institute

Vlasiator is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3
as published by the Free Software Foundation.

Vlasiator is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Vlasiator. If not, see <http://www.gnu.org/licenses/>.
*/

#include <cstdlib>
#include <iostream>
#include <cmath>

#include "spatial_cell.hpp"
#include "common.h"
#include "project.h"
#include "parameters.h"
#include "readparameters.h"
#include "vlasovmover.h"

enum cases {BXCASE,BYCASE,BZCASE};

// SET THIS TO BXCASE,BYCASE,OR BZCASE TO SELECT ONE OF THE THREE CASES

static int CASE = BZCASE;

using namespace std;

bool initializeProject(void) {return true;}
bool addProjectParameters(){return true;}
bool getProjectParameters(){return true;}

void setProjectCell(SpatialCell* cell) {
   // Set up cell parameters:
   calcCellParameters(&((*cell).parameters[0]), 0.0);
   
   cell->parameters[CellParams::RHOLOSSADJUST] = 0.0;
   cell->parameters[CellParams::RHOLOSSVELBOUNDARY] = 0.0;
   
   // Go through each velocity block in the velocity phase space grid.
   // Set the initial state and block parameters:
   creal dvx_block = SpatialCell::block_dvx; // Size of a block in vx-direction
   creal dvy_block = SpatialCell::block_dvy; //                    vy
   creal dvz_block = SpatialCell::block_dvz; //                    vz
   creal dvx_blockCell = SpatialCell::cell_dvx; // Size of one cell in a block in vx-direction
   creal dvy_blockCell = SpatialCell::cell_dvy; //                                vy
   creal dvz_blockCell = SpatialCell::cell_dvz; //                                vz
   
   for (uint kv=0; kv<P::vzblocks_ini; ++kv) 
      for (uint jv=0; jv<P::vyblocks_ini; ++jv)
         for (uint iv=0; iv<P::vxblocks_ini; ++iv) {
            creal vx_block = P::vxmin + iv*dvx_block; // vx-coordinate of the lower left corner
            creal vy_block = P::vymin + jv*dvy_block; // vy-
            creal vz_block = P::vzmin + kv*dvz_block; // vz-
            
            // Calculate volume average of distrib. function for each cell in the block.
            for (uint kc=0; kc<WID; ++kc) 
               for (uint jc=0; jc<WID; ++jc) 
                  for (uint ic=0; ic<WID; ++ic) {
                     creal vx_cell = vx_block + ic*dvx_blockCell;
                     creal vy_cell = vy_block + jc*dvy_blockCell;
                     creal vz_cell = vz_block + kc*dvz_blockCell;
                     Real average = 
                     calcPhaseSpaceDensity(cell->parameters[CellParams::XCRD],
                                           cell->parameters[CellParams::YCRD],
                                           cell->parameters[CellParams::ZCRD],
                                           cell->parameters[CellParams::DX],
                                           cell->parameters[CellParams::DY],
                                           cell->parameters[CellParams::DZ],
                                           vx_cell,vy_cell,vz_cell,
                                           dvx_blockCell,dvy_blockCell,dvz_blockCell);
                     
                     if(average!=0.0){
                        creal vx_cell_center = vx_block + (ic+convert<Real>(0.5))*dvx_blockCell;
                        creal vy_cell_center = vy_block + (jc+convert<Real>(0.5))*dvy_blockCell;
                        creal vz_cell_center = vz_block + (kc+convert<Real>(0.5))*dvz_blockCell;
                        cell->set_value(vx_cell_center,vy_cell_center,vz_cell_center,average);
                     }
                  }
         }
         calculateCellVelocityMoments(cell);
         
         //let's get rid of blocks not fulfilling the criteria here to save memory.
         cell->adjustSingleCellVelocityBlocks();
}

bool cellParametersChanged(creal& t) {return false;}

Real sign(creal value)
{
   if(abs(value) < 1e-5) return 0.0;
   else return value / abs(value);
}

Real calcPhaseSpaceDensity(creal& x,creal& y,creal& z,creal& dx,creal& dy,creal& dz,creal& vx,creal& vy,creal& vz,creal& dvx,creal& dvy,creal& dvz) {
   const Real T  = 1e-6;
   const Real n  = 1.0e7;
   creal ALPHA = 7.0 * M_PI / 4.0;

   Real VX,VY,VZ,ksi,eta;
   switch (CASE) {
      case BXCASE:
      ksi = ((y + 0.5 * dy)  * cos(ALPHA) + (z + 0.5 * dz) * sin(ALPHA)) / (2.0 * sqrt(2.0));
      eta = (-(y + 0.5 * dy)  * sin(ALPHA) + (z + 0.5 * dz) * cos(ALPHA)) / (2.0 * sqrt(2.0));
      VX = 0.0;
      VY = sign(cos(ALPHA)) * 0.5 + 0.1*cos(ALPHA) * sin(2.0 * M_PI * eta);
      VZ = sign(sin(ALPHA)) * 0.5 + 0.1*sin(ALPHA) * sin(2.0 * M_PI * eta); 
      break;
    case BYCASE:
      ksi = ((z + 0.5 * dz)  * cos(ALPHA) + (x + 0.5 * dx) * sin(ALPHA)) / (2.0 * sqrt(2.0));
      eta = (-(z + 0.5 * dz)  * sin(ALPHA) + (x + 0.5 * dx) * cos(ALPHA)) / (2.0 * sqrt(2.0));
      VX = sign(sin(ALPHA)) * 0.5 + 0.1*sin(ALPHA) * sin(2.0 * M_PI * eta);
      VY = 0.0;
      VZ = sign(cos(ALPHA)) * 0.5 + 0.1*cos(ALPHA) * sin(2.0 * M_PI * eta);
      break;
    case BZCASE:
      ksi = ((x + 0.5 * dx)  * cos(ALPHA) + (y + 0.5 * dy) * sin(ALPHA)) / (2.0 * sqrt(2.0));
      eta = (-(x + 0.5 * dx)  * sin(ALPHA) + (y + 0.5 * dy) * cos(ALPHA)) / (2.0 * sqrt(2.0));
      VX = sign(cos(ALPHA)) * 0.5 + 0.1*cos(ALPHA) * sin(2.0 * M_PI * eta);
      VY = sign(sin(ALPHA)) * 0.5 + 0.1*sin(ALPHA) * sin(2.0 * M_PI * eta);
      VZ = 0.0;
      break;
   }
   
   creal VX2 = (vx+0.5*dvx-VX)*(vx+0.5*dvx-VX);
   creal VY2 = (vy+0.5*dvy-VY)*(vy+0.5*dvy-VY);
   creal VZ2 = (vz+0.5*dvz-VZ)*(vz+0.5*dvz-VZ);
   
   creal CONST = physicalconstants::MASS_PROTON / 2.0 / physicalconstants::K_B / T;
   Real NORM = (physicalconstants::MASS_PROTON / 2.0 / M_PI / physicalconstants::K_B / T);
   NORM = n * pow(NORM,1.5);
   
   return NORM*exp(-CONST*(VX2+VY2+VZ2));
}

void calcBlockParameters(Real* blockParams) { }

void calcCellParameters(Real* cellParams,creal& t) {
   cellParams[CellParams::EX   ] = 0.0;
   cellParams[CellParams::EY   ] = 0.0;
   cellParams[CellParams::EZ   ] = 0.0;
   cellParams[CellParams::PERBX   ] = 0.0;
   cellParams[CellParams::PERBY   ] = 0.0;
   cellParams[CellParams::PERBZ   ] = 0.0;
   cellParams[CellParams::BGBX] = 0;
   cellParams[CellParams::BGBY] = 0;
   cellParams[CellParams::BGBZ] = 0;
   
   typedef Parameters P;
   creal x = cellParams[CellParams::XCRD] + 0.5 * cellParams[CellParams::DX];
   creal y = cellParams[CellParams::YCRD] + 0.5 * cellParams[CellParams::DY];
   creal z = cellParams[CellParams::ZCRD] + 0.5 * cellParams[CellParams::DZ];

   creal B0 = 1e-10;
   creal B1 = 1e-9;
   
   switch (CASE) {
    case BXCASE:
      if (y >= -0.2 && y <= 0.2)
	if (z >= -0.2 && z <= 0.2)
	  cellParams[CellParams::PERBX] = B1;
	  cellParams[CellParams::BGBX] = B0;
      break;
    case BYCASE:
      if (x >= -0.2 && x <= 0.2)
	if (z >= -0.2 && z <= 0.2)
	  cellParams[CellParams::PERBY] = B1;
	  cellParams[CellParams::BGBY] = B0;
      break;
    case BZCASE:
      if (x >= -0.2 && x <= 0.2)
	if (y >= -0.2 && y <= 0.2)
	  cellParams[CellParams::PERBZ] = B1;
	  cellParams[CellParams::BGBZ] = B0;
      break;
   }
}



void calcSimParameters(dccrg::Dccrg<SpatialCell>& mpiGrid, creal& t, Real& /*dt*/) {
   std::vector<uint64_t> cells = mpiGrid.get_cells();
   for (uint i = 0; i < cells.size(); ++i) {
      calcCellParameters(mpiGrid[cells[i]]->parameters, t);
   }
}



