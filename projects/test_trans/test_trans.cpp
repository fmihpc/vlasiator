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

using namespace std;

bool initializeProject(void) {return true;}
bool addProjectParameters(){return true;}
bool getProjectParameters(){return true;}

bool cellParametersChanged(creal& t) {return false;}

Real calcPhaseSpaceDensity(creal& x,creal& y,creal& z,creal& dx,creal& dy,creal& dz,
			   creal& vx,creal& vy,creal& vz,creal& dvx,creal& dvy,creal& dvz) {

   //Please use even number of cells in velocity and real space
   Real xyz[3];
   Real vxyz[3];
   
//location of this cell
   vxyz[0]=(vx+0.5*dvx)/dvx;
   vxyz[1]=(vy+0.5*dvy)/dvy;
   vxyz[2]=(vz+0.5*dvz)/dvz;

      
   xyz[0]=(x+0.5*dx)/dx;
   xyz[1]=(y+0.5*dy)/dy;
   xyz[2]=(z+0.5*dz)/dz;


   //real space coordinates of boxes
   //Assume an even number of spatial cells per grid dimension
   const Real box_real[8][3] = { { 1.5,1.5,1.5},
                                 {-1.5,1.5,1.5},
                                 {1.5,-1.5,1.5},
                                 {1.5,1.5,-1.5},
                                 {-1.5,-1.5,1.5},
                                 {-1.5,1.5,-1.5},
                                 {1.5,-1.5,-1.5},
                                 {-1.5,-1.5,-1.5}};
   
   //velocity space coordinates of boxes in reduced units
   //there is always an even amount of velocity cells per dimension (assuming WID is even) 
   const Real box_vel[8][3] = { { 1.5,1.5,1.5},
                                {-1.5,1.5,1.5},
                                {1.5,-1.5,1.5},
                                {1.5,1.5,-1.5},
                                {-1.5,-1.5,1.5},
                                {-1.5,1.5,-1.5},
                                {1.5,-1.5,-1.5},
                                {-1.5,-1.5,-1.5}};
   
   
   for(int box=0;box<8;box++){
      bool outsideBox=false;
      for(int i=0;i<3;i++){
         if(xyz[i]<(box_real[box][i]-0.1) ||
            xyz[i]>(box_real[box][i]+0.1) ||
            vxyz[i]<(box_vel[box][i]-0.1) ||
            vxyz[i]>(box_vel[box][i]+0.1)){
            outsideBox=true;
            break;
         }
      }
      
      if(!outsideBox) {
         return 1.0;
      }
   }
   return 0.0;
}

   
void calcBlockParameters(Real* blockParams) { }

void calcCellParameters(Real* cellParams,creal& t) {
   cellParams[CellParams::EX   ] = 0.0;
   cellParams[CellParams::EY   ] = 0.0;
   cellParams[CellParams::EZ   ] = 0.0;
   cellParams[CellParams::BX   ] = 0.0;
   cellParams[CellParams::BY   ] = 0.0;
   cellParams[CellParams::BZ   ] = 0.0;
   
   cellParams[CellParams::EXVOL] = 0.0;
   cellParams[CellParams::EYVOL] = 0.0;
   cellParams[CellParams::EZVOL] = 0.0;
   cellParams[CellParams::BXVOL] = 0.0;
   cellParams[CellParams::BYVOL] = 0.0;
   cellParams[CellParams::BZVOL] = 0.0;
}


void calcSimParameters(dccrg::Dccrg<SpatialCell>& mpiGrid, creal& t, Real& /*dt*/) {
   std::vector<uint64_t> cells = mpiGrid.get_cells();
   for (uint i = 0; i < cells.size(); ++i) {
      calcCellParameters(mpiGrid[cells[i]]->parameters, t);
   }
}

 
