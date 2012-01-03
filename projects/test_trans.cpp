#include <cstdlib>
#include <iostream>
#include <cmath>

#include "spatial_cell.hpp"
#include "common.h"
#include "project.h"
#include "parameters.h"

using namespace std;

bool initializeProject(void) {
   return true;
}

bool cellParametersChanged(creal& t) {return false;}

Real calcPhaseSpaceDensity(creal& x,creal& y,creal& z,creal& dx,creal& dy,creal& dz,
			   creal& vx,creal& vy,creal& vz,creal& dvx,creal& dvy,creal& dvz) {

   Real xyz[3];
   Real vxyz[3];
   Real dxyz[3];
   Real dvxyz[3];

   dxyz[0]=dx;
   dxyz[1]=dy;
   dxyz[2]=dz;
   dvxyz[0]=dvx;
   dvxyz[1]=dvy;
   dvxyz[2]=dvz;
   xyz[0]=x;
   xyz[1]=y;
   xyz[2]=z;
   vxyz[0]=vx;
   vxyz[1]=vy;
   vxyz[2]=vz;

   //in coordinate units of dx,dy,dz
   const Real box_real[8][3] = { { 2,2,2},
                                 {-2,2,2},
                                 {2,-2,2},
                                 {2,2,-2},
                                 {-2,-2,2},
                                 {-2,2,-2},
                                 {2,-2,-2},
                                 {-2,-2,-2}};
   
   //in coordinate units of dvx,dvy,dvz
   const Real box_vel[8][3] = { { 1,1,1},
                                {-1,1,1},
                                {1,-1,1},
                                {1,1,-1},
                                {-1,-1,1},
                                {-1,1,-1},
                                {1,-1,-1},
                                {-1,-1,-1}};
   
   
   for(int box=0;box<8;box++){
      bool outsideBox=false;
      for(int i=0;i<3;i++){
         if(xyz[i]<(box_real[box][i]-0.5)*dxyz[i] ||
            xyz[i]>(box_real[box][i]+0.5)*dxyz[i] ||
            vxyz[i]<(box_vel[box][i]-0.5)*dvxyz[i] ||
            vxyz[i]>(box_vel[box][i]+0.5)*dvxyz[i]){
            outsideBox=true;
            break;
         }
      }
      
      if(!outsideBox) return 1.0;
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

 
