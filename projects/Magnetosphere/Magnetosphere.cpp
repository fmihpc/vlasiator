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

using namespace std;

typedef magnetosphereParameters MP;
Real MP::T = {NAN};
Real MP::rho = NAN;
uint MP::nSpaceSamples = 0;
uint MP::nVelocitySamples = 0;



/** Integrate the distribution function over the given six-dimensional phase-space cell.
 * @param x Starting value of the x-coordinate of the cell.
 * @param y Starting value of the y-coordinate of the cell.
 * @param z Starting value of the z-coordinate of the cell.
 * @param dx The size of the cell in x-direction.
 * @param dy The size of the cell in y-direction.
 * @param dz The size of the cell in z-direction.
 * @param vx Starting value of the vx-coordinate of the cell.
 * @param vy Starting value of the vy-coordinate of the cell.
 * @param vz Starting value of the vz-coordinate of the cell.
 * @param dvx The size of the cell in vx-direction.
 * @param dvy The size of the cell in vy-direction.
 * @param dvz The size of the cell in vz-direction.
 * @return The volume average of the distribution function in the given phase space cell.
 * The physical unit of this quantity is 1 / (m^3 (m/s)^3).
 */
Real calcPhaseSpaceDensity(creal& x,creal& y,creal& z,creal& dx,creal& dy,creal& dz,
			   creal& vx,creal& vy,creal& vz,creal& dvx,creal& dvy,creal& dvz);

/** Calculate parameters for the given spatial cell at the given time.
 * Here you need to set values for the following array indices:
 * CellParams::EX, CellParams::EY, CellParams::EZ, CellParams::BX, CellParams::BY, and CellParams::BZ.
 * 
 * The following array indices contain the coordinates of the "lower left corner" of the cell: 
 * CellParams::XCRD, CellParams::YCRD, and CellParams::ZCRD.
 * The cell size is given in the following array indices: CellParams::DX, CellParams::DY, and CellParams::DZ.
 * @param cellParams Array containing cell parameters.
 * @param t The current value of time. This is passed as a convenience. If you need more detailed information 
 * of the state of the simulation, you can read it from Parameters.
 */
void calcCellParameters(Real* cellParams,creal& t);



bool initializeProject(void) {return true;}

bool addProjectParameters(){
   typedef Readparameters RP;
   RP::add("Magnetosphere.rho", "Number density (m^-3)", 0.0);
   RP::add("Magnetosphere.T", "Temperature (K)", 0.0);
   RP::add("Magnetosphere.nSpaceSamples", "Number of sampling points per spatial dimension", 2);
   RP::add("Magnetosphere.nVelocitySamples", "Number of sampling points per velocity dimension", 5);
   return true;
}

bool getProjectParameters(){
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
   typedef Readparameters RP;
   if(!RP::get("Magnetosphere.rho", MP::rho)) {
      if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
      exit(1);
   }
   if(!RP::get("Magnetosphere.T", MP::T)) {
      if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
      exit(1);
   }
   if(!RP::get("Magnetosphere.nSpaceSamples", MP::nSpaceSamples)) {
      if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
      exit(1);
   }
   if(!RP::get("Magnetosphere.nVelocitySamples", MP::nVelocitySamples)) {
      if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
      exit(1);
   }
   return true;
}

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

void dipole(creal x, creal y, creal z, Real& Bx, Real &By, Real& Bz) {
   creal k_0 = 8.0e15; // Wb m
   Real r = sqrt(x*x + y*y + z*z); // radial
   Real theta = atan2(sqrt(x*x + y*y), z); // polar
   Real phi = atan2(y, x); // azimuthal
   
   Bx = sin(theta) * cos(theta) * cos(phi);
   By = sin(theta) * cos(theta) * sin(phi);
   Bz = cos(theta)*cos(theta) - 1.0 / 3.0;
   
   Bx *= 3.0 * k_0 / (r*r*r);
   By *= 3.0 * k_0 / (r*r*r);
   Bz *= 3.0 * k_0 / (r*r*r);
}

Real getDistribValue(creal& x,creal& y, creal& z, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz) {
   return MP::rho * pow(physicalconstants::MASS_PROTON / (2.0 * M_PI * physicalconstants::K_B * MP::T), 1.5) *
   exp(- physicalconstants::MASS_PROTON * (vx*vx + vy*vy + vz*vz) / (2.0 * physicalconstants::K_B * MP::T));
}

Real calcPhaseSpaceDensity(creal& x, creal& y, creal& z, creal& dx, creal& dy, creal& dz, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz) {   
   if((MP::nSpaceSamples > 1) && (MP::nVelocitySamples > 1)) {
      creal d_x = dx / (MP::nSpaceSamples-1);
      creal d_y = dy / (MP::nSpaceSamples-1);
      creal d_z = dz / (MP::nSpaceSamples-1);
      creal d_vx = dvx / (MP::nVelocitySamples-1);
      creal d_vy = dvy / (MP::nVelocitySamples-1);
      creal d_vz = dvz / (MP::nVelocitySamples-1);
      
      Real avg = 0.0;
   // #pragma omp parallel for collapse(6) reduction(+:avg)
      // WARNING No threading here if calling functions are already threaded
      for (uint i=0; i<MP::nSpaceSamples; ++i)
         for (uint j=0; j<MP::nSpaceSamples; ++j)
            for (uint k=0; k<MP::nSpaceSamples; ++k)
               for (uint vi=0; vi<MP::nVelocitySamples; ++vi)
                  for (uint vj=0; vj<MP::nVelocitySamples; ++vj)
                     for (uint vk=0; vk<MP::nVelocitySamples; ++vk) {
                        avg += getDistribValue(x+i*d_x, y+j*d_y, z+k*d_z, vx+vi*d_vx, vy+vj*d_vy, vz+vk*d_vz, dvx, dvy, dvz);
                     }
                     return avg /
                            (MP::nSpaceSamples*MP::nSpaceSamples*MP::nSpaceSamples) /
                            (MP::nVelocitySamples*MP::nVelocitySamples*MP::nVelocitySamples);
   } else {
      return getDistribValue(x+0.5*dx, y+0.5*dy, z+0.5*dz, vx+0.5*dvx, vy+0.5*dvy, vz+0.5*dvz, dvx, dvy, dvz);
   }
   
//    CellID cellID = 1 + round((x - Parameters::xmin) / dx + 
//    (y - Parameters::ymin) / dy * Parameters::xcells_ini +
//    (z - Parameters::zmin) / dz * Parameters::ycells_ini * Parameters::xcells_ini);
   
//    return cellID * pow(physicalconstants::MASS_PROTON / (2.0 * M_PI * physicalconstants::K_B * MP::T), 1.5) *
//    exp(- physicalconstants::MASS_PROTON * (vx*vx + vy*vy + vz*vz) / (2.0 * physicalconstants::K_B * MP::T));
}


void calcCellParameters(Real* cellParams,creal& t) {
   creal x = cellParams[CellParams::XCRD];
   creal dx = cellParams[CellParams::DX];
   creal y = cellParams[CellParams::YCRD];
   creal dy = cellParams[CellParams::DY];
   creal z = cellParams[CellParams::ZCRD];
   creal dz = cellParams[CellParams::DZ];
   
   Real Bxavg, Byavg, Bzavg;
   if((MP::nSpaceSamples > 1) && (MP::nVelocitySamples > 1)) {
      Bxavg = Byavg = Bzavg = 0.0;
      Real d_x = dx / (MP::nSpaceSamples - 1);
      Real d_y = dy / (MP::nSpaceSamples - 1);
      Real d_z = dz / (MP::nSpaceSamples - 1);
   // #pragma omp parallel for collapse(3) reduction(+:Bxavg,Byavg,Bzavg)
      // WARNING No threading here if calling functions are already threaded
      for (uint i=0; i<MP::nSpaceSamples; ++i)
         for (uint j=0; j<MP::nSpaceSamples; ++j)
            for (uint k=0; k<MP::nSpaceSamples; ++k) {
               Real field[3];
               dipole(x+i*d_x, y+j*d_y, z+k*d_z, field[0], field[1], field[2]);
               Bxavg += field[0];
               Byavg += field[1];
               Bzavg += field[2];
            }
      cuint nPts = MP::nSpaceSamples*MP::nSpaceSamples*MP::nSpaceSamples;
      
      cellParams[CellParams::BGBX   ] = Bxavg / nPts;
      cellParams[CellParams::BGBY   ] = Byavg / nPts;
      cellParams[CellParams::BGBZ   ] = Bzavg / nPts;
   } else {
      dipole(x+0.5*dx, y+0.5*dy, z+0.5*dz, Bxavg, Byavg, Bzavg);
      cellParams[CellParams::BGBX   ] = Bxavg;
      cellParams[CellParams::BGBY   ] = Byavg;
      cellParams[CellParams::BGBZ   ] = Bzavg;
   }
}



