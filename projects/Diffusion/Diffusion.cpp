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
#include <iomanip>
#include <cmath>

#include "spatial_cell.hpp"
#include "common.h"
#include "project.h"
#include "parameters.h"
#include "readparameters.h"
#include "vlasovmover.h"

typedef diffusionParameters DiffP;
Real DiffP::B0 = NAN;
Real DiffP::DENSITY = NAN;
Real DiffP::SCA_X = NAN;
Real DiffP::SCA_Y = NAN;
Real DiffP::TEMPERATURE = NAN;
uint DiffP::nSpaceSamples = 0;
uint DiffP::nVelocitySamples = 0;

bool initializeProject(void) {return true;}

bool addProjectParameters(){
   typedef Readparameters RP;
   RP::add("Diffusion.B0", "Background field value (T)", 1.0e-9);
   RP::add("Diffusion.rho", "Number density (m^-3)", 1.0e7);
   RP::add("Diffusion.Temperature", "Temperature (K)", 2.0e6);
   RP::add("Diffusion.Scale_x", "Scale length in x (m)", 100000.0);
   RP::add("Diffusion.Scale_y", "Scale length in y (m)", 100000.0);
   RP::add("Diffusion.nSpaceSamples", "Number of sampling points per spatial dimension", 2);
   RP::add("Diffusion.nVelocitySamples", "Number of sampling points per velocity dimension", 5);
   
   return true;
}

bool getProjectParameters(){
   typedef Readparameters RP;
   RP::get("Diffusion.B0", DiffP::B0);
   RP::get("Diffusion.rho", DiffP::DENSITY);
   RP::get("Diffusion.Temperature", DiffP::TEMPERATURE);
   RP::get("Diffusion.Scale_x", DiffP::SCA_X);
   RP::get("Diffusion.Scale_y", DiffP::SCA_Y);
   RP::get("Diffusion.nSpaceSamples", DiffP::nSpaceSamples);
   RP::get("Diffusion.nVelocitySamples", DiffP::nVelocitySamples);
   
   return true;
}


bool cellParametersChanged(creal& t) {return false;}

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

Real getDistribValue(creal& x,creal& y,creal& z,
		     creal& vx,creal& vy,creal& vz) {
  
  creal mass = 1.67262171e-27; // m_p in kg
  creal k = 1.3806505e-23; // Boltzmann
  creal mu0 = 1.25663706144e-6; // mu_0
  creal q = 1.60217653e-19; // q_i
  
  return DiffP::DENSITY * pow(mass / (2.0 * M_PI * k * DiffP::TEMPERATURE), 1.5) * (
     5.0 * exp(- (pow(x, 2.0) / pow(DiffP::SCA_X, 2.0) +  pow(y, 2.0) / pow(DiffP::SCA_Y, 2.0))) * 
     exp(- mass * (pow(vx, 2.0) + pow(vy, 2.0) + pow(vz, 2.0)) / (2.0 * k * DiffP::TEMPERATURE))
    +
    exp(- mass * (pow(vx, 2.0) + pow(vy, 2.0) + pow(vz, 2.0)) / (2.0 * k * DiffP::TEMPERATURE)));
}

Real calcPhaseSpaceDensity(creal& x, creal& y, creal& z, creal& dx, creal& dy, creal& dz, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz) {
   creal d_x = dx / (DiffP::nSpaceSamples-1);
   creal d_y = dy / (DiffP::nSpaceSamples-1);
   creal d_z = dz / (DiffP::nSpaceSamples-1);
   creal d_vx = dvx / (DiffP::nVelocitySamples-1);
   creal d_vy = dvy / (DiffP::nVelocitySamples-1);
   creal d_vz = dvz / (DiffP::nVelocitySamples-1);
   Real avg = 0.0;
   for (uint i=0; i<DiffP::nSpaceSamples; ++i)
      for (uint j=0; j<DiffP::nSpaceSamples; ++j)
	 for (uint k=0; k<DiffP::nSpaceSamples; ++k)
	    for (uint vi=0; vi<DiffP::nVelocitySamples; ++vi)
	       for (uint vj=0; vj<DiffP::nVelocitySamples; ++vj)
		  for (uint vk=0; vk<DiffP::nVelocitySamples; ++vk)
		  {
		     avg += getDistribValue(x+i*d_x, y+j*d_y, z+k*d_z, vx+vi*d_vx, vy+vj*d_vy, vz+vk*d_vz);
		  }
   return avg / (DiffP::nSpaceSamples*DiffP::nSpaceSamples*DiffP::nSpaceSamples) / (DiffP::nVelocitySamples*DiffP::nVelocitySamples*DiffP::nVelocitySamples);
	       
}

void calcBlockParameters(Real* blockParams) {
   //blockParams[BlockParams::Q_PER_M] = 1.0;
}

void calcCellParameters(Real* cellParams,creal& t) {
  creal x = cellParams[CellParams::XCRD];
  creal dx = cellParams[CellParams::DX];
  
  cellParams[CellParams::EX   ] = 0.0;
  cellParams[CellParams::EY   ] = 0.0;
  cellParams[CellParams::EZ   ] = 0.0;
  cellParams[CellParams::PERBX   ] = 0.0;
  cellParams[CellParams::PERBY   ] = 0.0;
  cellParams[CellParams::PERBZ   ] = 0.0;
  cellParams[CellParams::BGBX   ] = 0.0;
  cellParams[CellParams::BGBY   ] = 0.0;
  cellParams[CellParams::BGBZ   ] = DiffP::B0;
}

// TODO use this instead: template <class Grid, class CellData> void calcSimParameters(Grid<CellData>& mpiGrid...
void calcSimParameters(dccrg::Dccrg<SpatialCell>& mpiGrid, creal& t, Real& /*dt*/) {
   std::vector<uint64_t> cells = mpiGrid.get_cells();
   for (uint i = 0; i < cells.size(); ++i) {
      calcCellParameters(mpiGrid[cells[i]]->parameters, t);
   }
}

