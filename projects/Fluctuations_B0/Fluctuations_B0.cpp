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

typedef fluctuationsParameters FlucP;
Real FlucP::BX0 = NAN;
Real FlucP::BY0 = NAN;
Real FlucP::BZ0 = NAN;
Real FlucP::DENSITY = NAN;
Real FlucP::TEMPERATURE = NAN;
Real FlucP::magPertAmp = NAN;
Real FlucP::densityPertAmp = NAN;
Real FlucP::velocityPertAmp = NAN;
Real FlucP::maxwCutoff = NAN;
uint FlucP::nSpaceSamples = 0;
uint FlucP::nVelocitySamples = 0;

bool initializeProject(void) { return true; }

bool addProjectParameters() {
   typedef Readparameters RP;
   RP::add("Fluctuations.BX0", "Background field value (T)", 1.0e-9);
   RP::add("Fluctuations.BY0", "Background field value (T)", 2.0e-9);
   RP::add("Fluctuations.BZ0", "Background field value (T)", 3.0e-9);
   RP::add("Fluctuations.rho", "Number density (m^-3)", 1.0e7);
   RP::add("Fluctuations.Temperature", "Temperature (K)", 2.0e6);
   RP::add("Fluctuations.magPertAmp", "Amplitude of the magnetic perturbation", 1.0e-9);
   RP::add("Fluctuations.densityPertAmp", "Amplitude factor of the density perturbation", 0.1);
   RP::add("Fluctuations.velocityPertAmp", "Amplitude of the velocity perturbation", 1.0e6);
   RP::add("Fluctuations.nSpaceSamples", "Number of sampling points per spatial dimension", 2);
   RP::add("Fluctuations.nVelocitySamples", "Number of sampling points per velocity dimension", 5);
   RP::add("Fluctuations.maxwCutoff", "Cutoff for the maxwellian distribution", 1e-12);
   return true;
}

bool getProjectParameters() {
   typedef Readparameters RP;
   RP::get("Fluctuations.BX0", FlucP::BX0);
   RP::get("Fluctuations.BY0", FlucP::BY0);
   RP::get("Fluctuations.BZ0", FlucP::BZ0);
   RP::get("Fluctuations.rho", FlucP::DENSITY);
   RP::get("Fluctuations.Temperature", FlucP::TEMPERATURE);
   RP::get("Fluctuations.magPertAmp", FlucP::magPertAmp);
   RP::get("Fluctuations.densityPertAmp", FlucP::densityPertAmp);
   RP::get("Fluctuations.velocityPertAmp", FlucP::velocityPertAmp);
   RP::get("Fluctuations.nSpaceSamples", FlucP::nSpaceSamples);
   RP::get("Fluctuations.nVelocitySamples", FlucP::nVelocitySamples);
   RP::get("Fluctuations.maxwCutoff", FlucP::maxwCutoff);
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

Real getDistribValue(creal& vx,creal& vy, creal& vz) {
   creal k = 1.3806505e-23; // Boltzmann
   creal mass = 1.67262171e-27; // m_p in kg
   return exp(- mass * (vx*vx + vy*vy + vz*vz) / (2.0 * k * FlucP::TEMPERATURE));
}

Real calcPhaseSpaceDensity(creal& x, creal& y, creal& z, creal& dx, creal& dy, creal& dz, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz) {
   if(vx < Parameters::vxmin + 0.5 * dvx ||
      vy < Parameters::vymin + 0.5 * dvy ||
      vz < Parameters::vzmin + 0.5 * dvz ||
      vx > Parameters::vxmax - 1.5 * dvx ||
      vy > Parameters::vymax - 1.5 * dvy ||
      vz > Parameters::vzmax - 1.5 * dvz
   ) return 0.0;
   
   creal mass = Parameters::m;
   creal q = Parameters::q;
   creal k = 1.3806505e-23; // Boltzmann
   creal mu0 = 1.25663706144e-6; // mu_0
   
   static int rndRho = 0;
   static int rndVel[3] = {0};
   #pragma omp threadprivate(rndRho,rndVel)
   int cellID = (int) (x / dx) +
                (int) (y / dy) * Parameters::xcells_ini +
                (int) (z / dz) * Parameters::xcells_ini * Parameters::ycells_ini;
#pragma omp critical
   {
      //critical region since srand not thread-safe. A nicer fix would be to use a thread-safe rng
      //e.g., http://linux.die.net/man/3/random_r (FIXME)
      srand(cellID);
      rndRho = rand();
      rndVel[0] =rand();
      rndVel[1] =rand();
      rndVel[2] =rand();
   } // end of omp critical region
   
   creal d_vx = dvx / (FlucP::nVelocitySamples-1);
   creal d_vy = dvy / (FlucP::nVelocitySamples-1);
   creal d_vz = dvz / (FlucP::nVelocitySamples-1);
   Real avg = 0.0;
   
   for (uint vi=0; vi<FlucP::nVelocitySamples; ++vi)
      for (uint vj=0; vj<FlucP::nVelocitySamples; ++vj)
	 for (uint vk=0; vk<FlucP::nVelocitySamples; ++vk)
         {
	    avg += getDistribValue(
	       vx+vi*d_vx - FlucP::velocityPertAmp * (0.5 - (double)rndVel[0] / (double)RAND_MAX),
	       vy+vj*d_vy - FlucP::velocityPertAmp * (0.5 - (double)rndVel[1] / (double)RAND_MAX),
	       vz+vk*d_vz - FlucP::velocityPertAmp * (0.5 - (double)rndVel[2] / (double)RAND_MAX));
         }
   
   creal result = avg *
	          FlucP::DENSITY * (1.0 + FlucP::densityPertAmp * (0.5 - (double)rndRho / (double)RAND_MAX)) *
	          pow(mass / (2.0 * M_PI * k * FlucP::TEMPERATURE), 1.5) /
//   	            (Parameters::vzmax - Parameters::vzmin) / 
                  (FlucP::nVelocitySamples*FlucP::nVelocitySamples*FlucP::nVelocitySamples);
   if(result < FlucP::maxwCutoff) {
      return 0.0;
   } else {
      return result;
   }
}
      
void calcBlockParameters(Real* blockParams) {}

void calcCellParameters(Real* cellParams,creal& t) {
   creal x = cellParams[CellParams::XCRD];
   creal dx = cellParams[CellParams::DX];
   creal y = cellParams[CellParams::YCRD];
   creal dy = cellParams[CellParams::DY];
   creal z = cellParams[CellParams::ZCRD];
   creal dz = cellParams[CellParams::DZ];
   
   int cellID = (int) (x / dx) +
   (int) (y / dy) * Parameters::xcells_ini +
   (int) (z / dz) * Parameters::xcells_ini * Parameters::ycells_ini;
   srand(cellID);

   // default initial values
   cellParams[CellParams::EX] = 0.0;
   cellParams[CellParams::EY] = 0.0;
   cellParams[CellParams::EZ] = 0.0;
   cellParams[CellParams::BX] = 0.0;
   cellParams[CellParams::BY] = 0.0;
   cellParams[CellParams::BZ] = 0.0;
   
   // user supplied initial values
   cellParams[CellParams::BX0] = FlucP::BX0;
   cellParams[CellParams::BY0] = FlucP::BY0;
   cellParams[CellParams::BZ0] = FlucP::BZ0;
#pragma omp critical
   {
      //critical region since srand not thread-safe. A nicer fix would be to use a thread-safe rng
      //e.g., http://linux.die.net/man/3/random_r (FIXME)
      cellParams[CellParams::BY   ] = FlucP::magPertAmp * (0.5 - (double)rand() / (double)RAND_MAX);
      cellParams[CellParams::BZ   ] = FlucP::magPertAmp * (0.5 - (double)rand() / (double)RAND_MAX);
   } // end of omp critical region
}

// TODO use this instead: template <class Grid, class CellData> void calcSimParameters(Grid<CellData>& mpiGrid...
void calcSimParameters(dccrg::Dccrg<SpatialCell>& mpiGrid, creal& t, Real& /*dt*/) {
   const std::vector<uint64_t> cells = mpiGrid.get_cells();
   for (uint i = 0; i < cells.size(); ++i) {
      calcCellParameters(mpiGrid[cells[i]]->parameters, t);
   }
}

