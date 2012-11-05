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


// TODO when projects are classes, the rnd values and state variables etc can be class members.
/*! Integrate the distribution function over the given six-dimensional phase-space cell.
 * \param x Starting value of the x-coordinate of the cell.
 * \param y Starting value of the y-coordinate of the cell.
 * \param z Starting value of the z-coordinate of the cell.
 * \param dx The size of the cell in x-direction.
 * \param dy The size of the cell in y-direction.
 * \param dz The size of the cell in z-direction.
 * \param vx Starting value of the vx-coordinate of the cell.
 * \param vy Starting value of the vy-coordinate of the cell.
 * \param vz Starting value of the vz-coordinate of the cell.
 * \param dvx The size of the cell in vx-direction.
 * \param dvy The size of the cell in vy-direction.
 * \param dvz The size of the cell in vz-direction.
 * \param rndRho Random number for the density perturbation.
 * \param rndVel Three random numbers for the velocity perturbation.
 * \return The volume average of the distribution function in the given phase space cell.
 * The physical unit of this quantity is 1 / (m^3 (m/s)^3).
 */
Real calcPhaseSpaceDensity(creal& x,creal& y,creal& z,
                           creal& dx,creal& dy,creal& dz,
                           creal& vx,creal& vy,creal& vz,
                           creal& dvx,creal& dvy,creal& dvz,
                           const int32_t& rndRho, const int32_t rndVel[3]);


typedef fluctuationsParameters FlucP;
Real FlucP::BX0 = NAN;
Real FlucP::BY0 = NAN;
Real FlucP::BZ0 = NAN;
Real FlucP::DENSITY = NAN;
Real FlucP::TEMPERATURE = NAN;
Real FlucP::magXPertAbsAmp = NAN;
Real FlucP::magYPertAbsAmp = NAN;
Real FlucP::magZPertAbsAmp = NAN;
Real FlucP::densityPertRelAmp = NAN;
Real FlucP::velocityPertAbsAmp = NAN;
Real FlucP::maxwCutoff = NAN;
uint FlucP::nSpaceSamples = 0;
uint FlucP::nVelocitySamples = 0;

bool initializeProject(void) {return true;}

bool addProjectParameters() {
   typedef Readparameters RP;
   RP::add("Fluctuations.BX0", "Background field value (T)", 1.0e-9);
   RP::add("Fluctuations.BY0", "Background field value (T)", 2.0e-9);
   RP::add("Fluctuations.BZ0", "Background field value (T)", 3.0e-9);
   RP::add("Fluctuations.rho", "Number density (m^-3)", 1.0e7);
   RP::add("Fluctuations.Temperature", "Temperature (K)", 2.0e6);
   RP::add("Fluctuations.magXPertAbsAmp", "Amplitude of the magnetic perturbation along x", 1.0e-9);
   RP::add("Fluctuations.magYPertAbsAmp", "Amplitude of the magnetic perturbation along y", 1.0e-9);
   RP::add("Fluctuations.magZPertAbsAmp", "Amplitude of the magnetic perturbation along z", 1.0e-9);
   RP::add("Fluctuations.densityPertRelAmp", "Amplitude factor of the density perturbation", 0.1);
   RP::add("Fluctuations.velocityPertAbsAmp", "Amplitude of the velocity perturbation", 1.0e6);
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
   RP::get("Fluctuations.magXPertAbsAmp", FlucP::magXPertAbsAmp);
   RP::get("Fluctuations.magYPertAbsAmp", FlucP::magYPertAbsAmp);
   RP::get("Fluctuations.magZPertAbsAmp", FlucP::magZPertAbsAmp);
   RP::get("Fluctuations.densityPertRelAmp", FlucP::densityPertRelAmp);
   RP::get("Fluctuations.velocityPertAbsAmp", FlucP::velocityPertAbsAmp);
   RP::get("Fluctuations.nSpaceSamples", FlucP::nSpaceSamples);
   RP::get("Fluctuations.nVelocitySamples", FlucP::nVelocitySamples);
   RP::get("Fluctuations.maxwCutoff", FlucP::maxwCutoff);
   return true;
}


void setProjectCell(SpatialCell* cell) {
   // Set up cell parameters:
   calcCellParameters(&((*cell).parameters[0]), 0.0);
   
   cell->parameters[CellParams::RHOLOSSADJUST] = 0.0;
   cell->parameters[CellParams::RHOLOSSVELBOUNDARY] = 0.0;
   
   creal x = cell->parameters[CellParams::XCRD];
   creal y = cell->parameters[CellParams::YCRD];
   creal z = cell->parameters[CellParams::ZCRD];
   creal dx = cell->parameters[CellParams::DX];
   creal dy = cell->parameters[CellParams::DY];
   creal dz = cell->parameters[CellParams::DZ];
   
   uint cellID = (int) ((x - Parameters::xmin) / dx) +
   (int) ((y - Parameters::ymin) / dy) * Parameters::xcells_ini +
   (int) ((z - Parameters::zmin) / dz) * Parameters::xcells_ini * Parameters::ycells_ini;
   
   // TODO when projects are classes, the rnd values and state variables etc can be class members.
   char rngStateBuffer[256];
   random_data rngDataBuffer;
   memset(&rngDataBuffer, 0, sizeof(rngDataBuffer));

   #ifndef _AIX
   initstate_r(cellID, &rngStateBuffer[0], 256, &rngDataBuffer);
   int32_t rndRho, rndVel[3];
   random_r(&rngDataBuffer, &rndRho);
   random_r(&rngDataBuffer, &rndVel[0]);
   random_r(&rngDataBuffer, &rndVel[1]);
   random_r(&rngDataBuffer, &rndVel[2]);
   #else
   initstate_r(cellID, &rngStateBuffer[0], 256, NULL, &rngDataBuffer);
   int64_t rndRho, rndVel[3];
   random_r(&rndRho, &rngDataBuffer);
   random_r(&rndVel[0], &rngDataBuffer);
   random_r(&rndVel[1], &rngDataBuffer);
   random_r(&rndVel[2], &rngDataBuffer);
   #endif
   
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
                     calcPhaseSpaceDensity(x, y, z, dx, dy, dz,
                                           vx_cell, vy_cell, vz_cell,
                                           dvx_blockCell, dvy_blockCell, dvz_blockCell,
                                           rndRho,
                                           rndVel);
                     
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

/* calcPhaseSpaceDensity needs to be thread-safe */
Real calcPhaseSpaceDensity(
   creal& x, creal& y, creal& z,
   creal& dx, creal& dy, creal& dz,
   creal& vx, creal& vy, creal& vz,
   creal& dvx, creal& dvy, creal& dvz,
   #ifndef _AIX
   const int32_t& rndRho,
   const int32_t rndVel[3]
   #else
   const int64_t& rndRho,
   const int64_t rndVel[3]
   #endif
) {
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
   
   creal d_vx = dvx / (FlucP::nVelocitySamples-1);
   creal d_vy = dvy / (FlucP::nVelocitySamples-1);
   creal d_vz = dvz / (FlucP::nVelocitySamples-1);
   Real avg = 0.0;
   
   for (uint vi=0; vi<FlucP::nVelocitySamples; ++vi)
      for (uint vj=0; vj<FlucP::nVelocitySamples; ++vj)
         for (uint vk=0; vk<FlucP::nVelocitySamples; ++vk)
         {
            avg += getDistribValue(
               vx+vi*d_vx - FlucP::velocityPertAbsAmp * (0.5 - (double)rndVel[0] / (double)RAND_MAX),
               vy+vj*d_vy - FlucP::velocityPertAbsAmp * (0.5 - (double)rndVel[1] / (double)RAND_MAX),
               vz+vk*d_vz - FlucP::velocityPertAbsAmp * (0.5 - (double)rndVel[2] / (double)RAND_MAX));
         }
   
   creal result = avg *
            FlucP::DENSITY * (1.0 + FlucP::densityPertRelAmp * (0.5 - (double)rndRho / (double)RAND_MAX)) *
            pow(mass / (2.0 * M_PI * k * FlucP::TEMPERATURE), 1.5) /
//            (Parameters::vzmax - Parameters::vzmin) / 
           (FlucP::nVelocitySamples*FlucP::nVelocitySamples*FlucP::nVelocitySamples);
   if(result < FlucP::maxwCutoff) {
      return 0.0;
   } else {
      return result;
   }
}
      


void calcCellParameters(Real* cellParams,creal& t) {
   creal x = cellParams[CellParams::XCRD];
   creal dx = cellParams[CellParams::DX];
   creal y = cellParams[CellParams::YCRD];
   creal dy = cellParams[CellParams::DY];
   creal z = cellParams[CellParams::ZCRD];
   creal dz = cellParams[CellParams::DZ];
   
   uint cellID = (int) ((x - Parameters::xmin) / dx) +
   (int) ((y - Parameters::ymin) / dy) * Parameters::xcells_ini +
   (int) ((z - Parameters::zmin) / dz) * Parameters::xcells_ini * Parameters::ycells_ini;
   
   cellParams[CellParams::EX   ] = 0.0;
   cellParams[CellParams::EY   ] = 0.0;
   cellParams[CellParams::EZ   ] = 0.0;
   
   // TODO when projects are classes, the rnd values and state variables etc can be class members.
   char rngStateBuffer[256];
   random_data rngDataBuffer;
   memset(&rngDataBuffer, 0, sizeof(rngDataBuffer));

   #ifndef _AIX
   initstate_r(
      cellID + (uint)(Parameters::xcells_ini*Parameters::ycells_ini*Parameters::zcells_ini),
      &rngStateBuffer[0],
      256,
      &rngDataBuffer
   );
   int32_t rndBuffer[3];
   random_r(&rngDataBuffer, &rndBuffer[0]);
   random_r(&rngDataBuffer, &rndBuffer[1]);
   random_r(&rngDataBuffer, &rndBuffer[2]);
   #else
   initstate_r(
      cellID + (uint)(Parameters::xcells_ini*Parameters::ycells_ini*Parameters::zcells_ini),
      &rngStateBuffer[0],
      256,
      NULL,
      &rngDataBuffer
   );
   int64_t rndBuffer[3];
   random_r(&rndBuffer[0], &rngDataBuffer);
   random_r(&rndBuffer[1], &rngDataBuffer);
   random_r(&rndBuffer[2], &rngDataBuffer);
   #endif
   
   cellParams[CellParams::BGBX]  = FlucP::BX0 ;
   cellParams[CellParams::PERBX] = FlucP::magXPertAbsAmp * (0.5 - (double)rndBuffer[0] / (double)RAND_MAX);
   cellParams[CellParams::BGBY]   = FlucP::BY0; 
   cellParams[CellParams::PERBY] = FlucP::magYPertAbsAmp * (0.5 - (double)rndBuffer[1] / (double)RAND_MAX);
   cellParams[CellParams::BGBZ] = FlucP::BZ0;
   cellParams[CellParams::PERBZ] = FlucP::magZPertAbsAmp * (0.5 - (double)rndBuffer[2] / (double)RAND_MAX);
}



