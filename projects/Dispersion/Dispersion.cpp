/*
 * This file is part of Vlasiator.
 * 
 * Copyright 2011, 2012 Finnish Meteorological Institute
 * 
 * Vlasiator is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 3
 * as published by the Free Software Foundation.
 * 
 * Vlasiator is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with Vlasiator. If not, see <http://www.gnu.org/licenses/>.
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

typedef dispersionParameters DispP;
Real DispP::B0 = NAN;
Real DispP::angleXY = NAN;
Real DispP::angleXZ = NAN;
Real DispP::DENSITY = NAN;
Real DispP::TEMPERATURE = NAN;
Real DispP::magXPertAbsAmp = NAN;
Real DispP::magYPertAbsAmp = NAN;
Real DispP::magZPertAbsAmp = NAN;
Real DispP::densityPertRelAmp = NAN;
Real DispP::velocityPertAbsAmp = NAN;
Real DispP::maxwCutoff = NAN;
uint DispP::nSpaceSamples = 0;
uint DispP::nVelocitySamples = 0;

bool initializeProject(void) {return true;}

bool addProjectParameters() {
   typedef Readparameters RP;
   RP::add("Dispersion.B0", "Guide magnetic field strength (T)", 1.0e-9);
   RP::add("Dispersion.angleXY", "Orientation of the guide magnetic field with respect to the x-axis in x-y plane (rad)", 0.001);
   RP::add("Dispersion.angleXZ", "Orientation of the guide magnetic field with respect to the x-axis in x-z plane (rad)", 0.001);
   RP::add("Dispersion.rho", "Number density (m^-3)", 1.0e7);
   RP::add("Dispersion.Temperature", "Temperature (K)", 2.0e6);
   RP::add("Dispersion.magXPertAbsAmp", "Absolute amplitude of the magnetic perturbation along x (T)", 1.0e-9);
   RP::add("Dispersion.magYPertAbsAmp", "Absolute amplitude of the magnetic perturbation along y (T)", 1.0e-9);
   RP::add("Dispersion.magZPertAbsAmp", "Absolute amplitude of the magnetic perturbation along z (T)", 1.0e-9);
   RP::add("Dispersion.densityPertRelAmp", "Relative amplitude of the density perturbation", 0.1);
   RP::add("Dispersion.velocityPertAbsAmp", "Absolute amplitude of the velocity perturbation", 1.0e6);
   RP::add("Dispersion.nSpaceSamples", "Number of sampling points per spatial dimension", 2);
   RP::add("Dispersion.nVelocitySamples", "Number of sampling points per velocity dimension", 5);
   RP::add("Dispersion.maxwCutoff", "Cutoff for the maxwellian distribution", 1e-12);
   return true;
}

bool getProjectParameters() {
   typedef Readparameters RP;
   RP::get("Dispersion.B0", DispP::B0);
   RP::get("Dispersion.angleXY", DispP::angleXY);
   RP::get("Dispersion.angleXZ", DispP::angleXZ);
   RP::get("Dispersion.rho", DispP::DENSITY);
   RP::get("Dispersion.Temperature", DispP::TEMPERATURE);
   RP::get("Dispersion.magXPertAbsAmp", DispP::magXPertAbsAmp);
   RP::get("Dispersion.magYPertAbsAmp", DispP::magYPertAbsAmp);
   RP::get("Dispersion.magZPertAbsAmp", DispP::magZPertAbsAmp);
   RP::get("Dispersion.densityPertRelAmp", DispP::densityPertRelAmp);
   RP::get("Dispersion.velocityPertAbsAmp", DispP::velocityPertAbsAmp);
   RP::get("Dispersion.nSpaceSamples", DispP::nSpaceSamples);
   RP::get("Dispersion.nVelocitySamples", DispP::nVelocitySamples);
   RP::get("Dispersion.maxwCutoff", DispP::maxwCutoff);
   return true;
}

bool cellParametersChanged(creal& t) {return false;}

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
   initstate_r(cellID, &rngStateBuffer[0], 256, &rngDataBuffer);
   
   int32_t rndRho, rndVel[3];
   random_r(&rngDataBuffer, &rndRho);
   random_r(&rngDataBuffer, &rndVel[0]);
   random_r(&rngDataBuffer, &rndVel[1]);
   random_r(&rngDataBuffer, &rndVel[2]);
   
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
   return exp(- mass * (vx*vx + vy*vy + vz*vz) / (2.0 * k * DispP::TEMPERATURE));
}

/* calcPhaseSpaceDensity needs to be thread-safe */
Real calcPhaseSpaceDensity(creal& x, creal& y, creal& z, creal& dx, creal& dy, creal& dz, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz, const int32_t& rndRho, const int32_t rndVel[3]) {
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
   
   creal d_vx = dvx / (DispP::nVelocitySamples-1);
   creal d_vy = dvy / (DispP::nVelocitySamples-1);
   creal d_vz = dvz / (DispP::nVelocitySamples-1);
   Real avg = 0.0;
   
   for (uint vi=0; vi<DispP::nVelocitySamples; ++vi)
      for (uint vj=0; vj<DispP::nVelocitySamples; ++vj)
         for (uint vk=0; vk<DispP::nVelocitySamples; ++vk)
         {
            avg += getDistribValue(
               vx+vi*d_vx - DispP::velocityPertAbsAmp * (0.5 - (double)rndVel[0] / (double)RAND_MAX),
                                   vy+vj*d_vy - DispP::velocityPertAbsAmp * (0.5 - (double)rndVel[1] / (double)RAND_MAX),
                                   vz+vk*d_vz - DispP::velocityPertAbsAmp * (0.5 - (double)rndVel[2] / (double)RAND_MAX));
         }
         
         creal result = avg *
         DispP::DENSITY * (1.0 + DispP::densityPertRelAmp * (0.5 - (double)rndRho / (double)RAND_MAX)) *
         pow(mass / (2.0 * M_PI * k * DispP::TEMPERATURE), 1.5) /
         //            (Parameters::vzmax - Parameters::vzmin) / 
         (DispP::nVelocitySamples*DispP::nVelocitySamples*DispP::nVelocitySamples);
         if(result < DispP::maxwCutoff) {
            return 0.0;
         } else {
            return result;
         }
}

void calcBlockParameters(Real* blockParams) {
   //blockParams[BlockParams::Q_PER_M] = 1.0;
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
   initstate_r(cellID + 
   (uint)(Parameters::xcells_ini*Parameters::ycells_ini*Parameters::zcells_ini), &rngStateBuffer[0], 256, &rngDataBuffer);
   
   int32_t rndBuffer[3];
   random_r(&rngDataBuffer, &rndBuffer[0]);
   random_r(&rngDataBuffer, &rndBuffer[1]);
   random_r(&rngDataBuffer, &rndBuffer[2]);

   cellParams[CellParams::PERBX] = 0.0;
   cellParams[CellParams::PERBY] = 0.0;
   cellParams[CellParams::PERBZ] = 0.0;
   
   cellParams[CellParams::BGBX] = DispP::B0 * cos(DispP::angleXY) * cos(DispP::angleXZ) +
      DispP::magXPertAbsAmp * (0.5 - (double)rndBuffer[0] / (double)RAND_MAX);
   cellParams[CellParams::BGBY] = DispP::B0 * sin(DispP::angleXY) * cos(DispP::angleXZ) + 
      DispP::magYPertAbsAmp * (0.5 - (double)rndBuffer[1] / (double)RAND_MAX);
   cellParams[CellParams::BGBZ] = DispP::B0 * sin(DispP::angleXZ) +
      DispP::magZPertAbsAmp * (0.5 - (double)rndBuffer[2] / (double)RAND_MAX);
}

// TODO use this instead: template <class Grid, class CellData> void calcSimParameters(Grid<CellData>& mpiGrid...
void calcSimParameters(dccrg::Dccrg<SpatialCell>& mpiGrid, creal& t, Real& /*dt*/) {
   const std::vector<uint64_t> cells = mpiGrid.get_cells();
   for (uint i = 0; i < cells.size(); ++i) {
      calcCellParameters(mpiGrid[cells[i]]->parameters, t);
   }
}

