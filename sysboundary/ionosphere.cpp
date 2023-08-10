/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
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

/*!\file ionosphere.cpp
 * \brief Implementation of the class SysBoundaryCondition::Ionosphere to handle cells classified as sysboundarytype::IONOSPHERE.
 */

#include <cstdlib>
#include <iostream>

#include "ionosphere.h"
#include "../projects/project.h"
#include "../projects/projects_common.h"
#include "../vlasovmover.h"
#include "../fieldsolver/fs_common.h"
#include "../fieldsolver/fs_limiters.h"
#include "../fieldsolver/ldz_magnetic_field.hpp"
#include "../common.h"
#include "../object_wrapper.h"

#ifndef NDEBUG
   #define DEBUG_IONOSPHERE
#endif
#ifdef DEBUG_SYSBOUNDARY
   #define DEBUG_IONOSPHERE
#endif

extern ARCH_MANAGED GridParameters meshParams;

namespace SBC {
   Ionosphere::Ionosphere(): SysBoundaryCondition() { }
   
   Ionosphere::~Ionosphere() { }
   
   void Ionosphere::addParameters() {
      Readparameters::add("ionosphere.centerX", "X coordinate of ionosphere center (m)", 0.0);
      Readparameters::add("ionosphere.centerY", "Y coordinate of ionosphere center (m)", 0.0);
      Readparameters::add("ionosphere.centerZ", "Z coordinate of ionosphere center (m)", 0.0);
      Readparameters::add("ionosphere.radius", "Radius of ionosphere (m).", 1.0e7);
      Readparameters::add("ionosphere.geometry", "Select the geometry of the ionosphere, 0: inf-norm (diamond), 1: 1-norm (square), 2: 2-norm (circle, DEFAULT), 3: 2-norm cylinder aligned with y-axis, use with polar plane/line dipole.", 2);
      Readparameters::add("ionosphere.precedence", "Precedence value of the ionosphere system boundary condition (integer), the higher the stronger.", 2);
      Readparameters::add("ionosphere.reapplyUponRestart", "If 0 (default), keep going with the state existing in the restart file. If 1, calls again applyInitialState. Can be used to change boundary condition behaviour during a run.", 0);

      // Per-population parameters
      for(uint i=0; i< getObjectWrapper().particleSpecies.size(); i++) {
         const std::string& pop = getObjectWrapper().particleSpecies[i].name;
         
         Readparameters::add(pop + "_ionosphere.rho", "Number density of the ionosphere (m^-3)", 0.0);
         Readparameters::add(pop + "_ionosphere.T", "Temperature of the ionosphere (K)", 0.0);
         Readparameters::add(pop + "_ionosphere.VX0", "Bulk velocity of ionospheric distribution function in X direction (m/s)", 0.0);
         Readparameters::add(pop + "_ionosphere.VY0", "Bulk velocity of ionospheric distribution function in X direction (m/s)", 0.0);
         Readparameters::add(pop + "_ionosphere.VZ0", "Bulk velocity of ionospheric distribution function in X direction (m/s)", 0.0);
         Readparameters::add(pop + "_ionosphere.fluffiness", "Inertia of boundary smoothing when copying neighbour's moments and velocity distributions (0=completely constant boundaries, 1=neighbours are interpolated immediately).", 0);
      }
   }
   
   void Ionosphere::getParameters() {
   
      Readparameters::get("ionosphere.centerX", this->center[0]);
      Readparameters::get("ionosphere.centerY", this->center[1]);
      Readparameters::get("ionosphere.centerZ", this->center[2]);
      Readparameters::get("ionosphere.radius", this->radius);
      Readparameters::get("ionosphere.geometry", this->geometry);
      Readparameters::get("ionosphere.precedence", this->precedence);
      uint reapply;
      Readparameters::get("ionosphere.reapplyUponRestart", reapply);
      this->applyUponRestart = false;
      if(reapply == 1) {
         this->applyUponRestart = true;
      }

      for(uint i=0; i< getObjectWrapper().particleSpecies.size(); i++) {
        const std::string& pop = getObjectWrapper().particleSpecies[i].name;
        IonosphereSpeciesParameters sP;

        Readparameters::get(pop + "_ionosphere.rho", sP.rho);
        Readparameters::get(pop + "_ionosphere.VX0", sP.V0[0]);
        Readparameters::get(pop + "_ionosphere.VY0", sP.V0[1]);
        Readparameters::get(pop + "_ionosphere.VZ0", sP.V0[2]);
        Readparameters::get(pop + "_ionosphere.fluffiness", sP.fluffiness);
        Readparameters::get(pop + "_ionosphere.T", sP.T);
        Readparameters::get(pop + "_Magnetosphere.nSpaceSamples", sP.nSpaceSamples);
        Readparameters::get(pop + "_Magnetosphere.nVelocitySamples", sP.nVelocitySamples);

        speciesParams.push_back(sP);
      }
   }

   bool Ionosphere::initFieldBoundary() {
      fieldBoundary = new IonosphereFieldBoundary(center, radius, geometry);
      return true;
   }
   
   bool Ionosphere::initSysBoundary(
      creal& t,
      Project &project
   ) {
      getParameters();
      isThisDynamic = false;

      // iniSysBoundary is only called once, generateTemplateCell must 
      // init all particle species
      generateTemplateCell(project);
      
      return true;
   }

   Real getR(creal x,creal y,creal z, uint geometry, Real center[3]) {

      Real r;
      
      switch(geometry) {
      case 0:
         // infinity-norm, result is a diamond/square with diagonals aligned on the axes in 2D
         r = fabs(x-center[0]) + fabs(y-center[1]) + fabs(z-center[2]);
         break;
      case 1:
         // 1-norm, result is is a grid-aligned square in 2D
         r = max(max(fabs(x-center[0]), fabs(y-center[1])), fabs(z-center[2]));
         break;
      case 2:
         // 2-norm (Cartesian), result is a circle in 2D
         r = sqrt((x-center[0])*(x-center[0]) + (y-center[1])*(y-center[1]) + (z-center[2])*(z-center[2]));
         break;
      case 3:
         // 2-norm (Cartesian) cylinder aligned on y-axis
         r = sqrt((x-center[0])*(x-center[0]) + (z-center[2])*(z-center[2]));
         break;
      default:
         std::cerr << __FILE__ << ":" << __LINE__ << ":" << "ionosphere.geometry has to be 0, 1 or 2." << std::endl;
         abort();
      }

      return r;
   }
   
   bool Ionosphere::assignSysBoundary(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                      FsGrid< fsgrids::technical, 1, FS_STENCIL_WIDTH> & technicalGrid) {
      const vector<CellID>& cells = getLocalCells();
      for(uint i=0; i<cells.size(); i++) {
         if(mpiGrid[cells[i]]->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
            continue;
         }
         
         creal* const cellParams = &(mpiGrid[cells[i]]->parameters[0]);
         creal dx = cellParams[CellParams::DX];
         creal dy = cellParams[CellParams::DY];
         creal dz = cellParams[CellParams::DZ];
         creal x = cellParams[CellParams::XCRD] + 0.5*dx;
         creal y = cellParams[CellParams::YCRD] + 0.5*dy;
         creal z = cellParams[CellParams::ZCRD] + 0.5*dz;
         
         if(getR(x,y,z,this->geometry,this->center) < this->radius) {
            mpiGrid[cells[i]]->sysBoundaryFlag = this->getIndex();
         }
      }

      // Assign boundary flags to local fsgrid cells
      auto gridDims(technicalGrid.getLocalSize());  
      for (int k=0; k<gridDims[2]; k++) {
         for (int j=0; j<gridDims[1]; j++) {
            for (int i=0; i<gridDims[0]; i++) {
               const auto& coords = technicalGrid.getPhysicalCoords(i,j,k);
               
               // Shift to the center of the fsgrid cell
               auto cellCenterCoords = coords;
               cellCenterCoords[0] += 0.5 * technicalGrid.DX;
               cellCenterCoords[1] += 0.5 * technicalGrid.DY;
               cellCenterCoords[2] += 0.5 * technicalGrid.DZ;

               if(getR(cellCenterCoords[0],cellCenterCoords[1],cellCenterCoords[2],this->geometry,this->center) < this->radius) {
                  technicalGrid.get(i,j,k,0).sysBoundaryFlag = this->getIndex();
               }

            }
         }
      }

      return true;
   }

   bool Ionosphere::applyInitialState(
      const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      FsGrid<Real, fsgrids::bfield::N_BFIELD, FS_STENCIL_WIDTH> & perBGrid,
      Project &project
   ) {
      const vector<CellID>& cells = getLocalCells();
      #pragma omp parallel for
      for (uint i=0; i<cells.size(); ++i) {
         SpatialCell* cell = mpiGrid[cells[i]];
         if (cell->sysBoundaryFlag != this->getIndex()) continue;
         
         for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID)
            setCellFromTemplate(cell,popID);
      }
      return true;
   }

   void Ionosphere::vlasovBoundaryCondition(
      const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const CellID& cellID,
      const uint popID,
      const bool calculate_V_moments
   ) {
      phiprof::start("vlasovBoundaryCondition (Ionosphere)");
      this->vlasovBoundaryFluffyCopyFromAllCloseNbrs(mpiGrid, cellID, popID, calculate_V_moments, this->speciesParams[popID].fluffiness);
      phiprof::stop("vlasovBoundaryCondition (Ionosphere)");
   }

   /**
    * NOTE: This function must initialize all particle species!
    * @param project
    */
   void Ionosphere::generateTemplateCell(Project &project) {
      // WARNING not 0.0 here or the dipole() function fails miserably.
      templateCell.sysBoundaryFlag = this->getIndex();
      templateCell.sysBoundaryLayer = 1;
      templateCell.parameters[CellParams::XCRD] = 1.0;
      templateCell.parameters[CellParams::YCRD] = 1.0;
      templateCell.parameters[CellParams::ZCRD] = 1.0;
      templateCell.parameters[CellParams::DX] = 1;
      templateCell.parameters[CellParams::DY] = 1;
      templateCell.parameters[CellParams::DZ] = 1;
      
      // Loop over particle species
      for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
         const IonosphereSpeciesParameters& sP = this->speciesParams[popID];
         const vector<vmesh::GlobalID> blocksToInitialize = findBlocksToInitialize(templateCell,popID);
         Realf* data = templateCell.get_data(popID);
         
         for (size_t i = 0; i < blocksToInitialize.size(); i++) {
            const vmesh::GlobalID blockGID = blocksToInitialize.at(i);
            const vmesh::LocalID blockLID = templateCell.get_velocity_block_local_id(blockGID,popID);
            const Real* block_parameters = templateCell.get_block_parameters(blockLID,popID);
            creal vxBlock = block_parameters[BlockParams::VXCRD];
            creal vyBlock = block_parameters[BlockParams::VYCRD];
            creal vzBlock = block_parameters[BlockParams::VZCRD];
            creal dvxCell = block_parameters[BlockParams::DVX];
            creal dvyCell = block_parameters[BlockParams::DVY];
            creal dvzCell = block_parameters[BlockParams::DVZ];

            creal x = templateCell.parameters[CellParams::XCRD];
            creal y = templateCell.parameters[CellParams::YCRD];
            creal z = templateCell.parameters[CellParams::ZCRD];
            creal dx = templateCell.parameters[CellParams::DX];
            creal dy = templateCell.parameters[CellParams::DY];
            creal dz = templateCell.parameters[CellParams::DZ];
         
            // Calculate volume average of distrib. function for each cell in the block.
            for (uint kc=0; kc<WID; ++kc) for (uint jc=0; jc<WID; ++jc) for (uint ic=0; ic<WID; ++ic) {
               creal vxCell = vxBlock + ic*dvxCell;
               creal vyCell = vyBlock + jc*dvyCell;
               creal vzCell = vzBlock + kc*dvzCell;
               Real average = 0.0;
               if(sP.nVelocitySamples > 1) {
                  creal d_vx = dvxCell / (sP.nVelocitySamples-1);
                  creal d_vy = dvyCell / (sP.nVelocitySamples-1);
                  creal d_vz = dvzCell / (sP.nVelocitySamples-1);
                  for (uint vi=0; vi<sP.nVelocitySamples; ++vi)
                     for (uint vj=0; vj<sP.nVelocitySamples; ++vj)
                        for (uint vk=0; vk<sP.nVelocitySamples; ++vk) {
                           average +=  shiftedMaxwellianDistribution(
                                                                     popID,
                                                                     vxCell + vi*d_vx,
                                                                     vyCell + vj*d_vy,
                                                                     vzCell + vk*d_vz
                                                                    );
                        }
                  average /= sP.nVelocitySamples * sP.nVelocitySamples * sP.nVelocitySamples;
               } else {
                  average = shiftedMaxwellianDistribution(
                                                          popID,
                                                          vxCell + 0.5*dvxCell,
                                                          vyCell + 0.5*dvyCell,
                                                          vzCell + 0.5*dvzCell
                                                         );
               }

               if (average !=0.0 ) {
                  data[blockLID*WID3+cellIndex(ic,jc,kc)] = average;
               }
            } // for-loop over cells in velocity block
         } // for-loop over velocity blocks

         // let's get rid of blocks not fulfilling the criteria here to save memory.
         templateCell.adjustSingleCellVelocityBlocks(popID);
      } // for-loop over particle species
      
      calculateCellMoments(&templateCell,true,true);

      // WARNING Time-independence assumed here. Normal moments computed in setProjectCell
      templateCell.parameters[CellParams::RHOM_R] = templateCell.parameters[CellParams::RHOM];
      templateCell.parameters[CellParams::VX_R] = templateCell.parameters[CellParams::VX];
      templateCell.parameters[CellParams::VY_R] = templateCell.parameters[CellParams::VY];
      templateCell.parameters[CellParams::VZ_R] = templateCell.parameters[CellParams::VZ];
      templateCell.parameters[CellParams::RHOQ_R] = templateCell.parameters[CellParams::RHOQ];
      templateCell.parameters[CellParams::P_11_R] = templateCell.parameters[CellParams::P_11];
      templateCell.parameters[CellParams::P_22_R] = templateCell.parameters[CellParams::P_22];
      templateCell.parameters[CellParams::P_33_R] = templateCell.parameters[CellParams::P_33];
      templateCell.parameters[CellParams::RHOM_V] = templateCell.parameters[CellParams::RHOM];
      templateCell.parameters[CellParams::VX_V] = templateCell.parameters[CellParams::VX];
      templateCell.parameters[CellParams::VY_V] = templateCell.parameters[CellParams::VY];
      templateCell.parameters[CellParams::VZ_V] = templateCell.parameters[CellParams::VZ];
      templateCell.parameters[CellParams::RHOQ_V] = templateCell.parameters[CellParams::RHOQ];
      templateCell.parameters[CellParams::P_11_V] = templateCell.parameters[CellParams::P_11];
      templateCell.parameters[CellParams::P_22_V] = templateCell.parameters[CellParams::P_22];
      templateCell.parameters[CellParams::P_33_V] = templateCell.parameters[CellParams::P_33];
   }
   
   Real Ionosphere::shiftedMaxwellianDistribution(
      const uint popID,
      creal& vx, creal& vy, creal& vz
   ) {
      
      const Real MASS = getObjectWrapper().particleSpecies[popID].mass;
      const IonosphereSpeciesParameters& sP = this->speciesParams[popID];

      return sP.rho * pow(MASS /
      (2.0 * M_PI * physicalconstants::K_B * sP.T), 1.5) *
      exp(-MASS * ((vx-sP.V0[0])*(vx-sP.V0[0]) + (vy-sP.V0[1])*(vy-sP.V0[1]) + (vz-sP.V0[2])*(vz-sP.V0[2])) /
      (2.0 * physicalconstants::K_B * sP.T));
   }

   std::vector<vmesh::GlobalID> Ionosphere::findBlocksToInitialize(spatial_cell::SpatialCell& cell,const uint popID) {
      vector<vmesh::GlobalID> blocksToInitialize;
      bool search = true;
      uint counter = 0;
      const uint8_t refLevel = 0;

      const vmesh::LocalID* vblocks_ini = cell.get_velocity_grid_length(popID,refLevel);

      while (search) {
         #warning TODO: add SpatialCell::getVelocityBlockMinValue() in place of sparseMinValue ? (if applicable)
         if (0.1 * getObjectWrapper().particleSpecies[popID].sparseMinValue > 
            shiftedMaxwellianDistribution(popID,counter*cell.get_velocity_grid_block_size(popID,refLevel)[0], 0.0, 0.0)
            || counter > vblocks_ini[0]) {
            search = false;
         }
         ++counter;
      }
      counter+=2;
      Real vRadiusSquared 
              = (Real)counter*(Real)counter
              * cell.get_velocity_grid_block_size(popID,refLevel)[0]
              * cell.get_velocity_grid_block_size(popID,refLevel)[0];

      for (uint kv=0; kv<vblocks_ini[2]; ++kv) 
         for (uint jv=0; jv<vblocks_ini[1]; ++jv)
            for (uint iv=0; iv<vblocks_ini[0]; ++iv) {
               vmesh::LocalID blockIndices[3];
               blockIndices[0] = iv;
               blockIndices[1] = jv;
               blockIndices[2] = kv;
               const vmesh::GlobalID blockGID = cell.get_velocity_block(popID,blockIndices,refLevel);
               Real blockCoords[3];
               cell.get_velocity_block_coordinates(popID,blockGID,blockCoords);
               Real blockSize[3];
               cell.get_velocity_block_size(popID,blockGID,blockSize);
               blockCoords[0] += 0.5*blockSize[0];
               blockCoords[1] += 0.5*blockSize[1];
               blockCoords[2] += 0.5*blockSize[2];
               //creal vx = P::vxmin + (iv+0.5) * cell.get_velocity_grid_block_size(popID)[0]; // vx-coordinate of the centre
               //creal vy = P::vymin + (jv+0.5) * cell.get_velocity_grid_block_size(popID)[1]; // vy-
               //creal vz = P::vzmin + (kv+0.5) * cell.get_velocity_grid_block_size(popID)[2]; // vz-
               
               if (blockCoords[0]*blockCoords[0] + blockCoords[1]*blockCoords[1] + blockCoords[2]*blockCoords[2] < vRadiusSquared) {
               //if (vx*vx + vy*vy + vz*vz < vRadiusSquared) {
                  // Adds velocity block to active population's velocity mesh
                  //const vmesh::GlobalID newBlockGID = cell.get_velocity_block(popID,vx,vy,vz);
                  cell.add_velocity_block(blockGID,popID);
                  blocksToInitialize.push_back(blockGID);
               }
            }

      return blocksToInitialize;
   }

   void Ionosphere::setCellFromTemplate(SpatialCell* cell,const uint popID) {
      copyCellData(&templateCell,cell,false,popID,true); // copy also vdf, _V
      copyCellData(&templateCell,cell,true,popID,false); // don't copy vdf again but copy _R now
   }

   std::string Ionosphere::getName() const {return "Ionosphere";}
   
   uint Ionosphere::getIndex() const {return sysboundarytype::IONOSPHERE;}
}
