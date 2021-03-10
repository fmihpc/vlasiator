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

/*!\file setmaxwellian.cpp
 * \brief Implementation of the class SysBoundaryCondition::SetMaxwellian to handle cells classified as sysboundarytype::MAXWELLIAN.
 */

#include <cstdlib>
#include <iostream>

#include "setmaxwellian.h"
#include "../vlasovmover.h"
#include "../object_wrapper.h"

namespace SBC {
   SetMaxwellian::SetMaxwellian(): SetByUser() {
   }
   SetMaxwellian::~SetMaxwellian() { }
   
   void SetMaxwellian::addParameters() {
      Readparameters::addComposing("maxwellian.face", "List of faces on which set Maxwellian boundary conditions are to be applied ([xyz][+-]).");
      Readparameters::add("maxwellian.precedence", "Precedence value of the set Maxwellian system boundary condition (integer), the higher the stronger.", 3);
      Readparameters::add("maxwellian.reapplyUponRestart", "If 0 (default), keep going with the state existing in the restart file. If 1, calls again applyInitialState. Can be used to change boundary condition behaviour during a run.", 0);

      // Per-population parameters
      for(uint i=0; i< getObjectWrapper().particleSpecies.size(); i++) {
        const std::string& pop = getObjectWrapper().particleSpecies[i].name;

        Readparameters::add(pop + "_maxwellian.file_x+", "Input files for the set Maxwellian inflow parameters on face x+. Data format per line: time (s) density (p/m^3) Temperature (K) Vx Vy Vz (m/s) Bx By Bz (T).", "");
        Readparameters::add(pop + "_maxwellian.file_x-", "Input files for the set Maxwellian inflow parameters on face x-. Data format per line: time (s) density (p/m^3) Temperature (K) Vx Vy Vz (m/s) Bx By Bz (T).", "");
        Readparameters::add(pop + "_maxwellian.file_y+", "Input files for the set Maxwellian inflow parameters on face y+. Data format per line: time (s) density (p/m^3) Temperature (K) Vx Vy Vz (m/s) Bx By Bz (T).", "");
        Readparameters::add(pop + "_maxwellian.file_y-", "Input files for the set Maxwellian inflow parameters on face y-. Data format per line: time (s) density (p/m^3) Temperature (K) Vx Vy Vz (m/s) Bx By Bz (T).", "");
        Readparameters::add(pop + "_maxwellian.file_z+", "Input files for the set Maxwellian inflow parameters on face z+. Data format per line: time (s) density (p/m^3) Temperature (K) Vx Vy Vz (m/s) Bx By Bz (T).", "");
        Readparameters::add(pop + "_maxwellian.file_z-", "Input files for the set Maxwellian inflow parameters on face z-. Data format per line: time (s) density (p/m^3) Temperature (K) Vx Vy Vz (m/s) Bx By Bz (T).", "");
        Readparameters::add(pop + "_maxwellian.nVelocitySamples", "Number of sampling points per velocity dimension (template cells)", 5);
        Readparameters::add(pop + "_maxwellian.dynamic", "Boolean value, is the set Maxwellian inflow dynamic in time or not.", 0);
      }
   }
   
   void SetMaxwellian::getParameters() {
      int myRank;
      MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
      if(!Readparameters::get("maxwellian.face", faceList)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!Readparameters::get("maxwellian.precedence", precedence)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      uint reapply;
      if(!Readparameters::get("maxwellian.reapplyUponRestart",reapply)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      };
      this->applyUponRestart = false;
      if(reapply == 1) {
         this->applyUponRestart = true;
      }

      // Per-population parameters
      for(uint i=0; i< getObjectWrapper().particleSpecies.size(); i++) {
         const std::string& pop = getObjectWrapper().particleSpecies[i].name;

         UserSpeciesParameters sP;
         sP.nParams = 9;

         if(!Readparameters::get(pop + "_maxwellian.dynamic", isThisDynamic)) {
            if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added for population " << pop << "!" << endl;
            exit(1);
         }
         if(!Readparameters::get(pop + "_maxwellian.file_x+", sP.files[0])) {
            if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added for population " << pop << "!" << endl;
            exit(1);
         }
         if(!Readparameters::get(pop + "_maxwellian.file_x-", sP.files[1])) {
            if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added for population " << pop << "!" << endl;
            exit(1);
         }
         if(!Readparameters::get(pop + "_maxwellian.file_y+", sP.files[2])) {
            if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added for population " << pop << "!" << endl;
            exit(1);
         }
         if(!Readparameters::get(pop + "_maxwellian.file_y-", sP.files[3])) {
            if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added for population " << pop << "!" << endl;
            exit(1);
         }
         if(!Readparameters::get(pop + "_maxwellian.file_z+", sP.files[4])) {
            if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added for population " << pop << "!" << endl;
            exit(1);
         }
         if(!Readparameters::get(pop + "_maxwellian.file_z-", sP.files[5])) {
            if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added for population " << pop << "!" << endl;
            exit(1);
         }
         if(!Readparameters::get(pop + "_maxwellian.nVelocitySamples", sP.nVelocitySamples)) {
            if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added for population " << pop << "!" << endl;
            exit(1);
         }

         speciesParams.push_back(sP);
      }
   }
   
   Real SetMaxwellian::maxwellianDistribution(
            const uint popID,
            creal& rho,
            creal& T,
            creal& vx, creal& vy, creal& vz
   ) {
      const Real MASS = getObjectWrapper().particleSpecies[popID].mass;
      return rho * pow(MASS /
      (2.0 * M_PI * physicalconstants::K_B * T), 1.5) *
      exp(-MASS * (vx*vx + vy*vy + vz*vz) /
      (2.0 * physicalconstants::K_B * T));
   }
   
   std::vector<vmesh::GlobalID> SetMaxwellian::findBlocksToInitialize(
      const uint popID,
      spatial_cell::SpatialCell& cell,
      creal& rho,
      creal& T,
      creal& VX0,
      creal& VY0,
      creal& VZ0) {
      vector<vmesh::GlobalID> blocksToInitialize;
      bool search = true;
      uint counter = 0;
      const uint8_t refLevel = 0;
      
      const vmesh::LocalID* vblocks_ini = cell.get_velocity_grid_length(popID,refLevel);

      while (search) {
         #warning TODO: add SpatialCell::getVelocityBlockMinValue() in place of sparseMinValue?
         if (0.1 * getObjectWrapper().particleSpecies[popID].sparseMinValue > 
             maxwellianDistribution(
                                    popID,
                                    rho,
                                    T,
                                    VX0 + counter*cell.get_velocity_grid_block_size(popID,refLevel)[0], VY0, VZ0
                                   )
             ||
             counter > vblocks_ini[0]
            ) {
            search = false;
         }
         counter++;
      }
      counter+=2;

      Real vRadiusSquared 
              = (Real)counter*(Real)counter
              * cell.get_velocity_grid_block_size(popID,refLevel)[0]
              * cell.get_velocity_grid_block_size(popID,refLevel)[0];
      
      for (uint kv=0; kv<vblocks_ini[2]; ++kv) 
         for (uint jv=0; jv<vblocks_ini[1]; ++jv)
            for (uint iv=0; iv<vblocks_ini[0]; ++iv) {
               vmesh::GlobalID blockIndices[3];
               blockIndices[0] = iv;
               blockIndices[1] = jv;
               blockIndices[2] = kv;
               const vmesh::GlobalID blockGID = cell.get_velocity_block(popID,blockIndices,refLevel);
               
               Real V_crds[3];
               cell.get_velocity_block_coordinates(popID,blockGID,V_crds);
               Real dV[3];
               cell.get_velocity_block_size(popID,blockGID,dV);
               V_crds[0] += 0.5*dV[0];
               V_crds[1] += 0.5*dV[1];
               V_crds[2] += 0.5*dV[2];
               Real R2 = ((V_crds[0]-VX0)*(V_crds[0]-VX0)
                       +  (V_crds[1]-VY0)*(V_crds[1]-VY0)
                       +  (V_crds[2]-VZ0)*(V_crds[2]-VZ0));

               if (R2 < vRadiusSquared) {
                  cell.add_velocity_block(blockGID,popID);
                  blocksToInitialize.push_back(blockGID);
               }
            }

      return blocksToInitialize;
   }
   
   /*!\brief Generate the template cell for the face corresponding to the index passed.
    * This function generates a spatial cell which is to be used as a template for the
    * system boundary condition.
    * \param templateCell Address of the template cell to be generated.
    * \param inputDataIndex Index used for the location of the input data.
    * \param t Current simulation time.
    */
   void SetMaxwellian::generateTemplateCell(
      spatial_cell::SpatialCell& templateCell,
      Real B[3],
      int inputDataIndex,
      creal& t
   ) {
      Real rho, T, Vx, Vy, Vz, Bx=0.0, By=0.0, Bz=0.0, buffer[8];
      
      
      templateCell.sysBoundaryFlag = this->getIndex();
      templateCell.sysBoundaryLayer = 1;
      
      templateCell.parameters[CellParams::XCRD] = 0.0;
      templateCell.parameters[CellParams::YCRD] = 0.0;
      templateCell.parameters[CellParams::ZCRD] = 0.0;
      templateCell.parameters[CellParams::DX] = 1;
      templateCell.parameters[CellParams::DY] = 1;
      templateCell.parameters[CellParams::DZ] = 1;
      
      // Init all particle species
      for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
         interpolate(inputDataIndex, popID, t, &buffer[0]);
         rho = buffer[0];
         T = buffer[1];
         Vx = buffer[2];
         Vy = buffer[3];
         Vz = buffer[4];
         Bx = buffer[5];
         By = buffer[6];
         Bz = buffer[7];

         vector<vmesh::GlobalID> blocksToInitialize = this->findBlocksToInitialize(popID,templateCell, rho, T, Vx, Vy, Vz);
         Realf* data = templateCell.get_data(popID);

         for(vmesh::GlobalID i=0; i<blocksToInitialize.size(); ++i) {
            const vmesh::GlobalID blockGID = blocksToInitialize[i];
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
               if(speciesParams[popID].nVelocitySamples > 1) {
                  creal d_vx = dvxCell / (speciesParams[popID].nVelocitySamples-1);
                  creal d_vy = dvyCell / (speciesParams[popID].nVelocitySamples-1);
                  creal d_vz = dvzCell / (speciesParams[popID].nVelocitySamples-1);
                  for (uint vi=0; vi<speciesParams[popID].nVelocitySamples; ++vi)
                    for (uint vj=0; vj<speciesParams[popID].nVelocitySamples; ++vj)
                      for (uint vk=0; vk<speciesParams[popID].nVelocitySamples; ++vk) {
                         average +=  maxwellianDistribution(
                                                            popID,
                                                            rho,
                                                            T,
                                                            vxCell + vi*d_vx - Vx,
                                                            vyCell + vj*d_vy - Vy,
                                                            vzCell + vk*d_vz - Vz
                                                           );
                      }
                  average /= speciesParams[popID].nVelocitySamples * speciesParams[popID].nVelocitySamples * speciesParams[popID].nVelocitySamples;
               } else {
                  average =   maxwellianDistribution(
                                                     popID,
                                                     rho,
                                                     T,
                                                     vxCell + 0.5*dvxCell - Vx,
                                                     vyCell + 0.5*dvyCell - Vy,
                                                     vzCell + 0.5*dvzCell - Vz
                                                    );
               }
               
               if (average != 0.0) {
                  data[blockLID*WID3+cellIndex(ic,jc,kc)] = average;
               } 
            } // for-loop over cells in velocity block
         } // for-loop over velocity blocks
         
         //let's get rid of blocks not fulfilling the criteria here to save
         //memory.
         templateCell.adjustSingleCellVelocityBlocks(popID);
      } // for-loop over particle species
      
      B[0] = Bx;
      B[1] = By;
      B[2] = Bz;
      
      calculateCellMoments(&templateCell,true,true);
      
      if(!this->isThisDynamic) {
         // WARNING Time-independence assumed here.
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
      } else {
         cerr << "ERROR: this is not dynamic in time, please code it!" << endl;
         abort();
      }

   }
   
   string SetMaxwellian::getName() const {return "SetMaxwellian";}
   uint SetMaxwellian::getIndex() const {return sysboundarytype::SET_MAXWELLIAN;}
   
} // namespace SBC
