/*
 * This file is part of Vlasiator.
 * 
 * Copyright 2010, 2011, 2012, 2013 Finnish Meteorological Institute
 */

/*!\file setmaxwellian.cpp
 * \brief Implementation of the class SysBoundaryCondition::SetMaxwellian to handle cells classified as sysboundarytype::MAXWELLIAN.
 */

#include <cstdlib>
#include <iostream>

#include "setmaxwellian.h"
#include "../vlasovmover.h"

namespace SBC {
   SetMaxwellian::SetMaxwellian(): SetByUser() {
      nParams = 9;
   }
   SetMaxwellian::~SetMaxwellian() { }
   
   void SetMaxwellian::addParameters() {
      Readparameters::addComposing("maxwellian.face", "List of faces on which set Maxwellian boundary conditions are to be applied ([xyz][+-]).");
      Readparameters::add("maxwellian.file_x+", "Input files for the set Maxwellian inflow parameters on face x+. Data format per line: time (s) density (p/m^3) Temperature (K) Vx Vy Vz (m/s) Bx By Bz (T).", "");
      Readparameters::add("maxwellian.file_x-", "Input files for the set Maxwellian inflow parameters on face x-. Data format per line: time (s) density (p/m^3) Temperature (K) Vx Vy Vz (m/s) Bx By Bz (T).", "");
      Readparameters::add("maxwellian.file_y+", "Input files for the set Maxwellian inflow parameters on face y+. Data format per line: time (s) density (p/m^3) Temperature (K) Vx Vy Vz (m/s) Bx By Bz (T).", "");
      Readparameters::add("maxwellian.file_y-", "Input files for the set Maxwellian inflow parameters on face y-. Data format per line: time (s) density (p/m^3) Temperature (K) Vx Vy Vz (m/s) Bx By Bz (T).", "");
      Readparameters::add("maxwellian.file_z+", "Input files for the set Maxwellian inflow parameters on face z+. Data format per line: time (s) density (p/m^3) Temperature (K) Vx Vy Vz (m/s) Bx By Bz (T).", "");
      Readparameters::add("maxwellian.file_z-", "Input files for the set Maxwellian inflow parameters on face z-. Data format per line: time (s) density (p/m^3) Temperature (K) Vx Vy Vz (m/s) Bx By Bz (T).", "");
      Readparameters::add("maxwellian.dynamic", "Boolean value, is the set Maxwellian inflow dynamic in time or not.", 0);
      Readparameters::add("maxwellian.precedence", "Precedence value of the set Maxwellian system boundary condition (integer), the higher the stronger.", 3);
      Readparameters::add("maxwellian.nSpaceSamples", "Number of sampling points per spatial dimension (template cells)", 2);
      Readparameters::add("maxwellian.nVelocitySamples", "Number of sampling points per velocity dimension (template cells)", 5);
   }
   
   void SetMaxwellian::getParameters() {
      int myRank;
      MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
      if(!Readparameters::get("maxwellian.face", faceList)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!Readparameters::get("maxwellian.dynamic", isThisDynamic)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!Readparameters::get("maxwellian.file_x+", files[0])) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!Readparameters::get("maxwellian.file_x-", files[1])) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!Readparameters::get("maxwellian.file_y+", files[2])) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!Readparameters::get("maxwellian.file_y-", files[3])) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!Readparameters::get("maxwellian.file_z+", files[4])) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!Readparameters::get("maxwellian.file_z-", files[5])) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!Readparameters::get("maxwellian.precedence", precedence)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!Readparameters::get("maxwellian.nSpaceSamples", nSpaceSamples)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!Readparameters::get("maxwellian.nVelocitySamples", nVelocitySamples)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
   }
   
   Real SetMaxwellian::maxwellianDistribution(
            const int& popID,
            creal& rho,
            creal& T,
            creal& vx, creal& vy, creal& vz
   ) {
      #warning All populations assumed to have the same T
      const Real MASS = getObjectWrapper().particleSpecies[popID].mass;
      return rho * pow(MASS /
      (2.0 * M_PI * physicalconstants::K_B * T), 1.5) *
      exp(-MASS * (vx*vx + vy*vy + vz*vz) /
      (2.0 * physicalconstants::K_B * T));
   }
   
   vector<vmesh::GlobalID> SetMaxwellian::findBlocksToInitialize(
                                                      const int& popID,
                                                      SpatialCell& cell,
                                                      creal& rho,
                                                      creal& T,
                                                      creal& VX0,
                                                      creal& VY0,
                                                      creal& VZ0
                                                     ) {
      vector<vmesh::GlobalID> blocksToInitialize;
      bool search = true;
      uint counter = 0;
      while (search) {
         if (0.1 * P::sparseMinValue >
             maxwellianDistribution(
                                    popID,
                                    rho,
                                    T,
                                    counter*SpatialCell::get_velocity_grid_block_size()[0], 0.0, 0.0
                                   )
             ||
             counter > P::vxblocks_ini
            ) {
            search = false;
         }
         counter++;
      }
      counter+=2;

      Real vRadiusSquared = (Real)counter*(Real)counter*SpatialCell::get_velocity_grid_block_size()[0]*SpatialCell::get_velocity_grid_block_size()[0];
      
      for (uint kv=0; kv<P::vzblocks_ini; ++kv) 
         for (uint jv=0; jv<P::vyblocks_ini; ++jv)
            for (uint iv=0; iv<P::vxblocks_ini; ++iv) {
               creal vx = P::vxmin + (iv+0.5) * SpatialCell::get_velocity_grid_block_size()[0]; // vx-coordinate of the centre
               creal vy = P::vymin + (jv+0.5) * SpatialCell::get_velocity_grid_block_size()[1]; // vy-
               creal vz = P::vzmin + (kv+0.5) * SpatialCell::get_velocity_grid_block_size()[2]; // vz-

               if ((vx-VX0)*(vx-VX0) + (vy-VY0)*(vy-VY0) + (vz-VZ0)*(vz-VZ0) < vRadiusSquared) {
                  cell.add_velocity_block(SpatialCell::get_velocity_block(vx, vy, vz),popID);
                  blocksToInitialize.push_back(SpatialCell::get_velocity_block(vx, vy, vz));
               }
            }

      return blocksToInitialize;
   }
   
   /*!\brief Generate the template cell for the face corresponding to the index passed.
    * This function generates a spatial cell which is to be used as a template for the
    * system boundary condition.
    * \param templateCell Addressof the template cell to be generated.
    * \param inputDataIndex Index used for the location of the input data.
    * \param t Current simulation time.
    */
   void SetMaxwellian::generateTemplateCell(spatial_cell::SpatialCell& templateCell, int inputDataIndex, creal& t) {
      Real rho, T, Vx, Vy, Vz, Bx, By, Bz, buffer[8];
      
      interpolate(inputDataIndex, t, &buffer[0]);
      rho = buffer[0];
      T = buffer[1];
      Vx = buffer[2];
      Vy = buffer[3];
      Vz = buffer[4];
      Bx = buffer[5];
      By = buffer[6];
      Bz = buffer[7];
      
      templateCell.sysBoundaryFlag = this->getIndex();
      templateCell.sysBoundaryLayer = 1;
      
      templateCell.parameters[CellParams::XCRD] = 0.0;
      templateCell.parameters[CellParams::YCRD] = 0.0;
      templateCell.parameters[CellParams::ZCRD] = 0.0;
      templateCell.parameters[CellParams::DX] = 1;
      templateCell.parameters[CellParams::DY] = 1;
      templateCell.parameters[CellParams::DZ] = 1;
      templateCell.parameters[CellParams::PERBX] = Bx;
      templateCell.parameters[CellParams::PERBY] = By;
      templateCell.parameters[CellParams::PERBZ] = Bz;
      
      templateCell.parameters[CellParams::RHOLOSSADJUST] = 0.0;
      templateCell.parameters[CellParams::RHOLOSSVELBOUNDARY] = 0.0;
      
      // Init all particle species
      for (int popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
      
         vector<vmesh::GlobalID> blocksToInitialize = this->findBlocksToInitialize(popID,templateCell, rho, T, Vx, Vy, Vz);
         for(vmesh::GlobalID i = 0; i < blocksToInitialize.size(); i++) {
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
               if(this->nVelocitySamples > 1) {
                  creal d_vx = dvxCell / (nVelocitySamples-1);
                  creal d_vy = dvyCell / (nVelocitySamples-1);
                  creal d_vz = dvzCell / (nVelocitySamples-1);
                  for (uint vi=0; vi<nVelocitySamples; ++vi)
                     for (uint vj=0; vj<nVelocitySamples; ++vj)
                        for (uint vk=0; vk<nVelocitySamples; ++vk) {
                           average +=  maxwellianDistribution(
                                          popID,
                                          rho,
                                          T,
                                          vxCell + vi*d_vx - Vx,
                                          vyCell + vj*d_vy - Vy,
                                          vzCell + vk*d_vz - Vz
                                       );
                        }
                  average /= this->nVelocitySamples * this->nVelocitySamples * this->nVelocitySamples;
               } else {
                  average =   maxwellianDistribution(
                                 popID,
                                 rho,
                                 T,
                                 vxCell + 0.5*dvxCell,
                                 vyCell + 0.5*dvyCell,
                                 vzCell + 0.5*dvzCell
                              );
               }
                  
               if(average!=0.0){
                  creal vxCellCenter = vxBlock + (ic+convert<Real>(0.5))*dvxCell;
                  creal vyCellCenter = vyBlock + (jc+convert<Real>(0.5))*dvyCell;
                  creal vzCellCenter = vzBlock + (kc+convert<Real>(0.5))*dvzCell;
                  templateCell.set_value(vxCellCenter,vyCellCenter,vzCellCenter,average);
               }
            } // for-loop over cells in velocity block
         } // for-loop over velocity blocks

         //let's get rid of blocks not fulfilling the criteria here to save memory.
         templateCell.adjustSingleCellVelocityBlocks(popID);
      } // for-loop over particle species

      calculateCellVelocityMoments(&templateCell, true);
      
      if(!this->isThisDynamic) {
         // WARNING Time-independence assumed here.
         templateCell.parameters[CellParams::RHO_DT2] = templateCell.parameters[CellParams::RHO];
         templateCell.parameters[CellParams::RHOVX_DT2] = templateCell.parameters[CellParams::RHOVX];
         templateCell.parameters[CellParams::RHOVY_DT2] = templateCell.parameters[CellParams::RHOVY];
         templateCell.parameters[CellParams::RHOVZ_DT2] = templateCell.parameters[CellParams::RHOVZ];
         templateCell.parameters[CellParams::PERBX_DT2] = templateCell.parameters[CellParams::PERBX];
         templateCell.parameters[CellParams::PERBY_DT2] = templateCell.parameters[CellParams::PERBY];
         templateCell.parameters[CellParams::PERBZ_DT2] = templateCell.parameters[CellParams::PERBZ];
      } else {
         cerr << "ERROR: this is not dynamic in time, please code it!" << endl;
         abort();
      }
   }
   
   string SetMaxwellian::getName() const {return "SetMaxwellian";}
   uint SetMaxwellian::getIndex() const {return sysboundarytype::SET_MAXWELLIAN;}
}
