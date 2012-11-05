/*
 * This file is part of Vlasiator.
 * 
 * Copyright 2010, 2011, 2012 Finnish Meteorological Institute
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
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

/*!\file setmaxwellian.cpp
 * \brief Implementation of the class SysBoundaryCondition::SetMaxwellian to handle cells classified as sysboundarytype::MAXWELLIAN.
 */

#include <cstdlib>
#include <iostream>

#include "setmaxwellian.h"
#include "../vlasovmover.h"


using namespace std;

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
      
      // Go through each velocity block in the velocity phase space grid.
      // Set the initial state and block parameters:
      creal dvx_block = spatial_cell::SpatialCell::block_dvx; // Size of a block in vx-direction
      creal dvy_block = spatial_cell::SpatialCell::block_dvy; //                    vy
      creal dvz_block = spatial_cell::SpatialCell::block_dvz; //                    vz
      creal dvx_blockCell = spatial_cell::SpatialCell::cell_dvx; // Size of one cell in a block in vx-direction
      creal dvy_blockCell = spatial_cell::SpatialCell::cell_dvy; //                                vy
      creal dvz_blockCell = spatial_cell::SpatialCell::cell_dvz; //                                vz
      
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
                        
                        creal d_vx = dvx_blockCell / (nVelocitySamples-1);
                        creal d_vy = dvy_blockCell / (nVelocitySamples-1);
                        creal d_vz = dvz_blockCell / (nVelocitySamples-1);
                        
                        Real average = 0.0;
                        for (uint vi=0; vi<nVelocitySamples; ++vi)
                           for (uint vj=0; vj<nVelocitySamples; ++vj)
                              for (uint vk=0; vk<nVelocitySamples; ++vk) {
                                 average += rho * pow(physicalconstants::MASS_PROTON / 
                                       (2.0 * M_PI * physicalconstants::K_B * T), 1.5) *
                                       exp(-physicalconstants::MASS_PROTON *
                                       (pow(vx_cell + vi*d_vx - Vx, 2.0) + 
                                       pow(vy_cell + vj*d_vy - Vy, 2.0) + 
                                       pow(vz_cell + vk*d_vz - Vz, 2.0)) / 
                                       (2.0 * physicalconstants::K_B * T));
                              }
                        average /= nVelocitySamples * nVelocitySamples * nVelocitySamples;
                        
                        if(average!=0.0){
                           creal vx_cell_center = vx_block + (ic+convert<Real>(0.5))*dvx_blockCell;
                           creal vy_cell_center = vy_block + (jc+convert<Real>(0.5))*dvy_blockCell;
                           creal vz_cell_center = vz_block + (kc+convert<Real>(0.5))*dvz_blockCell;
                           templateCell.set_value(vx_cell_center,vy_cell_center,vz_cell_center,average);
                        }
               }
      }
      calculateCellVelocityMoments(&templateCell);
      
      // WARNING Time-independence assumed here.
      templateCell.parameters[CellParams::RHO_DT2] = templateCell.parameters[CellParams::RHO];
      templateCell.parameters[CellParams::RHOVX_DT2] = templateCell.parameters[CellParams::RHOVX];
      templateCell.parameters[CellParams::RHOVY_DT2] = templateCell.parameters[CellParams::RHOVY];
      templateCell.parameters[CellParams::RHOVZ_DT2] = templateCell.parameters[CellParams::RHOVZ];
      templateCell.parameters[CellParams::RHO_R] = templateCell.parameters[CellParams::RHO];
      templateCell.parameters[CellParams::RHOVX_R] = templateCell.parameters[CellParams::RHOVX];
      templateCell.parameters[CellParams::RHOVY_R] = templateCell.parameters[CellParams::RHOVY];
      templateCell.parameters[CellParams::RHOVZ_R] = templateCell.parameters[CellParams::RHOVZ];
      templateCell.parameters[CellParams::RHO_V] = templateCell.parameters[CellParams::RHO];
      templateCell.parameters[CellParams::RHOVX_V] = templateCell.parameters[CellParams::RHOVX];
      templateCell.parameters[CellParams::RHOVY_V] = templateCell.parameters[CellParams::RHOVY];
      templateCell.parameters[CellParams::RHOVZ_V] = templateCell.parameters[CellParams::RHOVZ];
      
      //let's get rid of blocks not fulfilling the criteria here to save
      //memory.
      templateCell.adjustSingleCellVelocityBlocks();
   }
   
   string SetMaxwellian::getName() const {return "SetMaxwellian";}
   uint SetMaxwellian::getIndex() const {return sysboundarytype::SET_MAXWELLIAN;}
}
