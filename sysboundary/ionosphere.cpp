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

/*!\file ionosphere.cpp
 * \brief Implementation of the class SysBoundaryCondition::Ionosphere to handle cells classified as sysboundarytype::IONOSPHERE.
 */

#include <cstdlib>
#include <iostream>

#include "ionosphere.h"
#include "../projects/project.h"
#include "../projects/projects_common.h"
#include "../vlasovmover.h"
#include "../common.h"

namespace SBC {
   Ionosphere::Ionosphere(): SysBoundaryCondition() { }
   
   Ionosphere::~Ionosphere() { }
   
   void Ionosphere::addParameters() {
      Readparameters::add("ionosphere.centerX", "X coordinate of ionosphere center (m)", 0.0);
      Readparameters::add("ionosphere.centerY", "Y coordinate of ionosphere center (m)", 0.0);
      Readparameters::add("ionosphere.centerZ", "Z coordinate of ionosphere center (m)", 0.0);
      Readparameters::add("ionosphere.radius", "Radius of ionosphere (m).", 1.0e7);
      Readparameters::add("ionosphere.rho", "Number density of the ionosphere (m^-3)", 1.0e6);
      Readparameters::add("ionosphere.depth", "Depth in cells of ionosphere layer.", 1);
      Readparameters::add("ionosphere.taperRadius", "Width of the zone with a density tapering from the ionospheric value to the background (m)", 0.0);
      Readparameters::add("ionosphere.precedence", "Precedence value of the ionosphere system boundary condition (integer), the higher the stronger.", 2);
   }
   
   void Ionosphere::getParameters() {
      int myRank;
      MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
      if(!Readparameters::get("ionosphere.centerX", this->center[0])) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!Readparameters::get("ionosphere.centerY", this->center[1])) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!Readparameters::get("ionosphere.centerZ", this->center[2])) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!Readparameters::get("ionosphere.radius", this->radius)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!Readparameters::get("ionosphere.rho", this->rho)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!Readparameters::get("Magnetosphere.T", this->T)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!Readparameters::get("ionosphere.depth", this->depth)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!Readparameters::get("ionosphere.precedence", this->precedence)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!Readparameters::get("Magnetosphere.nSpaceSamples", this->nSpaceSamples)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!Readparameters::get("Magnetosphere.nVelocitySamples", this->nVelocitySamples)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
   }
   
   bool Ionosphere::initSysBoundary(
      creal& t,
      Project &project
   ) {
      getParameters();
      isThisDynamic = false;
      
      generateTemplateCell(project);
      
      return true;
   }
   
   bool Ionosphere::assignSysBoundary(dccrg::Dccrg<SpatialCell>& mpiGrid) {
      vector<CellID> cells = mpiGrid.get_cells();
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
         creal r = sqrt((x-center[0])*(x-center[0]) + (y-center[1])*(y-center[1]) + (z-center[2])*(z-center[2]));
         
         if(r < radius) {
            mpiGrid[cells[i]]->sysBoundaryFlag = this->getIndex();
         }
      }
      
      SpatialCell::set_mpi_transfer_type(Transfer::CELL_SYSBOUNDARYFLAG);
      mpiGrid.update_remote_neighbor_data(SYSBOUNDARIES_NEIGHBORHOOD_ID);
      
      vector<bool> iCanHasDoNotCompute(cells.size(), true);
      for(uint i=0; i<cells.size(); i++) {
         if(mpiGrid[cells[i]]->sysBoundaryFlag != this->getIndex()) {
            continue;
         }
         
         for(int ci=-depth; ci<=depth; ci++)
            for(int cj=-depth; cj<=depth; cj++)
               for(int ck=-depth; ck<=depth; ck++) {
                  if(ci == 0 && cj == 0 && ck == 0) continue;
                  CellID tmpCellID = getNeighbour(mpiGrid, cells[i], ci, cj, ck);
                  if((tmpCellID != INVALID_CELLID) &&
                     (mpiGrid[tmpCellID]->sysBoundaryFlag != this->getIndex())) {
                     iCanHasDoNotCompute[i] = false;
                  }
         }
      }
      for(uint i=0; i<cells.size(); i++) {
         if(mpiGrid[cells[i]]->sysBoundaryFlag != this->getIndex()) {
            continue;
         }
         if(iCanHasDoNotCompute[i]) {
            mpiGrid[cells[i]]->sysBoundaryFlag = sysboundarytype::DO_NOT_COMPUTE;
         }
      }
      
      return true;
   }
   
   bool Ionosphere::applyInitialState(
      const dccrg::Dccrg<SpatialCell>& mpiGrid,
      Project &project
   ) {
      vector<uint64_t> cells = mpiGrid.get_cells();
#pragma omp parallel for
      for (uint i=0; i<cells.size(); ++i) {
         SpatialCell* cell = mpiGrid[cells[i]];
         if(cell->sysBoundaryFlag != this->getIndex()) continue;
         setCellFromTemplate(cell);
      }
      return true;
   }
   
//    bool Ionosphere::applySysBoundaryCondition(
//       const dccrg::Dccrg<SpatialCell>& mpiGrid,
//       creal& t
//    ) {
//       return true;
//    }
   
   Real Ionosphere::fieldSolverBoundaryCondMagneticField(
      const dccrg::Dccrg<SpatialCell>& mpiGrid,
      const CellID& cellID,
      creal& dt,
      cuint& component
   ) {
      // The perturbed magnetic field is reset to 0.0, the dipole field is in the background component.
      return 0.0;
   }
   
   void Ionosphere::fieldSolverBoundaryCondElectricField(
      dccrg::Dccrg<SpatialCell>& mpiGrid,
      const CellID& cellID,
      cuint RKCase,
      cuint component
   ) {
      if((RKCase == RK_ORDER1) || (RKCase == RK_ORDER2_STEP2)) {
         mpiGrid[cellID]->parameters[CellParams::EX+component] = 0.0;
      } else {// RKCase == RK_ORDER2_STEP1
         mpiGrid[cellID]->parameters[CellParams::EX_DT2+component] = 0.0;
      }
   }
   
   void Ionosphere::fieldSolverBoundaryCondDerivatives(
      const dccrg::Dccrg<SpatialCell>& mpiGrid,
      const CellID& cellID,
      cuint& component
   ) {
      // WARNING This is crap!!
      this->setCellDerivativesToZero(mpiGrid, cellID, component);
   }
   
   void Ionosphere::fieldSolverBoundaryCondBVOLDerivatives(
      const dccrg::Dccrg<SpatialCell>& mpiGrid,
      const CellID& cellID,
      cuint& component
   ) {
      // WARNING This is crap!!
      this->setCellBVOLDerivativesToZero(mpiGrid, cellID, component);
   }
   
   void Ionosphere::vlasovBoundaryCondition(
      const dccrg::Dccrg<SpatialCell>& mpiGrid,
      const CellID& cellID
   ) {
      phiprof::start("vlasovBoundaryCondition (Ionosphere)");
      copyCellData(&templateCell, mpiGrid[cellID]);
      phiprof::stop("vlasovBoundaryCondition (Ionosphere)");
   }
   
   void Ionosphere::generateTemplateCell(Project &project) {
      // WARNING not 0.0 here or the dipole() function fails miserably.
      templateCell.sysBoundaryFlag = this->getIndex();
      templateCell.parameters[CellParams::XCRD] = 1.0;
      templateCell.parameters[CellParams::YCRD] = 1.0;
      templateCell.parameters[CellParams::ZCRD] = 1.0;
      templateCell.parameters[CellParams::DX] = 1;
      templateCell.parameters[CellParams::DY] = 1;
      templateCell.parameters[CellParams::DZ] = 1;
      
      vector<uint> blocksToInitialize = this->findBlocksToInitialize(templateCell);
      
      for(uint i = 0; i < blocksToInitialize.size(); i++) {
         Velocity_Block* blockPtr = templateCell.at(blocksToInitialize.at(i));
         creal vxBlock = blockPtr->parameters[BlockParams::VXCRD];
         creal vyBlock = blockPtr->parameters[BlockParams::VYCRD];
         creal vzBlock = blockPtr->parameters[BlockParams::VZCRD];
         creal dvxCell = SpatialCell::cell_dvx; // Size of one cell in a block in vx-direction
         creal dvyCell = SpatialCell::cell_dvy; //                                vy
         creal dvzCell = SpatialCell::cell_dvz; //                                vz
         
         creal x = templateCell.parameters[CellParams::XCRD];
         creal y = templateCell.parameters[CellParams::YCRD];
         creal z = templateCell.parameters[CellParams::ZCRD];
         creal dx = templateCell.parameters[CellParams::DX];
         creal dy = templateCell.parameters[CellParams::DY];
         creal dz = templateCell.parameters[CellParams::DZ];
         
         // Calculate volume average of distrib. function for each cell in the block.
         for (uint kc=0; kc<WID; ++kc) 
            for (uint jc=0; jc<WID; ++jc) 
               for (uint ic=0; ic<WID; ++ic) {
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
                              average += maxwellianDistribution(
                                 vxCell + vi*d_vx,
                                 vyCell + vj*d_vy,
                                 vzCell + vk*d_vz
                              );
                           }
                           average /= this->nVelocitySamples * this->nVelocitySamples * this->nVelocitySamples;
                  } else {
                     average = maxwellianDistribution(
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
         }
      }
      //let's get rid of blocks not fulfilling the criteria here to save
      //memory.
      templateCell.adjustSingleCellVelocityBlocks();
      
      calculateCellVelocityMoments(&templateCell);
      
      // WARNING Time-independence assumed here. Normal momentes computed in setProjectCell
      templateCell.parameters[CellParams::RHO_DT2] = templateCell.parameters[CellParams::RHO];
      templateCell.parameters[CellParams::RHOVX_DT2] = templateCell.parameters[CellParams::RHOVX];
      templateCell.parameters[CellParams::RHOVY_DT2] = templateCell.parameters[CellParams::RHOVY];
      templateCell.parameters[CellParams::RHOVZ_DT2] = templateCell.parameters[CellParams::RHOVZ];
   }
   
   Real Ionosphere::maxwellianDistribution(
      creal& vx, creal& vy, creal& vz
   ) {
      return this->rho * pow(physicalconstants::MASS_PROTON /
      (2.0 * M_PI * physicalconstants::K_B * this->T), 1.5) *
      exp(-physicalconstants::MASS_PROTON * (vx*vx + vy*vy + vz*vz) /
      (2.0 * physicalconstants::K_B * this->T));
   }
   
   vector<uint> Ionosphere::findBlocksToInitialize(
      SpatialCell& cell
   ) {
      vector<uint> blocksToInitialize;
      bool search = true;
      int counter = 0;
      
      while(search) {
         if(0.1 * P::sparseMinValue >
            maxwellianDistribution(counter*SpatialCell::block_dvx, 0.0, 0.0)
         ) {
            search = false;
         }
         counter++;
      }
      counter+=2;
      Real vRadiusSquared = (Real)counter*(Real)counter*SpatialCell::block_dvx*SpatialCell::block_dvx;
      
      for (uint kv=0; kv<P::vzblocks_ini; ++kv) 
         for (uint jv=0; jv<P::vyblocks_ini; ++jv)
            for (uint iv=0; iv<P::vxblocks_ini; ++iv) {
               creal vx = P::vxmin + (iv+0.5) * SpatialCell::block_dvx; // vx-coordinate of the centre
               creal vy = P::vymin + (jv+0.5) * SpatialCell::block_dvy; // vy-
               creal vz = P::vzmin + (kv+0.5) * SpatialCell::block_dvz; // vz-
               
               if(vx*vx + vy*vy + vz*vz < vRadiusSquared) {
                  cell.add_velocity_block(cell.get_velocity_block(vx, vy, vz));
                  blocksToInitialize.push_back(cell.get_velocity_block(vx, vy, vz));
               }
            }
            
            return blocksToInitialize;
   }
   
   void Ionosphere::setCellFromTemplate(SpatialCell *cell) {
      // The ionospheric cell has the same state as the initial state of non-system boundary cells so far.
      cell->parameters[CellParams::RHOLOSSADJUST] = 0.0;
      cell->parameters[CellParams::RHOLOSSVELBOUNDARY] = 0.0;
      
      copyCellData(&templateCell, cell);
   }
   
   std::string Ionosphere::getName() const {return "Ionosphere";}
   
   uint Ionosphere::getIndex() const {return sysboundarytype::IONOSPHERE;}
}
