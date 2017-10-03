/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://vlasiator.fmi.fi/
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

/*!\file project_boundary.cpp
 * \brief Implementation of the class SysBoundaryCondition::ProjectBoundary. 
 */

#include <cstdlib>
#include <iostream>

#include "project_boundary.h"
#include "../vlasovmover.h"
#include "../fieldsolver/fs_common.h"
#include "../object_wrapper.h"

using namespace std;

#warning Project boundaries do not yet support multipop, not checked whether any changes are needed.

namespace SBC {
   ProjectBoundary::ProjectBoundary(): SysBoundaryCondition() { 
      project = NULL;
   }
   ProjectBoundary::~ProjectBoundary() { 
      project = NULL;
   }
   
   void ProjectBoundary::addParameters() {
      Readparameters::addComposing("projectboundary.face", "List of faces on which outflow boundary conditions are to be applied ([xyz][+-]).");
      Readparameters::add("projectboundary.precedence", "Precedence value of the outflow system boundary condition (integer), the higher the stronger.", 4);
      return;
   }
   
   void ProjectBoundary::getParameters() {
      int myRank;
      MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
      if (!Readparameters::get("projectboundary.face", faceList)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if (!Readparameters::get("projectboundary.precedence", precedence)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      return;
   }
   
   bool ProjectBoundary::initSysBoundary(creal& t,Project& project) {
      this->project = &project;
      
      /* The array of bool describes which of the x+, x-, y+, y-, z+, z- faces are to have user-set system boundary conditions.
       * A true indicates the corresponding face will have user-set system boundary conditions.
       * The 6 elements correspond to x+, x-, y+, y-, z+, z- respectively.
       */
      bool success = true;
      for(uint i=0; i<6; i++) facesToProcess[i] = false;
      
      this->getParameters();
      
      vector<string>::const_iterator it;
      for (it = faceList.begin();
           it != faceList.end();
      it++) {
         if(*it == "x+") facesToProcess[0] = true;
         if(*it == "x-") facesToProcess[1] = true;
         if(*it == "y+") facesToProcess[2] = true;
         if(*it == "y-") facesToProcess[3] = true;
         if(*it == "z+") facesToProcess[4] = true;
         if(*it == "z-") facesToProcess[5] = true;
      }

      success = success & generateTemplateCell();
      return success;
   }
   
   bool ProjectBoundary::assignSysBoundary(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid) {
      vector<CellID> cells = mpiGrid.get_cells();
      for(uint i = 0; i < cells.size(); i++) {
         if(mpiGrid[cells[i]]->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) continue;
         creal* const cellParams = &(mpiGrid[cells[i]]->parameters[0]);
         creal dx = cellParams[CellParams::DX];
         creal dy = cellParams[CellParams::DY];
         creal dz = cellParams[CellParams::DZ];
         creal x = cellParams[CellParams::XCRD] + 0.5*dx;
         creal y = cellParams[CellParams::YCRD] + 0.5*dy;
         creal z = cellParams[CellParams::ZCRD] + 0.5*dz;
         
         bool isThisCellOnAFace[6];
         determineFace(&isThisCellOnAFace[0], x, y, z, dx, dy, dz);
         // Comparison of the array defining which faces to use and the array telling on which faces this cell is
         bool doAssign = false;
         for(uint j=0; j<6; j++) doAssign = doAssign || (facesToProcess[j] && isThisCellOnAFace[j]);
         if(doAssign) {
            mpiGrid[cells[i]]->sysBoundaryFlag = this->getIndex();
         }
      }
      return true;
   }
   
   bool ProjectBoundary::applyInitialState(
                                           const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                           Project &project
                                          ) {
      bool success = true;
      const vector<CellID>& cells = getLocalCells();

      for (size_t c=0; c<cells.size(); ++c) {
         SpatialCell* cell = mpiGrid[cells[c]];
         if (cell->sysBoundaryFlag != getIndex()) continue;
         this->project->setCell(cell);
      }

      return success;
   }
   
   Real ProjectBoundary::fieldSolverBoundaryCondMagneticField(
                                                              const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                                              const CellID& cellID,
                                                              creal& dt,
                                                              cuint& component
                                                             ) {
      return 0.0;
   }

   void ProjectBoundary::fieldSolverBoundaryCondElectricField(
                                                              dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                                              const CellID& cellID,
                                                              cuint RKCase,
                                                              cuint component
                                                             ) {
      dccrg::Types<3>::indices_t indices = mpiGrid.mapping.get_indices(cellID);
      SpatialCell* nbr = NULL;
      SpatialCell* cell = mpiGrid[cellID];
      
      int offset=0;
      if ((RKCase == RK_ORDER1) || (RKCase == RK_ORDER2_STEP2)) offset = CellParams::EX;
      else offset = CellParams::EX_DT2;

      creal dx = cell->parameters[CellParams::DX];
      creal dy = cell->parameters[CellParams::DY];
      creal dz = cell->parameters[CellParams::DZ];
      creal x = cell->parameters[CellParams::XCRD] + 0.5*dx;
      creal y = cell->parameters[CellParams::YCRD] + 0.5*dy;
      creal z = cell->parameters[CellParams::ZCRD] + 0.5*dz;
      bool isThisCellOnAFace[6];
      determineFace(&isThisCellOnAFace[0], x, y, z, dx, dy, dz);
      
      for (uint i=0; i<6; i++) if (facesToProcess[i] && isThisCellOnAFace[i]) {
         switch (i) {
          case 0:
            --indices[0];
            nbr = mpiGrid[ mpiGrid.mapping.get_cell_from_indices(indices,0) ];
            cell->parameters[offset+component] = 0;
            break;
          case 1:
            ++indices[0];
            nbr = mpiGrid[ mpiGrid.mapping.get_cell_from_indices(indices,0) ];
            cell->parameters[offset+component] = 0;
            break;
          case 2:
            --indices[1];
            nbr = mpiGrid[ mpiGrid.mapping.get_cell_from_indices(indices,0) ];
            cell->parameters[offset+component] = 0;
            break;
          case 3:
            ++indices[1];
            nbr = mpiGrid[ mpiGrid.mapping.get_cell_from_indices(indices,0) ];
            cell->parameters[offset+component] = 0;
            break;
          case 4:
            --indices[2];
            nbr = mpiGrid[ mpiGrid.mapping.get_cell_from_indices(indices,0) ];
            cell->parameters[offset+component] = 0;
            break;
          case 5:
            ++indices[2];
            nbr = mpiGrid[ mpiGrid.mapping.get_cell_from_indices(indices,0) ];
            cell->parameters[offset+component] = 0;
            break;
         }         
      }
   }

   void ProjectBoundary::fieldSolverBoundaryCondHallElectricField(
                                                                  dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid,
                                                                  const CellID& cellID,
                                                                  cuint RKCase,
                                                                  cuint component
                                                                 ) {
      switch (component) {
       case 0:
         mpiGrid[cellID]->parameters[CellParams::EXHALL_000_100] = 0.0;
         mpiGrid[cellID]->parameters[CellParams::EXHALL_010_110] = 0.0;
         mpiGrid[cellID]->parameters[CellParams::EXHALL_001_101] = 0.0;
         mpiGrid[cellID]->parameters[CellParams::EXHALL_011_111] = 0.0;
         break;
       case 1:
         mpiGrid[cellID]->parameters[CellParams::EYHALL_000_010] = 0.0;
         mpiGrid[cellID]->parameters[CellParams::EYHALL_100_110] = 0.0;
         mpiGrid[cellID]->parameters[CellParams::EYHALL_001_011] = 0.0;
         mpiGrid[cellID]->parameters[CellParams::EYHALL_101_111] = 0.0;
         break;
       case 2:
         mpiGrid[cellID]->parameters[CellParams::EZHALL_000_001] = 0.0;
         mpiGrid[cellID]->parameters[CellParams::EZHALL_100_101] = 0.0;
         mpiGrid[cellID]->parameters[CellParams::EZHALL_010_011] = 0.0;
         mpiGrid[cellID]->parameters[CellParams::EZHALL_110_111] = 0.0;
         break;
       default:
         cerr << __FILE__ << ":" << __LINE__ << ":" << " Invalid component" << endl;
      }
   }
   
   void ProjectBoundary::fieldSolverBoundaryCondDerivatives(
                                                            dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                                            const CellID& cellID,
                                                            cuint& RKCase,
                                                            cuint& component
                                                           ) {
      this->setCellDerivativesToZero(mpiGrid, cellID, component);
   }
   
   void ProjectBoundary::fieldSolverBoundaryCondBVOLDerivatives(
                                                                const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                                                const CellID& cellID,
                                                                cuint& component
                                                               ) {
      this->setCellBVOLDerivativesToZero(mpiGrid, cellID, component);
   }
   
   void ProjectBoundary::vlasovBoundaryCondition(
                                                 const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                                 const CellID& cellID,
                                                 const uint popID
                                                ) {
      SpatialCell* cell = mpiGrid[cellID];
      cell->get_velocity_mesh(popID)   = templateCell.get_velocity_mesh(popID);
      cell->get_velocity_blocks(popID) = templateCell.get_velocity_blocks(popID);
   }

   void ProjectBoundary::getFaces(bool* faces) {
      for(uint i=0; i<6; i++) faces[i] = facesToProcess[i];
   }
   
   bool ProjectBoundary::generateTemplateCell() {
      if (project == NULL) return false;
      templateCell.parameters[CellParams::XCRD] = -2*Parameters::xmin;
      templateCell.parameters[CellParams::YCRD] = -2*Parameters::ymin;
      templateCell.parameters[CellParams::ZCRD] = -2*Parameters::zmin;
      templateCell.parameters[CellParams::DX] = -2*Parameters::dx_ini;
      templateCell.parameters[CellParams::DY] = -2*Parameters::dy_ini;
      templateCell.parameters[CellParams::DZ] = -2*Parameters::dz_ini;
      project->setCell(&templateCell);
      return true;
   }
   
   string ProjectBoundary::getName() const {
      return "ProjectBoundary";
   }

   uint ProjectBoundary::getIndex() const {
      return sysboundarytype::PROJECT;
   }
}
