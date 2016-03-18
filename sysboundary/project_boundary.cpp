/*
 * This file is part of Vlasiator.
 * 
 * Copyright 2015 Finnish Meteorological Institute
 * 
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
      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 3, 2> & perBGrid,
      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 3, 2> & perBDt2Grid,
                                                              const CellID& cellID,
                                                              creal& dt,
                                                              cuint& component
                                                             ) {
      return 0.0;
   }

   void ProjectBoundary::fieldSolverBoundaryCondElectricField(
      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 3, 2> & EGrid,
      cint i,
      cint j,
      cint k,
      cuint component
   ) {
      SpatialCell* nbr = NULL;
      SpatialCell* cell = mpiGrid[cellID];
      
      creal dx = EGrid.DX;
      creal dy = EGrid.DY;
      creal dz = EGrid.DZ;
#warning the positions have to be got from the new grid
      creal x = cell->parameters[CellParams::XCRD] + 0.5*dx;
      creal y = cell->parameters[CellParams::YCRD] + 0.5*dy;
      creal z = cell->parameters[CellParams::ZCRD] + 0.5*dz;
      bool isThisCellOnAFace[6];
      determineFace(&isThisCellOnAFace[0], x, y, z, dx, dy, dz);
      
#warning this function makes little sense as is. It made little before, too.
      for (uint fi=0; fi<6; fi++) if (facesToProcess[fi] && isThisCellOnAFace[fi]) {
         switch (fi) {
          case 0:
            EGrid.get(i,j,k)[fsgrids::efield::EX+component] = 0;
            break;
          case 1:
            EGrid.get(i,j,k)[fsgrids::efield::EX+component] = 0;
            break;
          case 2:
            EGrid.get(i,j,k)[fsgrids::efield::EX+component] = 0;
            break;
          case 3:
            EGrid.get(i,j,k)[fsgrids::efield::EX+component] = 0;
            break;
          case 4:
            EGrid.get(i,j,k)[fsgrids::efield::EX+component] = 0;
            break;
          case 5:
            EGrid.get(i,j,k)[fsgrids::efield::EX+component] = 0;
            break;
         }
      }
   }
   
   void ProjectBoundary::fieldSolverBoundaryCondGradPeElectricField(
      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, 3, 2> & EGradPeGrid,
      cint i,
      cint j,
      cint k,
      cuint component
   ) {
      EGradPeGrid.get(i,j,k)[fsgrids::egradpe::EXGRADPE+component] = 0.0;
   }
   
   void ProjectBoundary::fieldSolverBoundaryCondHallElectricField(
      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, 3, 2> & EHallGrid,
      cint i,
      cint j,
      cint k,
      cuint component
   ) {
      const std::array<Real, fsgrids::ehall::N_EHALL> * cp = EHallGrid.get(i,j,k);
      switch (component) {
         case 0:
            cp[fsgrids::ehall::EXHALL_000_100] = 0.0;
            cp[fsgrids::ehall::EXHALL_010_110] = 0.0;
            cp[fsgrids::ehall::EXHALL_001_101] = 0.0;
            cp[fsgrids::ehall::EXHALL_011_111] = 0.0;
            break;
         case 1:
            cp[fsgrids::ehall::EYHALL_000_010] = 0.0;
            cp[fsgrids::ehall::EYHALL_100_110] = 0.0;
            cp[fsgrids::ehall::EYHALL_001_011] = 0.0;
            cp[fsgrids::ehall::EYHALL_101_111] = 0.0;
            break;
         case 2:
            cp[fsgrids::ehall::EZHALL_000_001] = 0.0;
            cp[fsgrids::ehall::EZHALL_100_101] = 0.0;
            cp[fsgrids::ehall::EZHALL_010_011] = 0.0;
            cp[fsgrids::ehall::EZHALL_110_111] = 0.0;
            break;
         default:
            cerr << __FILE__ << ":" << __LINE__ << ":" << " Invalid component" << endl;
   }
   
   void ProjectBoundary::fieldSolverBoundaryCondDerivatives(
      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, 3, 2> & dPerBGrid,
      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 3, 2> & dMomentsGrid,
                                                            const CellID& cellID,
                                                            cuint& RKCase,
                                                            cuint& component
                                                           ) {
      this->setCellDerivativesToZero(dPerBGrid, dMomentsGrid, cellID, component);
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
                                                 const int& popID
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
