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
   ProjectBoundary::ProjectBoundary(): SysBoundaryCondition() { }
   ProjectBoundary::~ProjectBoundary() { }
   
   void ProjectBoundary::addParameters() {
      return;
   }
   
   void ProjectBoundary::getParameters() {
      return;
   }
   
   bool ProjectBoundary::initSysBoundary(
      creal& t,
      Project &project
   ) {
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
      
      success = loadInputData();
      success = success & generateTemplateCells(t);
      
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
      bool success;
      
      success = setCellsFromTemplate(mpiGrid);
      
      return true;
   }
   
   Real ProjectBoundary::fieldSolverBoundaryCondMagneticField(
      const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const CellID& cellID,
      creal& dt,
      cuint& component
   ) {
      Real result = 0.0;
      const SpatialCell* cell = mpiGrid[cellID];
      creal dx = cell->parameters[CellParams::DX];
      creal dy = cell->parameters[CellParams::DY];
      creal dz = cell->parameters[CellParams::DZ];
      creal x = cell->parameters[CellParams::XCRD] + 0.5*dx;
      creal y = cell->parameters[CellParams::YCRD] + 0.5*dy;
      creal z = cell->parameters[CellParams::ZCRD] + 0.5*dz;
      
      bool isThisCellOnAFace[6];
      determineFace(&isThisCellOnAFace[0], x, y, z, dx, dy, dz);
      
      for(uint i=0; i<6; i++) {
         if(isThisCellOnAFace[i]) {
            if(dt == 0.0) {
               result = templateCells[i].parameters[CellParams::PERBX + component];
            } else {
               result = templateCells[i].parameters[CellParams::PERBX_DT2 + component];
            }
            break; // This effectively sets the precedence of faces through the order of faces.
         }
      }
      return result;
   }
   
   void ProjectBoundary::fieldSolverBoundaryCondElectricField(
      dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
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
   
   void ProjectBoundary::fieldSolverBoundaryCondHallElectricField(
      dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid,
      const CellID& cellID,
      cuint RKCase,
      cuint component
   ) {
      switch(component) {
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
         const CellID& cellID
      ) {
      SpatialCell* cell = mpiGrid[cellID];
      
      for (int popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
         cell->get_velocity_mesh(popID) = templateCell->get_velocity_mesh(popID);
         cell->get_velocity_blocks(popID) = templateCell->get_velocity_blocks(popID);
      }
   }

   bool ProjectBoundary::setCellsFromTemplate(const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid) {
      vector<uint64_t> cells = mpiGrid.get_cells();
#pragma omp parallel for
      for (uint i=0; i<cells.size(); i++) {
         SpatialCell* cell = mpiGrid[cells[i]];
         if(cell->sysBoundaryFlag != this->getIndex()) continue;
         
         creal dx = cell->parameters[CellParams::DX];
         creal dy = cell->parameters[CellParams::DY];
         creal dz = cell->parameters[CellParams::DZ];
         creal x = cell->parameters[CellParams::XCRD] + 0.5*dx;
         creal y = cell->parameters[CellParams::YCRD] + 0.5*dy;
         creal z = cell->parameters[CellParams::ZCRD] + 0.5*dz;
         
         bool isThisCellOnAFace[6];
         determineFace(&isThisCellOnAFace[0], x, y, z, dx, dy, dz);
         
         for(uint i=0; i<6; i++) {
            if(facesToProcess[i] && isThisCellOnAFace[i]) {
               cell->parameters[CellParams::PERBX] = templateCells[i].parameters[CellParams::PERBX];
               cell->parameters[CellParams::PERBY] = templateCells[i].parameters[CellParams::PERBY];
               cell->parameters[CellParams::PERBZ] = templateCells[i].parameters[CellParams::PERBZ];
               
               cell->parameters[CellParams::RHOLOSSADJUST] = 0.0;
               cell->parameters[CellParams::RHOLOSSVELBOUNDARY] = 0.0;
               
               copyCellData(&templateCells[i], cell,true);
               break; // This effectively sets the precedence of faces through the order of faces.
            }
         }
      }
      return true;
   }
   
   void ProjectBoundary::getFaces(bool* faces) {
      for(uint i=0; i<6; i++) faces[i] = facesToProcess[i];
   }
   
   bool ProjectBoundary::loadInputData() {
      for(uint i=0; i<6; i++) {
         if(facesToProcess[i]) {
            inputData[i] = loadFile(&(files[i][0]));
         } else {
            vector<Real> tmp1;
            vector<vector<Real> > tmp2;
            for(uint j=0; j<nParams; j++) {
               tmp1.push_back(-1.0);
            }
            tmp2.push_back(tmp1);
            inputData[i] = tmp2;
         }
      }
      return true;
   }
   
   /*! Loops through the array of template cells and generates the ones needed. The function
    * generateTemplateCell is defined in the inheriting class such as to have the specific
    * condition needed.
    * \param t Simulation time.
    * \sa generateTemplateCell
    */
   bool ProjectBoundary::generateTemplateCells(creal& t) {
# pragma omp parallel for
      for(uint i=0; i<6; i++) {
         int index;
         if(facesToProcess[i]) {
            generateTemplateCell(templateCells[i], i, t);
         }
      }
      return true;
   }

   void ProjectBoundary::generateTemplateCell(spatial_cell::SpatialCell& templateCell, int inputDataIndex, creal& t) {
      cerr << "Base class ProjectBoundary::generateTemplateCell() called instead of derived class function!" << endl;
   }
   
   string ProjectBoundary::getName() const {
      return "ProjectBoundary";
   }

   uint ProjectBoundary::getIndex() const {
      return sysboundarytype::PROJECT;
   }
}
