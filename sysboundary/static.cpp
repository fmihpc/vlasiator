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

/*!\file static.cpp
 * \brief Implementation of the class SysBoundaryCondition::Static to handle cells classified as sysboundarytype::STATIC.
 */

#include <cstdlib>
#include <iostream>

#include "../object_wrapper.h"
#include "static.h"
#include "../projects/projects_common.h"
#include "../fieldsolver/fs_common.h"
#include "../fieldsolver/ldz_magnetic_field.hpp"
#include "../vlasovmover.h"

#ifndef NDEBUG
   #define DEBUG_STATIC
#endif
#ifdef DEBUG_SYSBOUNDARY
   #define DEBUG_STATIC
#endif

using namespace std;

namespace SBC {
   Static::Static(): SysBoundaryCondition() { }
   Static::~Static() { }
   
   void Static::addParameters() {
      Readparameters::add("static.precedence", "Precedence value of the static system boundary condition (integer), the higher the stronger.", 4);
      //Readparameters::add("static.reapplyUponRestart", "If 0 (default), keep going with the state existing in the restart file. If 1, calls again applyInitialState. Can be used to change boundary condition behaviour during a run.", 0);

      // Per-population parameters
      for(uint i=0; i< getObjectWrapper().particleSpecies.size(); i++) {
        const std::string& pop = getObjectWrapper().particleSpecies[i].name;

        //Readparameters::addComposing(pop + "_static.reapplyFaceUponRestart", "List of faces on which static boundary conditions are to be reapplied upon restart ([xyz][+-]).");
        Readparameters::addComposing(pop + "_static.face", "List of faces on which static boundary conditions are to be applied ([xyz][+-]).");

      }
   }
   
   void Static::getParameters() {
      int myRank;
      MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
      if(!Readparameters::get("static.precedence", precedence)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }

      // Per-species parameters
      for(uint i=0; i< getObjectWrapper().particleSpecies.size(); i++) {
        const std::string& pop = getObjectWrapper().particleSpecies[i].name;
        StaticSpeciesParameters sP;

        // Unless we find out otherwise, we assume that this species will not be treated at any boundary
        for(int j=0; j<6; j++) {
          sP.facesToSkipVlasov[j] = true;
        }

        std::vector<std::string> thisSpeciesFaceList;
        if(!Readparameters::get(pop + "_static.face", thisSpeciesFaceList)) {
           if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added for population " << pop << "!" << endl;
           exit(1);
        }

        for(auto& face : thisSpeciesFaceList) {
          if(face == "x+") { facesToProcess[0] = true; sP.facesToSkipVlasov[0] = false; }
          if(face == "x-") { facesToProcess[1] = true; sP.facesToSkipVlasov[1] = false; }
          if(face == "y+") { facesToProcess[2] = true; sP.facesToSkipVlasov[2] = false; }
          if(face == "y-") { facesToProcess[3] = true; sP.facesToSkipVlasov[3] = false; }
          if(face == "z+") { facesToProcess[4] = true; sP.facesToSkipVlasov[4] = false; }
          if(face == "z-") { facesToProcess[5] = true; sP.facesToSkipVlasov[5] = false; }
        }
        speciesParams.push_back(sP);
      }
   }
   
   bool Static::initSysBoundary(
      creal& t,
      Project &project
   ) {
      /* The array of bool describes which of the x+, x-, y+, y-, z+, z- faces are to have static system boundary conditions.
       * A true indicates the corresponding face will have static.
       * The 6 elements correspond to x+, x-, y+, y-, z+, z- respectively.
       */
      for(uint i=0; i<6; i++) {
         facesToProcess[i] = false;
      }      
      this->getParameters();

      return true;
   }
   
   bool Static::assignSysBoundary(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                   FsGrid< fsgrids::technical, 2> & technicalGrid) {

      bool doAssign;
      std::array<bool,6> isThisCellOnAFace;
      
      // Assign boundary flags to local DCCRG cells
      vector<CellID> cells = mpiGrid.get_cells();
      for(const auto& dccrgId : cells) {
         if(mpiGrid[dccrgId]->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) continue;
         creal* const cellParams = &(mpiGrid[dccrgId]->parameters[0]);
         creal dx = cellParams[CellParams::DX];
         creal dy = cellParams[CellParams::DY];
         creal dz = cellParams[CellParams::DZ];
         creal x = cellParams[CellParams::XCRD] + 0.5*dx;
         creal y = cellParams[CellParams::YCRD] + 0.5*dy;
         creal z = cellParams[CellParams::ZCRD] + 0.5*dz;
         
         isThisCellOnAFace.fill(false);
         determineFace(isThisCellOnAFace.data(), x, y, z, dx, dy, dz);
         
         // Comparison of the array defining which faces to use and the array telling on which faces this cell is
         doAssign = false;
         for(int j=0; j<6; j++) doAssign = doAssign || (facesToProcess[j] && isThisCellOnAFace[j]);
         if(doAssign) {
            mpiGrid[dccrgId]->sysBoundaryFlag = this->getIndex();
         }         
      }
      
      // Assign boundary flags to local fsgrid cells
      const std::array<int, 3> gridDims(technicalGrid.getLocalSize());  
      for (int k=0; k<gridDims[2]; k++) {
         for (int j=0; j<gridDims[1]; j++) {
            for (int i=0; i<gridDims[0]; i++) {
               const auto& coords = technicalGrid.getPhysicalCoords(i,j,k);

               // Shift to the center of the fsgrid cell
               auto cellCenterCoords = coords;
               cellCenterCoords[0] += 0.5 * technicalGrid.DX;
               cellCenterCoords[1] += 0.5 * technicalGrid.DY;
               cellCenterCoords[2] += 0.5 * technicalGrid.DZ;
               const auto refLvl = mpiGrid.get_refinement_level(mpiGrid.get_existing_cell(cellCenterCoords));

               if(refLvl == -1) {
                  cerr << "Error, could not get refinement level of remote DCCRG cell " << __FILE__ << " " << __LINE__ << endl;
               }
               
               creal dx = P::dx_ini * pow(2,-refLvl);
               creal dy = P::dy_ini * pow(2,-refLvl);
               creal dz = P::dz_ini * pow(2,-refLvl);

               isThisCellOnAFace.fill(false);
               doAssign = false;

               determineFace(isThisCellOnAFace.data(), cellCenterCoords[0], cellCenterCoords[1], cellCenterCoords[2], dx, dy, dz);
               for(int iface=0; iface<6; iface++) doAssign = doAssign || (facesToProcess[iface] && isThisCellOnAFace[iface]);
               if(doAssign) {
                  technicalGrid.get(i,j,k)->sysBoundaryFlag = this->getIndex();
               }
            }
         }
      }
      
      return true;
   }
   
   bool Static::applyInitialState(
      const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2> & perBGrid,
      Project &project
   ) {
      const vector<CellID>& cells = getLocalCells();
      #pragma omp parallel for
      for (uint i=0; i<cells.size(); ++i) {
         SpatialCell* cell = mpiGrid[cells[i]];
         if (cell->sysBoundaryFlag != this->getIndex()) continue;
         
         if(!Parameters::isRestart) {
            // Defined in project.cpp, used here as the static cell has the same state 
            // as the initial state of non-system boundary cells.
            project.setCell(cell);
            // WARNING Time-independence assumed here.
            cell->parameters[CellParams::RHOM_DT2] = cell->parameters[CellParams::RHOM];
            cell->parameters[CellParams::RHOQ_DT2] = cell->parameters[CellParams::RHOQ];
            cell->parameters[CellParams::VX_DT2] = cell->parameters[CellParams::VX];
            cell->parameters[CellParams::VY_DT2] = cell->parameters[CellParams::VY];
            cell->parameters[CellParams::VZ_DT2] = cell->parameters[CellParams::VZ];
            cell->parameters[CellParams::P_11_DT2] = cell->parameters[CellParams::P_11];
            cell->parameters[CellParams::P_22_DT2] = cell->parameters[CellParams::P_22];
            cell->parameters[CellParams::P_33_DT2] = cell->parameters[CellParams::P_33];
         }
      }
      return true;
   }

   Real Static::fieldSolverBoundaryCondMagneticField(
      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2> & bGrid,
      FsGrid< fsgrids::technical, 2> & technicalGrid,
      cint i,
      cint j,
      cint k,
      creal& dt,
      cuint& component
   ) {
      return bGrid.get(i,j,k)->at(fsgrids::bfield::PERBX+component);
   }

   void Static::fieldSolverBoundaryCondElectricField(
      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 2> & EGrid,
      cint i,
      cint j,
      cint k,
      cuint component
   ) { 
     //EGrid.get(i,j,k)->at(fsgrids::efield::EX+component) = 0.0;
   }
   
   void Static::fieldSolverBoundaryCondHallElectricField(
      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, 2> & EHallGrid,
      cint i,
      cint j,
      cint k,
      cuint component
   ) { 
     /*
      std::array<Real, fsgrids::ehall::N_EHALL> * cp = EHallGrid.get(i,j,k);
      switch (component) {
         case 0:
            cp->at(fsgrids::ehall::EXHALL_000_100) = 0.0;
            cp->at(fsgrids::ehall::EXHALL_010_110) = 0.0;
            cp->at(fsgrids::ehall::EXHALL_001_101) = 0.0;
            cp->at(fsgrids::ehall::EXHALL_011_111) = 0.0;
            break;
         case 1:
            cp->at(fsgrids::ehall::EYHALL_000_010) = 0.0;
            cp->at(fsgrids::ehall::EYHALL_100_110) = 0.0;
            cp->at(fsgrids::ehall::EYHALL_001_011) = 0.0;
            cp->at(fsgrids::ehall::EYHALL_101_111) = 0.0;
            break;
         case 2:
            cp->at(fsgrids::ehall::EZHALL_000_001) = 0.0;
            cp->at(fsgrids::ehall::EZHALL_100_101) = 0.0;
            cp->at(fsgrids::ehall::EZHALL_010_011) = 0.0;
            cp->at(fsgrids::ehall::EZHALL_110_111) = 0.0;
            break;
         default:
            cerr << __FILE__ << ":" << __LINE__ << ":" << " Invalid component" << endl;
      }
     */
   }
   
   void Static::fieldSolverBoundaryCondGradPeElectricField(
      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, 2> & EGradPeGrid,
      cint i,
      cint j,
      cint k,
      cuint component
   ) {  
     //EGradPeGrid.get(i,j,k)->at(fsgrids::egradpe::EXGRADPE+component) = 0.0;
   }
   
   void Static::fieldSolverBoundaryCondDerivatives(
      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, 2> & dPerBGrid,
      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 2> & dMomentsGrid,
      cint i,
      cint j,
      cint k,
      cuint& RKCase,
      cuint& component
   ) {  
     //this->setCellDerivativesToZero(dPerBGrid, dMomentsGrid, i, j, k, component);
   }
   
   void Static::fieldSolverBoundaryCondBVOLDerivatives(
      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, 2> & volGrid,
      cint i,
      cint j,
      cint k,
      cuint& component
   ) {  
     //this->setCellBVOLDerivativesToZero(volGrid, i, j, k, component);
   }
   
   /**
    * NOTE that this is called once for each particle species!
    */
   void Static::vlasovBoundaryCondition(
      const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const CellID& cellID,
      const uint popID,
      const bool calculate_V_moments
   ) { 
     return;
   }

   void Static::getFaces(bool* faces) {
      for(uint i=0; i<6; i++) faces[i] = facesToProcess[i];
   }
      
   std::string Static::getName() const {return "Static";}
   uint Static::getIndex() const {return sysboundarytype::STATIC;}
      
}
