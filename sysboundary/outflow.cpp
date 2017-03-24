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

/*!\file outflow.cpp
 * \brief Implementation of the class SysBoundaryCondition::Outflow to handle cells classified as sysboundarytype::OUTFLOW.
 */

#include <cstdlib>
#include <iostream>

#include "../object_wrapper.h"
#include "outflow.h"
#include "../projects/projects_common.h"
#include "../fieldsolver/fs_common.h"
#include "../fieldsolver/ldz_magnetic_field.hpp"
#include "../vlasovmover.h"

#ifndef NDEBUG
   #define DEBUG_OUTFLOW
#endif
#ifdef DEBUG_SYSBOUNDARY
   #define DEBUG_OUTFLOW
#endif

using namespace std;

namespace SBC {
   Outflow::Outflow(): SysBoundaryCondition() { }
   Outflow::~Outflow() { }
   
   void Outflow::addParameters() {
      const std::string defStr = "Copy";
      Readparameters::addComposing("outflow.face", "List of faces on which outflow boundary conditions are to be applied ([xyz][+-]).");
      Readparameters::addComposing("outflow.faceNoFields", "List of faces on which no field outflow boundary conditions are to be applied ([xyz][+-]).");
      Readparameters::add("outflow.vlasovScheme_face_x+", "Scheme to use on the face x+ (Copy, Limit, None)", defStr);
      Readparameters::add("outflow.vlasovScheme_face_x-", "Scheme to use on the face x- (Copy, Limit, None)", defStr);
      Readparameters::add("outflow.vlasovScheme_face_y+", "Scheme to use on the face y+ (Copy, Limit, None)", defStr);
      Readparameters::add("outflow.vlasovScheme_face_y-", "Scheme to use on the face y- (Copy, Limit, None)", defStr);
      Readparameters::add("outflow.vlasovScheme_face_z+", "Scheme to use on the face z+ (Copy, Limit, None)", defStr);
      Readparameters::add("outflow.vlasovScheme_face_z-", "Scheme to use on the face z- (Copy, Limit, None)", defStr);
      Readparameters::add("outflow.precedence", "Precedence value of the outflow system boundary condition (integer), the higher the stronger.", 4);
      Readparameters::add("outflow.quench", "Factor by which to quench the inflowing parts of the velocity distribution function.", 1.0);
      Readparameters::add("outflow.reapplyUponRestart", "If 0 (default), keep going with the state existing in the restart file. If 1, calls again applyInitialState. Can be used to change boundary condition behaviour during a run.", 0);
   }
   
   void Outflow::getParameters() {
      std::array<std::string, 6> vlasovSysBoundarySchemeName;
      int myRank;
      MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
      if(!Readparameters::get("outflow.face", this->faceList)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!Readparameters::get("outflow.faceNoFields", this->faceNoFieldsList)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!Readparameters::get("outflow.vlasovScheme_face_x+", vlasovSysBoundarySchemeName[0])) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!Readparameters::get("outflow.vlasovScheme_face_x-", vlasovSysBoundarySchemeName[1])) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!Readparameters::get("outflow.vlasovScheme_face_y+", vlasovSysBoundarySchemeName[2])) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!Readparameters::get("outflow.vlasovScheme_face_y-", vlasovSysBoundarySchemeName[3])) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!Readparameters::get("outflow.vlasovScheme_face_z+", vlasovSysBoundarySchemeName[4])) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!Readparameters::get("outflow.vlasovScheme_face_z-", vlasovSysBoundarySchemeName[5])) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      for(uint i=0; i<6 ; i++) {
         if(vlasovSysBoundarySchemeName[i] == "None") {
            faceVlasovScheme[i] = vlasovscheme::NONE;
         } else if (vlasovSysBoundarySchemeName[i] == "Copy") {
            faceVlasovScheme[i] = vlasovscheme::COPY;
         } else if(vlasovSysBoundarySchemeName[i] == "Limit") {
            faceVlasovScheme[i] = vlasovscheme::LIMIT;
         } else {
            if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: " << vlasovSysBoundarySchemeName[i] << " is an invalid Outflow Vlasov scheme!" << endl;
            exit(1);
         }
      }
      if(!Readparameters::get("outflow.precedence", precedence)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!Readparameters::get("outflow.quench", this->quenchFactor)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      uint reapply;
      if(!Readparameters::get("outflow.reapplyUponRestart",reapply)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      };
      this->applyUponRestart = false;
      if(reapply == 1) {
         this->applyUponRestart = true;
      }
   }
   
   bool Outflow::initSysBoundary(
      creal& t,
      Project &project
   ) {
      /* The array of bool describes which of the x+, x-, y+, y-, z+, z- faces are to have outflow system boundary conditions.
       * A true indicates the corresponding face will have outflow.
       * The 6 elements correspond to x+, x-, y+, y-, z+, z- respectively.
       */
      for(uint i=0; i<6; i++) {
         facesToProcess[i] = false;
         facesToSkipFields[i] = false;
      }
      
      this->getParameters();
      
      isThisDynamic = false;
      
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
      for (it = faceNoFieldsList.begin();
           it != faceNoFieldsList.end();
      it++) {
         if(*it == "x+") facesToSkipFields[0] = true;
         if(*it == "x-") facesToSkipFields[1] = true;
         if(*it == "y+") facesToSkipFields[2] = true;
         if(*it == "y-") facesToSkipFields[3] = true;
         if(*it == "z+") facesToSkipFields[4] = true;
         if(*it == "z-") facesToSkipFields[5] = true;
      }
      return true;
   }
   
   bool Outflow::assignSysBoundary(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid) {
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
   
   bool Outflow::applyInitialState(
      const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      Project &project
   ) {
      const vector<CellID>& cells = getLocalCells();
      #pragma omp parallel for
      for (uint i=0; i<cells.size(); ++i) {
         SpatialCell* cell = mpiGrid[cells[i]];
         if (cell->sysBoundaryFlag != this->getIndex()) continue;
         
         // Defined in project.cpp, used here as the outflow cell has the same state 
         // as the initial state of non-system boundary cells.
         project.setCell(cell);
         // WARNING Time-independence assumed here.
         cell->parameters[CellParams::RHO_DT2] = cell->parameters[CellParams::RHO];
         cell->parameters[CellParams::RHOVX_DT2] = cell->parameters[CellParams::RHOVX];
         cell->parameters[CellParams::RHOVY_DT2] = cell->parameters[CellParams::RHOVY];
         cell->parameters[CellParams::RHOVZ_DT2] = cell->parameters[CellParams::RHOVZ];
      }

      return true;
   }

   Real Outflow::fieldSolverBoundaryCondMagneticField(
      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2> & perBGrid,
      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2> & perBDt2Grid,
      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 2> & EGrid,
      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 2> & EDt2Grid,
      FsGrid< fsgrids::technical, 2> & technicalGrid,
      cint i,
      cint j,
      cint k,
      creal& dt,
      cuint& RKCase,
      cuint& component
   ) {
      Real fieldValue = -1.0;
      
      creal dx =technicalGrid.DX;
      creal dy =technicalGrid.DY;
      creal dz =technicalGrid.DZ;
      const std::array<int, 3> globalIndices = technicalGrid.getGlobalIndices(i,j,k);
      creal x = (convert<Real>(globalIndices[0])+0.5)*dx + Parameters::xmin;
      creal y = (convert<Real>(globalIndices[1])+0.5)*dy + Parameters::ymin;
      creal z = (convert<Real>(globalIndices[2])+0.5)*dz + Parameters::zmin;
      
      bool isThisCellOnAFace[6];
      determineFace(&isThisCellOnAFace[0], x, y, z, dx, dy, dz, true);
      
      cuint sysBoundaryLayer = technicalGrid.get(i,j,k)->sysBoundaryLayer;
      
      if (sysBoundaryLayer == 1
         && isThisCellOnAFace[0]                            // we are on the face
         && this->facesToProcess[0]                         // we are supposed to do this face
         && !this->facesToSkipFields[0]                     // we are not supposed to skip fields on this face
         && component == 0                                  // we do the component normal to this face
         && !(isThisCellOnAFace[2] || isThisCellOnAFace[3] || isThisCellOnAFace[4] || isThisCellOnAFace[5]) // we are not in a corner
      ) {
         propagateMagneticField(perBGrid, perBDt2Grid, EGrid, EDt2Grid, i, j, k, dt, RKCase, true, false, false);
         if(RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
            fieldValue = perBGrid.get(i,j,k)->at(fsgrids::bfield::PERBX + component);
         } else {
            fieldValue = perBDt2Grid.get(i,j,k)->at(fsgrids::bfield::PERBX + component);
         }
      } else if (sysBoundaryLayer == 1
         && isThisCellOnAFace[2]                            // we are on the face
         && this->facesToProcess[2]                         // we are supposed to do this face
         && !this->facesToSkipFields[2]                     // we are not supposed to skip fields on this face
         && component == 1                                  // we do the component normal to this face
         && !(isThisCellOnAFace[0] || isThisCellOnAFace[1] || isThisCellOnAFace[4] || isThisCellOnAFace[5]) // we are not in a corner
      ) {
         propagateMagneticField(perBGrid, perBDt2Grid, EGrid, EDt2Grid, i, j, k, dt, RKCase, false, true, false);
         if(RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
            fieldValue = perBGrid.get(i,j,k)->at(fsgrids::bfield::PERBX + component);
         } else {
            fieldValue = perBDt2Grid.get(i,j,k)->at(fsgrids::bfield::PERBX + component);
         }
      } else if (sysBoundaryLayer == 1
         && isThisCellOnAFace[4]                            // we are on the face
         && this->facesToProcess[4]                         // we are supposed to do this face
         && !this->facesToSkipFields[4]                     // we are not supposed to skip fields on this face
         && component == 2                                  // we do the component normal to this face
         && !(isThisCellOnAFace[0] || isThisCellOnAFace[1] || isThisCellOnAFace[2] || isThisCellOnAFace[3]) // we are not in a corner
      ) {
         propagateMagneticField(perBGrid, perBDt2Grid, EGrid, EDt2Grid, i, j, k, dt, RKCase, false, false, true);
         if(RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
            fieldValue = perBGrid.get(i,j,k)->at(fsgrids::bfield::PERBX + component);
         } else {
            fieldValue = perBDt2Grid.get(i,j,k)->at(fsgrids::bfield::PERBX + component);
         }
      } else {
         if(RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
            fieldValue = fieldBoundaryCopyFromExistingFaceNbrMagneticField(perBGrid, technicalGrid, i, j, k, component);
         } else {
            fieldValue = fieldBoundaryCopyFromExistingFaceNbrMagneticField(perBDt2Grid, technicalGrid, i, j, k, component);
         }
      }
      
      return fieldValue;
   }

   void Outflow::fieldSolverBoundaryCondElectricField(
      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 2> & EGrid,
      cint i,
      cint j,
      cint k,
      cuint component
   ) {
      EGrid.get(i,j,k)->at(fsgrids::efield::EX+component) = 0.0;
   }
   
   void Outflow::fieldSolverBoundaryCondHallElectricField(
      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, 2> & EHallGrid,
      cint i,
      cint j,
      cint k,
      cuint component
   ) {
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
   }
   
   void Outflow::fieldSolverBoundaryCondGradPeElectricField(
      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, 2> & EGradPeGrid,
      cint i,
      cint j,
      cint k,
      cuint component
   ) {
      EGradPeGrid.get(i,j,k)->at(fsgrids::egradpe::EXGRADPE+component) = 0.0;
   }
   
   void Outflow::fieldSolverBoundaryCondDerivatives(
      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, 2> & dPerBGrid,
      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 2> & dMomentsGrid,
      cint i,
      cint j,
      cint k,
      cuint& RKCase,
      cuint& component
   ) {
      this->setCellDerivativesToZero(dPerBGrid, dMomentsGrid, i, j, k, component);
   }
   
   void Outflow::fieldSolverBoundaryCondBVOLDerivatives(
      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, 2> & volGrid,
      cint i,
      cint j,
      cint k,
      cuint& component
   ) {
      this->setCellBVOLDerivativesToZero(volGrid, i, j, k, component);
   }
   
   /**
    * NOTE that this is called once for each particle species!
    * @param mpiGrid
    * @param cellID
    */
   void Outflow::vlasovBoundaryCondition(
      const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const CellID& cellID,
      const int& popID
   ) {
//      phiprof::start("vlasovBoundaryCondition (Outflow)");
      
      SpatialCell* cell = mpiGrid[cellID];
      creal* const cellParams = cell->parameters;
      creal dx = cellParams[CellParams::DX];
      creal dy = cellParams[CellParams::DY];
      creal dz = cellParams[CellParams::DZ];
      creal x = cellParams[CellParams::XCRD] + 0.5*dx;
      creal y = cellParams[CellParams::YCRD] + 0.5*dy;
      creal z = cellParams[CellParams::ZCRD] + 0.5*dz;
      
      bool isThisCellOnAFace[6];
      determineFace(&isThisCellOnAFace[0], x, y, z, dx, dy, dz, true);
      
      for(uint i=0; i<6; i++) {
         if(isThisCellOnAFace[i] && facesToProcess[i]) {
            switch(this->faceVlasovScheme[i]) {
               case vlasovscheme::NONE:
                  break;
               case vlasovscheme::COPY:
                  if (cell->sysBoundaryLayer == 1) {
                     vlasovBoundaryCopyFromTheClosestNbr(mpiGrid,cellID,false,popID);
                  } else {
                     vlasovBoundaryCopyFromTheClosestNbr(mpiGrid,cellID,true,popID);
                  }
                  break;
               case vlasovscheme::LIMIT:
                  if (cell->sysBoundaryLayer == 1) {
                     vlasovBoundaryCopyFromTheClosestNbrAndLimit(mpiGrid,cellID,popID);
                  } else {
                     vlasovBoundaryCopyFromTheClosestNbr(mpiGrid,cellID,true,popID);
                  }
                  break;
               default:
                  std::cerr << __FILE__ << ":" << __LINE__ << "ERROR: invalid Outflow Vlasov scheme!" << std::endl;
                  exit(1);
                  break;
            }
         }
      }
      
//      phiprof::stop("vlasovBoundaryCondition (Outflow)");
   }
   
   void Outflow::getFaces(bool* faces) {
      for(uint i=0; i<6; i++) faces[i] = facesToProcess[i];
   }
   
   std::string Outflow::getName() const {return "Outflow";}
   uint Outflow::getIndex() const {return sysboundarytype::OUTFLOW;}
      
}
