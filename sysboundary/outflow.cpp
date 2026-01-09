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
#include "../vlasovsolver/vlasovmover.h"
#include "../grid.h"

#ifdef DEBUG_VLASIATOR
   #define DEBUG_OUTFLOW
#endif
#ifdef DEBUG_SYSBOUNDARY
   #define DEBUG_OUTFLOW
#endif

using namespace std;

namespace SBC {
   Outflow::Outflow(): OuterBoundaryCondition() { }
   Outflow::~Outflow() { }

   void Outflow::addParameters() {
      const string defStr = "Copy";
      Readparameters::addComposing("outflow.faceNoFields", "List of faces on which no field outflow boundary conditions are to be applied ([xyz][+-]).");
      Readparameters::add("outflow.precedence", "Precedence value of the outflow system boundary condition (integer), the higher the stronger.", 4);
      Readparameters::add("outflow.reapplyUponRestart", "If 0 (default), keep going with the state existing in the restart file. If 1, calls again applyInitialState. Can be used to change boundary condition behaviour during a run.", 0);

      // Per-population parameters
      for(uint i=0; i< getObjectWrapper().particleSpecies.size(); i++) {
        const string& pop = getObjectWrapper().particleSpecies[i].name;

        Readparameters::addComposing(pop + "_outflow.reapplyFaceUponRestart", "List of faces on which outflow boundary conditions are to be reapplied upon restart ([xyz][+-]).");
        Readparameters::addComposing(pop + "_outflow.face", "List of faces on which outflow boundary conditions are to be applied ([xyz][+-]).");
        Readparameters::add(pop + "_outflow.vlasovScheme_face_x+", "Scheme to use on the face x+ (Copy, None)", defStr);
        Readparameters::add(pop + "_outflow.vlasovScheme_face_x-", "Scheme to use on the face x- (Copy, None)", defStr);
        Readparameters::add(pop + "_outflow.vlasovScheme_face_y+", "Scheme to use on the face y+ (Copy, None)", defStr);
        Readparameters::add(pop + "_outflow.vlasovScheme_face_y-", "Scheme to use on the face y- (Copy, None)", defStr);
        Readparameters::add(pop + "_outflow.vlasovScheme_face_z+", "Scheme to use on the face z+ (Copy, None)", defStr);
        Readparameters::add(pop + "_outflow.vlasovScheme_face_z-", "Scheme to use on the face z- (Copy, None)", defStr);

        Readparameters::add(pop + "_outflow.quench", "Factor by which to quench the inflowing parts of the velocity distribution function.", 1.0);
      }
   }

   void Outflow::getParameters() {
      int myRank;
      MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
      Readparameters::get("outflow.faceNoFields", this->faceNoFieldsList);
      Readparameters::get("outflow.precedence", precedence);
      uint reapply;
      Readparameters::get("outflow.reapplyUponRestart", reapply);
      this->applyUponRestart = false;
      if(reapply == 1) {
         this->applyUponRestart = true;
      }

      // Per-species parameters
      for(uint i=0; i< getObjectWrapper().particleSpecies.size(); i++) {
        const string& pop = getObjectWrapper().particleSpecies[i].name;
        OutflowSpeciesParameters sP;

        // Unless we find out otherwise, we assume that this species will not be treated at any boundary
        for(int j=0; j<6; j++) {
          sP.facesToSkipVlasov[j] = true;
        }

        vector<string> thisSpeciesFaceList;
        Readparameters::get(pop + "_outflow.face", thisSpeciesFaceList);

        for(auto& face : thisSpeciesFaceList) {
          if(face == "x+") { facesToProcess[0] = true; sP.facesToSkipVlasov[0] = false; }
          if(face == "x-") { facesToProcess[1] = true; sP.facesToSkipVlasov[1] = false; }
          if(face == "y+") { facesToProcess[2] = true; sP.facesToSkipVlasov[2] = false; }
          if(face == "y-") { facesToProcess[3] = true; sP.facesToSkipVlasov[3] = false; }
          if(face == "z+") { facesToProcess[4] = true; sP.facesToSkipVlasov[4] = false; }
          if(face == "z-") { facesToProcess[5] = true; sP.facesToSkipVlasov[5] = false; }
        }

        Readparameters::get(pop + "_outflow.reapplyFaceUponRestart", sP.faceToReapplyUponRestartList);
        array<string, 6> vlasovSysBoundarySchemeName;
        Readparameters::get(pop + "_outflow.vlasovScheme_face_x+", vlasovSysBoundarySchemeName[0]);
        Readparameters::get(pop + "_outflow.vlasovScheme_face_x-", vlasovSysBoundarySchemeName[1]);
        Readparameters::get(pop + "_outflow.vlasovScheme_face_y+", vlasovSysBoundarySchemeName[2]);

        Readparameters::get(pop + "_outflow.vlasovScheme_face_y-", vlasovSysBoundarySchemeName[3]);
        Readparameters::get(pop + "_outflow.vlasovScheme_face_z+", vlasovSysBoundarySchemeName[4]);
        Readparameters::get(pop + "_outflow.vlasovScheme_face_z-", vlasovSysBoundarySchemeName[5]);
        for(uint j=0; j<6 ; j++) {
           if(vlasovSysBoundarySchemeName[j] == "None") {
              sP.faceVlasovScheme[j] = vlasovscheme::NONE;
           } else if (vlasovSysBoundarySchemeName[j] == "Copy") {
              sP.faceVlasovScheme[j] = vlasovscheme::COPY;
           } else {
              abort_mpi("ERROR: " + vlasovSysBoundarySchemeName[j] + " is an invalid Outflow Vlasov scheme!");
           }
        }

        Readparameters::get(pop + "_outflow.quench", sP.quenchFactor);

        speciesParams.push_back(sP);
      }
   }

   void Outflow::initSysBoundary(
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
         facesToReapply[i] = false;
      }

      this->getParameters();

      dynamic = false;

      vector<string>::const_iterator it;
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

      for(uint i=0; i< getObjectWrapper().particleSpecies.size(); i++) {
         OutflowSpeciesParameters& sP = this->speciesParams[i];
         for (it = sP.faceToReapplyUponRestartList.begin();
              it != sP.faceToReapplyUponRestartList.end();
         it++) {
            if(*it == "x+") facesToReapply[0] = true;
            if(*it == "x-") facesToReapply[1] = true;
            if(*it == "y+") facesToReapply[2] = true;
            if(*it == "y-") facesToReapply[3] = true;
            if(*it == "z+") facesToReapply[4] = true;
            if(*it == "z-") facesToReapply[5] = true;
         }
      }
   }


   void Outflow::assignSysBoundary(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid,
                                   FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid) {

      bool doAssign;
      array<bool,6> isThisCellOnAFace;

      // Assign boundary flags to local DCCRG cells
      const vector<CellID>& cells = getLocalCells();
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
         for(int j=0; j<6; j++) {
            doAssign = doAssign || (facesToProcess[j] && isThisCellOnAFace[j]);
         }
         if(doAssign) {
            mpiGrid[dccrgId]->sysBoundaryFlag = this->getIndex();
         }
      }

      // Assign boundary flags to local fsgrid cells
      const array<FsGridTools::FsIndex_t, 3> gridDims(technicalGrid.getLocalSize());
      for (FsGridTools::FsIndex_t k=0; k<gridDims[2]; k++) {
         for (FsGridTools::FsIndex_t j=0; j<gridDims[1]; j++) {
            for (FsGridTools::FsIndex_t i=0; i<gridDims[0]; i++) {
               const auto& coords = technicalGrid.getPhysicalCoords(i,j,k);

               // Shift to the center of the fsgrid cell
               auto cellCenterCoords = coords;
               cellCenterCoords[0] += 0.5 * technicalGrid.DX;
               cellCenterCoords[1] += 0.5 * technicalGrid.DY;
               cellCenterCoords[2] += 0.5 * technicalGrid.DZ;
               const auto refLvl = mpiGrid.get_refinement_level(mpiGrid.get_existing_cell(cellCenterCoords));

               if (refLvl == -1) {
                  abort_mpi("ERROR: Could not get refinement level of remote DCCRG cell!", 1);
               }

               creal dx = P::dx_ini / pow(2, refLvl);
               creal dy = P::dy_ini / pow(2, refLvl);
               creal dz = P::dz_ini / pow(2, refLvl);

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
   }

   void Outflow::applyInitialState(
      dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
      FsGrid< array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
      FsGrid<std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH>& BgBGrid,
      Project &project
   ) {
      const vector<CellID>& cells = getLocalCells();
      #pragma omp parallel for schedule(static)
      for (uint i=0; i<cells.size(); ++i) {
         CellID id = cells[i];
         SpatialCell* cell = mpiGrid[id];
         if (cell->sysBoundaryFlag != this->getIndex()) {
            continue;
         }

         bool doApply = true;

         if(Parameters::isRestart) {
            std::array<bool, 6> isThisCellOnAFace;
            determineFace(isThisCellOnAFace, mpiGrid, id);

            doApply=false;
            // Comparison of the array defining which faces to use and the array telling on which faces this cell is
            for (uint j=0; j<6; j++) {
               doApply = doApply || (facesToReapply[j] && isThisCellOnAFace[j]);
            }
         }

         if (doApply) {
            project.setCell(cell); // We set everything including VDF even in L2 cells to avoid a pile of spaghetti. Won't get communicated.
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
   }

   void Outflow::updateState(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid,
                             FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
                             FsGrid<array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH>& perBGrid,
                             FsGrid<std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH>& BgBGrid,
                             creal t) {}

   Real Outflow::fieldSolverBoundaryCondMagneticField(
      FsGrid< array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & bGrid,
      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & bgbGrid,
      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
      cint i,
      cint j,
      cint k,
      creal dt,
      cuint component
   ) {
      switch(component) {
      case 0:
         return fieldBoundaryCopyFromSolvingNbrMagneticField(bGrid, technicalGrid, i, j, k, component, compute::BX);
      case 1:
         return fieldBoundaryCopyFromSolvingNbrMagneticField(bGrid, technicalGrid, i, j, k, component, compute::BY);
      case 2:
         return fieldBoundaryCopyFromSolvingNbrMagneticField(bGrid, technicalGrid, i, j, k, component, compute::BZ);
      default:
         return 0.0;
      }
   }

   void Outflow::fieldSolverBoundaryCondElectricField(
      FsGrid< array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
      cint i,
      cint j,
      cint k,
      cuint component
   ) {
      EGrid.get(i,j,k)->at(fsgrids::efield::EX+component) = 0.0;
   }

   void Outflow::fieldSolverBoundaryCondHallElectricField(
      FsGrid< array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
      cint i,
      cint j,
      cint k,
      cuint component
   ) {
      array<Real, fsgrids::ehall::N_EHALL> * cp = EHallGrid.get(i,j,k);
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
      FsGrid< array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
      cint i,
      cint j,
      cint k,
      cuint component
   ) {
      EGradPeGrid.get(i,j,k)->at(fsgrids::egradpe::EXGRADPE+component) = 0.0;
   }

   void Outflow::fieldSolverBoundaryCondDerivatives(
      FsGrid< array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
      FsGrid< array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
      cint i,
      cint j,
      cint k,
      cuint RKCase,
      cuint component
   ) {
      this->setCellDerivativesToZero(dPerBGrid, dMomentsGrid, i, j, k, component);
   }

   void Outflow::fieldSolverBoundaryCondBVOLDerivatives(
      FsGrid< array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
      cint i,
      cint j,
      cint k,
      cuint component
   ) {
      this->setCellBVOLDerivativesToZero(volGrid, i, j, k, component);
   }

   /**
    * NOTE that this is called once for each particle species!
    * @param mpiGrid
    * @param cellID
    */
   void Outflow::vlasovBoundaryCondition(
      dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const CellID& cellID,
      const uint popID,
      const bool calculate_V_moments
   ) {

      const OutflowSpeciesParameters& sP = this->speciesParams[popID];
      if (mpiGrid[cellID]->sysBoundaryFlag != this->getIndex()) {
         return;
      }

      std::array<bool, 6> isThisCellOnAFace;
      determineFace(isThisCellOnAFace, mpiGrid, cellID, true);

      for(uint i=0; i<6; i++) {
         if(isThisCellOnAFace[i] && facesToProcess[i] && !sP.facesToSkipVlasov[i]) {
            switch(sP.faceVlasovScheme[i]) {
               case vlasovscheme::NONE:
                  break;
               case vlasovscheme::COPY:
                  if(mpiGrid[cellID]->sysBoundaryLayer == 1) {
                     vlasovBoundaryCopyFromTheClosestNbr(mpiGrid,cellID,false,popID,calculate_V_moments); // copies VDF too
                  } else {
                     vlasovBoundaryCopyFromTheClosestNbr(mpiGrid,cellID,true,popID,calculate_V_moments); // no VDF copy
                  }
                  break;
               default:
                  abort_mpi("ERROR: invalid Outflow Vlasov scheme", 1);
            }
         }
      }
   }

   void Outflow::setupL2OutflowAtRestart(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid) {
      // These updates emulated from SysBoundary::applySysBoundaryVlasovConditions()
      // Needs a call to updateRemoteVelocityBlockLists(), done in the SysBoundary class.
      SpatialCell::set_mpi_transfer_type(Transfer::CELL_PARAMETERS | Transfer::POP_METADATA | Transfer::CELL_SYSBOUNDARYFLAG, true);
      mpiGrid.update_copies_of_remote_neighbors(Neighborhoods::SYSBOUNDARIES_EXTENDED);

      for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
         SpatialCell::setCommunicatedSpecies(popID);
         updateRemoteVelocityBlockLists(mpiGrid, popID, Neighborhoods::SYSBOUNDARIES);
         SpatialCell::set_mpi_transfer_type(Transfer::VEL_BLOCK_DATA, true);
         mpiGrid.update_copies_of_remote_neighbors(Neighborhoods::SYSBOUNDARIES);
      }

      const vector<CellID>& cells = getLocalCells();
      #pragma omp parallel
      {
         #pragma omp for schedule(guided,1)
         for(uint i=0; i<cells.size(); i++) {
            const CellID cellID = cells[i];
            // As of 20250505 the loop only does something for layer 1 in COPY mode so the check for layer 1 was moved here for earlier loop continuation.
            if(mpiGrid[cellID]->sysBoundaryFlag != this->getIndex() || mpiGrid[cellID]->sysBoundaryLayer != 1) {
               continue;
            }

            std::array<bool, 6> isThisCellOnAFace;
            determineFace(isThisCellOnAFace, mpiGrid, cellID, true);

            for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
               const OutflowSpeciesParameters& sP = this->speciesParams[popID];
               for(uint i=0; i<6; i++) {
                  if(isThisCellOnAFace[i] && facesToProcess[i] && !sP.facesToSkipVlasov[i]) {
                     switch(sP.faceVlasovScheme[i]) {
                        case vlasovscheme::NONE:
                           break;
                        case vlasovscheme::COPY:
                           // if(mpiGrid[cellID]->sysBoundaryLayer == 1) { // This is actually now (20250505) moved up a few lines for earlier loop continuation. Reinstate if other cases change!
                              vlasovBoundaryCopyFromTheClosestNbr(mpiGrid,cellID,false,popID,true); // first false means copy VDF too, second true means V moments
                              vlasovBoundaryCopyFromTheClosestNbr(mpiGrid,cellID,false,popID,false); // first false means copy VDF too, second false means R moments
                           // } // see comment above
                           break;
                        default:
                           abort_mpi("ERROR: invalid Outflow Vlasov scheme", 1);
                     } // switch
                  } // if on face
               } // faces
            } // populations
         } // cells

         #pragma omp barrier
         #pragma omp single
         {
            SpatialCell::set_mpi_transfer_type(Transfer::ALL_SPATIAL_DATA); // No need to update VDFs, we only copy the VDFs from L1 to L2 below.
            mpiGrid.update_copies_of_remote_neighbors(Neighborhoods::SYSBOUNDARIES);
         }
         #pragma omp barrier // maybe useless

         // then 2nd pass and copy from the closest L1 outflow neighbor
         #pragma omp for schedule(guided,1)
         for(uint i=0; i<cells.size(); i++) {
            const CellID cellID = cells[i];
            // As of 20250505 the loop only does something for layer 1 in COPY mode so the check for layer 2 was moved here for earlier loop continuation.
            if(mpiGrid[cellID]->sysBoundaryFlag != this->getIndex() || mpiGrid[cellID]->sysBoundaryLayer != 2) {
               continue;
            }

            std::array<bool, 6> isThisCellOnAFace;
            determineFace(isThisCellOnAFace, mpiGrid, cellID, true);

            for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
               const OutflowSpeciesParameters& sP = this->speciesParams[popID];
               for(uint i=0; i<6; i++) {
                  if(isThisCellOnAFace[i] && facesToProcess[i] && !sP.facesToSkipVlasov[i]) {
                     switch(sP.faceVlasovScheme[i]) {
                        case vlasovscheme::NONE:
                           break;
                        case vlasovscheme::COPY:
                           // if(mpiGrid[cellID]->sysBoundaryLayer == 2) { // This is actually now (20250505) moved up a few lines for earlier loop continuation. Reinstate if other cases change!
                              vlasovBoundaryCopyFromTheClosestL1OutflowNbr(mpiGrid,cellID,true,popID,true); // first true means copy moments only, second true means V moments
                              vlasovBoundaryCopyFromTheClosestL1OutflowNbr(mpiGrid,cellID,true,popID,false); // first true means copy moments only, second false means R moments
                           // } // see comment above
                           break;
                        default:
                           abort_mpi("ERROR: invalid Outflow Vlasov scheme", 1);
                     } // switch
                  } // if on face
               } // faces
            } // populations
         } // cells
      } // pragma omp parallel
   }

   void Outflow::getFaces(bool* faces) {
      for(uint i=0; i<6; i++) {
         faces[i] = facesToProcess[i];
      }
   }

   string Outflow::getName() const {return "Outflow";}
   uint Outflow::getIndex() const {return sysboundarytype::OUTFLOW;}

}
