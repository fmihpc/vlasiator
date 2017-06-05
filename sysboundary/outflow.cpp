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
      Readparameters::addComposing("outflow.faceNoFields", "List of faces on which no field outflow boundary conditions are to be applied ([xyz][+-]).");
      Readparameters::add("outflow.precedence", "Precedence value of the outflow system boundary condition (integer), the higher the stronger.", 4);
      Readparameters::add("outflow.reapplyUponRestart", "If 0 (default), keep going with the state existing in the restart file. If 1, calls again applyInitialState. Can be used to change boundary condition behaviour during a run.", 0);

      // Per-population parameters
      for(uint i=0; i< getObjectWrapper().particleSpecies.size(); i++) {
        const std::string& pop = getObjectWrapper().particleSpecies[i].name;

        Readparameters::addComposing(pop + "_outflow.face", "List of faces on which outflow boundary conditions are to be applied ([xyz][+-]).");
        Readparameters::add(pop + "_outflow.vlasovScheme_face_x+", "Scheme to use on the face x+ (Copy, Limit, None)", defStr);
        Readparameters::add(pop + "_outflow.vlasovScheme_face_x-", "Scheme to use on the face x- (Copy, Limit, None)", defStr);
        Readparameters::add(pop + "_outflow.vlasovScheme_face_y+", "Scheme to use on the face y+ (Copy, Limit, None)", defStr);
        Readparameters::add(pop + "_outflow.vlasovScheme_face_y-", "Scheme to use on the face y- (Copy, Limit, None)", defStr);
        Readparameters::add(pop + "_outflow.vlasovScheme_face_z+", "Scheme to use on the face z+ (Copy, Limit, None)", defStr);
        Readparameters::add(pop + "_outflow.vlasovScheme_face_z-", "Scheme to use on the face z- (Copy, Limit, None)", defStr);

        Readparameters::add(pop + "_outflow.quench", "Factor by which to quench the inflowing parts of the velocity distribution function.", 1.0);
      }
   }
   
   void Outflow::getParameters() {
      int myRank;
      MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
      if(!Readparameters::get("outflow.faceNoFields", this->faceNoFieldsList)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!Readparameters::get("outflow.precedence", precedence)) {
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

      // Per-species parameters
      for(uint i=0; i< getObjectWrapper().particleSpecies.size(); i++) {
        const std::string& pop = getObjectWrapper().particleSpecies[i].name;
        OutflowSpeciesParameters sP;

        // Unless we find out otherwise, we assume that this species will not be treated at any boundary
        for(int i=0; i<6; i++) {
          sP.facesToSkipVlasov[i] = true;
        }

        std::vector<std::string> thisSpeciesFaceList;
        if(!Readparameters::get(pop + "_outflow.face", thisSpeciesFaceList)) {
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

        std::array<std::string, 6> vlasovSysBoundarySchemeName;
        if(!Readparameters::get(pop + "_outflow.vlasovScheme_face_x+", vlasovSysBoundarySchemeName[0])) {
           if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added for population " << pop << "!" << endl;
           exit(1);
        }
        if(!Readparameters::get(pop + "_outflow.vlasovScheme_face_x-", vlasovSysBoundarySchemeName[1])) {
           if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added for population " << pop << "!" << endl;
           exit(1);
        }
        if(!Readparameters::get(pop + "_outflow.vlasovScheme_face_y+", vlasovSysBoundarySchemeName[2])) {
           if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added for population " << pop << "!" << endl;
           exit(1);
        }
        if(!Readparameters::get(pop + "_outflow.vlasovScheme_face_y-", vlasovSysBoundarySchemeName[3])) {
           if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added for population " << pop << "!" << endl;
           exit(1);
        }
        if(!Readparameters::get(pop + "_outflow.vlasovScheme_face_z+", vlasovSysBoundarySchemeName[4])) {
           if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added for population " << pop << "!" << endl;
           exit(1);
        }
        if(!Readparameters::get(pop + "_outflow.vlasovScheme_face_z-", vlasovSysBoundarySchemeName[5])) {
           if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added for population " << pop << "!" << endl;
           exit(1);
        }
        for(uint i=0; i<6 ; i++) {
           if(vlasovSysBoundarySchemeName[i] == "None") {
              sP.faceVlasovScheme[i] = vlasovscheme::NONE;
           } else if (vlasovSysBoundarySchemeName[i] == "Copy") {
              sP.faceVlasovScheme[i] = vlasovscheme::COPY;
           } else if(vlasovSysBoundarySchemeName[i] == "Limit") {
              sP.faceVlasovScheme[i] = vlasovscheme::LIMIT;
           } else {
              if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: " << vlasovSysBoundarySchemeName[i] << " is an invalid Outflow Vlasov scheme!" << endl;
              exit(1);
           }
        }


        if(!Readparameters::get(pop + "_outflow.quench", sP.quenchFactor)) {
           if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
           exit(1);
        }

        speciesParams.push_back(sP);
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
         cell->parameters[CellParams::RHOM_DT2] = cell->parameters[CellParams::RHOM];
         cell->parameters[CellParams::RHOMVX_DT2] = cell->parameters[CellParams::RHOMVX];
         cell->parameters[CellParams::RHOMVY_DT2] = cell->parameters[CellParams::RHOMVY];
         cell->parameters[CellParams::RHOMVZ_DT2] = cell->parameters[CellParams::RHOMVZ];
         cell->parameters[CellParams::RHOQ_DT2] = cell->parameters[CellParams::RHOQ];
      }

      return true;
   }

   Real Outflow::fieldSolverBoundaryCondMagneticField(
      const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const std::vector<fs_cache::CellCache>& cellCache,
      const uint16_t& localID,
      creal& dt,
      cuint& RKCase,
      cint& offset,
      cuint& component
   ) {
      const CellID cellID = cellCache[localID].cellID;
      Real fieldValue = -1.0;
      
      creal* const cellParams = cellCache[localID].cells[fs_cache::calculateNbrID(1  ,1  ,1  )]->parameters;
      creal dx = cellParams[CellParams::DX];
      creal dy = cellParams[CellParams::DY];
      creal dz = cellParams[CellParams::DZ];
      creal x = cellParams[CellParams::XCRD] + 0.5*dx;
      creal y = cellParams[CellParams::YCRD] + 0.5*dy;
      creal z = cellParams[CellParams::ZCRD] + 0.5*dz;
      
      bool isThisCellOnAFace[6];
      determineFace(&isThisCellOnAFace[0], x, y, z, dx, dy, dz, true);
      
      cuint sysBoundaryLayer = cellCache[localID].cells[fs_cache::calculateNbrID(1  ,1  ,1  )]->sysBoundaryLayer;
      
      if (sysBoundaryLayer == 1
         && isThisCellOnAFace[0]                            // we are on the face
         && this->facesToProcess[0]                         // we are supposed to do this face
         && !this->facesToSkipFields[0]                     // we are not supposed to skip fields on this face
         && component == 0                                  // we do the component normal to this face
         && !(isThisCellOnAFace[2] || isThisCellOnAFace[3] || isThisCellOnAFace[4] || isThisCellOnAFace[5]) // we are not in a corner
      ) {
         const vector<uint16_t> cellVector = {{localID}};
         propagateMagneticField(cellCache, cellVector, dt, RKCase, true, false, false);
         fieldValue = cellParams[CellParams::PERBX + component + offset];
      } else if (sysBoundaryLayer == 1
         && isThisCellOnAFace[2]                            // we are on the face
         && this->facesToProcess[2]                         // we are supposed to do this face
         && !this->facesToSkipFields[2]                     // we are not supposed to skip fields on this face
         && component == 1                                  // we do the component normal to this face
         && !(isThisCellOnAFace[0] || isThisCellOnAFace[1] || isThisCellOnAFace[4] || isThisCellOnAFace[5]) // we are not in a corner
      ) {
         const vector<uint16_t> cellVector = {{localID}};
         propagateMagneticField(cellCache, cellVector, dt, RKCase, false, true, false);
         fieldValue = cellParams[CellParams::PERBX + component + offset];
      } else if (sysBoundaryLayer == 1
         && isThisCellOnAFace[4]                            // we are on the face
         && this->facesToProcess[4]                         // we are supposed to do this face
         && !this->facesToSkipFields[4]                     // we are not supposed to skip fields on this face
         && component == 2                                  // we do the component normal to this face
         && !(isThisCellOnAFace[0] || isThisCellOnAFace[1] || isThisCellOnAFace[2] || isThisCellOnAFace[3]) // we are not in a corner
      ) {
         const vector<uint16_t> cellVector = {{localID}};
         propagateMagneticField(cellCache, cellVector, dt, RKCase, false, false, true);
         fieldValue = cellParams[CellParams::PERBX + component + offset];
      } else {
         bool skipThisOne = false;
         for(uint i=0; i<6; i++) {
            if(isThisCellOnAFace[i] && this->facesToProcess[i] && this->facesToSkipFields[i]) {
               skipThisOne = true;
               break;
            }
         }
         if(skipThisOne) {
            fieldValue = cellParams[CellParams::PERBX + component + offset]; // copy the existing value, we do not touch it
         } else {
            fieldValue = fieldBoundaryCopyFromExistingFaceNbrMagneticField(mpiGrid, cellID, component + offset);
         }
      }
      
      return fieldValue;
   }

   void Outflow::fieldSolverBoundaryCondElectricField(
      dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const CellID& cellID,
      cuint RKCase,
      cuint component
   ) {
      creal* const cellParams = mpiGrid[cellID]->parameters;
      creal dx = cellParams[CellParams::DX];
      creal dy = cellParams[CellParams::DY];
      creal dz = cellParams[CellParams::DZ];
      creal x = cellParams[CellParams::XCRD] + 0.5*dx;
      creal y = cellParams[CellParams::YCRD] + 0.5*dy;
      creal z = cellParams[CellParams::ZCRD] + 0.5*dz;
      
      bool isThisCellOnAFace[6];
      determineFace(&isThisCellOnAFace[0], x, y, z, dx, dy, dz, true);
      
      for(uint i=0; i<6; i++) {
         if(isThisCellOnAFace[i] && this->facesToSkipFields[i]) {
            return;
         }
      }
      
      if((RKCase == RK_ORDER1) || (RKCase == RK_ORDER2_STEP2)) {
         mpiGrid[cellID]->parameters[CellParams::EX+component] = 0.0;
      } else {// RKCase == RK_ORDER2_STEP1
         mpiGrid[cellID]->parameters[CellParams::EX_DT2+component] = 0.0;
      }
   }
   
   void Outflow::fieldSolverBoundaryCondHallElectricField(
      fs_cache::CellCache& cache,
      cuint RKCase,
      cuint component
   ) {

      Real* cp = cache.cells[fs_cache::calculateNbrID(1,1,1)]->parameters;
      
      creal dx = cp[CellParams::DX];
      creal dy = cp[CellParams::DY];
      creal dz = cp[CellParams::DZ];
      creal x = cp[CellParams::XCRD] + 0.5*dx;
      creal y = cp[CellParams::YCRD] + 0.5*dy;
      creal z = cp[CellParams::ZCRD] + 0.5*dz;
      
      bool isThisCellOnAFace[6];
      determineFace(&isThisCellOnAFace[0], x, y, z, dx, dy, dz, true);
      
      for(uint i=0; i<6; i++) {
         if(isThisCellOnAFace[i] && this->facesToSkipFields[i]) {
            return;
         }
      }
      
      switch (component) {
         case 0:
            cp[CellParams::EXHALL_000_100] = 0.0;
            cp[CellParams::EXHALL_010_110] = 0.0;
            cp[CellParams::EXHALL_001_101] = 0.0;
            cp[CellParams::EXHALL_011_111] = 0.0;
            break;
         case 1:
            cp[CellParams::EYHALL_000_010] = 0.0;
            cp[CellParams::EYHALL_100_110] = 0.0;
            cp[CellParams::EYHALL_001_011] = 0.0;
            cp[CellParams::EYHALL_101_111] = 0.0;
            break;
         case 2:
            cp[CellParams::EZHALL_000_001] = 0.0;
            cp[CellParams::EZHALL_100_101] = 0.0;
            cp[CellParams::EZHALL_010_011] = 0.0;
            cp[CellParams::EZHALL_110_111] = 0.0;
            break;
         default:
            cerr << __FILE__ << ":" << __LINE__ << ":" << " Invalid component" << endl;
      }
   }
   
   void Outflow::fieldSolverBoundaryCondGradPeElectricField(
      fs_cache::CellCache& cache,
      cuint RKCase,
      cuint component
   ) {
      
      Real* cp = cache.cells[fs_cache::calculateNbrID(1,1,1)]->parameters;
      
      switch (component) {
         case 0:
            //             cp[CellParams::EXGRADPE_000_100] = 0.0;
            //             cp[CellParams::EXGRADPE_010_110] = 0.0;
            //             cp[CellParams::EXGRADPE_001_101] = 0.0;
            //             cp[CellParams::EXGRADPE_011_111] = 0.0;
            cp[CellParams::EXGRADPE] = 0.0;
            break;
         case 1:
            //             cp[CellParams::EYGRADPE_000_010] = 0.0;
            //             cp[CellParams::EYGRADPE_100_110] = 0.0;
            //             cp[CellParams::EYGRADPE_001_011] = 0.0;
            //             cp[CellParams::EYGRADPE_101_111] = 0.0;
            cp[CellParams::EYGRADPE] = 0.0;
            break;
         case 2:
            //             cp[CellParams::EZGRADPE_000_001] = 0.0;
            //             cp[CellParams::EZGRADPE_100_101] = 0.0;
            //             cp[CellParams::EZGRADPE_010_011] = 0.0;
            //             cp[CellParams::EZGRADPE_110_111] = 0.0;
            cp[CellParams::EZGRADPE] = 0.0;
            break;
         default:
            cerr << __FILE__ << ":" << __LINE__ << ":" << " Invalid component" << endl;
      }
   }
   
   Real Outflow::fieldBoundaryCopyFromExistingFaceNbrMagneticField(
      const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const CellID& cellID,
      cuint& component
   ) {
      const CellID closestCell = getTheClosestNonsysboundaryCell(cellID);
      
      #ifdef DEBUG_OUTFLOW
      if (mpiGrid[closestCell]->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY) {
         stringstream ss;
         ss << "ERROR, outflow cell " << cellID << " uses value from sysboundary nbr " << closestCell;
         ss << " in " << __FILE__ << ":" << __LINE__ << endl;
         cerr << ss.str();
         exit(1);
      }
      
      if (closestCell == INVALID_CELLID) {
         cerr << cellID << " " << __FILE__ << ":" << __LINE__ << ": No closest cell found!" << endl;
         abort();
      }
      #endif
      
      return mpiGrid[closestCell]->parameters[CellParams::PERBX+component];
   }
   
   void Outflow::fieldSolverBoundaryCondDerivatives(
      dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const CellID& cellID,
      cuint& RKCase,
      cuint& component
   ) {
      this->setCellDerivativesToZero(mpiGrid, cellID, component);
   }
   
   void Outflow::fieldSolverBoundaryCondBVOLDerivatives(
      const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const CellID& cellID,
      cuint& component
   ) {
      this->setCellBVOLDerivativesToZero(mpiGrid, cellID, component);
   }
   
   /**
    * NOTE that this is called once for each particle species!
    * @param mpiGrid
    * @param cellID
    */
   void Outflow::vlasovBoundaryCondition(
      const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const CellID& cellID,
      const uint popID
   ) {
//      phiprof::start("vlasovBoundaryCondition (Outflow)");
      
      const OutflowSpeciesParameters& sP = this->speciesParams[popID];
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
         if(isThisCellOnAFace[i] && facesToProcess[i] && !sP.facesToSkipVlasov[i]) {
            switch(sP.faceVlasovScheme[i]) {
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
