/*
 * This file is part of Vlasiator.
 * 
 * Copyright 2010, 2011, 2012, 2013 Finnish Meteorological Institute
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
      Readparameters::addComposing("outflow.face", "List of faces on which outflow boundary conditions are to be applied ([xyz][+-]).");
      Readparameters::add("outflow.precedence", "Precedence value of the outflow system boundary condition (integer), the higher the stronger.", 4);
      Readparameters::add("outflow.reapplyUponRestart", "If 0 (default), keep going with the state existing in the restart file. If 1, calls again applyInitialState. Can be used to change boundary condition behaviour during a run.", 0);
   }
   
   void Outflow::getParameters() {
      int myRank;
      MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
      if(!Readparameters::get("outflow.face", faceList)) {
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
   }
   
   bool Outflow::initSysBoundary(
      creal& t,
      Project &project
   ) {
      /* The array of bool describes which of the x+, x-, y+, y-, z+, z- faces are to have outflow system boundary conditions.
       * A true indicates the corresponding face will have outflow.
       * The 6 elements correspond to x+, x-, y+, y-, z+, z- respectively.
       */
      for(uint i=0; i<6; i++) facesToProcess[i] = false;
      
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
         && component == 0                                  // we do the component normal to this face
         && !(isThisCellOnAFace[2] || isThisCellOnAFace[3] || isThisCellOnAFace[4] || isThisCellOnAFace[5]) // we are not in a corner
      ) {
         const vector<uint16_t> cellVector = {{localID}};
         propagateMagneticField(cellCache, cellVector, dt, RKCase, true, false, false);
         fieldValue = cellCache[localID].cells[fs_cache::calculateNbrID(1  ,1  ,1  )]->parameters[CellParams::PERBX + component + offset];
      } else if (sysBoundaryLayer == 1
         && isThisCellOnAFace[2]                            // we are on the face
         && this->facesToProcess[2]                         // we are supposed to do this face
         && component == 1                                  // we do the component normal to this face
         && !(isThisCellOnAFace[0] || isThisCellOnAFace[1] || isThisCellOnAFace[4] || isThisCellOnAFace[5]) // we are not in a corner
      ) {
         const vector<uint16_t> cellVector = {{localID}};
         propagateMagneticField(cellCache, cellVector, dt, RKCase, false, true, false);
         fieldValue = cellCache[localID].cells[fs_cache::calculateNbrID(1  ,1  ,1  )]->parameters[CellParams::PERBX + component + offset];
      } else if (sysBoundaryLayer == 1
         && isThisCellOnAFace[4]                            // we are on the face
         && this->facesToProcess[4]                         // we are supposed to do this face
         && component == 2                                  // we do the component normal to this face
         && !(isThisCellOnAFace[0] || isThisCellOnAFace[1] || isThisCellOnAFace[2] || isThisCellOnAFace[3]) // we are not in a corner
      ) {
         const vector<uint16_t> cellVector = {{localID}};
         propagateMagneticField(cellCache, cellVector, dt, RKCase, false, false, true);
         fieldValue = cellCache[localID].cells[fs_cache::calculateNbrID(1  ,1  ,1  )]->parameters[CellParams::PERBX + component + offset];
      } else {
         fieldValue = fieldBoundaryCopyFromExistingFaceNbrMagneticField(mpiGrid, cellID, component + offset);
      }
      
      return fieldValue;
   }

   void Outflow::fieldSolverBoundaryCondElectricField(
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
   
   void Outflow::fieldSolverBoundaryCondHallElectricField(
      fs_cache::CellCache& cache,
      cuint RKCase,
      cuint component
   ) {

      Real* cp = cache.cells[fs_cache::calculateNbrID(1,1,1)]->parameters;
      
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
      const int& popID
   ) {
//      phiprof::start("vlasovBoundaryCondition (Outflow)");
      vlasovBoundaryCopyFromTheClosestNbr(mpiGrid,cellID,popID);
//      phiprof::stop("vlasovBoundaryCondition (Outflow)");
   }
   
   void Outflow::getFaces(bool* faces) {
      for(uint i=0; i<6; i++) faces[i] = facesToProcess[i];
   }
   
   std::string Outflow::getName() const {return "Outflow";}
   uint Outflow::getIndex() const {return sysboundarytype::OUTFLOW;}
}
