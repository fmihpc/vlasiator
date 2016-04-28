/*
 * This file is part of Vlasiator.
 * 
 * Copyright 2015 Finnish Meteorological Institute
 * 
 */

/*!\file antisymmetric.cpp
 * \brief Implementation of the class SysBoundaryCondition::Antisymmetric to 
 * handle cells classified as sysboundarytype::ANTISYMMETRIC.
 */

#include <cstdlib>
#include <iostream>

#include "antisymmetric.h"
#include "../projects/projects_common.h"
#include "../fieldsolver/fs_common.h"

#ifndef NDEBUG
   #define DEBUG_ANTISYMMETRIC
#endif
#ifdef DEBUG_SYSBOUNDARY
   #define DEBUG_ANTISYMMETRIC
#endif

using namespace std;

namespace SBC {
   Antisymmetric::Antisymmetric(): SysBoundaryCondition() { }
   Antisymmetric::~Antisymmetric() { }
   
   void Antisymmetric::addParameters() {
      Readparameters::addComposing("antisymmetric.face", "List of faces on which outflow boundary conditions are to be applied ([xyz][+-]).");
      Readparameters::add("antisymmetric.precedence", "Precedence value of the outflow system boundary condition (integer), the higher the stronger.", 4);
   }
   
   void Antisymmetric::getParameters() {
      int myRank;
      MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
      if(!Readparameters::get("antisymmetric.face", faceList)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!Readparameters::get("antisymmetric.precedence", precedence)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
   }
   
   bool Antisymmetric::initSysBoundary(
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
   
   bool Antisymmetric::assignSysBoundary(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid) {
      const vector<CellID>& cells = getLocalCells();
      for (size_t c=0; c<cells.size(); ++c) {
         if (mpiGrid[cells[c]]->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) continue;
         creal* const cellParams = mpiGrid[cells[c]]->parameters;
         creal dx = cellParams[CellParams::DX];
         creal dy = cellParams[CellParams::DY];
         creal dz = cellParams[CellParams::DZ];
         creal x = cellParams[CellParams::XCRD] + 0.5*dx;
         creal y = cellParams[CellParams::YCRD] + 0.5*dy;
         creal z = cellParams[CellParams::ZCRD] + 0.5*dz;

         bool isThisCellOnAFace[6];
         determineFace(&isThisCellOnAFace[0], x, y, z, dx, dy, dz);

         // Comparison of the array defining which faces to use and the 
         // array telling on which faces this cell is
         bool doAssign = false;
         for (uint j=0; j<6; j++) doAssign = doAssign || (facesToProcess[j] && isThisCellOnAFace[j]);
         if (doAssign) {
            uint flag = getIndex();
            //if (x <  Parameters::xmin + 2*Parameters::dx_ini) flag = sysboundarytype::DO_NOT_COMPUTE;
            //if (x >= Parameters::xmax - 2*Parameters::dx_ini) flag = sysboundarytype::DO_NOT_COMPUTE;
            if (y < Parameters::ymin+Parameters::dy_ini)      flag = sysboundarytype::DO_NOT_COMPUTE;

            mpiGrid[cells[c]]->sysBoundaryFlag = flag;
         }
      }
      return true;
   }
   
   bool Antisymmetric::applyInitialState(
      const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      Project &project
   ) {
      vector<uint64_t> cells = mpiGrid.get_cells();
      #pragma omp parallel for
      for (uint i=0; i<cells.size(); ++i) {
         SpatialCell* cell = mpiGrid[cells[i]];
         if(cell->sysBoundaryFlag != this->getIndex()) continue;
         
         // Defined in project.cpp, used here as the outflow cell has the same state as the initial state of non-system boundary cells.
         project.setCell(cell);
         // WARNING Time-independence assumed here.
         cell->parameters[CellParams::RHO_DT2] = cell->parameters[CellParams::RHO];
         cell->parameters[CellParams::RHOVX_DT2] = cell->parameters[CellParams::RHOVX];
         cell->parameters[CellParams::RHOVY_DT2] = cell->parameters[CellParams::RHOVY];
         cell->parameters[CellParams::RHOVZ_DT2] = cell->parameters[CellParams::RHOVZ];
      }
      
      return true;
   }
   
   Real Antisymmetric::fieldSolverBoundaryCondMagneticField(
                                                      const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                                      const CellID& cellID,
                                                      creal& dt,
                                                      cuint& component
   ) {
      if (dt == 0.0) {
         return fieldBoundaryCopyFromExistingFaceNbrMagneticField(mpiGrid, cellID, component);
      } else { // Return PERB[XYZ]_DT2
         cint offset = CellParams::PERBX_DT2 - CellParams::PERBX;
         return fieldBoundaryCopyFromExistingFaceNbrMagneticField(mpiGrid, cellID, component+offset);
      }
   }

   void Antisymmetric::fieldSolverBoundaryCondElectricField(
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
   
   void Antisymmetric::fieldSolverBoundaryCondHallElectricField(
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
   
   Real Antisymmetric::fieldBoundaryCopyFromExistingFaceNbrMagneticField(
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
   
   void Antisymmetric::fieldSolverBoundaryCondDerivatives(
      dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const CellID& cellID,
      cuint& RKCase,
      cuint& component
   ) {
      this->setCellDerivativesToZero(mpiGrid, cellID, component);
   }
   
   void Antisymmetric::fieldSolverBoundaryCondBVOLDerivatives(
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
   void Antisymmetric::vlasovBoundaryCondition(
      const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const CellID& cellID,
      const int& popID
   ) {
      //cerr << "AS vlasovBoundaryCondition cell " << cellID << " called " << endl;
      
//      phiprof::start("vlasovBoundaryCondition (Antisymmetric)");
//      vlasovBoundaryCopyFromTheClosestNbr(mpiGrid,cellID,popID);
//      phiprof::stop("vlasovBoundaryCondition (Antisymmetric)");

      dccrg::Types<3>::indices_t indices = mpiGrid.mapping.get_indices(cellID);      
      if (indices[1] == 1 || indices[1] == Parameters::ycells_ini-2) {
         if (indices[1] == 1) ++indices[1];
         else                 --indices[1];
      } else {
         cerr << "ERROR in antisymmetric boundaries in " << __FILE__ << ":" << __LINE__ << endl;
         return;
      }
      
      // Clear data in this cell
      vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& vmesh    = mpiGrid[cellID]->get_velocity_mesh(popID);
      vmesh::VelocityBlockContainer<vmesh::LocalID>& blockContainer = mpiGrid[cellID]->get_velocity_blocks(popID);
      vmesh.clear();
      blockContainer.clear();
      
      // Get data from neighbor cell
      const CellID nbrID = mpiGrid.mapping.get_cell_from_indices(indices,0);
      const vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& nbrVmesh    = mpiGrid[nbrID]->get_velocity_mesh(popID);
      const vmesh::VelocityBlockContainer<vmesh::LocalID>& nbrBlockContainer = mpiGrid[nbrID]->get_velocity_blocks(popID);

      //cerr << "BOUNDARY CELL " << cellID << " HAS " << vmesh.size() << " BLOCKS BEFORE COPY NBR " << nbrID << " HAS ";
      //cerr << nbrVmesh.size() << endl;
      
      // vx,vz components are untouched, vy is both mirrored (vy -> -vy)
      // and copied (vy -> vy).
      for (vmesh::LocalID srcLID=0; srcLID<nbrVmesh.size(); ++srcLID) {
         Real V[3];
         const vmesh::GlobalID srcGID = nbrVmesh.getGlobalID(srcLID);
         nbrVmesh.getBlockCoordinates(srcGID,V);
         
         Real dV[3];
         nbrVmesh.getBlockSize(srcGID,dV);
         for (int t=0; t<3; ++t) dV[t] /= WID;

         // Pointer to source data
         const Realf* srcData = nbrBlockContainer.getData(srcLID);

         Real V_trgt[3];
         for (unsigned int k=0; k<WID; ++k) for (unsigned int j=0; j<WID; ++j) for (unsigned int i=0; i<WID; ++i) {
//            if (V[1] + (j+0.5)*dV[1] > 0) continue;

            for (int dir=-1; dir<0; dir+=2) {
//            for (int dir=-1; dir<2; dir+=2) {
               V_trgt[0] =  V[0] + (i+0.5)*dV[0];
               V_trgt[1] = (V[1] + (j+0.5)*dV[1]) * dir;
               V_trgt[2] =  V[2] + (k+0.5)*dV[2];
            
               const vmesh::GlobalID trgtGID = vmesh.getGlobalID(0,V_trgt);
            
               Real V_trgt_block[3];
               Real dV_trgt_block[3];
               vmesh.getBlockCoordinates(trgtGID,V_trgt_block);
               vmesh.getBlockSize(trgtGID,dV_trgt_block);
               for (int t=0; t<3; ++t) dV_trgt_block[t] /= WID;

               // Make sure target block exists
               if (vmesh.getLocalID(trgtGID) == vmesh::INVALID_LOCALID) {
                  mpiGrid[cellID]->add_velocity_block(trgtGID,popID);
               }
            
               // Get target block local ID
               const vmesh::LocalID trgtLID = vmesh.getLocalID(trgtGID);
               if (trgtLID == vmesh::INVALID_LOCALID) {
                  cerr << "ERROR, got invalid local ID in antisymmetric" << endl;
                  continue;
               }

               // Calculate target cell indices
               const int ii = static_cast<int>((V_trgt[0] - V_trgt_block[0]) / dV_trgt_block[0]);
               const int jj = static_cast<int>((V_trgt[1] - V_trgt_block[1]) / dV_trgt_block[1]);
               const int kk = static_cast<int>((V_trgt[2] - V_trgt_block[2]) / dV_trgt_block[2]);
            
               /*bool ok = true;
               if (ii < 0 || ii >= WID) ok = false;
               if (jj < 0 || jj >= WID) ok = false;
               if (kk < 0 || kk >= WID) ok = false;
               if (ok == false) {
                  cerr << "ERROR " << ii << ' ' << jj << ' ' << kk << endl;
                  exit(1);               
               }*/

               /*uint8_t ref=0;
               vmesh::LocalID i_srcBlock[3];
               nbrVmesh.getIndices(srcGID,ref,i_srcBlock[0],i_srcBlock[1],i_srcBlock[2]);
            
               vmesh::LocalID i_trgtBlock[3];
               vmesh.getIndices(trgtGID,ref,i_trgtBlock[0],i_trgtBlock[1],i_trgtBlock[2]);
            
               cerr << "\t src indices: " << i_srcBlock[0] << ' ' << i_srcBlock[1] << ' ' << i_srcBlock[2] << " (";
               cerr << i << ' ' << j << ' ' << k << ") trgt ";
               cerr << i_trgtBlock[0] << ' ' << i_trgtBlock[1] << ' ' << i_trgtBlock[2] << " (";
               cerr << ii << ' ' << jj << ' ' << kk << ")" << endl;
               */
               Realf* data = blockContainer.getData(trgtLID);
               data[cellIndex(ii,jj,kk)] += srcData[cellIndex(i,j,k)];
            }
         } // for-loops over phase-space cells in source block
      } // for-loop over velocity blocks in neighbor cell
      
      //cerr << "BOUNDARY CELL HAS " << vmesh.size() << " velocity blocks" << endl;
   }

   void Antisymmetric::getFaces(bool* faces) {
      for(uint i=0; i<6; i++) faces[i] = facesToProcess[i];
   }
   
   std::string Antisymmetric::getName() const {return "Antisymmetric";}

   uint Antisymmetric::getIndex() const {return sysboundarytype::ANTISYMMETRIC;}
}
