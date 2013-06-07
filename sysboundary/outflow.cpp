/*
 * This file is part of Vlasiator.
 * 
 * Copyright 2010, 2011, 2012, 2013 Finnish Meteorological Institute
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

/*!\file outflow.cpp
 * \brief Implementation of the class SysBoundaryCondition::Outflow to handle cells classified as sysboundarytype::OUTFLOW.
 */

#include <cstdlib>
#include <iostream>

#include "outflow.h"
#include "../projects/projects_common.h"
#include "../fieldsolver.h"

using namespace std;

namespace SBC {
   Outflow::Outflow(): SysBoundaryCondition() { }
   Outflow::~Outflow() { }
   
   void Outflow::addParameters() {
      Readparameters::addComposing("outflow.face", "List of faces on which outflow boundary conditions are to be applied ([xyz][+-]).");
      Readparameters::add("outflow.precedence", "Precedence value of the outflow system boundary condition (integer), the higher the stronger.", 4);
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
   
   bool Outflow::assignSysBoundary(dccrg::Dccrg<SpatialCell>& mpiGrid) {
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
      const dccrg::Dccrg<SpatialCell>& mpiGrid,
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
   
   Real Outflow::fieldSolverBoundaryCondMagneticField(
      const dccrg::Dccrg<SpatialCell>& mpiGrid,
      const CellID& cellID,
      creal& dt,
      cuint& component
   ) {
      if(dt == 0.0) {
         return fieldBoundaryCopyFromExistingFaceNbrMagneticField(mpiGrid, cellID, component);
      } else { // Return PERB[XYZ]_DT2
         cint offset = CellParams::PERBX_DT2 - CellParams::PERBX;
         return fieldBoundaryCopyFromExistingFaceNbrMagneticField(mpiGrid, cellID, component+offset);
      }
   }
   
   void Outflow::fieldSolverBoundaryCondElectricField(
      dccrg::Dccrg<SpatialCell>& mpiGrid,
      const CellID& cellID,
      cuint RKCase,
      cuint component
   ) {
      switch(component) {
         case 0:
            calculateEdgeElectricFieldX(mpiGrid, cellID, RKCase);
            break;
         case 1:
            calculateEdgeElectricFieldY(mpiGrid, cellID, RKCase);
            break;
         case 2:
            calculateEdgeElectricFieldZ(mpiGrid, cellID, RKCase);
            break;
         default:
            cerr << "ERROR: Reached end of switch in Outflow::fieldSolverBoundaryCondElectricField." << endl;
            abort();
      }
   }
   
   Real Outflow::fieldBoundaryCopyFromExistingFaceNbrMagneticField(
      const dccrg::Dccrg<SpatialCell>& mpiGrid,
      const CellID& cellID,
      cuint& component
   ) {
      const CellID closestCell = getClosestNonsysboundaryCell(mpiGrid, cellID);
      
      if(closestCell == INVALID_CELLID) {
         cerr << cellID << " " << __FILE__ << ":" << __LINE__ << ": No closest cell found!" << endl;
         abort();
      }
      
      return mpiGrid[closestCell]->parameters[CellParams::PERBX+component];
      
/*
 // get_neighbors                                                *
 // if all neighbors are of the same size then they are in z order, e.g. with a neighborhood size of 2 the first neighbor is at offset (-2, -2, -2) from the given cell, the second one is at (-1, -2, -2), etc, in size units of the given cell.
 
 const vector<CellID>* neighbors = mpiGrid.get_neighbors(cellID);

      // CORRECT
      // -1  0  0 62
      //  0 -1  0 58
      //  0  0 -1 38
      //  1  0  0 63
      //  0  1  0 67
      //  0  0  1 87
      // WARNING minus one, C++ array!!
      
      // If -x/-y/-z neighbours are missing, copy from +x/+y/+z
      if(((*neighbors)[61] == dccrg::error_cell) &&
         ((*neighbors)[57] == dccrg::error_cell) &&
         ((*neighbors)[37] == dccrg::error_cell)) {
         return mpiGrid[getNeighbour(mpiGrid,cellID,+1,+1,+1)]->parameters[CellParams::PERBX+component];
      }
      // If -x/-y/+z neighbours are missing, copy from +x/+y/-z
      if(((*neighbors)[61] == dccrg::error_cell) &&
         ((*neighbors)[57] == dccrg::error_cell) &&
         ((*neighbors)[86] == dccrg::error_cell)) {
         return mpiGrid[getNeighbour(mpiGrid,cellID,+1,+1,-1)]->parameters[CellParams::PERBX+component];
      }
      // If -x/+y/-z neighbours are missing, copy from +x/-y/+z
      if(((*neighbors)[61] == dccrg::error_cell) &&
         ((*neighbors)[66] == dccrg::error_cell) &&
         ((*neighbors)[37] == dccrg::error_cell)) {
         return mpiGrid[getNeighbour(mpiGrid,cellID,+1,-1,+1)]->parameters[CellParams::PERBX+component];
      }
      // If +x/-y/-z neighbours are missing, copy from -x/+y/+z
      if(((*neighbors)[62] == dccrg::error_cell) &&
         ((*neighbors)[57] == dccrg::error_cell) &&
         ((*neighbors)[37] == dccrg::error_cell)) {
         return mpiGrid[getNeighbour(mpiGrid,cellID,-1,+1,+1)]->parameters[CellParams::PERBX+component];
      }
      // If -x/+y/+z neighbours are missing, copy from +x/-y/-z
      if(((*neighbors)[61] == dccrg::error_cell) &&
         ((*neighbors)[66] == dccrg::error_cell) &&
         ((*neighbors)[86] == dccrg::error_cell)) {
         return mpiGrid[getNeighbour(mpiGrid,cellID,+1,-1,-1)]->parameters[CellParams::PERBX+component];
      }
      // If +x/-y/+z neighbours are missing, copy from -x/+y/-z
      if(((*neighbors)[62] == dccrg::error_cell) &&
         ((*neighbors)[57] == dccrg::error_cell) &&
         ((*neighbors)[86] == dccrg::error_cell)) {
         return mpiGrid[getNeighbour(mpiGrid,cellID,-1,+1,-1)]->parameters[CellParams::PERBX+component];
      }
      // If +x/+y/-z neighbours are missing, copy from -x/-y/+z
      if(((*neighbors)[62] == dccrg::error_cell) &&
         ((*neighbors)[66] == dccrg::error_cell) &&
         ((*neighbors)[37] == dccrg::error_cell)) {
         return mpiGrid[getNeighbour(mpiGrid,cellID,-1,-1,+1)]->parameters[CellParams::PERBX+component];
      }
      // If +x/+y/+z neighbours are missing, copy from -x/-y/-z
      if(((*neighbors)[62] == dccrg::error_cell) &&
         ((*neighbors)[66] == dccrg::error_cell) &&
         ((*neighbors)[86] == dccrg::error_cell)) {
         return mpiGrid[getNeighbour(mpiGrid,cellID,-1,-1,-1)]->parameters[CellParams::PERBX+component];
      }
      
      // If -z/-y neighbours are missing, copy from +z/+y
      if(((*neighbors)[37] == dccrg::error_cell) &&
         ((*neighbors)[57] == dccrg::error_cell)) {
         return mpiGrid[getNeighbour(mpiGrid,cellID,0,+1,+1)]->parameters[CellParams::PERBX+component];
      }
      // If -z/-x neighbours are missing, copy from +z/+x
      if(((*neighbors)[37] == dccrg::error_cell) &&
         ((*neighbors)[61] == dccrg::error_cell)) {
         return mpiGrid[getNeighbour(mpiGrid,cellID,+1,0,+1)]->parameters[CellParams::PERBX+component];
      }
      // If -z/+y neighbours are missing, copy from +z/-y
      if(((*neighbors)[37] == dccrg::error_cell) &&
         ((*neighbors)[66] == dccrg::error_cell)) {
         return mpiGrid[getNeighbour(mpiGrid,cellID,0,-1,+1)]->parameters[CellParams::PERBX+component];
      }
      // If -z/+x neighbours are missing, copy from +z/-x
      if(((*neighbors)[37] == dccrg::error_cell) &&
         ((*neighbors)[62] == dccrg::error_cell)) {
         return mpiGrid[getNeighbour(mpiGrid,cellID,-1,0,+1)]->parameters[CellParams::PERBX+component];
      }
      // If +z/-y neighbours are missing, copy from -z/+y
      if(((*neighbors)[86] == dccrg::error_cell) &&
         ((*neighbors)[57] == dccrg::error_cell)) {
         return mpiGrid[getNeighbour(mpiGrid,cellID,0,+1,-1)]->parameters[CellParams::PERBX+component];
      }
      // If +z/+x neighbours are missing, copy from -z/-x
      if(((*neighbors)[62] == dccrg::error_cell) &&
         ((*neighbors)[86] == dccrg::error_cell)) {
         return mpiGrid[getNeighbour(mpiGrid,cellID,-1,0,-1)]->parameters[CellParams::PERBX+component];
      }
      // If +z/+y neighbours are missing, copy from -z/-y
      if(((*neighbors)[86] == dccrg::error_cell) &&
         ((*neighbors)[66] == dccrg::error_cell)) {
         return mpiGrid[getNeighbour(mpiGrid,cellID,0,-1,-1)]->parameters[CellParams::PERBX+component];
      }
      // If +z/-x neighbours are missing, copy from -z/+x
      if(((*neighbors)[61] == dccrg::error_cell) &&
         ((*neighbors)[86] == dccrg::error_cell)) {
         return mpiGrid[getNeighbour(mpiGrid,cellID,+1,0,-1)]->parameters[CellParams::PERBX+component];
      }
      // If -x/-y neighbours are missing, copy from +x/+y
      if(((*neighbors)[61] == dccrg::error_cell) &&
         ((*neighbors)[57] == dccrg::error_cell)) {
         return mpiGrid[getNeighbour(mpiGrid,cellID,+1,+1,0)]->parameters[CellParams::PERBX+component];
      }
      // If -x/+y neighbours are missing, copy from +x/-y
      if(((*neighbors)[61] == dccrg::error_cell) &&
         ((*neighbors)[66] == dccrg::error_cell)) {
         return mpiGrid[getNeighbour(mpiGrid,cellID,+1,-1,0)]->parameters[CellParams::PERBX+component];
      }
      // If +x/-y neighbours are missing, copy from -x/+y
      if(((*neighbors)[62] == dccrg::error_cell) &&
         ((*neighbors)[57] == dccrg::error_cell)) {
         return mpiGrid[getNeighbour(mpiGrid,cellID,-1,+1,0)]->parameters[CellParams::PERBX+component];
      }
      // If +x/+y neighbours are missing, copy from -x/-y
      if(((*neighbors)[62] == dccrg::error_cell) &&
         ((*neighbors)[66] == dccrg::error_cell)) {
         return mpiGrid[getNeighbour(mpiGrid,cellID,-1,-1,0)]->parameters[CellParams::PERBX+component];
      }
      // If -z neighbour is missing, copy from +z
      if((*neighbors)[37] == dccrg::error_cell) {
         return mpiGrid[getNeighbour(mpiGrid,cellID,0,0,+1)]->parameters[CellParams::PERBX+component];
      }
      // If -y neighbour is missing, copy from +y
      if((*neighbors)[57] == dccrg::error_cell) {
         return mpiGrid[getNeighbour(mpiGrid,cellID,0,+1,0)]->parameters[CellParams::PERBX+component];
      }
      // If -x neighbour is missing, copy from +x
      if((*neighbors)[61] == dccrg::error_cell) {
         return mpiGrid[getNeighbour(mpiGrid,cellID,+1,0,0)]->parameters[CellParams::PERBX+component];
      }
      // If +z neighbour is missing, copy from -z
      if((*neighbors)[86] == dccrg::error_cell) {
         return mpiGrid[getNeighbour(mpiGrid,cellID,0,0,-1)]->parameters[CellParams::PERBX+component];
      }
      // If +y neighbour is missing, copy from -y
      if((*neighbors)[66] == dccrg::error_cell) {
         return mpiGrid[getNeighbour(mpiGrid,cellID,0,-1,0)]->parameters[CellParams::PERBX+component];
      }
      // If +x neighbour is missing, copy from -x
      if((*neighbors)[62] == dccrg::error_cell) {
         return mpiGrid[getNeighbour(mpiGrid,cellID,-1,0,0)]->parameters[CellParams::PERBX+component];
      }
      cerr << cellID << " " << __FILE__  << ":" << __LINE__
           << ": Got to the end of fieldBoundaryCopyFromExistingFaceNbrMagneticField(), this should not happen!"
           << endl;
      return 0;*/
   }
   
   void Outflow::fieldSolverBoundaryCondDerivatives(
      dccrg::Dccrg<SpatialCell>& mpiGrid,
      const CellID& cellID,
      cuint& RKCase,
      cuint& component
   ) {
      this->setCellDerivativesToZero(mpiGrid, cellID, component);
   }
   
   void Outflow::fieldSolverBoundaryCondBVOLDerivatives(
      const dccrg::Dccrg<SpatialCell>& mpiGrid,
      const CellID& cellID,
      cuint& component
   ) {
      this->setCellBVOLDerivativesToZero(mpiGrid, cellID, component);
   }
   
   void Outflow::vlasovBoundaryCondition(
      const dccrg::Dccrg<SpatialCell>& mpiGrid,
      const CellID& cellID
   ) {
//      phiprof::start("vlasovBoundaryCondition (Outflow)");
      vlasovBoundaryCopyFromExistingFaceNbr(mpiGrid, cellID);
//      phiprof::stop("vlasovBoundaryCondition (Outflow)");
   }
   
   void Outflow::vlasovBoundaryCopyFromExistingFaceNbr(
      const dccrg::Dccrg<SpatialCell>& mpiGrid,
      const CellID& cellID
   ) {
      const CellID closestCell = getClosestNonsysboundaryCell(mpiGrid, cellID);
      
      if(closestCell == INVALID_CELLID) {
         cerr << __FILE__ << ":" << __LINE__ << ": No closest cell found!" << endl;
         abort();
      }
      //Do not allow block adjustment, the block structure when calling vlasovBoundaryCondition should be static
      copyCellData(mpiGrid[closestCell], mpiGrid[cellID],false);
   }
   
   void Outflow::getFaces(bool* faces) {
      for(uint i=0; i<6; i++) faces[i] = facesToProcess[i];
   }
   
   std::string Outflow::getName() const {return "Outflow";}
   uint Outflow::getIndex() const {return sysboundarytype::OUTFLOW;}
}
