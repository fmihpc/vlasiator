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

/*!\file sysboundarycondition.cpp
 * \brief Implementation of the base class SysBoundaryCondition to handle system boundary cells.
 * 
 * \sa donotcompute.cpp ionosphere.cpp outflow.cpp setbyuser.cpp setmaxwellian.cpp
 * 
 */

#include <cstdlib>
#include <iostream>

#include "../parameters.h"
#include "../vlasovmover.h"
#include "sysboundarycondition.h"
#include "../projects/projects_common.h"

using namespace std;

namespace SBC {
   // ************************************************************
   // ***** DEFINITIONS FOR BOUNDARYCONDITION BASE CLASS *****
   // ************************************************************
   
   /*!\brief Function used to determine on which face(s) if any the cell at given coordinates is.
    * 
    * This function is used by some of the classes inheriting from this base class.
    * 
    * Depth is hard-coded to be 2 as other parts of the code (field solver especially) rely on that.
    * 
    * \param isThisCellOnAFace Pointer to an array of 6 bool returning of each face whether the cell is on that face. Order: 0 x+; 1 x-; 2 y+; 3 y-; 4 z+; 5 z-
    * \param x Cell x coordinate
    * \param y Cell y coordinate
    * \param z Cell z coordinate
    * \param dx Cell dx size
    * \param dy Cell dy size
    * \param dz Cell dz size
    * \param excludeSlicesAndPeriodicDimensions If true, do not consider a cell to be part of the face if that face has a depth of 1 only (single-cell thick slices/columns) or if that direciton is periodic..
    */
   void SysBoundaryCondition::determineFace(
      bool* isThisCellOnAFace,
      creal x, creal y, creal z,
      creal dx, creal dy, creal dz,
      const bool excludeSlicesAndPeriodicDimensions //=false (default)
   ) {
      for(uint i=0; i<6; i++) {
         isThisCellOnAFace[i] = false;
      }
      if(x > Parameters::xmax - 2.0*dx) {
         isThisCellOnAFace[0] = true;
      }
      if(x < Parameters::xmin + 2.0*dx) {
         isThisCellOnAFace[1] = true;
      }
      if(y > Parameters::ymax - 2.0*dy) {
         isThisCellOnAFace[2] = true;
      }
      if(y < Parameters::ymin + 2.0*dy) {
         isThisCellOnAFace[3] = true;
      }
      if(z > Parameters::zmax - 2.0*dz) {
         isThisCellOnAFace[4] = true;
      }
      if(z < Parameters::zmin + 2.0*dz) {
         isThisCellOnAFace[5] = true;
      }
      if(excludeSlicesAndPeriodicDimensions == true) {
         if(Parameters::xcells_ini == 1 || this->isPeriodic[0]) {
            isThisCellOnAFace[0] = false;
            isThisCellOnAFace[1] = false;
         }
         if(Parameters::ycells_ini == 1 || this->isPeriodic[1]) {
            isThisCellOnAFace[2] = false;
            isThisCellOnAFace[3] = false;
         }
         if(Parameters::zcells_ini == 1 || this->isPeriodic[2]) {
            isThisCellOnAFace[4] = false;
            isThisCellOnAFace[5] = false;
         }
      }
   }
   
   /*! SysBoundaryCondition base class constructor. The constructor is empty.*/
   SysBoundaryCondition::SysBoundaryCondition() { }
   
   /*! SysBoundaryCondition base class virtual destructor. The destructor is empty.*/
   SysBoundaryCondition::~SysBoundaryCondition() { }
   
   /*! SysBoundaryCondition base class instance of the addParameters function. Should not be used, each derived class should have its own.*/
   void SysBoundaryCondition::addParameters() {
      cerr << "ERROR: SysBoundaryCondition::addParameters called instead of derived class function!" << endl;
   }
   
   /*! Function used to set the system boundary condition cell's derivatives to 0.
    * \param mpiGrid Grid
    * \param cellID The cell's ID.
    * \param component 0: x-derivatives, 1: y-derivatives, 2: z-derivatives, 3: xy-derivatives, 4: xz-derivatives, 5: yz-derivatives.
    */
   void SysBoundaryCondition::setCellDerivativesToZero(
      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, 2> & dPerBGrid,
      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 2> & dMomentsGrid,
      cint i,
      cint j,
      cint k,
      cuint& component
   ) {
      std::array<Real, fsgrids::dperb::N_DPERB> * dPerBGrid0 = dPerBGrid.get(i,j,k);
      std::array<Real, fsgrids::dmoments::N_DMOMENTS> * dMomentsGrid0 = dMomentsGrid.get(i,j,k);
      switch(component) {
         case 0: // x, xx
            dMomentsGrid0->at(fsgrids::dmoments::drhomdx) = 0.0;
            dMomentsGrid0->at(fsgrids::dmoments::drhoqdx) = 0.0;
            dMomentsGrid0->at(fsgrids::dmoments::dp11dx) = 0.0;
            dMomentsGrid0->at(fsgrids::dmoments::dp22dx) = 0.0;
            dMomentsGrid0->at(fsgrids::dmoments::dp33dx) = 0.0;
            dPerBGrid0->at(fsgrids::dperb::dPERBydx)  = 0.0;
            dPerBGrid0->at(fsgrids::dperb::dPERBzdx)  = 0.0;
            dMomentsGrid0->at(fsgrids::dmoments::dVxdx)  = 0.0;
            dMomentsGrid0->at(fsgrids::dmoments::dVydx)  = 0.0;
            dMomentsGrid0->at(fsgrids::dmoments::dVzdx)  = 0.0;
            dPerBGrid0->at(fsgrids::dperb::dPERBydxx) = 0.0;
            dPerBGrid0->at(fsgrids::dperb::dPERBzdxx) = 0.0;
            break;
         case 1: // y, yy
            dMomentsGrid0->at(fsgrids::dmoments::drhomdy) = 0.0;
            dMomentsGrid0->at(fsgrids::dmoments::drhoqdy) = 0.0;
            dMomentsGrid0->at(fsgrids::dmoments::dp11dy) = 0.0;
            dMomentsGrid0->at(fsgrids::dmoments::dp22dy) = 0.0;
            dMomentsGrid0->at(fsgrids::dmoments::dp33dy) = 0.0;
            dPerBGrid0->at(fsgrids::dperb::dPERBxdy)  = 0.0;
            dPerBGrid0->at(fsgrids::dperb::dPERBzdy)  = 0.0;
            dMomentsGrid0->at(fsgrids::dmoments::dVxdy)  = 0.0;
            dMomentsGrid0->at(fsgrids::dmoments::dVydy)  = 0.0;
            dMomentsGrid0->at(fsgrids::dmoments::dVzdy)  = 0.0;
            dPerBGrid0->at(fsgrids::dperb::dPERBxdyy) = 0.0;
            dPerBGrid0->at(fsgrids::dperb::dPERBzdyy) = 0.0;
            break;
         case 2: // z, zz
            dMomentsGrid0->at(fsgrids::dmoments::drhomdz) = 0.0;
            dMomentsGrid0->at(fsgrids::dmoments::drhoqdz) = 0.0;
            dMomentsGrid0->at(fsgrids::dmoments::dp11dz) = 0.0;
            dMomentsGrid0->at(fsgrids::dmoments::dp22dz) = 0.0;
            dMomentsGrid0->at(fsgrids::dmoments::dp33dz) = 0.0;
            dPerBGrid0->at(fsgrids::dperb::dPERBxdz)  = 0.0;
            dPerBGrid0->at(fsgrids::dperb::dPERBydz)  = 0.0;
            dMomentsGrid0->at(fsgrids::dmoments::dVxdz)  = 0.0;
            dMomentsGrid0->at(fsgrids::dmoments::dVydz)  = 0.0;
            dMomentsGrid0->at(fsgrids::dmoments::dVzdz)  = 0.0;
            dPerBGrid0->at(fsgrids::dperb::dPERBxdzz) = 0.0;
            dPerBGrid0->at(fsgrids::dperb::dPERBydzz) = 0.0;
            break;
         case 3: // xy
            dPerBGrid0->at(fsgrids::dperb::dPERBzdxy) = 0.0;
            break;
         case 4: // xz
            dPerBGrid0->at(fsgrids::dperb::dPERBydxz) = 0.0;
            break;
         case 5: // yz
            dPerBGrid0->at(fsgrids::dperb::dPERBxdyz) = 0.0;
            break;
         default:
            cerr << __FILE__ << ":" << __LINE__ << ":" << " Invalid component" << endl;
      }
   }
   
   /*! Function used to set the system boundary condition cell's BVOL derivatives to 0.
    * \param mpiGrid Grid
    * \param cellID The cell's ID.
    * \param component 0: x-derivatives, 1: y-derivatives, 2: z-derivatives.
    */
   void SysBoundaryCondition::setCellBVOLDerivativesToZero(
      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, 2> & volGrid,
      cint i,
      cint j,
      cint k,
      cuint& component
   ) {
      std::array<Real, fsgrids::volfields::N_VOL> * volGrid0 = volGrid.get(i,j,k);
      switch(component) {
         case 0:
            volGrid0->at(fsgrids::volfields::dPERBYVOLdx) = 0.0;
            volGrid0->at(fsgrids::volfields::dPERBZVOLdx) = 0.0;
            break;
         case 1:
            volGrid0->at(fsgrids::volfields::dPERBXVOLdy) = 0.0;
            volGrid0->at(fsgrids::volfields::dPERBZVOLdy) = 0.0;
            break;
         case 2:
            volGrid0->at(fsgrids::volfields::dPERBXVOLdz) = 0.0;
            volGrid0->at(fsgrids::volfields::dPERBYVOLdz) = 0.0;
            break;
         default:
            cerr << __FILE__ << ":" << __LINE__ << ":" << " Invalid component" << endl;
      }
   }
   
   /*! Function used to copy the distribution and moments from (one of) the closest sysboundarytype::NOT_SYSBOUNDARY cell.
    * \param mpiGrid Grid
    * \param cellID The cell's ID.
    * \param copyMomentsOnly If true, do not touch velocity space.
    */
   void SysBoundaryCondition::vlasovBoundaryCopyFromTheClosestNbr(
         const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
         const CellID& cellID,
         const bool& copyMomentsOnly,
         const uint popID,
         const bool calculate_V_moments
   ) {
      const CellID closestCell = getTheClosestNonsysboundaryCell(cellID);
      
      if(closestCell == INVALID_CELLID) {
         cerr << __FILE__ << ":" << __LINE__ << ": No closest cell found!" << endl;
         abort();
      }
      
      copyCellData(mpiGrid[closestCell],mpiGrid[cellID], copyMomentsOnly, popID, calculate_V_moments);
   }
   
   /*! Function used to average and copy the distribution and moments from all the closest sysboundarytype::NOT_SYSBOUNDARY cells.
    * \param mpiGrid Grid
    * \param cellID The cell's ID.
    */
   void SysBoundaryCondition::vlasovBoundaryCopyFromAllClosestNbrs(
      const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const CellID& cellID,const uint popID, const bool calculate_V_moments
   ) {
      const std::vector<CellID> closestCells = getAllClosestNonsysboundaryCells(cellID);
      
      if(closestCells[0] == INVALID_CELLID) {
         cerr << __FILE__ << ":" << __LINE__ << ": No closest cell found!" << endl;
         abort();
      }
      averageCellData(mpiGrid, closestCells, mpiGrid[cellID], popID, calculate_V_moments);
   }
   
   /*! Function used to average and copy the distribution and moments from all the close sysboundarytype::NOT_SYSBOUNDARY cells.
    * \param mpiGrid Grid
    * \param cellID The cell's ID.
    */
   void SysBoundaryCondition::vlasovBoundaryFluffyCopyFromAllCloseNbrs(
      const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const CellID& cellID,const uint popID,const bool calculate_V_moments, creal fluffiness
   ) {
      const std::vector<CellID> closeCells = getAllCloseNonsysboundaryCells(cellID);
      
      if(closeCells[0] == INVALID_CELLID) {
         cerr << __FILE__ << ":" << __LINE__ << ": No close cell found!" << endl;
         abort();
      }
      averageCellData(mpiGrid, closeCells, mpiGrid[cellID], popID, calculate_V_moments, fluffiness);
   }
   
   /*! Function used to copy the distribution from (one of) the closest sysboundarytype::NOT_SYSBOUNDARY cell but limiting to values no higher than where it can flow into. Moments are recomputed.
    * \param mpiGrid Grid
    * \param cellID The cell's ID.
    */
   void SysBoundaryCondition::vlasovBoundaryCopyFromTheClosestNbrAndLimit(
      const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const CellID& cellID,
      const uint popID
      ) {
      const CellID closestCell = getTheClosestNonsysboundaryCell(cellID);
      SpatialCell * from = mpiGrid[closestCell];
      SpatialCell * to = mpiGrid[cellID];
      
      if(closestCell == INVALID_CELLID) {
         cerr << __FILE__ << ":" << __LINE__ << ": No closest cell found!" << endl;
         abort();
      }
      
      const std::array<SpatialCell*,27> flowtoCells = getFlowtoCells(cellID);
      //Do not allow block adjustment, the block structure when calling vlasovBoundaryCondition should be static
      //just copy data to existing blocks, no modification of to blocks allowed
      for (vmesh::LocalID blockLID=0; blockLID<to->get_number_of_velocity_blocks(popID); ++blockLID) {
         const vmesh::GlobalID blockGID = to->get_velocity_block_global_id(blockLID,popID);
//          const Realf* fromBlock_data = from->get_data(from->get_velocity_block_local_id(blockGID) );
         Realf* toBlock_data = to->get_data(blockLID,popID);
         if (from->get_velocity_block_local_id(blockGID,popID) == from->invalid_local_id()) {
            for (unsigned int i = 0; i < VELOCITY_BLOCK_LENGTH; i++) {
               toBlock_data[i] = 0.0; //block did not exist in from cell, fill with zeros.
            }
         } else {
            const Real* blockParameters = to->get_block_parameters(blockLID, popID);
            // check where cells are
            creal vxBlock = blockParameters[BlockParams::VXCRD];
            creal vyBlock = blockParameters[BlockParams::VYCRD];
            creal vzBlock = blockParameters[BlockParams::VZCRD];
            creal dvxCell = blockParameters[BlockParams::DVX];
            creal dvyCell = blockParameters[BlockParams::DVY];
            creal dvzCell = blockParameters[BlockParams::DVZ];
            
            std::array<Realf*,27> flowtoCellsBlockCache = getFlowtoCellsBlock(flowtoCells, blockGID, popID);
            
            for (uint kc=0; kc<WID; ++kc) {
               for (uint jc=0; jc<WID; ++jc) {
                  for (uint ic=0; ic<WID; ++ic) {
                     velocity_cell_indices_t indices = {ic, jc, kc};
                     const uint cell = from->get_velocity_cell(indices);
                     
                     creal vxCellCenter = vxBlock + (ic+convert<Real>(0.5))*dvxCell;
                     creal vyCellCenter = vyBlock + (jc+convert<Real>(0.5))*dvyCell;
                     creal vzCellCenter = vzBlock + (kc+convert<Real>(0.5))*dvzCell;
                     const int vxCellSign = vxCellCenter < 0 ? -1 : 1;
                     const int vyCellSign = vyCellCenter < 0 ? -1 : 1;
                     const int vzCellSign = vzCellCenter < 0 ? -1 : 1;
                     Realf value = from->get_value(blockGID, cell, popID);
                     //loop over spatial cells in quadrant of influence
                     for(int dvx = 0 ; dvx <= 1; dvx++) {
                        for(int dvy = 0 ; dvy <= 1; dvy++) {
                           for(int dvz = 0 ; dvz <= 1; dvz++) {
                              const int flowToId = nbrID(dvx * vxCellSign, dvy * vyCellSign, dvz * vzCellSign);
                              if(flowtoCells.at(flowToId)){
                                 value = min(value, flowtoCellsBlockCache.at(flowToId)[cell]);
                              }
                           }
                        }
                     }
                     to->set_value(blockGID, cell,  value, popID);
                  }
               }
            }
         }
      }
   }
   
   /*! Function used to copy the distribution and moments from one cell to another.
    * \param from Pointer to parent cell to copy from.
    * \param to Pointer to destination cell.
    */
   void SysBoundaryCondition::copyCellData(
            SpatialCell* from,
            SpatialCell* to,
            const bool copyMomentsOnly,
            const uint popID,
            const bool calculate_V_moments
   ) {
      if (popID == 0) {
         if (calculate_V_moments) {
            to->parameters[CellParams::RHOM_V] = from->parameters[CellParams::RHOM_V];
            to->parameters[CellParams::VX_V] = from->parameters[CellParams::VX_V];
            to->parameters[CellParams::VY_V] = from->parameters[CellParams::VY_V];
            to->parameters[CellParams::VZ_V] = from->parameters[CellParams::VZ_V];
            to->parameters[CellParams::RHOQ_V] = from->parameters[CellParams::RHOQ_V];
            to->parameters[CellParams::P_11_V] = from->parameters[CellParams::P_11_V];
            to->parameters[CellParams::P_22_V] = from->parameters[CellParams::P_22_V];
            to->parameters[CellParams::P_33_V] = from->parameters[CellParams::P_33_V];
         } else {
            to->parameters[CellParams::RHOM_R] = from->parameters[CellParams::RHOM_R];
            to->parameters[CellParams::VX_R] = from->parameters[CellParams::VX_R];
            to->parameters[CellParams::VY_R] = from->parameters[CellParams::VY_R];
            to->parameters[CellParams::VZ_R] = from->parameters[CellParams::VZ_R];
            to->parameters[CellParams::RHOQ_R] = from->parameters[CellParams::RHOQ_R];
            to->parameters[CellParams::P_11_R] = from->parameters[CellParams::P_11_R];
            to->parameters[CellParams::P_22_R] = from->parameters[CellParams::P_22_R];
            to->parameters[CellParams::P_33_R] = from->parameters[CellParams::P_33_R];
         }
      }
      
       if(!copyMomentsOnly) { // Do this only if copyMomentsOnly is false.
         to->set_population(from->get_population(popID), popID);
      } else {
         if (calculate_V_moments) {
            to->get_population(popID).RHO_V = from->get_population(popID).RHO_V;
         } else {
            to->get_population(popID).RHO_R = from->get_population(popID).RHO_R;
         }
         
         for (uint i=0; i<3; i++) {
            if (calculate_V_moments) {
               to->get_population(popID).V_V[i] = from->get_population(popID).V_V[i];
               to->get_population(popID).P_V[i] = from->get_population(popID).P_V[i];
            } else {
               to->get_population(popID).V_R[i] = from->get_population(popID).V_R[i];
               to->get_population(popID).P_R[i] = from->get_population(popID).P_R[i];
            }
         }
      }
   }
   
   /*! Take a list of cells and set the destination cell distribution function to the average of the list's cells'.
    * \param mpiGrid Grid
    * \param cellList Vector of cells to copy from.
    * \param to Pointer to cell in which to set the averaged distribution.
    */
   void SysBoundaryCondition::averageCellData(
         const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
         const std::vector<CellID> cellList,
         SpatialCell *to,
         const uint popID,
         const bool calculate_V_moments,
         creal fluffiness /* default =0.0*/
   ) {
      const size_t numberOfCells = cellList.size();
      creal factor = fluffiness / convert<Real>(numberOfCells);
      

      // Rescale own vspace
      const Realf* toData = to->get_data(popID);
      for (vmesh::LocalID toBlockLID=0; toBlockLID<to->get_number_of_velocity_blocks(popID); ++toBlockLID) {
         // Pointer to target block data
         Realf* toData = to->get_data(toBlockLID,popID);
         
         // Add values from source cells
         for (uint kc=0; kc<WID; ++kc) for (uint jc=0; jc<WID; ++jc) for (uint ic=0; ic<WID; ++ic) {
            toData[cellIndex(ic,jc,kc)] *= 1.0 - fluffiness;
         }
         toData += SIZE_VELBLOCK;
      } // for-loop over velocity blocks

      
      for (size_t i=0; i<numberOfCells; i++) {
         const SpatialCell* incomingCell = mpiGrid[cellList[i]];

         const Realf* fromData = incomingCell->get_data(popID);
         for (vmesh::LocalID incBlockLID=0; incBlockLID<incomingCell->get_number_of_velocity_blocks(popID); ++incBlockLID) {
            // Global ID of the block containing incoming data
            vmesh::GlobalID incBlockGID = incomingCell->get_velocity_block_global_id(incBlockLID,popID);
            
            // Get local ID of the target block. If the block doesn't exist, create it.
            vmesh::GlobalID toBlockLID = to->get_velocity_block_local_id(incBlockGID,popID);
            if (toBlockLID == SpatialCell::invalid_local_id()) {
               to->add_velocity_block(incBlockGID,popID);
               toBlockLID = to->get_velocity_block_local_id(incBlockGID,popID);
            }
            
            // Pointer to target block data
            Realf* toData = to->get_data(toBlockLID,popID);

            // Add values from source cells
            for (uint kc=0; kc<WID; ++kc) for (uint jc=0; jc<WID; ++jc) for (uint ic=0; ic<WID; ++ic) {
               toData[cellIndex(ic,jc,kc)] += factor*fromData[cellIndex(ic,jc,kc)];
            }
            fromData += SIZE_VELBLOCK;
         } // for-loop over velocity blocks
      }
   }

   /*! Take neighboring distribution and reflect all parts going in the direction opposite to the normal vector given in.
    * \param mpiGrid Grid
    * \param cellID Cell in which to set the distribution where incoming velocity cells have been reflected/bounced.
    * \param nx Unit vector x component normal to the bounce/reflection plane.
    * \param ny Unit vector y component normal to the bounce/reflection plane.
    * \param nz Unit vector z component normal to the bounce/reflection plane.
    */
   void SysBoundaryCondition::vlasovBoundaryReflect(
         const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
         const CellID& cellID,
         creal& nx,
         creal& ny,
         creal& nz,
         const uint popID
   ) {
      SpatialCell * cell = mpiGrid[cellID];
      const std::vector<CellID> cellList = this->getAllClosestNonsysboundaryCells(cellID);
      const size_t numberOfCells = cellList.size();

      creal factor = 1.0 / convert<Real>(numberOfCells);
      
      cell->clear(popID);
      
      for (size_t i=0; i<numberOfCells; i++) {
         SpatialCell* incomingCell = mpiGrid[cellList[i]];
         const Real* blockParameters = incomingCell->get_block_parameters(popID);

         // add blocks
         for (vmesh::LocalID blockLID=0; blockLID<incomingCell->get_number_of_velocity_blocks(popID); ++blockLID) {
            // check where cells are
            creal vxBlock = blockParameters[BlockParams::VXCRD];
            creal vyBlock = blockParameters[BlockParams::VYCRD];
            creal vzBlock = blockParameters[BlockParams::VZCRD];
            creal dvxCell = blockParameters[BlockParams::DVX];
            creal dvyCell = blockParameters[BlockParams::DVY];
            creal dvzCell = blockParameters[BlockParams::DVZ];
            for (uint kc=0; kc<WID; ++kc) for (uint jc=0; jc<WID; ++jc) for (uint ic=0; ic<WID; ++ic) {
               creal vxCellCenter = vxBlock + (ic+convert<Real>(0.5))*dvxCell;
               creal vyCellCenter = vyBlock + (jc+convert<Real>(0.5))*dvyCell;
               creal vzCellCenter = vzBlock + (kc+convert<Real>(0.5))*dvzCell;
               // scalar product v.n
               creal vNormal = vxCellCenter*nx + vyCellCenter*ny + vzCellCenter*nz;
               if (vNormal >= 0.0) {
                  // Not flowing in, leave as is.
                  cell->increment_value(
                     vxCellCenter,
                     vyCellCenter,
                     vzCellCenter,
                     factor*incomingCell->get_value(vxCellCenter, vyCellCenter, vzCellCenter,popID),
                     popID
                  );
               } else {
                  // Flowing in, bounce off.
                  cell->increment_value(
                     vxCellCenter - 2.0*vNormal*nx,
                     vyCellCenter - 2.0*vNormal*ny,
                     vzCellCenter - 2.0*vNormal*nz,
                     factor*incomingCell->get_value(vxCellCenter, vyCellCenter, vzCellCenter,popID),
                     popID
                  );
               }
            } // for-loop over cells in velocity block
         } // for-loop over velocity blocks
         blockParameters += BlockParams::N_VELOCITY_BLOCK_PARAMS;
      } // for-loop over spatial cells
   }
   
   /*! Take neighboring distribution and absorb all parts going in the direction opposite to the normal vector given in.
    * \param mpiGrid Grid
    * \param cellID Cell in which to set the distribution where incoming velocity cells have been kept or swallowed.
    * \param nx Unit vector x component normal to the absorption plane.
    * \param ny Unit vector y component normal to the absorption plane.
    * \param nz Unit vector z component normal to the absorption plane.
    * \param quenchingFactor Multiplicative factor by which to scale the distribution function values. 0: absorb. ]0;1[: quench.
    */
   void SysBoundaryCondition::vlasovBoundaryAbsorb(
      const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const CellID& cellID,
      creal& nx,
      creal& ny,
      creal& nz,
      creal& quenchingFactor,
      const uint popID
   ) {
      SpatialCell* cell = mpiGrid[cellID];
      const std::vector<CellID> cellList = this->getAllClosestNonsysboundaryCells(cellID);
      const size_t numberOfCells = cellList.size();

      creal factor = 1.0 / convert<Real>(numberOfCells);
      
      cell->clear(popID);
      
      for (size_t i=0; i<numberOfCells; i++) {
         SpatialCell* incomingCell = mpiGrid[cellList[i]];
         const Real* blockParameters = incomingCell->get_block_parameters(popID);
         
         // add blocks
         for (vmesh::LocalID blockLID=0; blockLID<incomingCell->get_number_of_velocity_blocks(popID); ++blockLID) {
            // check where cells are
            creal vxBlock = blockParameters[BlockParams::VXCRD];
            creal vyBlock = blockParameters[BlockParams::VYCRD];
            creal vzBlock = blockParameters[BlockParams::VZCRD];
            creal dvxCell = blockParameters[BlockParams::DVX];
            creal dvyCell = blockParameters[BlockParams::DVY];
            creal dvzCell = blockParameters[BlockParams::DVZ];
            for (uint kc=0; kc<WID; ++kc) 
               for (uint jc=0; jc<WID; ++jc) 
                  for (uint ic=0; ic<WID; ++ic) {
                     creal vxCellCenter = vxBlock + (ic+convert<Real>(0.5))*dvxCell;
                     creal vyCellCenter = vyBlock + (jc+convert<Real>(0.5))*dvyCell;
                     creal vzCellCenter = vzBlock + (kc+convert<Real>(0.5))*dvzCell;
                     // scalar product v.n
                     creal vNormal = vxCellCenter*nx + vyCellCenter*ny + vzCellCenter*nz;
                     if(vNormal >= 0.0) {
                        // Not flowing in, leave as is.
                        cell->increment_value(
                           vxCellCenter,
                           vyCellCenter,
                           vzCellCenter,
                           factor*incomingCell->get_value(vxCellCenter, vyCellCenter, vzCellCenter),
                           popID
                        );
                     } else {
                        // Flowing in, bounce off.
                        cell->increment_value(
                           vxCellCenter,
                           vyCellCenter,
                           vzCellCenter,
                           factor*quenchingFactor*incomingCell->get_value(vxCellCenter, vyCellCenter, vzCellCenter),
                           popID
                        );
                     }
            }
         } // for-loop over velocity blocks
         blockParameters += BlockParams::N_VELOCITY_BLOCK_PARAMS;
      } // for-loop over spatial cells
   }


   /*! Updates the system boundary conditions after load balancing. This is called from e.g. the class SysBoundary.
    * \param mpiGrid Grid
    * \param local_cells_on_boundary Cells within this process
    * \retval success Returns true if the operation is successful
    */
   bool SysBoundaryCondition::updateSysBoundaryConditionsAfterLoadBalance(
      dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const vector<CellID> & local_cells_on_boundary
   ) {
      // Loop over cellids
      for( vector<CellID>::const_iterator it = local_cells_on_boundary.begin(); it != local_cells_on_boundary.end(); ++it ) {
         const CellID cellId = *it;
         std::vector<CellID> & closestCells = allClosestNonsysboundaryCells[cellId];
         closestCells.clear();
         std::vector<CellID> & closeCells = allCloseNonsysboundaryCells[cellId];
         closeCells.clear();
         std::array<SpatialCell*,27> & flowtoCells = allFlowtoCells[cellId];
         flowtoCells.fill(NULL);
         uint dist = numeric_limits<uint>::max();
      
         // First iteration of search to determine closest distance
         for(int i=-2; i<3; i++)
            for(int j=-2; j<3; j++)
               for(int k=-2; k<3; k++) {
                  const CellID cell = getNeighbour(mpiGrid,cellId,i,j,k);
                  if(cell != INVALID_CELLID) {
                     if(mpiGrid[cell]->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
                        cuint d2 = i*i+j*j+k*k;
                        if(d2 < dist) {
                           dist = d2;
                        }
                        // Flowto neighbours have distances of 1, 2 or 3 at a distance of 1 layer, 4, 5 or 6 at a distance of 2 layers.
                        // Furthermore one does not want to have the cell itself in this list.
                        if(d2 < 4 && i != 0 && j != 0 && k != 0) {
                           flowtoCells.at(i + 3*j + 9*k + 13) = mpiGrid[cell];
                        }
                        if(mpiGrid[cellId]->sysBoundaryLayer == 1 && abs(i) < 2 && abs(j) < 2 && abs(k) < 2) {
                           closeCells.push_back(cell);
                        }
                        if(mpiGrid[cellId]->sysBoundaryLayer == 2) {
                           closeCells.push_back(cell);
                        }
                     }
                  }
               }
         // Second iteration to record the cellIds of all cells at closest distance
         for(int i=-2; i<3; i++)
            for(int j=-2; j<3; j++)
               for(int k=-2; k<3; k++) {
                  const CellID cell = getNeighbour(mpiGrid,cellId,i,j,k);
                  if(cell != INVALID_CELLID) {
                     if(mpiGrid[cell]->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
                        cuint d2 = i*i+j*j+k*k;
                        if(d2 == dist) {
                           closestCells.push_back(cell);
                        }
                     }
                  }
               }
         if(closestCells.size() == 0) closestCells.push_back(INVALID_CELLID);
         if(closeCells.size() == 0) closeCells.push_back(INVALID_CELLID);
      }
      return true;
   }
   
   /*! Get the cellID of the first closest cell of type NOT_SYSBOUNDARY found.
    * \param i,j,k Coordinates of the cell to start looking from
    * \return The cell index of that cell
    * \sa getAllClosestNonsysboundaryCells
    */
   std::array<int, 3> SysBoundaryCondition::getTheClosestNonsysboundaryCell(
      FsGrid< fsgrids::technical, 2> & technicalGrid,
      cint i,
      cint j,
      cint k
   ) {
      const std::vector< std::array<int, 3> > closestCells = getAllClosestNonsysboundaryCells(technicalGrid, i, j, k);
      return closestCells.at(0);
   }
   
   /*! Get the cellIDs of all the closest cells of type NOT_SYSBOUNDARY.
    * \param i,j,k Coordinates of the cell to start looking from
    * \return The vector of cell indices of those cells
    * \sa getTheClosestNonsysboundaryCell
    */
   std::vector< std::array<int, 3> > SysBoundaryCondition::getAllClosestNonsysboundaryCells(
      FsGrid< fsgrids::technical, 2> & technicalGrid,
      cint i,
      cint j,
      cint k
   ) {
      int distance = std::numeric_limits<int>::max();
      std::vector< std::array<int,3> > closestCells;
      
      for (int kk=-2; kk<3; kk++) {
         for (int jj=-2; jj<3; jj++) {
            for (int ii=-2; ii<3 ; ii++) {
               if( technicalGrid.get(i+ii,j+jj,k+kk) && technicalGrid.get(i+ii,j+jj,k+kk)->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
                  distance = min(distance, ii*ii + jj*jj + kk*kk);
               }
            }
         }
      }
      
      for (int kk=-2; kk<3; kk++) {
         for (int jj=-2; jj<3; jj++) {
            for (int ii=-2; ii<3 ; ii++) {
               if( technicalGrid.get(i+ii,j+jj,k+kk) && technicalGrid.get(i+ii,j+jj,k+kk)->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
                  int d = ii*ii + jj*jj + kk*kk;
                  if( d == distance ) {
                     std::array<int, 3> cell = {i+ii, j+jj, k+kk};
                     closestCells.push_back(cell);
                  }
               }
            }
         }
      }
      
      if(closestCells.size() == 0) {
         std::array<int, 3> dummy  = {std::numeric_limits<int>::min()};
         closestCells.push_back(dummy);
      }
      
      return closestCells;
   }
   
   /*! Get the cellID of the first closest cell of type NOT_SYSBOUNDARY found.
    * \param cellID ID of the cell to start look from.
    * \return The cell index of that cell
    * \sa getAllClosestNonsysboundaryCells
    */
   CellID & SysBoundaryCondition::getTheClosestNonsysboundaryCell(
      const CellID& cellID
   ) {
      std::vector<CellID> & closestCells = allClosestNonsysboundaryCells.at(cellID);
      return closestCells.at(0);
   }
   
   /*! Get the cellIDs of all the closest cells of type NOT_SYSBOUNDARY.
    * \param cellID ID of the cell to start look from.
    * \return The vector of cell indices of those cells
    * \sa getTheClosestNonsysboundaryCell
    */
   std::vector<CellID> & SysBoundaryCondition::getAllClosestNonsysboundaryCells(
      const CellID& cellID
   ) {
      std::vector<CellID> & closestCells = allClosestNonsysboundaryCells.at(cellID);
      return closestCells;
   }
   
   /*! Get the cellIDs of all the close cells of type NOT_SYSBOUNDARY.
    * \param cellID ID of the cell to start look from.
    * \return The vector of cell indices of those cells
    */
   std::vector<CellID> & SysBoundaryCondition::getAllCloseNonsysboundaryCells(
      const CellID& cellID
   ) {
      std::vector<CellID> & closeCells = allCloseNonsysboundaryCells.at(cellID);
      return closeCells;
   }
   
   /*! Get the cellIDs of all flowto cells (cells into which the velocity distribution can flow and which is of type NOT_SYSBOUNDARY).
    * \param cellID ID of the cell to start look from.
    * \return The vector of cell indices of those cells
    */
   std::array<SpatialCell*,27> & SysBoundaryCondition::getFlowtoCells(
      const CellID& cellID
   ) {
      phiprof::start("getFlowtoCells");
      std::array<SpatialCell*,27> & flowtoCells = allFlowtoCells.at(cellID);
      phiprof::stop("getFlowtoCells");
      return flowtoCells;
   }
   
   std::array<Realf*,27> SysBoundaryCondition::getFlowtoCellsBlock(
      const std::array<SpatialCell*,27> flowtoCells,
      const vmesh::GlobalID blockGID,
      const uint popID
   ) {
      phiprof::start("getFlowtoCellsBlock");
      std::array<Realf*,27> flowtoCellsBlock;
      flowtoCellsBlock.fill(NULL);
      for (uint i=0; i<27; i++) {
         if(flowtoCells.at(i)) {
            flowtoCellsBlock.at(i) = flowtoCells.at(i)->get_data(flowtoCells.at(i)->get_velocity_block_local_id(blockGID,popID), popID);
         }
      }
      phiprof::stop("getFlowtoCellsBlock");
      return flowtoCellsBlock;
   }
   
   Real SysBoundaryCondition::fieldBoundaryCopyFromSolvingNbrMagneticField(
      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2> & bGrid,
      FsGrid< fsgrids::technical, 2> & technicalGrid,
      cint i,
      cint j,
      cint k,
      cuint component,
      cuint mask
   ) {

      int distance = std::numeric_limits<int>::max();
      std::vector< std::array<int,3> > closestCells;

      for (int kk=-2; kk<3; kk++) {
         for (int jj=-2; jj<3; jj++) {
            for (int ii=-2; ii<3 ; ii++) {
               if( technicalGrid.get(i+ii,j+jj,k+kk) // skip invalid cells returning NULL
                   && (technicalGrid.get(i+ii,j+jj,k+kk)->SOLVE & mask) == mask // Did that guy solve this component?
                   && technicalGrid.get(i+ii,j+jj,k+kk)->sysBoundaryFlag != sysboundarytype::DO_NOT_COMPUTE // Do not copy from there
               ) {
                  distance = min(distance, ii*ii + jj*jj + kk*kk);
               }
            }
         }
      }

      for (int kk=-2; kk<3; kk++) {
         for (int jj=-2; jj<3; jj++) {
            for (int ii=-2; ii<3 ; ii++) {
               if( technicalGrid.get(i+ii,j+jj,k+kk) // skip invalid cells returning NULL
                   && (technicalGrid.get(i+ii,j+jj,k+kk)->SOLVE & mask) == mask // Did that guy solve this component?
                   && technicalGrid.get(i+ii,j+jj,k+kk)->sysBoundaryFlag != sysboundarytype::DO_NOT_COMPUTE // Do not copy from there
               ) {
                  int d = ii*ii + jj*jj + kk*kk;
                  if( d == distance ) {
                     std::array<int, 3> cell = {i+ii, j+jj, k+kk};
                     closestCells.push_back(cell);
                  }
               }
            }
         }
      }

      if(closestCells.size() == 0) {
         cerr << __FILE__ << ":" << __LINE__ << ": No closest cell found!" << endl;
         abort();
      }

      return bGrid.get(closestCells[0][0], closestCells[0][1], closestCells[0][2])->at(fsgrids::bfield::PERBX+component);
   }
   
   /*! Function used in some cases to know which faces the system boundary condition is being applied to.
    * \param faces Pointer to array of 6 bool in which the values are returned whether the corresponding face is of that type. Order: 0 x+; 1 x-; 2 y+; 3 y-; 4 z+; 5 z-
    */
   void SysBoundaryCondition::getFaces(bool* faces) {
      cerr << "ERROR: SysBoundaryCondition::getFaces called instead of derived class function!" << endl;
      for(int i=0; i<6; i++) {
        faces[i]=false;
      }
   }
   
   /*! Get the precedence value of the system boundary condition.
    * \return The precedence value of the system boundary condition as set by parameter.
    */
   uint SysBoundaryCondition::getPrecedence() const {return precedence;}
   
   /*! Returns whether the boundary condition is dynamic in time.
    * \return Boolean value.
    */
   bool SysBoundaryCondition::isDynamic() const {return isThisDynamic;}
   
   void SysBoundaryCondition::setPeriodicity(
      bool isFacePeriodic[3]
   ) {
      for (uint i=0; i<3; i++) {
         this->isPeriodic[i] = isFacePeriodic[i];
      }
   }
   
   /*! Get a bool telling whether to call again applyInitialState upon restarting the simulation. */
   bool SysBoundaryCondition::doApplyUponRestart() const {return this->applyUponRestart;}
} // namespace SBC
