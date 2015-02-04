/*
This file is part of Vlasiator.

Copyright 2010, 2011, 2012, 2013 Finnish Meteorological Institute












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
    */
   void SysBoundaryCondition::determineFace(
      bool* isThisCellOnAFace,
      creal x, creal y, creal z,
      creal dx, creal dy, creal dz
   ) {
      for(uint i=0; i<6; i++) {
         isThisCellOnAFace[i] = false;
      }
      if(x > Parameters::xmax - 2.0*dx) isThisCellOnAFace[0] = true;
      if(x < Parameters::xmin + 2.0*dx) isThisCellOnAFace[1] = true;
      if(y > Parameters::ymax - 2.0*dy) isThisCellOnAFace[2] = true;
      if(y < Parameters::ymin + 2.0*dy) isThisCellOnAFace[3] = true;
      if(z > Parameters::zmax - 2.0*dz) isThisCellOnAFace[4] = true;
      if(z < Parameters::zmin + 2.0*dz) isThisCellOnAFace[5] = true;
   }
   
   /*! SysBoundaryCondition base class constructor. The constructor is empty.*/
   SysBoundaryCondition::SysBoundaryCondition() { }
   
   /*! SysBoundaryCondition base class virtual destructor. The destructor is empty.*/
   SysBoundaryCondition::~SysBoundaryCondition() { }
   
   void SysBoundaryCondition::addParameters() {
      cerr << "ERROR: SysBoundaryCondition::addParameters called instead of derived class function!" << endl;
   }
   
   void SysBoundaryCondition::getParameters() {
      cerr << "ERROR: SysBoundaryCondition::getParameters called instead of derived class function!" << endl;
   }
   
   /*! Function called during initialisation to set the system boundary condition's parameters.
    * This function must initialize the boundary condition function for all particle species,
    * i.e., its behavior is not allowed to depend on SysBoundaryCondition::activePopID.
    */
   bool SysBoundaryCondition::initSysBoundary(
      creal& t,
      Project &project
   ) {
      cerr << "ERROR: SysBoundaryCondition::initSysBoundary called instead of derived class function!" << endl;
      return false;
   }

   /*! Function used to assign the system boundary condition type to a cell.
    * \return The system boundary condition type's index
    */
   bool SysBoundaryCondition::assignSysBoundary(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid) {
      cerr << "ERROR: SysBoundaryCondition::assignSysBoundary called instead of derived class function!" << endl;
      return false;
   }
   
   /*! Function used to apply the system boundary condition initial state to a cell. 
    * \return Boolean true on success, false on failure.
    */
   bool SysBoundaryCondition::applyInitialState(
      const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      Project &project
   ) {
      cerr << "ERROR: SysBoundaryCondition::applyInitialState called instead of derived class function!" << endl;
      return false;
   }
   
   /*! Function used to return the system boundary condition cell's magnetic field component.
    * \param mpiGrid The grid.
    * \param cellID The cell's ID.
    * \param dt 0 when at the beginning of a time step, non-zero for the _DT2 (Runge-Kutta stepping) value.
    * \param component 0: x-component, 1: y-component, 2: z-component.
    * 
    * \return The requested component value.
    */
   Real SysBoundaryCondition::fieldSolverBoundaryCondMagneticField(
      const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const CellID& cellID,
      creal& dt,
      cuint& component
   ) {
      cerr << "ERROR: SysBoundaryCondition::fieldSolverBoundaryCondMagneticField called instead of derived class function!" << endl;
      exit(1);
   }
   
   /*! Function used to return the system boundary condition cell's main electric field component (without Hall term).
    * \param mpiGrid The grid.
    * \param cellID The cell's ID.
    * \param RKCase The step in Runge-Kutta (use values from enum defined in common.h).
    * \param component 0: x-component, 1: y-component, 2: z-component.
    * 
    * \return The requested component value.
    * 
    * \sa fieldSolverBoundaryCondHallElectricField
    */
   void SysBoundaryCondition::fieldSolverBoundaryCondElectricField(
      dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const CellID& cellID,
      cuint RKCase,
      cuint component
   ) {
      cerr << "ERROR: SysBoundaryCondition::fieldSolverBoundaryCondElectricField called instead of derived class function!" << endl;
      exit(1);
   }
   
   /*! Function used to compute the system boundary condition cell's Hall electric field components and save them into the cell parameters.
    * \param mpiGrid The grid.
    * \param cellID The cell's ID.
    * \param RKCase The step in Runge-Kutta (use values from enum defined in common.h).
    * \param component 0: x-component, 1: y-component, 2: z-component.
    */
   void SysBoundaryCondition::fieldSolverBoundaryCondHallElectricField(
      dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid,
      const CellID& cellID,
      cuint RKCase,
      cuint component
   ) {
      cerr << "ERROR: SysBoundaryCondition::fieldSolverBoundaryCondHallElectricField called instead of derived class function!" << endl;
      exit(1);
   }
   
   /*! Function used to compute the system boundary condition cell's derivatives and save them into the cell derivatives.
    * \param mpiGrid The grid.
    * \param cellID The cell's ID.
    * \param RKCase The step in Runge-Kutta (use values from enum defined in common.h).
    * \param component 0: x-derivatives, 1: y-derivatives, 2: z-derivatives, 3: xy-derivatives, 4: xz-derivatives, 5: yz-derivatives.
    */
   void SysBoundaryCondition::fieldSolverBoundaryCondDerivatives(
      dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const CellID& cellID,
      cuint& RKCase,
      cuint& component
   ) {
      cerr << "ERROR: SysBoundaryCondition::fieldSolverBoundaryCondDerivatives called instead of derived class function!" << endl;
      exit(1);
   }
   
   /*! Function used to compute the system boundary condition cell's BVOL derivatives and save them into the cell derivatives.
    * \param mpiGrid The grid.
    * \param cellID The cell's ID.
    * \param component 0: x-derivatives, 1: y-derivatives, 2: z-derivatives.
    */
   void SysBoundaryCondition::fieldSolverBoundaryCondBVOLDerivatives(
      const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const CellID& cellID,
      cuint& component
   ) {
      cerr << "ERROR: SysBoundaryCondition::fieldSolverBoundaryCondBVOLDerivatives called instead of derived class function!" << endl;
      exit(1);
   }
   
   /*! Function used to set the system boundary condition cell's derivatives to 0.
    * \param mpiGrid The grid.
    * \param cellID The cell's ID.
    * \param component 0: x-derivatives, 1: y-derivatives, 2: z-derivatives, 3: xy-derivatives, 4: xz-derivatives, 5: yz-derivatives.
    */
   void SysBoundaryCondition::setCellDerivativesToZero(
      const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const CellID& cellID,
      cuint& component
   ) {
      Real* const derivs = &(mpiGrid[cellID]->derivatives[0]);
      switch(component) {
         case 0: // x, xx
            derivs[fieldsolver::drhodx] = 0.0;
            derivs[fieldsolver::dp11dx] = 0.0;
            derivs[fieldsolver::dp22dx] = 0.0;
            derivs[fieldsolver::dp33dx] = 0.0;
            derivs[fieldsolver::dPERBydx]  = 0.0;
            derivs[fieldsolver::dPERBzdx]  = 0.0;
            derivs[fieldsolver::dVxdx]  = 0.0;
            derivs[fieldsolver::dVydx]  = 0.0;
            derivs[fieldsolver::dVzdx]  = 0.0;
            derivs[fieldsolver::dPERBydxx] = 0.0;
            derivs[fieldsolver::dPERBzdxx] = 0.0;
            break;
         case 1: // y, yy
            derivs[fieldsolver::drhody] = 0.0;
            derivs[fieldsolver::dp11dy] = 0.0;
            derivs[fieldsolver::dp22dy] = 0.0;
            derivs[fieldsolver::dp33dy] = 0.0;
            derivs[fieldsolver::dPERBxdy]  = 0.0;
            derivs[fieldsolver::dPERBzdy]  = 0.0;
            derivs[fieldsolver::dVxdy]  = 0.0;
            derivs[fieldsolver::dVydy]  = 0.0;
            derivs[fieldsolver::dVzdy]  = 0.0;
            derivs[fieldsolver::dPERBxdyy] = 0.0;
            derivs[fieldsolver::dPERBzdyy] = 0.0;
            break;
         case 2: // z, zz
            derivs[fieldsolver::drhodz] = 0.0;
            derivs[fieldsolver::dp11dz] = 0.0;
            derivs[fieldsolver::dp22dz] = 0.0;
            derivs[fieldsolver::dp33dz] = 0.0;
            derivs[fieldsolver::dPERBxdz]  = 0.0;
            derivs[fieldsolver::dPERBydz]  = 0.0;
            derivs[fieldsolver::dVxdz]  = 0.0;
            derivs[fieldsolver::dVydz]  = 0.0;
            derivs[fieldsolver::dVzdz]  = 0.0;
            derivs[fieldsolver::dPERBxdzz] = 0.0;
            derivs[fieldsolver::dPERBydzz] = 0.0;
            break;
         case 3: // xy
            derivs[fieldsolver::dPERBzdxy] = 0.0;
            break;
         case 4: // xz
            derivs[fieldsolver::dPERBydxz] = 0.0;
            break;
         case 5: // yz
            derivs[fieldsolver::dPERBxdyz] = 0.0;
            break;
         default:
            cerr << __FILE__ << ":" << __LINE__ << ":" << " Invalid component" << endl;
      }
   }
   
   /*! Function used to set the system boundary condition cell's BVOL derivatives to 0.
    * \param mpiGrid The grid.
    * \param cellID The cell's ID.
    * \param component 0: x-derivatives, 1: y-derivatives, 2: z-derivatives.
    */
   void SysBoundaryCondition::setCellBVOLDerivativesToZero(
      const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const CellID& cellID,
      cuint& component
   ) {
      Real* const derivs = &(mpiGrid[cellID]->derivativesBVOL[0]);
      switch(component) {
         case 0:
            derivs[bvolderivatives::dPERBYVOLdx] = 0.0;
            derivs[bvolderivatives::dPERBZVOLdx] = 0.0;
            break;
         case 1:
            derivs[bvolderivatives::dPERBXVOLdy] = 0.0;
            derivs[bvolderivatives::dPERBZVOLdy] = 0.0;
            break;
         case 2:
            derivs[bvolderivatives::dPERBXVOLdz] = 0.0;
            derivs[bvolderivatives::dPERBYVOLdz] = 0.0;
            break;
         default:
            cerr << __FILE__ << ":" << __LINE__ << ":" << " Invalid component" << endl;
      }
   }
   
   /*! Function used to compute the system boundary condition cell's distribution function and moments.
    * \param mpiGrid The grid.
    * \param cellID The cell's ID.
    */
   void SysBoundaryCondition::vlasovBoundaryCondition(
      const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const CellID& cellID
   ) {
      cerr << "ERROR: SysBoundaryCondition::vlasovBoundaryCondition called instead of derived class function!" << endl;
      exit(1);
   }
   
   /*! Function used to copy the distribution and moments from (one of) the closest sysboundarytype::NOT_SYSBOUNDARY cell.
    * \param mpiGrid The grid.
    * \param cellID The cell's ID.
    */
   void SysBoundaryCondition::vlasovBoundaryCopyFromTheClosestNbr(
      const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const CellID& cellID
   ) {
      const CellID closestCell = getTheClosestNonsysboundaryCell(mpiGrid, cellID);
      
      if(closestCell == INVALID_CELLID) {
         cerr << __FILE__ << ":" << __LINE__ << ": No closest cell found!" << endl;
         abort();
      }
      //Do not allow block adjustment, the block structure when calling vlasovBoundaryCondition should be static
      copyCellData(mpiGrid[closestCell], mpiGrid[cellID],false);
   }
   
   /*! Function used to average and copy the distribution and moments from all the closest sysboundarytype::NOT_SYSBOUNDARY cells.
    * \param mpiGrid The grid.
    * \param cellID The cell's ID.
    */
   void SysBoundaryCondition::vlasovBoundaryCopyFromAllClosestNbrs(
      const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const CellID& cellID
   ) {
      const std::vector<CellID> closestCells = getAllClosestNonsysboundaryCells(mpiGrid, cellID);
      
      if(closestCells[0] == INVALID_CELLID) {
         cerr << __FILE__ << ":" << __LINE__ << ": No closest cell found!" << endl;
         abort();
      }
      averageCellData(mpiGrid, closestCells, mpiGrid[cellID]);
   }
   
   /*! Function used to copy the distribution and moments from one cell to another. In layer 2, copy only the moments.
    * \param from Pointer to parent cell to copy from.
    * \param to Pointer to destination cell.
    * \param allowBlockAdjustment If true, blocks can be created or destroyed. If false, only blocks existing in the destination cell are copied.
    */
   void SysBoundaryCondition::copyCellData(
      SpatialCell *from,
      SpatialCell *to,
      bool allowBlockAdjustment
   ) {      
      // WARNING Time-independence assumed here. _R and _V not copied, 
      // as boundary conditions cells should not set/use them.
      to->parameters[CellParams::RHO_DT2] = from->parameters[CellParams::RHO_DT2];
      to->parameters[CellParams::RHOVX_DT2] = from->parameters[CellParams::RHOVX_DT2];
      to->parameters[CellParams::RHOVY_DT2] = from->parameters[CellParams::RHOVY_DT2];
      to->parameters[CellParams::RHOVZ_DT2] = from->parameters[CellParams::RHOVZ_DT2];
      to->parameters[CellParams::P_11_DT2] = from->parameters[CellParams::P_11_DT2];
      to->parameters[CellParams::P_22_DT2] = from->parameters[CellParams::P_22_DT2];
      to->parameters[CellParams::P_33_DT2] = from->parameters[CellParams::P_33_DT2];
      to->parameters[CellParams::RHO] = from->parameters[CellParams::RHO];
      to->parameters[CellParams::RHOVX] = from->parameters[CellParams::RHOVX];
      to->parameters[CellParams::RHOVY] = from->parameters[CellParams::RHOVY];
      to->parameters[CellParams::RHOVZ] = from->parameters[CellParams::RHOVZ];
      to->parameters[CellParams::P_11] = from->parameters[CellParams::P_11];
      to->parameters[CellParams::P_22] = from->parameters[CellParams::P_22];
      to->parameters[CellParams::P_33] = from->parameters[CellParams::P_33];
      
      // Do this only for the first layer, the other layers do not need this.
      if (to->sysBoundaryLayer != 1) return;

      if (allowBlockAdjustment) {
         // prepare list of blocks to remove. It is not safe to loop over velocity_block_list while adding/removing blocks
         std::vector<vmesh::GlobalID> blocksToRemove;
         for (vmesh::LocalID block_i=0; block_i<to->get_number_of_velocity_blocks(); ++block_i) {
            const vmesh::GlobalID blockGID = to->get_velocity_block_global_id(block_i);
            
            // If this block does not exist in from, mark it for removal.
            if (from->get_velocity_block_local_id(blockGID) == from->invalid_local_id()) {
               blocksToRemove.push_back(blockGID);
            }
         }

         // remove blocks
         for (size_t b=0; b<blocksToRemove.size(); ++b) {
            cuint blockID=blocksToRemove[b];
            to->remove_velocity_block(blockID);
         }

         // add blocks
         const Realf* fromBlock_data = from->get_data();
         Realf* toBlock_data = to->get_data();
         for (vmesh::LocalID block_i=0; block_i<from->get_number_of_velocity_blocks(); ++block_i) {
            const vmesh::GlobalID blockGID = from->get_velocity_block_global_id(block_i);
            to->add_velocity_block(blockGID);
            const vmesh::LocalID toBlockLID = to->get_velocity_block_local_id(blockGID);

            for (unsigned int i = 0; i < SIZE_VELBLOCK; i++) {
               toBlock_data[toBlockLID*SIZE_VELBLOCK+i] = fromBlock_data[block_i*SIZE_VELBLOCK+i];
            }
         }
      } else {
         //just copy data to existing blocks, no modification of to blocks allowed
         const Realf* fromBlock_data = from->get_data();
         Realf* toBlock_data = to->get_data();
         for (vmesh::LocalID block_i=0; block_i<to->get_number_of_velocity_blocks(); ++block_i) {
            const vmesh::GlobalID blockGID = to->get_velocity_block_global_id(block_i);
            const vmesh::LocalID fromBlockLID = from->get_velocity_block_local_id(blockGID);
            if (from->get_velocity_block_local_id(blockGID) == from->invalid_local_id()) {
               for (unsigned int i = 0; i < VELOCITY_BLOCK_LENGTH; i++) {
                  toBlock_data[block_i*SIZE_VELBLOCK+i] = 0.0; //block did not exist in from cell, fill with zeros.
               }
            } else {
               for (unsigned int i = 0; i < VELOCITY_BLOCK_LENGTH; i++) {
                  toBlock_data[block_i*SIZE_VELBLOCK+i] = fromBlock_data[fromBlockLID*SIZE_VELBLOCK+i];
               }
            }
         }
      }
   }
   
   /*! Take a list of cells and set the destination cell distribution function to the average of the list's cells'.
    *  For layer 1 the whole distribution function is copied.
    *  For layer >1, only moments are copied
    * \param mpiGrid Grid
    * \param cellList List of cells to copy from.
    * \param to Cell in which to set the averaged distribution.
    */
   void SysBoundaryCondition::averageCellData(
      const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const std::vector<CellID> cellList,
      SpatialCell *to
   ) {
      const size_t numberOfCells = cellList.size();
      if(numberOfCells == 1) {
         copyCellData(mpiGrid[cellList[0]], to, true);
      } else {
         creal factor = 1.0 / convert<Real>(numberOfCells);
         
         to->parameters[CellParams::RHO_DT2] = 0.0;
         to->parameters[CellParams::RHOVX_DT2] = 0.0;
         to->parameters[CellParams::RHOVY_DT2] = 0.0;
         to->parameters[CellParams::RHOVZ_DT2] = 0.0;
         to->parameters[CellParams::P_11_DT2] = 0.0;
         to->parameters[CellParams::P_22_DT2] = 0.0;
         to->parameters[CellParams::P_33_DT2] = 0.0;
         to->parameters[CellParams::RHO] = 0.0;
         to->parameters[CellParams::RHOVX] = 0.0;
         to->parameters[CellParams::RHOVY] = 0.0;
         to->parameters[CellParams::RHOVZ] = 0.0;
         to->parameters[CellParams::P_11] = 0.0;
         to->parameters[CellParams::P_22] = 0.0;
         to->parameters[CellParams::P_33] = 0.0;
         
         to->clear();
         
         for (size_t i=0; i<numberOfCells; i++) {
            const SpatialCell * incomingCell = mpiGrid[cellList[i]];
            
            // WARNING Time-independence assumed here. _R and _V not copied, as boundary conditions cells should not set/use them
            to->parameters[CellParams::RHO_DT2] += factor*incomingCell->parameters[CellParams::RHO_DT2];
            to->parameters[CellParams::RHOVX_DT2] += factor*incomingCell->parameters[CellParams::RHOVX_DT2];
            to->parameters[CellParams::RHOVY_DT2] += factor*incomingCell->parameters[CellParams::RHOVY_DT2];
            to->parameters[CellParams::RHOVZ_DT2] += factor*incomingCell->parameters[CellParams::RHOVZ_DT2];
            to->parameters[CellParams::P_11_DT2] += factor*incomingCell->parameters[CellParams::P_11_DT2];
            to->parameters[CellParams::P_22_DT2] += factor*incomingCell->parameters[CellParams::P_22_DT2];
            to->parameters[CellParams::P_33_DT2] += factor*incomingCell->parameters[CellParams::P_33_DT2];
            to->parameters[CellParams::RHO] += factor*incomingCell->parameters[CellParams::RHO];
            to->parameters[CellParams::RHOVX] += factor*incomingCell->parameters[CellParams::RHOVX];
            to->parameters[CellParams::RHOVY] += factor*incomingCell->parameters[CellParams::RHOVY];
            to->parameters[CellParams::RHOVZ] += factor*incomingCell->parameters[CellParams::RHOVZ];
            to->parameters[CellParams::P_11] += factor*incomingCell->parameters[CellParams::P_11];
            to->parameters[CellParams::P_22] += factor*incomingCell->parameters[CellParams::P_22];
            to->parameters[CellParams::P_33] += factor*incomingCell->parameters[CellParams::P_33];

            // Do this only for the first layer, the other layers do not need this.
            if (to->sysBoundaryLayer != 1) continue;
            
            const Real* blockParameters = incomingCell->get_block_parameters();
            for (vmesh::LocalID blockLID=0; blockLID<incomingCell->get_number_of_velocity_blocks(); ++blockLID) {
               //const Real* blockParameters = incomingCell->get_block_parameters(blockLID);
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
                  to->increment_value(
                     vxCellCenter,
                     vyCellCenter,
                     vzCellCenter,
                     factor*incomingCell->get_value(vxCellCenter, vyCellCenter, vzCellCenter)
                  );
               }
               blockParameters += BlockParams::N_VELOCITY_BLOCK_PARAMS;
            } // for-loop over velocity blocks
         }
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
      creal& nz
   ) {
      SpatialCell * cell = mpiGrid[cellID];
      const std::vector<CellID> cellList = this->getAllClosestNonsysboundaryCells(mpiGrid, cellID);
      const size_t numberOfCells = cellList.size();
      creal factor = 1.0 / convert<Real>(numberOfCells);
      
      cell->clear();
      
      for (size_t i=0; i<numberOfCells; i++) {
         SpatialCell* incomingCell = mpiGrid[cellList[i]];
         const Real* blockParameters = incomingCell->get_block_parameters();

         // add blocks
         for (vmesh::LocalID blockLID=0; blockLID<incomingCell->get_number_of_velocity_blocks(); ++blockLID) {
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
               if(vNormal >= 0.0) {
                  // Not flowing in, leave as is.
                  cell->increment_value(
                     vxCellCenter,
                     vyCellCenter,
                     vzCellCenter,
                     factor*incomingCell->get_value(vxCellCenter, vyCellCenter, vzCellCenter)
                  );
               } else {
                  // Flowing in, bounce off.
                  cell->increment_value(
                     vxCellCenter - 2.0*vNormal*nx,
                     vyCellCenter - 2.0*vNormal*ny,
                     vzCellCenter - 2.0*vNormal*nz,
                     factor*incomingCell->get_value(vxCellCenter, vyCellCenter, vzCellCenter)
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
    * \param nx Unit vector x component normal to the bounce/reflection plane.
    * \param ny Unit vector y component normal to the bounce/reflection plane.
    * \param nz Unit vector z component normal to the bounce/reflection plane.
    */
   void SysBoundaryCondition::vlasovBoundaryAbsorb(
      const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const CellID& cellID,
      creal& nx,
      creal& ny,
      creal& nz
   ) {
      SpatialCell * cell = mpiGrid[cellID];
      const std::vector<CellID> cellList = this->getAllClosestNonsysboundaryCells(mpiGrid, cellID);
      const size_t numberOfCells = cellList.size();
      creal factor = 1.0 / convert<Real>(numberOfCells);
      
      cell->clear();
      
      for (size_t i=0; i<numberOfCells; i++) {
         SpatialCell * incomingCell = mpiGrid[cellList[i]];
         const Real* blockParameters = incomingCell->get_block_parameters();   
         
         // add blocks
         for (vmesh::LocalID blockLID=0; blockLID<incomingCell->get_number_of_velocity_blocks(); ++blockLID) {
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
               if(vNormal >= 0.0) {
                  // Not flowing in, leave as is.
                  cell->increment_value(
                     vxCellCenter,
                     vyCellCenter,
                     vzCellCenter,
                     factor*incomingCell->get_value(vxCellCenter, vyCellCenter, vzCellCenter)
                  );
               } else {
                  // Flowing in, bounce off.
                  cell->increment_value(
                     vxCellCenter,
                     vyCellCenter,
                     vzCellCenter,
                     0.0
                  );
               }
            }
         } // for-loop over velocity blocks
         blockParameters += BlockParams::N_VELOCITY_BLOCK_PARAMS;
      } // for-loop over spatial cells
   }
   
   /*! Get the cellID of the first closest cell of type NOT_SYSBOUNDARY found.
    * \return The cell index of that cell
    * \sa getAllClosestNonsysboundaryCells
    */
   CellID SysBoundaryCondition::getTheClosestNonsysboundaryCell(
      const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const CellID& cellID
   ) {
      phiprof::start("getTheClosestNonsysboundaryCell");
      CellID closestCell = INVALID_CELLID;
      uint dist = numeric_limits<uint>::max();
      
      for(int i=-2; i<3; i++)
         for(int j=-2; j<3; j++)
            for(int k=-2; k<3; k++) {
               const CellID cell = getNeighbour(mpiGrid,cellID,i,j,k);
               if(cell != INVALID_CELLID) {
                  if(mpiGrid[cell]->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
                     cuint d2 =  i*i+j*j+k*k;
                     if(d2 < dist) {
                        dist = d2;
                        closestCell = cell;
                     }
                  }
               }
            }
      phiprof::stop("getTheClosestNonsysboundaryCell");
      return closestCell;
   }
   
   /*! Get the cellIDs of all the closest cells of type NOT_SYSBOUNDARY.
    * \return The vector of cell indices of those cells
    * \sa getTheClosestNonsysboundaryCell
    */
   std::vector<CellID> SysBoundaryCondition::getAllClosestNonsysboundaryCells(
      const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const CellID& cellID
   ) {
      phiprof::start("getAllClosestNonsysboundaryCells");
      std::vector<CellID> closestCells;
      uint dist = numeric_limits<uint>::max();
      
      // First iteration of search to determine closest distance
      for(int i=-2; i<3; i++)
         for(int j=-2; j<3; j++)
            for(int k=-2; k<3; k++) {
               const CellID cell = getNeighbour(mpiGrid,cellID,i,j,k);
               if(cell != INVALID_CELLID) {
                  if(mpiGrid[cell]->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
                     cuint d2 = i*i+j*j+k*k;
                     if(d2 < dist) {
                        dist = d2;
                     }
                  }
               }
            }
      // Second iteration to record the cellIDs of all cells at closest distance
      for(int i=-2; i<3; i++)
         for(int j=-2; j<3; j++)
            for(int k=-2; k<3; k++) {
               const CellID cell = getNeighbour(mpiGrid,cellID,i,j,k);
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
      phiprof::stop("getAllClosestNonsysboundaryCells");
      return closestCells;
   }
   
   /*! Function used in some cases to know which faces the system boundary condition is being applied to.*/
   void SysBoundaryCondition::getFaces(bool* faces) {
      cerr << "ERROR: SysBoundaryCondition::getFaces called instead of derived class function!" << endl;
   }
   
   /*! Get the name of the system boundary condition.
    * \return The name of the system boundary. The base class function returns an empty string.
    */
   string SysBoundaryCondition::getName() const {
      cerr << "ERROR: SysBoundaryCondition::getName called instead of derived class function!" << endl;
      return string("");
   }
   
   /*! Get the enum index of the system boundary condition.
    * \return The index of the system boundary condition as enumerated in namespace sysboundarytype. The base class function returns 0.
    */
   uint SysBoundaryCondition::getIndex() const {
      cerr << "ERROR: SysBoundaryCondition::getIndex called instead of derived class function!" << endl;
      return 0;
   }
   
   /*! Get the precedence value of the system boundary condition.
    * \return The precedence value of the system boundary condition as set by parameter.
    */
   uint SysBoundaryCondition::getPrecedence() const {return precedence;}
   
   /*! Returns whether the boundary condition is dynamic in time.
    * \return Boolean value.
    */
   bool SysBoundaryCondition::isDynamic() const {return isThisDynamic;}
   
} // namespace SBC
