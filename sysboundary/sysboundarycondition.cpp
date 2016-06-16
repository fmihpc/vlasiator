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
   
   /*! SysBoundaryCondition base class instance of the getParameters function. Should not be used, each derived class should have its own.*/
   void SysBoundaryCondition::getParameters() {
      cerr << "ERROR: SysBoundaryCondition::getParameters called instead of derived class function!" << endl;
   }
   
   /*! Function called during initialisation to set the system boundary condition's parameters.
    * This function must initialize the boundary condition function for all particle species,
    * i.e., its behavior is not allowed to depend on SysBoundaryCondition::activePopID.
    * Base version should not be used, each derived class should have its own.
    */
   bool SysBoundaryCondition::initSysBoundary(
      creal& t,
      Project &project
   ) {
      cerr << "ERROR: SysBoundaryCondition::initSysBoundary called instead of derived class function!" << endl;
      return false;
   }

   /*! Function used to assign the system boundary condition type to a cell.
    * \param mpiGrid Grid
    * \return The system boundary condition type's index
    */
   bool SysBoundaryCondition::assignSysBoundary(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid) {
      cerr << "ERROR: SysBoundaryCondition::assignSysBoundary called instead of derived class function!" << endl;
      return false;
   }
   
   /*! Function used to apply the system boundary condition initial state to a cell.
    * \param mpiGrid Grid
    * \param project Project object
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
    * \param mpiGrid Grid
    * \param cellCache Field solver cell cache
    * \param localID Field solver cache local cell ID
    * \param dt Time step length
    * \param RKCase The step in Runge-Kutta (use values from enum defined in common.h).
    * \param offset Offset of 0-index to the PERBX component, usually 0 or CellParams::PERBX_DT2 - CellParams::PERBX
    * \param component 0: x-component, 1: y-component, 2: z-component.
    * 
    * \return The requested component value.
    */

   Real SysBoundaryCondition::fieldSolverBoundaryCondMagneticField(
      const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const std::vector<fs_cache::CellCache>& cellCache,
      const uint16_t& localID,
      creal& dt,
      cuint& RKCase,
      cint& offset,
      cuint& component
   ) {
      cerr << "ERROR: SysBoundaryCondition::fieldSolverBoundaryCondMagneticField called instead of derived class function!" << endl;
      exit(1);
   }
   
   /*! Function used to return the system boundary condition cell's main electric field component (without Hall term).
    * \param mpiGrid Grid
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
    * \param cache Field solver cell cache
    * \param RKCase The step in Runge-Kutta (use values from enum defined in common.h).
    * \param component 0: x-component, 1: y-component, 2: z-component.
    */
   void SysBoundaryCondition::fieldSolverBoundaryCondHallElectricField(
      fs_cache::CellCache& cache,
      cuint RKCase,
      cuint component
   ) {
      cerr << "ERROR: SysBoundaryCondition::fieldSolverBoundaryCondHallElectricField called instead of derived class function!" << endl;
      exit(1);
   }
   
   /*! Function used to compute the system boundary condition cell's electron pressure gradient electric field components and save them into the cell parameters.
    * \param cache Field solver cell cache.
    * \param RKCase The step in Runge-Kutta (use values from enum defined in common.h).
    * \param component 0: x-component, 1: y-component, 2: z-component.
    */
   void SysBoundaryCondition::fieldSolverBoundaryCondGradPeElectricField(
      fs_cache::CellCache& cache,
      cuint RKCase,
      cuint component
   ) {
      cerr << "ERROR: SysBoundaryCondition::fieldSolverBoundaryCondGradPeElectricField called instead of derived class function!" << endl;
      exit(1);
   }
   
   /*! Function used to compute the system boundary condition cell's derivatives and save them into the cell derivatives.
    * \param mpiGrid Grid
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
    * \param mpiGrid Grid
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
    * \param mpiGrid Grid
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
    * \param mpiGrid Grid
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
   

   /*! Function used to compute the system boundary condition for 
    * the distribution function and moments in the given spatial cell.
    * @param mpiGrid The grid.
    * @param cellID The cell's ID.
    * @param popID Particle species ID.*/
   void SysBoundaryCondition::vlasovBoundaryCondition(
         const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
         const CellID& cellID,
         const int& popID
   ) {
      cerr << "ERROR: SysBoundaryCondition::vlasovBoundaryCondition called instead of derived class function!" << endl;
      exit(1);
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
         const int& popID
   ) {
      const CellID closestCell = getTheClosestNonsysboundaryCell(cellID);
      
      if(closestCell == INVALID_CELLID) {
         cerr << __FILE__ << ":" << __LINE__ << ": No closest cell found!" << endl;
         abort();
      }
      //Do not allow block adjustment, the block structure when calling vlasovBoundaryCondition should be static
      copyCellData(mpiGrid[closestCell],mpiGrid[cellID],false, copyMomentsOnly, popID);
   }
   
   /*! Function used to average and copy the distribution and moments from all the closest sysboundarytype::NOT_SYSBOUNDARY cells.
    * \param mpiGrid Grid
    * \param cellID The cell's ID.
    */
   void SysBoundaryCondition::vlasovBoundaryCopyFromAllClosestNbrs(
      const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const CellID& cellID,const int& popID
   ) {
      const std::vector<CellID> closestCells = getAllClosestNonsysboundaryCells(cellID);
      
      if(closestCells[0] == INVALID_CELLID) {
         cerr << __FILE__ << ":" << __LINE__ << ": No closest cell found!" << endl;
         abort();
      }
      averageCellData(mpiGrid, closestCells, mpiGrid[cellID],popID);
   }
   
   /*! Function used to copy the distribution from (one of) the closest sysboundarytype::NOT_SYSBOUNDARY cell but limiting to values no higher than where it can flow into. Moments are recomputed.
    * \param mpiGrid Grid
    * \param cellID The cell's ID.
    */
   void SysBoundaryCondition::vlasovBoundaryCopyFromTheClosestNbrAndLimit(
      const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const CellID& cellID,
      const int& popID
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
      calculateCellMoments(to,true,true);
   }
   
   /*! Function used to copy the distribution and moments from one cell to another. In layer 2, copy only the moments.
    * \param from Pointer to parent cell to copy from.
    * \param to Pointer to destination cell.
    * \param allowBlockAdjustment If true, blocks can be created or destroyed. If false, only blocks existing in the destination cell are copied.
    */
   void SysBoundaryCondition::copyCellData(
            SpatialCell* from,
            SpatialCell* to,
            bool allowBlockAdjustment,
            const bool& copyMomentsOnly,
            const int& popID
   ) {
      // WARNING Time-independence assumed here. _R and _V not copied, 
      // as boundary conditions cells should not set/use them.
      if (popID == 0) {
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
      }
      if(to->sysBoundaryLayer == 1 && !copyMomentsOnly) { // Do this only for the first layer, the other layers do not need this. Do only if copyMomentsOnly is false.

      // Do this only for the first layer, the other layers do not need this.
      if (to->sysBoundaryLayer != 1) return;

      if (allowBlockAdjustment) {
         // prepare list of blocks to remove. It is not safe to loop 
         // over velocity_block_list while adding/removing blocks
         std::vector<vmesh::GlobalID> blocksToRemove;
         for (vmesh::LocalID block_i=0; block_i<to->get_number_of_velocity_blocks(popID); ++block_i) {
            const vmesh::GlobalID blockGID = to->get_velocity_block_global_id(block_i,popID);

            // If this block does not exist in from, mark it for removal.
            if (from->get_velocity_block_local_id(blockGID,popID) == from->invalid_local_id()) {
               blocksToRemove.push_back(blockGID);
            }
         }

         // remove blocks
         for (size_t b=0; b<blocksToRemove.size(); ++b) {
            cuint blockID=blocksToRemove[b];
            to->remove_velocity_block(blockID,popID);
         }

         // add blocks
         const Realf* fromBlock_data = from->get_data(popID);
         for (vmesh::LocalID block_i=0; block_i<from->get_number_of_velocity_blocks(popID); ++block_i) {
            const vmesh::GlobalID blockGID = from->get_velocity_block_global_id(block_i,popID);
            
            // Ensure that target block exists in 'to' cell. 
            // We must get a new pointer to the 'to' data array here 
            // because add_velocity_block may reallocate it.
            to->add_velocity_block(blockGID,popID);
            Realf* toBlock_data = to->get_data(popID);
            const vmesh::LocalID toBlockLID = to->get_velocity_block_local_id(blockGID,popID);

            for (unsigned int i = 0; i < SIZE_VELBLOCK; i++) {
               toBlock_data[toBlockLID*SIZE_VELBLOCK+i] = fromBlock_data[block_i*SIZE_VELBLOCK+i];
            }
         }
      } else {         
         //just copy data to existing blocks, no modification of to blocks allowed
         const Realf* fromBlock_data = from->get_data(popID);
         Realf* toBlock_data = to->get_data(popID);
         for (vmesh::LocalID block_i=0; block_i<to->get_number_of_velocity_blocks(popID); ++block_i) {
            const vmesh::GlobalID blockGID = to->get_velocity_block_global_id(block_i,popID);
            const vmesh::LocalID fromBlockLID = from->get_velocity_block_local_id(blockGID,popID);
            if (from->get_velocity_block_local_id(blockGID,popID) == from->invalid_local_id()) {
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
   }
   
   /*! Take a list of cells and set the destination cell distribution function to the average of the list's cells'.
    *  For layer 1 the whole distribution function is copied.
    *  For layer >1, only moments are copied
    * \param mpiGrid Grid
    * \param cellList Vector of cells to copy from.
    * \param to Pointer to cell in which to set the averaged distribution.
    */
   void SysBoundaryCondition::averageCellData(
         const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
         const std::vector<CellID> cellList,
         SpatialCell *to,
         const int& popID
   ) {
      const size_t numberOfCells = cellList.size();
      if(numberOfCells == 1) {
         copyCellData(mpiGrid[cellList[0]], to, true, false, popID);
      } else {
         creal factor = 1.0 / convert<Real>(numberOfCells);

         if (popID == 0) {
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
         }
         to->clear(popID);
         
         for (size_t i=0; i<numberOfCells; i++) {
            const SpatialCell* incomingCell = mpiGrid[cellList[i]];
            
            // WARNING Time-independence assumed here. _R and _V not copied, as boundary conditions cells should not set/use them
            if (popID == 0) {
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
            }

            // Do this only for the first layer, the other layers do not need this.
            if (to->sysBoundaryLayer != 1) continue;

            const Real* blockParameters = incomingCell->get_block_parameters(popID);
            const Realf* fromData = incomingCell->get_data(popID);
            for (vmesh::LocalID incBlockLID=0; incBlockLID<incomingCell->get_number_of_velocity_blocks(popID); ++incBlockLID) {
               // Check where cells are
               creal vxBlock = blockParameters[BlockParams::VXCRD];
               creal vyBlock = blockParameters[BlockParams::VYCRD];
               creal vzBlock = blockParameters[BlockParams::VZCRD];
               creal dvxCell = blockParameters[BlockParams::DVX];
               creal dvyCell = blockParameters[BlockParams::DVY];
               creal dvzCell = blockParameters[BlockParams::DVZ];
               
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
         creal& nz,
         const int& popID
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
      const int& popID
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
      }
      return true;
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
      phiprof::start("getAllClosestNonsysboundaryCells");
      std::vector<CellID> & closestCells = allClosestNonsysboundaryCells.at(cellID);
      phiprof::stop("getAllClosestNonsysboundaryCells");
      return closestCells;
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
      const int& popID
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
   
   /*! Function used in some cases to know which faces the system boundary condition is being applied to.
    * \param faces Pointer to array of 6 bool in which the values are returned whether the corresponding face is of that type. Order: 0 x+; 1 x-; 2 y+; 3 y-; 4 z+; 5 z-
    */
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
