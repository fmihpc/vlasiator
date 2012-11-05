/*
This file is part of Vlasiator.

Copyright 2010, 2011, 2012 Finnish Meteorological Institute

Vlasiator is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3
as published by the Free Software Foundation.

Vlasiator is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
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
    */
   void SysBoundaryCondition::determineFace(bool* isThisCellOnAFace,
                                            creal x, creal y, creal z,
                                            creal dx, creal dy, creal dz) {
      for(uint i=0; i<6; i++) {
         isThisCellOnAFace[i] = false;
      }
      if(x > Parameters::xmax - dx) isThisCellOnAFace[0] = true;
      if(x < Parameters::xmin + dx) isThisCellOnAFace[1] = true;
      if(y > Parameters::ymax - dy) isThisCellOnAFace[2] = true;
      if(y < Parameters::ymin + dy) isThisCellOnAFace[3] = true;
      if(z > Parameters::zmax - dz) isThisCellOnAFace[4] = true;
      if(z < Parameters::zmin + dz) isThisCellOnAFace[5] = true;
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
   
   /*! Function called at initialisation to set the system boundary condition's parameters.
    */
   bool SysBoundaryCondition::initSysBoundary(creal& t) {
      cerr << "ERROR: SysBoundaryCondition::initSysBoundary called instead of derived class function!" << endl;
      return false;
   }
   
   /*! Function used to assign the system boundary condition type to a cell.
    * \param cellParams Pointer to the cell's parameters array.
    * \return The system boundary condition type's index
    */
   bool SysBoundaryCondition::assignSysBoundary(dccrg::Dccrg<SpatialCell>& mpiGrid) {
      cerr << "ERROR: SysBoundaryCondition::assignSysBoundary called instead of derived class function!" << endl;
      return false;
   }
   
   /*! Function used to apply the system boundary condition initial state to a cell. */
   bool SysBoundaryCondition::applyInitialState(const dccrg::Dccrg<SpatialCell>& mpiGrid) {
      cerr << "ERROR: SysBoundaryCondition::applyInitialState called instead of derived class function!" << endl;
      return false;
   }
   
//    /*! Function used to apply the system boundary condition state to a cell at given time t. */
//    bool SysBoundaryCondition::applySysBoundaryCondition(
//       const dccrg::Dccrg<SpatialCell>& mpiGrid,
//       creal& t
//    ) {
//       cerr << "ERROR: SysBoundaryCondition::applySysBoundaryCondition called instead of derived class function!" << endl;
//       return false;
//    }
   
   Real SysBoundaryCondition::fieldSolverBoundaryCondMagneticField(
      const dccrg::Dccrg<SpatialCell>& mpiGrid,
      const CellID& cellID,
      creal& dt,
      cuint& component
   ) {
      cerr << "ERROR: SysBoundaryCondition::fieldSolverBoundaryCondMagneticField called instead of derived class function!" << endl;
      exit(1);
   }
   
   void SysBoundaryCondition::fieldSolverBoundaryCondElectricField(
      dccrg::Dccrg<SpatialCell>& mpiGrid,
      const CellID& cellID,
      cuint RKCase,
      cuint component
   ) {
      cerr << "ERROR: SysBoundaryCondition::fieldSolverBoundaryCondElectricField called instead of derived class function!" << endl;
      exit(1);
   }
   
   void SysBoundaryCondition::fieldSolverBoundaryCondDerivatives(
      const dccrg::Dccrg<SpatialCell>& mpiGrid,
      const CellID& cellID,
      cuint& component
   ) {
      cerr << "ERROR: SysBoundaryCondition::fieldSolverBoundaryCondDerivatives called instead of derived class function!" << endl;
      exit(1);
   }
   
   void SysBoundaryCondition::fieldSolverBoundaryCondBVOLDerivatives(
      const dccrg::Dccrg<SpatialCell>& mpiGrid,
      const CellID& cellID,
      cuint& component
   ) {
      cerr << "ERROR: SysBoundaryCondition::fieldSolverBoundaryCondBVOLDerivatives called instead of derived class function!" << endl;
      exit(1);
   }
   
   void SysBoundaryCondition::setCellDerivativesToZero(
      const dccrg::Dccrg<SpatialCell>& mpiGrid,
      const CellID& cellID,
      cuint& component
   ) {
      Real* const derivs = &(mpiGrid[cellID]->derivatives[0]);
      switch(component) {
         case 0:
            derivs[fieldsolver::drhodx] = 0.0;
            derivs[fieldsolver::dBydx]  = 0.0;
            derivs[fieldsolver::dBzdx]  = 0.0;
            derivs[fieldsolver::dVxdx]  = 0.0;
            derivs[fieldsolver::dVydx]  = 0.0;
            derivs[fieldsolver::dVzdx]  = 0.0;
            break;
         case 1:
            derivs[fieldsolver::drhody] = 0.0;
            derivs[fieldsolver::dBxdy]  = 0.0;
            derivs[fieldsolver::dBzdy]  = 0.0;
            derivs[fieldsolver::dVxdy]  = 0.0;
            derivs[fieldsolver::dVydy]  = 0.0;
            derivs[fieldsolver::dVzdy]  = 0.0;
            break;
         case 2:
            derivs[fieldsolver::drhodz] = 0.0;
            derivs[fieldsolver::dBxdz]  = 0.0;
            derivs[fieldsolver::dBydz]  = 0.0;
            derivs[fieldsolver::dVxdz]  = 0.0;
            derivs[fieldsolver::dVydz]  = 0.0;
            derivs[fieldsolver::dVzdz]  = 0.0;
            break;
         default:
            cerr << "Invalid component" << endl;
      }
   }
   
   void SysBoundaryCondition::setCellBVOLDerivativesToZero(
      const dccrg::Dccrg<SpatialCell>& mpiGrid,
      const CellID& cellID,
      cuint& component
   ) {
      Real* const derivs = &(mpiGrid[cellID]->derivativesBVOL[0]);
      switch(component) {
         case 0:
            derivs[bvolderivatives::dBYVOLdx] = 0.0;
            derivs[bvolderivatives::dBZVOLdx] = 0.0;
            break;
         case 1:
            derivs[bvolderivatives::dBXVOLdy] = 0.0;
            derivs[bvolderivatives::dBZVOLdy] = 0.0;
            break;
         case 2:
            derivs[bvolderivatives::dBXVOLdz] = 0.0;
            derivs[bvolderivatives::dBYVOLdz] = 0.0;
            break;
         default:
            cerr << "Invalid component" << endl;
      }
   }
   
   void SysBoundaryCondition::vlasovBoundaryCondition(
      const dccrg::Dccrg<SpatialCell>& mpiGrid,
      const CellID& cellID
   ) {
      cerr << "ERROR: SysBoundaryCondition::vlasovBoundaryCondition called instead of derived class function!" << endl;
      exit(1);
   }
   
   //if the spatialcells are neighbors
   void SysBoundaryCondition::copyCellData(SpatialCell *from, SpatialCell *to)
   {
      for(uint block_i=0;
          block_i<to->number_of_blocks;
          block_i++) {
         cuint blockID=to->velocity_block_list[block_i];
         if(from->get_block_has_content(blockID) == false) {
            to->remove_velocity_block(blockID);
         }
      }
      
      for(uint block_i=0;
          block_i<from->number_of_blocks;
          block_i++) {
         cuint blockID = from->velocity_block_list[block_i];
         
         to->add_velocity_block(blockID);
         
         const Velocity_Block* block = to->at(blockID);
         
         for (uint kc=0; kc<WID; ++kc)
            for (uint jc=0; jc<WID; ++jc)
               for (uint ic=0; ic<WID; ++ic) {
                  creal vx_cell_center = block->parameters[BlockParams::VXCRD] + (ic+0.5)*block->parameters[BlockParams::DVX];
                  creal vy_cell_center = block->parameters[BlockParams::VYCRD] + (jc+0.5)*block->parameters[BlockParams::DVY];
                  creal vz_cell_center = block->parameters[BlockParams::VZCRD] + (kc+0.5)*block->parameters[BlockParams::DVZ];
                  creal value = from->get_value(vx_cell_center, vy_cell_center, vz_cell_center);
                  to->set_value(vx_cell_center, vy_cell_center, vz_cell_center, value);
               }
      }
      
      // WARNING Time-independence assumed here. _R and _V not copied, as boundary conditions cells should not set/use them
      to->parameters[CellParams::RHO_DT2] = from->parameters[CellParams::RHO_DT2];
      to->parameters[CellParams::RHOVX_DT2] = from->parameters[CellParams::RHOVX_DT2];
      to->parameters[CellParams::RHOVY_DT2] = from->parameters[CellParams::RHOVY_DT2];
      to->parameters[CellParams::RHOVZ_DT2] = from->parameters[CellParams::RHOVZ_DT2];
      to->parameters[CellParams::RHO] = from->parameters[CellParams::RHO];
      to->parameters[CellParams::RHOVX] = from->parameters[CellParams::RHOVX];
      to->parameters[CellParams::RHOVY] = from->parameters[CellParams::RHOVY];
      to->parameters[CellParams::RHOVZ] = from->parameters[CellParams::RHOVZ];

      

      //let's get rid of blocks not fulfilling the criteria here to save memory.
      to->adjustSingleCellVelocityBlocks();
   }
   
//    void SysBoundaryCondition::zeroCellData(SpatialCell *to)
//    {
// //       cout << "WARNING zeroing out a cell in zeroCellData !" << endl;
//       for(unsigned int block_i=0; block_i<to->number_of_blocks;block_i++){
//          unsigned int block = to->velocity_block_list[block_i];         
//          Velocity_Block* to_block_ptr = to->at(block);
//          //if block does not exist in from cell, then the empty null_block will be returned        
//          for (uint i=0; i<SIZE_VELBLOCK; i++) to_block_ptr->data[i]=0.0;
//       }
//    }
   
   /*! Function used in some cases to know which faces the system boundary condition is being applied to.*/
   void SysBoundaryCondition::getFaces(bool* faces) {
      cerr << "ERROR: SysBoundaryCondition::getFaces called instead of derived class function!" << endl;
   }
   
   /*! Get the name of the system boundary condition.
    * @return The name of the data. The base class function returns an empty string.
    */
   string SysBoundaryCondition::getName() const {
      cerr << "ERROR: SysBoundaryCondition::getName called instead of derived class function!" << endl;
      return string("");
   }
   
   /*! Get the enum index of the system boundary condition.
    * @return The index of the system boundary condition as enumerated in namespace sysboundarytype. The base class function returns 0.
    */
   uint SysBoundaryCondition::getIndex() const {
      cerr << "ERROR: SysBoundaryCondition::getIndex called instead of derived class function!" << endl;
      return 0;
   }
   
   /*! Get the precedence value of the system boundary condition.
    * @return The precedence value of the system boundary condition as set by parameter.
    */
   uint SysBoundaryCondition::getPrecedence() const {return precedence;}
   
   /*! Returns whether the boundary condition is dynamic in time.
    * @return Boolean value.
    */
   bool SysBoundaryCondition::isDynamic() const {return isThisDynamic;}
   
} // namespace SBC
