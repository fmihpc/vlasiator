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
   
   /*! Function called at initialisation to set the system boundary condition's parameters.
    */
   bool SysBoundaryCondition::initSysBoundary(
      creal& t,
      Project &project
   ) {
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
   bool SysBoundaryCondition::applyInitialState(
      const dccrg::Dccrg<SpatialCell>& mpiGrid,
      Project &project
   ) {
      cerr << "ERROR: SysBoundaryCondition::applyInitialState called instead of derived class function!" << endl;
      return false;
   }
   
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
      dccrg::Dccrg<SpatialCell>& mpiGrid,
      const CellID& cellID,
      cuint& RKCase,
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
         case 0: // x, xx
            derivs[fieldsolver::drhodx] = 0.0;
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
   void SysBoundaryCondition::copyCellData(SpatialCell *from, SpatialCell *to,bool allowBlockAdjustment)
   {
      if(to->sysBoundaryLayer == 1) { // Do this only for the first layer, the inner layer does not need this.

         if(allowBlockAdjustment) {
/*prepare list of blocks to remove. It is not safe to loop over
 * velocity_block_list while adding/removing blocks*/
            std::vector<uint> blocksToRemove;
            for(uint block_i=0;
                block_i<to->number_of_blocks;
                block_i++) {
               cuint blockID=to->velocity_block_list[block_i];
               if(from->is_null_block(from->at(blockID))) {
                  //this block does not exist in from -> mark for removal.
                  blocksToRemove.push_back(blockID);
               }
            }
         
            /*remove blocks*/
            for(uint block_i=0;
                block_i<blocksToRemove.size();
                block_i++) {
               cuint blockID=blocksToRemove[block_i];
               to->remove_velocity_block(blockID);
            }
         
            /*add blocks*/
            for(uint block_i=0;
                block_i<from->number_of_blocks;
                block_i++) {
               cuint blockID = from->velocity_block_list[block_i];

               to->add_velocity_block(blockID);          
               const Velocity_Block* toBlock = to->at(blockID);
               const Velocity_Block* fromBlock = from->at(blockID);
               for (unsigned int i = 0; i < VELOCITY_BLOCK_LENGTH; i++) {
                  toBlock->data[i] = fromBlock->data[i];
               }
            }
         }
         else{
            //just copy data to existing blocks, no modification of to blocks allowed
            for(uint block_i=0; block_i<to->number_of_blocks;block_i++) {
               cuint blockID=to->velocity_block_list[block_i];
               const Velocity_Block* toBlock = to->at(blockID);
               const Velocity_Block* fromBlock = from->at(blockID);
               if(from->is_null_block(from->at(blockID))) {
                  for (unsigned int i = 0; i < VELOCITY_BLOCK_LENGTH; i++) {
                     toBlock->data[i] = 0.0; //block did not exist in from cell, fill with zeros.
                  }
               }
               else {
                  for (unsigned int i = 0; i < VELOCITY_BLOCK_LENGTH; i++) {
                     toBlock->data[i] = fromBlock->data[i];
                  }     
               }
            }
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
   }
   
   CellID SysBoundaryCondition::getClosestNonsysboundaryCell(
      const dccrg::Dccrg<SpatialCell>& mpiGrid,
      const CellID& cellID
   ) {
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
                        dist = min(dist, d2);
                        closestCell = cell;
                     }
                  }
               }
            }
      return closestCell;
   }
   
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
