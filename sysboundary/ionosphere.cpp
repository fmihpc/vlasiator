/*
 * This file is part of Vlasiator.
 * 
 * Copyright 2010, 2011, 2012, 2013 Finnish Meteorological Institute
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 */

/*!\file ionosphere.cpp
 * \brief Implementation of the class SysBoundaryCondition::Ionosphere to handle cells classified as sysboundarytype::IONOSPHERE.
 */

#include <cstdlib>
#include <iostream>

#include "ionosphere.h"
#include "../projects/project.h"
#include "../projects/projects_common.h"
#include "../vlasovmover.h"
#include "../fieldsolver/fs_common.h"
#include "../fieldsolver/fs_limiters.h"
#include "../common.h"

namespace SBC {
   Ionosphere::Ionosphere(): SysBoundaryCondition() { }
   
   Ionosphere::~Ionosphere() { }
   
   void Ionosphere::addParameters() {
      Readparameters::add("ionosphere.centerX", "X coordinate of ionosphere center (m)", 0.0);
      Readparameters::add("ionosphere.centerY", "Y coordinate of ionosphere center (m)", 0.0);
      Readparameters::add("ionosphere.centerZ", "Z coordinate of ionosphere center (m)", 0.0);
      Readparameters::add("ionosphere.radius", "Radius of ionosphere (m).", 1.0e7);
      Readparameters::add("ionosphere.rho", "Number density of the ionosphere (m^-3)", 1.0e6);
      Readparameters::add("ionosphere.VX0", "Bulk velocity of ionospheric distribution function in X direction (m/s)", 0.0);
      Readparameters::add("ionosphere.VY0", "Bulk velocity of ionospheric distribution function in X direction (m/s)", 0.0);
      Readparameters::add("ionosphere.VZ0", "Bulk velocity of ionospheric distribution function in X direction (m/s)", 0.0);
      Readparameters::add("ionosphere.taperRadius", "Width of the zone with a density tapering from the ionospheric value to the background (m)", 0.0);
      Readparameters::add("ionosphere.precedence", "Precedence value of the ionosphere system boundary condition (integer), the higher the stronger.", 2);
   }
   
   void Ionosphere::getParameters() {
      int myRank;
      MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
      if(!Readparameters::get("ionosphere.centerX", this->center[0])) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!Readparameters::get("ionosphere.centerY", this->center[1])) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!Readparameters::get("ionosphere.centerZ", this->center[2])) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!Readparameters::get("ionosphere.radius", this->radius)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!Readparameters::get("ionosphere.rho", this->rho)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!Readparameters::get("ionosphere.VX0", this->VX0)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!Readparameters::get("ionosphere.VY0", this->VY0)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!Readparameters::get("ionosphere.VZ0", this->VZ0)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!Readparameters::get("Magnetosphere.T", this->T)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!Readparameters::get("ionosphere.precedence", this->precedence)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!Readparameters::get("Magnetosphere.nSpaceSamples", this->nSpaceSamples)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
      if(!Readparameters::get("Magnetosphere.nVelocitySamples", this->nVelocitySamples)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      }
   }
   
   bool Ionosphere::initSysBoundary(
      creal& t,
      Project &project
   ) {
      getParameters();
      isThisDynamic = false;
      
      generateTemplateCell(project);
      
      return true;
   }
   
   bool Ionosphere::assignSysBoundary(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid) {
      vector<CellID> cells = mpiGrid.get_cells();
      for(uint i=0; i<cells.size(); i++) {
         if(mpiGrid[cells[i]]->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
            continue;
         }
         creal* const cellParams = &(mpiGrid[cells[i]]->parameters[0]);
         creal dx = cellParams[CellParams::DX];
         creal dy = cellParams[CellParams::DY];
         creal dz = cellParams[CellParams::DZ];
         creal x = cellParams[CellParams::XCRD] + 0.5*dx;
         creal y = cellParams[CellParams::YCRD] + 0.5*dy;
         creal z = cellParams[CellParams::ZCRD] + 0.5*dz;
         creal r = sqrt((x-center[0])*(x-center[0]) + (y-center[1])*(y-center[1]) + (z-center[2])*(z-center[2]));
         
         if(r < radius) {
            mpiGrid[cells[i]]->sysBoundaryFlag = this->getIndex();
         }
      }
      return true;
   }
   
   bool Ionosphere::applyInitialState(
      const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      Project &project
   ) {
      vector<uint64_t> cells = mpiGrid.get_cells();
#pragma omp parallel for
      for (uint i=0; i<cells.size(); ++i) {
         SpatialCell* cell = mpiGrid[cells[i]];
         if(cell->sysBoundaryFlag != this->getIndex()) continue;
         setCellFromTemplate(cell);
      }
      return true;
   }
   
//    bool Ionosphere::applySysBoundaryCondition(
//       const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
//       creal& t
//    ) {
//       return true;
//    }
   
   Real Ionosphere::fieldSolverBoundaryCondMagneticField(
      const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const CellID& cellID,
      creal& dt,
      cuint& component
   ) {
      // The perturbed magnetic field is reset to 0.0, the dipole field is in the background component.
      return 0.0;
   }
   
   void Ionosphere::fieldSolverBoundaryCondElectricField(
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
      return;
   }
   
   void Ionosphere::fieldSolverBoundaryCondHallElectricField(
      dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid,
      const CellID& cellID,
      cuint RKCase,
      cuint component
   ) {
      switch(component) {
         case 0:
            mpiGrid[cellID]->parameters[CellParams::EXHALL_000_100] = 0.0;
            mpiGrid[cellID]->parameters[CellParams::EXHALL_010_110] = 0.0;
            mpiGrid[cellID]->parameters[CellParams::EXHALL_001_101] = 0.0;
            mpiGrid[cellID]->parameters[CellParams::EXHALL_011_111] = 0.0;
            break;
         case 1:
            mpiGrid[cellID]->parameters[CellParams::EYHALL_000_010] = 0.0;
            mpiGrid[cellID]->parameters[CellParams::EYHALL_100_110] = 0.0;
            mpiGrid[cellID]->parameters[CellParams::EYHALL_001_011] = 0.0;
            mpiGrid[cellID]->parameters[CellParams::EYHALL_101_111] = 0.0;
            break;
         case 2:
            mpiGrid[cellID]->parameters[CellParams::EZHALL_000_001] = 0.0;
            mpiGrid[cellID]->parameters[CellParams::EZHALL_100_101] = 0.0;
            mpiGrid[cellID]->parameters[CellParams::EZHALL_010_011] = 0.0;
            mpiGrid[cellID]->parameters[CellParams::EZHALL_110_111] = 0.0;
            break;
         default:
            cerr << __FILE__ << ":" << __LINE__ << ":" << " Invalid component" << endl;
      }
   }
   
   void Ionosphere::fieldSolverBoundaryCondDerivatives(
      dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const CellID& cellID,
      cuint& RKCase,
      cuint& component
   ) {
      // For B: use background B + perturbed B in normal cells, only background B in the ionosphere and DO_NOT_COMPUTE cells
      // For RHO and V: use self. One could also use the DO_NOT_COMPUTE cells and give them the ionospheric values too but that means more code changes than just here.
      Real* const array       = mpiGrid[cellID]->derivatives;
      CellID leftNbrID,rghtNbrID;
      creal* rhovLeft = NULL;
      creal* left = NULL;
      creal* cent = mpiGrid[cellID]->parameters;
      creal* rhovRght = NULL;
      creal* rght = NULL;
      CellID botLeftNbrID, botRghtNbrID, topLeftNbrID, topRghtNbrID;
      creal* botLeft = NULL;
      creal* botRght = NULL;
      creal* topLeft = NULL;
      creal* topRght = NULL;
      switch(component) {
         namespace cp = CellParams;
         namespace fs = fieldsolver;
         case 0: // x,xx
            leftNbrID = getNeighbourID(mpiGrid,cellID,2-1,2  ,2  );
            rghtNbrID = getNeighbourID(mpiGrid,cellID,2+1,2  ,2  );
            left = mpiGrid[leftNbrID]->parameters;
            rght = mpiGrid[rghtNbrID]->parameters;
            if(mpiGrid[leftNbrID]->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
               rhovLeft = mpiGrid[cellID]->parameters;
            } else{
               rhovLeft = mpiGrid[leftNbrID]->parameters;
            }
            if(mpiGrid[rghtNbrID]->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
               rhovRght = mpiGrid[cellID]->parameters;
            } else{
               rhovRght = mpiGrid[rghtNbrID]->parameters;
            }
            if(RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
               array[fs::drhodx] = limiter(rhovLeft[cp::RHO],cent[cp::RHO],rhovRght[cp::RHO]);
               array[fs::dp11dx] = limiter(rhovLeft[cp::P_11],cent[cp::P_11],rhovRght[cp::P_11]);
               array[fs::dp22dx] = limiter(rhovLeft[cp::P_22],cent[cp::P_22],rhovRght[cp::P_22]);
               array[fs::dp33dx] = limiter(rhovLeft[cp::P_33],cent[cp::P_33],rhovRght[cp::P_33]);
               array[fs::dVxdx]  = limiter(divideIfNonZero(rhovLeft[cp::RHOVX], rhovLeft[cp::RHO]),
                                           divideIfNonZero(    cent[cp::RHOVX],     cent[cp::RHO]),
                                           divideIfNonZero(rhovRght[cp::RHOVX], rhovRght[cp::RHO]));
               array[fs::dVydx]  = limiter(divideIfNonZero(rhovLeft[cp::RHOVY], rhovLeft[cp::RHO]),
                                           divideIfNonZero(    cent[cp::RHOVY],     cent[cp::RHO]),
                                           divideIfNonZero(rhovRght[cp::RHOVY], rhovRght[cp::RHO]));
               array[fs::dVzdx]  = limiter(divideIfNonZero(rhovLeft[cp::RHOVZ], rhovLeft[cp::RHO]),
                                           divideIfNonZero(    cent[cp::RHOVZ],     cent[cp::RHO]),
                                           divideIfNonZero(rhovRght[cp::RHOVZ], rhovRght[cp::RHO]));
               array[fs::dPERBydx]  = limiter(left[cp::PERBY],cent[cp::PERBY],rght[cp::PERBY]);
               array[fs::dBGBydx]  = limiter(left[cp::BGBY],cent[cp::BGBY],rght[cp::BGBY]);
               array[fs::dPERBzdx]  = limiter(left[cp::PERBZ],cent[cp::PERBZ],rght[cp::PERBZ]);
               array[fs::dBGBzdx]  = limiter(left[cp::BGBZ],cent[cp::BGBZ],rght[cp::BGBZ]);
               if(Parameters::ohmHallTerm < 2) {
                  array[fs::dPERBydxx] = 0.0;
                  array[fs::dPERBzdxx] = 0.0;
               } else {
                  array[fs::dPERBydxx] = left[cp::PERBY] + rght[cp::PERBY] - 2.0*cent[cp::PERBY];
                  array[fs::dPERBzdxx] = left[cp::PERBZ] + rght[cp::PERBZ] - 2.0*cent[cp::PERBZ];
               }
            }
            if (RKCase == RK_ORDER2_STEP1) {
               array[fs::drhodx] = limiter(rhovLeft[cp::RHO_DT2],cent[cp::RHO_DT2],rhovRght[cp::RHO_DT2]);
               array[fs::dp11dx] = limiter(rhovLeft[cp::P_11_DT2],cent[cp::P_11_DT2],rhovRght[cp::P_11_DT2]);
               array[fs::dp22dx] = limiter(rhovLeft[cp::P_22_DT2],cent[cp::P_22_DT2],rhovRght[cp::P_22_DT2]);
               array[fs::dp33dx] = limiter(rhovLeft[cp::P_33_DT2],cent[cp::P_33_DT2],rhovRght[cp::P_33_DT2]);
               array[fs::dVxdx]  = limiter(divideIfNonZero(rhovLeft[cp::RHOVX_DT2], rhovLeft[cp::RHO_DT2]),
                                           divideIfNonZero(    cent[cp::RHOVX_DT2],     cent[cp::RHO_DT2]),
                                           divideIfNonZero(rhovRght[cp::RHOVX_DT2], rhovRght[cp::RHO_DT2]));
               array[fs::dVydx]  = limiter(divideIfNonZero(rhovLeft[cp::RHOVY_DT2], rhovLeft[cp::RHO_DT2]),
                                           divideIfNonZero(    cent[cp::RHOVY_DT2],     cent[cp::RHO_DT2]),
                                           divideIfNonZero(rhovRght[cp::RHOVY_DT2], rhovRght[cp::RHO_DT2]));
               array[fs::dVzdx]  = limiter(divideIfNonZero(rhovLeft[cp::RHOVZ_DT2], rhovLeft[cp::RHO_DT2]),
                                           divideIfNonZero(    cent[cp::RHOVZ_DT2],     cent[cp::RHO_DT2]),
                                           divideIfNonZero(rhovRght[cp::RHOVZ_DT2], rhovRght[cp::RHO_DT2]));
               array[fs::dPERBydx]  = limiter(left[cp::PERBY_DT2],cent[cp::PERBY_DT2],rght[cp::PERBY_DT2]);
               array[fs::dBGBydx]  = limiter(left[cp::BGBY],cent[cp::BGBY],rght[cp::BGBY]);
               array[fs::dPERBzdx]  = limiter(left[cp::PERBZ_DT2],cent[cp::PERBZ_DT2],rght[cp::PERBZ_DT2]);
               array[fs::dBGBzdx]  = limiter(left[cp::BGBZ],cent[cp::BGBZ],rght[cp::BGBZ]);
               if(Parameters::ohmHallTerm < 2) {
                  array[fs::dPERBydxx] = 0.0;
                  array[fs::dPERBzdxx] = 0.0;
               } else {
                  array[fs::dPERBydxx] = left[cp::PERBY_DT2] + rght[cp::PERBY_DT2] - 2.0*cent[cp::PERBY_DT2];
                  array[fs::dPERBzdxx] = left[cp::PERBZ_DT2] + rght[cp::PERBZ_DT2] - 2.0*cent[cp::PERBZ_DT2];
               }
            }
            break;
         case 1: // y,yy
            leftNbrID = getNeighbourID(mpiGrid,cellID,2  ,2-1,2  );
            rghtNbrID = getNeighbourID(mpiGrid,cellID,2  ,2+1,2  );
            left = mpiGrid[leftNbrID]->parameters;
            rght = mpiGrid[rghtNbrID]->parameters;
            if(mpiGrid[leftNbrID]->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
               rhovLeft = mpiGrid[cellID]->parameters;
            } else{
               rhovLeft = mpiGrid[leftNbrID]->parameters;
            }
            if(mpiGrid[rghtNbrID]->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
               rhovRght = mpiGrid[cellID]->parameters;
            } else{
               rhovRght = mpiGrid[rghtNbrID]->parameters;
            }
            if(RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
               array[fs::drhody] = limiter(rhovLeft[cp::RHO],cent[cp::RHO],rhovRght[cp::RHO]);
               array[fs::dp11dy] = limiter(rhovLeft[cp::P_11],cent[cp::P_11],rhovRght[cp::P_11]);
               array[fs::dp22dy] = limiter(rhovLeft[cp::P_22],cent[cp::P_22],rhovRght[cp::P_22]);
               array[fs::dp33dy] = limiter(rhovLeft[cp::P_33],cent[cp::P_33],rhovRght[cp::P_33]);
               array[fs::dVxdy]  = limiter(divideIfNonZero(rhovLeft[cp::RHOVX], rhovLeft[cp::RHO]),
                                           divideIfNonZero(    cent[cp::RHOVX],     cent[cp::RHO]),
                                           divideIfNonZero(rhovRght[cp::RHOVX], rhovRght[cp::RHO]));
               array[fs::dVydy]  = limiter(divideIfNonZero(rhovLeft[cp::RHOVY], rhovLeft[cp::RHO]),
                                           divideIfNonZero(    cent[cp::RHOVY],     cent[cp::RHO]),
                                           divideIfNonZero(rhovRght[cp::RHOVY], rhovRght[cp::RHO]));
               array[fs::dVzdy]  = limiter(divideIfNonZero(rhovLeft[cp::RHOVZ], rhovLeft[cp::RHO]),
                                           divideIfNonZero(    cent[cp::RHOVZ],     cent[cp::RHO]),
                                           divideIfNonZero(rhovRght[cp::RHOVZ], rhovRght[cp::RHO]));
               array[fs::dPERBxdy]  = limiter(left[cp::PERBX],cent[cp::PERBX],rght[cp::PERBX]);
               array[fs::dBGBxdy]  = limiter(left[cp::BGBX],cent[cp::BGBX],rght[cp::BGBX]);
               array[fs::dPERBzdy]  = limiter(left[cp::PERBZ],cent[cp::PERBZ],rght[cp::PERBZ]);
               array[fs::dBGBzdy]  = limiter(left[cp::BGBZ],cent[cp::BGBZ],rght[cp::BGBZ]);
               if(Parameters::ohmHallTerm < 2) {
                  array[fs::dPERBxdyy] = 0.0;
                  array[fs::dPERBzdyy] = 0.0;
               } else {
                  array[fs::dPERBxdyy] = left[cp::PERBX] + rght[cp::PERBX] - 2.0*cent[cp::PERBX];
                  array[fs::dPERBzdyy] = left[cp::PERBZ] + rght[cp::PERBZ] - 2.0*cent[cp::PERBZ];
               }
            }
            if (RKCase == RK_ORDER2_STEP1) {
               array[fs::drhody] = limiter(rhovLeft[cp::RHO_DT2],cent[cp::RHO_DT2],rhovRght[cp::RHO_DT2]);
               array[fs::dp11dy] = limiter(rhovLeft[cp::P_11_DT2],cent[cp::P_11_DT2],rhovRght[cp::P_11_DT2]);
               array[fs::dp22dy] = limiter(rhovLeft[cp::P_22_DT2],cent[cp::P_22_DT2],rhovRght[cp::P_22_DT2]);
               array[fs::dp33dy] = limiter(rhovLeft[cp::P_33_DT2],cent[cp::P_33_DT2],rhovRght[cp::P_33_DT2]);
               array[fs::dVxdy]  = limiter(divideIfNonZero(rhovLeft[cp::RHOVX_DT2], rhovLeft[cp::RHO_DT2]),
                                           divideIfNonZero(    cent[cp::RHOVX_DT2],     cent[cp::RHO_DT2]),
                                           divideIfNonZero(rhovRght[cp::RHOVX_DT2], rhovRght[cp::RHO_DT2]));
               array[fs::dVydy]  = limiter(divideIfNonZero(rhovLeft[cp::RHOVY_DT2], rhovLeft[cp::RHO_DT2]),
                                           divideIfNonZero(    cent[cp::RHOVY_DT2],     cent[cp::RHO_DT2]),
                                           divideIfNonZero(rhovRght[cp::RHOVY_DT2], rhovRght[cp::RHO_DT2]));
               array[fs::dVzdy]  = limiter(divideIfNonZero(rhovLeft[cp::RHOVZ_DT2], rhovLeft[cp::RHO_DT2]),
                                           divideIfNonZero(    cent[cp::RHOVZ_DT2],     cent[cp::RHO_DT2]),
                                           divideIfNonZero(rhovRght[cp::RHOVZ_DT2], rhovRght[cp::RHO_DT2]));
               array[fs::dPERBxdy]  = limiter(left[cp::PERBX_DT2],cent[cp::PERBX_DT2],rght[cp::PERBX_DT2]);
               array[fs::dBGBxdy]  = limiter(left[cp::BGBX],cent[cp::BGBX],rght[cp::BGBX]);
               array[fs::dPERBzdy]  = limiter(left[cp::PERBZ_DT2],cent[cp::PERBZ_DT2],rght[cp::PERBZ_DT2]);
               array[fs::dBGBzdy]  = limiter(left[cp::BGBZ],cent[cp::BGBZ],rght[cp::BGBZ]);
               if(Parameters::ohmHallTerm < 2) {
                  array[fs::dPERBxdyy] = 0.0;
                  array[fs::dPERBzdyy] = 0.0;
               } else {
                  array[fs::dPERBxdyy] = left[cp::PERBX_DT2] + rght[cp::PERBX_DT2] - 2.0*cent[cp::PERBX_DT2];
                  array[fs::dPERBzdyy] = left[cp::PERBZ_DT2] + rght[cp::PERBZ_DT2] - 2.0*cent[cp::PERBZ_DT2];
               }
            }
            break;
         case 2: // z, zz
            leftNbrID = getNeighbourID(mpiGrid,cellID,2  ,2  ,2-1);
            rghtNbrID = getNeighbourID(mpiGrid,cellID,2  ,2  ,2+1);
            left = mpiGrid[leftNbrID]->parameters;
            rght = mpiGrid[rghtNbrID]->parameters;
            if(mpiGrid[leftNbrID]->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
               rhovLeft = mpiGrid[cellID]->parameters;
            } else{
               rhovLeft = mpiGrid[leftNbrID]->parameters;
            }
            if(mpiGrid[rghtNbrID]->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
               rhovRght = mpiGrid[cellID]->parameters;
            } else{
               rhovRght = mpiGrid[rghtNbrID]->parameters;
            }
            if(RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
               array[fs::drhodz] = limiter(rhovLeft[cp::RHO],cent[cp::RHO],rhovRght[cp::RHO]);
               array[fs::dp11dz] = limiter(rhovLeft[cp::P_11],cent[cp::P_11],rhovRght[cp::P_11]);
               array[fs::dp22dz] = limiter(rhovLeft[cp::P_22],cent[cp::P_22],rhovRght[cp::P_22]);
               array[fs::dp33dz] = limiter(rhovLeft[cp::P_33],cent[cp::P_33],rhovRght[cp::P_33]);
               array[fs::dVxdz]  = limiter(divideIfNonZero(rhovLeft[cp::RHOVX], rhovLeft[cp::RHO]),
                                           divideIfNonZero(    cent[cp::RHOVX],     cent[cp::RHO]),
                                           divideIfNonZero(rhovRght[cp::RHOVX], rhovRght[cp::RHO]));
               array[fs::dVydz]  = limiter(divideIfNonZero(rhovLeft[cp::RHOVY], rhovLeft[cp::RHO]),
                                           divideIfNonZero(    cent[cp::RHOVY],     cent[cp::RHO]),
                                           divideIfNonZero(rhovRght[cp::RHOVY], rhovRght[cp::RHO]));
               array[fs::dVzdz]  = limiter(divideIfNonZero(rhovLeft[cp::RHOVZ], rhovLeft[cp::RHO]),
                                           divideIfNonZero(    cent[cp::RHOVZ],     cent[cp::RHO]),
                                           divideIfNonZero(rhovRght[cp::RHOVZ], rhovRght[cp::RHO]));
               array[fs::dPERBxdz]  = limiter(left[cp::PERBX],cent[cp::PERBX],rght[cp::PERBX]);
               array[fs::dBGBxdz]  = limiter(left[cp::BGBX],cent[cp::BGBX],rght[cp::BGBX]);
               array[fs::dPERBydz]  = limiter(left[cp::PERBY],cent[cp::PERBY],rght[cp::PERBY]);
               array[fs::dBGBydz]  = limiter(left[cp::BGBY],cent[cp::BGBY],rght[cp::BGBY]);
               if(Parameters::ohmHallTerm < 2) {
                  array[fs::dPERBxdzz] = 0.0;
                  array[fs::dPERBydzz] = 0.0;
               } else {
                  array[fs::dPERBxdzz] = left[cp::PERBX] + rght[cp::PERBX] - 2.0*cent[cp::PERBX];
                  array[fs::dPERBydzz] = left[cp::PERBY] + rght[cp::PERBY] - 2.0*cent[cp::PERBY];
               }
            }
            if (RKCase == RK_ORDER2_STEP1) {
               array[fs::drhodz] = limiter(rhovLeft[cp::RHO_DT2],cent[cp::RHO_DT2],rhovRght[cp::RHO_DT2]);
               array[fs::dp11dz] = limiter(rhovLeft[cp::P_11_DT2],cent[cp::P_11_DT2],rhovRght[cp::P_11_DT2]);
               array[fs::dp22dz] = limiter(rhovLeft[cp::P_22_DT2],cent[cp::P_22_DT2],rhovRght[cp::P_22_DT2]);
               array[fs::dp33dz] = limiter(rhovLeft[cp::P_33_DT2],cent[cp::P_33_DT2],rhovRght[cp::P_33_DT2]);
               array[fs::dVxdz]  = limiter(divideIfNonZero(rhovLeft[cp::RHOVX_DT2], rhovLeft[cp::RHO_DT2]),
                                           divideIfNonZero(    cent[cp::RHOVX_DT2],     cent[cp::RHO_DT2]),
                                           divideIfNonZero(rhovRght[cp::RHOVX_DT2], rhovRght[cp::RHO_DT2]));
               array[fs::dVydz]  = limiter(divideIfNonZero(rhovLeft[cp::RHOVY_DT2], rhovLeft[cp::RHO_DT2]),
                                           divideIfNonZero(    cent[cp::RHOVY_DT2],     cent[cp::RHO_DT2]),
                                           divideIfNonZero(rhovRght[cp::RHOVY_DT2], rhovRght[cp::RHO_DT2]));
               array[fs::dVzdz]  = limiter(divideIfNonZero(rhovLeft[cp::RHOVZ_DT2], rhovLeft[cp::RHO_DT2]),
                                           divideIfNonZero(    cent[cp::RHOVZ_DT2],     cent[cp::RHO_DT2]),
                                           divideIfNonZero(rhovRght[cp::RHOVZ_DT2], rhovRght[cp::RHO_DT2]));
               array[fs::dPERBxdz]  = limiter(left[cp::PERBX_DT2],cent[cp::PERBX_DT2],rght[cp::PERBX_DT2]);
               array[fs::dBGBxdz]  = limiter(left[cp::BGBX],cent[cp::BGBX],rght[cp::BGBX]);
               array[fs::dPERBydz]  = limiter(left[cp::PERBY_DT2],cent[cp::PERBY_DT2],rght[cp::PERBY_DT2]);
               array[fs::dBGBydz]  = limiter(left[cp::BGBY],cent[cp::BGBY],rght[cp::BGBY]);
               if(Parameters::ohmHallTerm < 2) {
                  array[fs::dPERBxdzz] = 0.0;
                  array[fs::dPERBydzz] = 0.0;
               } else {
                  array[fs::dPERBxdzz] = left[cp::PERBX_DT2] + rght[cp::PERBX_DT2] - 2.0*cent[cp::PERBX_DT2];
                  array[fs::dPERBydzz] = left[cp::PERBY_DT2] + rght[cp::PERBY_DT2] - 2.0*cent[cp::PERBY_DT2];
               }
            }
            break;
         case 3: // xy
            if(Parameters::ohmHallTerm < 2) {
               array[fs::dPERBzdxy] = 0.0;
            } else {
               botLeftNbrID = getNeighbourID(mpiGrid,cellID,2-1,2-1,2  );
               botRghtNbrID = getNeighbourID(mpiGrid,cellID,2+1,2-1,2  );
               topLeftNbrID = getNeighbourID(mpiGrid,cellID,2-1,2+1,2  );
               topRghtNbrID = getNeighbourID(mpiGrid,cellID,2+1,2+1,2  );
               botLeft = mpiGrid[botLeftNbrID]->parameters;
               botRght = mpiGrid[botRghtNbrID]->parameters;
               topLeft = mpiGrid[topLeftNbrID]->parameters;
               topRght = mpiGrid[topRghtNbrID]->parameters;
               
               if(RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
                  array[fs::dPERBzdxy] = FOURTH * (botLeft[cp::PERBZ] + topRght[cp::PERBZ] - botRght[cp::PERBZ] - topLeft[cp::PERBZ]);
               }
               if (RKCase == RK_ORDER2_STEP1) {
                  array[fs::dPERBzdxy] = FOURTH * (botLeft[cp::PERBZ_DT2] + topRght[cp::PERBZ_DT2] - botRght[cp::PERBZ_DT2] - topLeft[cp::PERBZ_DT2]);
               }
            }
            break;
         case 4: // xz
            if(Parameters::ohmHallTerm < 2) {
               array[fs::dPERBydxz] = 0.0;
            } else {
               botLeftNbrID = getNeighbourID(mpiGrid,cellID,2-1,2  ,2-1);
               botRghtNbrID = getNeighbourID(mpiGrid,cellID,2+1,2  ,2-1);
               topLeftNbrID = getNeighbourID(mpiGrid,cellID,2-1,2  ,2+1);
               topRghtNbrID = getNeighbourID(mpiGrid,cellID,2+1,2  ,2+1);
               botLeft = mpiGrid[botLeftNbrID]->parameters;
               botRght = mpiGrid[botRghtNbrID]->parameters;
               topLeft = mpiGrid[topLeftNbrID]->parameters;
               topRght = mpiGrid[topRghtNbrID]->parameters;
               
               if(RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
                  array[fs::dPERBydxz] = FOURTH * (botLeft[cp::PERBY] + topRght[cp::PERBY] - botRght[cp::PERBY] - topLeft[cp::PERBY]);
               }
               if (RKCase == RK_ORDER2_STEP1) {
                  array[fs::dPERBydxz] = FOURTH * (botLeft[cp::PERBY_DT2] + topRght[cp::PERBY_DT2] - botRght[cp::PERBY_DT2] - topLeft[cp::PERBY_DT2]);
               }
            }
            break;
         case 5: // yz
            if(Parameters::ohmHallTerm < 2) {
               array[fs::dPERBxdyz] = 0.0;
            } else {
               botLeftNbrID = getNeighbourID(mpiGrid,cellID,2  ,2-1,2-1);
               botRghtNbrID = getNeighbourID(mpiGrid,cellID,2  ,2+1,2-1);
               topLeftNbrID = getNeighbourID(mpiGrid,cellID,2  ,2-1,2+1);
               topRghtNbrID = getNeighbourID(mpiGrid,cellID,2  ,2+1,2+1);
               botLeft = mpiGrid[botLeftNbrID]->parameters;
               botRght = mpiGrid[botRghtNbrID]->parameters;
               topLeft = mpiGrid[topLeftNbrID]->parameters;
               topRght = mpiGrid[topRghtNbrID]->parameters;
               
               if(RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
                  array[fs::dPERBxdyz] = FOURTH * (botLeft[cp::PERBX] + topRght[cp::PERBX] - botRght[cp::PERBX] - topLeft[cp::PERBX]);
               }
               if (RKCase == RK_ORDER2_STEP1) {
                  array[fs::dPERBxdyz] = FOURTH * (botLeft[cp::PERBX_DT2] + topRght[cp::PERBX_DT2] - botRght[cp::PERBX_DT2] - topLeft[cp::PERBX_DT2]);
               }
            }
            break;
         default:
            cerr << __FILE__ << ":" << __LINE__ << ":" << " Invalid component" << endl;
      }
   }
   
   void Ionosphere::fieldSolverBoundaryCondBVOLDerivatives(
      const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const CellID& cellID,
      cuint& component
   ) {
      // FIXME This should be OK as the BVOL derivatives are only used for Lorentz force JXB, which is not applied on the ionosphere cells.
      this->setCellBVOLDerivativesToZero(mpiGrid, cellID, component);
   }
   
   void Ionosphere::vlasovBoundaryCondition(
      const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const CellID& cellID
   ) {
      //phiprof::start("vlasovBoundaryCondition (Ionosphere)");
      //No need to copy distribution function, will be constant throughout simulation
      //copyCellData(&templateCell, mpiGrid[cellID],false);
      //phiprof::stop("vlasovBoundaryCondition (Ionosphere)");
   }
   
   void Ionosphere::generateTemplateCell(Project &project) {
      // WARNING not 0.0 here or the dipole() function fails miserably.
      templateCell.sysBoundaryFlag = this->getIndex();
      templateCell.sysBoundaryLayer = 1;
      templateCell.parameters[CellParams::XCRD] = 1.0;
      templateCell.parameters[CellParams::YCRD] = 1.0;
      templateCell.parameters[CellParams::ZCRD] = 1.0;
      templateCell.parameters[CellParams::DX] = 1;
      templateCell.parameters[CellParams::DY] = 1;
      templateCell.parameters[CellParams::DZ] = 1;
      
      vector<uint> blocksToInitialize = this->findBlocksToInitialize(templateCell);
      
      for (uint i = 0; i < blocksToInitialize.size(); i++) {
         const vmesh::GlobalID blockGID = blocksToInitialize.at(i);
         const vmesh::LocalID blockLID = templateCell.get_velocity_block_local_id(blockGID);
         const Real* block_parameters = templateCell.get_block_parameters(blockLID);
         creal vxBlock = block_parameters[BlockParams::VXCRD];
         creal vyBlock = block_parameters[BlockParams::VYCRD];
         creal vzBlock = block_parameters[BlockParams::VZCRD];
         creal dvxCell = block_parameters[BlockParams::DVX];
         creal dvyCell = block_parameters[BlockParams::DVY];
         creal dvzCell = block_parameters[BlockParams::DVZ];

         creal x = templateCell.parameters[CellParams::XCRD];
         creal y = templateCell.parameters[CellParams::YCRD];
         creal z = templateCell.parameters[CellParams::ZCRD];
         creal dx = templateCell.parameters[CellParams::DX];
         creal dy = templateCell.parameters[CellParams::DY];
         creal dz = templateCell.parameters[CellParams::DZ];
         
         // Calculate volume average of distrib. function for each cell in the block.
         for (uint kc=0; kc<WID; ++kc) 
            for (uint jc=0; jc<WID; ++jc) 
               for (uint ic=0; ic<WID; ++ic) {
                  creal vxCell = vxBlock + ic*dvxCell;
                  creal vyCell = vyBlock + jc*dvyCell;
                  creal vzCell = vzBlock + kc*dvzCell;
                  Real average = 0.0;
                  if(this->nVelocitySamples > 1) {
                     creal d_vx = dvxCell / (nVelocitySamples-1);
                     creal d_vy = dvyCell / (nVelocitySamples-1);
                     creal d_vz = dvzCell / (nVelocitySamples-1);
                     for (uint vi=0; vi<nVelocitySamples; ++vi)
                        for (uint vj=0; vj<nVelocitySamples; ++vj)
                           for (uint vk=0; vk<nVelocitySamples; ++vk) {
                              average += shiftedMaxwellianDistribution(
                                 vxCell + vi*d_vx,
                                 vyCell + vj*d_vy,
                                 vzCell + vk*d_vz
                              );
                           }
                           average /= this->nVelocitySamples * this->nVelocitySamples * this->nVelocitySamples;
                  } else {
                     average = shiftedMaxwellianDistribution(
                        vxCell + 0.5*dvxCell,
                        vyCell + 0.5*dvyCell,
                        vzCell + 0.5*dvzCell
                     );
                  }
                  
                  if(average!=0.0){
                     creal vxCellCenter = vxBlock + (ic+convert<Real>(0.5))*dvxCell;
                     creal vyCellCenter = vyBlock + (jc+convert<Real>(0.5))*dvyCell;
                     creal vzCellCenter = vzBlock + (kc+convert<Real>(0.5))*dvzCell;
                     templateCell.set_value(vxCellCenter,vyCellCenter,vzCellCenter,average);
                  }
         }
      }
      //let's get rid of blocks not fulfilling the criteria here to save
      //memory.
      templateCell.adjustSingleCellVelocityBlocks();
      
      calculateCellVelocityMoments(&templateCell, true);
      
      // WARNING Time-independence assumed here. Normal moments computed in setProjectCell
      templateCell.parameters[CellParams::RHO_DT2] = templateCell.parameters[CellParams::RHO];
      templateCell.parameters[CellParams::RHOVX_DT2] = templateCell.parameters[CellParams::RHOVX];
      templateCell.parameters[CellParams::RHOVY_DT2] = templateCell.parameters[CellParams::RHOVY];
      templateCell.parameters[CellParams::RHOVZ_DT2] = templateCell.parameters[CellParams::RHOVZ];
   }
   
   Real Ionosphere::shiftedMaxwellianDistribution(
      creal& vx, creal& vy, creal& vz
   ) {
      return this->rho * pow(physicalconstants::MASS_PROTON /
      (2.0 * M_PI * physicalconstants::K_B * this->T), 1.5) *
      exp(-physicalconstants::MASS_PROTON * ((vx-this->VX0)*(vx-this->VX0) + (vy-this->VY0)*(vy-this->VY0) + (vz-this->VZ0)*(vz-this->VZ0)) /
      (2.0 * physicalconstants::K_B * this->T));
   }
   
   vector<uint> Ionosphere::findBlocksToInitialize(SpatialCell& cell) {
      vector<uint> blocksToInitialize;
      bool search = true;
      uint counter = 0;
      
      while (search) {
         if (0.1 * P::sparseMinValue >
            //shiftedMaxwellianDistribution(counter*blockSize[0], 0.0, 0.0) || counter > P::vxblocks_ini) {
            shiftedMaxwellianDistribution(counter*SpatialCell::get_velocity_base_grid_block_size()[0], 0.0, 0.0) || counter > P::vxblocks_ini) {
            search = false;
         }
         ++counter;
      }
      counter+=2;
      //Real vRadiusSquared = (Real)counter*(Real)counter*SpatialCell::block_dvx*SpatialCell::block_dvx;
      Real vRadiusSquared = (Real)counter*(Real)counter*SpatialCell::get_velocity_base_grid_block_size()[0]*SpatialCell::get_velocity_base_grid_block_size()[0];

      #warning FIXME this should loop over blocks in base grid      
      for (uint kv=0; kv<P::vzblocks_ini; ++kv) 
         for (uint jv=0; jv<P::vyblocks_ini; ++jv)
            for (uint iv=0; iv<P::vxblocks_ini; ++iv) {
               creal vx = P::vxmin + (iv+0.5) * SpatialCell::get_velocity_base_grid_block_size()[0]; // vx-coordinate of the centre
               creal vy = P::vymin + (jv+0.5) * SpatialCell::get_velocity_base_grid_block_size()[1]; // vy-
               creal vz = P::vzmin + (kv+0.5) * SpatialCell::get_velocity_base_grid_block_size()[2]; // vz-
               
               if (vx*vx + vy*vy + vz*vz < vRadiusSquared) {
                  cell.add_velocity_block(cell.get_velocity_block(vx, vy, vz));
                  blocksToInitialize.push_back(cell.get_velocity_block(vx, vy, vz));
               }
            }
            
      return blocksToInitialize;
   }
   
   void Ionosphere::setCellFromTemplate(SpatialCell *cell) {
      // The ionospheric cell has the same state as the initial state of non-system boundary cells so far.
      cell->parameters[CellParams::RHOLOSSADJUST] = 0.0;
      cell->parameters[CellParams::RHOLOSSVELBOUNDARY] = 0.0;
      //Copy, and allow to change blocks
      copyCellData(&templateCell, cell,true);
   }
   
   std::string Ionosphere::getName() const {return "Ionosphere";}
   
   uint Ionosphere::getIndex() const {return sysboundarytype::IONOSPHERE;}
}
