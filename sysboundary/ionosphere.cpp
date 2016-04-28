/*
 * This file is part of Vlasiator.
 * 
 * Copyright 2010-2015 Finnish Meteorological Institute
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
#include "../object_wrapper.h"

#ifndef NDEBUG
   #define DEBUG_IONOSPHERE
#endif
#ifdef DEBUG_SYSBOUNDARY
   #define DEBUG_IONOSPHERE
#endif

namespace SBC {
   Ionosphere::Ionosphere(): SysBoundaryCondition() { }
   
   Ionosphere::~Ionosphere() { }
   
   void Ionosphere::addParameters() {
      Readparameters::add("ionosphere.centerX", "X coordinate of ionosphere center (m)", 0.0);
      Readparameters::add("ionosphere.centerY", "Y coordinate of ionosphere center (m)", 0.0);
      Readparameters::add("ionosphere.centerZ", "Z coordinate of ionosphere center (m)", 0.0);
      Readparameters::add("ionosphere.radius", "Radius of ionosphere (m).", 1.0e7);
      Readparameters::add("ionosphere.geometry", "Select the geometry of the ionosphere, 0: inf-norm (diamond), 1: 1-norm (square), 2: 2-norm (circle, DEFAULT), 3: 2-norm cylinder aligned with y-axis, use with polar plane/line dipole.", 2);
      Readparameters::add("ionosphere.rho", "Number density of the ionosphere (m^-3)", 1.0e6);
      Readparameters::add("ionosphere.VX0", "Bulk velocity of ionospheric distribution function in X direction (m/s)", 0.0);
      Readparameters::add("ionosphere.VY0", "Bulk velocity of ionospheric distribution function in X direction (m/s)", 0.0);
      Readparameters::add("ionosphere.VZ0", "Bulk velocity of ionospheric distribution function in X direction (m/s)", 0.0);
      Readparameters::add("ionosphere.taperRadius", "Width of the zone with a density tapering from the ionospheric value to the background (m)", 0.0);
      Readparameters::add("ionosphere.precedence", "Precedence value of the ionosphere system boundary condition (integer), the higher the stronger.", 2);
      Readparameters::add("ionosphere.reapplyUponRestart", "If 0 (default), keep going with the state existing in the restart file. If 1, calls again applyInitialState. Can be used to change boundary condition behaviour during a run.", 0);
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
      if(!Readparameters::get("ionosphere.geometry", this->geometry)) {
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
      uint reapply;
      if(!Readparameters::get("ionosphere.reapplyUponRestart",reapply)) {
         if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
         exit(1);
      };
      this->applyUponRestart = false;
      if(reapply == 1) {
         this->applyUponRestart = true;
      }
   }
   
   bool Ionosphere::initSysBoundary(
      creal& t,
      Project &project
   ) {
      getParameters();
      isThisDynamic = false;

      // iniSysBoundary is only called once, generateTemplateCell must 
      // init all particle species
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
         Real r;
         
         switch(this->geometry) {
            case 0:
               // infinity-norm, result is a diamond/square with diagonals aligned on the axes in 2D
               r = fabs(x-center[0]) + fabs(y-center[1]) + fabs(z-center[2]);
               break;
            case 1:
               // 1-norm, result is is a grid-aligned square in 2D
               r = max(max(fabs(x-center[0]), fabs(y-center[1])), fabs(z-center[2]));
               break;
            case 2:
               // 2-norm (Cartesian), result is a circle in 2D
               r = sqrt((x-center[0])*(x-center[0]) + (y-center[1])*(y-center[1]) + (z-center[2])*(z-center[2]));
               break;
            case 3:
               // 2-norm (Cartesian) cylinder aligned on y-axis
               r = sqrt((x-center[0])*(x-center[0]) + (z-center[2])*(z-center[2]));
               break;
            default:
               std::cerr << __FILE__ << ":" << __LINE__ << ":" << "ionosphere.geometry has to be 0, 1 or 2." << std::endl;
               abort();
         }
         
         if(r < this->radius) {
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
         if (cell->sysBoundaryFlag != this->getIndex()) continue;
         
         for (unsigned int popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID)
            setCellFromTemplate(cell,popID);
      }
      return true;
   }

   std::array<Real, 3> Ionosphere::fieldSolverGetNormalDirection(
      const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const CellID& cellID
   ) {
      phiprof::start("Ionosphere::fieldSolverGetNormalDirection");
      std::array<Real, 3> normalDirection{{ 0.0, 0.0, 0.0 }};
      
      static creal DIAG2 = 1.0 / sqrt(2.0);
      static creal DIAG3 = 1.0 / sqrt(3.0);
      
      creal dx = mpiGrid[cellID]->parameters[CellParams::DX];
      creal dy = mpiGrid[cellID]->parameters[CellParams::DY];
      creal dz = mpiGrid[cellID]->parameters[CellParams::DZ];
      creal x = mpiGrid[cellID]->parameters[CellParams::XCRD] + 0.5*dx;
      creal y = mpiGrid[cellID]->parameters[CellParams::YCRD] + 0.5*dy;
      creal z = mpiGrid[cellID]->parameters[CellParams::ZCRD] + 0.5*dz;
      creal xsign = divideIfNonZero(x, fabs(x));
      creal ysign = divideIfNonZero(y, fabs(y));
      creal zsign = divideIfNonZero(z, fabs(z));
      
      Real length = 0.0;
      
      if (Parameters::xcells_ini == 1) {
         if (Parameters::ycells_ini == 1) {
            if (Parameters::zcells_ini == 1) {
               // X,Y,Z
               std::cerr << __FILE__ << ":" << __LINE__ << ":" << "What do you expect to do with a single-cell simulation of ionosphere boundary type? Stop kidding." << std::endl;
               abort();
               // end of X,Y,Z
            } else {
               // X,Y
               normalDirection[2] = zsign;
               // end of X,Y
            }
         } else if (Parameters::zcells_ini == 1) {
            // X,Z
            normalDirection[1] = ysign;
            // end of X,Z
         } else {
            // X
            switch(this->geometry) {
               case 0:
                  normalDirection[1] = DIAG2*ysign;
                  normalDirection[2] = DIAG2*zsign;
                  break;
               case 1:
                  if(fabs(y) == fabs(z)) {
                     normalDirection[1] = ysign*DIAG2;
                     normalDirection[2] = zsign*DIAG2;
                     break;
                  }
                  if(fabs(y) > (this->radius - dy)) {
                     normalDirection[1] = ysign;
                     break;
                  }
                  if(fabs(z) > (this->radius - dz)) {
                     normalDirection[2] = zsign;
                     break;
                  }
                  if(fabs(y) > (this->radius - 2.0*dy)) {
                     normalDirection[1] = ysign;
                     break;
                  }
                  if(fabs(z) > (this->radius - 2.0*dz)) {
                     normalDirection[2] = zsign;
                     break;
                  }
                  break;
               case 2:
                  length = sqrt(y*y + z*z);
                  normalDirection[1] = y / length;
                  normalDirection[2] = z / length;
                  break;
               default:
                  std::cerr << __FILE__ << ":" << __LINE__ << ":" << "ionosphere.geometry has to be 0, 1 or 2 with this grid shape." << std::endl;
                  abort();
            }
            // end of X
         }
      } else if (Parameters::ycells_ini == 1) {
         if (Parameters::zcells_ini == 1) {
            // Y,Z
            normalDirection[0] = xsign;
            // end of Y,Z
         } else {
            // Y
            switch(this->geometry) {
               case 0:
                  normalDirection[0] = DIAG2*xsign;
                  normalDirection[2] = DIAG2*zsign;
                  break;
               case 1:
                  if(fabs(x) == fabs(z)) {
                     normalDirection[0] = xsign*DIAG2;
                     normalDirection[2] = zsign*DIAG2;
                     break;
                  }
                  if(fabs(x) > (this->radius - dx)) {
                     normalDirection[0] = xsign;
                     break;
                  }
                  if(fabs(z) > (this->radius - dz)) {
                     normalDirection[2] = zsign;
                     break;
                  }
                  if(fabs(x) > (this->radius - 2.0*dx)) {
                     normalDirection[0] = xsign;
                     break;
                  }
                  if(fabs(z) > (this->radius - 2.0*dz)) {
                     normalDirection[2] = zsign;
                     break;
                  }
                  break;
               case 2:
               case 3:
                  length = sqrt(x*x + z*z);
                  normalDirection[0] = x / length;
                  normalDirection[2] = z / length;
                  break;
               default:
                  std::cerr << __FILE__ << ":" << __LINE__ << ":" << "ionosphere.geometry has to be 0, 1, 2 or 3 with this grid shape." << std::endl;
                  abort();
            }
            // end of Y
         }
      } else if (Parameters::zcells_ini == 1) {
         // Z
         switch(this->geometry) {
            case 0:
               normalDirection[0] = DIAG2*xsign;
               normalDirection[1] = DIAG2*ysign;
               break;
            case 1:
               if(fabs(x) == fabs(y)) {
                  normalDirection[0] = xsign*DIAG2;
                  normalDirection[1] = ysign*DIAG2;
                  break;
               }
               if(fabs(x) > (this->radius - dx)) {
                  normalDirection[0] = xsign;
                  break;
               }
               if(fabs(y) > (this->radius - dy)) {
                  normalDirection[1] = ysign;
                  break;
               }
               if(fabs(x) > (this->radius - 2.0*dx)) {
                  normalDirection[0] = xsign;
                  break;
               }
               if(fabs(y) > (this->radius - 2.0*dy)) {
                  normalDirection[1] = ysign;
                  break;
               }
               break;
            case 2:
               length = sqrt(x*x + y*y);
               normalDirection[0] = x / length;
               normalDirection[1] = y / length;
               break;
            default:
               std::cerr << __FILE__ << ":" << __LINE__ << ":" << "ionosphere.geometry has to be 0, 1 or 2 with this grid shape." << std::endl;
               abort();
         }
         // end of Z
      } else {
         // 3D
         switch(this->geometry) {
            case 0:
               normalDirection[0] = DIAG3*xsign;
               normalDirection[1] = DIAG3*ysign;
               normalDirection[2] = DIAG3*zsign;
               break;
            case 1:
               if(fabs(x) == fabs(y) && fabs(x) == fabs(z) && fabs(x) > this->radius - dx) {
                  normalDirection[0] = xsign*DIAG3;
                  normalDirection[1] = ysign*DIAG3;
                  normalDirection[2] = zsign*DIAG3;
                  break;
               }
               if(fabs(x) == fabs(y) && fabs(x) == fabs(z) && fabs(x) > this->radius - 2.0*dx) {
                  normalDirection[0] = xsign*DIAG3;
                  normalDirection[1] = ysign*DIAG3;
                  normalDirection[2] = zsign*DIAG3;
                  break;
               }
               if(fabs(x) == fabs(y) && fabs(x) > this->radius - dx && fabs(z) < this->radius - dz) {
                  normalDirection[0] = xsign*DIAG2;
                  normalDirection[1] = ysign*DIAG2;
                  normalDirection[2] = 0.0;
                  break;
               }
               if(fabs(y) == fabs(z) && fabs(y) > this->radius - dy && fabs(x) < this->radius - dx) {
                  normalDirection[0] = 0.0;
                  normalDirection[1] = ysign*DIAG2;
                  normalDirection[2] = zsign*DIAG2;
                  break;
               }
               if(fabs(x) == fabs(z) && fabs(x) > this->radius - dx && fabs(y) < this->radius - dy) {
                  normalDirection[0] = xsign*DIAG2;
                  normalDirection[1] = 0.0;
                  normalDirection[2] = zsign*DIAG2;
                  break;
               }
               if(fabs(x) == fabs(y) && fabs(x) > this->radius - 2.0*dx && fabs(z) < this->radius - 2.0*dz) {
                  normalDirection[0] = xsign*DIAG2;
                  normalDirection[1] = ysign*DIAG2;
                  normalDirection[2] = 0.0;
                  break;
               }
               if(fabs(y) == fabs(z) && fabs(y) > this->radius - 2.0*dy && fabs(x) < this->radius - 2.0*dx) {
                  normalDirection[0] = 0.0;
                  normalDirection[1] = ysign*DIAG2;
                  normalDirection[2] = zsign*DIAG2;
                  break;
               }
               if(fabs(x) == fabs(z) && fabs(x) > this->radius - 2.0*dx && fabs(y) < this->radius - 2.0*dy) {
                  normalDirection[0] = xsign*DIAG2;
                  normalDirection[1] = 0.0;
                  normalDirection[2] = zsign*DIAG2;
                  break;
               }
               if(fabs(x) > (this->radius - dx)) {
                  normalDirection[0] = xsign;
                  break;
               }
               if(fabs(y) > (this->radius - dy)) {
                  normalDirection[1] = ysign;
                  break;
               }
               if(fabs(z) > (this->radius - dz)) {
                  normalDirection[2] = zsign;
                  break;
               }
               if(fabs(x) > (this->radius - 2.0*dx)) {
                  normalDirection[0] = xsign;
                  break;
               }
               if(fabs(y) > (this->radius - 2.0*dy)) {
                  normalDirection[1] = ysign;
                  break;
               }
               if(fabs(z) > (this->radius - 2.0*dz)) {
                  normalDirection[2] = zsign;
                  break;
               }
               break;
            case 2:
               length = sqrt(x*x + y*y + z*z);
               normalDirection[0] = x / length;
               normalDirection[1] = y / length;
               normalDirection[2] = z / length;
               break;
            case 3:
               length = sqrt(x*x + z*z);
               normalDirection[0] = x / length;
               normalDirection[2] = z / length;
               break;
            default:
               std::cerr << __FILE__ << ":" << __LINE__ << ":" << "ionosphere.geometry has to be 0, 1, 2 or 3 with this grid shape." << std::endl;
               abort();
         }
         // end of 3D
      }
      
#ifdef DEBUG_IONOSPHERE
      // Uncomment one of the following line for debugging output to evaluate the correctness of the results. Best used with a single process and single thread.
//       if (mpiGrid[cellID]->sysBoundaryLayer == 1) std::cerr << x << " " << y << " " << z << " " << normalDirection[0] << " " << normalDirection[1] << " " << normalDirection[2] << std::endl;
//       if (mpiGrid[cellID]->sysBoundaryLayer == 2) std::cerr << x << " " << y << " " << z << " " << normalDirection[0] << " " << normalDirection[1] << " " << normalDirection[2] << std::endl;
//      std::cerr << x << " " << y << " " << z << " " << normalDirection[0] << " " << normalDirection[1] << " " << normalDirection[2] << std::endl;
#endif
      phiprof::stop("Ionosphere::fieldSolverGetNormalDirection");
      return normalDirection;
   }
   
   /*! We want here to
    * 
    * -- Average perturbed face B from the nearest neighbours
    * 
    * -- Retain only the normal components of perturbed face B
    */
   Real Ionosphere::fieldSolverBoundaryCondMagneticField(
      const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const std::vector<fs_cache::CellCache>& cellCache,
      const uint16_t& localID,
      creal& dt,
      cuint& RKCase,
      cint& offset,
      cuint& component
   ) {
      const CellID cellID = cellCache[localID].cellID;
      std::vector<CellID> closestCells = getAllClosestNonsysboundaryCells(cellID);
      if (closestCells.size() == 1 && closestCells[0] == INVALID_CELLID) {
         std::cerr << __FILE__ << ":" << __LINE__ << ":" << "No closest cells found!" << std::endl;
         abort();
      }

      // Sum perturbed B component over all nearest NOT_SYSBOUNDARY neighbours
      std::array<Real, 3> averageB = {{ 0.0 }};
      for (uint i=0; i<closestCells.size(); i++) {
         #ifdef DEBUG_IONOSPHERE
         if (mpiGrid[closestCells[i]]->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY) {
            stringstream ss;
            ss << "ERROR, ionosphere cell " << cellID << " uses value from sysboundary nbr " << closestCells[i];
            ss << " in " << __FILE__ << ":" << __LINE__ << endl;
            cerr << ss.str();
            exit(1);
         }
         #endif
         averageB[0] += mpiGrid[closestCells[i]]->parameters[CellParams::PERBX+offset];
         averageB[1] += mpiGrid[closestCells[i]]->parameters[CellParams::PERBY+offset];
         averageB[2] += mpiGrid[closestCells[i]]->parameters[CellParams::PERBZ+offset];
      }

      // Average and project to normal direction
      std::array<Real, 3> normalDirection = fieldSolverGetNormalDirection(mpiGrid, cellID);
      for(uint i=0; i<3; i++) {
         averageB[i] *= normalDirection[i] / closestCells.size();
      }

      // Return (B.n)*normalVector[component]
      return (averageB[0]+averageB[1]+averageB[2])*normalDirection[component];
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
      fs_cache::CellCache& cache,
      cuint RKCase,
      cuint component
   ) {
      Real* cp = cache.cells[fs_cache::calculateNbrID(1,1,1)]->parameters;
      
      switch(component) {
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
   
   void Ionosphere::fieldSolverBoundaryCondGradPeElectricField(
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
   
   void Ionosphere::fieldSolverBoundaryCondDerivatives(
      dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const CellID& cellID,
      cuint& RKCase,
      cuint& component
   ) {
      this->setCellDerivativesToZero(mpiGrid, cellID, component);
      return;
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
      const CellID& cellID,
      const int& popID
   ) {
//       phiprof::start("vlasovBoundaryCondition (Ionosphere)");
//       const SpatialCell * cell = mpiGrid[cellID];
//       this->vlasovBoundaryCopyFromAllClosestNbrs(mpiGrid, cellID);
//       phiprof::stop("vlasovBoundaryCondition (Ionosphere)");
   }

   /**
    * NOTE: This function must initialize all particle species!
    * @param project
    */
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
      
      // Loop over particle species
      for (unsigned int popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
         const vector<vmesh::GlobalID> blocksToInitialize = findBlocksToInitialize(templateCell,popID);
         Realf* data = templateCell.get_data(popID);
         
         for (size_t i = 0; i < blocksToInitialize.size(); i++) {
            const vmesh::GlobalID blockGID = blocksToInitialize.at(i);
            const vmesh::LocalID blockLID = templateCell.get_velocity_block_local_id(blockGID,popID);
            const Real* block_parameters = templateCell.get_block_parameters(blockLID,popID);
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
            for (uint kc=0; kc<WID; ++kc) for (uint jc=0; jc<WID; ++jc) for (uint ic=0; ic<WID; ++ic) {
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
                           average +=  shiftedMaxwellianDistribution(
                                                                     popID,
                                                                     vxCell + vi*d_vx,
                                                                     vyCell + vj*d_vy,
                                                                     vzCell + vk*d_vz
                                                                    );
                        }
                  average /= this->nVelocitySamples * this->nVelocitySamples * this->nVelocitySamples;
               } else {
                  average = shiftedMaxwellianDistribution(
                                                          popID,
                                                          vxCell + 0.5*dvxCell,
                                                          vyCell + 0.5*dvyCell,
                                                          vzCell + 0.5*dvzCell
                                                         );
               }

               if (average !=0.0 ) {
                  data[blockLID*WID3+cellIndex(ic,jc,kc)] = average;
               }
            } // for-loop over cells in velocity block
         } // for-loop over velocity blocks

         // let's get rid of blocks not fulfilling the criteria here to save memory.
         templateCell.adjustSingleCellVelocityBlocks(popID);
      } // for-loop over particle species
      
      calculateCellMoments(&templateCell,true,true);

      // WARNING Time-independence assumed here. Normal moments computed in setProjectCell
      templateCell.parameters[CellParams::RHO_DT2] = templateCell.parameters[CellParams::RHO];
      templateCell.parameters[CellParams::RHOVX_DT2] = templateCell.parameters[CellParams::RHOVX];
      templateCell.parameters[CellParams::RHOVY_DT2] = templateCell.parameters[CellParams::RHOVY];
      templateCell.parameters[CellParams::RHOVZ_DT2] = templateCell.parameters[CellParams::RHOVZ];
   }
   
   Real Ionosphere::shiftedMaxwellianDistribution(
      const int& popID,
      creal& vx, creal& vy, creal& vz
   ) {
      
      #warning All species have the same VX0,VY0,VZ0,T
      const Real MASS = getObjectWrapper().particleSpecies[popID].mass;

      return this->rho * pow(MASS /
      (2.0 * M_PI * physicalconstants::K_B * this->T), 1.5) *
      exp(-MASS * ((vx-this->VX0)*(vx-this->VX0) + (vy-this->VY0)*(vy-this->VY0) + (vz-this->VZ0)*(vz-this->VZ0)) /
      (2.0 * physicalconstants::K_B * this->T));
   }

   std::vector<vmesh::GlobalID> Ionosphere::findBlocksToInitialize(spatial_cell::SpatialCell& cell,const int& popID) {
      vector<vmesh::GlobalID> blocksToInitialize;
      bool search = true;
      uint counter = 0;
      const uint8_t refLevel = 0;

      const vmesh::LocalID* vblocks_ini = cell.get_velocity_grid_length(popID,refLevel);

      while (search) {
         #warning TODO: add SpatialCell::getVelocityBlockMinValue() in place of sparseMinValue ? (if applicable)
         if (0.1 * getObjectWrapper().particleSpecies[popID].sparseMinValue > 
            shiftedMaxwellianDistribution(popID,counter*cell.get_velocity_grid_block_size(popID,refLevel)[0], 0.0, 0.0)
            || counter > vblocks_ini[0]) {
            search = false;
         }
         ++counter;
      }
      counter+=2;
      Real vRadiusSquared 
              = (Real)counter*(Real)counter
              * cell.get_velocity_grid_block_size(popID,refLevel)[0]
              * cell.get_velocity_grid_block_size(popID,refLevel)[0];

      for (uint kv=0; kv<vblocks_ini[2]; ++kv) 
         for (uint jv=0; jv<vblocks_ini[1]; ++jv)
            for (uint iv=0; iv<vblocks_ini[0]; ++iv) {
               vmesh::LocalID blockIndices[3];
               blockIndices[0] = iv;
               blockIndices[1] = jv;
               blockIndices[2] = kv;
               const vmesh::GlobalID blockGID = cell.get_velocity_block(popID,blockIndices,refLevel);
               Real blockCoords[3];
               cell.get_velocity_block_coordinates(popID,blockGID,blockCoords);
               Real blockSize[3];
               cell.get_velocity_block_size(popID,blockGID,blockSize);
               blockCoords[0] += 0.5*blockSize[0];
               blockCoords[1] += 0.5*blockSize[1];
               blockCoords[2] += 0.5*blockSize[2];
               //creal vx = P::vxmin + (iv+0.5) * cell.get_velocity_grid_block_size(popID)[0]; // vx-coordinate of the centre
               //creal vy = P::vymin + (jv+0.5) * cell.get_velocity_grid_block_size(popID)[1]; // vy-
               //creal vz = P::vzmin + (kv+0.5) * cell.get_velocity_grid_block_size(popID)[2]; // vz-
               
               if (blockCoords[0]*blockCoords[0] + blockCoords[1]*blockCoords[1] + blockCoords[2]*blockCoords[2] < vRadiusSquared) {
               //if (vx*vx + vy*vy + vz*vz < vRadiusSquared) {
                  // Adds velocity block to active population's velocity mesh
                  //const vmesh::GlobalID newBlockGID = cell.get_velocity_block(popID,vx,vy,vz);
                  cell.add_velocity_block(blockGID,popID);
                  blocksToInitialize.push_back(blockGID);
               }
            }

      return blocksToInitialize;
   }

   void Ionosphere::setCellFromTemplate(SpatialCell* cell,const int& popID) {
      // The ionospheric cell has the same state as the initial state of non-system boundary cells so far.
      if (popID == 0) {
         cell->parameters[CellParams::RHOLOSSADJUST] = 0.0;
         cell->parameters[CellParams::RHOLOSSVELBOUNDARY] = 0.0;
      }

      //Copy, and allow to change blocks
      copyCellData(&templateCell,cell,true,false,popID);
   }

   std::string Ionosphere::getName() const {return "Ionosphere";}
   
   uint Ionosphere::getIndex() const {return sysboundarytype::IONOSPHERE;}
}
