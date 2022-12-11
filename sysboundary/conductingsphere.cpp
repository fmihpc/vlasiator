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

/*!\file conductingsphere.cpp
 * \brief Implementation of the class SysBoundaryCondition::Conductingsphere to handle cells classified as sysboundarytype::CONDUCTINGSPHERE.
 */

#include <cstdlib>
#include <iostream>

#include "conductingsphere.h"
#include "../projects/project.h"
#include "../projects/projects_common.h"
#include "../vlasovmover.h"
#include "../fieldsolver/fs_common.h"
#include "../fieldsolver/fs_limiters.h"
#include "../fieldsolver/ldz_magnetic_field.hpp"
#include "../common.h"
#include "../object_wrapper.h"

#ifndef NDEBUG
   #define DEBUG_CONDUCTINGSPHERE
#endif
#ifdef DEBUG_SYSBOUNDARY
   #define DEBUG_CONDUCTINGSPHERE
#endif

namespace SBC {
   Conductingsphere::Conductingsphere(): SysBoundaryCondition() { }
   
   Conductingsphere::~Conductingsphere() { }
   
   void Conductingsphere::addParameters() {
      Readparameters::add("conductingsphere.centerX", "X coordinate of conductingsphere center (m)", 0.0);
      Readparameters::add("conductingsphere.centerY", "Y coordinate of conductingsphere center (m)", 0.0);
      Readparameters::add("conductingsphere.centerZ", "Z coordinate of conductingsphere center (m)", 0.0);
      Readparameters::add("conductingsphere.radius", "Radius of conductingsphere (m).", 1.0e7);
      Readparameters::add("conductingsphere.geometry", "Select the geometry of the conductingsphere, 0: inf-norm (diamond), 1: 1-norm (square), 2: 2-norm (circle, DEFAULT), 3: 2-norm cylinder aligned with y-axis, use with polar plane/line dipole.", 2);
      Readparameters::add("conductingsphere.precedence", "Precedence value of the conductingsphere system boundary condition (integer), the higher the stronger.", 2);
      Readparameters::add("conductingsphere.reapplyUponRestart", "If 0 (default), keep going with the state existing in the restart file. If 1, calls again applyInitialState. Can be used to change boundary condition behaviour during a run.", 0);

      // Per-population parameters
      for(uint i=0; i< getObjectWrapper().particleSpecies.size(); i++) {
         const std::string& pop = getObjectWrapper().particleSpecies[i].name;
         
         Readparameters::add(pop + "_conductingsphere.taperRadius", "Width of the zone with a density tapering from the conducting sphere value to the background (m)", 0.0);
         Readparameters::add(pop + "_conductingsphere.rho", "Number density of the conductingsphere (m^-3)", 0.0);
         Readparameters::add(pop + "_conductingsphere.T", "Temperature of the conductingsphere (K)", 0.0);
         Readparameters::add(pop + "_conductingsphere.VX0", "Bulk velocity of conductospheric distribution function in X direction (m/s)", 0.0);
         Readparameters::add(pop + "_conductingsphere.VY0", "Bulk velocity of conductospheric distribution function in X direction (m/s)", 0.0);
         Readparameters::add(pop + "_conductingsphere.VZ0", "Bulk velocity of conductospheric distribution function in X direction (m/s)", 0.0);
         Readparameters::add(pop + "_conductingsphere.fluffiness", "Inertia of boundary smoothing when copying neighbour's moments and velocity distributions (0=completely constant boundaries, 1=neighbours are interpolated immediately).", 0);
      }
   }
   
   void Conductingsphere::getParameters() {

      Readparameters::get("conductingsphere.centerX", this->center[0]);
      Readparameters::get("conductingsphere.centerY", this->center[1]);
      Readparameters::get("conductingsphere.centerZ", this->center[2]);
      Readparameters::get("conductingsphere.radius", this->radius);
      Readparameters::get("conductingsphere.geometry", this->geometry);
      Readparameters::get("conductingsphere.precedence", this->precedence);
      uint reapply;
      Readparameters::get("conductingsphere.reapplyUponRestart",reapply);
      this->applyUponRestart = false;
      if(reapply == 1) {
         this->applyUponRestart = true;
      }

      for(uint i=0; i< getObjectWrapper().particleSpecies.size(); i++) {
        const std::string& pop = getObjectWrapper().particleSpecies[i].name;
        ConductingsphereSpeciesParameters sP;

        Readparameters::get(pop + "_conductingsphere.rho", sP.rho);
        Readparameters::get(pop + "_conductingsphere.VX0", sP.V0[0]);
        Readparameters::get(pop + "_conductingsphere.VY0", sP.V0[1]);
        Readparameters::get(pop + "_conductingsphere.VZ0", sP.V0[2]);
        Readparameters::get(pop + "_conductingsphere.fluffiness", sP.fluffiness);
        Readparameters::get(pop + "_conductingsphere.T", sP.T);
        Readparameters::get(pop + "_Magnetosphere.nSpaceSamples", sP.nSpaceSamples);
        Readparameters::get(pop + "_Magnetosphere.nVelocitySamples", sP.nVelocitySamples);

        // Failsafe, if density or temperature is zero, read from Magnetosphere
        // (compare the corresponding verbose handling in projects/Magnetosphere/Magnetosphere.cpp)
        if(sP.T == 0) {
            Readparameters::get(pop + "_Magnetosphere.T", sP.T);
         }
         if(sP.rho == 0) {
            Readparameters::get(pop + "_Magnetosphere.rho", sP.rho);
         }

         speciesParams.push_back(sP);
      }
   }
   
   bool Conductingsphere::initSysBoundary(
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

   Real getR(creal x,creal y,creal z, uint geometry, Real center[3]) {

      Real r;
      
      switch(geometry) {
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
         std::cerr << __FILE__ << ":" << __LINE__ << ":" << "conductingsphere.geometry has to be 0, 1 or 2." << std::endl;
         abort();
      }

      return r;
   }
   
   bool Conductingsphere::assignSysBoundary(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid) {
      const vector<CellID>& cells = getLocalCells();
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
         
         if(getR(x,y,z,this->geometry,this->center) < this->radius) {
            mpiGrid[cells[i]]->sysBoundaryFlag = this->getIndex();
         }
      }

      return true;
   }

   bool Conductingsphere::applyInitialState(
      const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
      Project &project
   ) {
      const vector<CellID>& cells = getLocalCells();
      #pragma omp parallel for
      for (uint i=0; i<cells.size(); ++i) {
         SpatialCell* cell = mpiGrid[cells[i]];
         if (cell->sysBoundaryFlag != this->getIndex()) continue;
         
         for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID)
            setCellFromTemplate(cell,popID);
      }
      return true;
   }

   std::array<Real, 3> Conductingsphere::fieldSolverGetNormalDirection(
      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
      cint i,
      cint j,
      cint k
   ) {
      phiprof::start("Conductingsphere::fieldSolverGetNormalDirection");
      std::array<Real, 3> normalDirection{{ 0.0, 0.0, 0.0 }};
      
      static creal DIAG2 = 1.0 / sqrt(2.0);
      static creal DIAG3 = 1.0 / sqrt(3.0);
      
      creal dx = technicalGrid.DX;
      creal dy = technicalGrid.DY;
      creal dz = technicalGrid.DZ;
      const std::array<int, 3> globalIndices = technicalGrid.getGlobalIndices(i,j,k);
      creal x = P::xmin + (convert<Real>(globalIndices[0])+0.5)*dx;
      creal y = P::ymin + (convert<Real>(globalIndices[1])+0.5)*dy;
      creal z = P::zmin + (convert<Real>(globalIndices[2])+0.5)*dz;
      creal xsign = divideIfNonZero(x, fabs(x));
      creal ysign = divideIfNonZero(y, fabs(y));
      creal zsign = divideIfNonZero(z, fabs(z));
      
      Real length = 0.0;
      
      if (Parameters::xcells_ini == 1) {
         if (Parameters::ycells_ini == 1) {
            if (Parameters::zcells_ini == 1) {
               // X,Y,Z
               std::cerr << __FILE__ << ":" << __LINE__ << ":" << "What do you expect to do with a single-cell simulation of conductingsphere boundary type? Stop kidding." << std::endl;
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
                  std::cerr << __FILE__ << ":" << __LINE__ << ":" << "conductingsphere.geometry has to be 0, 1 or 2 with this grid shape." << std::endl;
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
                  std::cerr << __FILE__ << ":" << __LINE__ << ":" << "conductingsphere.geometry has to be 0, 1, 2 or 3 with this grid shape." << std::endl;
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
               std::cerr << __FILE__ << ":" << __LINE__ << ":" << "conductingsphere.geometry has to be 0, 1 or 2 with this grid shape." << std::endl;
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
               std::cerr << __FILE__ << ":" << __LINE__ << ":" << "conductingsphere.geometry has to be 0, 1, 2 or 3 with this grid shape." << std::endl;
               abort();
         }
         // end of 3D
      }
      
      phiprof::stop("Conductingsphere::fieldSolverGetNormalDirection");
      return normalDirection;
   }
   
   /*! We want here to
    * 
    * -- Average perturbed face B from the nearest neighbours
    * 
    * -- Retain only the normal components of perturbed face B
    */
   Real Conductingsphere::fieldSolverBoundaryCondMagneticField(
      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & bGrid,
      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
      cint i,
      cint j,
      cint k,
      creal& dt,
      cuint& component
   ) {
      if (technicalGrid.get(i,j,k)->sysBoundaryLayer == 1) {
         switch(component) {
            case 0:
               if (  ((technicalGrid.get(i-1,j,k)->SOLVE & compute::BX) == compute::BX)
                  && ((technicalGrid.get(i+1,j,k)->SOLVE & compute::BX) == compute::BX)
               ) {
                  return 0.5 * (bGrid.get(i-1,j,k)->at(fsgrids::bfield::PERBX) + bGrid.get(i+1,j,k)->at(fsgrids::bfield::PERBX));
               } else if ((technicalGrid.get(i-1,j,k)->SOLVE & compute::BX) == compute::BX) {
                  return bGrid.get(i-1,j,k)->at(fsgrids::bfield::PERBX);
               } else if ((technicalGrid.get(i+1,j,k)->SOLVE & compute::BX) == compute::BX) {
                  return bGrid.get(i+1,j,k)->at(fsgrids::bfield::PERBX);
               } else {
                  Real retval = 0.0;
                  uint nCells = 0;
                  if ((technicalGrid.get(i,j-1,k)->SOLVE & compute::BX) == compute::BX) {
                     retval += bGrid.get(i,j-1,k)->at(fsgrids::bfield::PERBX);
                     nCells++;
                  }
                  if ((technicalGrid.get(i,j+1,k)->SOLVE & compute::BX) == compute::BX) {
                     retval += bGrid.get(i,j+1,k)->at(fsgrids::bfield::PERBX);
                     nCells++;
                  }
                  if ((technicalGrid.get(i,j,k-1)->SOLVE & compute::BX) == compute::BX) {
                     retval += bGrid.get(i,j,k-1)->at(fsgrids::bfield::PERBX);
                     nCells++;
                  }
                  if ((technicalGrid.get(i,j,k+1)->SOLVE & compute::BX) == compute::BX) {
                     retval += bGrid.get(i,j,k+1)->at(fsgrids::bfield::PERBX);
                     nCells++;
                  }
                  if (nCells == 0) {
                     for (int a=i-1; a<i+2; a++) {
                        for (int b=j-1; b<j+2; b++) {
                           for (int c=k-1; c<k+2; c++) {
                              if ((technicalGrid.get(a,b,c)->SOLVE & compute::BX) == compute::BX) {
                                 retval += bGrid.get(a,b,c)->at(fsgrids::bfield::PERBX);
                                 nCells++;
                              }
                           }
                        }
                     }
                  }
                  if (nCells == 0) {
                     cerr << __FILE__ << ":" << __LINE__ << ": ERROR: this should not have fallen through." << endl;
                     return 0.0;
                  }
                  return retval / nCells;
               }
            case 1:
               if (  (technicalGrid.get(i,j-1,k)->SOLVE & compute::BY) == compute::BY
                  && (technicalGrid.get(i,j+1,k)->SOLVE & compute::BY) == compute::BY
               ) {
                  return 0.5 * (bGrid.get(i,j-1,k)->at(fsgrids::bfield::PERBY) + bGrid.get(i,j+1,k)->at(fsgrids::bfield::PERBY));
               } else if ((technicalGrid.get(i,j-1,k)->SOLVE & compute::BY) == compute::BY) {
                  return bGrid.get(i,j-1,k)->at(fsgrids::bfield::PERBY);
               } else if ((technicalGrid.get(i,j+1,k)->SOLVE & compute::BY) == compute::BY) {
                  return bGrid.get(i,j+1,k)->at(fsgrids::bfield::PERBY);
               } else {
                  Real retval = 0.0;
                  uint nCells = 0;
                  if ((technicalGrid.get(i-1,j,k)->SOLVE & compute::BY) == compute::BY) {
                     retval += bGrid.get(i-1,j,k)->at(fsgrids::bfield::PERBY);
                     nCells++;
                  }
                  if ((technicalGrid.get(i+1,j,k)->SOLVE & compute::BY) == compute::BY) {
                     retval += bGrid.get(i+1,j,k)->at(fsgrids::bfield::PERBY);
                     nCells++;
                  }
                  if ((technicalGrid.get(i,j,k-1)->SOLVE & compute::BY) == compute::BY) {
                     retval += bGrid.get(i,j,k-1)->at(fsgrids::bfield::PERBY);
                     nCells++;
                  }
                  if ((technicalGrid.get(i,j,k+1)->SOLVE & compute::BY) == compute::BY) {
                     retval += bGrid.get(i,j,k+1)->at(fsgrids::bfield::PERBY);
                     nCells++;
                  }
                  if (nCells == 0) {
                     for (int a=i-1; a<i+2; a++) {
                        for (int b=j-1; b<j+2; b++) {
                           for (int c=k-1; c<k+2; c++) {
                              if ((technicalGrid.get(a,b,c)->SOLVE & compute::BY) == compute::BY) {
                                 retval += bGrid.get(a,b,c)->at(fsgrids::bfield::PERBY);
                                 nCells++;
                              }
                           }
                        }
                     }
                  }
                  if (nCells == 0) {
                     cerr << __FILE__ << ":" << __LINE__ << ": ERROR: this should not have fallen through." << endl;
                     return 0.0;
                  }
                  return retval / nCells;
               }
            case 2:
               if (  (technicalGrid.get(i,j,k-1)->SOLVE & compute::BZ) == compute::BZ
                  && (technicalGrid.get(i,j,k+1)->SOLVE & compute::BZ) == compute::BZ
               ) {
                  return 0.5 * (bGrid.get(i,j,k-1)->at(fsgrids::bfield::PERBZ) + bGrid.get(i,j,k+1)->at(fsgrids::bfield::PERBZ));
               } else if ((technicalGrid.get(i,j,k-1)->SOLVE & compute::BZ) == compute::BZ) {
                  return bGrid.get(i,j,k-1)->at(fsgrids::bfield::PERBZ);
               } else if ((technicalGrid.get(i,j,k+1)->SOLVE & compute::BZ) == compute::BZ) {
                  return bGrid.get(i,j,k+1)->at(fsgrids::bfield::PERBZ);
               } else {
                  Real retval = 0.0;
                  uint nCells = 0;
                  if ((technicalGrid.get(i-1,j,k)->SOLVE & compute::BZ) == compute::BZ) {
                     retval += bGrid.get(i-1,j,k)->at(fsgrids::bfield::PERBZ);
                     nCells++;
                  }
                  if ((technicalGrid.get(i+1,j,k)->SOLVE & compute::BZ) == compute::BZ) {
                     retval += bGrid.get(i+1,j,k)->at(fsgrids::bfield::PERBZ);
                     nCells++;
                  }
                  if ((technicalGrid.get(i,j-1,k)->SOLVE & compute::BZ) == compute::BZ) {
                     retval += bGrid.get(i,j-1,k)->at(fsgrids::bfield::PERBZ);
                     nCells++;
                  }
                  if ((technicalGrid.get(i,j+1,k)->SOLVE & compute::BZ) == compute::BZ) {
                     retval += bGrid.get(i,j+1,k)->at(fsgrids::bfield::PERBZ);
                     nCells++;
                  }
                  if (nCells == 0) {
                     for (int a=i-1; a<i+2; a++) {
                        for (int b=j-1; b<j+2; b++) {
                           for (int c=k-1; c<k+2; c++) {
                              if ((technicalGrid.get(a,b,c)->SOLVE & compute::BZ) == compute::BZ) {
                                 retval += bGrid.get(a,b,c)->at(fsgrids::bfield::PERBZ);
                                 nCells++;
                              }
                           }
                        }
                     }
                  }
                  if (nCells == 0) {
                     cerr << __FILE__ << ":" << __LINE__ << ": ERROR: this should not have fallen through." << endl;
                     return 0.0;
                  }
                  return retval / nCells;
               }
            default:
               cerr << "ERROR: conductingsphere boundary tried to copy nonsensical magnetic field component " << component << endl;
               return 0.0;
         }
      } else { // L2 cells
         Real retval = 0.0;
         uint nCells = 0;
         for (int a=i-1; a<i+2; a++) {
            for (int b=j-1; b<j+2; b++) {
               for (int c=k-1; c<k+2; c++) {
                  if (technicalGrid.get(a,b,c)->sysBoundaryLayer == 1) {
                     retval += bGrid.get(a,b,c)->at(fsgrids::bfield::PERBX + component);
                     nCells++;
                  }
               }
            }
         }
         if (nCells == 0) {
            cerr << __FILE__ << ":" << __LINE__ << ": ERROR: this should not have fallen through." << endl;
            return 0.0;
         }
         return retval / nCells;
      }
   }

   /*! We want here to
    *
    * -- Retain only the boundary-normal projection of perturbed face B
    */
   void Conductingsphere::fieldSolverBoundaryCondMagneticFieldProjection(
      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & bGrid,
      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
      cint i,
      cint j,
      cint k
   ) {
      // Projection of B-field to normal direction
      Real BdotN = 0;
      std::array<Real, 3> normalDirection = fieldSolverGetNormalDirection(technicalGrid, i, j, k);
      for(uint component=0; component<3; component++) {
         BdotN += bGrid.get(i,j,k)->at(fsgrids::bfield::PERBX+component) * normalDirection[component];
      }
      // Apply to any components that were not solved
      if ((technicalGrid.get(i,j,k)->sysBoundaryLayer == 2) ||
          ((technicalGrid.get(i,j,k)->sysBoundaryLayer == 1) && ((technicalGrid.get(i,j,k)->SOLVE & compute::BX) != compute::BX))
         ) {
         bGrid.get(i,j,k)->at(fsgrids::bfield::PERBX) = BdotN*normalDirection[0];
      }
      if ((technicalGrid.get(i,j,k)->sysBoundaryLayer == 2) ||
          ((technicalGrid.get(i,j,k)->sysBoundaryLayer == 1) && ((technicalGrid.get(i,j,k)->SOLVE & compute::BY) != compute::BY))
         ) {
         bGrid.get(i,j,k)->at(fsgrids::bfield::PERBY) = BdotN*normalDirection[1];
      }
      if ((technicalGrid.get(i,j,k)->sysBoundaryLayer == 2) ||
          ((technicalGrid.get(i,j,k)->sysBoundaryLayer == 1) && ((technicalGrid.get(i,j,k)->SOLVE & compute::BZ) != compute::BZ))
         ) {
         bGrid.get(i,j,k)->at(fsgrids::bfield::PERBZ) = BdotN*normalDirection[2];
      }
   }

   void Conductingsphere::fieldSolverBoundaryCondElectricField(
      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
      cint i,
      cint j,
      cint k,
      cuint component
   ) {
      EGrid.get(i,j,k)->at(fsgrids::efield::EX+component) = 0.0;
   }
   
   void Conductingsphere::fieldSolverBoundaryCondHallElectricField(
      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
      cint i,
      cint j,
      cint k,
      cuint component
   ) {
      std::array<Real, fsgrids::ehall::N_EHALL> * cp = EHallGrid.get(i,j,k);
      switch (component) {
         case 0:
            cp->at(fsgrids::ehall::EXHALL_000_100) = 0.0;
            cp->at(fsgrids::ehall::EXHALL_010_110) = 0.0;
            cp->at(fsgrids::ehall::EXHALL_001_101) = 0.0;
            cp->at(fsgrids::ehall::EXHALL_011_111) = 0.0;
            break;
         case 1:
            cp->at(fsgrids::ehall::EYHALL_000_010) = 0.0;
            cp->at(fsgrids::ehall::EYHALL_100_110) = 0.0;
            cp->at(fsgrids::ehall::EYHALL_001_011) = 0.0;
            cp->at(fsgrids::ehall::EYHALL_101_111) = 0.0;
            break;
         case 2:
            cp->at(fsgrids::ehall::EZHALL_000_001) = 0.0;
            cp->at(fsgrids::ehall::EZHALL_100_101) = 0.0;
            cp->at(fsgrids::ehall::EZHALL_010_011) = 0.0;
            cp->at(fsgrids::ehall::EZHALL_110_111) = 0.0;
            break;
         default:
            cerr << __FILE__ << ":" << __LINE__ << ":" << " Invalid component" << endl;
      }
   }
   
   void Conductingsphere::fieldSolverBoundaryCondGradPeElectricField(
      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
      cint i,
      cint j,
      cint k,
      cuint component
   ) {
      EGradPeGrid.get(i,j,k)->at(fsgrids::egradpe::EXGRADPE+component) = 0.0;
   }
   
   void Conductingsphere::fieldSolverBoundaryCondDerivatives(
      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
      cint i,
      cint j,
      cint k,
      cuint& RKCase,
      cuint& component
   ) {
      this->setCellDerivativesToZero(dPerBGrid, dMomentsGrid, i, j, k, component);
      return;
   }
   
   void Conductingsphere::fieldSolverBoundaryCondBVOLDerivatives(
      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, 2> & volGrid,
      cint i,
      cint j,
      cint k,
      cuint& component
   ) {
      // FIXME This should be OK as the BVOL derivatives are only used for Lorentz force JXB, which is not applied on the conducting sphere cells.
      this->setCellBVOLDerivativesToZero(volGrid, i, j, k, component);
   }
   
   void Conductingsphere::vlasovBoundaryCondition(
      const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const CellID& cellID,
      const uint popID,
      const bool calculate_V_moments
   ) {
      phiprof::start("vlasovBoundaryCondition (Conductingsphere)");
      this->vlasovBoundaryFluffyCopyFromAllCloseNbrs(mpiGrid, cellID, popID, calculate_V_moments, this->speciesParams[popID].fluffiness);
      phiprof::stop("vlasovBoundaryCondition (Conductingsphere)");
   }

   /**
    * NOTE: This function must initialize all particle species!
    * @param project
    */
   void Conductingsphere::generateTemplateCell(Project &project) {
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
      for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
         const ConductingsphereSpeciesParameters& sP = this->speciesParams[popID];
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
               if(sP.nVelocitySamples > 1) {
                  creal d_vx = dvxCell / (sP.nVelocitySamples-1);
                  creal d_vy = dvyCell / (sP.nVelocitySamples-1);
                  creal d_vz = dvzCell / (sP.nVelocitySamples-1);
                  for (uint vi=0; vi<sP.nVelocitySamples; ++vi)
                     for (uint vj=0; vj<sP.nVelocitySamples; ++vj)
                        for (uint vk=0; vk<sP.nVelocitySamples; ++vk) {
                           average +=  shiftedMaxwellianDistribution(
                                                                     popID,
                                                                     vxCell + vi*d_vx,
                                                                     vyCell + vj*d_vy,
                                                                     vzCell + vk*d_vz
                                                                    );
                        }
                  average /= sP.nVelocitySamples * sP.nVelocitySamples * sP.nVelocitySamples;
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
      templateCell.parameters[CellParams::RHOM_R] = templateCell.parameters[CellParams::RHOM];
      templateCell.parameters[CellParams::VX_R] = templateCell.parameters[CellParams::VX];
      templateCell.parameters[CellParams::VY_R] = templateCell.parameters[CellParams::VY];
      templateCell.parameters[CellParams::VZ_R] = templateCell.parameters[CellParams::VZ];
      templateCell.parameters[CellParams::RHOQ_R] = templateCell.parameters[CellParams::RHOQ];
      templateCell.parameters[CellParams::P_11_R] = templateCell.parameters[CellParams::P_11];
      templateCell.parameters[CellParams::P_22_R] = templateCell.parameters[CellParams::P_22];
      templateCell.parameters[CellParams::P_33_R] = templateCell.parameters[CellParams::P_33];
      templateCell.parameters[CellParams::RHOM_V] = templateCell.parameters[CellParams::RHOM];
      templateCell.parameters[CellParams::VX_V] = templateCell.parameters[CellParams::VX];
      templateCell.parameters[CellParams::VY_V] = templateCell.parameters[CellParams::VY];
      templateCell.parameters[CellParams::VZ_V] = templateCell.parameters[CellParams::VZ];
      templateCell.parameters[CellParams::RHOQ_V] = templateCell.parameters[CellParams::RHOQ];
      templateCell.parameters[CellParams::P_11_V] = templateCell.parameters[CellParams::P_11];
      templateCell.parameters[CellParams::P_22_V] = templateCell.parameters[CellParams::P_22];
      templateCell.parameters[CellParams::P_33_V] = templateCell.parameters[CellParams::P_33];
   }
   
   Real Conductingsphere::shiftedMaxwellianDistribution(
      const uint popID,
      creal& vx, creal& vy, creal& vz
   ) {
      
      const Real MASS = getObjectWrapper().particleSpecies[popID].mass;
      const ConductingsphereSpeciesParameters& sP = this->speciesParams[popID];

      return sP.rho * pow(MASS /
      (2.0 * M_PI * physicalconstants::K_B * sP.T), 1.5) *
      exp(-MASS * ((vx-sP.V0[0])*(vx-sP.V0[0]) + (vy-sP.V0[1])*(vy-sP.V0[1]) + (vz-sP.V0[2])*(vz-sP.V0[2])) /
      (2.0 * physicalconstants::K_B * sP.T));
   }

   std::vector<vmesh::GlobalID> Conductingsphere::findBlocksToInitialize(spatial_cell::SpatialCell& cell,const uint popID) {
      vector<vmesh::GlobalID> blocksToInitialize;
      bool search = true;
      uint counter = 0;
      const uint8_t refLevel = 0;

      const vmesh::LocalID* vblocks_ini = cell.get_velocity_grid_length(popID,refLevel);

      while (search) {
         if (0.1 * cell.getVelocityBlockMinValue(popID) > shiftedMaxwellianDistribution(popID,counter*cell.get_velocity_grid_block_size(popID,refLevel)[0], 0.0, 0.0) || counter > vblocks_ini[0]) {
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

   void Conductingsphere::setCellFromTemplate(SpatialCell* cell,const uint popID) {
      copyCellData(&templateCell,cell,false,popID,true); // copy also vdf, _V
      copyCellData(&templateCell,cell,true,popID,false); // don't copy vdf again but copy _R now
   }

   std::string Conductingsphere::getName() const {return "Conductingsphere";}
   
   uint Conductingsphere::getIndex() const {return sysboundarytype::CONDUCTINGSPHERE;}
}
