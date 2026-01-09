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

/*!\file copysphere.cpp
 * \brief Implementation of the class SysBoundaryCondition::Copysphere to handle cells classified as sysboundarytype::COPYSPHERE.
 */

#include <cstdlib>
#include <iostream>

#include "copysphere.h"
#include "../projects/project.h"
#include "../projects/projects_common.h"
#include "../vlasovsolver/vlasovmover.h"
#include "../fieldsolver/fs_common.h"
#include "../fieldsolver/fs_limiters.h"
#include "../fieldsolver/ldz_magnetic_field.hpp"
#include "../fieldtracing/fieldtracing.h"
#include "../common.h"
#include "../object_wrapper.h"


#ifdef DEBUG_VLASIATOR
   #define DEBUG_COPYSPHERE
#endif
#ifdef DEBUG_SYSBOUNDARY
   #define DEBUG_COPYSPHERE
#endif

namespace SBC {
   Copysphere::Copysphere(): SysBoundaryCondition() { }

   Copysphere::~Copysphere() { }

   void Copysphere::addParameters() {
      Readparameters::add("copysphere.centerX", "X coordinate of copysphere center (m)", 0.0);
      Readparameters::add("copysphere.centerY", "Y coordinate of copysphere center (m)", 0.0);
      Readparameters::add("copysphere.centerZ", "Z coordinate of copysphere center (m)", 0.0);
      Readparameters::add("copysphere.radius", "Radius of copysphere (m).", 1.0e7);
      Readparameters::add("copysphere.geometry", "Select the geometry of the copysphere, 0: inf-norm (diamond), 1: 1-norm (square), 2: 2-norm (circle, DEFAULT), 3: 2-norm cylinder aligned with y-axis, use with polar plane/line dipole.", 2);
      Readparameters::add("copysphere.precedence", "Precedence value of the copysphere system boundary condition (integer), the higher the stronger.", 2);
      Readparameters::add("copysphere.reapplyUponRestart", "If 0 (default), keep going with the state existing in the restart file. If 1, calls again applyInitialState. Can be used to change boundary condition behaviour during a run.", 0);
      Readparameters::add("copysphere.zeroPerB","If 0 (default), normal copysphere behaviour of magnetic field at inner boundary. If 1, keep magnetic field static at the inner boundary",0);

      // Per-population parameters
      for(uint i=0; i< getObjectWrapper().particleSpecies.size(); i++) {
         const std::string& pop = getObjectWrapper().particleSpecies[i].name;

         Readparameters::add(pop + "_copysphere.rho", "Number density of the copysphere (m^-3)", 0.0);
         Readparameters::add(pop + "_copysphere.T", "Temperature of the copysphere (K)", 0.0);
         Readparameters::add(pop + "_copysphere.VX0", "Bulk velocity of copyspheric distribution function in X direction (m/s)", 0.0);
         Readparameters::add(pop + "_copysphere.VY0", "Bulk velocity of copyspheric distribution function in X direction (m/s)", 0.0);
         Readparameters::add(pop + "_copysphere.VZ0", "Bulk velocity of copyspheric distribution function in X direction (m/s)", 0.0);
         Readparameters::add(pop + "_copysphere.fluffiness", "Inertia of boundary smoothing when copying neighbour's moments and velocity distributions (0=completely constant boundaries, 1=neighbours are interpolated immediately).", 0);
      }
   }

   void Copysphere::getParameters() {

      Readparameters::get("copysphere.centerX", this->center[0]);
      Readparameters::get("copysphere.centerY", this->center[1]);
      Readparameters::get("copysphere.centerZ", this->center[2]);
      Readparameters::get("copysphere.radius", this->radius);
      FieldTracing::fieldTracingParameters.innerBoundaryRadius = this->radius;
      Readparameters::get("copysphere.geometry", this->geometry);
      Readparameters::get("copysphere.precedence", this->precedence);
      uint reapply;
      Readparameters::get("copysphere.reapplyUponRestart",reapply);
      this->applyUponRestart = false;
      if(reapply == 1) {
         this->applyUponRestart = true;
      }
      uint noperb;
      Readparameters::get("copysphere.zeroPerB",noperb);
      this->zeroPerB = false;
      if(noperb == 1) {
         this->zeroPerB = true;
      }

      for(uint i=0; i< getObjectWrapper().particleSpecies.size(); i++) {
        const std::string& pop = getObjectWrapper().particleSpecies[i].name;
        CopysphereSpeciesParameters sP;

        Readparameters::get(pop + "_copysphere.rho", sP.rho);
        Readparameters::get(pop + "_copysphere.VX0", sP.V0[0]);
        Readparameters::get(pop + "_copysphere.VY0", sP.V0[1]);
        Readparameters::get(pop + "_copysphere.VZ0", sP.V0[2]);
        Readparameters::get(pop + "_copysphere.fluffiness", sP.fluffiness);
        Readparameters::get(pop + "_copysphere.T", sP.T);

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

   void Copysphere::initSysBoundary(
      creal& t,
      Project &project
   ) {
      getParameters();
      dynamic = false;

      // iniSysBoundary is only called once, generateTemplateCell must
      // init all particle species
      generateTemplateCell(project);
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
         abort_mpi("copysphere.geometry has to be 0, 1 or 2.", 1);
      }

      return r;
   }

   void Copysphere::assignSysBoundary(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
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

   }

   void Copysphere::applyInitialState(
      dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
      FsGrid<std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH>& BgBGrid,
      Project &project
   ) {
      const vector<CellID>& cells = getLocalCells();
      #pragma omp parallel for
      for (uint i=0; i<cells.size(); ++i) {
         SpatialCell* cell = mpiGrid[cells[i]];
         if (cell->sysBoundaryFlag != this->getIndex()) continue;
         for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
            setCellFromTemplate(cell,popID);
            #ifdef DEBUG_VLASIATOR
            // Verify current mesh and blocks
            if (!cell->checkMesh(popID)) {
               printf("ERROR in vmesh check: %s at %d\n",__FILE__,__LINE__);
            }
            #endif
         }
      }
   }

   std::array<Real, 3> Copysphere::fieldSolverGetNormalDirection(
      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
      cint i,
      cint j,
      cint k
   ) {
      std::array<Real, 3> normalDirection{{ 0.0, 0.0, 0.0 }};

      static creal DIAG2 = 1.0 / sqrt(2.0);
      static creal DIAG3 = 1.0 / sqrt(3.0);

      creal dx = technicalGrid.DX;
      creal dy = technicalGrid.DY;
      creal dz = technicalGrid.DZ;
      const std::array<FsGridTools::FsIndex_t, 3> globalIndices = technicalGrid.getGlobalIndices(i,j,k);
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
               abort_mpi("What do you expect to do with a single-cell simulation of copysphere boundary type? Stop kidding.", 1);
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
                  std::cerr << __FILE__ << ":" << __LINE__ << ":" << "copysphere.geometry has to be 0, 1 or 2 with this grid shape." << std::endl;
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
                  std::cerr << __FILE__ << ":" << __LINE__ << ":" << "copysphere.geometry has to be 0, 1, 2 or 3 with this grid shape." << std::endl;
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
               abort_mpi("copysphere.geometry has to be 0, 1 or 2 with this grid shape.", 1);
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
               abort_mpi("copysphere.geometry has to be 0, 1, 2 or 3 with this grid shape.", 1);
         }
         // end of 3D
      }
      return normalDirection;
   }

   /*! We want here to
    *
    * -- Average perturbed face B from the nearest neighbours
    *
    * -- Retain only the normal components of perturbed face B
    */
   Real Copysphere::fieldSolverBoundaryCondMagneticField(
      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & bGrid,
      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & bgbGrid,
      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
      cint i,
      cint j,
      cint k,
      creal dt,
      cuint component
   ) {
      if(this->zeroPerB == true){
         return bGrid.get(i,j,k)->at(fsgrids::bfield::PERBX+component);
      } else {
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
                  cerr << "ERROR: copysphere boundary tried to copy nonsensical magnetic field component " << component << endl;
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

   }

   void Copysphere::fieldSolverBoundaryCondElectricField(
      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
      cint i,
      cint j,
      cint k,
      cuint component
   ) {
      EGrid.get(i,j,k)->at(fsgrids::efield::EX+component) = 0.0;
   }

   void Copysphere::fieldSolverBoundaryCondHallElectricField(
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

   void Copysphere::fieldSolverBoundaryCondGradPeElectricField(
      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
      cint i,
      cint j,
      cint k,
      cuint component
   ) {
      EGradPeGrid.get(i,j,k)->at(fsgrids::egradpe::EXGRADPE+component) = 0.0;
   }

   void Copysphere::fieldSolverBoundaryCondDerivatives(
      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
      cint i,
      cint j,
      cint k,
      cuint RKCase,
      cuint component
   ) {
      this->setCellDerivativesToZero(dPerBGrid, dMomentsGrid, i, j, k, component);
      return;
   }

   void Copysphere::fieldSolverBoundaryCondBVOLDerivatives(
      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
      cint i,
      cint j,
      cint k,
      cuint component
   ) {
      // FIXME This should be OK as the BVOL derivatives are only used for Lorentz force JXB, which is not applied on the copy sphere cells.
      this->setCellBVOLDerivativesToZero(volGrid, i, j, k, component);
   }

   void Copysphere::vlasovBoundaryCondition(
      dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const CellID& cellID,
      const uint popID,
      const bool calculate_V_moments
   ) {
      this->vlasovBoundaryFluffyCopyFromAllCloseNbrs(mpiGrid, cellID, popID, calculate_V_moments, this->speciesParams[popID].fluffiness);
   }

   /**
    * NOTE: This function must initialize all particle species!
    * @param project
    */
   void Copysphere::generateTemplateCell(Project &project) {
      // WARNING not 0.0 here or the dipole() function fails miserably.
      templateCell.sysBoundaryFlag = this->getIndex();
      templateCell.sysBoundaryLayer = 1;
      templateCell.parameters[CellParams::XCRD] = 1.0;
      templateCell.parameters[CellParams::YCRD] = 1.0;
      templateCell.parameters[CellParams::ZCRD] = 1.0;
      templateCell.parameters[CellParams::DX] = 1;
      templateCell.parameters[CellParams::DY] = 1;
      templateCell.parameters[CellParams::DZ] = 1;

      Real initRho, initT, initV0X, initV0Y, initV0Z;
      // Loop over particle species
      for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
         templateCell.clear(popID,false); //clear, do not de-allocate memory
         const CopysphereSpeciesParameters& sP = this->speciesParams[popID];
         const Real mass = getObjectWrapper().particleSpecies[popID].mass;
         initRho = sP.rho;
         initT = sP.T;
         initV0X = sP.V0[0];
         initV0Y = sP.V0[1];
         initV0Z = sP.V0[2];

         // Find list of blocks to initialize.
         const uint nRequested = SBC::findMaxwellianBlocksToInitialize(popID,templateCell, initRho, initT, initV0X, initV0Y, initV0Z);
         // stores in vmesh->getGrid() (localToGlobalMap)
         // with count in cell.get_population(popID).N_blocks

         // Resize and populate mesh
         templateCell.prepare_to_receive_blocks(popID);

         // Set the reservation value (capacity is increased in add_velocity_blocks
         const Realf minValue = templateCell.getVelocityBlockMinValue(popID);

         // fills v-space into target

         #ifdef USE_GPU
         vmesh::VelocityMesh *vmesh = templateCell.dev_get_velocity_mesh(popID);
         vmesh::VelocityBlockContainer* VBC = templateCell.dev_get_velocity_blocks(popID);
         #else
         vmesh::VelocityMesh *vmesh = templateCell.get_velocity_mesh(popID);
         vmesh::VelocityBlockContainer* VBC = templateCell.get_velocity_blocks(popID);
         #endif
         // Loop over blocks
         Realf rhosum = 0;
         arch::parallel_reduce<arch::null>(
            {WID, WID, WID, nRequested},
            ARCH_LOOP_LAMBDA (const uint i, const uint j, const uint k, const uint initIndex, Realf *lsum ) {
               vmesh::GlobalID *GIDlist = vmesh->getGrid()->data();
               Realf* bufferData = VBC->getData();
               const vmesh::GlobalID blockGID = GIDlist[initIndex];
               // Calculate parameters for new block
               Real blockCoords[6];
               vmesh->getBlockInfo(blockGID,&blockCoords[0]);
               creal vxBlock = blockCoords[0];
               creal vyBlock = blockCoords[1];
               creal vzBlock = blockCoords[2];
               creal dvxCell = blockCoords[3];
               creal dvyCell = blockCoords[4];
               creal dvzCell = blockCoords[5];
               ARCH_INNER_BODY(i, j, k, initIndex, lsum) {
                  creal vx = vxBlock + (i+0.5)*dvxCell - initV0X;
                  creal vy = vyBlock + (j+0.5)*dvyCell - initV0Y;
                  creal vz = vzBlock + (k+0.5)*dvzCell - initV0Z;
                  const Realf value = projects::MaxwellianPhaseSpaceDensity(vx,vy,vz,initT,initRho,mass);
                  bufferData[initIndex*WID3 + k*WID2 + j*WID + i] = value;
                  //lsum[0] += value;
               };
            }, rhosum);

         #ifdef USE_GPU
         // Set and apply the reservation value
         templateCell.setReservation(popID,nRequested,true); // Force to this value
         templateCell.applyReservation(popID);
         #endif

         //let's get rid of blocks not fulfilling the criteria here to save memory.
         templateCell.adjustSingleCellVelocityBlocks(popID,true);

      } // for-loop over particle species

      calculateCellMoments(&templateCell,true,false,true);

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

   void Copysphere::setCellFromTemplate(SpatialCell* cell,const uint popID) {
      copyCellData(&templateCell,cell,false,popID,true); // copy also vdf, _V
      copyCellData(&templateCell,cell,true,popID,false); // don't copy vdf again but copy _R now
      #ifdef USE_GPU
      cell->setReservation(popID,templateCell.getReservation(popID));
      #endif
   }

   std::string Copysphere::getName() const {return "Copysphere";}
   void Copysphere::getFaces(bool *faces) {}

   void Copysphere::updateState(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid,
                                FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
                                FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH>& perBGrid,
                                FsGrid<std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH>& BgBGrid,
                                creal t) {}

   uint Copysphere::getIndex() const {return sysboundarytype::COPYSPHERE;}
}
