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
 * \brief Implementation of the class SysBoundaryCondition::Copysphere to handle cells classified as
 * sysboundarytype::COPYSPHERE.
 */

#include <cstdlib>
#include <iostream>

#include "../common.h"
#include "../fieldsolver/fs_common.h"
#include "../fieldsolver/fs_limiters.h"
#include "../fieldsolver/ldz_magnetic_field.hpp"
#include "../fieldtracing/fieldtracing.h"
#include "../object_wrapper.h"
#include "../projects/project.h"
#include "../projects/projects_common.h"
#include "../vlasovsolver/vlasovmover.h"
#include "copysphere.h"

#ifdef DEBUG_VLASIATOR
#define DEBUG_COPYSPHERE
#endif
#ifdef DEBUG_SYSBOUNDARY
#define DEBUG_COPYSPHERE
#endif

namespace SBC {
   Copysphere::Copysphere() : SysBoundaryCondition() {}
   
   Copysphere::~Copysphere() {}
   
   void Copysphere::addParameters() {
      Readparameters::add("copysphere.centerX", "X coordinate of copysphere center (m)", 0.0);
      Readparameters::add("copysphere.centerY", "Y coordinate of copysphere center (m)", 0.0);
      Readparameters::add("copysphere.centerZ", "Z coordinate of copysphere center (m)", 0.0);
      Readparameters::add("copysphere.radius", "Radius of copysphere (m).", 1.0e7);
      Readparameters::add("copysphere.geometry",
                          "Select the geometry of the copysphere, 0: inf-norm (diamond), 1: 1-norm (square), 2: 2-norm "
                          "(circle, DEFAULT), 3: 2-norm cylinder aligned with y-axis, use with polar plane/line dipole.",
                          2);
      Readparameters::add(
          "copysphere.precedence",
          "Precedence value of the copysphere system boundary condition (integer), the higher the stronger.", 2);
      Readparameters::add("copysphere.reapplyUponRestart",
                          "If 0 (default), keep going with the state existing in the restart file. If 1, calls again "
                          "applyInitialState. Can be used to change boundary condition behaviour during a run.",
                          0);
      Readparameters::add("copysphere.zeroPerB",
                          "If 0 (default), normal copysphere behaviour of magnetic field at inner boundary. If 1, keep "
                          "magnetic field static at the inner boundary",
                          0);
   
      // Per-population parameters
      for (uint i = 0; i < getObjectWrapper().particleSpecies.size(); i++) {
         const std::string& pop = getObjectWrapper().particleSpecies[i].name;
   
         Readparameters::add(pop + "_copysphere.rho", "Number density of the copysphere (m^-3)", 0.0);
         Readparameters::add(pop + "_copysphere.T", "Temperature of the copysphere (K)", 0.0);
         Readparameters::add(pop + "_copysphere.VX0",
                             "Bulk velocity of copyspheric distribution function in X direction (m/s)", 0.0);
         Readparameters::add(pop + "_copysphere.VY0",
                             "Bulk velocity of copyspheric distribution function in X direction (m/s)", 0.0);
         Readparameters::add(pop + "_copysphere.VZ0",
                             "Bulk velocity of copyspheric distribution function in X direction (m/s)", 0.0);
         Readparameters::add(pop + "_copysphere.fluffiness",
                             "Inertia of boundary smoothing when copying neighbour's moments and velocity distributions "
                             "(0=completely constant boundaries, 1=neighbours are interpolated immediately).",
                             0);
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
      Readparameters::get("copysphere.reapplyUponRestart", reapply);
      this->applyUponRestart = false;
      if (reapply == 1) {
         this->applyUponRestart = true;
      }
      uint noperb;
      Readparameters::get("copysphere.zeroPerB", noperb);
      this->zeroPerB = false;
      if (noperb == 1) {
         this->zeroPerB = true;
      }
   
      for (uint i = 0; i < getObjectWrapper().particleSpecies.size(); i++) {
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
         if (sP.T == 0) {
            Readparameters::get(pop + "_Magnetosphere.T", sP.T);
         }
         if (sP.rho == 0) {
            Readparameters::get(pop + "_Magnetosphere.rho", sP.rho);
         }
   
         speciesParams.push_back(sP);
      }
   }
   
   void Copysphere::initSysBoundary(creal& t, Project& project) {
      getParameters();
      dynamic = false;
   
      // iniSysBoundary is only called once, generateTemplateCell must
      // init all particle species
      generateTemplateCell(project);
   }
   
   Real getR(creal x, creal y, creal z, uint geometry, Real center[3]) {
   
      Real r;
   
      switch (geometry) {
      case 0:
         // infinity-norm, result is a diamond/square with diagonals aligned on the axes in 2D
         r = fabs(x - center[0]) + fabs(y - center[1]) + fabs(z - center[2]);
         break;
      case 1:
         // 1-norm, result is is a grid-aligned square in 2D
         r = max(max(fabs(x - center[0]), fabs(y - center[1])), fabs(z - center[2]));
         break;
      case 2:
         // 2-norm (Cartesian), result is a circle in 2D
         r = sqrt((x - center[0]) * (x - center[0]) + (y - center[1]) * (y - center[1]) +
                  (z - center[2]) * (z - center[2]));
         break;
      case 3:
         // 2-norm (Cartesian) cylinder aligned on y-axis
         r = sqrt((x - center[0]) * (x - center[0]) + (z - center[2]) * (z - center[2]));
         break;
      default:
         abort_mpi("copysphere.geometry has to be 0, 1 or 2.", 1);
      }
   
      return r;
   }
   
   void Copysphere::assignSysBoundary(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid,
                                      std::span<fsgrids::technical> technical, FieldSolverGrid &fsgrid) {
      const vector<CellID>& cells = getLocalCells();
      for (uint i = 0; i < cells.size(); i++) {
         if (mpiGrid[cells[i]]->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
            continue;
         }
   
         creal* const cellParams = &(mpiGrid[cells[i]]->parameters[0]);
         creal dx = cellParams[CellParams::DX];
         creal dy = cellParams[CellParams::DY];
         creal dz = cellParams[CellParams::DZ];
         creal x = cellParams[CellParams::XCRD] + 0.5 * dx;
         creal y = cellParams[CellParams::YCRD] + 0.5 * dy;
         creal z = cellParams[CellParams::ZCRD] + 0.5 * dz;
   
         if (getR(x, y, z, this->geometry, this->center) < this->radius) {
            mpiGrid[cells[i]]->sysBoundaryFlag = this->getIndex();
         }
      }
   }
   
   void Copysphere::applyInitialState(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid,
                                      std::span<fsgrids::technical> technical, FieldSolverGrid &fsgrid,
                                      std::span<std::array<Real, fsgrids::bfield::N_BFIELD>> perb,
                                      fsgrids::bgbspan bgb, Project& project) {
      const vector<CellID>& cells = getLocalCells();
   #pragma omp parallel for
      for (uint i = 0; i < cells.size(); ++i) {
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
      std::span<fsgrids::technical> technical,
      FieldSolverGrid &fsgrid,
      cint i,
      cint j,
      cint k
   ) {
      std::array<Real, 3> normalDirection{{ 0.0, 0.0, 0.0 }};

      static creal DIAG2 = 1.0 / sqrt(2.0);
      static creal DIAG3 = 1.0 / sqrt(3.0);
   
      const auto& gridSpacing = fsgrid.getGridSpacing();
   
      creal dx = gridSpacing[0];
      creal dy = gridSpacing[1];
      creal dz = gridSpacing[2];
      const std::array<fsgrid::FsSize_t, 3> globalIndices = fsgrid.localToGlobal(i, j, k);
      creal x = P::xmin + (convert<Real>(globalIndices[0]) + 0.5) * dx;
      creal y = P::ymin + (convert<Real>(globalIndices[1]) + 0.5) * dy;
      creal z = P::zmin + (convert<Real>(globalIndices[2]) + 0.5) * dz;
      creal xsign = divideIfNonZero(x, fabs(x));
      creal ysign = divideIfNonZero(y, fabs(y));
      creal zsign = divideIfNonZero(z, fabs(z));
   
      Real length = 0.0;
   
      if (Parameters::xcells_ini == 1) {
         if (Parameters::ycells_ini == 1) {
            if (Parameters::zcells_ini == 1) {
               // X,Y,Z
               abort_mpi(
                   "What do you expect to do with a single-cell simulation of copysphere boundary type? Stop kidding.", 1);
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
            switch (this->geometry) {
            case 0:
               normalDirection[1] = DIAG2 * ysign;
               normalDirection[2] = DIAG2 * zsign;
               break;
            case 1:
               if (fabs(y) == fabs(z)) {
                  normalDirection[1] = ysign * DIAG2;
                  normalDirection[2] = zsign * DIAG2;
                  break;
               }
               if (fabs(y) > (this->radius - dy)) {
                  normalDirection[1] = ysign;
                  break;
               }
               if (fabs(z) > (this->radius - dz)) {
                  normalDirection[2] = zsign;
                  break;
               }
               if (fabs(y) > (this->radius - 2.0 * dy)) {
                  normalDirection[1] = ysign;
                  break;
               }
               if (fabs(z) > (this->radius - 2.0 * dz)) {
                  normalDirection[2] = zsign;
                  break;
               }
               break;
            case 2:
               length = sqrt(y * y + z * z);
               normalDirection[1] = y / length;
               normalDirection[2] = z / length;
               break;
            default:
               std::cerr << __FILE__ << ":" << __LINE__ << ":"
                         << "copysphere.geometry has to be 0, 1 or 2 with this grid shape." << std::endl;
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
            switch (this->geometry) {
            case 0:
               normalDirection[0] = DIAG2 * xsign;
               normalDirection[2] = DIAG2 * zsign;
               break;
            case 1:
               if (fabs(x) == fabs(z)) {
                  normalDirection[0] = xsign * DIAG2;
                  normalDirection[2] = zsign * DIAG2;
                  break;
               }
               if (fabs(x) > (this->radius - dx)) {
                  normalDirection[0] = xsign;
                  break;
               }
               if (fabs(z) > (this->radius - dz)) {
                  normalDirection[2] = zsign;
                  break;
               }
               if (fabs(x) > (this->radius - 2.0 * dx)) {
                  normalDirection[0] = xsign;
                  break;
               }
               if (fabs(z) > (this->radius - 2.0 * dz)) {
                  normalDirection[2] = zsign;
                  break;
               }
               break;
            case 2:
            case 3:
               length = sqrt(x * x + z * z);
               normalDirection[0] = x / length;
               normalDirection[2] = z / length;
               break;
            default:
               std::cerr << __FILE__ << ":" << __LINE__ << ":"
                         << "copysphere.geometry has to be 0, 1, 2 or 3 with this grid shape." << std::endl;
               abort();
            }
            // end of Y
         }
      } else if (Parameters::zcells_ini == 1) {
         // Z
         switch (this->geometry) {
         case 0:
            normalDirection[0] = DIAG2 * xsign;
            normalDirection[1] = DIAG2 * ysign;
            break;
         case 1:
            if (fabs(x) == fabs(y)) {
               normalDirection[0] = xsign * DIAG2;
               normalDirection[1] = ysign * DIAG2;
               break;
            }
            if (fabs(x) > (this->radius - dx)) {
               normalDirection[0] = xsign;
               break;
            }
            if (fabs(y) > (this->radius - dy)) {
               normalDirection[1] = ysign;
               break;
            }
            if (fabs(x) > (this->radius - 2.0 * dx)) {
               normalDirection[0] = xsign;
               break;
            }
            if (fabs(y) > (this->radius - 2.0 * dy)) {
               normalDirection[1] = ysign;
               break;
            }
            break;
         case 2:
            length = sqrt(x * x + y * y);
            normalDirection[0] = x / length;
            normalDirection[1] = y / length;
            break;
         default:
            abort_mpi("copysphere.geometry has to be 0, 1 or 2 with this grid shape.", 1);
         }
         // end of Z
      } else {
         // 3D
         switch (this->geometry) {
         case 0:
            normalDirection[0] = DIAG3 * xsign;
            normalDirection[1] = DIAG3 * ysign;
            normalDirection[2] = DIAG3 * zsign;
            break;
         case 1:
            if (fabs(x) == fabs(y) && fabs(x) == fabs(z) && fabs(x) > this->radius - dx) {
               normalDirection[0] = xsign * DIAG3;
               normalDirection[1] = ysign * DIAG3;
               normalDirection[2] = zsign * DIAG3;
               break;
            }
            if (fabs(x) == fabs(y) && fabs(x) == fabs(z) && fabs(x) > this->radius - 2.0 * dx) {
               normalDirection[0] = xsign * DIAG3;
               normalDirection[1] = ysign * DIAG3;
               normalDirection[2] = zsign * DIAG3;
               break;
            }
            if (fabs(x) == fabs(y) && fabs(x) > this->radius - dx && fabs(z) < this->radius - dz) {
               normalDirection[0] = xsign * DIAG2;
               normalDirection[1] = ysign * DIAG2;
               normalDirection[2] = 0.0;
               break;
            }
            if (fabs(y) == fabs(z) && fabs(y) > this->radius - dy && fabs(x) < this->radius - dx) {
               normalDirection[0] = 0.0;
               normalDirection[1] = ysign * DIAG2;
               normalDirection[2] = zsign * DIAG2;
               break;
            }
            if (fabs(x) == fabs(z) && fabs(x) > this->radius - dx && fabs(y) < this->radius - dy) {
               normalDirection[0] = xsign * DIAG2;
               normalDirection[1] = 0.0;
               normalDirection[2] = zsign * DIAG2;
               break;
            }
            if (fabs(x) == fabs(y) && fabs(x) > this->radius - 2.0 * dx && fabs(z) < this->radius - 2.0 * dz) {
               normalDirection[0] = xsign * DIAG2;
               normalDirection[1] = ysign * DIAG2;
               normalDirection[2] = 0.0;
               break;
            }
            if (fabs(y) == fabs(z) && fabs(y) > this->radius - 2.0 * dy && fabs(x) < this->radius - 2.0 * dx) {
               normalDirection[0] = 0.0;
               normalDirection[1] = ysign * DIAG2;
               normalDirection[2] = zsign * DIAG2;
               break;
            }
            if (fabs(x) == fabs(z) && fabs(x) > this->radius - 2.0 * dx && fabs(y) < this->radius - 2.0 * dy) {
               normalDirection[0] = xsign * DIAG2;
               normalDirection[1] = 0.0;
               normalDirection[2] = zsign * DIAG2;
               break;
            }
            if (fabs(x) > (this->radius - dx)) {
               normalDirection[0] = xsign;
               break;
            }
            if (fabs(y) > (this->radius - dy)) {
               normalDirection[1] = ysign;
               break;
            }
            if (fabs(z) > (this->radius - dz)) {
               normalDirection[2] = zsign;
               break;
            }
            if (fabs(x) > (this->radius - 2.0 * dx)) {
               normalDirection[0] = xsign;
               break;
            }
            if (fabs(y) > (this->radius - 2.0 * dy)) {
               normalDirection[1] = ysign;
               break;
            }
            if (fabs(z) > (this->radius - 2.0 * dz)) {
               normalDirection[2] = zsign;
               break;
            }
            break;
         case 2:
            length = sqrt(x * x + y * y + z * z);
            normalDirection[0] = x / length;
            normalDirection[1] = y / length;
            normalDirection[2] = z / length;
            break;
         case 3:
            length = sqrt(x * x + z * z);
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
   Real Copysphere::fieldSolverBoundaryCondMagneticField(std::span<const std::array<Real, fsgrids::bfield::N_BFIELD>> b,
                                                         std::span<const std::array<Real, fsgrids::bgbfield::N_BGB>> bgb,
                                                         std::span<const fsgrids::technical> technical,
                                                         const std::array<Real, 3>& gridSpacing,
                                                         const std::array<fsgrid::FsSize_t, 3>& globalCoordinates,
                                                         const fsgrid::FsStencil& stencil, cuint component) {
      const uint32_t perbComponent = fsgrids::bfield::PERBX + component;
      const uint32_t bitfield = 1 << component;
   
      // clang-format off
      static constexpr std::array permutations = {
          std::array {
              0, 1, 2, 3, 4, 5,
          },
          std::array {
              2, 3, 0, 1, 4, 5,
          },
          std::array {
              4, 5, 0, 1, 2, 3,
          },
      };
   
      const std::array permutation = permutations[component];
   
      const std::array<size_t, 6> inds = {
          stencil.moo(),
          stencil.poo(),
          stencil.omo(),
          stencil.opo(),
          stencil.oom(),
          stencil.oop(),
      };
      // clang-format on
   
      auto bitFieldSet = [&bitfield](auto& tech) { return (tech.SOLVE & bitfield) == bitfield; };
      auto sbLayerIsOne = [](auto& tech) { return tech.sysBoundaryLayer == 1; };
      auto averageNeigbours = [&technical, &b, &inds, &permutation, &perbComponent,
                               &bitFieldSet](auto begin, auto end, auto& sum, auto& nCells) {
         for (size_t i = begin; i < end; i++) {
            const auto j = inds[permutation[i]];
            if (bitFieldSet(technical[j])) {
               sum += b[j][perbComponent];
               nCells++;
            }
         }
      };
   
      auto averageAllNeighbours = [&stencil, &technical, &b, &perbComponent](auto predicateLambda, auto& sum,
                                                                             auto& nCells) {
         for (const auto& i : stencil.indices()) {
            if (predicateLambda(technical[i])) {
               sum += b[i][perbComponent];
               nCells++;
            }
         }
      };
   
      Real sum = 0.0;
      uint nCells = 0;
      if (this->zeroPerB == true) {
         sum = b[stencil.ooo()][perbComponent];
         nCells = 1;
      } else {
         if (sbLayerIsOne(technical[stencil.ooo()])) {
            averageNeigbours(0ul, 2ul, sum, nCells);
   
            if (nCells == 0) {
               averageNeigbours(2ul, 6ul, sum, nCells);
            }
   
            if (nCells == 0) {
               averageAllNeighbours(bitFieldSet, sum, nCells);
            }
         } else {
            // L2 cells
            averageAllNeighbours(sbLayerIsOne, sum, nCells);
         }
      }
   
      if (nCells == 0) {
         cerr << __FILE__ << ":" << __LINE__ << ": ERROR: this should not have fallen through." << endl;
         sum = 0.0;
         nCells = 1;
      }
   
      return sum / nCells;
   }
   
   void Copysphere::fieldSolverBoundaryCondElectricField(std::span<std::array<Real, fsgrids::efield::N_EFIELD>> e,
                                                         const fsgrid::FsStencil& stencil, cuint component) {
      e[stencil.ooo()][fsgrids::efield::EX + component] = 0.0;
   }
   
   void Copysphere::fieldSolverBoundaryCondHallElectricField(std::span<std::array<Real, fsgrids::ehall::N_EHALL>> ehall,
                                                             const fsgrid::FsStencil& stencil, cuint component) {
      std::array<Real, fsgrids::ehall::N_EHALL>& cp = ehall[stencil.ooo()];
      switch (component) {
      case 0:
         cp[fsgrids::ehall::EXHALL_000_100] = 0.0;
         cp[fsgrids::ehall::EXHALL_010_110] = 0.0;
         cp[fsgrids::ehall::EXHALL_001_101] = 0.0;
         cp[fsgrids::ehall::EXHALL_011_111] = 0.0;
         break;
      case 1:
         cp[fsgrids::ehall::EYHALL_000_010] = 0.0;
         cp[fsgrids::ehall::EYHALL_100_110] = 0.0;
         cp[fsgrids::ehall::EYHALL_001_011] = 0.0;
         cp[fsgrids::ehall::EYHALL_101_111] = 0.0;
         break;
      case 2:
         cp[fsgrids::ehall::EZHALL_000_001] = 0.0;
         cp[fsgrids::ehall::EZHALL_100_101] = 0.0;
         cp[fsgrids::ehall::EZHALL_010_011] = 0.0;
         cp[fsgrids::ehall::EZHALL_110_111] = 0.0;
         break;
      default:
         cerr << __FILE__ << ":" << __LINE__ << ":"
              << " Invalid component" << endl;
      }
   }
   
   void Copysphere::fieldSolverBoundaryCondGradPeElectricField(
       std::span<std::array<Real, fsgrids::egradpe::N_EGRADPE>> EGradPe, const fsgrid::FsStencil& stencil,
       cuint component) {
      EGradPe[stencil.ooo()][fsgrids::egradpe::EXGRADPE + component] = 0.0;
   }
   
   void Copysphere::fieldSolverBoundaryCondDerivatives(std::span<std::array<Real, fsgrids::dperb::N_DPERB>> dperb,
                                                       fsgrids::dmomentsspan dmoments,
                                                       const fsgrid::FsStencil& stencil, cuint RKCase, cuint component) {
      this->setCellDerivativesToZero(dperb, dmoments, stencil, component);
   }
   
   void Copysphere::fieldSolverBoundaryCondBVOLDerivatives(std::span<std::array<Real, fsgrids::volfields::N_VOL>> vols,
                                                           const fsgrid::FsStencil& stencil, cuint component) {
      // FIXME This should be OK as the BVOL derivatives are only used for Lorentz force JXB, which is not applied on the
      // copy sphere cells.
      this->setCellBVOLDerivativesToZero(vols, stencil, component);
   }
   
   void Copysphere::vlasovBoundaryCondition(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid,
                                            const CellID& cellID, const uint popID, const bool calculate_V_moments) {
      this->vlasovBoundaryFluffyCopyFromAllCloseNbrs(mpiGrid, cellID, popID, calculate_V_moments,
                                                     this->speciesParams[popID].fluffiness);
   }
   
   /**
    * NOTE: This function must initialize all particle species!
    * @param project
    */
   void Copysphere::generateTemplateCell(Project& project) {
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
      for (uint popID = 0; popID < getObjectWrapper().particleSpecies.size(); ++popID) {
         templateCell.clear(popID, false); // clear, do not de-allocate memory
         const CopysphereSpeciesParameters& sP = this->speciesParams[popID];
         const Real mass = getObjectWrapper().particleSpecies[popID].mass;
         initRho = sP.rho;
         initT = sP.T;
         initV0X = sP.V0[0];
         initV0Y = sP.V0[1];
         initV0Z = sP.V0[2];

         // Find list of blocks to initialize.
         const uint nRequested =
             SBC::findMaxwellianBlocksToInitialize(popID, templateCell, initRho, initT, initV0X, initV0Y, initV0Z);
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
         vmesh::VelocityMesh* vmesh = templateCell.get_velocity_mesh(popID);
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
   
      calculateCellMoments(&templateCell, true, false, true);
   
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
   
   void Copysphere::setCellFromTemplate(SpatialCell* cell, const uint popID) {
      copyCellData(&templateCell, cell, false, popID, true); // copy also vdf, _V
      copyCellData(&templateCell, cell, true, popID, false); // don't copy vdf again but copy _R now
   #ifdef USE_GPU
      cell->setReservation(popID, templateCell.getReservation(popID));
   #endif
   }
   
   std::string Copysphere::getName() const { return "Copysphere"; }
   void Copysphere::getFaces(bool* faces) {}
   
   void Copysphere::updateState(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid,
                                std::span<fsgrids::technical> technical, FieldSolverGrid &fsgrid,
                                std::span<std::array<Real, fsgrids::bfield::N_BFIELD>> perb,
                                fsgrids::bgbspan bgb, creal t) {}
   
   uint Copysphere::getIndex() const { return sysboundarytype::COPYSPHERE; }
} // namespace SBC
