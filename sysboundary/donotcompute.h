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

#ifndef DONOTCOMPUTE_H
#define DONOTCOMPUTE_H

#include <vector>
#include "../definitions.h"
#include "../readparameters.h"
#include "../spatial_cells/spatial_cell_wrapper.hpp"
#include "sysboundarycondition.h"

using namespace projects;

namespace SBC {
   /*!\brief DoNotCompute is a class handling cells not to be computed.
    * 
    * DoNotCompute is a class handling cells tagged as sysboundarytype::DO_NOT_COMPUTE by a system boundary condition (e.g. SysBoundaryCondition::Ionosphere).
    */
   class DoNotCompute: public SysBoundaryCondition {
   public:
      DoNotCompute();
      virtual ~DoNotCompute();
      
      static void addParameters();
      virtual void getParameters() override;
      
      void initSysBoundary(
         creal& t,
         Project &project
      ) override ;
      void assignSysBoundary(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                     FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid) override;
      void applyInitialState(
         dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
         FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
         FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
         FsGrid<std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH>& BgBGrid,
         Project &project
      ) override;
      void updateState(
         dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid,
         FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
         FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> &perBGrid,
         FsGrid<std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH>& BgBGrid,
         creal t
      ) override;
      void getFaces(bool *faces) override;
      std::string getName() const override;
      uint getIndex() const override;

      // Explicit warning functions to inform the user if a doNotCompute cell gets computed
      Real fieldSolverBoundaryCondMagneticField(
         FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
         FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & bgbGrid,
         FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
         cint i,
         cint j,
         cint k,
         creal dt,
         cuint component
      ) override { std::cerr << "ERROR: DoNotCompute::fieldSolverBoundaryCondMagneticField called!" << std::endl; return 0.;}
      void fieldSolverBoundaryCondElectricField(
         FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
         cint i,
         cint j,
         cint k,
         cuint component
      ) override { std::cerr << "ERROR: DoNotCompute::fieldSolverBoundaryCondElectricField called!" << std::endl;}
      void fieldSolverBoundaryCondHallElectricField(
         FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
         cint i,
         cint j,
         cint k,
         cuint component
      ) override { std::cerr << "ERROR: DoNotCompute::fieldSolverBoundaryCondHallElectricField called!" << std::endl;}
      void fieldSolverBoundaryCondGradPeElectricField(
         FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
         cint i,
         cint j,
         cint k,
         cuint component
      ) override { std::cerr << "ERROR: DoNotCompute::fieldSolverBoundaryCondGradPeElectricField called!" << std::endl;}
      void fieldSolverBoundaryCondDerivatives(
         FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
         FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
         cint i,
         cint j,
         cint k,
         cuint RKCase,
         cuint component
      ) override{ std::cerr << "ERROR: DoNotCompute::fieldSolverBoundaryCondDerivatives called!" << std::endl;}
      void fieldSolverBoundaryCondBVOLDerivatives(
         FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
         cint i,
         cint j,
         cint k,
         cuint component
      ) override { std::cerr << "ERROR: DoNotCompute::fieldSolverBoundaryCondBVOLDerivatives called!" << std::endl;}
      void vlasovBoundaryCondition(
          dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
          const CellID& cellID,
          const uint popID,
          const bool calculate_V_moments
      ) override { std::cerr << "ERROR: DoNotCompute::vlasovBoundaryCondition called!" << std::endl;}
   };
}

#endif
