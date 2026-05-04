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

#ifndef COPYSPHERE_H
#define COPYSPHERE_H

#include <vector>
#include "../definitions.h"
#include "../readparameters.h"
#include "../spatial_cells/spatial_cell_wrapper.hpp"
#include "sysboundarycondition.h"

using namespace projects;
using namespace std;

namespace SBC {

   struct CopysphereSpeciesParameters {
      Real rho;
      Real V0[3];
      Real T;
      Real fluffiness;
   };

   /*!\brief Copysphere is a class applying an ionosphere-ish boundary conditions.
    * 
    * Copysphere is a class handling cells tagged as sysboundarytype::COPYSPHERE by this system boundary condition. It applies copy boundary conditions to perturbed magnetic field.
    * 
    * These consist in:
    * - Do nothing for the distribution (keep the initial state constant in time);
    * - Copy the closest neighbors' perturbed B and average it;
    * - Null out the electric fields.
    *
    * For 3D magnetospheric simulations, you might be interesting in trying the ionosphere boundary instead!
    */
   class Copysphere: public SysBoundaryCondition {
   public:
      Copysphere();
      virtual ~Copysphere();
      
      static void addParameters();
      virtual void getParameters() override;
      
      virtual void initSysBoundary(
         creal& t,
         Project &project
      ) override;
      virtual void assignSysBoundary(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                     FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid) override;
      virtual void applyInitialState(
         dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
         FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
         FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
         FsGrid<std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH>& BgBGrid,
         Project &project
      ) override;
      virtual void updateState(
         dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid,
         FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
         FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> &perBGrid,
         FsGrid<std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> &BgBGrid,
         creal t
      ) override;
      virtual Real fieldSolverBoundaryCondMagneticField(
         FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & bGrid,
         FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & bgbGrid,
         FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
         cint i,
         cint j,
         cint k,
         creal dt,
         cuint component
      ) override;
      virtual void fieldSolverBoundaryCondElectricField(
         FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
         cint i,
         cint j,
         cint k,
         cuint component
      ) override;
      virtual void fieldSolverBoundaryCondHallElectricField(
         FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
         cint i,
         cint j,
         cint k,
         cuint component
      ) override;
      virtual void fieldSolverBoundaryCondGradPeElectricField(
         FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
         cint i,
         cint j,
         cint k,
         cuint component
      ) override;
      virtual void fieldSolverBoundaryCondDerivatives(
         FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
         FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
         cint i,
         cint j,
         cint k,
         cuint RKCase,
         cuint component
      ) override;
      virtual void fieldSolverBoundaryCondBVOLDerivatives(
         FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
         cint i,
         cint j,
         cint k,
         cuint component
      ) override;
      virtual void vlasovBoundaryCondition(
         dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
         const CellID& cellID,
         const uint popID,
         const bool calculate_V_moments
      ) override;
      
      void getFaces(bool *faces) override;
      virtual std::string getName() const override;
      virtual uint getIndex() const override;
      
      void generateTemplateCell(Project &project);
      void setCellFromTemplate(SpatialCell* cell,const uint popID);
      
      std::array<Real, 3> fieldSolverGetNormalDirection(
         FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
         cint i,
         cint j,
         cint k
      );
      
      Real center[3]; /*!< Coordinates of the centre of the copy sphere. */
      Real radius; /*!< Radius of the copy sphere. */
      uint geometry; /*!< Geometry of the copy sphere, 0: inf-norm (diamond), 1: 1-norm (square), 2: 2-norm (circle, DEFAULT), 3: polar-plane cylinder with line dipole. */

      std::vector<CopysphereSpeciesParameters> speciesParams;
      bool zeroPerB;
      
      spatial_cell::SpatialCell templateCell;

   };
}

#endif
