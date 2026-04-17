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

      virtual void initSysBoundary(creal& t, Project& project) override;
      virtual void assignSysBoundary(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid,
                                     fsgrids::technicalspan technical, FieldSolverGrid& fsgrid) override;
      virtual void applyInitialState(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid,
                                     fsgrids::technicalspan technical, FieldSolverGrid& fsgrid,
                                     fsgrids::perbspan perb,
                                     fsgrids::bgbspan bgb,
                                     Project& project) override;
      virtual void updateState(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid,
                               fsgrids::technicalspan technical, FieldSolverGrid& fsgrid,
                               fsgrids::perbspan perb,
                               fsgrids::bgbspan bgb, creal t) override;
      virtual Real fieldSolverBoundaryCondMagneticField(fsgrids::perbspan b,
                                                        fsgrids::constbgbspan bgb,
                                                        fsgrids::consttechnicalspan technical,
                                                        const std::array<Real, 3>& gridSpacing,
                                                        const std::array<fsgrid::FsSize_t, 3>& globalCoordinates,
                                                        const fsgrid::FsStencil& stencil, cuint component) override;
      virtual void fieldSolverBoundaryCondElectricField(fsgrids::efieldspan e,
                                                        const fsgrid::FsStencil& stencil, cuint component) override;
      virtual void fieldSolverBoundaryCondHallElectricField(fsgrids::ehallspan ehall,
                                                            const fsgrid::FsStencil& stencil, cuint component) override;
      virtual void
      fieldSolverBoundaryCondGradPeElectricField(fsgrids::egradpespan EGradPe,
                                                 const fsgrid::FsStencil& stencil, cuint component) override;
      virtual void
      fieldSolverBoundaryCondDerivatives(fsgrids::dperbspan dperb,
                                         fsgrids::dmomentsspan dmoments,
                                         const fsgrid::FsStencil& stencil, cuint RKCase, cuint component) override;
      virtual void fieldSolverBoundaryCondBVOLDerivatives(fsgrids::volspan vols,
                                                          const fsgrid::FsStencil& stencil, cuint component) override;
      virtual void vlasovBoundaryCondition(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid,
                                           const CellID& cellID, const uint popID,
                                           const bool calculate_V_moments) override;

      void getFaces(bool* faces) override;
      virtual std::string getName() const override;
      virtual uint getIndex() const override;

      void generateTemplateCell(Project &project);
      void setCellFromTemplate(SpatialCell* cell,const uint popID);

      std::array<Real, 3> fieldSolverGetNormalDirection(
         fsgrids::technicalspan technical, FieldSolverGrid &fsgrid,
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
