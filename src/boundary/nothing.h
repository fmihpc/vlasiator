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

#ifndef NOTHING_H
#define NOTHING_H

#include "../definitions.h"
#include "../readparameters.h"
#include "../spatial_cell.hpp"
#include "boundarycondition.h"
#include <vector>

using namespace projects;

namespace BC {
/*!\brief Nothing is a class handling cells not to be computed.*/
class Nothing : public BoundaryCondition {
public:
   Nothing();
   ~Nothing() override;

   static void addParameters();
   void getParameters() override;

   void initBoundary(creal t, Project &project) override;
   void assignBoundary(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid,
                       FsGrid<fsgrids::technical, FS_STENCIL_WIDTH> &technicalGrid) override;
   void applyInitialState(const dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid,
                          FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> &perBGrid,
                          Project &project) override;
   void updateState(const dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid,
                    FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> &perBGrid, creal t) override;
   std::string getName() const override;
   uint getIndex() const override;

   // Explicit warning functions to inform the user if a Nothing cell gets computed
   Real
   fieldSolverBoundaryCondMagneticField(FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> &perBGrid,
                                        FsGrid<fsgrids::technical, FS_STENCIL_WIDTH> &technicalGrid, cint i, cint j,
                                        cint k, creal dt, cuint component) override {
      std::string errmsg = "ERROR: calling Nothing::";
      errmsg += __func__;
      abort_mpi(errmsg);
      return 0.;
   }
   void
   fieldSolverBoundaryCondElectricField(FsGrid<std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> &EGrid,
                                        cint i, cint j, cint k, cuint component) override {
      std::string errmsg = "ERROR: calling Nothing::";
      errmsg += __func__;
      abort_mpi(errmsg);
   }
   void fieldSolverBoundaryCondHallElectricField(
       FsGrid<std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> &EHallGrid, cint i, cint j, cint k,
       cuint component) override {
      std::string errmsg = "ERROR: calling Nothing::";
      errmsg += __func__;
      abort_mpi(errmsg);
   }
   void fieldSolverBoundaryCondGradPeElectricField(
       FsGrid<std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> &EGradPeGrid, cint i, cint j, cint k,
       cuint component) override {
      std::string errmsg = "ERROR: calling Nothing::";
      errmsg += __func__;
      abort_mpi(errmsg);
   }
   void fieldSolverBoundaryCondDerivatives(
       FsGrid<std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> &dPerBGrid,
       FsGrid<std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> &dMomentsGrid, cint i, cint j, cint k,
       cuint RKCase, cuint component) override {
      std::string errmsg = "ERROR: calling Nothing::";
      errmsg += __func__;
      abort_mpi(errmsg);
   }
   void fieldSolverBoundaryCondBVOLDerivatives(
       FsGrid<std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> &volGrid, cint i, cint j, cint k,
       cuint component) override {
      std::string errmsg = "ERROR: calling Nothing::";
      errmsg += __func__;
      abort_mpi(errmsg);
   }
   void vlasovBoundaryCondition(const dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid,
                                const CellID &cellID, const uint popID, const bool calculate_V_moments) override {
      std::string errmsg = "ERROR: calling Nothing::";
      errmsg += __func__;
      abort_mpi(errmsg);
   }
};
} // namespace BC

#endif
