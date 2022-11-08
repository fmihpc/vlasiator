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

#ifndef SETBYUSER_H
#define SETBYUSER_H

#include <vector>
#include "../definitions.h"
#include "../readparameters.h"
#include "../spatial_cell.hpp"
#include "sysboundarycondition.h"

namespace SBC {

/*!\brief Class for boundary conditions with user-set settings.
 * To be implemented by the user.
 */
class User : public OuterBoundaryCondition {
public:
   User();
   ~User() override;

   static void addParameters();
   void getParameters() override;

   virtual void initSysBoundary(creal t, Project &project);
   virtual void assignSysBoundary(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid,
                       FsGrid<fsgrids::technical, FS_STENCIL_WIDTH> &technicalGrid);
   virtual void applyInitialState(const dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid,
                          FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> &perBGrid,
                          Project &project);
   virtual void updateState(const dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid,
                    FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> &perBGrid, creal t);
   virtual Real fieldSolverBoundaryCondMagneticField(FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> &bGrid,
       FsGrid<fsgrids::technical, FS_STENCIL_WIDTH> &technicalGrid,
       cint i,
       cint j,
       cint k,
       creal dt,
       cuint component);
   virtual void fieldSolverBoundaryCondMagneticFieldProjection(
       FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> &bGrid,
       FsGrid<fsgrids::technical, FS_STENCIL_WIDTH> &technicalGrid,
       cint i,
       cint j,
       cint k);
   virtual void fieldSolverBoundaryCondElectricField(FsGrid<std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> &EGrid,
       cint i,
       cint j,
       cint k,
       cuint component);
   virtual void fieldSolverBoundaryCondHallElectricField(
       FsGrid<std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> &EHallGrid,
       cint i,
       cint j,
       cint k,
       cuint component);
   virtual void fieldSolverBoundaryCondGradPeElectricField(
       FsGrid<std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> &EGradPeGrid,
       cint i,
       cint j,
       cint k,
       cuint component);
   virtual void fieldSolverBoundaryCondDerivatives(
       FsGrid<std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> &dPerBGrid,
       FsGrid<std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> &dMomentsGrid,
       cint i,
       cint j,
       cint k,
       cuint RKCase, cuint component);
   virtual void fieldSolverBoundaryCondBVOLDerivatives(
       FsGrid<std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> &volGrid,
       cint i,
       cint j,
       cint k,
       cuint component) ;
   virtual void vlasovBoundaryCondition(const dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid,
                                const CellID &cellID, const uint popID, const bool doCalcMomentsV) ;

   virtual void getFaces(bool *faces) ;

   virtual std::string getName() const ;
   virtual uint getIndex() const ;
};
} // namespace SBC

#endif
