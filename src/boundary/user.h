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

#ifndef USER_H
#define USER_H

#include "../definitions.h"
#include "../readparameters.h"
#include "../spatial_cell.hpp"
#include "boundarycondition.h"
#include <vector>

namespace BC
{

/*!\brief Class for boundary conditions with user-set settings.
 * To be implemented by the user.
 */
class User : public BoundaryCondition
{
public:
   User();
   ~User() override;

   static void addParameters();
   void getParameters() override;

   bool initBoundary(creal &t, Project &project) override;
   bool assignBoundary(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid,
                       FsGrid<fsgrids::technical, 2> &technicalGrid) override;
   bool applyInitialState(const dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid,
                          FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, 2> &perBGrid, Project &project) override;
   Real fieldSolverBoundaryCondMagneticField(FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, 2> &bGrid,
                                             FsGrid<fsgrids::technical, 2> &technicalGrid, cint i, cint j, cint k,
                                             creal &dt, cuint &component) override;
   void fieldSolverBoundaryCondElectricField(FsGrid<std::array<Real, fsgrids::efield::N_EFIELD>, 2> &EGrid, cint i,
                                             cint j, cint k, cuint component) override;
   void fieldSolverBoundaryCondHallElectricField(FsGrid<std::array<Real, fsgrids::ehall::N_EHALL>, 2> &EHallGrid,
                                                 cint i, cint j, cint k, cuint component) override;
   void
   fieldSolverBoundaryCondGradPeElectricField(FsGrid<std::array<Real, fsgrids::egradpe::N_EGRADPE>, 2> &EGradPeGrid,
                                              cint i, cint j, cint k, cuint component) override;
   void fieldSolverBoundaryCondDerivatives(FsGrid<std::array<Real, fsgrids::dperb::N_DPERB>, 2> &dPerBGrid,
                                           FsGrid<std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 2> &dMomentsGrid,
                                           cint i, cint j, cint k, cuint &RKCase, cuint &component) override;
   void fieldSolverBoundaryCondBVOLDerivatives(FsGrid<std::array<Real, fsgrids::volfields::N_VOL>, 2> &volGrid, cint i,
                                               cint j, cint k, cuint &component) override;
   void vlasovBoundaryCondition(const dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid,
                                const CellID &cellID, const uint popID, const bool calculate_V_moments) override{
       // No need to do anything in this function, as the propagators do not touch the distribution function
   };

   void getFaces(bool *faces) override;

   std::string getName() const override;
   uint getIndex() const override;
};
} // namespace BC

#endif
