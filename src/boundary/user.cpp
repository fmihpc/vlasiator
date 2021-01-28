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

/*!\file user.cpp
 * \brief Implementation of the class BoundaryCondition::User.
 * This is the interface for customized boundary implementation.
 */

#include <cstdlib>
#include <fstream>
#include <iostream>

#include "../fieldsolver/fs_common.h"
#include "../object_wrapper.h"
#include "../vlasovmover.h"
#include "user.h"

#ifndef NDEBUG
#define DEBUG_USER
#endif
#ifdef DEBUG_BOUNDARY
#define DEBUG_USER
#endif

using namespace std;

namespace BC
{
User::User() : BoundaryCondition() {}
User::~User() {}

void User::addParameters() {}

void User::initBoundary(creal t, Project &project) {}

void User::getParameters() {}

void User::assignBoundary(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid,
                          FsGrid<fsgrids::technical, 2> &technicalGrid)
{
}

void User::applyInitialState(const dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid,
                             FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, 2> &perBGrid, Project &project)
{
}

void User::vlasovBoundaryCondition(const dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid,
                                   const CellID &cellID, const uint popID, const bool doCalcMomentsV, creal t)
{
}

Real User::fieldSolverBoundaryCondMagneticField(FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, 2> &bGrid,
                                                FsGrid<fsgrids::technical, 2> &technicalGrid,
                                                cint i, cint j, cint k, creal dt, cuint component)
{
}

void User::fieldSolverBoundaryCondElectricField(FsGrid<std::array<Real, fsgrids::efield::N_EFIELD>, 2> &EGrid, cint i,
                                                cint j, cint k, cuint component)
{
}

void User::fieldSolverBoundaryCondHallElectricField(FsGrid<std::array<Real, fsgrids::ehall::N_EHALL>, 2> &EHallGrid,
                                                    cint i, cint j, cint k, cuint component)
{
}

void User::fieldSolverBoundaryCondGradPeElectricField(
    FsGrid<std::array<Real, fsgrids::egradpe::N_EGRADPE>, 2> &EGradPeGrid, cint i, cint j, cint k, cuint component)
{
}

void User::fieldSolverBoundaryCondDerivatives(FsGrid<std::array<Real, fsgrids::dperb::N_DPERB>, 2> &dPerBGrid,
                                              FsGrid<std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 2> &dMomentsGrid,
                                              cint i, cint j, cint k, cuint RKCase, cuint component)
{
}

void User::fieldSolverBoundaryCondBVOLDerivatives(FsGrid<std::array<Real, fsgrids::volfields::N_VOL>, 2> &volGrid,
                                                  cint i, cint j, cint k, cuint component)
{
}

void User::getFaces(bool *faces) {}

string User::getName() const { return "User"; }
uint User::getIndex() const { return boundarytype::USER; }

} // namespace BC
