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

#ifndef OUTFLOW_H
#define OUTFLOW_H

#include "../definitions.h"
#include "../readparameters.h"
#include "../spatial_cell.hpp"
#include "boundarycondition.h"
#include <vector>

namespace BC {

struct OutflowSpeciesParameters {
   /*! Array of bool telling which faces are going to be skipped by the Vlasov boundary condition.*/
   std::array<bool, 6> facesToSkipVlasov;
   /*! List of schemes to use for the Vlasov outflow boundary conditions on each face ([xyz][+-]). */
   std::array<uint, 6> faceVlasovScheme;
   /*! List of faces on which outflow boundary conditions are to be reapplied upon restart ([xyz][+-]). */
   std::vector<std::string> faceToReapplyUponRestartList;

   /*! Factor by which to quench the inflowing parts of the velocity distribution function.*/
   Real quenchFactor;
};

/*!\brief Outflow is a class applying copy/outflow boundary conditions.
 *
 * Outflow is a class handling cells tagged as boundarytype::OUTFLOW by this
 * boundary condition. It applies copy/outflow boundary conditions.
 *
 * These consist in:
 * - Copy the distribution and moments from the nearest NOT_BOUNDARY cell;
 * - Copy the perturbed B components from the nearest NOT_BOUNDARY cell.
 * EXCEPTION: the face components adjacent to the simulation domain at the
 * +x/+y/+z faces are propagated still.
 */
class Outflow : public BoundaryCondition {
public:
   Outflow();
   ~Outflow() override;

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
   Real
   fieldSolverBoundaryCondMagneticField(FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> &bGrid,
                                        FsGrid<fsgrids::technical, FS_STENCIL_WIDTH> &technicalGrid, cint i, cint j,
                                        cint k, creal dt, cuint component) override;
   void
   fieldSolverBoundaryCondElectricField(FsGrid<std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> &EGrid,
                                        cint i, cint j, cint k, cuint component) override;
   void fieldSolverBoundaryCondHallElectricField(
       FsGrid<std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> &EHallGrid, cint i, cint j, cint k,
       cuint component) override;
   void fieldSolverBoundaryCondGradPeElectricField(
       FsGrid<std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> &EGradPeGrid, cint i, cint j, cint k,
       cuint component) override;
   void fieldSolverBoundaryCondDerivatives(
       FsGrid<std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> &dPerBGrid,
       FsGrid<std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> &dMomentsGrid, cint i, cint j, cint k,
       cuint RKCase, cuint component) override;
   void fieldSolverBoundaryCondBVOLDerivatives(
       FsGrid<std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> &volGrid, cint i, cint j, cint k,
       cuint component) override;
   void vlasovBoundaryCondition(const dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid,
                                const CellID &cellID, const uint popID, const bool doCalcMomentsV) override;

   void getFaces(bool *faces) override;
   std::string getName() const override;
   uint getIndex() const override;

protected:
   /*! Array of bool telling which faces are going to be processed by the boundary condition.*/
   bool facesToProcess[6];
   /*! Array of bool telling which faces are going to be processed by the fields boundary condition.*/
   bool facesToSkipFields[6];
   /*! Array of bool telling which faces are going to be reapplied upon restart.*/
   bool facesToReapply[6];
   /*! List of faces on which outflow boundary conditions are to be applied ([xyz][+-]). */
   std::vector<std::string> faceList;
   /*! List of faces on which no fields outflow boundary conditions are to be applied ([xyz][+-]). */
   std::vector<std::string> faceNoFieldsList;
   std::vector<OutflowSpeciesParameters> speciesParams;

   /*! Factor by which to quench the inflowing parts of the velocity distribution function.*/
   Real quenchFactor;

   enum vlasovscheme { NONE, COPY, LIMIT, N_SCHEMES };

}; // class Outflow
} // namespace BC

#endif
