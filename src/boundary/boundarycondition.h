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

#ifndef BOUNDARYCONDITION_H
#define BOUNDARYCONDITION_H

#include <dccrg.hpp>
#include <dccrg_cartesian_geometry.hpp>
#include <fsgrid.hpp>

#include "../common.h"
#include "../definitions.h"
#include "../projects/project.h"
#include "../spatial_cell.hpp"
#include <vector>

using namespace spatial_cell;
using namespace projects;

namespace BC {
/*!\brief BC::BoundaryCondition is the base class for boundary conditions.
 *
 * BC::BoundaryCondition defines a base class for applying boundary
 * conditions. Specific boundary conditions inherit from this base class.
 *
 * The initBoundary function is used to initialise the internal workings
 * needed by the boundary condition to run (e.g. importing parameters,
 * initialising class members). assignBoundary is used to determine whether a
 * given cell is within the domain of boundary condition. applyInitialState
 * is called to initialise a boundary cell's parameters and velocity space.
 *
 * If needed, a user can write his or her own BC::BoundaryConditions, which
 * are loaded when the simulation initializes.
 */
class BoundaryCondition {
public:
   BoundaryCondition();
   virtual ~BoundaryCondition();

   static void addParameters();
   virtual void getParameters() = 0;

   virtual void initBoundary(creal t, Project &project) = 0;
   virtual void assignBoundary(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid,
                               FsGrid<fsgrids::technical, FS_STENCIL_WIDTH> &technicalGrid) = 0;
   virtual void applyInitialState(const dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid,
                                  FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> &perBGrid,
                                  Project &project) = 0;
   virtual void updateState(const dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid,
                            FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> &perBGrid,
                            creal t) = 0;
   virtual Real
   fieldSolverBoundaryCondMagneticField(FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> &bGrid,
                                        FsGrid<fsgrids::technical, FS_STENCIL_WIDTH> &technicalGrid, cint i, cint j,
                                        cint k, creal dt, cuint component) = 0;
   virtual void
   fieldSolverBoundaryCondElectricField(FsGrid<std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> &EGrid,
                                        cint i, cint j, cint k, cuint component) = 0;
   virtual void fieldSolverBoundaryCondHallElectricField(
       FsGrid<std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> &EHallGrid, cint i, cint j, cint k,
       cuint component) = 0;
   virtual void fieldSolverBoundaryCondGradPeElectricField(
       FsGrid<std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> &EGradPeGrid, cint i, cint j, cint k,
       cuint component) = 0;
   virtual void fieldSolverBoundaryCondDerivatives(
       FsGrid<std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> &dPerBGrid,
       FsGrid<std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> &dMomentsGrid, cint i, cint j, cint k,
       cuint RKCase, cuint component) = 0;
   virtual void fieldSolverBoundaryCondBVOLDerivatives(
       FsGrid<std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> &volGrid, cint i, cint j, cint k,
       cuint component) = 0;
   static void
   setCellDerivativesToZero(FsGrid<std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> &dPerBGrid,
                            FsGrid<std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> &dMomentsGrid,
                            cint i, cint j, cint k, cuint component);
   static void
   setCellBVOLDerivativesToZero(FsGrid<std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> &volGrid, cint i,
                                cint j, cint k, cuint component);

   /** This function computes the Vlasov (distribution function)
    * boundary condition for the given particle species only.
    * It is NOT allowed to change block structure in cell.
    * @param mpiGrid Parallel grid.
    * @param cellID Spatial cell ID.
    * @param popID Particle species ID.*/
   virtual void vlasovBoundaryCondition(const dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid,
                                        const CellID &cellID, const uint popID, const bool doCalcMomentsV) = 0;

   virtual void getFaces(bool *faces);
   virtual std::string getName() const = 0;
   virtual uint getIndex() const = 0;
   uint getPrecedence() const;
   bool isDynamic() const;

   void updateBoundaryConditionsAfterLoadBalance(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid,
                                                 const std::vector<CellID> &local_cells_on_boundary);
   bool doApplyUponRestart() const;
   void setPeriodicity(bool isFacePeriodic[3]);

protected:
   /*! Precedence value of the boundary condition. */
   uint precedence;
   /*! Is the boundary condition dynamic in time or not. */
   bool dynamic;
   /*! Array of bool telling whether the system is periodic in any direction. */
   bool periodic[3];
   /*! Map of closest nonboundarycells. Used in getAllClosestNonboundaryCells. */
   std::unordered_map<CellID, std::vector<CellID>> allClosestNonboundaryCells;
   /*! Map of close nonboundarycells. Used in getAllCloseNonboundaryCells. */
   std::unordered_map<CellID, std::vector<CellID>> allCloseNonboundaryCells;

   /*! Array of cells into which the distribution function can flow. Used in getAllFlowtoCells. Cells into which one
    * cannot flow are set to INVALID_CELLID. */
   std::unordered_map<CellID, std::array<SpatialCell *, 27>> allFlowtoCells;
   /*! bool telling whether to call again applyInitialState upon restarting the simulation. */
   bool applyUponRestart;

   void determineFace(bool *isThisCellOnAFace, creal x, creal y, creal z, creal dx, creal dy, creal dz,
                      const bool excludeSlicesAndPeriodicDimensions = false);
   void copyCellData(SpatialCell *from, SpatialCell *to, const bool copyMomentsOnly, const uint popID,
                     const bool doCalcMomentsV);
   void averageCellData(const dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid,
                        std::vector<CellID> cellList, SpatialCell *to, const uint popID, const bool doCalcMomentsV,
                        creal fluffiness = 0);
   std::array<SpatialCell *, 27> &getFlowtoCells(const CellID &cellID);

   std::array<Realf *, 27> getFlowtoCellsBlock(const std::array<SpatialCell *, 27> flowtoCells,
                                               const vmesh::GlobalID blockGID, const uint popID);

   /*! Helper function to get the index of a neighboring cell in the arrays in allFlowtoCells.
    * \param i Offset in x direction (-1, 0 or 1)
    * \param j Offset in y direction (-1, 0 or 1)
    * \param k Offset in z direction (-1, 0 or 1)
    * \retval int Index in the flowto cell array (0 to 26, indexed from - to + x, y, z.
    */
   inline int nbrID(const int i, const int j, const int k) { return (k + 1) * 9 + (j + 1) * 3 + i + 1; }

   void vlasovBoundaryCopyFromClosestNbr(const dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid,
                                         const CellID &cellID, const bool &copyMomentsOnly, const uint popID,
                                         const bool doCalcMomentsV);
   void vlasovBoundaryCopyFromClosestNbrAndLimit(const dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid,
                                                 const CellID &cellID, const uint popID);
   void vlasovBoundaryCopyFromAllClosestNbrs(const dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid,
                                             const CellID &cellID, const uint popID, const bool doCalcMomentsV);
   void vlasovBoundaryFluffyCopyFromAllCloseNbrs(const dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid,
                                                 const CellID &cellID, const uint popID, const bool doCalcMomentsV,
                                                 creal fluffiness);
   void vlasovBoundaryReflect(const dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid, const CellID &cellID,
                              creal &nx, creal &ny, creal &nz, const uint popID);
   void vlasovBoundaryAbsorb(const dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid, const CellID &cellID,
                             creal &nx, creal &ny, creal &nz, creal &quenchingFactor, const uint popID);
   std::array<int, 3> getTheClosestNonboundaryCell(FsGrid<fsgrids::technical, FS_STENCIL_WIDTH> &technicalGrid, cint i,
                                                   cint j, cint k);
   std::vector<std::array<int, 3>>
   getAllClosestNonboundaryCells(FsGrid<fsgrids::technical, FS_STENCIL_WIDTH> &technicalGrid, cint i, cint j, cint k);
   CellID &getTheClosestNonboundaryCell(const CellID &cellID);
   std::vector<CellID> &getAllClosestNonboundaryCells(const CellID &cellID);
   std::vector<CellID> &getAllCloseNonboundaryCells(const CellID &cellID);
   Real fieldBoundaryCopyFromSolvingNbrMagneticField(
       FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> &bGrid,
       FsGrid<fsgrids::technical, FS_STENCIL_WIDTH> &technicalGrid, cint i, cint j, cint k, cuint component,
       cuint mask);
};

void abort_mpi(const std::string str, const int err_type = 0);

} // namespace BC

#endif
