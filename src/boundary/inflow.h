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

#ifndef INFLOW_H
#define INFLOW_H

#include "../definitions.h"
#include "../readparameters.h"
#include "../spatial_cell.hpp"
#include "boundarycondition.h"
#include <vector>

using namespace std;

namespace BC {

struct InflowSpeciesParameters {
   /*! Vector containing a vector for each face which has the current boundary
    * condition. Each of these vectors has one line per input data line
    * (time point). The length of the lines is nParams.*/
   vector<vector<Real>> inputData[6];
   /*! Input files for the inflow boundary conditions. */
   string files[6];

   /*! Number of space- and velocityspace samples used when creating phase space
    * densities */
   uint nSpaceSamples;
   uint nVelocitySamples;

   /*! Number of parameters per input file line. */
   uint nParams;
};

/*!\brief Base class for boundary conditions with settings and parameters read
 * from file.
 *
 * Inflow is a base class for e.g. BoundaryConditon::SetMaxwellian. It defines
 * the managing functions to set boundary conditions on the faces of the
 * simulation domain.
 *
 * This class handles the import and interpolation in time of the input
 * parameters read from file as well as the assignment of the state from the
 * template cells.
 *
 * The daughter classes have then to handle parameters and generate the
 * template cells as wished from the data returned.
 */
class Inflow : public BoundaryCondition {
public:
   Inflow();
   virtual ~Inflow();

   virtual void getParameters() = 0;

   virtual void initBoundary(creal t, Project &project);
   virtual void assignBoundary(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid,
                               FsGrid<fsgrids::technical, FS_STENCIL_WIDTH> &technicalGrid);
   virtual void applyInitialState(const dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid,
                                  FsGrid<array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> &perBGrid,
                                  Project &project);
   virtual void updateState(const dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid,
                            FsGrid<array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> &perBGrid, creal t);
   virtual Real
   fieldSolverBoundaryCondMagneticField(FsGrid<array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> &bGrid,
                                        FsGrid<fsgrids::technical, FS_STENCIL_WIDTH> &technicalGrid, cint i, cint j,
                                        cint k, creal dt, cuint component);
   virtual void
   fieldSolverBoundaryCondElectricField(FsGrid<array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> &EGrid, cint i,
                                        cint j, cint k, cuint component);
   virtual void
   fieldSolverBoundaryCondHallElectricField(FsGrid<array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> &EHallGrid,
                                            cint i, cint j, cint k, cuint component);
   virtual void fieldSolverBoundaryCondGradPeElectricField(
       FsGrid<array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> &EGradPeGrid, cint i, cint j, cint k,
       cuint component);
   virtual void fieldSolverBoundaryCondDerivatives(
       FsGrid<array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> &dPerBGrid,
       FsGrid<array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> &dMomentsGrid, cint i, cint j, cint k,
       cuint RKCase, cuint component);
   virtual void
   fieldSolverBoundaryCondBVOLDerivatives(FsGrid<array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> &volGrid,
                                          cint i, cint j, cint k, cuint component);
   virtual void vlasovBoundaryCondition(const dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid,
                                        const CellID &cellID, const uint popID, const bool doCalcMomentsV);

   virtual void getFaces(bool *faces);

   virtual string getName() const = 0;
   virtual uint getIndex() const = 0;

protected:
   /*! Array of bool telling which faces are going to be processed by the boundary condition.*/
   bool facesToProcess[6];
   /*! Array of template spatial cells replicated over the corresponding
    * simulation volume face. Only the template for an active face is actually
    * being touched at all by the code. */
   spatial_cell::SpatialCell templateCells[6];
   Real templateB[6][3];
   /*! Time interval for applying the dynamic BC. */
   Real tInterval;
   /*! Last simulation time the dynamic BC is applied. */
   Real tLastApply;
   /*! List of faces on which inflow boundary conditions are to be applied ([xyz][+-]). */
   vector<string> faceList;
   vector<InflowSpeciesParameters> speciesParams;
   vector<vector<Real>> loadFile(const char *file, unsigned int nParams);

   void loadInputData(const uint popID);
   void interpolate(const int inputDataIndex, const uint popID, creal t, Real *outputData);
   void generateTemplateCells(creal t);
   virtual void generateTemplateCell(spatial_cell::SpatialCell &templateCell, Real (&B)[3], int inputDataIndex,
                                     creal t) = 0;
   void setCellsFromTemplate(const dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid, const uint popID);
   void setBFromTemplate(const dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid,
                         FsGrid<array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> &perBGrid);
};
} // namespace BC

#endif
