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

#include <vector>
#include "../definitions.h"
#include "../readparameters.h"
#include "../spatial_cells/spatial_cell_wrapper.hpp"
#include "sysboundarycondition.h"

namespace SBC {

struct InflowSpeciesParameters {
   /*! Vector containing a vector for each face which has the current boundary condition. Each of these vectors has one
    * line per input data line (time point). The length of the lines is nParams.*/
   std::vector<std::vector<Real>> inputData[6];
   /*! Input files for the inflow boundary conditions. */
   std::string files[6];

   /*! Number of parameters per input file line. */
   uint nParams;
};

/*!\brief Base class for boundary conditions with settings and parameters read from file.
 *
 * Inflow is a base class for e.g. SysBoundaryConditon::Maxwellian. It defines methods to set boundary conditions on the
 * faces of the simulation domain.
 *
 * This class handles the import and interpolation in time of the input parameters read from file as well as the
 * assignment of the state from the template cells.
 *
 * The daughter classes have then to handle parameters and generate the template cells as wished from the data returned.
 */
class Inflow : public OuterBoundaryCondition {
public:
   Inflow();
   virtual ~Inflow();

   //virtual void getParameters() = 0;

   virtual void initSysBoundary(creal& t, Project& project) override;
   virtual void applyInitialState(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid,
                                  std::span<fsgrids::technical> technical, FieldSolverGrid& fsgrid,
                                  fsgrids::perbspan perb,
                                  fsgrids::bgbspan bgb, Project& project) override;
   virtual void updateState(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid,
                            std::span<fsgrids::technical> technical, FieldSolverGrid& fsgrid,
                            fsgrids::perbspan perb,
                            fsgrids::bgbspan bgb, creal t) override;
   virtual Real fieldSolverBoundaryCondMagneticField(fsgrids::perbspan b,
                                                     std::span<const std::array<Real, fsgrids::bgbfield::N_BGB>> bgb,
                                                     std::span<const fsgrids::technical> technical,
                                                     const std::array<Real, 3>& gridSpacing,
                                                     const std::array<fsgrid::FsSize_t, 3>& globalCoordinates,
                                                     const fsgrid::FsStencil& stencil, cuint component) override;
   virtual void fieldSolverBoundaryCondElectricField(std::span<std::array<Real, fsgrids::efield::N_EFIELD>> e,
                                                     const fsgrid::FsStencil& stencil, cuint component) override;
   virtual void fieldSolverBoundaryCondHallElectricField(std::span<std::array<Real, fsgrids::ehall::N_EHALL>> ehall,
                                                         const fsgrid::FsStencil& stencil, cuint component) override;
   virtual void
   fieldSolverBoundaryCondGradPeElectricField(std::span<std::array<Real, fsgrids::egradpe::N_EGRADPE>> EGradPe,
                                              const fsgrid::FsStencil& stencil, cuint component) override;
   virtual void fieldSolverBoundaryCondDerivatives(fsgrids::dperbspan dperb,
                                                   fsgrids::dmomentsspan dmoments,
                                                   const fsgrid::FsStencil& stencil, cuint RKCase,
                                                   cuint component) override;
   virtual void fieldSolverBoundaryCondBVOLDerivatives(fsgrids::volspan vols,
                                                       const fsgrid::FsStencil& stencil, cuint component) override;
   virtual void vlasovBoundaryCondition(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid,
                                        const CellID& cellID, const uint popID, const bool doCalcMomentsV) override;
   virtual void getFaces(bool* faces) override;
   virtual std::string getName() const override {
      std::cerr << "ERROR: base class Inflow::getParameters called!" << std::endl;
      return "ERROR";
   }
   virtual uint getIndex() const = 0; // {
   //   std::cerr << "ERROR: base class Inflow::getIndex called!" << std::endl;
   //   return sysboundarytype::N_SYSBOUNDARY_CONDITIONS;
//   }

protected:
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
   std::vector<std::string> faceList;
   std::vector<InflowSpeciesParameters> speciesParams;
   std::vector<std::vector<Real>> loadFile(const char* file, unsigned int nParams);

   void loadInputData(const uint popID);
   void interpolate(const int inputDataIndex, const uint popID, creal t, Real* outputData);
   void generateTemplateCells(creal t);
   virtual void generateTemplateCell(spatial_cell::SpatialCell& templateCell, Real (&B)[3], int inputDataIndex,
                                     creal t) {
      std::cerr << "ERROR: base class Inflow::generateTemplateCell called!" << std::endl;
   }
   void setCellsFromTemplate(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid, const uint popID);
   void setBFromTemplate(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid,
                         fsgrids::perbspan perb,
                         fsgrids::bgbspan bgb,
                         std::span<fsgrids::technical> technical, FieldSolverGrid &fsgrid);
};
} // namespace SBC

#endif
