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

/*!\file outflow.cpp
 * \brief Implementation of the class BoundaryCondition::Outflow to handle cells classified as boundarytype::OUTFLOW.
 */

#include <cstdlib>
#include <iostream>

#include "../fieldsolver/fs_common.h"
#include "../fieldsolver/ldz_magnetic_field.hpp"
#include "../object_wrapper.h"
#include "../projects/projects_common.h"
#include "../vlasovsolver/vlasovmover.h"
#include "outflow.h"

#ifndef NDEBUG
#define DEBUG_OUTFLOW
#endif
#ifdef DEBUG_BOUNDARY
#define DEBUG_OUTFLOW
#endif

using namespace std;

namespace BC
{
Outflow::Outflow() : BoundaryCondition() {}
Outflow::~Outflow() {}

void Outflow::addParameters()
{
   const std::string defStr = "Copy";
   Readparameters::addComposing(
       "outflow.faceNoFields",
       "List of faces on which no field outflow boundary conditions are to be applied ([xyz][+-]).");
   Readparameters::add("outflow.precedence",
                       "Precedence value of the outflow boundary condition (integer), the higher the stronger.", 4);
   Readparameters::add("outflow.reapplyUponRestart",
                       "If 0 (default), keep going with the state existing in the restart file. If 1, calls again "
                       "applyInitialState. Can be used to change boundary condition behaviour during a run.",
                       0);

   // Per-population parameters
   for (uint i = 0; i < getObjectWrapper().particleSpecies.size(); i++)
   {
      const std::string &pop = getObjectWrapper().particleSpecies[i].name;

      Readparameters::addComposing(
          pop + "_outflow.reapplyFaceUponRestart",
          "List of faces on which outflow boundary conditions are to be reapplied upon restart ([xyz][+-]).");
      Readparameters::addComposing(pop + "_outflow.face",
                                   "List of faces on which outflow boundary conditions are to be applied ([xyz][+-]).");
      Readparameters::add(pop + "_outflow.vlasovScheme_face_x+", "Scheme to use on the face x+ (Copy, Limit, None)",
                          defStr);
      Readparameters::add(pop + "_outflow.vlasovScheme_face_x-", "Scheme to use on the face x- (Copy, Limit, None)",
                          defStr);
      Readparameters::add(pop + "_outflow.vlasovScheme_face_y+", "Scheme to use on the face y+ (Copy, Limit, None)",
                          defStr);
      Readparameters::add(pop + "_outflow.vlasovScheme_face_y-", "Scheme to use on the face y- (Copy, Limit, None)",
                          defStr);
      Readparameters::add(pop + "_outflow.vlasovScheme_face_z+", "Scheme to use on the face z+ (Copy, Limit, None)",
                          defStr);
      Readparameters::add(pop + "_outflow.vlasovScheme_face_z-", "Scheme to use on the face z- (Copy, Limit, None)",
                          defStr);

      Readparameters::add(pop + "_outflow.quench",
                          "Factor by which to quench the inflowing parts of the velocity distribution function.", 1.0);
   }
}

void Outflow::getParameters()
{

   Readparameters::get("outflow.faceNoFields", this->faceNoFieldsList);
   Readparameters::get("outflow.precedence", precedence);
   uint reapply;
   Readparameters::get("outflow.reapplyUponRestart", reapply);
   this->applyUponRestart = false;
   if (reapply == 1)
   {
      this->applyUponRestart = true;
   }

   // Per-species parameters
   for (uint i = 0; i < getObjectWrapper().particleSpecies.size(); i++)
   {
      const std::string &pop = getObjectWrapper().particleSpecies[i].name;
      OutflowSpeciesParameters sP;

      // Unless we find out otherwise, we assume that this species will not be treated at any boundary
      for (int j = 0; j < 6; j++)
      {
         sP.facesToSkipVlasov[j] = true;
      }

      std::vector<std::string> thisSpeciesFaceList;
      Readparameters::get(pop + "_outflow.face", thisSpeciesFaceList);

      for (auto &face : thisSpeciesFaceList)
      {
         if (face == "x+")
         {
            facesToProcess[0] = true;
            sP.facesToSkipVlasov[0] = false;
         }
         if (face == "x-")
         {
            facesToProcess[1] = true;
            sP.facesToSkipVlasov[1] = false;
         }
         if (face == "y+")
         {
            facesToProcess[2] = true;
            sP.facesToSkipVlasov[2] = false;
         }
         if (face == "y-")
         {
            facesToProcess[3] = true;
            sP.facesToSkipVlasov[3] = false;
         }
         if (face == "z+")
         {
            facesToProcess[4] = true;
            sP.facesToSkipVlasov[4] = false;
         }
         if (face == "z-")
         {
            facesToProcess[5] = true;
            sP.facesToSkipVlasov[5] = false;
         }
      }

      Readparameters::get(pop + "_outflow.reapplyFaceUponRestart", sP.faceToReapplyUponRestartList);
      std::array<std::string, 6> vlasovBoundarySchemeName;
      Readparameters::get(pop + "_outflow.vlasovScheme_face_x+", vlasovBoundarySchemeName[0]);
      Readparameters::get(pop + "_outflow.vlasovScheme_face_x-", vlasovBoundarySchemeName[1]);
      Readparameters::get(pop + "_outflow.vlasovScheme_face_y+", vlasovBoundarySchemeName[2]);

      Readparameters::get(pop + "_outflow.vlasovScheme_face_y-", vlasovBoundarySchemeName[3]);
      Readparameters::get(pop + "_outflow.vlasovScheme_face_z+", vlasovBoundarySchemeName[4]);
      Readparameters::get(pop + "_outflow.vlasovScheme_face_z-", vlasovBoundarySchemeName[5]);
      for (uint j = 0; j < 6; j++)
      {
         if (vlasovBoundarySchemeName[j] == "None")
         {
            sP.faceVlasovScheme[j] = vlasovscheme::NONE;
         }
         else if (vlasovBoundarySchemeName[j] == "Copy")
         {
            sP.faceVlasovScheme[j] = vlasovscheme::COPY;
         }
         else if (vlasovBoundarySchemeName[j] == "Limit")
         {
            sP.faceVlasovScheme[j] = vlasovscheme::LIMIT;
         }
         else
         {
            abort_mpi(" ERROR: " + vlasovBoundarySchemeName[j] + " is an invalid Outflow Vlasov scheme!");
         }
      }

      Readparameters::get(pop + "_outflow.quench", sP.quenchFactor);

      speciesParams.push_back(sP);
   }
}

void Outflow::initBoundary(creal t, Project &project)
{
   // The array of bool describes which of the faces are to have outflow
   // boundary conditions, in the order of x+, x-, y+, y-, z+, z-.
   for (uint i = 0; i < 6; i++)
   {
      facesToProcess[i] = false;
      facesToSkipFields[i] = false;
      facesToReapply[i] = false;
   }

   this->getParameters();

   dynamic = false;

   vector<string>::const_iterator it;
   for (it = faceNoFieldsList.begin(); it != faceNoFieldsList.end(); it++)
   {
      if (*it == "x+") facesToSkipFields[0] = true;
      if (*it == "x-") facesToSkipFields[1] = true;
      if (*it == "y+") facesToSkipFields[2] = true;
      if (*it == "y-") facesToSkipFields[3] = true;
      if (*it == "z+") facesToSkipFields[4] = true;
      if (*it == "z-") facesToSkipFields[5] = true;
   }

   for (uint i = 0; i < getObjectWrapper().particleSpecies.size(); i++)
   {
      OutflowSpeciesParameters &sP = this->speciesParams[i];
      for (it = sP.faceToReapplyUponRestartList.begin(); it != sP.faceToReapplyUponRestartList.end(); it++)
      {
         if (*it == "x+") facesToReapply[0] = true;
         if (*it == "x-") facesToReapply[1] = true;
         if (*it == "y+") facesToReapply[2] = true;
         if (*it == "y-") facesToReapply[3] = true;
         if (*it == "z+") facesToReapply[4] = true;
         if (*it == "z-") facesToReapply[5] = true;
      }
   }
}

void Outflow::assignBoundary(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid,
                             FsGrid<fsgrids::technical, FS_STENCIL_WIDTH> &technicalGrid)
{

   bool doAssign;
   std::array<bool, 6> isThisCellOnAFace;

   // Assign boundary flags to local DCCRG cells
   vector<CellID> cells = mpiGrid.get_cells();
   for (const auto &dccrgId : cells)
   {
      if (mpiGrid[dccrgId]->boundaryFlag == boundarytype::NOTHING) continue;
      creal *const cellParams = &(mpiGrid[dccrgId]->parameters[0]);
      creal dx = cellParams[CellParams::DX];
      creal dy = cellParams[CellParams::DY];
      creal dz = cellParams[CellParams::DZ];
      creal x = cellParams[CellParams::XCRD] + 0.5 * dx;
      creal y = cellParams[CellParams::YCRD] + 0.5 * dy;
      creal z = cellParams[CellParams::ZCRD] + 0.5 * dz;

      isThisCellOnAFace.fill(false);
      determineFace(isThisCellOnAFace.data(), x, y, z, dx, dy, dz);

      // Comparison of the array defining which faces to use and the array telling on which faces this cell is
      doAssign = false;
      for (int j = 0; j < 6; j++)
         doAssign = doAssign || (facesToProcess[j] && isThisCellOnAFace[j]);
      if (doAssign)
      {
         mpiGrid[dccrgId]->boundaryFlag = this->getIndex();
      }
   }

   // Assign boundary flags to local fsgrid cells
   const std::array<int, 3> gridDims(technicalGrid.getLocalSize());
   for (int k = 0; k < gridDims[2]; k++)
   {
      for (int j = 0; j < gridDims[1]; j++)
      {
         for (int i = 0; i < gridDims[0]; i++)
         {
            const auto &coords = technicalGrid.getPhysicalCoords(i, j, k);

            // Shift to the center of the fsgrid cell
            auto cellCenterCoords = coords;
            cellCenterCoords[0] += 0.5 * technicalGrid.DX;
            cellCenterCoords[1] += 0.5 * technicalGrid.DY;
            cellCenterCoords[2] += 0.5 * technicalGrid.DZ;
            const auto refLvl = mpiGrid.get_refinement_level(mpiGrid.get_existing_cell(cellCenterCoords));

            if (refLvl == -1) abort_mpi("ERROR: Could not get refinement level of remote DCCRG cell!", 1);

            creal dx = P::dx_ini * pow(2, -refLvl);
            creal dy = P::dy_ini * pow(2, -refLvl);
            creal dz = P::dz_ini * pow(2, -refLvl);

            isThisCellOnAFace.fill(false);
            doAssign = false;

            determineFace(isThisCellOnAFace.data(), cellCenterCoords[0], cellCenterCoords[1], cellCenterCoords[2], dx,
                          dy, dz);
            for (int iface = 0; iface < 6; iface++)
               doAssign = doAssign || (facesToProcess[iface] && isThisCellOnAFace[iface]);
            if (doAssign)
            {
               technicalGrid.get(i, j, k)->boundaryFlag = this->getIndex();
            }
         }
      }
   }
}

void Outflow::applyInitialState(const dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid,
                                FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> &perBGrid,
                                Project &project) {
   const vector<CellID> &cells = getLocalCells();
#pragma omp parallel for
   for (uint iC = 0; iC < cells.size(); ++iC) {
      SpatialCell *cell = mpiGrid[cells[iC]];
      if (cell->boundaryFlag != this->getIndex())
         continue;

      bool doApply = true;

      if (Parameters::isRestart) {
         creal *const cellParams = &(mpiGrid[cells[iC]]->parameters[0]);
         creal dx = cellParams[CellParams::DX];
         creal dy = cellParams[CellParams::DY];
         creal dz = cellParams[CellParams::DZ];
         creal x = cellParams[CellParams::XCRD] + 0.5 * dx;
         creal y = cellParams[CellParams::YCRD] + 0.5 * dy;
         creal z = cellParams[CellParams::ZCRD] + 0.5 * dz;

         bool isThisCellOnAFace[6];
         determineFace(&isThisCellOnAFace[0], x, y, z, dx, dy, dz);

         doApply = false;
         // Comparison of the array defining which faces to use and the array telling on which faces this cell is
         for (uint i = 0; i < getObjectWrapper().particleSpecies.size(); i++) {
            for (uint j = 0; j < 6; j++) {
               doApply = doApply || (facesToReapply[j] && isThisCellOnAFace[j]);
            }
         }
      }

      if (doApply) {
         // Defined in project.cpp, used here as the outflow cell has the same state
         // as the initial state of non-boundary cells.
         project.setCell(cell);

         cell->parameters[CellParams::RHOM_DT2] = cell->parameters[CellParams::RHOM];
         cell->parameters[CellParams::RHOQ_DT2] = cell->parameters[CellParams::RHOQ];
         cell->parameters[CellParams::VX_DT2] = cell->parameters[CellParams::VX];
         cell->parameters[CellParams::VY_DT2] = cell->parameters[CellParams::VY];
         cell->parameters[CellParams::VZ_DT2] = cell->parameters[CellParams::VZ];
         cell->parameters[CellParams::P_11_DT2] = cell->parameters[CellParams::P_11];
         cell->parameters[CellParams::P_22_DT2] = cell->parameters[CellParams::P_22];
         cell->parameters[CellParams::P_33_DT2] = cell->parameters[CellParams::P_33];
      }
   }
}

void Outflow::updateState(const dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid,
                          FsGrid<array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> &perBGrid, creal t)
{
}

Real Outflow::fieldSolverBoundaryCondMagneticField(
    FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> &bGrid,
    FsGrid<fsgrids::technical, FS_STENCIL_WIDTH> &technicalGrid, cint i, cint j, cint k, creal dt, cuint component)
{
   switch (component)
   {
   case 0:
      return fieldBoundaryCopyFromSolvingNbrMagneticField(bGrid, technicalGrid, i, j, k, component, compute::BX);
      break;
   case 1:
      return fieldBoundaryCopyFromSolvingNbrMagneticField(bGrid, technicalGrid, i, j, k, component, compute::BY);
      break;
   case 2:
      return fieldBoundaryCopyFromSolvingNbrMagneticField(bGrid, technicalGrid, i, j, k, component, compute::BZ);
      break;
   default:
      return 0.0;
      break;
   }
}

void Outflow::fieldSolverBoundaryCondElectricField(
    FsGrid<std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> &EGrid, cint i, cint j, cint k,
    cuint component)
{
   EGrid.get(i, j, k)->at(fsgrids::efield::EX + component) = 0.0;
}

void Outflow::fieldSolverBoundaryCondHallElectricField(
    FsGrid<std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> &EHallGrid, cint i, cint j, cint k,
    cuint component)
{
   std::array<Real, fsgrids::ehall::N_EHALL> *cp = EHallGrid.get(i, j, k);
   switch (component)
   {
   case 0:
      cp->at(fsgrids::ehall::EXHALL_000_100) = 0.0;
      cp->at(fsgrids::ehall::EXHALL_010_110) = 0.0;
      cp->at(fsgrids::ehall::EXHALL_001_101) = 0.0;
      cp->at(fsgrids::ehall::EXHALL_011_111) = 0.0;
      break;
   case 1:
      cp->at(fsgrids::ehall::EYHALL_000_010) = 0.0;
      cp->at(fsgrids::ehall::EYHALL_100_110) = 0.0;
      cp->at(fsgrids::ehall::EYHALL_001_011) = 0.0;
      cp->at(fsgrids::ehall::EYHALL_101_111) = 0.0;
      break;
   case 2:
      cp->at(fsgrids::ehall::EZHALL_000_001) = 0.0;
      cp->at(fsgrids::ehall::EZHALL_100_101) = 0.0;
      cp->at(fsgrids::ehall::EZHALL_010_011) = 0.0;
      cp->at(fsgrids::ehall::EZHALL_110_111) = 0.0;
      break;
   default:
      abort_mpi("Invalid component");
   }
}

void Outflow::fieldSolverBoundaryCondGradPeElectricField(
    FsGrid<std::array<Real, fsgrids::egradpe::N_EGRADPE>, 2> &EGradPeGrid, cint i, cint j, cint k, cuint component)
{
   EGradPeGrid.get(i, j, k)->at(fsgrids::egradpe::EXGRADPE + component) = 0.0;
}

void Outflow::fieldSolverBoundaryCondDerivatives(
    FsGrid<std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> &dPerBGrid,
    FsGrid<std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> &dMomentsGrid, cint i, cint j, cint k,
    cuint RKCase, cuint component)
{
   this->setCellDerivativesToZero(dPerBGrid, dMomentsGrid, i, j, k, component);
}

void Outflow::fieldSolverBoundaryCondBVOLDerivatives(
    FsGrid<std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> &volGrid, cint i, cint j, cint k,
    cuint component)
{
   this->setCellBVOLDerivativesToZero(volGrid, i, j, k, component);
}

/**
 * NOTE that this is called once for each particle species!
 * @param mpiGrid
 * @param cellID
 */
void Outflow::vlasovBoundaryCondition(const dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid,
                                      const CellID &cellID, const uint popID, const bool doCalcMoments_V)
{
   const OutflowSpeciesParameters &sP = this->speciesParams[popID];
   SpatialCell *cell = mpiGrid[cellID];
   creal *const cellParams = cell->parameters.data();
   creal dx = cellParams[CellParams::DX];
   creal dy = cellParams[CellParams::DY];
   creal dz = cellParams[CellParams::DZ];
   creal x = cellParams[CellParams::XCRD] + 0.5 * dx;
   creal y = cellParams[CellParams::YCRD] + 0.5 * dy;
   creal z = cellParams[CellParams::ZCRD] + 0.5 * dz;

   bool isThisCellOnAFace[6];
   determineFace(&isThisCellOnAFace[0], x, y, z, dx, dy, dz, true);

   for (uint i = 0; i < 6; i++)
   {
      if (isThisCellOnAFace[i] && facesToProcess[i] && !sP.facesToSkipVlasov[i])
      {
         switch (sP.faceVlasovScheme[i])
         {
         case vlasovscheme::NONE:
            break;
         case vlasovscheme::COPY:
            vlasovBoundaryCopyFromClosestNbr(mpiGrid, cellID, false, popID, doCalcMoments_V);
            break;
         case vlasovscheme::LIMIT:
            vlasovBoundaryCopyFromClosestNbrAndLimit(mpiGrid, cellID, popID);
            break;
         default:
            abort_mpi("ERROR: invalid Outflow Vlasov scheme", 1);
            break;
         }
      }
   }
}

void Outflow::getFaces(bool *faces)
{
   for (uint i = 0; i < 6; i++)
      faces[i] = facesToProcess[i];
}

std::string Outflow::getName() const { return "Outflow"; }
uint Outflow::getIndex() const { return boundarytype::OUTFLOW; }

} // namespace BC
