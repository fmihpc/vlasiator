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

/*!\file inflow.cpp
 * \brief Implementation of the class BoundaryCondition::Inflow.
 * This serves as the base class for further classes like BoundaryCondition::Maxwellian.
 */

#include <cstdlib>
#include <fstream>
#include <iostream>

#include "../fieldsolver/fs_common.h"
#include "../object_wrapper.h"
#include "../vlasovsolver/vlasovmover.h"
#include "inflow.h"

#ifndef NDEBUG
#define DEBUG_INFLOW
#endif
#ifdef DEBUG_BOUNDARY
#define DEBUG_INFLOW
#endif

using namespace std;

namespace BC {
Inflow::Inflow() : BoundaryCondition() {}
Inflow::~Inflow() {}

void Inflow::initBoundary(creal t, Project &project) {
   // The array of bool describes which of the faces are to have inflow boundary
   // conditions, in the order of x+, x-, y+, y-, z+, z-.
   for (uint i = 0; i < 6; i++)
      facesToProcess[i] = false;

   this->getParameters();

   vector<string>::const_iterator it;
   for (it = faceList.begin(); it != faceList.end(); ++it) {
      if (*it == "x+")
         facesToProcess[0] = true;
      if (*it == "x-")
         facesToProcess[1] = true;
      if (*it == "y+")
         facesToProcess[2] = true;
      if (*it == "y-")
         facesToProcess[3] = true;
      if (*it == "z+")
         facesToProcess[4] = true;
      if (*it == "z-")
         facesToProcess[5] = true;
   }

   for (unsigned int i = 0; i < speciesParams.size(); i++)
      loadInputData(i);

   generateTemplateCells(t);
   tLastApply = t;
}

void Inflow::assignBoundary(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid,
                            FsGrid<fsgrids::technical, FS_STENCIL_WIDTH> &technicalGrid) {
   bool doAssign;
   array<bool, 6> isThisCellOnAFace;

   vector<CellID> cells = mpiGrid.get_cells();
   for (uint i = 0; i < cells.size(); i++) {
      if (mpiGrid[cells[i]]->boundaryFlag == boundarytype::NOTHING)
         continue;
      creal *const cellParams = &(mpiGrid[cells[i]]->parameters[0]);
      creal dx = cellParams[CellParams::DX];
      creal dy = cellParams[CellParams::DY];
      creal dz = cellParams[CellParams::DZ];
      creal x = cellParams[CellParams::XCRD] + 0.5 * dx;
      creal y = cellParams[CellParams::YCRD] + 0.5 * dy;
      creal z = cellParams[CellParams::ZCRD] + 0.5 * dz;

      isThisCellOnAFace.fill(false);
      determineFace(isThisCellOnAFace.data(), x, y, z, dx, dy, dz);
      // Comparison of the array defining which faces to use and the array
      // telling on which faces this cell is
      doAssign = false;
      for (int j = 0; j < 6; j++)
         doAssign = doAssign || (facesToProcess[j] && isThisCellOnAFace[j]);
      if (doAssign) {
         mpiGrid[cells[i]]->boundaryFlag = this->getIndex();
      }
   }

   // Assign boundary flags to local fsgrid cells
   const array<int, 3> gridDims(technicalGrid.getLocalSize());
   for (int k = 0; k < gridDims[2]; k++) {
      for (int j = 0; j < gridDims[1]; j++) {
         for (int i = 0; i < gridDims[0]; i++) {
            const auto coords = technicalGrid.getPhysicalCoords(i, j, k);

            // Shift to the center of the fsgrid cell
            auto cellCenterCoords = coords;
            cellCenterCoords[0] += 0.5 * technicalGrid.DX;
            cellCenterCoords[1] += 0.5 * technicalGrid.DY;
            cellCenterCoords[2] += 0.5 * technicalGrid.DZ;
            const auto refLvl = mpiGrid.get_refinement_level(mpiGrid.get_existing_cell(cellCenterCoords));

            if (refLvl == -1)
               abort_mpi("Error, could not get refinement level of remote DCCRG cell!", 1);

            creal dx = P::dx_ini * pow(2, -refLvl);
            creal dy = P::dy_ini * pow(2, -refLvl);
            creal dz = P::dz_ini * pow(2, -refLvl);

            isThisCellOnAFace.fill(false);
            doAssign = false;

            determineFace(isThisCellOnAFace.data(), cellCenterCoords[0], cellCenterCoords[1], cellCenterCoords[2], dx,
                          dy, dz);
            for (int iface = 0; iface < 6; iface++)
               doAssign = doAssign || (facesToProcess[iface] && isThisCellOnAFace[iface]);
            if (doAssign) {
               technicalGrid.get(i, j, k)->boundaryFlag = this->getIndex();
            }
         }
      }
   }
}

void Inflow::applyInitialState(const dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid,
                               FsGrid<array<Real, fsgrids::bfield::N_BFIELD>, 2> &perBGrid, Project &project) {
   for (uint popID = 0; popID < getObjectWrapper().particleSpecies.size(); ++popID)
      setCellsFromTemplate(mpiGrid, popID);

   setBFromTemplate(mpiGrid, perBGrid);
}

void Inflow::updateState(const dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid,
                         FsGrid<array<Real, fsgrids::bfield::N_BFIELD>, 2> &perBGrid, creal t) {
   if (t - tLastApply < tInterval) {
      return;
   } else {
      tLastApply = t;
   }
#pragma omp parallel for
   for (uint i = 0; i < 6; i++) {
      if (facesToProcess[i])
         generateTemplateCell(templateCells[i], templateB[i], i, t);
   }

   for (uint popID = 0; popID < getObjectWrapper().particleSpecies.size(); ++popID)
      setCellsFromTemplate(mpiGrid, popID);

   setBFromTemplate(mpiGrid, perBGrid);
}

Real Inflow::fieldSolverBoundaryCondMagneticField(
    FsGrid<array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> &bGrid,
    FsGrid<fsgrids::technical, FS_STENCIL_WIDTH> &technicalGrid, cint i, cint j, cint k, creal dt, cuint component) {
   Real result = 0.0;
   creal dx = Parameters::dx_ini;
   creal dy = Parameters::dy_ini;
   creal dz = Parameters::dz_ini;
   const array<int, 3> globalIndices = technicalGrid.getGlobalIndices(i, j, k);
   creal x = (convert<Real>(globalIndices[0]) + 0.5) * technicalGrid.DX + Parameters::xmin;
   creal y = (convert<Real>(globalIndices[1]) + 0.5) * technicalGrid.DY + Parameters::ymin;
   creal z = (convert<Real>(globalIndices[2]) + 0.5) * technicalGrid.DZ + Parameters::zmin;

   bool isThisCellOnAFace[6];
   determineFace(&isThisCellOnAFace[0], x, y, z, dx, dy, dz, true);

   for (uint i = 0; i < 6; i++) {
      if (isThisCellOnAFace[i]) {
         result = templateB[i][component];
         break; // This effectively sets the precedence of faces through the order of faces.
      }
   }
   return result;
}

void Inflow::fieldSolverBoundaryCondElectricField(
    FsGrid<array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> &EGrid, cint i, cint j, cint k, cuint component) {
   EGrid.get(i, j, k)->at(fsgrids::efield::EX + component) = 0.0;
}

void Inflow::fieldSolverBoundaryCondHallElectricField(
    FsGrid<array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> &EHallGrid, cint i, cint j, cint k,
    cuint component) {
   array<Real, fsgrids::ehall::N_EHALL> *cp = EHallGrid.get(i, j, k);
   switch (component) {
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
      abort_mpi("Invalid component", 1);
   }
}

void Inflow::fieldSolverBoundaryCondGradPeElectricField(
    FsGrid<array<Real, fsgrids::egradpe::N_EGRADPE>, 2> &EGradPeGrid, cint i, cint j, cint k, cuint component) {
   EGradPeGrid.get(i, j, k)->at(fsgrids::egradpe::EXGRADPE + component) = 0.0;
}

void Inflow::fieldSolverBoundaryCondDerivatives(
    FsGrid<array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> &dPerBGrid,
    FsGrid<array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> &dMomentsGrid, cint i, cint j, cint k,
    cuint RKCase, cuint component) {
   this->setCellDerivativesToZero(dPerBGrid, dMomentsGrid, i, j, k, component);
}

void Inflow::fieldSolverBoundaryCondBVOLDerivatives(
    FsGrid<array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> &volGrid, cint i, cint j, cint k,
    cuint component) {
   this->setCellBVOLDerivativesToZero(volGrid, i, j, k, component);
}

void Inflow::vlasovBoundaryCondition(const dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid,
                                     const CellID &cellID, const uint popID, const bool doCalcMomentsV) {
   if (dynamic) {
      SpatialCell *cell = mpiGrid[cellID];

      creal dx = cell->parameters[CellParams::DX];
      creal dy = cell->parameters[CellParams::DY];
      creal dz = cell->parameters[CellParams::DZ];
      creal x = cell->parameters[CellParams::XCRD] + 0.5 * dx;
      creal y = cell->parameters[CellParams::YCRD] + 0.5 * dy;
      creal z = cell->parameters[CellParams::ZCRD] + 0.5 * dz;

      bool isThisCellOnAFace[6];
      determineFace(&isThisCellOnAFace[0], x, y, z, dx, dy, dz, true);

      for (uint i = 0; i < 6; i++) {
         if (facesToProcess[i] && isThisCellOnAFace[i]) {
            copyCellData(&templateCells[i], cell, false, popID, true); // copy also vdf, _V
            copyCellData(&templateCells[i], cell, true, popID, false); // don't copy vdf again but copy _R now
            break; // Effectively sets the precedence of faces through the order of faces.
         }
      }
   }
}

void Inflow::setBFromTemplate(const dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid,
                              FsGrid<array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> &perBGrid) {
   array<bool, 6> isThisCellOnAFace;
   const array<int, 3> gridDims(perBGrid.getLocalSize());

   for (int k = 0; k < gridDims[2]; k++) {
      for (int j = 0; j < gridDims[1]; j++) {
         for (int i = 0; i < gridDims[0]; i++) {
            const auto coords = perBGrid.getPhysicalCoords(i, j, k);

            // TODO: This code up to determineFace() should be in a separate
            // function, it gets called in a lot of places.
            // Shift to the center of the fsgrid cell
            auto cellCenterCoords = coords;
            cellCenterCoords[0] += 0.5 * perBGrid.DX;
            cellCenterCoords[1] += 0.5 * perBGrid.DY;
            cellCenterCoords[2] += 0.5 * perBGrid.DZ;

            const auto refLvl = mpiGrid.get_refinement_level(mpiGrid.get_existing_cell(cellCenterCoords));

            if (refLvl == -1)
               abort_mpi("Error, could not get refinement level of remote DCCRG cell!", 1);

            creal dx = P::dx_ini * pow(2, -refLvl);
            creal dy = P::dy_ini * pow(2, -refLvl);
            creal dz = P::dz_ini * pow(2, -refLvl);

            determineFace(isThisCellOnAFace.data(), cellCenterCoords[0], cellCenterCoords[1], cellCenterCoords[2], dx,
                          dy, dz);

            for (uint iface = 0; iface < 6; iface++) {
               if (facesToProcess[iface] && isThisCellOnAFace[iface]) {
                  perBGrid.get(i, j, k)->at(fsgrids::bfield::PERBX) = templateB[iface][0];
                  perBGrid.get(i, j, k)->at(fsgrids::bfield::PERBY) = templateB[iface][1];
                  perBGrid.get(i, j, k)->at(fsgrids::bfield::PERBZ) = templateB[iface][2];
                  break;
               }
            }
         }
      }
   }
}

void Inflow::setCellsFromTemplate(const dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid,
                                  const uint popID) {
   vector<CellID> cells = mpiGrid.get_cells();
#pragma omp parallel for
   for (size_t c = 0; c < cells.size(); c++) {
      SpatialCell *cell = mpiGrid[cells[c]];
      if (cell->boundaryFlag != this->getIndex())
         continue;

      creal dx = cell->parameters[CellParams::DX];
      creal dy = cell->parameters[CellParams::DY];
      creal dz = cell->parameters[CellParams::DZ];
      creal x = cell->parameters[CellParams::XCRD] + 0.5 * dx;
      creal y = cell->parameters[CellParams::YCRD] + 0.5 * dy;
      creal z = cell->parameters[CellParams::ZCRD] + 0.5 * dz;

      bool isThisCellOnAFace[6];
      determineFace(&isThisCellOnAFace[0], x, y, z, dx, dy, dz, true);

      for (uint i = 0; i < 6; i++) {
         if (facesToProcess[i] && isThisCellOnAFace[i]) {
            copyCellData(&templateCells[i], cell, false, popID, true); // copy also vdf, _V
            copyCellData(&templateCells[i], cell, true, popID, false); // don't copy vdf again but copy _R now
            break; // Effectively sets the precedence of faces through the order of faces.
         }
      }
   }
}

void Inflow::getFaces(bool *faces) {
   for (uint i = 0; i < 6; i++)
      faces[i] = facesToProcess[i];
}

void Inflow::loadInputData(const uint popID) {
   InflowSpeciesParameters &sP = speciesParams[popID];

   for (uint i = 0; i < 6; i++) {
      if (facesToProcess[i])
         sP.inputData[i] = loadFile(sP.files[i].c_str(), sP.nParams);
   }
}

/*! Load inflow boundary data from given file.
 * The number of entries per line is given by nParams which is defined as a
 * parameter from the configuration file/command line.
 * \param fn Name of the data file.
 * \retval dataset Vector of Real vectors.
 */
vector<vector<Real>> Inflow::loadFile(const char *fn, const unsigned int nParams) {
   vector<vector<Real>> dataset(0, vector<Real>(nParams, 0));

   ifstream fi;
   fi.open(fn);

   uint nlines = 0;
   string line;
   while (getline(fi, line)) {
      vector<Real> vars(nParams, -7777);
      // Skip the comments
      if (line[0] == '#')
         continue;

      stringstream ss(line);

      int i = 0;
      Real num = 0;
      while (ss >> num) {
         if (i == nParams) {
            cerr << "Extra input values at line " << nlines + 1 << " in " << fn << endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
         }
         vars[i] = num;
         i++;
      }

      for (vector<Real>::iterator v = vars.begin(); v != vars.end(); ++v) {
         if (fabs(*v + 7777.) < numeric_limits<double>::epsilon()) {
            cerr << "Missing input values at line " << nlines + 1 << " in " << fn << endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
         }
      }

      dataset.push_back(vars);
      nlines++;
   }

   if (nlines < 1) {
      cerr << "Input file " << fn << " is empty!" << endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
   } else if (nlines > 1) {
      for (uint i = 1; i < nlines; ++i) {
         if (dataset[i][0] < dataset[i - 1][0]) {
            cerr << "Parameter data must be in ascending temporal order!" << endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
         }
      }
   }

   fi.close();

   return dataset;
}

/*! Loops through the array of template cells and generates the ones needed.
 * The function generateTemplateCell is defined in the inheriting class such as
 * to have the specific condition needed.
 * \param t Simulation time
 * \sa generateTemplateCell
 */
void Inflow::generateTemplateCells(creal t) {
#pragma omp parallel for
   for (uint i = 0; i < 6; i++)
      if (facesToProcess[i])
         generateTemplateCell(templateCells[i], templateB[i], i, t);
}

/*!Interpolate the input data to the given time.
 * The first entry of each line is assumed to be the time.
 * \param inputDataIndex Index used to get the correct face's input data.
 * \param t Current simulation time.
 * \param outputData Pointer to the location where to write the result.
 * Make sure from the calling side that nParams Real values can be written
 * there!
 */
void Inflow::interpolate(const int inputDataIndex, const uint popID, creal t, Real *outputData) {
   InflowSpeciesParameters &sP = speciesParams[popID];

   // Find first data[0] value which is >= t
   int i1 = 0, i2 = 0;
   bool found = false;
   Real s; // 0 <= s < 1

   // Use the first value of data if interpolating for time before data starts
   if (t < sP.inputData[inputDataIndex][0][0]) {
      i1 = i2 = 0;
      s = 0;
   } else {
      for (uint i = 0; i < sP.inputData[inputDataIndex].size(); i++) {
         if (sP.inputData[inputDataIndex][i][0] >= t) {
            found = true;
            i2 = (int)i;
            break;
         }
      }
      if (found) {
         // i2 is now "ceil(t)"
         i1 = i2 - 1;
         if (i1 < 0) {
            i1 = i2 = 0;
            s = 0.0;
         } else {
            // normal case, now both i1 and i2 are >= 0 and < nlines, and i1 = i2-1
            s = (t - sP.inputData[inputDataIndex][i1][0]) /
                (sP.inputData[inputDataIndex][i2][0] - sP.inputData[inputDataIndex][i1][0]);
         }
      } else {
         i1 = i2 = sP.inputData[inputDataIndex].size() - 1;
         s = 0.0;
      }
   }

   creal s1 = 1 - s;

   for (uint i = 0; i < sP.nParams - 1; i++) {
      outputData[i] = s1 * sP.inputData[inputDataIndex][i1][i + 1] + s * sP.inputData[inputDataIndex][i2][i + 1];
   }
}
} // namespace BC
