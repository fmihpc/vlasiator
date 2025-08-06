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
 * \brief Implementation of the class SysBoundaryCondition::Inflow.
 * This serves as the base class for inherited classes like SysBoundaryCondition::Maxwellian.
 */

#include <cstdlib>
#include <fstream>
#include <iostream>

#include "../fieldsolver/fs_common.h"
#include "../grid.h"
#include "../object_wrapper.h"
#include "../vlasovsolver/vlasovmover.h"
#include "inflow.h"

#ifdef DEBUG_VLASIATOR
#define DEBUG_INFLOW
#endif
#ifdef DEBUG_SYSBOUNDARY
#define DEBUG_INFLOW
#endif

namespace SBC {
   
   Inflow::Inflow() : OuterBoundaryCondition() {}
   Inflow::~Inflow() {}
   
   void Inflow::initSysBoundary(creal& t, Project& project) {
      // The array of bool describes which of the faces are to have inflow boundary
      // conditions, in the order of x+, x-, y+, y-, z+, z-.
      for(uint i=0; i<6; i++) {
         facesToProcess[i] = false;
      }
   
      this->getParameters();
   
      for (auto& it : faceList) {
         if (it == "x+") {
            facesToProcess[0] = true;
         } else if (it == "x-") {
            facesToProcess[1] = true;
         } else if (it == "y+") {
            facesToProcess[2] = true;
         } else if (it == "y-") {
            facesToProcess[3] = true;
         } else if (it == "z+") {
            facesToProcess[4] = true;
         } else if (it == "z-") {
            facesToProcess[5] = true;
         }
      }
   
      for (unsigned int i = 0; i < speciesParams.size(); i++) {
         loadInputData(i);
      }
   
      generateTemplateCells(t);
      tLastApply = t;
   }
   
   void Inflow::applyInitialState(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid,
                                  fsgrids::technicalspan technical, FieldSolverGrid &fsgrid,
                                  fsgrids::perbspan perb,
                                  fsgrids::bgbspan bgb, Project& project) {
      for (uint popID = 0; popID < getObjectWrapper().particleSpecies.size(); ++popID) {
         setCellsFromTemplate(mpiGrid, popID);
      }
      setBFromTemplate(mpiGrid, perb, bgb, technical, fsgrid);
   }
   
   void Inflow::updateState(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid,
                            fsgrids::technicalspan technical, FieldSolverGrid &fsgrid,
                            fsgrids::perbspan perb,
                            fsgrids::bgbspan bgb, creal t) {
      if (t - tLastApply < tInterval) {
         return;
      } else {
         tLastApply = t;
      }
      for (uint i = 0; i < 6; i++) {
         if (facesToProcess[i]) {
            generateTemplateCell(templateCells[i], templateB[i], i, t);
         }
      }
      for (uint popID = 0; popID < getObjectWrapper().particleSpecies.size(); ++popID) {
         setCellsFromTemplate(mpiGrid, popID);
      }
   
      setBFromTemplate(mpiGrid, perb, bgb, technical, fsgrid);
   
      // Ensure up-to-date velocity block counts for all neighbours
      for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
         updateRemoteVelocityBlockLists(mpiGrid,popID,Neighborhoods::FULL);
      }
   }
   
   Real Inflow::fieldSolverBoundaryCondMagneticField(fsgrids::perbspan b,
                                                     fsgrids::constbgbspan bgb,
                                                     fsgrids::consttechnicalspan technical,
                                                     const std::array<Real, 3>& gridSpacing,
                                                     const std::array<fsgrid::FsSize_t, 3>& globalCoordinates,
                                                     const fsgrid::FsStencil& stencil, cuint component) {
      Real result = 0.0;
      creal dx = Parameters::dx_ini;
      creal dy = Parameters::dy_ini;
      creal dz = Parameters::dz_ini;
      creal x = (convert<Real>(globalCoordinates[0]) + 0.5) * gridSpacing[0] + Parameters::xmin;
      creal y = (convert<Real>(globalCoordinates[1]) + 0.5) * gridSpacing[1] + Parameters::ymin;
      creal z = (convert<Real>(globalCoordinates[2]) + 0.5) * gridSpacing[2] + Parameters::zmin;
   
      bool isThisCellOnAFace[6];
      determineFace(&isThisCellOnAFace[0], x, y, z, dx, dy, dz, true);
   
      for (uint i = 0; i < 6; i++) {
         if (isThisCellOnAFace[i]) {
            result = templateB[i][component];
            break; // This effectively sets the precedence of faces through the order of faces.
         }
      }
   
      // There are projects that have non-uniform and non-zero perturbed B, e.g. Magnetosphere with dipole type 4.
      // We cannot jsut take the value from the templateCell, we also need a copy of the value from initialization.
      // This value is stored in the BgBGrid at fsgrids::bgbfield::BGBXVDCORR,BGBYVDCORR,BGBZVDCORR
      result += bgb[stencil.ooo()][fsgrids::bgbfield::BGBXVDCORR + component];
      return result;
   }
   
   void Inflow::fieldSolverBoundaryCondElectricField(fsgrids::efieldspan e,
                                                     const fsgrid::FsStencil& stencil, cuint component) {
      e[stencil.ooo()][fsgrids::efield::EX + component] = 0.0;
   }
   
   void Inflow::fieldSolverBoundaryCondHallElectricField(fsgrids::ehallspan ehall,
                                                         const fsgrid::FsStencil& stencil, cuint component) {
      std::array<Real, fsgrids::ehall::N_EHALL>& cp = ehall[stencil.ooo()];
      switch (component) {
      case 0:
         cp[fsgrids::ehall::EXHALL_000_100] = 0.0;
         cp[fsgrids::ehall::EXHALL_010_110] = 0.0;
         cp[fsgrids::ehall::EXHALL_001_101] = 0.0;
         cp[fsgrids::ehall::EXHALL_011_111] = 0.0;
         break;
      case 1:
         cp[fsgrids::ehall::EYHALL_000_010] = 0.0;
         cp[fsgrids::ehall::EYHALL_100_110] = 0.0;
         cp[fsgrids::ehall::EYHALL_001_011] = 0.0;
         cp[fsgrids::ehall::EYHALL_101_111] = 0.0;
         break;
      case 2:
         cp[fsgrids::ehall::EZHALL_000_001] = 0.0;
         cp[fsgrids::ehall::EZHALL_100_101] = 0.0;
         cp[fsgrids::ehall::EZHALL_010_011] = 0.0;
         cp[fsgrids::ehall::EZHALL_110_111] = 0.0;
         break;
      default:
         abort_mpi("Invalid component", 1);
      }
   }
   
   void Inflow::fieldSolverBoundaryCondGradPeElectricField(
      fsgrids::egradpespan EGradPe, const fsgrid::FsStencil& stencil,
      cuint component) {
      EGradPe[stencil.ooo()][fsgrids::egradpe::EXGRADPE + component] = 0.0;
   }
   
   void Inflow::fieldSolverBoundaryCondDerivatives(fsgrids::dperbspan dperb,
                                                   fsgrids::dmomentsspan dmoments,
                                                   const fsgrid::FsStencil& stencil, cuint RKCase, cuint component) {
      this->setCellDerivativesToZero(dperb, dmoments, stencil, component);
   }
   
   void Inflow::fieldSolverBoundaryCondBVOLDerivatives(fsgrids::volspan vols,
                                                       const fsgrid::FsStencil& stencil, cuint component) {
      this->setCellBVOLDerivativesToZero(vols, stencil, component);
   }
   
   void Inflow::vlasovBoundaryCondition(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid,
                                        const CellID& cellID, const uint popID, const bool doCalcMomentsV) {
      // This is a no-op because both template cell generation and block data copying takes place in
      // updateState() (at pre-set intervals only)
   }
   
   void determineFaceNoClassMembersInflow(
      bool* isThisCellOnAFace,
      creal x, creal y, creal z,
      creal dx, creal dy, creal dz,
      const std::array<bool, 3> periodicity,
      const bool excludeSlicesAndPeriodicDimensions=false // (default)
   ) {   
      for(uint i=0; i<6; i++) {
         isThisCellOnAFace[i] = false;
      }     
      if(x > Parameters::xmax - dx * 2) {
         isThisCellOnAFace[0] = true;
      }  
      if(x < Parameters::xmin + dx * 2) {
         isThisCellOnAFace[1] = true;
      }     
      if(y > Parameters::ymax - dy * 2) {
         isThisCellOnAFace[2] = true;
      }     
      if(y < Parameters::ymin + dy * 2) {
         isThisCellOnAFace[3] = true;
      }
      if(z > Parameters::zmax - dz * 2) {
         isThisCellOnAFace[4] = true;
      }
      if(z < Parameters::zmin + dz * 2) {
         isThisCellOnAFace[5] = true;
      }
      if(excludeSlicesAndPeriodicDimensions == true) {
         if (Parameters::xcells_ini == 1 || periodicity[0]) {
            isThisCellOnAFace[0] = false;
            isThisCellOnAFace[1] = false; 
         }     
         if (Parameters::ycells_ini == 1 || periodicity[1]) {
            isThisCellOnAFace[2] = false;
            isThisCellOnAFace[3] = false;
         }     
         if (Parameters::zcells_ini == 1 || periodicity[2]) {
            isThisCellOnAFace[4] = false; 
            isThisCellOnAFace[5] = false; 
         }     
      }
   }

   void Inflow::setBFromTemplate(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid,
                                 fsgrids::perbspan perb,
                                 fsgrids::bgbspan bgb,
                                 fsgrids::technicalspan technical, FieldSolverGrid &fsgrid) {
      const auto& gridSpacing = fsgrid.getGridSpacing();
      const auto facesToProcess_local = this->facesToProcess;
      std::array<bool, 3> periodic_local = this->periodic;
      const auto templateB_local = this->templateB;
      fsgrid.parallel_for_coords([](int timerId) -> phiprof::Timer { return phiprof::Timer{timerId}; },
                                 phiprof::initializeTimer("setBFromTemplate"), technical,
                                 [=](const fsgrid::FsStencil& stencil, cuint sysBoundaryFlag, cuint sysBoundaryLayer, std::array<Real, 3> coords) {
         creal dx = P::dx_ini * pow(2, -technical[stencil.ooo()].refLevel);
         creal dy = P::dy_ini * pow(2, -technical[stencil.ooo()].refLevel);
         creal dz = P::dz_ini * pow(2, -technical[stencil.ooo()].refLevel);

         std::array<bool, 6> isThisCellOnAFace = {{false}};
 
         determineFaceNoClassMembersInflow(isThisCellOnAFace.data(), coords[0] + 0.5 * gridSpacing[0], coords[1] + 0.5 * gridSpacing[1], coords[2] + 0.5 * gridSpacing[2], dx, dy, dz, periodic_local);
   
         for (uint iface = 0; iface < 6; iface++) {
            if (facesToProcess_local[iface] && isThisCellOnAFace[iface]) {
               perb[stencil.ooo()][fsgrids::bfield::PERBX] = templateB_local[iface][0] + bgb[stencil.ooo()][fsgrids::bgbfield::BGBXVDCORR];
               perb[stencil.ooo()][fsgrids::bfield::PERBY] = templateB_local[iface][1] + bgb[stencil.ooo()][fsgrids::bgbfield::BGBYVDCORR];
               perb[stencil.ooo()][fsgrids::bfield::PERBZ] = templateB_local[iface][2] + bgb[stencil.ooo()][fsgrids::bgbfield::BGBZVDCORR];
               break;
            }
         }
      });
   }
   
   void Inflow::setCellsFromTemplate(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid, const uint popID) {
      // Assign boundary flags to local DCCRG cells
      const std::vector<CellID>& cells = getLocalCells();
      #pragma omp parallel for schedule(dynamic, 1)
      for (size_t c = 0; c < cells.size(); c++) {
         SpatialCell* cell = mpiGrid[cells[c]];
         if (cell->sysBoundaryFlag != this->getIndex()) {
            continue;
         }
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
   #ifdef USE_GPU
               cell->setReservation(popID, templateCells[i].getReservation(popID));
   #endif
               break; // Effectively sets the precedence of faces through the order of faces.
            }
         }
      }
      // Verify current mesh and blocks
      // if (!cell->checkMesh(popID)) {
      //    printf("ERROR in vmesh check: %s at %d\n",__FILE__,__LINE__);
      // }
   }
   
   void Inflow::getFaces(bool* faces) {
      for (uint i = 0; i < 6; i++) {
         faces[i] = facesToProcess[i];
      }
   }
   
   void Inflow::loadInputData(const uint popID) {
      InflowSpeciesParameters& sP = speciesParams[popID];
   
      for (uint i = 0; i < 6; i++) {
         if (facesToProcess[i])
            sP.inputData[i] = loadFile(sP.files[i].c_str(), sP.nParams);
      }
   }
   
   /*! Load inflow boundary data from given file.
    * The number of entries per line is given by nParams which is defined as a parameter
    * from the configuration file/command line.
    * \param fn Name of the data file.
    * \retval dataset Vector of Real vectors.
    */
   vector<std::vector<Real>> Inflow::loadFile(const char* fn, const unsigned int nParams) {
      vector<std::vector<Real>> dataset(0, std::vector<Real>(nParams, 0));
   
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
            if (i == int(nParams)) {
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
    * The function generateTemplateCell is defined in the child class to have the specific condition needed.
    * \param t Simulation time
    * \sa generateTemplateCell
    */
   void Inflow::generateTemplateCells(creal t) {
      for (uint i = 0; i < 6; i++) {
         if (facesToProcess[i]) {
            generateTemplateCell(templateCells[i], templateB[i], i, t);
         }
      }
   }
   
   /*!Interpolate the input data to the given time.
    * The first entry of each line is assumed to be the time.
    * \param inputDataIndex Index used to get the correct face's input data.
    * \param t Current simulation time.
    * \param outputData Pointer to the location where to write the result.
    * Make sure from the calling side that nParams Real values can be written there!
    */
   void Inflow::interpolate(const int inputDataIndex, const uint popID, creal t, Real* outputData) {
      InflowSpeciesParameters& sP = speciesParams[popID];
   
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
   
} // namespace SBC
