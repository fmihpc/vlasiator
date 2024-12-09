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

#include <cstdlib>

#include "fs_common.h"
#include "ldz_volume.hpp"

#ifdef DEBUG_VLASIATOR
   #define DEBUG_FSOLVER
#endif

void calculateVolumeAveragedFields(
    fsgrid::FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH>& perBGrid,
    fsgrid::FsGrid<std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH>& EGrid,
    fsgrid::FsGrid<std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH>& dPerBGrid,
    fsgrid::FsGrid<std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH>& volGrid,
    fsgrid::FsGrid<fsgrids::technical, FS_STENCIL_WIDTH>& technicalGrid, int32_t i, int32_t j, int32_t k) {
   std::array<Real, fsgrids::volfields::N_VOL>& vol = *volGrid.get(i, j, k);

   // Calculate reconstruction coefficients for this cell:
   // This handles domain edges so no need to skip DO_NOT_COMPUTE or OUTER_BOUNDARY_PADDING cells.
   std::array<Real, Rec::N_REC_COEFFICIENTS> perturbedCoefficients;
   reconstructionCoefficients(perBGrid, dPerBGrid, perturbedCoefficients, i, j, k, 2);

   // Calculate volume average of B:
   vol[fsgrids::volfields::PERBXVOL] = perturbedCoefficients[Rec::a_0];
   vol[fsgrids::volfields::PERBYVOL] = perturbedCoefficients[Rec::b_0];
   vol[fsgrids::volfields::PERBZVOL] = perturbedCoefficients[Rec::c_0];

   // This avoids out of domain accesses below.
   if (technicalGrid.get(i, j, k)->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE ||
       technicalGrid.get(i, j, k)->sysBoundaryFlag == sysboundarytype::OUTER_BOUNDARY_PADDING) {
      return;
   }
   // Calculate volume average of E (FIXME NEEDS IMPROVEMENT):
   std::array<Real, fsgrids::efield::N_EFIELD>* EGrid_i1j1k1 = EGrid.get(i, j, k);
   if (technicalGrid.get(i, j, k)->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY ||
       (technicalGrid.get(i, j, k)->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY &&
        technicalGrid.get(i, j, k)->sysBoundaryLayer == 1)) {
#ifdef DEBUG_FSOLVER
      bool ok = true;
      if (technicalGrid.get(i, j + 1, k) == NULL)
         ok = false;
      if (technicalGrid.get(i, j, k + 1) == NULL)
         ok = false;
      if (technicalGrid.get(i, j + 1, k + 1) == NULL)
         ok = false;
      if (ok == false) {
         stringstream ss;
         ss << "ERROR, got NULL neighbor in " << __FILE__ << ":" << __LINE__ << endl;
         cerr << ss.str();
         exit(1);
      }
#endif

      std::array<Real, fsgrids::efield::N_EFIELD>* EGrid_i1j2k1 = EGrid.get(i, j + 1, k);
      std::array<Real, fsgrids::efield::N_EFIELD>* EGrid_i1j1k2 = EGrid.get(i, j, k + 1);
      std::array<Real, fsgrids::efield::N_EFIELD>* EGrid_i1j2k2 = EGrid.get(i, j + 1, k + 1);

      CHECK_FLOAT(EGrid_i1j1k1->at(fsgrids::efield::EX));
      CHECK_FLOAT(EGrid_i1j2k1->at(fsgrids::efield::EX));
      CHECK_FLOAT(EGrid_i1j1k2->at(fsgrids::efield::EX));
      CHECK_FLOAT(EGrid_i1j2k2->at(fsgrids::efield::EX));
      vol[fsgrids::volfields::EXVOL] =
          FOURTH * (EGrid_i1j1k1->at(fsgrids::efield::EX) + EGrid_i1j2k1->at(fsgrids::efield::EX) +
                    EGrid_i1j1k2->at(fsgrids::efield::EX) + EGrid_i1j2k2->at(fsgrids::efield::EX));
      CHECK_FLOAT(vol[fsgrids::volfields::EXVOL]);
   } else {
      vol[fsgrids::volfields::EXVOL] = 0.0;
   }

   if (technicalGrid.get(i, j, k)->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY ||
       (technicalGrid.get(i, j, k)->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY &&
        technicalGrid.get(i, j, k)->sysBoundaryLayer == 1)) {
#ifdef DEBUG_FSOLVER
      bool ok = true;
      if (technicalGrid.get(i + 1, j, k) == NULL)
         ok = false;
      if (technicalGrid.get(i, j, k + 1) == NULL)
         ok = false;
      if (technicalGrid.get(i + 1, j, k + 1) == NULL)
         ok = false;
      if (ok == false) {
         stringstream ss;
         ss << "ERROR, got NULL neighbor in " << __FILE__ << ":" << __LINE__ << endl;
         cerr << ss.str();
         exit(1);
      }
#endif

      std::array<Real, fsgrids::efield::N_EFIELD>* EGrid_i2j1k1 = EGrid.get(i + 1, j, k);
      std::array<Real, fsgrids::efield::N_EFIELD>* EGrid_i1j1k2 = EGrid.get(i, j, k + 1);
      std::array<Real, fsgrids::efield::N_EFIELD>* EGrid_i2j1k2 = EGrid.get(i + 1, j, k + 1);

      CHECK_FLOAT(EGrid_i1j1k1->at(fsgrids::efield::EY));
      CHECK_FLOAT(EGrid_i2j1k1->at(fsgrids::efield::EY));
      CHECK_FLOAT(EGrid_i1j1k2->at(fsgrids::efield::EY));
      CHECK_FLOAT(EGrid_i2j1k2->at(fsgrids::efield::EY));
      vol[fsgrids::volfields::EYVOL] =
          FOURTH * (EGrid_i1j1k1->at(fsgrids::efield::EY) + EGrid_i2j1k1->at(fsgrids::efield::EY) +
                    EGrid_i1j1k2->at(fsgrids::efield::EY) + EGrid_i2j1k2->at(fsgrids::efield::EY));
      CHECK_FLOAT(vol[fsgrids::volfields::EYVOL]);
   } else {
      vol[fsgrids::volfields::EYVOL] = 0.0;
   }

   if (technicalGrid.get(i, j, k)->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY ||
       (technicalGrid.get(i, j, k)->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY &&
        technicalGrid.get(i, j, k)->sysBoundaryLayer == 1)) {
#ifdef DEBUG_FSOLVER
      bool ok = true;
      if (technicalGrid.get(i + 1, j, k) == NULL)
         ok = false;
      if (technicalGrid.get(i, j + 1, k) == NULL)
         ok = false;
      if (technicalGrid.get(i + 1, j + 1, k) == NULL)
         ok = false;
      if (ok == false) {
         stringstream ss;
         ss << "ERROR, got NULL neighbor in " << __FILE__ << ":" << __LINE__ << endl;
         cerr << ss.str();
         exit(1);
      }
#endif

      std::array<Real, fsgrids::efield::N_EFIELD>* EGrid_i2j1k1 = EGrid.get(i + 1, j, k);
      std::array<Real, fsgrids::efield::N_EFIELD>* EGrid_i1j2k1 = EGrid.get(i, j + 1, k);
      std::array<Real, fsgrids::efield::N_EFIELD>* EGrid_i2j2k1 = EGrid.get(i + 1, j + 1, k);

      CHECK_FLOAT(EGrid_i1j1k1->at(fsgrids::efield::EZ));
      CHECK_FLOAT(EGrid_i2j1k1->at(fsgrids::efield::EZ));
      CHECK_FLOAT(EGrid_i1j2k1->at(fsgrids::efield::EZ));
      CHECK_FLOAT(EGrid_i2j2k1->at(fsgrids::efield::EZ));
      vol[fsgrids::volfields::EZVOL] =
          FOURTH * (EGrid_i1j1k1->at(fsgrids::efield::EZ) + EGrid_i2j1k1->at(fsgrids::efield::EZ) +
                    EGrid_i1j2k1->at(fsgrids::efield::EZ) + EGrid_i2j2k1->at(fsgrids::efield::EZ));
      CHECK_FLOAT(vol[fsgrids::volfields::EZVOL]);
   } else {
      vol[fsgrids::volfields::EZVOL] = 0.0;
   }
}

void calculateVolumeAveragedFieldsSimple(
    fsgrid::FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH>& perBGrid,
    fsgrid::FsGrid<std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH>& EGrid,
    fsgrid::FsGrid<std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH>& dPerBGrid,
    fsgrid::FsGrid<std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH>& volGrid,
    fsgrid::FsGrid<fsgrids::technical, FS_STENCIL_WIDTH>& technicalGrid) {
   const auto& localSize = technicalGrid.getLocalSize();
   const size_t N_cells = localSize[0] * localSize[1] * localSize[2];

   phiprof::Timer timer{"Calculate volume averaged fields"};
   int parallelTimerId{phiprof::initializeTimer("volume averaged fields compute cells")};
#pragma omp parallel
   {
      phiprof::Timer parallelTimer{parallelTimerId};
#pragma omp for collapse(2)
      for (fsgrid::FsIndex_t k = 0; k < localSize[2]; k++) {
         for (fsgrid::FsIndex_t j = 0; j < localSize[1]; j++) {
            for (fsgrid::FsIndex_t i = 0; i < localSize[0]; i++) {
               calculateVolumeAveragedFields(perBGrid, EGrid, dPerBGrid, volGrid, technicalGrid, i, j, k);
            }
         }
      }
   }

   timer.stop(N_cells, "Spatial Cells");
}
