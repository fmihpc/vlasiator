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
   const auto& stencil = technicalGrid.makeStencil(i, j, k);
   std::array<Real, fsgrids::volfields::N_VOL>& vol = *volGrid.get(i, j, k);
   std::span<const std::array<Real, fsgrids::efield::N_EFIELD>> E = EGrid.getData();
   const auto sbflag = technicalGrid.get(i, j, k)->sysBoundaryFlag;

#ifdef DEBUG_FSOLVER
   const bool ok = stencil.cellExists(i, j + 1, k) && stencil.cellExists(i, j, k + 1) &&
                   stencil.cellExists(i, j + 1, k + 1) && stencil.cellExists(i + 1, j, k) &&
                   stencil.cellExists(i + 1, j + 1, k) && stencil.cellExists(i + 1, j, k + 1);

   if (ok == false) {
      std::cerr << "Out-of-bounds access in " << __FILE__ << ":" << __LINE__ << endl;
      exit(1);
   }
#endif

   // Calculate reconstruction coefficients for this cell:
   // This handles domain edges so no need to skip DO_NOT_COMPUTE or OUTER_BOUNDARY_PADDING cells.
   std::array<Real, Rec::N_REC_COEFFICIENTS> perturbedCoefficients;
   reconstructionCoefficients(perBGrid, dPerBGrid, perturbedCoefficients, i, j, k, 2);

   // Calculate volume average of B:
   vol[fsgrids::volfields::PERBXVOL] = perturbedCoefficients[Rec::a_0];
   vol[fsgrids::volfields::PERBYVOL] = perturbedCoefficients[Rec::b_0];
   vol[fsgrids::volfields::PERBZVOL] = perturbedCoefficients[Rec::c_0];

   // This avoids out of domain accesses below.
   if (sbflag == sysboundarytype::DO_NOT_COMPUTE || sbflag == sysboundarytype::OUTER_BOUNDARY_PADDING) {
      return;
   }

   if (sbflag == sysboundarytype::NOT_SYSBOUNDARY || sbflag == 1) {
      const std::array<Real, fsgrids::efield::N_EFIELD>& E_i1j1k1 = E[stencil.center()];
      const std::array<Real, fsgrids::efield::N_EFIELD>& E_i1j2k1 = E[stencil.up()];
      const std::array<Real, fsgrids::efield::N_EFIELD>& E_i1j1k2 = E[stencil.near()];
      const std::array<Real, fsgrids::efield::N_EFIELD>& E_i1j2k2 = E[stencil.upnear()];
      const std::array<Real, fsgrids::efield::N_EFIELD>& E_i2j1k1 = E[stencil.right()];
      const std::array<Real, fsgrids::efield::N_EFIELD>& E_i2j1k2 = E[stencil.rightnear()];
      const std::array<Real, fsgrids::efield::N_EFIELD>& E_i2j2k1 = E[stencil.rightup()];

      CHECK_FLOAT(E_i1j1k1[fsgrids::efield::EX]);
      CHECK_FLOAT(E_i1j1k1[fsgrids::efield::EY]);
      CHECK_FLOAT(E_i1j1k1[fsgrids::efield::EZ]);
      CHECK_FLOAT(E_i1j2k1[fsgrids::efield::EX]);
      CHECK_FLOAT(E_i1j2k1[fsgrids::efield::EZ]);
      CHECK_FLOAT(E_i1j1k2[fsgrids::efield::EX]);
      CHECK_FLOAT(E_i1j1k2[fsgrids::efield::EY]);
      CHECK_FLOAT(E_i1j2k2[fsgrids::efield::EX]);
      CHECK_FLOAT(E_i2j1k1[fsgrids::efield::EY]);
      CHECK_FLOAT(E_i2j1k1[fsgrids::efield::EZ]);
      CHECK_FLOAT(E_i2j1k2[fsgrids::efield::EY]);
      CHECK_FLOAT(E_i2j2k1[fsgrids::efield::EZ]);

      vol[fsgrids::volfields::EXVOL] = FOURTH * (E_i1j1k1[fsgrids::efield::EX] + E_i1j2k1[fsgrids::efield::EX] +
                                                 E_i1j1k2[fsgrids::efield::EX] + E_i1j2k2[fsgrids::efield::EX]);
      vol[fsgrids::volfields::EYVOL] = FOURTH * (E_i1j1k1[fsgrids::efield::EY] + E_i2j1k1[fsgrids::efield::EY] +
                                                 E_i1j1k2[fsgrids::efield::EY] + E_i2j1k2[fsgrids::efield::EY]);
      vol[fsgrids::volfields::EZVOL] = FOURTH * (E_i1j1k1[fsgrids::efield::EZ] + E_i2j1k1[fsgrids::efield::EZ] +
                                                 E_i1j2k1[fsgrids::efield::EZ] + E_i2j2k1[fsgrids::efield::EZ]);
   } else {
      vol[fsgrids::volfields::EXVOL] = 0.0;
      vol[fsgrids::volfields::EYVOL] = 0.0;
      vol[fsgrids::volfields::EZVOL] = 0.0;
   }

   CHECK_FLOAT(vol[fsgrids::volfields::EXVOL]);
   CHECK_FLOAT(vol[fsgrids::volfields::EYVOL]);
   CHECK_FLOAT(vol[fsgrids::volfields::EZVOL]);
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
