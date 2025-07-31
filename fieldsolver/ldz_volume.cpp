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

void calculateVolumeAveragedFields(fsgrids::perbspan perb,
                                   std::span<const std::array<Real, fsgrids::efield::N_EFIELD>> e,
                                   std::span<const std::array<Real, fsgrids::dperb::N_DPERB>> dperb,
                                   fsgrids::volspan vols,
                                   const fsgrid::FsStencil& stencil, cuint sysBoundaryFlag, cuint sysBoundaryLayer) {
   const auto center = stencil.ooo();
   std::array<Real, fsgrids::volfields::N_VOL>& vol = vols[center];
   const auto sbflag = sysBoundaryFlag;

#ifdef DEBUG_FSOLVER
   const bool ok = stencil.cellExists(0, 1, 0) && stencil.cellExists(0, 0, 1) && stencil.cellExists(0, 1, 1) &&
                   stencil.cellExists(1, 0, 0) && stencil.cellExists(1, 1, 0) && stencil.cellExists(1, 0, 1);

   if (ok == false) {
      std::cerr << "Out-of-bounds access in " << __FILE__ << ":" << __LINE__ << endl;
      exit(1);
   }
#endif

   // Calculate reconstruction coefficients for this cell:
   // This handles domain edges so no need to skip DO_NOT_COMPUTE or OUTER_BOUNDARY_PADDING cells.
   const auto perturbedCoefficients = reconstructionCoefficients(perb, dperb, stencil, 2);

   // Calculate volume average of B:
   vol[fsgrids::volfields::PERBXVOL] = perturbedCoefficients[Rec::a_0];
   vol[fsgrids::volfields::PERBYVOL] = perturbedCoefficients[Rec::b_0];
   vol[fsgrids::volfields::PERBZVOL] = perturbedCoefficients[Rec::c_0];

   // This avoids out of domain accesses below.
   if (sbflag == sysboundarytype::DO_NOT_COMPUTE || sbflag == sysboundarytype::OUTER_BOUNDARY_PADDING) {
      return;
   }

   if (sbflag == sysboundarytype::NOT_SYSBOUNDARY || sbflag == 1) {
      const std::array<Real, fsgrids::efield::N_EFIELD>& E_i1j1k1 = e[center];
      const std::array<Real, fsgrids::efield::N_EFIELD>& E_i1j2k1 = e[stencil.opo()];
      const std::array<Real, fsgrids::efield::N_EFIELD>& E_i1j1k2 = e[stencil.oop()];
      const std::array<Real, fsgrids::efield::N_EFIELD>& E_i1j2k2 = e[stencil.opp()];
      const std::array<Real, fsgrids::efield::N_EFIELD>& E_i2j1k1 = e[stencil.poo()];
      const std::array<Real, fsgrids::efield::N_EFIELD>& E_i2j1k2 = e[stencil.pop()];
      const std::array<Real, fsgrids::efield::N_EFIELD>& E_i2j2k1 = e[stencil.ppo()];

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

void calculateVolumeAveragedFieldsSimple(fsgrids::perbspan perb,
                                         std::span<std::array<Real, fsgrids::efield::N_EFIELD>> e,
                                         fsgrids::dperbspan dperb,
                                         fsgrids::volspan vol,
                                         fsgrids::technicalspan technical, FieldSolverGrid &fsgrid) {
   phiprof::Timer timer{"Calculate volume averaged fields"};
   const size_t numCells = fsgrid.getNumCells();
   fsgrid.parallel_for([](int timerId) -> phiprof::Timer { return phiprof::Timer{timerId}; },
                       phiprof::initializeTimer("volume averaged fields compute cells"), technical,
                       [=](const fsgrid::FsStencil& stencil, cuint sysBoundaryFlag, cuint sysBoundaryLayer) {
                          calculateVolumeAveragedFields(perb, e, dperb, vol, stencil, sysBoundaryFlag,
                                                        sysBoundaryLayer);
                       });
   timer.stop(numCells, "Spatial Cells");
}
