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

#ifndef NDEBUG
   #define DEBUG_FSOLVER
#endif

using namespace std;

void calculateVolumeAveragedFields(
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2> & perBGrid,
   FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 2> & EGrid,
   FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, 2> & dPerBGrid,
   FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, 2> & volGrid,
   FsGrid< fsgrids::technical, 2> & technicalGrid
) {
   //const std::array<int, 3> gridDims = technicalGrid.getLocalSize();
   const int* gridDims = &technicalGrid.getLocalSize()[0];
   const size_t N_cells = gridDims[0]*gridDims[1]*gridDims[2];
   phiprof::start("Calculate volume averaged fields");
   
   #pragma omp parallel for collapse(3)
   for (int k=0; k<gridDims[2]; k++) {
      for (int j=0; j<gridDims[1]; j++) {
         for (int i=0; i<gridDims[0]; i++) {
            if(technicalGrid.get(i,j,k)->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) continue;
            
            Real perturbedCoefficients[Rec::N_REC_COEFFICIENTS];
            std::array<Real, fsgrids::volfields::N_VOL> * volGrid0 = volGrid.get(i,j,k);
            
            // Calculate reconstruction coefficients for this cell:
            reconstructionCoefficients(
               perBGrid,
               dPerBGrid,
               perturbedCoefficients,
               i,
               j,
               k,
               2
            );
            
            // Calculate volume average of B:
            volGrid0->at(fsgrids::volfields::PERBXVOL) = perturbedCoefficients[Rec::a_0];
            volGrid0->at(fsgrids::volfields::PERBYVOL) = perturbedCoefficients[Rec::b_0];
            volGrid0->at(fsgrids::volfields::PERBZVOL) = perturbedCoefficients[Rec::c_0];

            // Calculate volume average of E (FIXME NEEDS IMPROVEMENT):
            std::array<Real, fsgrids::efield::N_EFIELD> * EGrid_i1j1k1 = EGrid.get(i,j,k);
            if ( technicalGrid.get(i,j,k)->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY ||
                (technicalGrid.get(i,j,k)->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY && technicalGrid.get(i,j,k)->sysBoundaryLayer == 1)
            ) {
               #ifdef DEBUG_FSOLVER
               bool ok = true;
               if (technicalGrid.get(i  ,j+1,k  ) == NULL) ok = false;
               if (technicalGrid.get(i  ,j  ,k+1) == NULL) ok = false;
               if (technicalGrid.get(i  ,j+1,k+1) == NULL) ok = false;
               if (ok == false) {
                  stringstream ss;
                  ss << "ERROR, got NULL neighbor in " << __FILE__ << ":" << __LINE__ << endl;
                  cerr << ss.str(); exit(1);
               }
               #endif

               std::array<Real, fsgrids::efield::N_EFIELD> * EGrid_i1j2k1 = EGrid.get(i  ,j+1,k  );
               std::array<Real, fsgrids::efield::N_EFIELD> * EGrid_i1j1k2 = EGrid.get(i  ,j  ,k+1);
               std::array<Real, fsgrids::efield::N_EFIELD> * EGrid_i1j2k2 = EGrid.get(i  ,j+1,k+1);

               CHECK_FLOAT(EGrid_i1j1k1->at(fsgrids::efield::EX))
               CHECK_FLOAT(EGrid_i1j2k1->at(fsgrids::efield::EX))
               CHECK_FLOAT(EGrid_i1j1k2->at(fsgrids::efield::EX))
               CHECK_FLOAT(EGrid_i1j2k2->at(fsgrids::efield::EX))
               volGrid0->at(fsgrids::volfields::EXVOL) = FOURTH*(EGrid_i1j1k1->at(fsgrids::efield::EX) + EGrid_i1j2k1->at(fsgrids::efield::EX) + EGrid_i1j1k2->at(fsgrids::efield::EX) + EGrid_i1j2k2->at(fsgrids::efield::EX));
               CHECK_FLOAT(volGrid0->at(fsgrids::volfields::EXVOL))
            } else {
               volGrid0->at(fsgrids::volfields::EXVOL) = 0.0;
            }

            if ( technicalGrid.get(i,j,k)->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY ||
                 (technicalGrid.get(i,j,k)->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY && technicalGrid.get(i,j,k)->sysBoundaryLayer == 1)
            ) {
               #ifdef DEBUG_FSOLVER
               bool ok = true;
               if (technicalGrid.get(i+1,j  ,k  ) == NULL) ok = false;
               if (technicalGrid.get(i  ,j  ,k+1) == NULL) ok = false;
               if (technicalGrid.get(i+1,j  ,k+1) == NULL) ok = false;
               if (ok == false) {
                  stringstream ss;
                  ss << "ERROR, got NULL neighbor in " << __FILE__ << ":" << __LINE__ << endl;
                  cerr << ss.str(); exit(1);
               }
               #endif

               std::array<Real, fsgrids::efield::N_EFIELD> * EGrid_i2j1k1 = EGrid.get(i+1,j  ,k  );
               std::array<Real, fsgrids::efield::N_EFIELD> * EGrid_i1j1k2 = EGrid.get(i  ,j  ,k+1);
               std::array<Real, fsgrids::efield::N_EFIELD> * EGrid_i2j1k2 = EGrid.get(i+1,j  ,k+1);

               CHECK_FLOAT(EGrid_i1j1k1->at(fsgrids::efield::EY))
               CHECK_FLOAT(EGrid_i2j1k1->at(fsgrids::efield::EY))
               CHECK_FLOAT(EGrid_i1j1k2->at(fsgrids::efield::EY))
               CHECK_FLOAT(EGrid_i2j1k2->at(fsgrids::efield::EY))
               volGrid0->at(fsgrids::volfields::EYVOL) = FOURTH*(EGrid_i1j1k1->at(fsgrids::efield::EY) + EGrid_i2j1k1->at(fsgrids::efield::EY) + EGrid_i1j1k2->at(fsgrids::efield::EY) + EGrid_i2j1k2->at(fsgrids::efield::EY));
               CHECK_FLOAT(volGrid0->at(fsgrids::volfields::EYVOL))
            } else {
               volGrid0->at(fsgrids::volfields::EYVOL) = 0.0;
            }

            if ( technicalGrid.get(i,j,k)->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY ||
                (technicalGrid.get(i,j,k)->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY && technicalGrid.get(i,j,k)->sysBoundaryLayer == 1)
            ) {
               #ifdef DEBUG_FSOLVER
               bool ok = true;
               if (technicalGrid.get(i+1,j  ,k  ) == NULL) ok = false;
               if (technicalGrid.get(i  ,j+1,k  ) == NULL) ok = false;
               if (technicalGrid.get(i+1,j+1,k  ) == NULL) ok = false;
               if (ok == false) {
                  stringstream ss;
                  ss << "ERROR, got NULL neighbor in " << __FILE__ << ":" << __LINE__ << endl;
                  cerr << ss.str(); exit(1);
               }
               #endif

               std::array<Real, fsgrids::efield::N_EFIELD> * EGrid_i2j1k1 = EGrid.get(i+1,j  ,k  );
               std::array<Real, fsgrids::efield::N_EFIELD> * EGrid_i1j2k1 = EGrid.get(i  ,j+1,k  );
               std::array<Real, fsgrids::efield::N_EFIELD> * EGrid_i2j2k1 = EGrid.get(i+1,j+1,k  );

               CHECK_FLOAT(EGrid_i1j1k1->at(fsgrids::efield::EZ))
               CHECK_FLOAT(EGrid_i2j1k1->at(fsgrids::efield::EZ))
               CHECK_FLOAT(EGrid_i1j2k1->at(fsgrids::efield::EZ))
               CHECK_FLOAT(EGrid_i2j2k1->at(fsgrids::efield::EZ))
               volGrid0->at(fsgrids::volfields::EZVOL) = FOURTH*(EGrid_i1j1k1->at(fsgrids::efield::EZ) + EGrid_i2j1k1->at(fsgrids::efield::EZ) + EGrid_i1j2k1->at(fsgrids::efield::EZ) + EGrid_i2j2k1->at(fsgrids::efield::EZ));
               CHECK_FLOAT(volGrid0->at(fsgrids::volfields::EZVOL))
            } else {
               volGrid0->at(fsgrids::volfields::EZVOL) = 0.0;
            }
         }
      }
   }
   
   phiprof::stop("Calculate volume averaged fields",N_cells,"Spatial Cells");
}
