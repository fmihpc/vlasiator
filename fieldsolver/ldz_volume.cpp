/*
 * This file is part of Vlasiator.
 * 
 * Copyright 2015 Finnish Meteorological Institute
 * 
 */

#include <cstdlib>

#include "fs_common.h"
#include "ldz_volume.hpp"

#ifndef NDEBUG
   #define DEBUG_FSOLVER
#endif

using namespace std;

void calculateVolumeAveragedFields(
   const FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 3, 2> & perBGrid,
   const FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 3, 2> & EGrid,
   const FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, 3, 2> & dPerBGrid,
   FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, 3, 2> & volGrid,
   const FsGrid< fsgrids::technical, 3, 2> & technicalGrid
) {
   const std::array<int, 3> gridDims = technicalGrid.getLocalSize();
   const size_t N_cells = gridDims[0]*gridDims[1]*gridDims[2];
   phiprof::start("Calculate volume averaged fields");
   
   #pragma omp parallel for collapse(3)
   for (uint k=0; k<gridDims[2]; k++) {
      for (uint j=0; j<gridDims[1]; j++) {
         for (uint i=0; i<gridDims[0]; i++) {
            if(technicalGrid.get(i,j,k)->sysBoundaryFlag == sysboundaries::DO_NOT_COMPUTE) continue;
            
            Real perturbedCoefficients[Rec::N_REC_COEFFICIENTS];
            const std::array<Real, fsgrids::volfields::N_VOL> * volGrid0 = volGrid.get(i,j,k);
            
            // Calculate reconstruction coefficients for this cell:
            reconstructionCoefficients(
               perBGrid,
               dPerBGrid,
               perturbedCoefficients,
               i,
               j,
               k,
               2,
            );
            
            // Calculate volume average of B:
            volGrid0[fsgrids::volfields::PERBXVOL] = perturbedCoefficients[Rec::a_0];
            volGrid0[fsgrids::volfields::PERBYVOL] = perturbedCoefficients[Rec::b_0];
            volGrid0[fsgrids::volfields::PERBZVOL] = perturbedCoefficients[Rec::c_0];

            // Calculate volume average of E (FIXME NEEDS IMPROVEMENT):
            const FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 3, 2>* const EGrid_i1j1k1 = EGrid.get(i,j,k);
            if ( technicalGrid.get(i,j,k)->sysBoundaryFlag == sysboundaries::NOT_SYSBOUNDARY ||
                (technicalGrid.get(i,j,k)->sysBoundaryFlag != sysboundaries::NOT_SYSBOUNDARY && technicalGrid.get(i,j,k)->sysBoundaryLayer == 1)
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
               
               const FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 3, 2>* const EGrid_i1j2k1 = EGrid.get(i  ,j+1,k  );
               const FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 3, 2>* const EGrid_i1j1k2 = EGrid.get(i  ,j  ,k+1);
               const FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 3, 2>* const EGrid_i1j2k2 = EGrid.get(i  ,j+1,k+1);
               
               CHECK_FLOAT(EGrid_i1j1k1[fsgrids::efield::EX])
               CHECK_FLOAT(EGrid_i1j2k1[fsgrids::efield::EX])
               CHECK_FLOAT(EGrid_i1j1k2[fsgrids::efield::EX])
               CHECK_FLOAT(EGrid_i1j2k2[fsgrids::efield::EX])
               volGrid0[fsgrids::volfields::EXVOL] = FOURTH*(EGrid_i1j1k1[fsgrids::efield::EX] + EGrid_i1j2k1[fsgrids::efield::EX] + EGrid_i1j1k2[fsgrids::efield::EX] + EGrid_i1j2k2[fsgrids::efield::EX]);
               CHECK_FLOAT(volGrid0[fsgrids::volfields::EXVOL])
            } else {
               volGrid0[fsgrids::volfields::EXVOL] = 0.0;
            }

            if ( technicalGrid.get(i,j,k)->sysBoundaryFlag == sysboundaries::NOT_SYSBOUNDARY ||
                (technicalGrid.get(i,j,k)->sysBoundaryFlag != sysboundaries::NOT_SYSBOUNDARY && technicalGrid.get(i,j,k)->sysBoundaryLayer == 1)
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
               
               const FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 3, 2>* const EGrid_i2j1k1 = EGrid.get(i+1,j  ,k  );
               const FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 3, 2>* const EGrid_i1j1k2 = EGrid.get(i  ,j  ,k+1);
               const FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 3, 2>* const EGrid_i2j1k2 = EGrid.get(i+1,j  ,k+1);
               
               CHECK_FLOAT(EGrid_i1j1k1[fsgrids::efield::EY])
               CHECK_FLOAT(EGrid_i2j1k1[fsgrids::efield::EY])
               CHECK_FLOAT(EGrid_i1j1k2[fsgrids::efield::EY])
               CHECK_FLOAT(EGrid_i2j1k2[fsgrids::efield::EY])
               volGrid0[fsgrids::volfields::EYVOL] = FOURTH*(EGrid_i1j1k1[fsgrids::efield::EY] + EGrid_i2j1k1[fsgrids::efield::EY] + EGrid_i1j1k2[fsgrids::efield::EY] + EGrid_i2j1k2[fsgrids::efield::EY]);
               CHECK_FLOAT(volGrid0[fsgrids::volfields::EYVOL])
            } else {
               volGrid0[fsgrids::volfields::EYVOL] = 0.0;
            }

            if ( technicalGrid.get(i,j,k)->sysBoundaryFlag == sysboundaries::NOT_SYSBOUNDARY ||
                (technicalGrid.get(i,j,k)->sysBoundaryFlag != sysboundaries::NOT_SYSBOUNDARY && technicalGrid.get(i,j,k)->sysBoundaryLayer == 1)
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
               
               const FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 3, 2>* const EGrid_i2j1k1 = EGrid.get(i+1,j  ,k  );
               const FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 3, 2>* const EGrid_i1j2k1 = EGrid.get(i  ,j+1,k  );
               const FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 3, 2>* const EGrid_i2j2k1 = EGrid.get(i+1,j+1,k  );
               
               CHECK_FLOAT(EGrid_i1j1k1[fsgrids::efield::EZ])
               CHECK_FLOAT(EGrid_i2j1k1[fsgrids::efield::EZ])
               CHECK_FLOAT(EGrid_i1j2k1[fsgrids::efield::EZ])
               CHECK_FLOAT(EGrid_i2j2k1[fsgrids::efield::EZ])
               volGrid0[fsgrids::volfields::EZVOL] = FOURTH*(EGrid_i1j1k1[fsgrids::efield::EZ] + EGrid_i2j1k1[fsgrids::efield::EZ] + EGrid_i1j2k1[fsgrids::efield::EZ] + EGrid_i2j2k1[fsgrids::efield::EZ]);
               CHECK_FLOAT(volGrid0[fsgrids::volfields::EZVOL])
            } else {
               volGrid0[fsgrids::volfields::EZVOL] = 0.0;
            }
         }
      }
   }
   
   phiprof::stop("Calculate volume averaged fields",N_cells,"Spatial Cells");
}
