/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://vlasiator.fmi.fi/
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
#include "fs_cache.h"
#include "ldz_volume.hpp"

#ifndef NDEBUG
   #define DEBUG_FSOLVER
#endif

using namespace std;

void calculateVolumeAveragedFields(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                   std::vector<fs_cache::CellCache>& cache,
                                   const std::vector<uint16_t>& cells) {
   phiprof::start("Calculate volume averaged fields");

   namespace fs = fieldsolver;
   namespace cp = CellParams;

   cuint EX_CELLS = (1 << calcNbrNumber(1,1,1))
                  | (1 << calcNbrNumber(1,2,1))
                  | (1 << calcNbrNumber(1,1,2))
                  | (1 << calcNbrNumber(1,2,2));

   cuint EY_CELLS = (1 << calcNbrNumber(1,1,1))
                  | (1 << calcNbrNumber(2,1,1))
                  | (1 << calcNbrNumber(1,1,2))
                  | (1 << calcNbrNumber(2,1,2));

   cuint EZ_CELLS = (1 << calcNbrNumber(1,1,1))
                  | (1 << calcNbrNumber(2,1,1))
                  | (1 << calcNbrNumber(1,2,1))
                  | (1 << calcNbrNumber(2,2,1));

   // NOTE: cache does not include DO_NOT_COMPUTE cells

   #pragma omp parallel for
   for (size_t c=0; c<cells.size(); ++c) {
      const uint16_t localID = cells[c];

      Real perturbedCoefficients[Rec::N_REC_COEFFICIENTS];
      cuint existingCells = cache[localID].existingCellsFlags;

      // Calculate reconstruction coefficients for this cell:
      reconstructionCoefficients(cache[localID],perturbedCoefficients,2,RK_ORDER1);

      // Calculate volume average of B:
      Real* const cellParams   = cache[localID].cells[fs_cache::calculateNbrID(1  ,1  ,1  )]->parameters;
      //Real* const cellParams   = cache[localID].parameters[fs_cache::C222];
      cellParams[cp::PERBXVOL] = perturbedCoefficients[Rec::a_0];
      cellParams[cp::PERBYVOL] = perturbedCoefficients[Rec::b_0];
      cellParams[cp::PERBZVOL] = perturbedCoefficients[Rec::c_0];

      // Calculate volume average of E (FIXME NEEDS IMPROVEMENT):
      creal* const cep_i1j1k1 = cellParams;

      if ((existingCells & EX_CELLS) == EX_CELLS) {
         #ifdef DEBUG_FSOLVER
         bool ok = true;
         if (cache[localID].cells[fs_cache::calculateNbrID(1  ,1+1,1  )] == NULL) ok = false;
         if (cache[localID].cells[fs_cache::calculateNbrID(1  ,1  ,1+1)] == NULL) ok = false;
         if (cache[localID].cells[fs_cache::calculateNbrID(1  ,1+1,1+1)] == NULL) ok = false;
         if (ok == false) {
            stringstream ss;
            ss << "ERROR, got NULL neighbor in " << __FILE__ << ":" << __LINE__ << endl;
            cerr << ss.str(); exit(1);
         }
         #endif
         
         creal* const cep_i1j2k1 = cache[localID].cells[fs_cache::calculateNbrID(1  ,1+1,1  )]->parameters;
         creal* const cep_i1j1k2 = cache[localID].cells[fs_cache::calculateNbrID(1  ,1  ,1+1)]->parameters;
         creal* const cep_i1j2k2 = cache[localID].cells[fs_cache::calculateNbrID(1  ,1+1,1+1)]->parameters;
         //creal* const cep_i1j2k1 = cache[localID].parameters[fs_cache::C232];
         //creal* const cep_i1j1k2 = cache[localID].parameters[fs_cache::C223];
         //creal* const cep_i1j2k2 = cache[localID].parameters[fs_cache::C233];
         
         CHECK_FLOAT(cep_i1j1k1[cp::EX])
         CHECK_FLOAT(cep_i1j2k1[cp::EX])
         CHECK_FLOAT(cep_i1j1k2[cp::EX])
         CHECK_FLOAT(cep_i1j2k2[cp::EX])
         cellParams[cp::EXVOL] = FOURTH*(cep_i1j1k1[cp::EX] + cep_i1j2k1[cp::EX] + cep_i1j1k2[cp::EX] + cep_i1j2k2[cp::EX]);
         CHECK_FLOAT(cellParams[cp::EXVOL])
      } else {
         cellParams[cp::EXVOL] = 0.0;
      }

      if ((existingCells & EY_CELLS) == EY_CELLS) {
         #ifdef DEBUG_FSOLVER
         bool ok = true;
         if (cache[localID].cells[fs_cache::calculateNbrID(1+1,1  ,1  )] == NULL) ok = false;
         if (cache[localID].cells[fs_cache::calculateNbrID(1  ,1  ,1+1)] == NULL) ok = false;
         if (cache[localID].cells[fs_cache::calculateNbrID(1+1,1  ,1+1)] == NULL) ok = false;
         if (ok == false) {
            stringstream ss;
            ss << "ERROR, got NULL neighbor in " << __FILE__ << ":" << __LINE__ << endl;
            cerr << ss.str(); exit(1);
         }
         #endif
         
         creal* const cep_i2j1k1 = cache[localID].cells[fs_cache::calculateNbrID(1+1,1  ,1  )]->parameters;
         creal* const cep_i1j1k2 = cache[localID].cells[fs_cache::calculateNbrID(1  ,1  ,1+1)]->parameters;
         creal* const cep_i2j1k2 = cache[localID].cells[fs_cache::calculateNbrID(1+1,1  ,1+1)]->parameters;
         //creal* const cep_i2j1k1 = cache[localID].parameters[fs_cache::C322];
         //creal* const cep_i1j1k2 = cache[localID].parameters[fs_cache::C223];
         //creal* const cep_i2j1k2 = cache[localID].parameters[fs_cache::C323];
         
         CHECK_FLOAT(cep_i1j1k1[cp::EY])
         CHECK_FLOAT(cep_i2j1k1[cp::EY])
         CHECK_FLOAT(cep_i1j1k2[cp::EY])
         CHECK_FLOAT(cep_i2j1k2[cp::EY])
         cellParams[cp::EYVOL] = FOURTH*(cep_i1j1k1[cp::EY] + cep_i2j1k1[cp::EY] + cep_i1j1k2[cp::EY] + cep_i2j1k2[cp::EY]);
         CHECK_FLOAT(cellParams[cp::EYVOL])
      } else {
         cellParams[cp::EYVOL] = 0.0;
      }

      if ((existingCells & EZ_CELLS) == EZ_CELLS) {
         #ifdef DEBUG_FSOLVER
         bool ok = true;
         if (cache[localID].cells[fs_cache::calculateNbrID(1+1,1  ,1  )] == NULL) ok = false;
         if (cache[localID].cells[fs_cache::calculateNbrID(1  ,1+1,1  )] == NULL) ok = false;
         if (cache[localID].cells[fs_cache::calculateNbrID(1+1,1+1,1  )] == NULL) ok = false;
         if (ok == false) {
            stringstream ss;
            ss << "ERROR, got NULL neighbor in " << __FILE__ << ":" << __LINE__ << endl;
            cerr << ss.str(); exit(1);
         }
         #endif
         
         creal* const cep_i2j1k1 = cache[localID].cells[fs_cache::calculateNbrID(1+1,1  ,1  )]->parameters;
         creal* const cep_i1j2k1 = cache[localID].cells[fs_cache::calculateNbrID(1  ,1+1,1  )]->parameters;
         creal* const cep_i2j2k1 = cache[localID].cells[fs_cache::calculateNbrID(1+1,1+1,1  )]->parameters;
         //creal* const cep_i2j1k1 = cache[localID].parameters[fs_cache::C322];
         //creal* const cep_i1j2k1 = cache[localID].parameters[fs_cache::C232];
         //creal* const cep_i2j2k1 = cache[localID].parameters[fs_cache::C332];

         CHECK_FLOAT(cep_i1j1k1[cp::EZ])
         CHECK_FLOAT(cep_i2j1k1[cp::EZ])
         CHECK_FLOAT(cep_i1j2k1[cp::EZ])
         CHECK_FLOAT(cep_i2j2k1[cp::EZ])
         cellParams[cp::EZVOL] = FOURTH*(cep_i1j1k1[cp::EZ] + cep_i2j1k1[cp::EZ] + cep_i1j2k1[cp::EZ] + cep_i2j2k1[cp::EZ]);
         CHECK_FLOAT(cellParams[cp::EZVOL])
      } else {
         cellParams[cp::EZVOL] = 0.0;
      }
   }

   phiprof::stop("Calculate volume averaged fields",cache.size(),"Spatial Cells");
}
         
         
         



