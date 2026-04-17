/*
 * This file is part of Vlasiator.
 * Copyright 2010-2025 Finnish Meteorological Institute and University of Helsinki
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

/*
   In this file we define functions and variables that are used for both
   CPU and GPU versions of pitchAngleDiffusion
*/

#include "../parameters.h"
#include "../object_wrapper.h"
#include <math.h>
//#include <cmath> // NaN Inf checks
#include <iostream>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <Eigen/Geometry>
#include "common_pitch_angle_diffusion.hpp"

/* Storage of Temperature anisotropy to beta parallel array for pitch-angle diffusion parametrization
   see: Parametrization of coefficients for sub-grid modeling of pitch-angle diffusion in global
   magnetospheric hybrid-Vlasov simulations, M. Dubart, M. Battarbee, U. Ganse, A. Osmane, F. Spanier,
   J. Suni, G. Cozzani, K. Horaites, K. Papadakis, Y. Pfau-Kempf, V. Tarvus, and M. Palmroth,
   Physics of Plasmas 30, 123903 (2023)
   https://doi.org/10.1063/5.0176376
 */
std::vector<Real> betaParaArray;
std::vector<Real> TanisoArray;
std::vector<Real> nu0Array;
size_t n_betaPara = 0;
size_t n_Taniso = 0;
bool nuArrayRead = false;

void readNuArrayFromFile() {
   if (nuArrayRead) {
      return;
   }

   // Read from NU0BOX.DAT (or other file if declared in parameters)
   std::string PATHfile = Parameters::PADnu0;
   std::ifstream FILEDmumu;
   FILEDmumu.open(PATHfile);

   // verify file access was successful
   if (!FILEDmumu.is_open()) {
      std::cerr<<"Error opening file "<<PATHfile<<"!"<<std::endl;
      if (FILEDmumu.fail()) {
         std::cerr<<strerror(errno)<<std::endl;
      }
      abort();
   }

   // Read betaPara strings from file
   std::string lineBeta;
   for (int i = 0; i < 2; i++) {
      std::getline(FILEDmumu,lineBeta);
   }
   std::istringstream issBeta(lineBeta);
   float numBeta;
   while ((issBeta >> numBeta)) { // Stream read from issBeta into numBeta
      betaParaArray.push_back(numBeta);
   }

   // Read Taniso strings from file
   std::string lineTaniso;
   for (int i = 0; i < 2; i++) {
      std::getline(FILEDmumu,lineTaniso);
   }
   std::istringstream issTaniso(lineTaniso);
   float numTaniso;
   while ((issTaniso >> numTaniso)) { // Stream read from issTaniso into numTaniso
      TanisoArray.push_back(numTaniso);
   }

   // Discard one line
   std::string lineDUMP;
   for (int i = 0; i < 1; i++) {
      std::getline(FILEDmumu,lineDUMP);
   }

   // Read values of nu0 from file
   std::string linenu0;
   n_betaPara = betaParaArray.size();
   n_Taniso = TanisoArray.size();
   nu0Array.resize(n_betaPara*n_Taniso);

   for (size_t i = 0; i < n_betaPara; i++) {
      std::getline(FILEDmumu,linenu0);
      std::istringstream issnu0(linenu0);
      std::vector<Real> tempLINE;
      float numTEMP;
      while((issnu0 >> numTEMP)) {
         tempLINE.push_back(numTEMP);
      }
      if (tempLINE.size() != n_Taniso) {
         std::cerr<<"ERROR! line "<<i<<" entry in "<<PATHfile<<" has "<<tempLINE.size()<<" entries instead of expected "<<n_Taniso<<"!"<<std::endl;
         abort();
      }
      for (size_t j = 0; j < n_Taniso; j++) {
         nu0Array[i*n_Taniso+j] = tempLINE[j];
      }
   }

   nuArrayRead = true;
   FILEDmumu.close();
}

/* Linear interpolation of diffusion coefficient from above array
 */
Realf interpolateNuFromArray(
   const Real Taniso_in,
   const Real betaParallel_in
   ) {
   Real Taniso = Taniso_in;
   Real betaParallel = betaParallel_in;
   int betaIndx = -1;
   int TanisoIndx = -1;
   for (size_t i = 0; i < betaParaArray.size(); i++) {
      if (betaParallel >= betaParaArray[i]) {
         betaIndx   = i;
      }
   }
   for (size_t i = 0; i < TanisoArray.size()  ; i++) {
      if (Taniso       >= TanisoArray[i]  ) {
         TanisoIndx = i;
      }
   }

   if ( (betaIndx < 0) || (TanisoIndx < 0) ) {
      // Values below table lower bounds; no diffusion required.
      return 0.0;
   } else {
      // Interpolate values from table; if values are above bounds, cap to maximum value.
      if (betaIndx >= (int)betaParaArray.size()-1) {
         betaIndx = (int)betaParaArray.size()-2; // force last bin
         betaParallel = betaParaArray[betaIndx+1]; // force interpolation to bin top
      }
      if (TanisoIndx >= (int)TanisoArray.size()-1) {
         TanisoIndx = (int)TanisoArray.size()-2; // force last bin
         Taniso = TanisoArray[TanisoIndx+1]; // force interpolation to bin top
      }
      // bi-linear interpolation with weighted mean to find nu0(betaParallel,Taniso)
      const Real beta1   = betaParaArray[betaIndx];
      const Real beta2   = betaParaArray[betaIndx+1];
      const Real Taniso1 = TanisoArray[TanisoIndx];
      const Real Taniso2 = TanisoArray[TanisoIndx+1];
      const Real nu011   = nu0Array[betaIndx*n_Taniso+TanisoIndx];
      const Real nu012   = nu0Array[betaIndx*n_Taniso+TanisoIndx+1];
      const Real nu021   = nu0Array[(betaIndx+1)*n_Taniso+TanisoIndx];
      const Real nu022   = nu0Array[(betaIndx+1)*n_Taniso+TanisoIndx+1];
      // Weights
      const Real w11 = (beta2 - betaParallel)*(Taniso2 - Taniso)  / ( (beta2 - beta1)*(Taniso2-Taniso1) );
      const Real w12 = (beta2 - betaParallel)*(Taniso  - Taniso1) / ( (beta2 - beta1)*(Taniso2-Taniso1) );
      const Real w21 = (betaParallel - beta1)*(Taniso2 - Taniso)  / ( (beta2 - beta1)*(Taniso2-Taniso1) );
      const Real w22 = (betaParallel - beta1)*(Taniso  - Taniso1) / ( (beta2 - beta1)*(Taniso2-Taniso1) );
      // Linear interpolation (with fudge factor divisor)
      return (w11*nu011 + w12*nu012 + w21*nu021 + w22*nu022)/Parameters::PADfudge;
   }
}

void computePitchAngleDiffusionParameters(
   SpatialCell& cell,
   const uint popID, size_t CellIdx, bool& currentSpatialLoopComplete,
   Realf& sparsity, std::array<Real,3>& b, Real& nu0
   ){

   sparsity   = 0.01 * cell.getVelocityBlockMinValue(popID);

   currentSpatialLoopComplete = false;

   // Diffusion coefficient to use in this cell
   nu0 = 0.0;

   // Compute b
   const std::array<Real,3> B = {cell.parameters[CellParams::PERBXVOL] +  cell.parameters[CellParams::BGBXVOL],
      cell.parameters[CellParams::PERBYVOL] +  cell.parameters[CellParams::BGBYVOL],
      cell.parameters[CellParams::PERBZVOL] +  cell.parameters[CellParams::BGBZVOL]};
   const Real Bnorm           = sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);
   b[0] = B[0]/Bnorm;
   b[1] = B[1]/Bnorm;
   b[2] = B[2]/Bnorm;

   if (P::PADcoefficient >= 0) {
      // User-provided single diffusion coefficient
      nu0 = P::PADcoefficient;
   } else {
      // Use nu0 values based on Taniso and betaPara read from file
      if (!nuArrayRead) {
         std::cerr<<" ERROR! Attempting to interpolate nu0 value but file has not been read."<<std::endl;
         abort();
      }

      // Perform Eigen rotation to find parallel and perpendicular pressure
      Eigen::Matrix3d rot = Eigen::Quaterniond::FromTwoVectors(Eigen::Vector3d{b[0], b[1], b[2]}, Eigen::Vector3d{0, 0, 1}).normalized().toRotationMatrix();
      Eigen::Matrix3d Ptensor {
         {cell.parameters[CellParams::P_11], cell.parameters[CellParams::P_12], cell.parameters[CellParams::P_13]},
         {cell.parameters[CellParams::P_12], cell.parameters[CellParams::P_22], cell.parameters[CellParams::P_23]},
         {cell.parameters[CellParams::P_13], cell.parameters[CellParams::P_23], cell.parameters[CellParams::P_33]},
      };
      Eigen::Matrix3d transposerot = rot.transpose();
      Eigen::Matrix3d Pprime = rot * Ptensor * transposerot;

      // Anisotropy
      Real Taniso = 0.0;
      if (Pprime(2, 2) > std::numeric_limits<Real>::min()) {
         Taniso = (Pprime(0, 0) + Pprime(1, 1)) / (2 * Pprime(2, 2));
      }
      // Beta Parallel
      Real betaParallel = 0.0;
      if (Bnorm > 0) {
         betaParallel = 2.0 * physicalconstants::MU_0 * Pprime(2, 2) / (Bnorm*Bnorm);
      }
      // Find anisotropy and beta parallel indexes from read table
      nu0 = interpolateNuFromArray(Taniso,betaParallel);
   }

   // Enable nu0 disk output; skip cells where diffusion is not required (or diffusion coefficient is very small).
   cell.parameters[CellParams::NU0] = nu0;
   if (nu0 <= 0.001) {
      currentSpatialLoopComplete = true;
   }
}
