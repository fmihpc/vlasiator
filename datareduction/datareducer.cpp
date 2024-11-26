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
#include <iostream>

#include "datareducer.h"
#include "../common.h"
#include "dro_populations.h"
#include "../sysboundary/ionosphere.h"
#include "../fieldtracing/fieldtracing.h"

using namespace std;

void initializeDataReducers(DataReducer * outputReducer, DataReducer * diagnosticReducer)
{
   typedef Parameters P;

   vector<string>::const_iterator it;
   for (it = P::outputVariableList.begin();
        it != P::outputVariableList.end();
        it++) {

      /* Note: Each data reducer generation should be followed by a call to setUnitMetaData
         with the following arguments:
         unit, unit in LaTeX formulation, variable in LaTeX formulation, conversion factor
      */

      // Sidestep mixed case errors
      std::string lowercase = *it;
      for(auto& c : lowercase) c = tolower(c);

      if(P::systemWriteAllDROs || lowercase == "fg_b" || lowercase == "b") { // Bulk magnetic field at Yee-Lattice locations
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_b",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]*3);

               // Iterate through fsgrid cells and extract total magnetic field
               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[3*(gridSize[1]*gridSize[0]*z + gridSize[0]*y + x)] =     (*BgBGrid.get(x,y,z))[fsgrids::BGBX]
                           + (*perBGrid.get(x,y,z))[fsgrids::PERBX];
                        retval[3*(gridSize[1]*gridSize[0]*z + gridSize[0]*y + x) + 1] = (*BgBGrid.get(x,y,z))[fsgrids::BGBY]
                           + (*perBGrid.get(x,y,z))[fsgrids::PERBY];
                        retval[3*(gridSize[1]*gridSize[0]*z + gridSize[0]*y + x) + 2] = (*BgBGrid.get(x,y,z))[fsgrids::BGBZ]
                           + (*perBGrid.get(x,y,z))[fsgrids::PERBZ];
                     }
                  }
               }
               return retval;
            }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"T","$\\mathrm{T}$","$B_\\mathrm{fg}$","1.0");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "fg_backgroundb" || lowercase == "backgroundb" || lowercase == "fg_b_background") { // Static (typically dipole) magnetic field part
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_b_background",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]*3);

               // Iterate through fsgrid cells and extract background B
               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[3*(gridSize[1]*gridSize[0]*z + gridSize[0]*y + x)] =     (*BgBGrid.get(x,y,z))[fsgrids::BGBX];
                        retval[3*(gridSize[1]*gridSize[0]*z + gridSize[0]*y + x) + 1] = (*BgBGrid.get(x,y,z))[fsgrids::BGBY];
                        retval[3*(gridSize[1]*gridSize[0]*z + gridSize[0]*y + x) + 2] = (*BgBGrid.get(x,y,z))[fsgrids::BGBZ];
                     }
                  }
               }
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size() - 1, "T", "$\\mathrm{T}$", "$B_\\mathrm{bg,fg}$", "1.0");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "fg_backgroundbvol" || lowercase == "backgroundbvol" || lowercase == "fg_b_background_vol") { // Static (typically dipole) magnetic field part, volume-averaged
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_b_background_vol",[](
            FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
            FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
            FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
            FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
            FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
            FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
            FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
            FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
            FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
            FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]*3);

               // Iterate through fsgrid cells and extract total BVOL
               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[3*(gridSize[1]*gridSize[0]*z + gridSize[0]*y + x)] =     (*BgBGrid.get(x,y,z))[fsgrids::BGBXVOL];
                        retval[3*(gridSize[1]*gridSize[0]*z + gridSize[0]*y + x) + 1] = (*BgBGrid.get(x,y,z))[fsgrids::BGBYVOL];
                        retval[3*(gridSize[1]*gridSize[0]*z + gridSize[0]*y + x) + 2] = (*BgBGrid.get(x,y,z))[fsgrids::BGBZVOL];
                     }
                  }
               }
               return retval;
            }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"T","$\\mathrm{T}$","$B_\\mathrm{bg,vol,fg}$","1.0");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }

      if(P::systemWriteAllDROs || lowercase == "fg_perturbedb" || lowercase == "perturbedb" || lowercase == "fg_b_perturbed") { // Fluctuating magnetic field part
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_b_perturbed",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]*3);

               // Iterate through fsgrid cells and extract values
               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[3*(gridSize[1]*gridSize[0]*z + gridSize[0]*y + x)] =     (*perBGrid.get(x,y,z))[fsgrids::PERBX];
                        retval[3*(gridSize[1]*gridSize[0]*z + gridSize[0]*y + x) + 1] = (*perBGrid.get(x,y,z))[fsgrids::PERBY];
                        retval[3*(gridSize[1]*gridSize[0]*z + gridSize[0]*y + x) + 2] = (*perBGrid.get(x,y,z))[fsgrids::PERBZ];
                     }
                  }
               }
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"T","$\\mathrm{T}$","$B_\\mathrm{per,fg}$","1.0");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "fg_e" || lowercase == "e") { // Bulk electric field at Yee-lattice locations
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_e",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]*3);

               // Iterate through fsgrid cells and extract E values
               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[3*(gridSize[1]*gridSize[0]*z + gridSize[0]*y + x)] =     (*EGrid.get(x,y,z))[fsgrids::EX];
                        retval[3*(gridSize[1]*gridSize[0]*z + gridSize[0]*y + x) + 1] = (*EGrid.get(x,y,z))[fsgrids::EY];
                        retval[3*(gridSize[1]*gridSize[0]*z + gridSize[0]*y + x) + 2] = (*EGrid.get(x,y,z))[fsgrids::EZ];
                     }
                  }
               }
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"V/m","$\\mathrm{V}\\,\\mathrm{m}^{-1}$","$E$","1.0");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "vg_rhom" || lowercase == "rhom") { // Overall mass density (summed over all populations)
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("vg_rhom",CellParams::RHOM,1));
         outputReducer->addMetadata(outputReducer->size()-1,"kg/m^3","$\\mathrm{kg}\\,\\mathrm{m}^{-3}$","$\\rho_\\mathrm{m}$","1.0");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "vg_drift") { // Nudge velocity drift near ionosphere
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("vg_drift",CellParams::BULKV_FORCING_X,3));
         outputReducer->addMetadata(outputReducer->size()-1,"m/s","$\\mathrm{m}\\,\\mathrm{s}^{-1}$","$V$","1.0");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "fg_rhom") { // Overall mass density (summed over all populations)
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_rhom",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               // Iterate through fsgrid cells and extract rho valuesg
               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = (*momentsGrid.get(x,y,z))[fsgrids::RHOM];
                     }
                  }
               }
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"kg/m^3","$\\mathrm{kg}\\,\\mathrm{m}^{-3}$","$\\rho_\\mathrm{m}$","1.0");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "vg_rhoq" || lowercase == "rhoq") { // Overall charge density (summed over all populations)
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("vg_rhoq",CellParams::RHOQ,1));
         outputReducer->addMetadata(outputReducer->size()-1,"C/m^3","$\\mathrm{C}\\,\\mathrm{m}^{-3}$","$\\rho_\\mathrm{q}$","1.0");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "fg_rhoq") { // Overall charge density (summed over all populations)
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_rhoq",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               // Iterate through fsgrid cells and extract charge density
               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = (*momentsGrid.get(x,y,z))[fsgrids::RHOQ];
                     }
                  }
               }
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"C/m^3","$\\mathrm{C}\\,\\mathrm{m}^{-3}$","$\\rho_\\mathrm{q}$","1.0");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "populations_rho" || lowercase == "populations_vg_rho") { // Per-population particle number density
         for(unsigned int i =0; i < getObjectWrapper().particleSpecies.size(); i++) {
            species::Species& species=getObjectWrapper().particleSpecies[i];
            const std::string& pop = species.name;
            outputReducer->addOperator(new DRO::DataReductionOperatorPopulations<Real>(pop + "/vg_rho", i, offsetof(spatial_cell::Population, RHO), 1));
            outputReducer->addMetadata(outputReducer->size()-1,"1/m^3","$\\mathrm{m}^{-3}$","$n_\\mathrm{"+pop+"}$","1.0");
         }
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }

      if(P::systemWriteAllDROs || lowercase == "v" || lowercase == "vg_v") { // Overall effective bulk density defining the center-of-mass frame from all populations
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("vg_v",CellParams::VX,3));
         outputReducer->addMetadata(outputReducer->size()-1,"m/s","$\\mathrm{m}\\,\\mathrm{s}^{-1}$","$V$","1.0");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "fg_v") { // Overall effective bulk density defining the center-of-mass frame from all populations
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_v",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]*3);

               // Iterate through fsgrid cells and extract bulk Velocity
               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[3*(gridSize[1]*gridSize[0]*z + gridSize[0]*y + x)] =     (*momentsGrid.get(x,y,z))[fsgrids::VX];
                        retval[3*(gridSize[1]*gridSize[0]*z + gridSize[0]*y + x) + 1] = (*momentsGrid.get(x,y,z))[fsgrids::VY];
                        retval[3*(gridSize[1]*gridSize[0]*z + gridSize[0]*y + x) + 2] = (*momentsGrid.get(x,y,z))[fsgrids::VZ];
                     }
                  }
               }
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"m/s","$\\mathrm{m}\\,\\mathrm{s}^{-1}$","$V$","1.0");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "populations_v" || lowercase == "populations_vg_v") { // Per population bulk velocities
         for(unsigned int i =0; i < getObjectWrapper().particleSpecies.size(); i++) {
            species::Species& species=getObjectWrapper().particleSpecies[i];
            const std::string& pop = species.name;
            outputReducer->addOperator(new DRO::DataReductionOperatorPopulations<Real>(pop + "/vg_v", i, offsetof(spatial_cell::Population, V), 3));
            outputReducer->addMetadata(outputReducer->size()-1,"m/s","$\\mathrm{m}\\,\\mathrm{s}^{-1}$","$V_\\mathrm{"+pop+"}$","1.0");
         }
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "populations_moments_backstream" || lowercase == "populations_moments_nonthermal" || lowercase == "populations_vg_moments_nonthermal") { // Per-population moments of the backstreaming part
         for(unsigned int i =0; i < getObjectWrapper().particleSpecies.size(); i++) {
            species::Species& species=getObjectWrapper().particleSpecies[i];
            const std::string& pop = species.name;
            outputReducer->addOperator(new DRO::VariableRhoNonthermal(i));
            outputReducer->addOperator(new DRO::VariableVNonthermal(i));
            outputReducer->addOperator(new DRO::VariablePTensorNonthermalDiagonal(i));
            outputReducer->addOperator(new DRO::VariablePTensorNonthermalOffDiagonal(i));
            outputReducer->addMetadata(outputReducer->size()-4,"1/m^3","$\\mathrm{m}^{-3}$","$n_\\mathrm{"+pop+",nt}$","1.0");
            outputReducer->addMetadata(outputReducer->size()-3,"m/s","$\\mathrm{m}\\,\\mathrm{s}^{-1}$","$V_\\mathrm{"+pop+",nt}$","1.0");
            outputReducer->addMetadata(outputReducer->size()-2,"Pa","$\\mathrm{Pa}$","$\\mathcal{P}_\\mathrm{"+pop+",nt}$","1.0");
            outputReducer->addMetadata(outputReducer->size()-1,"Pa","$\\mathrm{Pa}$","$\\mathcal{\\tilde{P}}_\\mathrm{"+pop+",nt}$","1.0");
         }
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "populations_moments_nonbackstream" || lowercase == "populations_moments_thermal" || lowercase == "populations_vg_moments_thermal") { // Per-population moments of the non-backstreaming (thermal?) part.
         for(unsigned int i =0; i < getObjectWrapper().particleSpecies.size(); i++) {
            species::Species& species=getObjectWrapper().particleSpecies[i];
            const std::string& pop = species.name;
            outputReducer->addOperator(new DRO::VariableRhoThermal(i));
            outputReducer->addOperator(new DRO::VariableVThermal(i));
            outputReducer->addOperator(new DRO::VariablePTensorThermalDiagonal(i));
            outputReducer->addOperator(new DRO::VariablePTensorThermalOffDiagonal(i));
            outputReducer->addMetadata(outputReducer->size()-4,"1/m^3","$\\mathrm{m}^{-3}$","$n_\\mathrm{"+pop+",th}$","1.0");
            outputReducer->addMetadata(outputReducer->size()-3,"m/s","$\\mathrm{m}\\,\\mathrm{s}^{-1}$","$V_\\mathrm{"+pop+",th}$","1.0");
            outputReducer->addMetadata(outputReducer->size()-2,"Pa","$\\mathrm{Pa}$","$\\mathcal{P}_\\mathrm{"+pop+",th}$","1.0");
            outputReducer->addMetadata(outputReducer->size()-1,"Pa","$\\mathrm{Pa}$","$\\mathcal{\\tilde{P}}_\\mathrm{"+pop+",th}$","1.0");
         }
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "populations_minvalue" || lowercase == "populations_effectivesparsitythreshold" || lowercase == "populations_vg_effectivesparsitythreshold") {
         // Effective sparsity threshold affecting each cell, if dynamic threshould algorithm is used
         for(unsigned int i =0; i < getObjectWrapper().particleSpecies.size(); i++) {
            species::Species& species=getObjectWrapper().particleSpecies[i];
            const std::string& pop = species.name;
            outputReducer->addOperator(new DRO::VariableEffectiveSparsityThreshold(i));
            outputReducer->addMetadata(outputReducer->size()-1,"s^3/m^6","$\\mathrm{m}^{-6}\\,\\mathrm{s}^{3}$","$f_\\mathrm{"+pop+",min}$","1.0");
         }
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "populations_rholossadjust" || lowercase == "populations_rho_loss_adjust" || lowercase == "populations_vg_rho_loss_adjust") {
         // Accumulated lost particle number, per population, in each cell, since last restart
         for(unsigned int i =0; i < getObjectWrapper().particleSpecies.size(); i++) {
            species::Species& species=getObjectWrapper().particleSpecies[i];
            const std::string& pop = species.name;
            outputReducer->addOperator(new DRO::DataReductionOperatorPopulations<Real>(pop + "/vg_rho_loss_adjust", i, offsetof(spatial_cell::Population, RHOLOSSADJUST), 1));
            outputReducer->addMetadata(outputReducer->size()-1,"1/m^3","$\\mathrm{m}^{-3}$","$\\Delta_\\mathrm{loss} n_\\mathrm{"+pop+"}$","1.0");
         }
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "lbweight" || lowercase == "vg_lbweight" || lowercase == "vg_loadbalanceweight" || lowercase == "vg_loadbalance_weight") {
         // Load balance metric for LB debugging
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("vg_loadbalance_weight",CellParams::LBWEIGHTCOUNTER,1));
         outputReducer->addMetadata(outputReducer->size()-1,"","","$\\mathrm{LB weight}$","");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "maxvdt" || lowercase == "vg_maxdt_acceleration") {
         // Overall maximum timestep constraint as calculated by the velocity space vlasov update
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("vg_maxdt_acceleration",CellParams::MAXVDT,1));
         outputReducer->addMetadata(outputReducer->size()-1,"s","$\\mathrm{s}$","$\\Delta t_\\mathrm{V,max}$","1.0");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "populations_maxvdt" || lowercase == "populations_vg_maxdt_acceleration" || lowercase == "populations_maxdt_acceleration") {
         // Per-population maximum timestep constraint as calculated by the velocity space vlasov update
         for(unsigned int i =0; i < getObjectWrapper().particleSpecies.size(); i++) {
            species::Species& species=getObjectWrapper().particleSpecies[i];
            const std::string& pop = species.name;
            outputReducer->addOperator(new DRO::DataReductionOperatorPopulations<Real>(pop + "/vg_maxdt_acceleration", i, offsetof(spatial_cell::Population, max_dt[1]), 1));
            outputReducer->addMetadata(outputReducer->size()-1,"s","$\\mathrm{s}$","$\\Delta t_\\mathrm{"+pop+",V,max}$","1.0");
         }
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "maxrdt" || lowercase == "vg_maxdt_translation") {
         // Overall maximum timestep constraint as calculated by the real space vlasov update
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("vg_maxdt_translation",CellParams::MAXRDT,1));
         outputReducer->addMetadata(outputReducer->size()-1,"s","$\\mathrm{s}$","$\\Delta t_\\mathrm{R,max}$","1.0");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "populations_maxrdt" || lowercase == "populations_vg_maxdt_translation" || lowercase == "populations_maxdt_translation") {
         // Per-population maximum timestep constraint as calculated by the real space vlasov update
         for(unsigned int i =0; i < getObjectWrapper().particleSpecies.size(); i++) {
            species::Species& species=getObjectWrapper().particleSpecies[i];
            const std::string& pop = species.name;
            outputReducer->addOperator(new DRO::DataReductionOperatorPopulations<Real>(pop + "/vg_maxdt_translation", i, offsetof(spatial_cell::Population, max_dt[0]), 1));
            outputReducer->addMetadata(outputReducer->size()-1,"s","$\\mathrm{s}$","$\\Delta t_\\mathrm{"+pop+",R,max}$","1.0");
         }
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "populations_energydensity" || lowercase == "populations_vg_energydensity") {
         // Per-population energy density in three energy ranges
         for(unsigned int i =0; i < getObjectWrapper().particleSpecies.size(); i++) {
            species::Species& species=getObjectWrapper().particleSpecies[i];
            const std::string& pop = species.name;
            outputReducer->addOperator(new DRO::VariableEnergyDensity(i));
            std::stringstream conversion;
            conversion << (1.0e-6)/physicalconstants::CHARGE;
            outputReducer->addMetadata(outputReducer->size()-1,"eV/cm^3","$\\mathrm{eV}\\,\\mathrm{cm}^{-3}$","$U_\\mathrm{"+pop+"}$",conversion.str());
         }
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "populations_precipitationflux" || lowercase == "populations_vg_precipitationdifferentialflux" || lowercase == "populations_precipitationdifferentialflux") {
         // Per-population precipitation differential flux (within loss cone)
         for(unsigned int i =0; i < getObjectWrapper().particleSpecies.size(); i++) {
            species::Species& species=getObjectWrapper().particleSpecies[i];
            const std::string& pop = species.name;
            outputReducer->addOperator(new DRO::VariablePrecipitationDiffFlux(i));
            std::stringstream conversion;
            conversion << (1.0e-4)*physicalconstants::CHARGE;
            outputReducer->addMetadata(outputReducer->size()-1,"1/(cm^2 sr s eV)","$\\mathrm{cm}^{-2}\\,\\mathrm{sr}^{-1}\\,\\mathrm{s}^{-1}\\,\\mathrm{eV}^{-1}$","$\\mathcal{F}_\\mathrm{"+pop+"}$",conversion.str());
         }
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "populations_precipitationlineflux" || lowercase == "populations_vg_precipitationlinedifferentialflux" || lowercase == "populations_precipitationlinedifferentialflux") {
         // Per-population precipitation differential flux (along line)
         for(unsigned int i =0; i < getObjectWrapper().particleSpecies.size(); i++) {
            species::Species& species=getObjectWrapper().particleSpecies[i];
            const std::string& pop = species.name;
            outputReducer->addOperator(new DRO::VariablePrecipitationLineDiffFlux(i));
            std::stringstream conversion;
            conversion << (1.0e-4)*physicalconstants::CHARGE;
            outputReducer->addMetadata(outputReducer->size()-1,"1/(cm^2 sr s eV)","$\\mathrm{cm}^{-2}\\,\\mathrm{sr}^{-1}\\,\\mathrm{s}^{-1}\\,\\mathrm{eV}^{-1}$","$\\mathcal{F}_\\mathrm{"+pop+"}$",conversion.str());
         }
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "populations_heatflux" || lowercase == "populations_vg_heatflux") {
         // Per-population heat flux vector
         for(unsigned int i =0; i < getObjectWrapper().particleSpecies.size(); i++) {
            species::Species& species=getObjectWrapper().particleSpecies[i];
            const std::string& pop = species.name;
            outputReducer->addOperator(new DRO::VariableHeatFluxVector(i));
            outputReducer->addMetadata(outputReducer->size()-1,"W/m^2","$\\mathrm{W}\\,\\mathrm{m}^{-2}$","$q_\\mathrm{"+pop+"}$","1.0");
         }
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if (P::systemWriteAllDROs || lowercase == "populations_nonmaxwellianity" || lowercase == "populations_vg_nonmaxwellianity") {
         // Per-population dimensionless non-maxwellianity parameter
         for (unsigned int i = 0; i < getObjectWrapper().particleSpecies.size(); i++) {
            species::Species& species = getObjectWrapper().particleSpecies[i];
            const std::string& pop = species.name;
            outputReducer->addOperator(new DRO::VariableNonMaxwellianity(i));
            outputReducer->addMetadata(outputReducer->size() - 1, "", "",
                                       "$\\tilde{\\epsilon}_\\mathrm{M," + pop + "}$", "1.0");
         }
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "maxfieldsdt" || lowercase == "fg_maxfieldsdt" || lowercase == "fg_maxdt_fieldsolver") {
         // Maximum timestep constraint as calculated by the fieldsolver
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_maxdt_fieldsolver",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               // Iterate through fsgrid cells and extract field solver timestep limit
               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = technicalGrid.get(x,y,z)->maxFsDt;
                     }
                  }
               }
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"s","$\\mathrm{s}$","$\\Delta t_\\mathrm{f,max}$","1.0");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "mpirank" || lowercase == "vg_rank") {
         // Map of spatial decomposition of the DCCRG grid into MPI ranks
         outputReducer->addOperator(new DRO::MPIrank);
         outputReducer->addMetadata(outputReducer->size()-1,"","","$\\mathrm{MPI rank}$","");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "fsgridrank" || lowercase == "fg_rank") {
         // Map of spatial decomposition of the FsGrid into MPI ranks
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_rank",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2],technicalGrid.getRank());
               return retval;
             }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"","","$\\mathrm{fGrid rank}$","");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "fg_amr_level") {
         // Map of spatial decomposition of the FsGrid into MPI ranks
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_amr_level",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH>& perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH>& EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH>& EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH>& EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH>& momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH>& dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH>& dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH>& BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH>& volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH>& technicalGrid)->std::vector<double> {


               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               // Iterate through fsgrid cells and extract corresponding AMR level
               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = technicalGrid.get(x,y,z)->refLevel;
                     }
                  }
               }
               return retval;
            }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"","","$\\mathrm{fGrid rank}$","");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "boundarytype" || lowercase == "vg_boundarytype") {
         // Type of boundarycells
         outputReducer->addOperator(new DRO::BoundaryType);
         outputReducer->addMetadata(outputReducer->size()-1,"","","$\\mathrm{vGrid Boundary type}$","");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "fsgridboundarytype" || lowercase == "fg_boundarytype") {
         // Type of boundarycells as stored in FSGrid
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_boundarytype",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               // Iterate through fsgrid cells and extract boundary flag
               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = technicalGrid.get(x,y,z)->sysBoundaryFlag;
                     }
                  }
               }
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"","","$\\mathrm{fGrid Boundary type}$","");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "boundarylayer" || lowercase == "vg_boundarylayer") {
         // For boundaries with multiple layers: layer count per cell
         outputReducer->addOperator(new DRO::BoundaryLayer);
         outputReducer->addMetadata(outputReducer->size()-1,"","","$\\mathrm{vGrid Boundary layer}$","");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "fsgridboundarylayer" || lowercase == "fg_boundarylayer") {
         // Type of boundarycells as stored in FSGrid
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_boundarylayer",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               // Iterate through fsgrid cells and extract boundary layer
               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = technicalGrid.get(x,y,z)->sysBoundaryLayer;
                     }
                  }
               }
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"","","$\\mathrm{fGrid Boundary layer}$","");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "populations_blocks" || lowercase == "populations_vg_blocks") {
         // Per-population velocity space block counts
         for(unsigned int i =0; i < getObjectWrapper().particleSpecies.size(); i++) {
            species::Species& species=getObjectWrapper().particleSpecies[i];
            const std::string& pop = species.name;
            outputReducer->addOperator(new DRO::Blocks(i));
            outputReducer->addMetadata(outputReducer->size()-1,"","","$\\mathrm{"+pop+" blocks}$","");
         }
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "fsaved" || lowercase == "vg_fsaved" || lowercase == "vg_f_saved") {
         // Boolean marker whether a velocity space is saved in a given spatial cell
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("vg_f_saved",CellParams::ISCELLSAVINGF,1));
         outputReducer->addMetadata(outputReducer->size()-1,"","","$f(v)_\\mathrm{saved}$","");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "populations_accsubcycles" || lowercase == "populations_acceleration_subcycles" || lowercase == "populations_vg_acceleration_subcycles") {
         // Per-population number of subcycles performed for velocity space update
         for(unsigned int i =0; i < getObjectWrapper().particleSpecies.size(); i++) {
            species::Species& species=getObjectWrapper().particleSpecies[i];
            const std::string& pop = species.name;
            outputReducer->addOperator(new DRO::DataReductionOperatorPopulations<uint>(pop + "/vg_acceleration_subcycles", i, offsetof(spatial_cell::Population, ACCSUBCYCLES), 1));
            outputReducer->addMetadata(outputReducer->size()-1,"","","$\\mathrm{"+pop+" Acc subcycles}$","");
         }
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "vole" || lowercase == "vg_vole" || lowercase == "evol" || lowercase == "vg_e_vol" || lowercase == "e_vol") {
         // Volume-averaged E field
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("vg_e_vol",CellParams::EXVOL,3));
         outputReducer->addMetadata(outputReducer->size()-1,"V/m","$\\mathrm{V}\\,\\mathrm{m}^{-1}$","$E_\\mathrm{vol,vg}$","1.0");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "fg_vole" || lowercase == "fg_e_vol" || lowercase == "fg_evol") {
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_e_vol",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]*3);

               // Iterate through fsgrid cells and extract EVOL
               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[3*(gridSize[1]*gridSize[0]*z + gridSize[0]*y + x)] =  (*volGrid.get(x,y,z))[fsgrids::volfields::EXVOL];
                        retval[3*(gridSize[1]*gridSize[0]*z + gridSize[0]*y + x) + 1] = (*volGrid.get(x,y,z))[fsgrids::volfields::EYVOL];
                        retval[3*(gridSize[1]*gridSize[0]*z + gridSize[0]*y + x) + 2] = (*volGrid.get(x,y,z))[fsgrids::volfields::EZVOL];
                     }
                  }
               }
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"V/m","$\\mathrm{V}\\,\\mathrm{m}^{-1}$","$E_\\mathrm{vol,fg}$","1.0");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "halle" || lowercase == "fg_halle" || lowercase == "fg_e_hall") {
         for(int index=0; index<fsgrids::N_EHALL; index++) {
            std::string reducer_name = "fg_e_hall_" + std::to_string(index);
            outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid(reducer_name,[index](
                         FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                         FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                         FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                         FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                         FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                         FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                         FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                         FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                         FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                         FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

                  std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
                  std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

                  // Iterate through fsgrid cells and extract EHall
                  for(int z=0; z<gridSize[2]; z++) {
                     for(int y=0; y<gridSize[1]; y++) {
                        for(int x=0; x<gridSize[0]; x++) {
                           retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = (*EHallGrid.get(x,y,z))[index];
                        }
                     }
                  }
                  return retval;
            }
            ));
            outputReducer->addMetadata(outputReducer->size()-1,"V/m","$\\mathrm{V}\\,\\mathrm{m}^{-1}$","$E_\\mathrm{Hall,"+std::to_string(index)+"}$","1.0");
         }
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase =="gradpee" || lowercase == "e_gradpe" || lowercase == "vg_e_gradpe") {
         // Electron pressure gradient contribution to the generalized ohm's law
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("vg_e_gradpe",CellParams::EXGRADPE,3));
         outputReducer->addMetadata(outputReducer->size()-1,"V/m","$\\mathrm{V}\\,\\mathrm{m}^{-1}$","$E_{\\nabla P_\\mathrm{e}}$","1.0");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "volb" || lowercase == "vg_volb" || lowercase == "b_vol" || lowercase == "bvol" || lowercase == "vg_bvol" || lowercase == "vg_b_vol") {
         // Volume-averaged magnetic field
         outputReducer->addOperator(new DRO::VariableBVol);
         outputReducer->addMetadata(outputReducer->size()-1,"T","$\\mathrm{T}$","$B_\\mathrm{vol,vg}$","1.0");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "fg_volb" || lowercase == "fg_bvol" || lowercase == "fg_b_vol") { // Static (typically dipole) magnetic field part
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_b_vol",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]*3);

               // Iterate through fsgrid cells and extract total BVOL
               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[3*(gridSize[1]*gridSize[0]*z + gridSize[0]*y + x)] =     (*BgBGrid.get(x,y,z))[fsgrids::BGBXVOL]
                           + (*volGrid.get(x,y,z))[fsgrids::PERBXVOL];
                        retval[3*(gridSize[1]*gridSize[0]*z + gridSize[0]*y + x) + 1] = (*BgBGrid.get(x,y,z))[fsgrids::BGBYVOL]
                           + (*volGrid.get(x,y,z))[fsgrids::PERBYVOL];
                        retval[3*(gridSize[1]*gridSize[0]*z + gridSize[0]*y + x) + 2] = (*BgBGrid.get(x,y,z))[fsgrids::BGBZVOL]
                           + (*volGrid.get(x,y,z))[fsgrids::PERBZVOL];
                     }
                  }
               }
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"T","$\\mathrm{T}$","$B_\\mathrm{vol,fg}$","1.0");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "backgroundvolb" || lowercase == "vg_b_background_vol") {
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("vg_b_background_vol",CellParams::BGBXVOL,3));
         outputReducer->addMetadata(outputReducer->size()-1,"T","$\\mathrm{T}$","$B_\\mathrm{vol,vg,bg}$","1.0");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "perturbedvolb" || lowercase == "vg_b_perturbed_vol") {
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("vg_b_perturbed_vol",CellParams::PERBXVOL,3));
         outputReducer->addMetadata(outputReducer->size()-1,"T","$\\mathrm{T}$","$B_\\mathrm{vol,vg,per}$","1.0");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "pressure" || lowercase == "vg_pressure") {
         // Overall scalar pressure from all populations
         outputReducer->addOperator(new DRO::VariablePressureSolver);
         outputReducer->addMetadata(outputReducer->size()-1,"Pa","$\\mathrm{Pa}$","$P_\\mathrm{solver}$","1.0");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "fg_pressure") {
         // Overall scalar pressure from all populations
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_pressure",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               // Iterate through fsgrid cells and extract boundary flag
               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        auto& moments=(*momentsGrid.get(x,y,z));
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = 1./3. * (moments[fsgrids::P_11] + moments[fsgrids::P_22] + moments[fsgrids::P_33]);
                     }
                  }
               }
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"Pa","$\\mathrm{Pa}$","$P_\\mathrm{fg}$","1.0");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "populations_ptensor" || lowercase == "populations_vg_ptensor") {
         // Per-population pressure tensor, stored as diagonal and offdiagonal components
         for(unsigned int i =0; i < getObjectWrapper().particleSpecies.size(); i++) {
            species::Species& species=getObjectWrapper().particleSpecies[i];
            const std::string& pop = species.name;
            outputReducer->addOperator(new DRO::VariablePTensorDiagonal(i));
            outputReducer->addMetadata(outputReducer->size()-1,"Pa","$\\mathrm{Pa}$","$\\mathcal{P}_\\mathrm{"+pop+"}$","1.0");
            outputReducer->addOperator(new DRO::VariablePTensorOffDiagonal(i));
            outputReducer->addMetadata(outputReducer->size()-1,"Pa","$\\mathrm{Pa}$","$\\mathcal{\\tilde{P}}_\\mathrm{"+pop+"}$","1.0");
         }
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "bvolderivs" || lowercase == "b_vol_derivs" || lowercase == "b_vol_derivatives" || lowercase == "vg_b_vol_derivatives" || lowercase == "derivs") {
         // Volume-averaged derivatives
         outputReducer->addOperator(new DRO::DataReductionOperatorBVOLDerivatives("vg_derivatives/vg_dperbxvoldx",bvolderivatives::dPERBXVOLdx,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorBVOLDerivatives("vg_derivatives/vg_dperbxvoldy",bvolderivatives::dPERBXVOLdy,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorBVOLDerivatives("vg_derivatives/vg_dperbxvoldz",bvolderivatives::dPERBXVOLdz,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorBVOLDerivatives("vg_derivatives/vg_dperbyvoldx",bvolderivatives::dPERBYVOLdx,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorBVOLDerivatives("vg_derivatives/vg_dperbyvoldy",bvolderivatives::dPERBYVOLdy,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorBVOLDerivatives("vg_derivatives/vg_dperbyvoldz",bvolderivatives::dPERBYVOLdz,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorBVOLDerivatives("vg_derivatives/vg_dperbzvoldx",bvolderivatives::dPERBZVOLdx,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorBVOLDerivatives("vg_derivatives/vg_dperbzvoldy",bvolderivatives::dPERBZVOLdy,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorBVOLDerivatives("vg_derivatives/vg_dperbzvoldz",bvolderivatives::dPERBZVOLdz,1));
         outputReducer->addMetadata(outputReducer->size()-9,"T/m","$\\mathrm{T}\\,\\mathrm{m}^{-1}$","$\\Delta B_{X,\\mathrm{per,vol,vg}} (\\Delta X)^{-1}$","1.0");
         outputReducer->addMetadata(outputReducer->size()-8,"T/m","$\\mathrm{T}\\,\\mathrm{m}^{-1}$","$\\Delta B_{X,\\mathrm{per,vol,vg}} (\\Delta Y)^{-1}$","1.0");
         outputReducer->addMetadata(outputReducer->size()-7,"T/m","$\\mathrm{T}\\,\\mathrm{m}^{-1}$","$\\Delta B_{X,\\mathrm{per,vol,vg}} (\\Delta Z)^{-1}$","1.0");
         outputReducer->addMetadata(outputReducer->size()-6,"T/m","$\\mathrm{T}\\,\\mathrm{m}^{-1}$","$\\Delta B_{Y,\\mathrm{per,vol,vg}} (\\Delta X)^{-1}$","1.0");
         outputReducer->addMetadata(outputReducer->size()-5,"T/m","$\\mathrm{T}\\,\\mathrm{m}^{-1}$","$\\Delta B_{Y,\\mathrm{per,vol,vg}} (\\Delta Y)^{-1}$","1.0");
         outputReducer->addMetadata(outputReducer->size()-4,"T/m","$\\mathrm{T}\\,\\mathrm{m}^{-1}$","$\\Delta B_{Y,\\mathrm{per,vol,vg}} (\\Delta Z)^{-1}$","1.0");
         outputReducer->addMetadata(outputReducer->size()-3,"T/m","$\\mathrm{T}\\,\\mathrm{m}^{-1}$","$\\Delta B_{Z,\\mathrm{per,vol,vg}} (\\Delta X)^{-1}$","1.0");
         outputReducer->addMetadata(outputReducer->size()-2,"T/m","$\\mathrm{T}\\,\\mathrm{m}^{-1}$","$\\Delta B_{Z,\\mathrm{per,vol,vg}} (\\Delta Y)^{-1}$","1.0");
         outputReducer->addMetadata(outputReducer->size()-1,"T/m","$\\mathrm{T}\\,\\mathrm{m}^{-1}$","$\\Delta B_{Z,\\mathrm{per,vol,vg}} (\\Delta Z)^{-1}$","1.0");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }


      // This long block is writing all the derivatives we are storing on the fsgrids
      // !! EXCEPT background b !!
      // that is, derivatives of perturbed b, perturbed bvol, rhom, rhoq, v, p11, p22, p33, and pe.
      // Note that so far we are not computing nor storing the fg_dperbidi components, hence we cannot write them out!
      // We do have the background ones, as well as the fg_dperbivoldi and their vg equivalent (elsewhere) as those are computed from the perbivol components.
      // As of summer 2023 they are proper derivatives in DROs, unlike in the code where they are differences.
      // Search for "fg_derivs" to find the end of this block.
      if(P::systemWriteAllDROs || lowercase == "fg_derivs") {
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_derivatives/fg_dperbxdy",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = dPerBGrid.get(x,y,z)->at(fsgrids::dperb::dPERBxdy) / dPerBGrid.DY;
                     }
                  }
               }
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"T/m","$\\mathrm{T}\\,\\mathrm{m}^{-1}$","$\\Delta B_{X,\\mathrm{per,fg}} (\\Delta Y)^{-1}$","1.0");
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_derivatives/fg_dperbxdz",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = dPerBGrid.get(x,y,z)->at(fsgrids::dperb::dPERBxdz) / dPerBGrid.DZ;
                     }
                  }
               }
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"T/m","$\\mathrm{T}\\,\\mathrm{m}^{-1}$","$\\Delta B_{X,\\mathrm{per,fg}} (\\Delta Z)^{-1}$","1.0");
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_derivatives/fg_dperbydx",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = dPerBGrid.get(x,y,z)->at(fsgrids::dperb::dPERBydx) / dPerBGrid.DX;
                     }
                  }
               }
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"T/m","$\\mathrm{T}\\,\\mathrm{m}^{-1}$","$\\Delta B_{Y,\\mathrm{per,fg}} (\\Delta X)^{-1}$","1.0");
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_derivatives/fg_dperbydz",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = dPerBGrid.get(x,y,z)->at(fsgrids::dperb::dPERBydz) / dPerBGrid.DZ;
                     }
                  }
               }
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"T/m","$\\mathrm{T}\\,\\mathrm{m}^{-1}$","$\\Delta B_{Y,\\mathrm{per,fg}} (\\Delta Z)^{-1}$","1.0");
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_derivatives/fg_dperbzdx",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = dPerBGrid.get(x,y,z)->at(fsgrids::dperb::dPERBzdx) / dPerBGrid.DX;
                     }
                  }
               }
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"T/m","$\\mathrm{T}\\,\\mathrm{m}^{-1}$","$\\Delta B_{Z,\\mathrm{per,fg}} (\\Delta X)^{-1}$","1.0");
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_derivatives/fg_dperbzdy",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = dPerBGrid.get(x,y,z)->at(fsgrids::dperb::dPERBzdy) / dPerBGrid.DY;
                     }
                  }
               }
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"T/m","$\\mathrm{T}\\,\\mathrm{m}^{-1}$","$\\Delta B_{Z,\\mathrm{per,fg}} (\\Delta Y)^{-1}$","1.0");
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_derivatives/fg_dperbxdyy",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = dPerBGrid.get(x,y,z)->at(fsgrids::dperb::dPERBxdyy) / dPerBGrid.DY / dPerBGrid.DY;
                     }
                  }
               }
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"T/m","$\\mathrm{T}\\,\\mathrm{m}^{-2}$","$\\Delta B_{X,\\mathrm{per,fg}} (\\Delta Y)^{-2}$","1.0");
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_derivatives/fg_dperbxdzz",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = dPerBGrid.get(x,y,z)->at(fsgrids::dperb::dPERBxdzz) / dPerBGrid.DZ / dPerBGrid.DZ;
                     }
                  }
               }
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"T/m","$\\mathrm{T}\\,\\mathrm{m}^{-2}$","$\\Delta B_{X,\\mathrm{per,fg}} (\\Delta Z)^{-2}$","1.0");
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_derivatives/fg_dperbxdyz",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = dPerBGrid.get(x,y,z)->at(fsgrids::dperb::dPERBxdyz) / dPerBGrid.DY / dPerBGrid.DZ;
                     }
                  }
               }
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"T/m","$\\mathrm{T}\\,\\mathrm{m}^{-2}$","$\\Delta B_{X,\\mathrm{per,fg}} (\\Delta Y \\Delta Z)^{-1}$","1.0");
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_derivatives/fg_dperbydxx",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = dPerBGrid.get(x,y,z)->at(fsgrids::dperb::dPERBydxx) / dPerBGrid.DX / dPerBGrid.DX;
                     }
                  }
               }
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"T/m","$\\mathrm{T}\\,\\mathrm{m}^{-2}$","$\\Delta B_{Y,\\mathrm{per,fg}} (\\Delta X)^{-2}$","1.0");
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_derivatives/fg_dperbydzz",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = dPerBGrid.get(x,y,z)->at(fsgrids::dperb::dPERBydzz) / dPerBGrid.DZ / dPerBGrid.DZ;
                     }
                  }
               }
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"T/m","$\\mathrm{T}\\,\\mathrm{m}^{-2}$","$\\Delta B_{Y,\\mathrm{per,fg}} (\\Delta Z)^{-2}$","1.0");
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_derivatives/fg_dperbydxz",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = dPerBGrid.get(x,y,z)->at(fsgrids::dperb::dPERBydxz) / dPerBGrid.DX / dPerBGrid.DZ;
                     }
                  }
               }
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"T/m","$\\mathrm{T}\\,\\mathrm{m}^{-2}$","$\\Delta B_{Y,\\mathrm{per,fg}} (\\Delta X \\Delta Z)^{-1}$","1.0");
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_derivatives/fg_dperbzdxx",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = dPerBGrid.get(x,y,z)->at(fsgrids::dperb::dPERBzdxx) / dPerBGrid.DX / dPerBGrid.DX;
                     }
                  }
               }
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"T/m","$\\mathrm{T}\\,\\mathrm{m}^{-2}$","$\\Delta B_{Z,\\mathrm{per,fg}} (\\Delta Z)^{-2}$","1.0");
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_derivatives/fg_dperbzdyy",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = dPerBGrid.get(x,y,z)->at(fsgrids::dperb::dPERBzdyy) / dPerBGrid.DY / dPerBGrid.DY;
                     }
                  }
               }
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"T/m","$\\mathrm{T}\\,\\mathrm{m}^{-2}$","$\\Delta B_{Z,\\mathrm{per,fg}} (\\Delta Y)^{-2}$","1.0");
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_derivatives/fg_dperbzdxy",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = dPerBGrid.get(x,y,z)->at(fsgrids::dperb::dPERBzdxy) / dPerBGrid.DX / dPerBGrid.DY;
                     }
                  }
               }
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"T/m","$\\mathrm{T}\\,\\mathrm{m}^{-2}$","$\\Delta B_{Z,\\mathrm{per,fg}} (\\Delta X \\Delta Y)^{-1}$","1.0");


         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_derivatives/fg_drhomdx",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = dMomentsGrid.get(x,y,z)->at(fsgrids::dmoments::drhomdx) / dMomentsGrid.DX;
                     }
                  }
               }
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"kg/m^4","$\\mathrm{kg}\\mathrm{m}^{-4}$","$\\Delta \\rho_{m,\\mathrm{fg}} (\\Delta X)^{-1}$","1.0");
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_derivatives/fg_drhomdy",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = dMomentsGrid.get(x,y,z)->at(fsgrids::dmoments::drhomdy) / dMomentsGrid.DY;
                     }
                  }
               }
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"kg/m^4","$\\mathrm{kg}\\mathrm{m}^{-4}$","$\\Delta \\rho_{m,\\mathrm{fg}} (\\Delta Y)^{-1}$","1.0");
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_derivatives/fg_drhomdz",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = dMomentsGrid.get(x,y,z)->at(fsgrids::dmoments::drhomdz) / dMomentsGrid.DZ;
                     }
                  }
               }
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"kg/m^4","$\\mathrm{kg}\\mathrm{m}^{-4}$","$\\Delta \\rho_{m,\\mathrm{fg}} (\\Delta Z)^{-1}$","1.0");
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_derivatives/fg_drhoqdx",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = dMomentsGrid.get(x,y,z)->at(fsgrids::dmoments::drhoqdx) / dMomentsGrid.DX;
                     }
                  }
               }
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"C/m^4","$\\mathrm{C}\\mathrm{m}^{-4}$","$\\Delta \\rho_{q,\\mathrm{fg}} (\\Delta X)^{-1}$","1.0");
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_derivatives/fg_drhoqdy",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = dMomentsGrid.get(x,y,z)->at(fsgrids::dmoments::drhoqdy) / dMomentsGrid.DY;
                     }
                  }
               }
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"C/m^4","$\\mathrm{C}\\mathrm{m}^{-4}$","$\\Delta \\rho_{q,\\mathrm{fg}} (\\Delta Y)^{-1}$","1.0");
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_derivatives/fg_drhoqdz",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = dMomentsGrid.get(x,y,z)->at(fsgrids::dmoments::drhoqdz) / dMomentsGrid.DZ;
                     }
                  }
               }
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"C/m^4","$\\mathrm{C}\\mathrm{m}^{-4}$","$\\Delta \\rho_{q,\\mathrm{fg}} (\\Delta Z)^{-1}$","1.0");
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_derivatives/fg_dp11dx",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = dMomentsGrid.get(x,y,z)->at(fsgrids::dmoments::dp11dx) / dMomentsGrid.DX;
                     }
                  }
               }
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"Pa/m","$\\mathrm{Pa}\\mathrm{m}^{-1}$","$\\Delta P_{11,\\mathrm{fg}} (\\Delta X)^{-1}$","1.0");
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_derivatives/fg_dp11dy",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = dMomentsGrid.get(x,y,z)->at(fsgrids::dmoments::dp11dy) / dMomentsGrid.DY;
                     }
                  }
               }
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"Pa/m","$\\mathrm{Pa}\\mathrm{m}^{-1}$","$\\Delta P_{11,\\mathrm{fg}} (\\Delta Y)^{-1}$","1.0");
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_derivatives/fg_dp11dz",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = dMomentsGrid.get(x,y,z)->at(fsgrids::dmoments::dp11dz) / dMomentsGrid.DZ;
                     }
                  }
               }
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"Pa/m","$\\mathrm{Pa}\\mathrm{m}^{-1}$","$\\Delta P_{11,\\mathrm{fg}} (\\Delta Z)^{-1}$","1.0");
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_derivatives/fg_dp22dx",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = dMomentsGrid.get(x,y,z)->at(fsgrids::dmoments::dp22dx) / dMomentsGrid.DX;
                     }
                  }
               }
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"Pa/m","$\\mathrm{Pa}\\mathrm{m}^{-1}$","$\\Delta P_{22,\\mathrm{fg}} (\\Delta X)^{-1}$","1.0");
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_derivatives/fg_dp22dy",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = dMomentsGrid.get(x,y,z)->at(fsgrids::dmoments::dp22dy) / dMomentsGrid.DY;
                     }
                  }
               }
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"Pa/m","$\\mathrm{Pa}\\mathrm{m}^{-1}$","$\\Delta P_{22,\\mathrm{fg}} (\\Delta Y)^{-1}$","1.0");
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_derivatives/fg_dp22dz",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = dMomentsGrid.get(x,y,z)->at(fsgrids::dmoments::dp22dz) / dMomentsGrid.DZ;
                     }
                  }
               }
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"Pa/m","$\\mathrm{Pa}\\mathrm{m}^{-1}$","$\\Delta P_{22,\\mathrm{fg}} (\\Delta Z)^{-1}$","1.0");
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_derivatives/fg_dp33dx",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = dMomentsGrid.get(x,y,z)->at(fsgrids::dmoments::dp33dx) / dMomentsGrid.DX;
                     }
                  }
               }
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"Pa/m","$\\mathrm{Pa}\\mathrm{m}^{-1}$","$\\Delta P_{33,\\mathrm{fg}} (\\Delta X)^{-1}$","1.0");
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_derivatives/fg_dp33dy",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = dMomentsGrid.get(x,y,z)->at(fsgrids::dmoments::dp33dy) / dMomentsGrid.DY;
                     }
                  }
               }
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"Pa/m","$\\mathrm{Pa}\\mathrm{m}^{-1}$","$\\Delta P_{33,\\mathrm{fg}} (\\Delta Y)^{-1}$","1.0");
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_derivatives/fg_dp33dz",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = dMomentsGrid.get(x,y,z)->at(fsgrids::dmoments::dp33dz) / dMomentsGrid.DZ;
                     }
                  }
               }
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"Pa/m","$\\mathrm{Pa}\\mathrm{m}^{-1}$","$\\Delta P_{33,\\mathrm{fg}} (\\Delta Z)^{-1}$","1.0");
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_derivatives/fg_dvxdx",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = dMomentsGrid.get(x,y,z)->at(fsgrids::dmoments::dVxdx) / dMomentsGrid.DX;
                     }
                  }
               }
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"1/s","$\\mathrm{s}^{-1}$","$\\Delta V_{X,\\mathrm{fg}} (\\Delta X)^{-1}$","1.0");
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_derivatives/fg_dvxdy",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = dMomentsGrid.get(x,y,z)->at(fsgrids::dmoments::dVxdy) / dMomentsGrid.DY;
                     }
                  }
               }
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"1/s","$\\mathrm{s}^{-1}$","$\\Delta V_{X,\\mathrm{fg}} (\\Delta Y)^{-1}$","1.0");
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_derivatives/fg_dvxdz",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = dMomentsGrid.get(x,y,z)->at(fsgrids::dmoments::dVxdz) / dMomentsGrid.DZ;
                     }
                  }
               }
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"1/s","$\\mathrm{s}^{-1}$","$\\Delta V_{X,\\mathrm{fg}} (\\Delta Z)^{-1}$","1.0");
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_derivatives/fg_dvydx",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = dMomentsGrid.get(x,y,z)->at(fsgrids::dmoments::dVydx) / dMomentsGrid.DX;
                     }
                  }
               }
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"1/s","$\\mathrm{s}^{-1}$","$\\Delta V_{Y,\\mathrm{fg}} (\\Delta X)^{-1}$","1.0");
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_derivatives/fg_dvydy",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = dMomentsGrid.get(x,y,z)->at(fsgrids::dmoments::dVydy) / dMomentsGrid.DY;
                     }
                  }
               }
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"1/s","$\\mathrm{s}^{-1}$","$\\Delta V_{Y,\\mathrm{fg}} (\\Delta Y)^{-1}$","1.0");
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_derivatives/fg_dvydz",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = dMomentsGrid.get(x,y,z)->at(fsgrids::dmoments::dVydz) / dMomentsGrid.DZ;
                     }
                  }
               }
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"1/s","$\\mathrm{s}^{-1}$","$\\Delta V_{Y,\\mathrm{fg}} (\\Delta Z)^{-1}$","1.0");
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_derivatives/fg_dvzdx",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = dMomentsGrid.get(x,y,z)->at(fsgrids::dmoments::dVzdx) / dMomentsGrid.DX;
                     }
                  }
               }
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"1/s","$\\mathrm{s}^{-1}$","$\\Delta V_{Z,\\mathrm{fg}} (\\Delta X)^{-1}$","1.0");
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_derivatives/fg_dvzdy",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = dMomentsGrid.get(x,y,z)->at(fsgrids::dmoments::dVzdy) / dMomentsGrid.DY;
                     }
                  }
               }
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"1/s","$\\mathrm{s}^{-1}$","$\\Delta V_{Z,\\mathrm{fg}} (\\Delta Y)^{-1}$","1.0");
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_derivatives/fg_dvzdz",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = dMomentsGrid.get(x,y,z)->at(fsgrids::dmoments::dVzdz) / dMomentsGrid.DZ;
                     }
                  }
               }
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"1/s","$\\mathrm{s}^{-1}$","$\\Delta V_{Z,\\mathrm{fg}} (\\Delta Z)^{-1}$","1.0");
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_derivatives/fg_dpedx",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = dMomentsGrid.get(x,y,z)->at(fsgrids::dmoments::dPedx) / dMomentsGrid.DX;
                     }
                  }
               }
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"Pa/m","$\\mathrm{Pa}\\mathrm{m}^{-1}$","$\\Delta P_\\mathrm{e,fg} (\\Delta X)^{-1}$","1.0");
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_derivatives/fg_dpedy",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = dMomentsGrid.get(x,y,z)->at(fsgrids::dmoments::dPedy) / dMomentsGrid.DY;
                     }
                  }
               }
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"Pa/m","$\\mathrm{Pa}\\mathrm{m}^{-1}$","$\\Delta P_\\mathrm{e,fg} (\\Delta Y)^{-1}$","1.0");
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_derivatives/fg_dpedz",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = dMomentsGrid.get(x,y,z)->at(fsgrids::dmoments::dPedz) / dMomentsGrid.DZ;
                     }
                  }
               }
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"Pa/m","$\\mathrm{Pa}\\mathrm{m}^{-1}$","$\\Delta P_\\mathrm{e,fg} (\\Delta Z)^{-1}$","1.0");

         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_derivatives/fg_dperbxvoldx",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = volGrid.get(x,y,z)->at(fsgrids::volfields::dPERBXVOLdx) / BgBGrid.DX;
                     }
                  }
               }
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"T/m","$\\mathrm{T}\\mathrm{m}^{-1}$","$\\Delta B_{X,\\mathrm{per,vol,fg}} (\\Delta X)^{-1}$","1.0");
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_derivatives/fg_dperbxvoldy",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = volGrid.get(x,y,z)->at(fsgrids::volfields::dPERBXVOLdy) / BgBGrid.DY;
                     }
                  }
               }
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"T/m","$\\mathrm{T}\\mathrm{m}^{-1}$","$\\Delta B_{X,\\mathrm{per,vol,fg}} (\\Delta Y)^{-1}$","1.0");
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_derivatives/fg_dperbxvoldz",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = volGrid.get(x,y,z)->at(fsgrids::volfields::dPERBXVOLdz) / BgBGrid.DZ;
                     }
                  }
               }
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"T/m","$\\mathrm{T}\\mathrm{m}^{-1}$","$\\Delta B_{X,\\mathrm{per,vol,fg}} (\\Delta Z)^{-1}$","1.0");
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_derivatives/fg_dperbyvoldx",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = volGrid.get(x,y,z)->at(fsgrids::volfields::dPERBYVOLdx) / BgBGrid.DX;
                     }
                  }
               }
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"T/m","$\\mathrm{T}\\mathrm{m}^{-1}$","$\\Delta B_{Y,\\mathrm{per,vol,fg}} (\\Delta X)^{-1}$","1.0");
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_derivatives/fg_dperbyvoldy",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = volGrid.get(x,y,z)->at(fsgrids::volfields::dPERBYVOLdy) / BgBGrid.DY;
                     }
                  }
               }
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"T/m","$\\mathrm{T}\\mathrm{m}^{-1}$","$\\Delta B_{Y,\\mathrm{per,vol,fg}} (\\Delta Y)^{-1}$","1.0");

         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_derivatives/fg_dperbyvoldz",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = volGrid.get(x,y,z)->at(fsgrids::volfields::dPERBYVOLdz) / BgBGrid.DZ;
                     }
                  }
               }
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"T/m","$\\mathrm{T}\\mathrm{m}^{-1}$","$\\Delta B_{Y,\\mathrm{per,vol,fg}} (\\Delta Z)^{-1}$","1.0");
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_derivatives/fg_dperbzvoldx",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = volGrid.get(x,y,z)->at(fsgrids::volfields::dPERBZVOLdx) / BgBGrid.DX;
                     }
                  }
               }
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"T/m","$\\mathrm{T}\\mathrm{m}^{-1}$","$\\Delta B_{Z,\\mathrm{per,vol,fg}} (\\Delta X)^{-1}$","1.0");
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_derivatives/fg_dperbzvoldy",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = volGrid.get(x,y,z)->at(fsgrids::volfields::dPERBZVOLdy) / BgBGrid.DY;
                     }
                  }
               }
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"T/m","$\\mathrm{T}\\mathrm{m}^{-1}$","$\\Delta B_{Z,\\mathrm{per,vol,fg}} (\\Delta Y)^{-1}$","1.0");
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_derivatives/fg_dperbzvoldz",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = volGrid.get(x,y,z)->at(fsgrids::volfields::dPERBZVOLdz) / BgBGrid.DZ;
                     }
                  }
               }
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"T/m","$\\mathrm{T}\\mathrm{m}^{-1}$","$\\Delta B_{Z,\\mathrm{per,vol,fg}} (\\Delta Z)^{-1}$","1.0");

         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      // End of the long block for fg_derivs
      // that writes all the fsgrid-stored derivatives.

      // The following long block writes all the background magnetic field derivatives
      // we store on fsgrid, that is fg_dbgbidj and fg_dbgbivoldj (also i==j).
      // They are derivatives in these DROs, not differences as in the code.
      // Search for "fg_derivs_b_background" to find the end of the block.
      if(P::systemWriteAllDROs || lowercase == "fg_derivs_b_background") { // includes all face and volume-averaged derivatives of BGB on fg
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_derivatives/fg_dbgbxdy",[](
            FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
            FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
            FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
            FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
            FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
            FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
            FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
            FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
            FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
            FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = BgBGrid.get(x,y,z)->at(fsgrids::bgbfield::dBGBxdy) / technicalGrid.DY;
                     }
                  }
               }
               return retval;
            }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"T/m","$\\mathrm{T}\\,\\mathrm{m}^{-1}$","$\\Delta B_{X,\\mathrm{bg,fg}} (\\Delta Y)^{-1}$","1.0");

         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_derivatives/fg_dbgbxdz",[](
            FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
            FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
            FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
            FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
            FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
            FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
            FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
            FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
            FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
            FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = BgBGrid.get(x,y,z)->at(fsgrids::bgbfield::dBGBxdz) / technicalGrid.DZ;
                     }
                  }
               }
               return retval;
            }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"T/m","$\\mathrm{T}\\,\\mathrm{m}^{-1}$","$\\Delta B_{X,\\mathrm{bg,fg}} (\\Delta Z)^{-1}$","1.0");

         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_derivatives/fg_dbgbydx",[](
            FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
            FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
            FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
            FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
            FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
            FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
            FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
            FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
            FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
            FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = BgBGrid.get(x,y,z)->at(fsgrids::bgbfield::dBGBydx) / technicalGrid.DX;
                     }
                  }
               }
               return retval;
            }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"T/m","$\\mathrm{T}\\,\\mathrm{m}^{-1}$","$\\Delta B_{Y,\\mathrm{bg,fg}} (\\Delta X)^{-1}$","1.0");

         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_derivatives/fg_dbgbydz",[](
            FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
            FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
            FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
            FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
            FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
            FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
            FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
            FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
            FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
            FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = BgBGrid.get(x,y,z)->at(fsgrids::bgbfield::dBGBydz) / technicalGrid.DZ;
                     }
                  }
               }
               return retval;
            }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"T/m","$\\mathrm{T}\\,\\mathrm{m}^{-1}$","$\\Delta B_{Y,\\mathrm{bg,fg}} (\\Delta Z)^{-1}$","1.0");

         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_derivatives/fg_dbgbzdx",[](
            FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
            FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
            FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
            FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
            FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
            FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
            FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
            FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
            FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
            FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = BgBGrid.get(x,y,z)->at(fsgrids::bgbfield::dBGBzdx) / technicalGrid.DX;
                     }
                  }
               }
               return retval;
            }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"T/m","$\\mathrm{T}\\,\\mathrm{m}^{-1}$","$\\Delta B_{Z,\\mathrm{bg,fg}} (\\Delta X)^{-1}$","1.0");

         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_derivatives/fg_dbgbzdy",[](
            FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
            FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
            FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
            FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
            FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
            FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
            FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
            FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
            FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
            FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = BgBGrid.get(x,y,z)->at(fsgrids::bgbfield::dBGBzdy) / technicalGrid.DY;
                     }
                  }
               }
               return retval;
            }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"T/m","$\\mathrm{T}\\,\\mathrm{m}^{-1}$","$\\Delta B_{Z,\\mathrm{bg,fg}} (\\Delta Y)^{-1}$","1.0");

         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_derivatives/fg_dbgbxvoldx",[](
            FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
            FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
            FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
            FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
            FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
            FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
            FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
            FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
            FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
            FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = BgBGrid.get(x,y,z)->at(fsgrids::bgbfield::dBGBXVOLdx) / technicalGrid.DX;
                     }
                  }
               }
               return retval;
            }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"T/m","$\\mathrm{T}\\,\\mathrm{m}^{-1}$","$\\Delta B_{X,\\mathrm{bg,vol,fg}} (\\Delta X)^{-1}$","1.0");

         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_derivatives/fg_dbgbxvoldy",[](
            FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
            FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
            FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
            FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
            FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
            FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
            FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
            FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
            FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
            FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = BgBGrid.get(x,y,z)->at(fsgrids::bgbfield::dBGBXVOLdy) / technicalGrid.DY;
                     }
                  }
               }
               return retval;
            }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"T/m","$\\mathrm{T}\\,\\mathrm{m}^{-1}$","$\\Delta B_{X,\\mathrm{bg,vol,fg}} (\\Delta Y)^{-1}$","1.0");

         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_derivatives/fg_dbgbxvoldz",[](
            FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
            FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
            FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
            FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
            FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
            FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
            FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
            FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
            FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
            FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = BgBGrid.get(x,y,z)->at(fsgrids::bgbfield::dBGBXVOLdz) / technicalGrid.DZ;
                     }
                  }
               }
               return retval;
            }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"T/m","$\\mathrm{T}\\,\\mathrm{m}^{-1}$","$\\Delta B_{X,\\mathrm{bg,vol,fg}} (\\Delta Z)^{-1}$","1.0");

         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_derivatives/fg_dbgbyvoldx",[](
            FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
            FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
            FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
            FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
            FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
            FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
            FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
            FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
            FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
            FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = BgBGrid.get(x,y,z)->at(fsgrids::bgbfield::dBGBYVOLdx) / technicalGrid.DX;
                     }
                  }
               }
               return retval;
            }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"T/m","$\\mathrm{T}\\,\\mathrm{m}^{-1}$","$\\Delta B_{Y,\\mathrm{bg,vol,fg}} (\\Delta X)^{-1}$","1.0");

         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_derivatives/fg_dbgbyvoldy",[](
            FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
            FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
            FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
            FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
            FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
            FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
            FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
            FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
            FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
            FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = BgBGrid.get(x,y,z)->at(fsgrids::bgbfield::dBGBYVOLdy) / technicalGrid.DY;
                     }
                  }
               }
               return retval;
            }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"T/m","$\\mathrm{T}\\,\\mathrm{m}^{-1}$","$\\Delta B_{Y,\\mathrm{bg,vol,fg}} (\\Delta Y)^{-1}$","1.0");

         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_derivatives/fg_dbgbyvoldz",[](
            FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
            FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
            FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
            FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
            FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
            FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
            FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
            FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
            FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
            FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = BgBGrid.get(x,y,z)->at(fsgrids::bgbfield::dBGBYVOLdz) / technicalGrid.DZ;
                     }
                  }
               }
               return retval;
            }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"T/m","$\\mathrm{T}\\,\\mathrm{m}^{-1}$","$\\Delta B_{Y,\\mathrm{bg,vol,fg}} (\\Delta Z)^{-1}$","1.0");

         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_derivatives/fg_dbgbzvoldx",[](
            FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
            FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
            FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
            FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
            FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
            FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
            FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
            FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
            FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
            FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = BgBGrid.get(x,y,z)->at(fsgrids::bgbfield::dBGBZVOLdx) / technicalGrid.DX;
                     }
                  }
               }
               return retval;
            }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"T/m","$\\mathrm{T}\\,\\mathrm{m}^{-1}$","$\\Delta B_{Z,\\mathrm{bg,vol,fg}} (\\Delta X)^{-1}$","1.0");

         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_derivatives/fg_dbgbzvoldy",[](
            FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
            FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
            FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
            FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
            FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
            FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
            FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
            FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
            FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
            FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = BgBGrid.get(x,y,z)->at(fsgrids::bgbfield::dBGBZVOLdy) / technicalGrid.DY;
                     }
                  }
               }
               return retval;
            }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"T/m","$\\mathrm{T}\\,\\mathrm{m}^{-1}$","$\\Delta B_{Z,\\mathrm{bg,vol,fg}} (\\Delta Y)^{-1}$","1.0");

         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_derivatives/fg_dbgbzvoldz",[](
            FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
            FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
            FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
            FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
            FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
            FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
            FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
            FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
            FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
            FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = BgBGrid.get(x,y,z)->at(fsgrids::bgbfield::dBGBZVOLdz) / technicalGrid.DZ;
                     }
                  }
               }
               return retval;
            }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"T/m","$\\mathrm{T}\\,\\mathrm{m}^{-1}$","$\\Delta B_{Z,\\mathrm{bg,vol,fg}} (\\Delta Z)^{-1}$","1.0");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      // fg_derivs_b_background
      // End fo the long block writing out all the background magnetic field derivatives from fsgrid.

      if(P::systemWriteAllDROs || lowercase == "vg_gridcoordinates") {
         // Spatial coordinates for each cell
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("vg_x",CellParams::XCRD,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("vg_y",CellParams::YCRD,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("vg_z",CellParams::ZCRD,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("vg_dx",CellParams::DX,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("vg_dy",CellParams::DY,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("vg_dz",CellParams::DZ,1));
         outputReducer->addMetadata(outputReducer->size()-6,"m","$\\mathrm{m}$","$X_\\mathrm{vg}$","1.0");
         outputReducer->addMetadata(outputReducer->size()-5,"m","$\\mathrm{m}$","$Y_\\mathrm{vg}$","1.0");
         outputReducer->addMetadata(outputReducer->size()-4,"m","$\\mathrm{m}$","$Z_\\mathrm{vg}$","1.0");
         outputReducer->addMetadata(outputReducer->size()-3,"m","$\\mathrm{m}$","$\\delta X_\\mathrm{vg}$","1.0");
         outputReducer->addMetadata(outputReducer->size()-2,"m","$\\mathrm{m}$","$\\delta Y_\\mathrm{vg}$","1.0");
         outputReducer->addMetadata(outputReducer->size()-1,"m","$\\mathrm{m}$","$\\delta Z_\\mathrm{vg}$","1.0");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "fg_gridcoordinates") {
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_x",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               // Iterate through fsgrid cells and extract X coordinate
               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = technicalGrid.getPhysicalCoords(x,y,z)[0];
                     }
                  }
               }
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"m","$\\mathrm{m}$","$X_\\mathrm{fg}$","1.0");
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_y",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               // Iterate through fsgrid cells and extract Y coordinate
               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = technicalGrid.getPhysicalCoords(x,y,z)[1];
                     }
                  }
               }
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"m","$\\mathrm{m}$","$Y_\\mathrm{fg}$","1.0");
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_z",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]);

               // Iterate through fsgrid cells and extract Z coordinate
               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[gridSize[1]*gridSize[0]*z + gridSize[0]*y + x] = technicalGrid.getPhysicalCoords(x,y,z)[2];
                     }
                  }
               }
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"m","$\\mathrm{m}$","$Z_\\mathrm{fg}$","1.0");
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_dx",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2], technicalGrid.DX);
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"m","$\\mathrm{m}$","$\\delta X_\\mathrm{fg}$","1.0");
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_dy",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2], technicalGrid.DY);
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"m","$\\mathrm{m}$","$\\delta Y_\\mathrm{fg}$","1.0");
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_dz",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2], technicalGrid.DZ);
               return retval;
         }
         ));
         outputReducer->addMetadata(outputReducer->size()-1,"m","$\\mathrm{m}$","$\\delta Z_\\mathrm{fg}$","1.0");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "vg_amr_drho") {
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("vg_amr_drho",CellParams::AMR_DRHO,1));
         outputReducer->addMetadata(outputReducer->size()-1,"","","$\\frac{\\Delta \\rho}{\\hat{rho}}$","");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "vg_amr_du") {
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("vg_amr_du",CellParams::AMR_DU,1));
         outputReducer->addMetadata(outputReducer->size()-1,"","","$\\frac{\\Delta U_1}{\\hat{U}_1}$","");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "vg_amr_dpsq") {
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("vg_amr_dpsq",CellParams::AMR_DPSQ,1));
         outputReducer->addMetadata(outputReducer->size()-1,"","","$\\frac{(\\Delta P)^2}{2 \\rho \\hat{U}_1}$","");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "vg_amr_dbsq") {
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("vg_amr_dbsq",CellParams::AMR_DBSQ,1));
         outputReducer->addMetadata(outputReducer->size()-1,"","","$\\frac{(\\Delta B_1)^2}{2 \\mu_0 \\hat{U}_1}$","");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "vg_amr_db") {
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("vg_amr_db",CellParams::AMR_DB,1));
         outputReducer->addMetadata(outputReducer->size()-1,"","","$\\frac{|\\Delta B_1|}{\\hat{B}_1}$","");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "vg_amr_alpha1") {
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("vg_amr_alpha1",CellParams::AMR_ALPHA1,1));
         outputReducer->addMetadata(outputReducer->size()-1,"","","$\\alpha_1$","");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "vg_amr_reflevel") {
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("vg_amr_reflevel",CellParams::REFINEMENT_LEVEL,1));
         outputReducer->addMetadata(outputReducer->size()-1,"","","ref","");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "vg_amr_alpha2") {
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("vg_amr_alpha2",CellParams::AMR_ALPHA2,1));
         outputReducer->addMetadata(outputReducer->size()-1,"","","$\\alpha_2$","");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "vg_pressure_anisotropy") {
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("vg_pressure_anisotropy",CellParams::P_ANISOTROPY,1));
         outputReducer->addMetadata(outputReducer->size()-1,"","","$P_\\perp / P_\\parallel$","");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "vg_amr_vorticity") {
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("vg_amr_vorticity",CellParams::AMR_VORTICITY,1));
         outputReducer->addMetadata(outputReducer->size()-1,"","","Vorticity","");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "ig_latitude") {
         outputReducer->addOperator(new DRO::DataReductionOperatorIonosphereNode("ig_latitude", [](
                     SBC::SphericalTriGrid& grid)->std::vector<Real> {

                     std::vector<Real> retval(grid.nodes.size());

                     for(uint i=0; i<grid.nodes.size(); i++) {

                        // TODO: This is geographic latitude. Should it be magnetic?
                        Real z = grid.nodes[i].x[2];
                        retval[i] = acos(z/SBC::Ionosphere::innerRadius) / M_PI * 180.;
                     }

                     return retval;
                  }));
         outputReducer->addMetadata(outputReducer->size()-1, "Degrees", "$^\\circ$", "L", "");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "ig_chi0") {
         outputReducer->addOperator(new DRO::DataReductionOperatorIonosphereNode("ig_chi0", [](
                     SBC::SphericalTriGrid& grid)->std::vector<Real> {

                     std::vector<Real> retval(grid.nodes.size());

                     for(uint n=0; n<grid.nodes.size(); n++) {
                        Real theta = acos(grid.nodes[n].x[2] / sqrt(grid.nodes[n].x[0]*grid.nodes[n].x[0] + grid.nodes[n].x[1]*grid.nodes[n].x[1] + grid.nodes[n].x[2]*grid.nodes[n].x[2])); // Latitude
                        if(theta > M_PI/2.) {
                           theta = M_PI - theta;
                        }
                        // Smoothstep with an edge at about 67 deg.
                        Real Chi0 = 0.01 + 0.99 * .5 * (1 + tanh((23. - theta * (180. / M_PI)) / 6));
                        retval[n] = Chi0;
                     }

                     return retval;
                  }));
         outputReducer->addMetadata(outputReducer->size()-1, "arb.unit.", "", "Chi0", "");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "ig_cellarea") {
         outputReducer->addOperator(new DRO::DataReductionOperatorIonosphereElement("ig_cellarea", [](
                     SBC::SphericalTriGrid& grid)->std::vector<Real> {

                     std::vector<Real> retval(grid.elements.size());

                     for(uint i=0; i<grid.elements.size(); i++) {
                        retval[i] = grid.elementArea(i);
                     }

                     return retval;
                     }));
         outputReducer->addMetadata(outputReducer->size()-1, "m^2", "$\\mathrm{m}^2$", "$A_m$", "1.0");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "ig_b") {
         outputReducer->addOperator(new DRO::DataReductionOperatorIonosphereNode("ig_b", [](
                     SBC::SphericalTriGrid& grid)->std::vector<Real> {

                     std::vector<Real> retval(grid.nodes.size()*3);

                     for(uint i=0; i<grid.nodes.size(); i++) {
                        retval[3*i] = grid.nodes[i].parameters[ionosphereParameters::NODE_BX];
                        retval[3*i+1] = grid.nodes[i].parameters[ionosphereParameters::NODE_BY];
                        retval[3*i+2] = grid.nodes[i].parameters[ionosphereParameters::NODE_BZ];
                     }

                     return retval;
                     }));
         outputReducer->addMetadata(outputReducer->size()-1, "T", "$\\mathrm{T}$", "$B$", "1.0");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }

      if(P::systemWriteAllDROs || lowercase == "ig_e") {
         outputReducer->addOperator(new DRO::DataReductionOperatorIonosphereElement("ig_e", [](
                     SBC::SphericalTriGrid& grid)->std::vector<Real> {

                     std::vector<Real> retval(grid.elements.size()*3);

                     for(uint i=0; i<grid.elements.size(); i++) {
                        // Calculate E from element basis functions
                        const std::array<Real, 3>& c1 = grid.nodes[grid.elements[i].corners[0]].x;
                        const std::array<Real, 3>& c2 = grid.nodes[grid.elements[i].corners[1]].x;
                        const std::array<Real, 3>& c3 = grid.nodes[grid.elements[i].corners[2]].x;

                        // ET contains the test function gradient vectors (normalized to potential 1)
                        std::array<std::array<Real,3>, 3> ET({grid.computeGradT(c2,c3,c1), grid.computeGradT(c3,c1,c2), grid.computeGradT(c1,c2,c3)});
                        for(int n=0; n<3; n++) {
                           // Multiply with the corresponding node potentials to get E vector
                           ET[0][n] *= -grid.nodes[grid.elements[i].corners[0]].parameters[ionosphereParameters::SOLUTION];
                           ET[1][n] *= -grid.nodes[grid.elements[i].corners[1]].parameters[ionosphereParameters::SOLUTION];
                           ET[2][n] *= -grid.nodes[grid.elements[i].corners[2]].parameters[ionosphereParameters::SOLUTION];
                        }
                        // Sum up element gradient functions to yield complete E inside this element.
                        for(int n=0; n<3; n++) {
                           retval[3*i + n] = ET[0][n] + ET[1][n] + ET[2][n];
                        }
                     }

                     return retval;
                     }));
         outputReducer->addMetadata(outputReducer->size()-1, "V/m", "$\\mathrm{V/m}$", "$E$", "1.0");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "ig_inplanecurrent") {
         outputReducer->addOperator(new DRO::DataReductionOperatorIonosphereElement("ig_inplanecurrent", [](
                     SBC::SphericalTriGrid& grid)->std::vector<Real> {

                     std::vector<Real> retval(grid.elements.size()*3);

                     for(uint i=0; i<grid.elements.size(); i++) {
                        // Get effective sigma tensor for this element
                        std::array<Real, 9> sigma = grid.sigmaAverage(i);

                        // Calculate E from element basis functions
                        std::array<Real, 3> E({0,0,0});
                        const std::array<Real, 3>& c1 = grid.nodes[grid.elements[i].corners[0]].x;
                        const std::array<Real, 3>& c2 = grid.nodes[grid.elements[i].corners[1]].x;
                        const std::array<Real, 3>& c3 = grid.nodes[grid.elements[i].corners[2]].x;

                        // ET contains the test function gradient vectors (normalized to potential 1)
                        std::array<std::array<Real,3>, 3> ET({grid.computeGradT(c2,c3,c1), grid.computeGradT(c3,c1,c2), grid.computeGradT(c1,c2,c3)});
                        for(int n=0; n<3; n++) {
                           // Multiply with the corresponding node potentials to get E vector
                           ET[0][n] *= -grid.nodes[grid.elements[i].corners[0]].parameters[ionosphereParameters::SOLUTION];
                           ET[1][n] *= -grid.nodes[grid.elements[i].corners[1]].parameters[ionosphereParameters::SOLUTION];
                           ET[2][n] *= -grid.nodes[grid.elements[i].corners[2]].parameters[ionosphereParameters::SOLUTION];
                        }
                        // Sum up element gradient functions to yield complete E inside this element.
                        for(int n=0; n<3; n++) {
                           E[n] = ET[0][n] + ET[1][n] + ET[2][n];
                        }

                        // Get J from Ohm's law (J = sigma * E)
                        for(int n=0; n<3; n++) {
                           for(int m=0; m<3; m++) {
                              retval[3*i + n] += sigma[3*n+m] * E[m];
                           }
                        }
                     }

                     return retval;
                     }));
         outputReducer->addMetadata(outputReducer->size()-1, "A/m^2", "$\\mathrm{A/m}^2$", "$J$", "1.0");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "ig_upmappedarea") {
         outputReducer->addOperator(new DRO::DataReductionOperatorIonosphereElement("ig_upmappedarea", [](
                     SBC::SphericalTriGrid& grid)->std::vector<Real> {

                     std::vector<Real> retval(grid.elements.size()*3);

                     for(uint i=0; i<grid.elements.size(); i++) {
                        std::array<Real, 3> area = grid.mappedElementArea(i);
                        retval[3*i] = area[0];
                        retval[3*i+1] = area[1];
                        retval[3*i+2] = area[2];
                     }

                     return retval;
                     }));
         outputReducer->addMetadata(outputReducer->size()-1, "m^2", "$\\mathrm{m}^2$", "$A_m$", "1.0");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "ig_sigmap") {
         outputReducer->addOperator(new DRO::DataReductionOperatorIonosphereNode("ig_sigmap", [](
                     SBC::SphericalTriGrid& grid)->std::vector<Real> {

                     std::vector<Real> retval(grid.nodes.size());

                     for(uint i=0; i<grid.nodes.size(); i++) {
                        retval[i] = grid.nodes[i].parameters[ionosphereParameters::SIGMAP];
                     }

                     return retval;
                     }));
         outputReducer->addMetadata(outputReducer->size()-1, "mho", "$\\mathrm{\\Omega^{-1}}$", "$\\Sigma_P$", "1.0");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "ig_sigmah") {
         outputReducer->addOperator(new DRO::DataReductionOperatorIonosphereNode("ig_sigmah", [](
                     SBC::SphericalTriGrid& grid)->std::vector<Real> {

                     std::vector<Real> retval(grid.nodes.size());

                     for(uint i=0; i<grid.nodes.size(); i++) {
                        retval[i] = grid.nodes[i].parameters[ionosphereParameters::SIGMAH];
                     }

                     return retval;
                     }));
         outputReducer->addMetadata(outputReducer->size()-1, "mho", "$\\mathrm{\\Omega^{-1}}$", "$\\Sigma_H$", "1.0");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "ig_sigmaparallel") {
         outputReducer->addOperator(new DRO::DataReductionOperatorIonosphereNode("ig_sigmaparallel", [](
                     SBC::SphericalTriGrid& grid)->std::vector<Real> {

                     std::vector<Real> retval(grid.nodes.size());

                     for(uint i=0; i<grid.nodes.size(); i++) {
                        retval[i] = grid.nodes[i].parameters[ionosphereParameters::SIGMAPARALLEL];
                     }

                     return retval;
                     }));
         outputReducer->addMetadata(outputReducer->size()-1, "mho", "$\\mathrm{\\Omega^{-1}}$", "$\\Sigma_\\parallel$", "1.0");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "ig_rhon") {
         outputReducer->addOperator(new DRO::DataReductionOperatorIonosphereNode("ig_rhon", [](
                     SBC::SphericalTriGrid& grid)->std::vector<Real> {

                     std::vector<Real> retval(grid.nodes.size());

                     for(uint i=0; i<grid.nodes.size(); i++) {
                        retval[i] = grid.nodes[i].parameters[ionosphereParameters::RHON];
                     }

                     return retval;
                     }));
         outputReducer->addMetadata(outputReducer->size()-1, "m^-3", "$\\mathrm{m^{-3}}$", "$\\n_e$", "1.0");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "ig_electrontemp") {
         outputReducer->addOperator(new DRO::DataReductionOperatorIonosphereNode("ig_electrontemp", [](
                     SBC::SphericalTriGrid& grid)->std::vector<Real> {

                     std::vector<Real> retval(grid.nodes.size());

                     for(uint i=0; i<grid.nodes.size(); i++) {
                        retval[i] = grid.nodes[i].electronTemperature();
                     }

                     return retval;
                     }));
         outputReducer->addMetadata(outputReducer->size()-1, "K", "$\\mathrm{K}$", "$T_e$", "1.0");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "ig_deltaphi") {
         outputReducer->addOperator(new DRO::DataReductionOperatorIonosphereNode("ig_deltaphi", [](
                     SBC::SphericalTriGrid& grid)->std::vector<Real> {

                     std::vector<Real> retval(grid.nodes.size());

                     for(uint i=0; i<grid.nodes.size(); i++) {
                        retval[i] = grid.nodes[i].deltaPhi();
                     }

                     return retval;
                     }));
         outputReducer->addMetadata(outputReducer->size()-1, "eV", "$\\mathrm{eV}$", "$\\Delta\\Phi$", "1.0");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "ig_precipitation") {
         outputReducer->addOperator(new DRO::DataReductionOperatorIonosphereNode("ig_precipitation", [](
                     SBC::SphericalTriGrid& grid)->std::vector<Real> {

                     std::vector<Real> retval(grid.nodes.size());

                     for(uint i=0; i<grid.nodes.size(); i++) {
                        retval[i] = grid.nodes[i].parameters[ionosphereParameters::PRECIP];
                     }

                     return retval;
                     }));
         outputReducer->addMetadata(outputReducer->size()-1, "W/m^2", "$\\mathrm{W m^{-2}}$", "$W_\\mathrm{precipitation}$", "1.0");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "ig_precipnumflux") {
         outputReducer->addOperator(new DRO::DataReductionOperatorIonosphereNode("ig_precipnumflux", [](SBC::SphericalTriGrid& grid)->std::vector<Real> {

                     std::array< Real, SBC::productionNumParticleEnergies+1 > particle_energy;
                     // Precalculate effective energy bins
                     // Make sure this stays in sync with sysboundary/ionosphere.cpp
                     for(int e=0; e<SBC::productionNumParticleEnergies; e++) {
                     particle_energy[e] = pow(10.0, -1.+e*(2.3+1.)/(SBC::productionNumParticleEnergies-1));
                     }
                     particle_energy[SBC::productionNumParticleEnergies] = 2*particle_energy[SBC::productionNumParticleEnergies-1] - particle_energy[SBC::productionNumParticleEnergies-2];

                     Real accenergy = SBC::productionMinAccEnergy;

                     std::vector<Real> retval(grid.nodes.size());
                     for(uint i=0; i<grid.nodes.size(); i++) {
                        Real temp_keV = physicalconstants::K_B * grid.nodes[i].electronTemperature() / physicalconstants::CHARGE / 1000;

                        for(int p=0; p<SBC::productionNumParticleEnergies; p++) {
                           Real energyparam = (particle_energy[p]-accenergy)/temp_keV; // = E_p / (kB T)
                           Real deltaE = (particle_energy[p+1] - particle_energy[p])* 1e3*physicalconstants::CHARGE;  // dE in J
                           retval[i] += grid.nodes[i].parameters[ionosphereParameters::RHON] * sqrt(1. / (2. * M_PI * physicalconstants::MASS_ELECTRON))
                           * particle_energy[p] / temp_keV / sqrt(temp_keV * 1e3 *physicalconstants::CHARGE)
                           * deltaE * exp(-energyparam); // Flux 1/m^2/s
                        }
                     }
                     return retval;
                     }));
         outputReducer->addMetadata(outputReducer->size()-1, "1/m^2/s", "$m^{-2} s^{-1}$", "$\\bar{F}_\\mathrm{precip}$", "1.0");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "ig_precipavgenergy") {
         outputReducer->addOperator(new DRO::DataReductionOperatorIonosphereNode("ig_precipavgenergy", [](SBC::SphericalTriGrid& grid)->std::vector<Real> {

                     std::array< Real, SBC::productionNumParticleEnergies+1 > particle_energy;
                     // Precalculate effective energy bins
                     // Make sure this stays in sync with sysboundary/ionosphere.cpp
                     for(int e=0; e<SBC::productionNumParticleEnergies; e++) {
                     particle_energy[e] = pow(10.0, -1.+e*(2.3+1.)/(SBC::productionNumParticleEnergies-1));
                     }
                     particle_energy[SBC::productionNumParticleEnergies] = 2*particle_energy[SBC::productionNumParticleEnergies-1] - particle_energy[SBC::productionNumParticleEnergies-2];

                     Real accenergy = SBC::productionMinAccEnergy;

                     std::vector<Real> retval(grid.nodes.size());
                     for(uint i=0; i<grid.nodes.size(); i++) {
                        Real numberFlux = 0;

                        // Calculate precipitating number flux at this node
                        // (TODO: this is completely copy'n'pasted from the
                        // ig_precipnumflux reducer above. Share code?)
                        Real temp_keV = physicalconstants::K_B * grid.nodes[i].electronTemperature() / physicalconstants::CHARGE / 1000;

                        for(int p=0; p<SBC::productionNumParticleEnergies; p++) {
                           Real energyparam = (particle_energy[p]-accenergy)/temp_keV; // = E_p / (kB T)
                           Real deltaE = (particle_energy[p+1] - particle_energy[p])* 1e3*physicalconstants::CHARGE;  // dE in J
                           numberFlux += grid.nodes[i].parameters[ionosphereParameters::RHON] * sqrt(1. / (2. * M_PI * physicalconstants::MASS_ELECTRON))
                           * particle_energy[p] / temp_keV / sqrt(temp_keV * 1e3 *physicalconstants::CHARGE)
                           * deltaE * exp(-energyparam); // Flux 1/m^2/s
                        }

                        // Average precipitating energy = energyFlux / numberFlux (in eV)
                        retval[i] = grid.nodes[i].parameters[ionosphereParameters::PRECIP] / numberFlux / physicalconstants::CHARGE;
                     }
                     return retval;
                     }));
         outputReducer->addMetadata(outputReducer->size()-1, "eV", "eV", "$\\bar{E}_\\mathrm{precip}$", "1.0");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "ig_potential") {
         outputReducer->addOperator(new DRO::DataReductionOperatorIonosphereNode("ig_potential", [](
                     SBC::SphericalTriGrid& grid)->std::vector<Real> {

                     std::vector<Real> retval(grid.nodes.size());

                     for(uint i=0; i<grid.nodes.size(); i++) {
                        retval[i] = grid.nodes[i].parameters[ionosphereParameters::SOLUTION];
                     }

                     return retval;
                     }));
         outputReducer->addMetadata(outputReducer->size()-1, "V", "$\\mathrm{V}$", "$\\phi_I$", "1.0");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "ig_solverinternals") {
         outputReducer->addOperator(new DRO::DataReductionOperatorIonosphereNode("ig_source", [](
                     SBC::SphericalTriGrid& grid)->std::vector<Real> {

                     std::vector<Real> retval(grid.nodes.size());

                     for(uint i=0; i<grid.nodes.size(); i++) {
                        retval[i] = grid.nodes[i].parameters[ionosphereParameters::SOURCE];
                     }

                     return retval;
                     }));
         outputReducer->addOperator(new DRO::DataReductionOperatorIonosphereNode("ig_residual", [](
                     SBC::SphericalTriGrid& grid)->std::vector<Real> {

                     std::vector<Real> retval(grid.nodes.size());

                     for(uint i=0; i<grid.nodes.size(); i++) {
                        retval[i] = grid.nodes[i].parameters[ionosphereParameters::RESIDUAL];
                     }

                     return retval;
                     }));
         outputReducer->addOperator(new DRO::DataReductionOperatorIonosphereNode("ig_p", [](
                     SBC::SphericalTriGrid& grid)->std::vector<Real> {

                     std::vector<Real> retval(grid.nodes.size());

                     for(uint i=0; i<grid.nodes.size(); i++) {
                        retval[i] = grid.nodes[i].parameters[ionosphereParameters::PPARAM];
                     }

                     return retval;
                     }));
         outputReducer->addOperator(new DRO::DataReductionOperatorIonosphereNode("ig_pp", [](
                     SBC::SphericalTriGrid& grid)->std::vector<Real> {

                     std::vector<Real> retval(grid.nodes.size());

                     for(uint i=0; i<grid.nodes.size(); i++) {
                        retval[i] = grid.nodes[i].parameters[ionosphereParameters::PPPARAM];
                     }

                     return retval;
                     }));
         outputReducer->addOperator(new DRO::DataReductionOperatorIonosphereNode("ig_z", [](
                     SBC::SphericalTriGrid& grid)->std::vector<Real> {

                     std::vector<Real> retval(grid.nodes.size());

                     for(uint i=0; i<grid.nodes.size(); i++) {
                        retval[i] = grid.nodes[i].parameters[ionosphereParameters::ZPARAM];
                     }

                     return retval;
                     }));
         outputReducer->addOperator(new DRO::DataReductionOperatorIonosphereNode("ig_zz", [](
                     SBC::SphericalTriGrid& grid)->std::vector<Real> {

                     std::vector<Real> retval(grid.nodes.size());

                     for(uint i=0; i<grid.nodes.size(); i++) {
                        retval[i] = grid.nodes[i].parameters[ionosphereParameters::ZZPARAM];
                     }

                     return retval;
                     }));

         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "ig_upmappednodecoords") {
         outputReducer->addOperator(new DRO::DataReductionOperatorIonosphereNode("ig_upmappednodecoords", [](
                     SBC::SphericalTriGrid& grid)->std::vector<Real> {

                     std::vector<Real> retval(grid.nodes.size()*3);

                     for(uint i=0; i<grid.nodes.size(); i++) {
                        retval[3*i] = grid.nodes[i].xMapped[0];
                        retval[3*i+1] = grid.nodes[i].xMapped[1];
                        retval[3*i+2] = grid.nodes[i].xMapped[2];
                     }

                     return retval;
                     }));
         outputReducer->addMetadata(outputReducer->size()-1, "m", "m", "$x_\\mathrm{mapped}$", "1.0");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "ig_upmappedb") {
         outputReducer->addOperator(new DRO::DataReductionOperatorIonosphereNode("ig_upmappedb", [](
                     SBC::SphericalTriGrid& grid)->std::vector<Real> {

                     std::vector<Real> retval(grid.nodes.size()*3);

                     for(uint i=0; i<grid.nodes.size(); i++) {
                        retval[3*i] = grid.nodes[i].parameters[ionosphereParameters::UPMAPPED_BX];
                        retval[3*i+1] = grid.nodes[i].parameters[ionosphereParameters::UPMAPPED_BY];
                        retval[3*i+2] = grid.nodes[i].parameters[ionosphereParameters::UPMAPPED_BZ];
                     }

                     return retval;
                     }));
         outputReducer->addMetadata(outputReducer->size()-1, "T", "T", "$B_\\mathrm{mapped}$", "1.0");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "ig_openclosed") {
         FieldTracing::fieldTracingParameters.doTraceOpenClosed = true;
         outputReducer->addOperator(new DRO::DataReductionOperatorIonosphereNode("ig_openclosed", [](
            SBC::SphericalTriGrid& grid)->std::vector<Real> {

               std::vector<Real> retval(grid.nodes.size());

               for(uint i=0; i<grid.nodes.size(); i++) {
                  retval[i] = (Real)grid.nodes[i].openFieldLine;
               }

               return retval;
            }));
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "ig_fac") {
         outputReducer->addOperator(new DRO::DataReductionOperatorIonosphereNode("ig_fac", [](
                     SBC::SphericalTriGrid& grid)->std::vector<Real> {

                     std::vector<Real> retval(grid.nodes.size());

                     for(uint i=0; i<grid.nodes.size(); i++) {
                        Real area = 0;
                        for(uint e=0; e<grid.nodes[i].numTouchingElements; e++) {
                           area += grid.elementArea(grid.nodes[i].touchingElements[e]);
                        }
                        area /= 3.; // As every element has 3 corners, don't double-count areas
                        retval[i] = grid.nodes[i].parameters[ionosphereParameters::SOURCE]/area;
                     }

                     return retval;
                     }));
         outputReducer->addMetadata(outputReducer->size()-1, "A/m^2", "$\\mathrm{A m}^{-2}$", "$I_\\mathrm{FAC}$", "1.0");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "vg_ionospherecoupling") {
         outputReducer->addOperator(new DRO::DataReductionOperatorMPIGridCell("vg_ionospherecoupling", 3, [](
                     const SpatialCell* cell)->std::vector<Real> {

                  std::vector<Real> retval(3);

                  // Just return a 0,0,0 vector for non-ionosphere cells
                  if(cell->sysBoundaryFlag != sysboundarytype::IONOSPHERE) {
                     retval[0]=0;
                     retval[1]=0;
                     retval[2]=0;
                     return retval;
                  }

                  std::array<Real, 3> x;
                  x[0] = cell->parameters[CellParams::XCRD] + cell->parameters[CellParams::DX];
                  x[1] = cell->parameters[CellParams::YCRD] + cell->parameters[CellParams::DY];
                  x[2] = cell->parameters[CellParams::ZCRD] + cell->parameters[CellParams::DZ];

                  std::array<std::pair<int, Real>, 3> coupling = FieldTracing::calculateIonosphereVlasovGridCoupling(x, SBC::ionosphereGrid.nodes, SBC::Ionosphere::radius);
                  for(int i=0; i<3; i++) {
                     uint coupledNode = coupling[i].first;
                     Real a = coupling[i].second;
                     retval[0] += a*SBC::ionosphereGrid.nodes[coupledNode].x[0] - (cell->parameters[CellParams::XCRD] + 0.5*cell->parameters[CellParams::DX]);
                     retval[1] = a*SBC::ionosphereGrid.nodes[coupledNode].x[1] - (cell->parameters[CellParams::YCRD] + 0.5*cell->parameters[CellParams::DY]);
                     retval[2] = a*SBC::ionosphereGrid.nodes[coupledNode].x[2] - (cell->parameters[CellParams::ZCRD] + 0.5*cell->parameters[CellParams::DZ]);
                  }

                  return retval;
               }));
         outputReducer->addMetadata(outputReducer->size()-1, "m", "m", "$x_\\mathrm{coupled}$", "1.0");
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "vg_connection") {
         FieldTracing::fieldTracingParameters.doTraceFullBox = true;
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("vg_connection",CellParams::CONNECTION,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("vg_connection_coordinates_fw",CellParams::CONNECTION_FW_X,3));
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("vg_connection_coordinates_bw",CellParams::CONNECTION_BW_X,3));
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "vg_fluxrope" || lowercase == "vg_curvature") {
         Parameters::computeCurvature = true;
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("vg_curvature",CellParams::CURVATUREX,3));
         if(P::systemWriteAllDROs || lowercase == "vg_fluxrope") {
            FieldTracing::fieldTracingParameters.doTraceFullBox = true;
            outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("vg_fluxrope",CellParams::FLUXROPE,1));
         }
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs || lowercase == "fg_curvature") {
         Parameters::computeCurvature = true;
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_curvature",[](
            FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
            FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
            FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
            FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
            FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
            FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
            FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
            FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
            FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
            FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid)->std::vector<double> {

               std::array<FsGridTools::FsIndex_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2]*3);

               for(int z=0; z<gridSize[2]; z++) {
                  for(int y=0; y<gridSize[1]; y++) {
                     for(int x=0; x<gridSize[0]; x++) {
                        retval[3*(gridSize[1]*gridSize[0]*z + gridSize[0]*y + x)] =     (*volGrid.get(x,y,z))[fsgrids::volfields::CURVATUREX];
                        retval[3*(gridSize[1]*gridSize[0]*z + gridSize[0]*y + x) + 1] = (*volGrid.get(x,y,z))[fsgrids::volfields::CURVATUREY];
                        retval[3*(gridSize[1]*gridSize[0]*z + gridSize[0]*y + x) + 2] = (*volGrid.get(x,y,z))[fsgrids::volfields::CURVATUREZ];
                     }
                  }
               }
               return retval;
            }
         ));
         if(!P::systemWriteAllDROs) {
            continue;
         }
      }
      if(P::systemWriteAllDROs) {
         break; // from the loop
      }
      // After all the continue; statements one should never land here.
      int myRank;
      MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
      if (myRank == MASTER_RANK) {
         std::cerr << __FILE__ << ":" << __LINE__ << ": The output variable " << *it << " is not defined." << std::endl;
      }
      MPI_Finalize();
      exit(1);
   }

   for (it = P::diagnosticVariableList.begin();
        it != P::diagnosticVariableList.end();
        it++) {

      // Sidestep mixed case errors
      std::string lowercase = *it;
      for(auto& c : lowercase) c = tolower(c);

      if(P::diagnosticWriteAllDROs || lowercase == "populations_blocks" || lowercase == "populations_vg_blocks") {
         // Per-population total block counts
         for(unsigned int i =0; i < getObjectWrapper().particleSpecies.size(); i++) {
            diagnosticReducer->addOperator(new DRO::Blocks(i));
         }
         if(!P::diagnosticWriteAllDROs) {
            continue;
         }
      }
      if(P::diagnosticWriteAllDROs || lowercase == "vg_rhom" || lowercase == "rhom") {
         // Overall mass density
         diagnosticReducer->addOperator(new DRO::DataReductionOperatorCellParams("vg_rhom",CellParams::RHOM,1));
         if(!P::diagnosticWriteAllDROs) {
            continue;
         }
      }
      if(P::diagnosticWriteAllDROs || lowercase == "populations_rholossadjust" || lowercase == "populations_rho_loss_adjust" || lowercase == "populations_vg_rho_loss_adjust") {
         // Per-particle overall lost particle number
         for(unsigned int i =0; i < getObjectWrapper().particleSpecies.size(); i++) {
            species::Species& species=getObjectWrapper().particleSpecies[i];
            const std::string& pop = species.name;
            diagnosticReducer->addOperator(new DRO::DataReductionOperatorPopulations<Real>(pop + "/vg_rho_loss_adjust", i, offsetof(spatial_cell::Population, RHOLOSSADJUST), 1));
         }
         if(!P::diagnosticWriteAllDROs) {
            continue;
         }
      }
      if(P::diagnosticWriteAllDROs || lowercase == "lbweight" || lowercase == "vg_lbweight" || lowercase == "vg_loadbalanceweight" || lowercase == "vg_loadbalance_weight" || lowercase == "loadbalance_weight") {
         diagnosticReducer->addOperator(new DRO::DataReductionOperatorCellParams("vg_loadbalance_weight",CellParams::LBWEIGHTCOUNTER,1));
         if(!P::diagnosticWriteAllDROs) {
            continue;
         }
      }
      if(P::diagnosticWriteAllDROs || lowercase == "maxvdt" || lowercase == "maxdt_acceleration" || lowercase == "vg_maxdt_acceleration") {
         diagnosticReducer->addOperator(new DRO::DataReductionOperatorCellParams("vg_maxdt_acceleration",CellParams::MAXVDT,1));
         if(!P::diagnosticWriteAllDROs) {
            continue;
         }
      }
      if(P::diagnosticWriteAllDROs || lowercase == "maxrdt" || lowercase == "maxdt_translation" || lowercase == "vg_maxdt_translation") {
         diagnosticReducer->addOperator(new DRO::DataReductionOperatorCellParams("vg_maxdt_translation",CellParams::MAXRDT,1));
         if(!P::diagnosticWriteAllDROs) {
            continue;
         }
      }
      if(P::diagnosticWriteAllDROs || lowercase == "maxfieldsdt" || lowercase == "maxdt_fieldsolver" || lowercase == "fg_maxfieldsdt" || lowercase == "fg_maxdt_fieldsolver") {
         diagnosticReducer->addOperator(new DRO::DataReductionOperatorCellParams("fg_maxdt_fieldsolver",CellParams::MAXFDT,1));
         if(!P::diagnosticWriteAllDROs) {
            continue;
         }
      }
      if(P::diagnosticWriteAllDROs || lowercase == "populations_maxdistributionfunction" || lowercase == "populations_vg_maxdistributionfunction") {
         for(unsigned int i =0; i < getObjectWrapper().particleSpecies.size(); i++) {
            diagnosticReducer->addOperator(new DRO::MaxDistributionFunction(i));
         }
         if(!P::diagnosticWriteAllDROs) {
            continue;
         }
      }
      if(P::diagnosticWriteAllDROs || lowercase == "populations_mindistributionfunction" || lowercase == "populations_vg_mindistributionfunction") {
         for(unsigned int i =0; i < getObjectWrapper().particleSpecies.size(); i++) {
            diagnosticReducer->addOperator(new DRO::MinDistributionFunction(i));
         }
         if(!P::diagnosticWriteAllDROs) {
            continue;
         }
      }
      if(P::diagnosticWriteAllDROs || lowercase == "populations_maxrdt" || lowercase == "populations_maxdt_translation" || lowercase == "populations_vg_maxdt_translation") {
         for(unsigned int i =0; i < getObjectWrapper().particleSpecies.size(); i++) {
            species::Species& species=getObjectWrapper().particleSpecies[i];
            const std::string& pop = species.name;
            diagnosticReducer->addOperator(new DRO::DataReductionOperatorPopulations<Real>(pop + "/vg_maxdt_translation", i, offsetof(spatial_cell::Population, max_dt[0]), 1));
         }
         if(!P::diagnosticWriteAllDROs) {
            continue;
         }
      }
      if(P::diagnosticWriteAllDROs || lowercase == "populations_maxvdt" || lowercase == "populations_maxdt_acceleration" || lowercase == "populations_vg_maxdt_acceleration") {
         for(unsigned int i =0; i < getObjectWrapper().particleSpecies.size(); i++) {
            species::Species& species=getObjectWrapper().particleSpecies[i];
            const std::string& pop = species.name;
            diagnosticReducer->addOperator(new DRO::DataReductionOperatorPopulations<Real>(pop + "/vg_maxdt_acceleration", i, offsetof(spatial_cell::Population, max_dt[1]), 1));
         }
         if(!P::diagnosticWriteAllDROs) {
            continue;
         }
      }
      if(P::diagnosticWriteAllDROs) {
         break; // from the loop
      }
      // After all the continue; statements one should never land here.
      int myRank;
      MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
      if (myRank == MASTER_RANK) {
         std::cerr << __FILE__ << ":" << __LINE__ << ": The diagnostic variable " << *it << " is not defined." << std::endl;
      }
      MPI_Finalize();
      exit(1);
   }
}

// ************************************************************
// ***** DEFINITIONS FOR DATAREDUCER CLASS *****
// ************************************************************

/** Constructor for class DataReducer.
 */
DataReducer::DataReducer() { }

/** Destructor for class DataReducer. All stored DataReductionOperators
 * are deleted.
 */
DataReducer::~DataReducer() {
   // Call delete for each DataReductionOperator:
   for (vector<DRO::DataReductionOperator*>::iterator it=operators.begin(); it!=operators.end(); ++it) {
      delete *it;
      *it = NULL;
   }
}

/** Add a new DRO::DataReductionOperator which has been created with new operation.
 * DataReducer will take care of deleting it.
 * @return If true, the given DRO::DataReductionOperator was added successfully.
 */
bool DataReducer::addOperator(DRO::DataReductionOperator* op) {
   operators.push_back(op);
   return true;
}

/** Get the name of a DataReductionOperator.
 * @param operatorID ID number of the operator whose name is requested.
 * @return Name of the operator.
 */
std::string DataReducer::getName(const unsigned int& operatorID) const {
   if (operatorID >= operators.size()) return "";
   return operators[operatorID]->getName();
}

/** Get info on the type of data calculated by the given DataReductionOperator.
 * A DataReductionOperator writes an array on disk. Each element of the array is a vector with n elements. Finally, each
 * vector element has a byte size, as given by the sizeof function.
 * @param operatorID ID number of the DataReductionOperator whose output data info is requested.
 * @param dataType Basic datatype, must be int, uint, or float.
 * @param dataSize Byte size of written datatype, for example double-precision floating points
 * have byte size of sizeof(double).
 * @param vectorSize How many elements are in the vector returned by the DataReductionOperator.
 * @return If true, DataReductionOperator was found and it returned sensible values.
 */
bool DataReducer::getDataVectorInfo(const unsigned int& operatorID,std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
   if (operatorID >= operators.size()) return false;
   return operators[operatorID]->getDataVectorInfo(dataType,dataSize,vectorSize);
}

/** Add a metadata to the specified DRO::DataReductionOperator.
 * @param operatorID ID number of the DataReductionOperator to add metadata to
 * @param unit string with the physical unit of the DRO result
 * @param unitLaTeX LaTeX-formatted string with the physical unit of the DRO result
 * @param variableLaTeX LaTeX-formatted string with a descriptive short name for the DRO result
 * @param conversionFactor floating point conversion factor between DRO result and SI units
 * @return If true, the given metadata  was added successfully.
 */
bool DataReducer::addMetadata(const unsigned int operatorID, std::string unit,std::string unitLaTeX,std::string variableLaTeX,std::string unitConversion) {
   if (operatorID >= operators.size()) return false;
   return operators[operatorID]->setUnitMetadata(unit,unitLaTeX,variableLaTeX,unitConversion);
}

/** Get metadata on the unit of data calculated by the given DataReductionOperator.
 * @param operatorID ID number of the DataReductionOperator whose output unit metadata is requested.
 * @param unit Physical unit of variable
 * @param unitLaTeX Physical unit of variable, written using LaTeX notation
 * @param unitConversion Floating point value of conversion factor to SI units
 * @return If true, DataReductionOperator was found and it returned sensible values.
 */
bool DataReducer::getMetadata(const unsigned int& operatorID,std::string& unit,std::string& unitLaTeX,std::string& variableLaTeX,std::string& unitConversion) const {
   if (operatorID >= operators.size()) return false;
   return operators[operatorID]->getUnitMetadata(unit, unitLaTeX, variableLaTeX, unitConversion);
}

/** Ask a DataReductionOperator if it wants to write parameters to the vlsv file header
 * @param operatorID ID number of the DataReductionOperator.
 * @return If true, then VLSVWriter should be passed to the DataReductionOperator.*/
bool DataReducer::hasParameters(const unsigned int& operatorID) const {
   if (operatorID >= operators.size()) return false;
   return dynamic_cast<DRO::DataReductionOperatorHasParameters*>(operators[operatorID]) != nullptr;
}

/** Request a DataReductionOperator to calculate its output data and to write it to the given buffer.
 * @param cell Pointer to spatial cell whose data is to be reduced.
 * @param operatorID ID number of the applied DataReductionOperator.
 * @param buffer Buffer in which DataReductionOperator should write its data.
 * @return If true, DataReductionOperator calculated and wrote data successfully.
 */
bool DataReducer::reduceData(const SpatialCell* cell,const unsigned int& operatorID,char* buffer) {
   // Tell the chosen operator which spatial cell we are counting:
   if (operatorID >= operators.size()) return false;
   if (operators[operatorID]->setSpatialCell(cell) == false) return false;

   if (operators[operatorID]->reduceData(cell,buffer) == false) return false;
   return true;
}

/** Request a DataReductionOperator to calculate its output data and to write it to the given variable.
 * @param cell Pointer to spatial cell whose data is to be reduced.
 * @param operatorID ID number of the applied DataReductionOperator.
 * @param result Real variable in which DataReductionOperator should write its result.
 * @return If true, DataReductionOperator calculated and wrote data successfully.
 */
bool DataReducer::reduceDiagnostic(const SpatialCell* cell,const unsigned int& operatorID,Real * result) {
   // Tell the chosen operator which spatial cell we are counting:
   if (operatorID >= operators.size()) return false;
   if (operators[operatorID]->setSpatialCell(cell) == false) return false;

   if (operators[operatorID]->reduceDiagnostic(cell,result) == false) return false;
   return true;
}

/** Get the number of DataReductionOperators stored in DataReducer.
 * @return Number of DataReductionOperators stored in DataReducer.
 */
unsigned int DataReducer::size() const {return operators.size();}

/** Write parameters related to given DataReductionOperator to the output file.
 * @param operatorID ID number of the selected DataReductionOperator.
 * @param vlsvWriter VLSV file writer that has output file open.
 * @return If true, DataReductionOperator wrote its parameters successfully.*/
bool DataReducer::writeParameters(const unsigned int& operatorID, vlsv::Writer& vlsvWriter) {
   if (operatorID >= operators.size()) return false;
   DRO::DataReductionOperatorHasParameters* parameterOperator = dynamic_cast<DRO::DataReductionOperatorHasParameters*>(operators[operatorID]);
   if(parameterOperator == nullptr) {
      return false;
   }
   return parameterOperator->writeParameters(vlsvWriter);
}
/** Write all data thet the given DataReductionOperator wants to obtain from fsgrid into the output file.
 */
bool DataReducer::writeFsGridData(
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
                      const std::string& meshName, const unsigned int operatorID,
                      vlsv::Writer& vlsvWriter,
                      const bool writeAsFloat) {

   if (operatorID >= operators.size()) return false;
   DRO::DataReductionOperatorFsGrid* DROf = dynamic_cast<DRO::DataReductionOperatorFsGrid*>(operators[operatorID]);
   if(!DROf) {
      return false;
   } else {
      return DROf->writeFsGridData(perBGrid, EGrid, EHallGrid, EGradPeGrid, momentsGrid, dPerBGrid, dMomentsGrid, BgBGrid, volGrid, technicalGrid, meshName, vlsvWriter, writeAsFloat);
   }
}

bool DataReducer::writeIonosphereGridData(
                     SBC::SphericalTriGrid& grid, const std::string& meshName,
                     const unsigned int operatorID, vlsv::Writer& vlsvWriter) {

   if (operatorID >= operators.size()) return false;
   DRO::DataReductionOperatorIonosphereElement* DROe = dynamic_cast<DRO::DataReductionOperatorIonosphereElement*>(operators[operatorID]);
   DRO::DataReductionOperatorIonosphereNode* DROn = dynamic_cast<DRO::DataReductionOperatorIonosphereNode*>(operators[operatorID]);
   if(DROe) {
      return DROe->writeIonosphereData(grid, vlsvWriter);
   } else if(DROn) {
      return DROn->writeIonosphereData(grid, vlsvWriter);
   } else {
      return false;
   }

}
