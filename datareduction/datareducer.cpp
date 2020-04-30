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
#include <boost/algorithm/string.hpp>    
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
      const std::string lowercase = boost::algorithm::to_lower_copy(*it);
      
      if(lowercase == "fg_b" || lowercase == "b") { // Bulk magnetic field at Yee-Lattice locations
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_b",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2>& perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 2>& EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, 2>& EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, 2>& EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 2>& momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, 2>& dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 2>& dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, 2>& BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, 2>& volGrid,
                      FsGrid< fsgrids::technical, 2>& technicalGrid)->std::vector<double> {

               std::array<int32_t,3>& gridSize = technicalGrid.getLocalSize();
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
         outputReducer->addMetadata(outputReducer->size()-1,"T","$\\mathrm{T}$","$B$","1.0");
         continue;	 
      }
      if(lowercase == "fg_backgroundb" || lowercase == "backgroundb" || lowercase == "fg_b_background") { // Static (typically dipole) magnetic field part
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_b_background",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2>& perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 2>& EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, 2>& EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, 2>& EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 2>& momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, 2>& dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 2>& dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, 2>& BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, 2>& volGrid,
                      FsGrid< fsgrids::technical, 2>& technicalGrid)->std::vector<double> {

               std::array<int32_t,3>& gridSize = technicalGrid.getLocalSize();
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
	 outputReducer->addMetadata(outputReducer->size()-1,"T","$\\mathrm{T}$","$B_\\mathrm{bg}$","1.0");
         continue;
      }
      if(lowercase == "fg_perturbedb" || lowercase == "perturbedb" || lowercase == "fg_b_perturbed") { // Fluctuating magnetic field part
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_b_perturbed",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2>& perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 2>& EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, 2>& EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, 2>& EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 2>& momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, 2>& dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 2>& dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, 2>& BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, 2>& volGrid,
                      FsGrid< fsgrids::technical, 2>& technicalGrid)->std::vector<double> {

               std::array<int32_t,3>& gridSize = technicalGrid.getLocalSize();
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
	 outputReducer->addMetadata(outputReducer->size()-1,"T","$\\mathrm{T}$","$B_\\mathrm{per}$)","1.0");
         continue;
      }
      if(lowercase == "fg_e" || lowercase == "e") { // Bulk electric field at Yee-lattice locations
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_e",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2>& perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 2>& EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, 2>& EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, 2>& EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 2>& momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, 2>& dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 2>& dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, 2>& BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, 2>& volGrid,
                      FsGrid< fsgrids::technical, 2>& technicalGrid)->std::vector<double> {

               std::array<int32_t,3>& gridSize = technicalGrid.getLocalSize();
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
         continue;
      }
      if(lowercase == "vg_rhom" || lowercase == "rhom") { // Overall mass density (summed over all populations)
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("vg_rhom",CellParams::RHOM,1));
	 outputReducer->addMetadata(outputReducer->size()-1,"kg/m^3","$\\mathrm{kg}\\,\\mathrm{m}^{-3}$","$\\rho_\\mathrm{m}$","1.0");
         continue;
      }
      if(lowercase == "fg_rhom") { // Overall mass density (summed over all populations)
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_rhom",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2>& perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 2>& EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, 2>& EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, 2>& EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 2>& momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, 2>& dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 2>& dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, 2>& BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, 2>& volGrid,
                      FsGrid< fsgrids::technical, 2>& technicalGrid)->std::vector<double> {

               std::array<int32_t,3>& gridSize = technicalGrid.getLocalSize();
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
         continue;
      }
      if(lowercase == "vg_rhoq" || lowercase == "rhoq") { // Overall charge density (summed over all populations)
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("vg_rhoq",CellParams::RHOQ,1));
	 outputReducer->addMetadata(outputReducer->size()-1,"C/m^3","$\\mathrm{C}\\,\\mathrm{m}^{-3}$","$\\rho_\\mathrm{q}$","1.0");
         continue;
      }
      if(lowercase == "fg_rhoq") { // Overall charge density (summed over all populations)
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_rhoq",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2>& perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 2>& EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, 2>& EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, 2>& EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 2>& momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, 2>& dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 2>& dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, 2>& BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, 2>& volGrid,
                      FsGrid< fsgrids::technical, 2>& technicalGrid)->std::vector<double> {

               std::array<int32_t,3>& gridSize = technicalGrid.getLocalSize();
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
         continue;
      }
      if(lowercase == "populations_rho" || lowercase == "populations_vg_rho") { // Per-population particle number density
         for(unsigned int i =0; i < getObjectWrapper().particleSpecies.size(); i++) {
            species::Species& species=getObjectWrapper().particleSpecies[i];
            const std::string& pop = species.name;
            outputReducer->addOperator(new DRO::DataReductionOperatorPopulations<Real>(pop + "/vg_rho", i, offsetof(spatial_cell::Population, RHO), 1));
	    outputReducer->addMetadata(outputReducer->size()-1,"1/m^3","$\\mathrm{m}^{-3}$","$n_\\mathrm{"+pop+"}$","1.0");
         }
         continue;
      }
      
      if(lowercase == "v" || lowercase == "vg_v") { // Overall effective bulk density defining the center-of-mass frame from all populations
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("vg_v",CellParams::VX,3));
	 outputReducer->addMetadata(outputReducer->size()-1,"m/s","$\\mathrm{m}\\,\\mathrm{s}^{-1}$","$V$","1.0");
         continue;
      }
      if(lowercase == "fg_v") { // Overall effective bulk density defining the center-of-mass frame from all populations
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_v",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2>& perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 2>& EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, 2>& EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, 2>& EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 2>& momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, 2>& dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 2>& dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, 2>& BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, 2>& volGrid,
                      FsGrid< fsgrids::technical, 2>& technicalGrid)->std::vector<double> {

               std::array<int32_t,3>& gridSize = technicalGrid.getLocalSize();
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
         continue;
      }
      if(lowercase == "populations_v" || lowercase == "populations_vg_v") { // Per population bulk velocities
         for(unsigned int i =0; i < getObjectWrapper().particleSpecies.size(); i++) {
            species::Species& species=getObjectWrapper().particleSpecies[i];
            const std::string& pop = species.name;
            outputReducer->addOperator(new DRO::DataReductionOperatorPopulations<Real>(pop + "/vg_v", i, offsetof(spatial_cell::Population, V), 3));
	    outputReducer->addMetadata(outputReducer->size()-1,"m/s","$\\mathrm{m}\\,\\mathrm{s}^{-1}$","$V_\\mathrm{"+pop+"}$","1.0");
         }
         continue;
      }
      if(lowercase == "populations_moments_backstream" || lowercase == "populations_moments_nonthermal" || lowercase == "populations_vg_moments_nonthermal") { // Per-population moments of the backstreaming part
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
         continue;
      }
      if(lowercase == "populations_moments_nonbackstream" || lowercase == "populations_moments_thermal" || lowercase == "populations_vg_moments_thermal") { // Per-population moments of the non-backstreaming (thermal?) part.
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
         continue;
      }
      if(lowercase == "populations_minvalue" || lowercase == "populations_effectivesparsitythreshold" || lowercase == "populations_vg_effectivesparsitythreshold") {
         // Effective sparsity threshold affecting each cell, if dynamic threshould algorithm is used
         for(unsigned int i =0; i < getObjectWrapper().particleSpecies.size(); i++) {
            species::Species& species=getObjectWrapper().particleSpecies[i];
            const std::string& pop = species.name;
            outputReducer->addOperator(new DRO::VariableEffectiveSparsityThreshold(i));
	    outputReducer->addMetadata(outputReducer->size()-1,"s^3/m^6","$\\mathrm{m}^{-6}\\,\\mathrm{s}^{3}$","$f_\\mathrm{"+pop+",min}$","1.0");
         }
         continue;
      }
      if(lowercase == "populations_rholossadjust" || lowercase == "populations_rho_loss_adjust" || lowercase == "populations_vg_rho_loss_adjust") {
         // Accumulated lost particle number, per population, in each cell, since last restart
         for(unsigned int i =0; i < getObjectWrapper().particleSpecies.size(); i++) {
            species::Species& species=getObjectWrapper().particleSpecies[i];
            const std::string& pop = species.name;
            outputReducer->addOperator(new DRO::DataReductionOperatorPopulations<Real>(pop + "/vg_rho_loss_adjust", i, offsetof(spatial_cell::Population, RHOLOSSADJUST), 1));
	    outputReducer->addMetadata(outputReducer->size()-1,"1/m^3","$\\mathrm{m}^{-3}$","$\\Delta_\\mathrm{loss} n_\\mathrm{"+pop+"}$","1.0");
         }
         continue;
      }
      if(lowercase == "lbweight" || lowercase == "vg_lbweight" || lowercase == "vg_loadbalanceweight" || lowercase == "vg_loadbalance_weight") {
         // Load balance metric for LB debugging
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("vg_loadbalance_weight",CellParams::LBWEIGHTCOUNTER,1));
	 outputReducer->addMetadata(outputReducer->size()-1,"","","$\\mathrm{LB weight}$","");
         continue;
      }
      if(lowercase == "maxvdt" || lowercase == "vg_maxdt_acceleration") {
         // Overall maximum timestep constraint as calculated by the velocity space vlasov update
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("vg_maxdt_acceleration",CellParams::MAXVDT,1));
	 outputReducer->addMetadata(outputReducer->size()-1,"s","$\\mathrm{s}$","$\\Delta t_\\mathrm{V,max}$","1.0");
         continue;
      }
      if(lowercase == "populations_maxvdt" || lowercase == "populations_vg_maxdt_acceleration" || lowercase == "populations_maxdt_acceleration") {
         // Per-population maximum timestep constraint as calculated by the velocity space vlasov update
         for(unsigned int i =0; i < getObjectWrapper().particleSpecies.size(); i++) {
            species::Species& species=getObjectWrapper().particleSpecies[i];
            const std::string& pop = species.name;
            outputReducer->addOperator(new DRO::DataReductionOperatorPopulations<Real>(pop + "/vg_maxdt_acceleration", i, offsetof(spatial_cell::Population, max_dt[1]), 1));
	    outputReducer->addMetadata(outputReducer->size()-1,"s","$\\mathrm{s}$","$\\Delta t_\\mathrm{"+pop+",V,max}$","1.0");
         }
         continue;
      }
      if(lowercase == "maxrdt" || lowercase == "vg_maxdt_translation") {
         // Overall maximum timestep constraint as calculated by the real space vlasov update
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("vg_maxdt_translation",CellParams::MAXRDT,1));
	 outputReducer->addMetadata(outputReducer->size()-1,"s","$\\mathrm{s}$","$\\Delta t_\\mathrm{R,max}$","1.0");
         continue;
      }
      if(lowercase == "populations_maxrdt" || lowercase == "populations_vg_maxdt_translation" || lowercase == "populations_maxdt_translation") {
         // Per-population maximum timestep constraint as calculated by the real space vlasov update
         for(unsigned int i =0; i < getObjectWrapper().particleSpecies.size(); i++) {
            species::Species& species=getObjectWrapper().particleSpecies[i];
            const std::string& pop = species.name;
            outputReducer->addOperator(new DRO::DataReductionOperatorPopulations<Real>(pop + "/vg_maxdt_translation", i, offsetof(spatial_cell::Population, max_dt[0]), 1));
	    outputReducer->addMetadata(outputReducer->size()-1,"s","$\\mathrm{s}$","$\\Delta t_\\mathrm{"+pop+",R,max}$","1.0");
         }
         continue;
      }
      if(lowercase == "populations_energydensity" || lowercase == "populations_vg_energydensity") {
         // Per-population energy density in three energy ranges
         for(unsigned int i =0; i < getObjectWrapper().particleSpecies.size(); i++) {
            species::Species& species=getObjectWrapper().particleSpecies[i];
            const std::string& pop = species.name;
            outputReducer->addOperator(new DRO::VariableEnergyDensity(i));
	    std::stringstream conversion;
	    conversion << (1.0e-6)/physicalconstants::CHARGE;
	    outputReducer->addMetadata(outputReducer->size()-1,"eV/cm^3","$\\mathrm{eV}\\,\\mathrm{cm}^{-3}$","$U_\\mathrm{"+pop+"}$",conversion.str());
         }
         continue;
      }
      if(lowercase == "populations_precipitationflux" || lowercase == "populations_vg_precipitationdifferentialflux" || lowercase == "populations_precipitationdifferentialflux") {
         // Per-population precipitation differential flux
         for(unsigned int i =0; i < getObjectWrapper().particleSpecies.size(); i++) {
            species::Species& species=getObjectWrapper().particleSpecies[i];
            const std::string& pop = species.name;
            outputReducer->addOperator(new DRO::VariablePrecipitationDiffFlux(i));
	    std::stringstream conversion;
	    conversion << (1.0e-4)*physicalconstants::CHARGE;
	    outputReducer->addMetadata(outputReducer->size()-1,"1/(cm^2 sr s eV)","$\\mathrm{cm}^{-2}\\,\\mathrm{sr}^{-1}\\,\\mathrm{s}^{-1}\\,\\mathrm{eV}^{-1}$","$\\mathcal{F}_\\mathrm{"+pop+"}$",conversion.str());
         }
         continue;
      }
      if(lowercase == "maxfieldsdt" || lowercase == "fg_maxfieldsdt" || lowercase == "fg_maxdt_fieldsolver") {
         // Maximum timestep constraint as calculated by the fieldsolver
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_maxdt_fieldsolver",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2>& perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 2>& EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, 2>& EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, 2>& EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 2>& momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, 2>& dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 2>& dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, 2>& BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, 2>& volGrid,
                      FsGrid< fsgrids::technical, 2>& technicalGrid)->std::vector<double> {

               std::array<int32_t,3>& gridSize = technicalGrid.getLocalSize();
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
         continue;
      }
      if(lowercase == "mpirank" || lowercase == "vg_rank") {
         // Map of spatial decomposition of the DCCRG grid into MPI ranks
         outputReducer->addOperator(new DRO::MPIrank);
	 outputReducer->addMetadata(outputReducer->size()-1,"","","$\\mathrm{MPI rank}$","");
         continue;
      }
      if(lowercase == "fsgridrank" || lowercase == "fg_rank") {
         // Map of spatial decomposition of the FsGrid into MPI ranks
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_rank",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2>& perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 2>& EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, 2>& EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, 2>& EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 2>& momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, 2>& dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 2>& dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, 2>& BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, 2>& volGrid,
                      FsGrid< fsgrids::technical, 2>& technicalGrid)->std::vector<double> {

               std::array<int32_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2],technicalGrid.getRank());
               return retval;
             }
         ));
	 outputReducer->addMetadata(outputReducer->size()-1,"","","$\\mathrm{fGrid rank}$","");
         continue;
      }
      if(lowercase == "boundarytype" || lowercase == "vg_boundarytype") {
         // Type of boundarycells
         outputReducer->addOperator(new DRO::BoundaryType);
	 outputReducer->addMetadata(outputReducer->size()-1,"","","$\\mathrm{vGrid Boundary type}$","");
         continue;
      }
      if(lowercase == "fsgridboundarytype" || lowercase == "fg_boundarytype") {
         // Type of boundarycells as stored in FSGrid
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_boundarytype",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2>& perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 2>& EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, 2>& EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, 2>& EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 2>& momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, 2>& dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 2>& dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, 2>& BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, 2>& volGrid,
                      FsGrid< fsgrids::technical, 2>& technicalGrid)->std::vector<double> {

               std::array<int32_t,3>& gridSize = technicalGrid.getLocalSize();
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
         continue;
      }
      if(lowercase == "boundarylayer" || lowercase == "vg_boundarylayer") {
         // For boundaries with multiple layers: layer count per cell
         outputReducer->addOperator(new DRO::BoundaryLayer);
	 outputReducer->addMetadata(outputReducer->size()-1,"","","$\\mathrm{vGrid Boundary layer}$","");
         continue;
      }
      if(lowercase == "fsgridboundarylayer" || lowercase == "fg_boundarylayer") {
         // Type of boundarycells as stored in FSGrid
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_boundarylayer",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2>& perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 2>& EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, 2>& EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, 2>& EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 2>& momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, 2>& dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 2>& dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, 2>& BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, 2>& volGrid,
                      FsGrid< fsgrids::technical, 2>& technicalGrid)->std::vector<double> {

               std::array<int32_t,3>& gridSize = technicalGrid.getLocalSize();
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
         continue;
      }
      if(lowercase == "populations_blocks" || lowercase == "populations_vg_blocks") {
         // Per-population velocity space block counts
         for(unsigned int i =0; i < getObjectWrapper().particleSpecies.size(); i++) {
            species::Species& species=getObjectWrapper().particleSpecies[i];
            const std::string& pop = species.name;
            outputReducer->addOperator(new DRO::Blocks(i));
	    outputReducer->addMetadata(outputReducer->size()-1,"","","$\\mathrm{"+pop+" blocks}$","");
         }
         continue;
      }
      if(lowercase == "fsaved" || lowercase == "vg_fsaved" || lowercase == "vg_f_saved") {
         // Boolean marker whether a velocity space is saved in a given spatial cell
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("vg_f_saved",CellParams::ISCELLSAVINGF,1));
	 outputReducer->addMetadata(outputReducer->size()-1,"","","$f(v)_\\mathrm{saved}$","");
         continue;
      }
      if(lowercase == "populations_accsubcycles" || lowercase == "populations_acceleration_subcycles" || lowercase == "populations_vg_acceleration_subcycles") {
         // Per-population number of subcycles performed for velocity space update
         for(unsigned int i =0; i < getObjectWrapper().particleSpecies.size(); i++) {
            species::Species& species=getObjectWrapper().particleSpecies[i];
            const std::string& pop = species.name;
            outputReducer->addOperator(new DRO::DataReductionOperatorPopulations<uint>(pop + "/vg_acceleration_subcycles", i, offsetof(spatial_cell::Population, ACCSUBCYCLES), 1));
	    outputReducer->addMetadata(outputReducer->size()-1,"","","$\\mathrm{"+pop+" Acc subcycles}$","");
         }
         continue;
      }
      if(lowercase == "vole" || lowercase == "vg_vole" || lowercase == "evol" || lowercase == "vg_e_vol" || lowercase == "e_vol") {
         // Volume-averaged E field
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("vg_e_vol",CellParams::EXVOL,3));
	 outputReducer->addMetadata(outputReducer->size()-1,"V/m","$\\mathrm{V}\\,\\mathrm{m}^{-1}$","$E_\\mathrm{vol,vg}$","1.0");
         continue;
      }
      if(lowercase == "fg_vole" || lowercase == "fg_e_vol" || lowercase == "fg_evol") {
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_e_vol",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2>& perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 2>& EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, 2>& EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, 2>& EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 2>& momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, 2>& dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 2>& dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, 2>& BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, 2>& volGrid,
                      FsGrid< fsgrids::technical, 2>& technicalGrid)->std::vector<double> {

               std::array<int32_t,3>& gridSize = technicalGrid.getLocalSize();
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
         continue;
      }
      if(lowercase == "halle" || lowercase == "fg_halle" || lowercase == "fg_e_hall") {
         for(int index=0; index<fsgrids::N_EHALL; index++) {
            std::string reducer_name = "fg_e_hall_" + std::to_string(index);
            outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid(reducer_name,[index](
                         FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2>& perBGrid,
                         FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 2>& EGrid,
                         FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, 2>& EHallGrid,
                         FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, 2>& EGradPeGrid,
                         FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 2>& momentsGrid,
                         FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, 2>& dPerBGrid,
                         FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 2>& dMomentsGrid,
                         FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, 2>& BgBGrid,
                         FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, 2>& volGrid,
                         FsGrid< fsgrids::technical, 2>& technicalGrid)->std::vector<double> {

                  std::array<int32_t,3>& gridSize = technicalGrid.getLocalSize();
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
         continue;
      }
      if(lowercase =="gradpee" || lowercase == "e_gradpe" || lowercase == "vg_e_gradpe") {
         // Electron pressure gradient contribution to the generalized ohm's law
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("vg_e_gradpe",CellParams::EXGRADPE,3));
	 outputReducer->addMetadata(outputReducer->size()-1,"V/m","$\\mathrm{V}\\,\\mathrm{m}^{-1}$","$E_{\\del P_\\mathrm{e}}$","1.0");
         continue;
      }
      if(lowercase == "volb" || lowercase == "vg_volb" || lowercase == "b_vol" || lowercase == "bvol" || lowercase == "vg_bvol" || lowercase == "vg_b_vol") {
         // Volume-averaged magnetic field
         outputReducer->addOperator(new DRO::VariableBVol);
	 outputReducer->addMetadata(outputReducer->size()-1,"T","$\\mathrm{T}$","$B_\\mathrm{vol,vg}$","1.0");
         continue;
      }
      if(lowercase == "fg_volb" || lowercase == "fg_bvol" || lowercase == "fg_b_vol") { // Static (typically dipole) magnetic field part
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_b_vol",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2>& perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 2>& EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, 2>& EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, 2>& EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 2>& momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, 2>& dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 2>& dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, 2>& BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, 2>& volGrid,
                      FsGrid< fsgrids::technical, 2>& technicalGrid)->std::vector<double> {

               std::array<int32_t,3>& gridSize = technicalGrid.getLocalSize();
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
         continue;
      }
      if(lowercase == "backgroundvolb" || lowercase == "vg_b_background_vol") {
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("vg_b_background_vol",CellParams::BGBXVOL,3));
	 outputReducer->addMetadata(outputReducer->size()-1,"T","$\\mathrm{T}$","$B_\\mathrm{vol,vg,bg}$","1.0");
         continue;
      }
      if(lowercase == "perturbedvolb" || lowercase == "vg_b_perturbed_vol") {
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("vg_b_perturbed_vol",CellParams::PERBXVOL,3));
	 outputReducer->addMetadata(outputReducer->size()-1,"T","$\\mathrm{T}$","$B_\\mathrm{vol,vg,per}$","1.0");
         continue;
      }
      if(lowercase == "pressure" || lowercase == "vg_pressure") {
         // Overall scalar pressure from all populations
         outputReducer->addOperator(new DRO::VariablePressureSolver);
	 outputReducer->addMetadata(outputReducer->size()-1,"Pa","$\\mathrm{Pa}$","$P_\\mathrm{solver}$","1.0");
         continue;
      }
      if(lowercase == "fg_pressure") {
         // Overall scalar pressure from all populations
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_pressure",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2>& perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 2>& EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, 2>& EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, 2>& EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 2>& momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, 2>& dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 2>& dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, 2>& BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, 2>& volGrid,
                      FsGrid< fsgrids::technical, 2>& technicalGrid)->std::vector<double> {

               std::array<int32_t,3>& gridSize = technicalGrid.getLocalSize();
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
         continue;
      }
      if(lowercase == "populations_ptensor" || lowercase == "populations_vg_ptensor") {
         // Per-population pressure tensor, stored as diagonal and offdiagonal components
         for(unsigned int i =0; i < getObjectWrapper().particleSpecies.size(); i++) {
            species::Species& species=getObjectWrapper().particleSpecies[i];
            const std::string& pop = species.name;
            outputReducer->addOperator(new DRO::VariablePTensorDiagonal(i));
	    outputReducer->addMetadata(outputReducer->size()-1,"Pa","$\\mathrm{Pa}$","$\\mathcal{P}_\\mathrm{"+pop+"}$","1.0");
            outputReducer->addOperator(new DRO::VariablePTensorOffDiagonal(i));
	    outputReducer->addMetadata(outputReducer->size()-1,"Pa","$\\mathrm{Pa}$","$\\mathcal{\\tilde{P}}_\\mathrm{"+pop+"}$","1.0");
         }
         continue;
      }
      if(lowercase == "bvolderivs" || lowercase == "b_vol_derivs" || lowercase == "b_vol_derivatives" || lowercase == "vg_b_vol_derivatives" || lowercase == "derivs" ) {
         // Volume-averaged derivatives
         outputReducer->addOperator(new DRO::DataReductionOperatorBVOLDerivatives("vg_dperbxvoldy",bvolderivatives::dPERBXVOLdy,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorBVOLDerivatives("vg_dperbxvoldz",bvolderivatives::dPERBXVOLdz,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorBVOLDerivatives("vg_dperbyvoldx",bvolderivatives::dPERBYVOLdx,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorBVOLDerivatives("vg_dperbyvoldz",bvolderivatives::dPERBYVOLdz,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorBVOLDerivatives("vg_dperbzvoldx",bvolderivatives::dPERBZVOLdx,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorBVOLDerivatives("vg_dperbzvoldy",bvolderivatives::dPERBZVOLdy,1));
	 outputReducer->addMetadata(outputReducer->size()-6,"T/m","$\\mathrm{T}\\,\\mathrm{m}^{-1}$","$\\Delta B_{X,\\mathrm{per,vol,vg}} (\\Delta Y)^{-1}$","1.0");
	 outputReducer->addMetadata(outputReducer->size()-5,"T/m","$\\mathrm{T}\\,\\mathrm{m}^{-1}$","$\\Delta B_{X,\\mathrm{per,vol,vg}} (\\Delta Z)^{-1}$","1.0");
	 outputReducer->addMetadata(outputReducer->size()-4,"T/m","$\\mathrm{T}\\,\\mathrm{m}^{-1}$","$\\Delta B_{Y,\\mathrm{per,vol,vg}} (\\Delta X)^{-1}$","1.0");
	 outputReducer->addMetadata(outputReducer->size()-3,"T/m","$\\mathrm{T}\\,\\mathrm{m}^{-1}$","$\\Delta B_{Y,\\mathrm{per,vol,vg}} (\\Delta Z)^{-1}$","1.0");
	 outputReducer->addMetadata(outputReducer->size()-2,"T/m","$\\mathrm{T}\\,\\mathrm{m}^{-1}$","$\\Delta B_{Z,\\mathrm{per,vol,vg}} (\\Delta X)^{-1}$","1.0");
	 outputReducer->addMetadata(outputReducer->size()-1,"T/m","$\\mathrm{T}\\,\\mathrm{m}^{-1}$","$\\Delta B_{Z,\\mathrm{per,vol,vg}} (\\Delta Y)^{-1}$","1.0");
         continue;
      }
      if(lowercase == "vg_gridcoordinates") {
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
         continue;
      }
      if(lowercase == "fg_gridcoordinates") { 
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_x",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2>& perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 2>& EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, 2>& EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, 2>& EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 2>& momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, 2>& dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 2>& dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, 2>& BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, 2>& volGrid,
                      FsGrid< fsgrids::technical, 2>& technicalGrid)->std::vector<double> {

               std::array<int32_t,3>& gridSize = technicalGrid.getLocalSize();
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
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2>& perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 2>& EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, 2>& EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, 2>& EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 2>& momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, 2>& dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 2>& dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, 2>& BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, 2>& volGrid,
                      FsGrid< fsgrids::technical, 2>& technicalGrid)->std::vector<double> {

               std::array<int32_t,3>& gridSize = technicalGrid.getLocalSize();
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
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2>& perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 2>& EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, 2>& EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, 2>& EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 2>& momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, 2>& dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 2>& dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, 2>& BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, 2>& volGrid,
                      FsGrid< fsgrids::technical, 2>& technicalGrid)->std::vector<double> {

               std::array<int32_t,3>& gridSize = technicalGrid.getLocalSize();
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
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2>& perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 2>& EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, 2>& EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, 2>& EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 2>& momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, 2>& dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 2>& dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, 2>& BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, 2>& volGrid,
                      FsGrid< fsgrids::technical, 2>& technicalGrid)->std::vector<double> {

               std::array<int32_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2], technicalGrid.DX);
               return retval;
         }
         ));
	 outputReducer->addMetadata(outputReducer->size()-1,"m","$\\mathrm{m}$","$\\delta X_\\mathrm{fg}$","1.0");
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_dy",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2>& perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 2>& EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, 2>& EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, 2>& EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 2>& momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, 2>& dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 2>& dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, 2>& BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, 2>& volGrid,
                      FsGrid< fsgrids::technical, 2>& technicalGrid)->std::vector<double> {

               std::array<int32_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2], technicalGrid.DY);
               return retval;
         }
         ));
	 outputReducer->addMetadata(outputReducer->size()-1,"m","$\\mathrm{m}$","$\\delta Y_\\mathrm{fg}$","1.0");
         outputReducer->addOperator(new DRO::DataReductionOperatorFsGrid("fg_dz",[](
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2>& perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 2>& EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, 2>& EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, 2>& EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 2>& momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, 2>& dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 2>& dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, 2>& BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, 2>& volGrid,
                      FsGrid< fsgrids::technical, 2>& technicalGrid)->std::vector<double> {

               std::array<int32_t,3>& gridSize = technicalGrid.getLocalSize();
               std::vector<double> retval(gridSize[0]*gridSize[1]*gridSize[2], technicalGrid.DZ);
               return retval;
         }
         ));
	 outputReducer->addMetadata(outputReducer->size()-1,"m","$\\mathrm{m}$","$\\delta Z_\\mathrm{fg}$","1.0");
         continue;
      }
      if(lowercase == "meshdata") {
         outputReducer->addOperator(new DRO::VariableMeshData);
	 outputReducer->addMetadata(outputReducer->size()-1,"","","\\mathrm{Mesh data}$","");
         continue;
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
      const std::string lowercase = boost::algorithm::to_lower_copy(*it);

      if(lowercase == "populations_blocks" || lowercase == "populations_vg_blocks") {
         // Per-population total block counts
         for(unsigned int i =0; i < getObjectWrapper().particleSpecies.size(); i++) {
            diagnosticReducer->addOperator(new DRO::Blocks(i));
         }
         continue;
      }
      if(lowercase == "vg_rhom" || lowercase == "rhom") {
         // Overall mass density
         diagnosticReducer->addOperator(new DRO::DataReductionOperatorCellParams("rhom",CellParams::RHOM,1));
         continue;
      }
      if(lowercase == "populations_rholossadjust" || lowercase == "populations_rho_loss_adjust" || lowercase == "populations_vg_rho_loss_adjust") {
         // Per-particle overall lost particle number
         for(unsigned int i =0; i < getObjectWrapper().particleSpecies.size(); i++) {
            species::Species& species=getObjectWrapper().particleSpecies[i];
            const std::string& pop = species.name;
            diagnosticReducer->addOperator(new DRO::DataReductionOperatorPopulations<Real>(pop + "/rho_loss_adjust", i, offsetof(spatial_cell::Population, RHOLOSSADJUST), 1));
         }
         continue;
      }
      //if(lowercase == "rholossvelboundary") {
      //   diagnosticReducer->addOperator(new DRO::DataReductionOperatorCellParams("rho_loss_velocity_boundary",CellParams::RHOLOSSVELBOUNDARY,1));
      //   continue;
      //}
      if(lowercase == "lbweight" || lowercase == "vg_lbweight" || lowercase == "vg_loadbalanceweight" || lowercase == "vg_loadbalance_weight" || lowercase == "loadbalance_weight") {
         diagnosticReducer->addOperator(new DRO::DataReductionOperatorCellParams("loadbalance_weight",CellParams::LBWEIGHTCOUNTER,1));
         continue;
      }
      if(lowercase == "maxvdt" || lowercase == "maxdt_acceleration" || lowercase == "vg_maxdt_acceleration") {
         diagnosticReducer->addOperator(new DRO::DataReductionOperatorCellParams("maxdt_acceleration",CellParams::MAXVDT,1));
         continue;
      }
      if(lowercase == "maxrdt" || lowercase == "maxdt_translation" || lowercase == "vg_maxdt_translation") {
         diagnosticReducer->addOperator(new DRO::DataReductionOperatorCellParams("maxdt_translation",CellParams::MAXRDT,1));
         continue;
      }
      if(lowercase == "maxfieldsdt" || lowercase == "maxdt_fieldsolver" || lowercase == "fg_maxfieldsdt" || lowercase == "fg_maxdt_fieldsolver") {
         diagnosticReducer->addOperator(new DRO::DataReductionOperatorCellParams("maxdt_fieldsolver",CellParams::MAXFDT,1));
         continue;
      }
      if(lowercase == "populations_maxdistributionfunction") {
         for(unsigned int i =0; i < getObjectWrapper().particleSpecies.size(); i++) {
            diagnosticReducer->addOperator(new DRO::MaxDistributionFunction(i));
         }
         continue;
      }
      if(lowercase == "populations_mindistributionfunction") {
         for(unsigned int i =0; i < getObjectWrapper().particleSpecies.size(); i++) {
            diagnosticReducer->addOperator(new DRO::MinDistributionFunction(i));
         }
         continue;
      }
      if(lowercase == "populations_maxrdt" || lowercase == "populations_maxdt_translation" || lowercase == "populations_vg_maxdt_translation") {
         for(unsigned int i =0; i < getObjectWrapper().particleSpecies.size(); i++) {
            species::Species& species=getObjectWrapper().particleSpecies[i];
            const std::string& pop = species.name;
            diagnosticReducer->addOperator(new DRO::DataReductionOperatorPopulations<Real>(pop + "/maxdt_translation", i, offsetof(spatial_cell::Population, max_dt[0]), 1));
         }
         continue;
      }
      if(lowercase == "populations_maxvdt" || lowercase == "populations_maxdt_acceleration" || lowercase == "populations_vg_maxdt_acceleration") {
         for(unsigned int i =0; i < getObjectWrapper().particleSpecies.size(); i++) {
            species::Species& species=getObjectWrapper().particleSpecies[i];
            const std::string& pop = species.name;
            diagnosticReducer->addOperator(new DRO::DataReductionOperatorPopulations<Real>(pop + "/maxdt_acceleration", i, offsetof(spatial_cell::Population, max_dt[1]), 1));
         }
         continue;
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

/** Ask a DataReductionOperator if it wants to take care of writing the data 
 * to output file instead of letting be handled in iowrite.cpp. 
 * @param operatorID ID number of the DataReductionOperator.
 * @return If true, then VLSVWriter should be passed to the DataReductionOperator.*/
bool DataReducer::handlesWriting(const unsigned int& operatorID) const {
   if (operatorID >= operators.size()) return false;
   return dynamic_cast<DRO::DataReductionOperatorHandlesWriting*>(operators[operatorID]) != nullptr;
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

/** Write all data from given DataReductionOperator to the output file.
 * @param operatorID ID number of the selected DataReductionOperator.
 * @param mpiGrid Parallel grid library.
 * @param cells Vector containing spatial cell IDs.
 * @param meshName Name of the spatial mesh in the output file.
 * @param vlsvWriter VLSV file writer that has output file open.
 * @return If true, DataReductionOperator wrote its data successfully.*/
bool DataReducer::writeData(const unsigned int& operatorID,
                  const dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                  const std::vector<CellID>& cells,const std::string& meshName,
                  vlsv::Writer& vlsvWriter) {
   if (operatorID >= operators.size()) return false;
   DRO::DataReductionOperatorHandlesWriting* writingOperator = dynamic_cast<DRO::DataReductionOperatorHandlesWriting*>(operators[operatorID]);
   if(writingOperator == nullptr) {
      return false;
   }
   return writingOperator->writeData(mpiGrid,cells,meshName,vlsvWriter);
}

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
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2>& perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, 2>& EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, 2>& EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, 2>& EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 2>& momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, 2>& dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, 2>& dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, 2>& BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, 2>& volGrid,
                      FsGrid< fsgrids::technical, 2>& technicalGrid,
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
