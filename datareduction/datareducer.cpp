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
using namespace std;

void initializeDataReducers(DataReducer * outputReducer, DataReducer * diagnosticReducer)
{
   typedef Parameters P;

   vector<string>::const_iterator it;
   for (it = P::outputVariableList.begin();
        it != P::outputVariableList.end();
        it++) {
      if(*it == "B") { // Bulk magnetic field at Yee-Lattice locations
         outputReducer->addOperator(new DRO::VariableB);
         continue;
      }
      if(*it == "BackgroundB") { // Static (typically dipole) magnetic field part
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("background_B",CellParams::BGBX,3));
         continue;
      }
      if(*it == "PerturbedB") { // Fluctuating magnetic field part
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("perturbed_B",CellParams::PERBX,3));
         continue;
      }
      if(*it == "E") { // Bulk electric field at Yee-lattice locations
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("E",CellParams::EX,3));
         continue;
      }
      if(*it == "Rhom") { // Overall mass density (summed over all populations)
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("rhom",CellParams::RHOM,1));
         continue;
      }
      if(*it == "Rhoq") { // Overall charge density (summed over all populations)
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("rhoq",CellParams::RHOQ,1));
         continue;
      }
      if(*it == "populations_Rho") { // Per-population particle number density
         for(unsigned int i =0; i < getObjectWrapper().particleSpecies.size(); i++) {
            species::Species& species=getObjectWrapper().particleSpecies[i];
            const std::string& pop = species.name;
            outputReducer->addOperator(new DRO::DataReductionOperatorPopulations<Real>(pop + "/rho", i, offsetof(spatial_cell::Population, RHO), 1));
         }
         continue;
      }
      
      if(*it == "V") { // Overall effective bulk density defining the center-of-mass frame from all populations
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("V",CellParams::VX,3));
         continue;
      }
      if(*it == "populations_V") { // Per population bulk velocities
         for(unsigned int i =0; i < getObjectWrapper().particleSpecies.size(); i++) {
            species::Species& species=getObjectWrapper().particleSpecies[i];
            const std::string& pop = species.name;
            outputReducer->addOperator(new DRO::DataReductionOperatorPopulations<Real>(pop + "/V", i, offsetof(spatial_cell::Population, V), 3));
         }
         continue;
      }
      if(*it == "populations_moments_Backstream") { // Per-population moments of the backstreaming part
         for(unsigned int i =0; i < getObjectWrapper().particleSpecies.size(); i++) {
            outputReducer->addOperator(new DRO::VariableRhoBackstream(i));
            outputReducer->addOperator(new DRO::VariableVBackstream(i));
            outputReducer->addOperator(new DRO::VariablePTensorBackstreamDiagonal(i));
            outputReducer->addOperator(new DRO::VariablePTensorBackstreamOffDiagonal(i));
         }
         continue;
      }
      if(*it == "populations_moments_NonBackstream") { // Per-population moments of the non-backstreaming (thermal?) part.
         for(unsigned int i =0; i < getObjectWrapper().particleSpecies.size(); i++) {
            outputReducer->addOperator(new DRO::VariableRhoNonBackstream(i));
            outputReducer->addOperator(new DRO::VariableVNonBackstream(i));
            outputReducer->addOperator(new DRO::VariablePTensorNonBackstreamDiagonal(i));
            outputReducer->addOperator(new DRO::VariablePTensorNonBackstreamOffDiagonal(i));
         }
         continue;
      }
      if(*it == "populations_MinValue" || *it == "populations_EffectiveSparsityThreshold") {
         // Effective sparsity threshold affecting each cell, if dynamic threshould algorithm is used
         for(unsigned int i =0; i < getObjectWrapper().particleSpecies.size(); i++) {
            outputReducer->addOperator(new DRO::VariableEffectiveSparsityThreshold(i));
         }
         continue;
      }
      if(*it == "populations_RhoLossAdjust") {
         // Accumulated lost particle number, per population, in each cell, since last restart
         for(unsigned int i =0; i < getObjectWrapper().particleSpecies.size(); i++) {
            species::Species& species=getObjectWrapper().particleSpecies[i];
            const std::string& pop = species.name;
            outputReducer->addOperator(new DRO::DataReductionOperatorPopulations<Real>(pop + "/rho_loss_adjust", i, offsetof(spatial_cell::Population, RHOLOSSADJUST), 1));
         }
         continue;
      }
      if(*it == "LBweight") {
         // Load balance metric for LB debugging
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("LB_weight",CellParams::LBWEIGHTCOUNTER,1));
         continue;
      }
      if(*it == "MaxVdt") {
         // Overall maximum timestep constraint as calculated by the velocity space vlasov update
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("max_v_dt",CellParams::MAXVDT,1));
         continue;
      }
      if(*it == "populations_MaxVdt") {
         // Per-population maximum timestep constraint as calculated by the velocity space vlasov update
         for(unsigned int i =0; i < getObjectWrapper().particleSpecies.size(); i++) {
            species::Species& species=getObjectWrapper().particleSpecies[i];
            const std::string& pop = species.name;
            outputReducer->addOperator(new DRO::DataReductionOperatorPopulations<Real>(pop + "/MaxVdt", i, offsetof(spatial_cell::Population, max_dt[1]), 1));
         }
         continue;
      }
      if(*it == "MaxRdt") {
         // Overall maximum timestep constraint as calculated by the real space vlasov update
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("max_r_dt",CellParams::MAXRDT,1));
         continue;
      }
      if(*it == "populations_MaxRdt") {
         // Per-population maximum timestep constraint as calculated by the real space vlasov update
         for(unsigned int i =0; i < getObjectWrapper().particleSpecies.size(); i++) {
            species::Species& species=getObjectWrapper().particleSpecies[i];
            const std::string& pop = species.name;
            outputReducer->addOperator(new DRO::DataReductionOperatorPopulations<Real>(pop + "/MaxRdt", i, offsetof(spatial_cell::Population, max_dt[0]), 1));
         }
         continue;
      }
      if(*it == "MaxFieldsdt") {
         // Maximum timestep constraint as calculated by the fieldsolver
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("max_fields_dt",CellParams::MAXFDT,1));
         continue;
      }
      if(*it == "MPIrank") {
         // Map of spatial decomposition of the DCCRG grid into MPI ranks
         outputReducer->addOperator(new DRO::MPIrank);
         continue;
      }
      if(*it == "FsGridRank") {
         // Map of spatial decomposition of the FsGrid into MPI ranks
         outputReducer->addOperator(new DRO::FsGridRank);
         continue;
      }
      if(*it == "BoundaryType") {
         // Type of boundarycells
         outputReducer->addOperator(new DRO::BoundaryType);
         continue;
      }
      if(*it == "FsGridBoundaryType") {
         // Type of boundarycells as stored in FSGrid
         outputReducer->addOperator(new DRO::FsGridBoundaryType);
         continue;
      }
      if(*it == "BoundaryLayer") {
         // For boundaries with multiple layers: layer count per cell
         outputReducer->addOperator(new DRO::BoundaryLayer);
         continue;
      }
      if (*it == "populations_Blocks") {
         // Per-population velocity space block counts
         for(unsigned int i =0; i < getObjectWrapper().particleSpecies.size(); i++) {
            outputReducer->addOperator(new DRO::Blocks(i));
         }
         continue;
      }
      if(*it == "fSaved") {
         // Boolean marker whether a velocity space is saved in a given spatial cell
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("fSaved",CellParams::ISCELLSAVINGF,1));
         continue;
      }
      if(*it == "populations_accSubcycles") {
         // Per-population number of subcycles performed for velocity space update
         for(unsigned int i =0; i < getObjectWrapper().particleSpecies.size(); i++) {
            species::Species& species=getObjectWrapper().particleSpecies[i];
            const std::string& pop = species.name;
            outputReducer->addOperator(new DRO::DataReductionOperatorPopulations<uint>(pop + "/acc_subcycles", i, offsetof(spatial_cell::Population, ACCSUBCYCLES), 1));
         }
         continue;
      }
      if(*it == "VolE") {
         // Volume-averaged E field
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("E_vol",CellParams::EXVOL,3));
         continue;
      }
      if(*it == "HallE") {
         // 12 corner components of the hall-effect contribution to the electric field
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("EXHALL_000_100",CellParams::EXHALL_000_100,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("EXHALL_001_101",CellParams::EXHALL_001_101,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("EXHALL_010_110",CellParams::EXHALL_010_110,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("EXHALL_011_111",CellParams::EXHALL_011_111,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("EYHALL_000_010",CellParams::EYHALL_000_010,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("EYHALL_001_011",CellParams::EYHALL_001_011,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("EYHALL_100_110",CellParams::EYHALL_100_110,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("EYHALL_101_111",CellParams::EYHALL_101_111,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("EZHALL_000_001",CellParams::EZHALL_000_001,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("EZHALL_010_011",CellParams::EZHALL_010_011,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("EZHALL_100_101",CellParams::EZHALL_100_101,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("EZHALL_110_111",CellParams::EZHALL_110_111,1));
         continue;
      }
      if(*it =="GradPeE") {
         // Electron pressure gradient contribution to the generalized ohm's law
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("EGRADPE",CellParams::EXGRADPE,3));
         continue;
      }
      if(*it == "VolB") {
         // Volume-averaged magnetic field
         outputReducer->addOperator(new DRO::VariableBVol);
         continue;
      }
      if(*it == "BackgroundVolB") {
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("BGB_vol",CellParams::BGBXVOL,3));
         continue;
      }
      if(*it == "PerturbedVolB") {
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("PERB_vol",CellParams::PERBXVOL,3));
         continue;
      }
      if(*it == "Pressure") {
         // Overall scalar pressure from all populations
         outputReducer->addOperator(new DRO::VariablePressureSolver);
         continue;
      }
      if(*it == "populations_PTensor") {
         // Per-population pressure tensor, stored as diagonal and offdiagonal components
         for(unsigned int i =0; i < getObjectWrapper().particleSpecies.size(); i++) {
            outputReducer->addOperator(new DRO::VariablePTensorDiagonal(i));
            outputReducer->addOperator(new DRO::VariablePTensorOffDiagonal(i));
         }
         continue;
      }
      if(*it == "derivs") {
         // Derivatives of all quantities that might be of interest
         outputReducer->addOperator(new DRO::DataReductionOperatorDerivatives("drhomdx",fieldsolver::drhomdx,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorDerivatives("drhomdy",fieldsolver::drhomdy,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorDerivatives("drhomdz",fieldsolver::drhomdz,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorDerivatives("drhoqdx",fieldsolver::drhoqdx,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorDerivatives("drhoqdy",fieldsolver::drhoqdy,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorDerivatives("drhoqdz",fieldsolver::drhoqdz,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorDerivatives("dp11dx",fieldsolver::dp11dx,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorDerivatives("dp22dx",fieldsolver::dp22dx,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorDerivatives("dp33dx",fieldsolver::dp33dx,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorDerivatives("dp11dy",fieldsolver::dp11dy,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorDerivatives("dp22dy",fieldsolver::dp22dy,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorDerivatives("dp33dy",fieldsolver::dp33dy,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorDerivatives("dp11dz",fieldsolver::dp11dz,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorDerivatives("dp22dz",fieldsolver::dp22dz,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorDerivatives("dp33dz",fieldsolver::dp33dz,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorDerivatives("dPERBxdy",fieldsolver::dPERBxdy,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorDerivatives("dBGBxdy",fieldsolver::dBGBxdy,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorDerivatives("dPERBxdz",fieldsolver::dPERBxdz,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorDerivatives("dBGBxdz",fieldsolver::dBGBxdz,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorDerivatives("dPERBydx",fieldsolver::dPERBydx,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorDerivatives("dBGBydx",fieldsolver::dBGBydx,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorDerivatives("dPERBydz",fieldsolver::dPERBydz,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorDerivatives("dBGBydz",fieldsolver::dBGBydz,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorDerivatives("dPERBzdx",fieldsolver::dPERBzdx,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorDerivatives("dBGBzdx",fieldsolver::dBGBzdx,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorDerivatives("dPERBzdy",fieldsolver::dPERBzdy,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorDerivatives("dBGBzdy",fieldsolver::dBGBzdy,1));
         if(Parameters::ohmHallTerm == 2) {
            outputReducer->addOperator(new DRO::DataReductionOperatorDerivatives("dPERBxdyy",fieldsolver::dPERBxdyy,1));
            outputReducer->addOperator(new DRO::DataReductionOperatorDerivatives("dPERBxdzz",fieldsolver::dPERBxdzz,1));
            outputReducer->addOperator(new DRO::DataReductionOperatorDerivatives("dPERBydxx",fieldsolver::dPERBydxx,1));
            outputReducer->addOperator(new DRO::DataReductionOperatorDerivatives("dPERBydzz",fieldsolver::dPERBydzz,1));
            outputReducer->addOperator(new DRO::DataReductionOperatorDerivatives("dPERBzdxx",fieldsolver::dPERBzdxx,1));
            outputReducer->addOperator(new DRO::DataReductionOperatorDerivatives("dPERBzdyy",fieldsolver::dPERBzdyy,1));
            outputReducer->addOperator(new DRO::DataReductionOperatorDerivatives("dPERBxdyz",fieldsolver::dPERBxdyz,1));
            outputReducer->addOperator(new DRO::DataReductionOperatorDerivatives("dPERBydxz",fieldsolver::dPERBydxz,1));
            outputReducer->addOperator(new DRO::DataReductionOperatorDerivatives("dPERBzdxy",fieldsolver::dPERBzdxy,1));
         }
         outputReducer->addOperator(new DRO::DataReductionOperatorDerivatives("dVxdx",fieldsolver::dVxdx,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorDerivatives("dVxdy",fieldsolver::dVxdy,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorDerivatives("dVxdz",fieldsolver::dVxdz,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorDerivatives("dVydx",fieldsolver::dVydx,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorDerivatives("dVydy",fieldsolver::dVydy,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorDerivatives("dVydz",fieldsolver::dVydz,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorDerivatives("dVzdx",fieldsolver::dVzdx,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorDerivatives("dVzdy",fieldsolver::dVzdy,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorDerivatives("dVzdz",fieldsolver::dVzdz,1));
         continue;
      }
      if(*it == "BVOLderivs") {
         // Volume-averaged derivatives
         outputReducer->addOperator(new DRO::DataReductionOperatorBVOLDerivatives("dPERBXVOLdy",bvolderivatives::dPERBXVOLdy,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorBVOLDerivatives("dBGBXVOLdy",bvolderivatives::dBGBXVOLdy,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorBVOLDerivatives("dPERBXVOLdz",bvolderivatives::dPERBXVOLdz,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorBVOLDerivatives("dBGBXVOLdz",bvolderivatives::dBGBXVOLdz,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorBVOLDerivatives("dPERBYVOLdx",bvolderivatives::dPERBYVOLdx,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorBVOLDerivatives("dBGBYVOLdx",bvolderivatives::dBGBYVOLdx,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorBVOLDerivatives("dPERBYVOLdz",bvolderivatives::dPERBYVOLdz,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorBVOLDerivatives("dBGBYVOLdz",bvolderivatives::dBGBYVOLdz,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorBVOLDerivatives("dPERBZVOLdx",bvolderivatives::dPERBZVOLdx,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorBVOLDerivatives("dBGBZVOLdx",bvolderivatives::dBGBZVOLdx,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorBVOLDerivatives("dPERBZVOLdy",bvolderivatives::dPERBZVOLdy,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorBVOLDerivatives("dBGBZVOLdy",bvolderivatives::dBGBZVOLdy,1));
         continue;
      }
      if(*it == "GridCoordinates") {
         // Spatial coordinates for each cell
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("X",CellParams::XCRD,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("Y",CellParams::YCRD,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("Z",CellParams::ZCRD,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("DX",CellParams::DX,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("DY",CellParams::DY,1));
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("DZ",CellParams::DZ,1));
         continue;
      }
      
      if (*it == "Potential") {
         // Poisson soler potential
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("poisson/potential",CellParams::PHI,1));
         continue;
      }
      if (*it == "BackgroundVolE") {
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("poisson/BGE_vol",CellParams::BGEXVOL,3));
         continue;
      }
      if (*it == "ChargeDensity") {
         // Poisson-solver charge density
         // TODO: This is redundant with Rhoq
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("poisson/rho_q",CellParams::RHOQ_TOT,1));
         continue;
      }
      if (*it == "PotentialError") {
         // Poisson solver convergence measure
         outputReducer->addOperator(new DRO::DataReductionOperatorCellParams("poisson/pot_error",CellParams::PHI_TMP,1));
         continue;
      }
      if (*it == "MeshData") {
         outputReducer->addOperator(new DRO::VariableMeshData);
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
      if(*it == "FluxB") {
         // Overall magnetic flux through the simulation plane
         diagnosticReducer->addOperator(new DRO::DiagnosticFluxB);
         continue;
      }
      if(*it == "FluxE") {
         // Overall electric flux through the simulation plane
         diagnosticReducer->addOperator(new DRO::DiagnosticFluxE);
         continue;
      }
      if (*it == "populations_Blocks") {
         // Per-population total block counts
         for(unsigned int i =0; i < getObjectWrapper().particleSpecies.size(); i++) {
            diagnosticReducer->addOperator(new DRO::Blocks(i));
         }
         continue;
      }
      if(*it == "Rhom") {
         // Overall mass density
         diagnosticReducer->addOperator(new DRO::DataReductionOperatorCellParams("rho",CellParams::RHOM,1));
         continue;
      }
      if(*it == "populations_RhoLossAdjust") {
         // Per-particle overall lost particle number
         for(unsigned int i =0; i < getObjectWrapper().particleSpecies.size(); i++) {
            species::Species& species=getObjectWrapper().particleSpecies[i];
            const std::string& pop = species.name;
            diagnosticReducer->addOperator(new DRO::DataReductionOperatorPopulations<Real>(pop + "/rho_loss_adjust", i, offsetof(spatial_cell::Population, RHOLOSSADJUST), 1));
         }
         continue;
      }
      //if(*it == "RhoLossVelBoundary") {
      //   diagnosticReducer->addOperator(new DRO::DataReductionOperatorCellParams("rho_loss_velocity_boundary",CellParams::RHOLOSSVELBOUNDARY,1));
      //   continue;
      //}
      if(*it == "LBweight") {
         diagnosticReducer->addOperator(new DRO::DataReductionOperatorCellParams("LB_weight",CellParams::LBWEIGHTCOUNTER,1));
         continue;
      }
      if(*it == "MaxVdt") {
         diagnosticReducer->addOperator(new DRO::DataReductionOperatorCellParams("max_v_dt",CellParams::MAXVDT,1));
         continue;
      }
      if(*it == "MaxRdt") {
         diagnosticReducer->addOperator(new DRO::DataReductionOperatorCellParams("max_r_dt",CellParams::MAXRDT,1));
         continue;
      }
      if(*it == "MaxFieldsdt") {
         diagnosticReducer->addOperator(new DRO::DataReductionOperatorCellParams("max_fields_dt",CellParams::MAXFDT,1));
         continue;
      }
      if(*it == "populations_MaxDistributionFunction") {
         for(unsigned int i =0; i < getObjectWrapper().particleSpecies.size(); i++) {
            diagnosticReducer->addOperator(new DRO::MaxDistributionFunction(i));
         }
         continue;
      }
      if(*it == "populations_MinDistributionFunction") {
         for(unsigned int i =0; i < getObjectWrapper().particleSpecies.size(); i++) {
            diagnosticReducer->addOperator(new DRO::MinDistributionFunction(i));
         }
         continue;
      }
      if(*it == "populations_MaxRdt") {
         for(unsigned int i =0; i < getObjectWrapper().particleSpecies.size(); i++) {
            species::Species& species=getObjectWrapper().particleSpecies[i];
            const std::string& pop = species.name;
            diagnosticReducer->addOperator(new DRO::DataReductionOperatorPopulations<Real>(pop + "/Blocks", i, offsetof(spatial_cell::Population, max_dt[0]), 1));
         }
         continue;
      }
      if(*it == "populations_MaxVdt") {
         for(unsigned int i =0; i < getObjectWrapper().particleSpecies.size(); i++) {
            species::Species& species=getObjectWrapper().particleSpecies[i];
            const std::string& pop = species.name;
            diagnosticReducer->addOperator(new DRO::DataReductionOperatorPopulations<Real>(pop + "/Blocks", i, offsetof(spatial_cell::Population, max_dt[1]), 1));
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

/** Ask a DataReductionOperator if it wants to take care of writing the data 
 * to output file instead of letting be handled in iowrite.cpp. 
 * @param operatorID ID number of the DataReductionOperator.
 * @return If true, then VLSVWriter should be passed to the DataReductionOperator.*/
bool DataReducer::handlesWriting(const unsigned int& operatorID) const {
   if (operatorID >= operators.size()) return false;
   return dynamic_cast<DRO::DataReductionOperatorHandlesWriting*>(operators[operatorID]) != nullptr;
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
