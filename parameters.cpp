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

#include "parameters.h"
#include "object_wrapper.h"
#include "particle_species.h"
#include "readparameters.h"
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <set>
#include <unistd.h>

#include "fieldtracing/fieldtracing.h"

#ifndef NAN
#define NAN 0
#endif

using namespace std;

typedef Parameters P;

extern Logger logFile;
// Using numeric_limits<Real>::max() leads to FP exceptions inside boost programoptions, use a slightly smaller value to
// avoid...

const Real LARGE_REAL = 1e20;
// Define static members:
Real P::xmin = NAN;
Real P::xmax = NAN;
Real P::ymin = NAN;
Real P::ymax = NAN;
Real P::zmin = NAN;
Real P::zmax = NAN;
Real P::dx_ini = NAN;
Real P::dy_ini = NAN;
Real P::dz_ini = NAN;

uint P::xcells_ini = numeric_limits<uint>::max();
uint P::ycells_ini = numeric_limits<uint>::max();
uint P::zcells_ini = numeric_limits<uint>::max();

Real P::t = 0;
Real P::t_min = 0;
Real P::t_max = LARGE_REAL;
Real P::dt_ceil = -1.0; 
Real P::dt = NAN;
Real P::vlasovSolverMaxCFL = NAN;
Real P::vlasovSolverMinCFL = NAN;
bool P::vlasovSolverGhostTranslate = false;
uint P::vlasovSolverGhostTranslateExtent = 0;
Real P::fieldSolverMaxCFL = NAN;
Real P::fieldSolverMinCFL = NAN;
uint P::fieldSolverSubcycles = 1;


uint P::tstep = 0;
uint P::tstep_min = 0;
uint P::tstep_max = 0;
uint P::diagnosticInterval = numeric_limits<uint>::max();
bool P::writeInitialState = true;
bool P::writeFullBGB = false;

bool P::meshRepartitioned = true;
bool P::prepareForRebalance = false;
vector<CellID> P::localCells;

bool P::adaptGPUWID = true;
uint P::GPUallocations = 128;

vector<string> P::systemWriteName;
vector<string> P::systemWritePath;
vector<Real> P::systemWriteTimeInterval;
vector<int> P::systemWriteDistributionWriteStride;
vector<int> P::systemWriteDistributionWriteXlineStride;
vector<int> P::systemWriteDistributionWriteYlineStride;
vector<int> P::systemWriteDistributionWriteZlineStride;
vector<Real> P::systemWriteDistributionWriteShellRadius;
vector<int> P::systemWriteDistributionWriteShellStride;
vector<bool> P::systemWriteFsGrid;
bool P::systemWriteAllDROs;
bool P::diagnosticWriteAllDROs;
vector<int> P::systemWrites;
vector<pair<string, string>> P::systemWriteHints;
vector<pair<string, string>> P::restartWriteHints;
vector<pair<string, string>> P::restartReadHints;

Real P::saveRestartWalltimeInterval = -1.0;
uint P::saveRecoverTstepInterval = 0;
uint P::exitAfterRestarts = numeric_limits<uint>::max();
uint P::recoverMaxFiles = 0;
uint64_t P::vlsvBufferSize = 0;
int P::restartStripeFactor = 0;
int P::systemStripeFactor = 0;
string P::restartWritePath = string("");
string P::recoverWritePath = string("");

bool P::recalculateStencils = true;
bool P::propagateVlasovAcceleration = true;
bool P::propagateVlasovTranslation = true;
bool P::propagateField = true;

bool P::dynamicTimestep = true;

Real P::maxWaveVelocity = 0.0;
uint P::maxFieldSolverSubcycles = 0.0;
int P::maxSlAccelerationSubcycles = 0.0;
Real P::resistivity = NAN;
bool P::fieldSolverDiffusiveEterms = true;
uint P::ohmHallTerm = 0;
uint P::ohmGradPeTerm = 0;
Real P::electronTemperature = 0.0;
Real P::electronDensity = 0.0;
Real P::electronPTindex = 1.0;

string P::restartFileName = string("");
bool P::isRestart = false;
int P::writeAsFloat = false;
int P::writeRestartAsFloat = false;
string P::loadBalanceAlgorithm = string("");
std::map<std::string, std::string> P::loadBalanceOptions;
uint P::rebalanceInterval = numeric_limits<uint>::max();

vector<string> P::outputVariableList;
vector<string> P::diagnosticVariableList;

string P::projectName = string("");

bool P::vlasovAccelerateMaxwellianBoundaries = false;
Real P::maxSlAccelerationRotation = 10.0;
Real P::hallMinimumRhom = physicalconstants::MASS_PROTON;
Real P::hallMinimumRhoq = physicalconstants::CHARGE;

bool P::bailout_write_restart = false;
Real P::bailout_min_dt = NAN;
Real P::bailout_max_memory = 1073741824.;
uint P::bailout_velocity_space_wall_margin = 0;

bool P::amrTransShortPencils = false;
int P::amrMaxSpatialRefLevel = 0;
int P::amrMaxAllowedSpatialRefLevel = -1;
bool P::adaptRefinement = false;
bool P::refineOnRestart = false;
bool P::forceRefinement = false;
bool P::shouldFilter = false;
bool P::useAlpha1 = true;
Real P::alpha1RefineThreshold = 0.5;
Real P::alpha1CoarsenThreshold = -1.0;
bool P::useAlpha2 = true;
Real P::alpha2RefineThreshold = 0.5;
Real P::alpha2CoarsenThreshold = -1.0;
bool P::useVorticity = false;
Real P::vorticityRefineThreshold = 0.5;
Real P::vorticityCoarsenThreshold = -1.0;
bool P::useAnisotropy = false;
Real P::anisotropyRefineThreshold = 0.5;
Real P::anisotropyCoarsenThreshold = -1.0;
int P::anisotropyMaxReflevel = 2;
Real P::alphaDRhoWeight = 1.0;
Real P::alphaDUWeight = 1.0;
Real P::alphaDPSqWeight = 1.0;
Real P::alphaDBSqWeight = 1.0;
Real P::alphaDBWeight = 1.0;

uint P::refineCadence = 5;
Real P::refineAfter = 0.0;
Real P::refineRadius = LARGE_REAL;
int P::refineBoxNumber = 0;
std::vector<Real> P::refinementMinX;
std::vector<Real> P::refinementMinY;
std::vector<Real> P::refinementMinZ;
std::vector<Real> P::refinementMaxX;
std::vector<Real> P::refinementMaxY;
std::vector<Real> P::refinementMaxZ;
int P::maxFilteringPasses = 0;
int P::amrBoxNumber = 0;
std::vector<uint> P::amrBoxHalfWidthX;
std::vector<uint> P::amrBoxHalfWidthY;
std::vector<uint> P::amrBoxHalfWidthZ;
std::vector<Realf> P::amrBoxCenterX;
std::vector<Realf> P::amrBoxCenterY;
std::vector<Realf> P::amrBoxCenterZ;
std::vector<int> P::amrBoxMaxLevel;
vector<string> P::blurPassString;
vector<int> P::numPasses;

bool P::artificialPADiff;
Realf P::PADcoefficient;
Realf P::PADCFL;
int P::PADvbins;
int P::PADmubins;
string P::PADnu0 = string("");
Realf P::PADfudge;

std::array<FsGridTools::Task_t,3> P::manualFsGridDecomposition = {0,0,0};
std::array<FsGridTools::Task_t,3> P::overrideReadFsGridDecomposition = {0,0,0};

std::string tracerString; /*!< Fieldline tracer to use for coupling ionosphere and magnetosphere */
bool P::computeCurvature;

bool P::addParameters() {
   typedef Readparameters RP;
   // the other default parameters we read through the add/get interface
   RP::add("io.diagnostic_write_interval", "Write diagnostic output every arg time steps", numeric_limits<uint>::max());

   RP::addComposing(
       "io.system_write_t_interval",
       "Save the simulation every arg simulated seconds. Negative values disable writes. [Define for all groups.]");
   RP::addComposing("io.system_write_file_name",
                    "Save the simulation to this file name series. [Define for all groups.]");
   RP::addComposing("io.system_write_path",
                    "Save this series in this location. Default is ./ [Define for all groups or none.]");
   RP::addComposing("io.system_write_distribution_stride",
                    "Every this many cells write out their velocity space. 0 is none. [Define for all groups.]");
   RP::addComposing("io.system_write_distribution_xline_stride",
                    "Every this many lines of cells along the x direction write out their velocity space. 0 is none. "
                    "[Define for all groups.]");
   RP::addComposing("io.system_write_distribution_yline_stride",
                    "Every this many lines of cells along the y direction write out their velocity space. 0 is none. "
                    "[Define for all groups.]");
   RP::addComposing("io.system_write_distribution_zline_stride",
                    "Every this many lines of cells along the z direction write out their velocity space. 0 is none. "
                    "[Define for all groups.]");
   RP::addComposing("io.system_write_distribution_shell_radius",
                    "At cells intersecting spheres with those radii centred at the origin write out their velocity "
                    "space. 0 is none.");
   RP::addComposing("io.system_write_distribution_shell_stride",
                    "Every this many cells for those on selected shells write out their velocity space. 0 is none.");
   RP::addComposing("io.system_write_fsgrid_variables", "If 0 don't write fsgrid DROs, if 1 do write them.");
   RP::addComposing(
       "io.system_write_mpiio_hint_key",
       "MPI-IO hint key passed to the non-restart IO. Has to be matched by io.system_write_mpiio_hint_value.");
   RP::addComposing(
       "io.system_write_mpiio_hint_value",
       "MPI-IO hint value passed to the non-restart IO. Has to be matched by io.system_write_mpiio_hint_key.");
   RP::addComposing(
       "io.restart_write_mpiio_hint_key",
       "MPI-IO hint key passed to the restart IO. Has to be matched by io.restart_write_mpiio_hint_value.");
   RP::addComposing(
       "io.restart_write_mpiio_hint_value",
       "MPI-IO hint value passed to the restart IO. Has to be matched by io.restart_write_mpiio_hint_key.");
   RP::addComposing(
       "io.restart_read_mpiio_hint_key",
       "MPI-IO hint key passed to the restart IO. Has to be matched by io.restart_read_mpiio_hint_value.");
   RP::addComposing(
       "io.restart_read_mpiio_hint_value",
       "MPI-IO hint value passed to the restart IO. Has to be matched by io.restart_read_mpiio_hint_key.");

   RP::add("io.write_initial_state",
           "Write initial state, not even the 0.5 dt propagation is done. Do not use for restarting. ", false);

   RP::add("io.write_full_bgb_data", "Write a dedicated file containing all BGB components and first derivatives, then exit.", false);

   RP::add("io.restart_walltime_interval",
           "Save the complete simulation in given walltime intervals. Negative values disable writes.", -1.0);
   RP::add("io.number_of_restarts", "Exit the simulation after certain number of walltime-based restarts.",
           numeric_limits<uint>::max());
   RP::add("io.recover_tstep_interval",
           "Save the complete simulation in given tstep intervals. 0 disables writes.", 0);
   RP::add("io.number_of_recovers", "Overwrite recovers cyclically after this number of recovers written.", 2);
   RP::add("io.vlsv_buffer_size",
           "Buffer size passed to VLSV writer (bytes, up to uint64_t), default 0 as this is sensible on sisu", 0);
   RP::add("io.write_restart_stripe_factor", "Stripe factor for restart and initial grid writing. Default 0 to inherit.", 0);
   RP::add("io.write_system_stripe_factor", "Stripe factor for bulk file writing. Default 0 to inherit.", 0);
   RP::add("io.write_as_float", "If true, write in floats instead of doubles", false);
   RP::add("io.restart_write_path",
           "Path to the location where restart files should be written. Defaults to the local directory, also if the "
           "specified destination is not writeable.",
           string("./"));
   RP::add("io.recover_write_path",
           "Path to the location where recover files should be written. Defaults to the local directory, also if the "
           "specified destination is not writeable.",
           string("./"));

   RP::add("propagate_field", "Propagate magnetic field during the simulation", true);
   RP::add("propagate_vlasov_acceleration",
           "Propagate distribution functions during the simulation in velocity space. If false, it is propagated with "
           "zero length timesteps.",
           true);
   RP::add("propagate_vlasov_translation",
           "Propagate distribution functions during the simulation in ordinary space. If false, it is propagated with "
           "zero length timesteps.",
           true);
   RP::add("dynamic_timestep", "If true,  timestep is set based on  CFL limits (default on)", true);
   RP::add("hallMinimumRho",
           "Minimum rho value used for the Hall and electron pressure gradient terms in the Lorentz force and in the "
           "field solver. Default is very low and has no effect in practice.",
           1.0);
   RP::add("project",
           "Specify the name of the project to use. Supported to date (20150610): Alfven Diffusion Dispersion "
           "Distributions Firehose Flowthrough Fluctuations Harris KHB Larmor Magnetosphere Multipeak Riemann1 Shock "
           "Shocktest Template test_fp testHall test_trans verificationLarmor",
           string(""));

   RP::add("restart.write_as_float", "If true, write restart fields in floats instead of doubles", false);
   RP::add("restart.filename", "Restart from this vlsv file. No restart if empty file.", string(""));

   RP::add(
       "restart.overrideReadFsGridDecompositionX",
       "Manual FsGridDecomposition for field solver grid stored in a restart file.", 0);
   RP::add(
       "restart.overrideReadFsGridDecompositionY",
       "Manual FsGridDecomposition for field solver grid stored in a restart file.", 0);
   RP::add(
       "restart.overrideReadFsGridDecompositionZ",
       "Manual FsGridDecomposition for field solver grid stored in a restart file.", 0);

   RP::add("gridbuilder.x_min", "Minimum value of the x-coordinate.", NAN);
   RP::add("gridbuilder.x_max", "Maximum value of the x-coordinate.", NAN);
   RP::add("gridbuilder.y_min", "Minimum value of the y-coordinate.", NAN);
   RP::add("gridbuilder.y_max", "Maximum value of the y-coordinate.", NAN);
   RP::add("gridbuilder.z_min", "Minimum value of the z-coordinate.", NAN);
   RP::add("gridbuilder.z_max", "Maximum value of the z-coordinate.", NAN);
   RP::add("gridbuilder.x_length", "Number of cells in x-direction in initial grid.", 0);
   RP::add("gridbuilder.y_length", "Number of cells in y-direction in initial grid.", 0);
   RP::add("gridbuilder.z_length", "Number of cells in z-direction in initial grid.", 0);

   RP::add("gridbuilder.dt", "Initial timestep in seconds.", 0.0);

   RP::add("gridbuilder.t_max",
           "Maximum simulation time, in seconds. If timestep_max limit is hit first this time will never be reached",
           LARGE_REAL);
   RP::add("gridbuilder.timestep_max",
           "Max. value for timesteps. If t_max limit is hit first, this step will never be reached",
           numeric_limits<uint>::max());
   RP::add("gridbuilder.dt_ceil",
           "Maximum simulation dt in seconds.",
           -1.0);

   // Field solver parameters
   RP::add("fieldsolver.maxWaveVelocity",
           "Maximum wave velocity allowed in the fastest velocity determination in m/s, default unlimited", LARGE_REAL);
   RP::add("fieldsolver.maxSubcycles", "Maximum allowed field solver subcycles", 1);
   RP::add("fieldsolver.resistivity", "Resistivity for the eta*J term in Ohm's law.", 0.0);
   RP::add("fieldsolver.diffusiveEterms", "Enable diffusive terms in the computation of E", true);
   RP::add(
       "fieldsolver.ohmHallTerm",
       "Enable/choose spatial order of the Hall term in Ohm's law. 0: off, 1: 1st spatial order, 2: 2nd spatial order",
       0);
   RP::add(
       "fieldsolver.ohmGradPeTerm",
       "Enable/choose spatial order of the electron pressure gradient term in Ohm's law. 0: off, 1: 1st spatial order.",
       0);
   RP::add("fieldsolver.electronTemperature",
           "Upstream (anchor point) electron temperature to be used for the electron pressure gradient term (K).", 0.0);
   RP::add("fieldsolver.electronDensity",
           "Upstream (anchor point) electron density to be used for the electron pressure gradient term (m^-3).", 0.0);
   RP::add("fieldsolver.electronPTindex",
           "Polytropic index for the equation of state to solve the electron pressure gradient term. 0 is isobaric, 1 is isothermal, 1.667 is adiabatic "
           "electrons, ",
           0.0);
   RP::add("fieldsolver.maxCFL",
           "The maximum CFL limit for field propagation. Used to set timestep if dynamic_timestep is true.", 0.5);
   RP::add("fieldsolver.minCFL",
           "The minimum CFL limit for field propagation. Used to set timestep if dynamic_timestep is true.", 0.4);

   RP::add(
       "fieldsolver.manualFsGridDecompositionX",
       "Manual FsGridDecomposition for field solver grid.", 0);
   RP::add(
       "fieldsolver.manualFsGridDecompositionY",
       "Manual FsGridDecomposition for field solver grid.", 0);
   RP::add(
       "fieldsolver.manualFsGridDecompositionZ",
       "Manual FsGridDecomposition for field solver grid.", 0);


   // Vlasov solver parameters
   RP::add("vlasovsolver.maxSlAccelerationRotation",
           "Maximum rotation angle (degrees) allowed by the Semi-Lagrangian solver (Use >25 values with care)", 25.0);
   RP::add("vlasovsolver.maxSlAccelerationSubcycles", "Maximum number of subcycles for acceleration", 1);
   RP::add("vlasovsolver.maxCFL",
           "The maximum CFL limit for vlasov propagation in ordinary space. Used to set timestep if dynamic_timestep "
           "is true.",
           0.99);
   RP::add("vlasovsolver.minCFL",
           "The minimum CFL limit for vlasov propagation in ordinary space. Used to set timestep if dynamic_timestep "
           "is true.",
           0.8);
   RP::add("vlasovsolver.accelerateMaxwellianBoundaries",
           "Propagate maxwellian boundary cell contents in velocity space. Default false.",
           false);
   RP::add("vlasovsolver.GhostTranslate","Boolean for activating all-local ghost translation",false);
   RP::add("vlasovsolver.GhostTranslateExtent","Stencil size in all-local ghost translation (default: VLASOV_STENCIL_WIDTH+1",0);

   // Load balancing parameters
   RP::add("loadBalance.algorithm", "Load balancing algorithm to be used", string("RCB"));
   RP::add("loadBalance.tolerance", "Load imbalance tolerance", string("1.05"));
   RP::add("loadBalance.rebalanceInterval", "Load rebalance interval (steps)", 10);

   RP::addComposing("loadBalance.optionKey", "Zoltan option key. Has to be matched by loadBalance.optionValue.");
   RP::addComposing("loadBalance.optionValue", "Zoltan option value. Has to be matched by loadBalance.optionKey.");

   // Output variable parameters
   RP::add("io.system_write_all_data_reducers", "If 0 don't write all DROs, if 1 do write them.", false);
   // NOTE Do not remove the : before the list of variable names as this is parsed by tools/check_vlasiator_cfg.sh
   RP::addComposing("variables.output",
                    string() +
                        "List of data reduction operators (DROs) to add to the grid file output.  Each variable to be "
                        "added has to be on a new line output = XXX. Names are case insensitive.  " +
                        "Available (20250413): " + "fg_b fg_b_background fg_b_perturbed fg_b_background_vol fg_derivs_b_background fg_e " +
                        "vg_rhom vg_rhoq populations_vg_rho " + "fg_rhom fg_rhoq " + "vg_v fg_v populations_vg_v " +
                        "populations_vg_moments_thermal populations_vg_moments_nonthermal " +
                        "populations_vg_effectivesparsitythreshold populations_vg_rho_loss_adjust " +
                        "populations_vg_energydensity populations_vg_precipitationdifferentialflux " +
                        "populations_vg_heatflux " + "populations_vg_1dmuspace " +
                        "populations_vg_nonmaxwellianity " +
                        "vg_maxdt_acceleration vg_maxdt_translation populations_vg_maxdt_acceleration " +
                        "populations_vg_maxdt_translation " +
                        "fg_maxdt_fieldsolver " + "vg_rank fg_rank fg_amr_level vg_loadbalance_weight " +
                        "vg_boundarytype fg_boundarytype vg_boundarylayer fg_boundarylayer " +
                        "populations_vg_blocks vg_f_saved " + "populations_vg_acceleration_subcycles " +
                        "vg_e_vol fg_e_vol " +
                        "fg_e_hall vg_e_gradpe fg_b_vol vg_b_vol vg_b_background_vol vg_b_perturbed_vol " +
                        "vg_pressure fg_pressure populations_vg_ptensor " + "vg_b_vol_derivatives fg_derivs " +
                        "ig_fac ig_latitude ig_chi0 ig_cellarea ig_upmappedarea ig_sigmap ig_sigmah ig_sigmaparallel ig_rhon " +
                        "ig_electrontemp ig_solverinternals ig_upmappednodecoords ig_upmappedb ig_openclosed ig_potential "+
                        "ig_precipitation ig_deltaphi "+
                        "ig_inplanecurrent ig_b ig_e vg_drift vg_ionospherecoupling vg_connection vg_fluxrope fg_curvature "+
                        "vg_amr_drho vg_amr_du vg_amr_dpsq vg_amr_dbsq vg_amr_db vg_amr_alpha1 vg_amr_reflevel vg_amr_alpha2 "+
                        "vg_gridcoordinates fg_gridcoordinates vg_pressure_anisotropy vg_amr_vorticity");

   RP::addComposing(
       "variables_deprecated.output",
       string() + "List of deprecated names for data reduction operators (DROs). Names are case insensitive. " +
           "Available (20250413): " + "B BackgroundB fg_BackgroundB PerturbedB fg_PerturbedB " + "E " +
           "Rhom Rhoq populations_Rho " + "V populations_V " +
           "populations_moments_Backstream populations_moments_NonBackstream " +
           "populations_moments_thermal populations_moments_nonthermal " +
           "populations_minvalue populations_EffectiveSparsityThreshold populations_RhoLossAdjust "
           "populations_rho_loss_adjust populations_1dmuspace " +
           "populations_EnergyDensity populations_PrecipitationFlux populations_precipitationdifferentialflux" +
           "LBweight vg_lbweight vg_loadbalanceweight MaxVdt MaxRdt populations_MaxVdt populations_MaxRdt " +
           "populations_maxdt_acceleration populations_maxdt_translation MaxFieldsdt fg_maxfieldsdt" +
           "MPIrank FsGridRank " + "FsGridBoundaryType BoundaryType FsGridBoundaryLayer BoundaryLayer " +
           "populations_Blocks fSaved vg_fsaved" + "populations_accSubcycles populations_acceleration_subcycles" +
           "VolE vg_VolE Evol E_vol fg_VolE fg_Evol " +
           "HallE fg_HallE GradPeE e_gradpe VolB vg_VolB fg_VolB B_vol Bvol vg_Bvol fg_volB fg_Bvol " +
           "BackgroundVolB PerturbedVolB " + "Pressure vg_Pressure fg_Pressure populations_PTensor " +
           "BVOLderivs b_vol_derivs");

   RP::add("io.diagnostic_write_all_data_reducers", "Write all available diagnostic reducers", false);
   // NOTE Do not remove the : before the list of variable names as this is parsed by tools/check_vlasiator_cfg.sh
   RP::addComposing("variables.diagnostic",
                    string() +
                        "List of data reduction operators (DROs) to add to the diagnostic runtime output. Each "
                        "variable to be added has to be on a new line diagnostic = XXX. Names are case insensitive. " +
                        "Available (20250130): " + "populations_vg_blocks " +
                        "vg_rhom populations_vg_rho_loss_adjust " + "vg_loadbalance_weight " +
                        "vg_maxdt_acceleration vg_maxdt_translation " + "fg_maxdt_fieldsolver " +
                        "populations_vg_maxdt_acceleration populations_vg_maxdt_translation ");

   RP::addComposing("variables_deprecated.diagnostic",
                    string() +
                        "List of deprecated data reduction operators (DROs) to add to the diagnostic runtime output. "
                        "Names are case insensitive. " +
                        "Available (20250130): " + "rhom populations_rholossadjust populations_rho_loss_adjust " +
                        "populations_blocks lbweight loadbalance_weight " + "vg_lbweight vg_loadbalanceweight " +
                        "maxvdt maxdt_acceleration " + "maxrdt maxdt_translation " +
                        "populations_maxvdt populations_maxrdt " +
                        "populations_maxdt_acceleration populations_maxdt_translation " +
                        "maxfieldsdt maxdt_fieldsolver fg_maxfieldsdt");

   // bailout parameters
   RP::add("bailout.write_restart",
           "If 1, write a restart file on bailout. Gets reset when sending a STOP (1) or a KILL (0).", true);
   RP::add("bailout.min_dt", "Minimum time step below which bailout occurs (s).", 1e-6);
   RP::add("bailout.max_memory", "Maximum amount of memory used per node (in GiB) over which bailout occurs.",
           1073741824.);
   RP::add("bailout.velocity_space_wall_block_margin", "Distance from the velocity space limits in blocks, if the distribution function reaches that distance from the wall we bail out to avoid hitting the wall.", 1);

   // Spatial Refinement parameters
   RP::add("AMR.max_spatial_level", "Maximum absolute spatial mesh refinement level", (uint)0);
   RP::add("AMR.max_allowed_spatial_level", "Maximum currently allowed spatial mesh refinement level", -1);
   RP::add("AMR.should_refine","If false, do not refine Vlasov grid regardless of max spatial level",true);
   RP::add("AMR.adapt_refinement","If true, re-refine vlasov grid every refine_cadence balance", false);
   RP::add("AMR.refine_on_restart","If true, re-refine vlasov grid on restart. DEPRECATED, consider using the DOMR command", false);
   RP::add("AMR.force_refinement","If true, refine/unrefine the vlasov grid to match the config on restart", false);
   RP::add("AMR.should_filter","If true, filter vlasov grid with boxcar filter on restart",false);
   RP::add("AMR.use_alpha1","Use the maximum of dimensionless gradients alpha_1 as a refinement index", true);
   RP::add("AMR.alpha1_refine_threshold","Determines the minimum value of alpha_1 to refine cells", 0.5);
   RP::add("AMR.alpha1_coarsen_threshold","Determines the maximum value of alpha_1 to unrefine cells, default half of the refine threshold", -1.0);
   RP::add("AMR.use_alpha2","Use J/B_perp as a refinement index", true);
   RP::add("AMR.alpha2_refine_threshold","Determines the minimum value of alpha_2 to refine cells", 0.5);
   RP::add("AMR.alpha2_coarsen_threshold","Determines the maximum value of alpha_2 to unrefine cells, default half of the refine threshold", -1.0);
   RP::add("AMR.use_vorticity","Use vorticity as a refinement index", false);
   RP::add("AMR.vorticity_refine_threshold","Determines the minimum value of vorticity to refine cells", 0.5);
   RP::add("AMR.vorticity_coarsen_threshold","Determines the maximum value of vorticity to unrefine cells, default half of the refine threshold", -1.0);
   RP::add("AMR.use_anisotropy","Use pressure anisotropy as a refinement index", false);
   RP::add("AMR.anisotropy_refine_threshold","Determines the maximum value of pressure anisotropy to refine cells", 0.5);
   RP::add("AMR.anisotropy_coarsen_threshold","Determines the minimum value of pressure anisotropy to unrefine cells, default twice the refine threshold", -1.0);
   RP::add("AMR.anisotropy_max_reflevel","When anisotropy is below the refine threshold, defines the maximum level to refine to", 2);
   RP::add("AMR.refine_cadence","Refine every nth load balance", 5);
   RP::add("AMR.refine_after","Start refinement after this many simulation seconds", 0.0);
   RP::add("AMR.refine_radius","Maximum distance from origin to allow refinement within. Only induced refinement allowed outside this radius.", LARGE_REAL);
   RP::add("AMR.number_of_refine_boxes", "How many boxes outside which to suppress refinement, that number of box edges have to then be defined as well. If more than 1 box is defined, refinement is suppressed outside the union of the volumes of all boxes.", 0);
   RP::addComposing("AMR.refinement_min_x", "Refinement minimum X coordinate, no refinement at x < this value (m) except induced refinement.");
   RP::addComposing("AMR.refinement_min_y", "Refinement minimum Y coordinate, no refinement at y < this value (m) except induced refinement.");
   RP::addComposing("AMR.refinement_min_z", "Refinement minimum Z coordinate, no refinement at z < this value (m) except induced refinement.");
   RP::addComposing("AMR.refinement_max_x", "Refinement maximum X coordinate, no refinement at x > this value (m) except induced refinement.");
   RP::addComposing("AMR.refinement_max_y", "Refinement maximum Y coordinate, no refinement at y > this value (m) except induced refinement.");
   RP::addComposing("AMR.refinement_max_z", "Refinement maximum Z coordinate, no refinement at z > this value (m) except induced refinement.");
   RP::add("AMR.alpha1_drho_weight","Multiplier for delta rho (plasma density) in alpha calculation", 1.0);
   RP::add("AMR.alpha1_du_weight","Multiplier for delta U (total kinetic + field energy density) in alpha calculation", 1.0);
   RP::add("AMR.alpha1_dpsq_weight","Multiplier for delta p squared (kinetic energy) in alpha calculation", 1.0);
   RP::add("AMR.alpha1_dbsq_weight","Multiplier for delta B squared (field energy) in alpha calculation", 1.0);
   RP::add("AMR.alpha1_db_weight","Multiplier for delta B (magnetic field strength) in alpha calculation", 1.0);
   RP::add("AMR.number_of_boxes", "How many boxes to be refined, that number of centers and sizes have to then be defined as well.", 0);
   RP::addComposing("AMR.box_half_width_x", "Half width in x of the box that is refined");
   RP::addComposing("AMR.box_half_width_y", "Half width in y of the box that is refined");
   RP::addComposing("AMR.box_half_width_z", "Half width in z of the box that is refined");
   RP::addComposing("AMR.box_center_x", "x coordinate of the center of the box that is refined");
   RP::addComposing("AMR.box_center_y", "y coordinate of the center of the box that is refined");
   RP::addComposing("AMR.box_center_z", "z coordinate of the center of the box that is refined");
   RP::addComposing("AMR.box_max_level", "max refinement level of the box that is refined");
   RP::add("AMR.transShortPencils", "if true, use one-cell pencils", false);
   RP::addComposing("AMR.filterpasses", string("AMR filter passes for each individual refinement level"));
   RP::add("adaptGPUWID", "if true, will halve velocity block counts if GPU is in use and WID==8", true);
   RP::add("GPUallocations", "How many parallel GPU vlasov allocations to make? (default 128)", 128);

   // Diffusion parameters
   RP::add("PAD.enable","Enable Artificial pitch-angle diffusion",0);
   RP::add("PAD.coefficient","Set artificial pitch-angle diffusion coefficient (overriding .DAT file)",-1);
   RP::add("PAD.CFL","Set CFL condition",0.1);
   RP::add("PAD.vbins","number of bins for velocity",200);
   RP::add("PAD.mubins","number of bins for mu",30);
   RP::add("PAD.file","Path of txt file for nu0", string("NU0BOX.DAT"));
   RP::add("PAD.fudge","Divide diffusion coefficient nu0 (read from file) by a fudge factor (see Dubart et al 2023)",4);
   
   // Fieldtracing
   RP::add("fieldtracing.fieldLineTracer", "Field line tracing method to use for coupling ionosphere and magnetosphere (options are: Euler, BS)", std::string("Euler"));
   RP::add("fieldtracing.tracer_max_allowed_error", "Maximum allowed error for the adaptive field line tracers ", 1000);
   RP::add("fieldtracing.tracer_max_attempts", "Maximum allowed attempts for the adaptive field line tracers", 100);
   RP::add("fieldtracing.tracer_min_dx", "Minimum allowed field line tracer step length for the adaptive field line tracers (m)", 100e3);
   RP::add("fieldtracing.fullbox_and_fluxrope_max_absolute_distance_to_trace", "Maximum absolute distance in m to trace along the field line before ending. Defaults to the sum of the simulation box edge lengths LX+LY+LZ if set <= 0.", -1);
   RP::add("fieldtracing.fullbox_max_incomplete_cells", "Maximum fraction of cells left incomplete when stopping tracing loop for full box tracing. Defaults to zero to process all, will be slow at scale! Both fluxrope_max_incomplete_cells and fullbox_max_incomplete_cells will be achieved.", 0);
   RP::add("fieldtracing.fluxrope_max_incomplete_cells", "Maximum fraction of cells left incomplete when stopping loop for flux rope tracing. Defaults to zero to process all, will be slow at scale! Both fluxrope_max_incomplete_cells and fullbox_max_incomplete_cells will be achieved.", 0);
   RP::add("fieldtracing.use_reconstruction_cache", "Use the cache to store reconstruction coefficients. (0: don't, 1: use)", 0);
   RP::add("fieldtracing.fluxrope_max_curvature_radii_to_trace", "Maximum number of seedpoint curvature radii to trace forward and backward from each DCCRG cell to find flux ropes", 10);
   RP::add("fieldtracing.fluxrope_max_curvature_radii_extent", "Maximum extent in seedpoint curvature radii from the seed a field line is allowed to extend to be counted as a flux rope", 2);
   RP::add("fieldtracing.min_allowed_x", "Trace for x coordinates larger than this limit (in m).", -LARGE_REAL);
   RP::add("fieldtracing.min_allowed_y", "Trace for y coordinates larger than this limit (in m).", -LARGE_REAL);
   RP::add("fieldtracing.min_allowed_z", "Trace for z coordinates larger than this limit (in m).", -LARGE_REAL);
   RP::add("fieldtracing.max_allowed_x", "Trace for x coordinates smaller than this limit (in m).", LARGE_REAL);
   RP::add("fieldtracing.max_allowed_y", "Trace for y coordinates smaller than this limit (in m).", LARGE_REAL);
   RP::add("fieldtracing.max_allowed_z", "Trace for z coordinates smaller than this limit (in m).", LARGE_REAL);

   return true;
}

void Parameters::getParameters() {
   typedef Readparameters RP;
   // get numerical values of the parameters
   RP::get("io.diagnostic_write_interval", P::diagnosticInterval);
   RP::get("io.diagnostic_write_all_data_reducers", P::diagnosticWriteAllDROs);
   RP::get("io.system_write_t_interval", P::systemWriteTimeInterval);
   RP::get("io.system_write_file_name", P::systemWriteName);
   RP::get("io.system_write_path", P::systemWritePath);
   RP::get("io.system_write_distribution_stride", P::systemWriteDistributionWriteStride);
   RP::get("io.system_write_distribution_xline_stride", P::systemWriteDistributionWriteXlineStride);
   RP::get("io.system_write_distribution_yline_stride", P::systemWriteDistributionWriteYlineStride);
   RP::get("io.system_write_distribution_zline_stride", P::systemWriteDistributionWriteZlineStride);
   RP::get("io.system_write_distribution_shell_radius", P::systemWriteDistributionWriteShellRadius);
   RP::get("io.system_write_distribution_shell_stride", P::systemWriteDistributionWriteShellStride);
   RP::get("io.system_write_fsgrid_variables", P::systemWriteFsGrid);
   RP::get("io.system_write_all_data_reducers", P::systemWriteAllDROs);
   RP::get("io.write_initial_state", P::writeInitialState);
   RP::get("io.write_full_bgb_data", P::writeFullBGB);
   RP::get("io.restart_walltime_interval", P::saveRestartWalltimeInterval);
   RP::get("io.recover_tstep_interval", P::saveRecoverTstepInterval);
   RP::get("io.number_of_restarts", P::exitAfterRestarts);
   RP::get("io.number_of_recovers", P::recoverMaxFiles);
   RP::get("io.vlsv_buffer_size", P::vlsvBufferSize);
   RP::get("io.write_restart_stripe_factor", P::restartStripeFactor);
   RP::get("io.write_system_stripe_factor", P::systemStripeFactor);
   RP::get("io.restart_write_path", P::restartWritePath);
   RP::get("io.recover_write_path", P::recoverWritePath);
   RP::get("io.write_as_float", P::writeAsFloat);

   // Checks for validity of io and restart parameters
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
   const string prefix = string("./");
   if (access(&(P::restartWritePath[0]), W_OK) != 0) {
      if (myRank == MASTER_RANK) {
         cerr << "ERROR restart write path " << P::restartWritePath << " not writeable, defaulting to local directory."
              << endl;
      }
      P::restartWritePath = prefix;
   }
   if (access(&(P::recoverWritePath[0]), W_OK) != 0) {
      if (myRank == MASTER_RANK) {
         cerr << "ERROR recover write path " << P::recoverWritePath << " not writeable, defaulting to local directory."
              << endl;
      }
      P::recoverWritePath = prefix;
   }
   if (P::recoverMaxFiles == 0) {
      P::recoverMaxFiles = 1; // If we leave it at zero a manual DORC will divide by zero when computing the index.
   }
   size_t maxSize = 0;
   maxSize = max(maxSize, P::systemWriteTimeInterval.size());
   maxSize = max(maxSize, P::systemWriteName.size());
   maxSize = max(maxSize, P::systemWritePath.size());
   maxSize = max(maxSize, P::systemWriteDistributionWriteStride.size());
   maxSize = max(maxSize, P::systemWriteDistributionWriteXlineStride.size());
   maxSize = max(maxSize, P::systemWriteDistributionWriteYlineStride.size());
   maxSize = max(maxSize, P::systemWriteDistributionWriteZlineStride.size());
   if (P::systemWriteTimeInterval.size() != maxSize) {
      if (myRank == MASTER_RANK) {
         cerr << "ERROR io.system_write_t_interval should be defined for all file types." << endl;
         MPI_Abort(MPI_COMM_WORLD, 1);
      }
   }
   if (P::systemWriteName.size() != maxSize) {
      if (myRank == MASTER_RANK) {
         cerr << "ERROR io.system_write_file_name should be defined for all file types." << endl;
         MPI_Abort(MPI_COMM_WORLD, 1);
      }
   }
   if (P::systemWritePath.size() != maxSize && P::systemWritePath.size() != 0) {
      if (myRank == MASTER_RANK) {
         cerr << "ERROR io.system_write_path should be defined for all file types or none at all." << endl;
         MPI_Abort(MPI_COMM_WORLD, 1);
      }
   }
   if (P::systemWriteDistributionWriteStride.size() != maxSize) {
      if (myRank == MASTER_RANK) {
         cerr << "ERROR io.system_write_distribution_stride should be defined for all file types." << endl;
         MPI_Abort(MPI_COMM_WORLD, 1);
      }
   }
   if (P::systemWriteDistributionWriteXlineStride.size() != maxSize) {
      if (myRank == MASTER_RANK) {
         cerr << "ERROR io.system_write_distribution_xline_stride should be defined for all file types." << endl;
         MPI_Abort(MPI_COMM_WORLD, 1);
      }
   }
   if (P::systemWriteDistributionWriteYlineStride.size() != maxSize) {
      if (myRank == MASTER_RANK) {
         cerr << "ERROR io.system_write_distribution_yline_stride should be defined for all file types." << endl;
         MPI_Abort(MPI_COMM_WORLD, 1);
      }
   }
   if (P::systemWriteDistributionWriteZlineStride.size() != maxSize) {
      if (myRank == MASTER_RANK) {
         cerr << "ERROR io.system_write_distribution_zline_stride should be defined for all file types." << endl;
         MPI_Abort(MPI_COMM_WORLD, 1);
      }
   }
   if (P::systemWriteDistributionWriteShellStride.size() != P::systemWriteDistributionWriteShellRadius.size()) {
      if (myRank == MASTER_RANK) {
         cerr << "ERROR You should set the same number of io.system_write_distribution_shell_stride "
              << "and io.system_write_distribution_shell_radius." << endl;
         MPI_Abort(MPI_COMM_WORLD, 1);
      }
   }
   if (P::systemWriteFsGrid.size() != maxSize) {
      if (P::systemWriteFsGrid.size() == 0) {
         for (uint i = 0; i < maxSize; i++) {
            P::systemWriteFsGrid.push_back(true);
         }
      } else {
         if (myRank == MASTER_RANK) {
            cerr << "ERROR io.system_write_fsgrid_variables should be defined for all file types (or none at all)." << endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
         }
      }
   }
   if (P::systemWritePath.size() == 0) {
      for (uint i = 0; i < P::systemWriteName.size(); i++) {
         P::systemWritePath.push_back(string("./"));
      }
   } else {
      for (uint i = 0; i < P::systemWritePath.size(); i++) {
         if (access(&(P::systemWritePath.at(i)[0]), W_OK) != 0) {
            if (myRank == MASTER_RANK) {
               cerr << "ERROR " << P::systemWriteName.at(i) << " write path " << P::systemWritePath.at(i)
                    << " not writeable, defaulting to local directory." << endl;
            }
            P::systemWritePath.at(i) = prefix;
         }
      }
   }

   bool includefSaved = false;
   for(uint i=0; i<maxSize; i++) {
      if(P::systemWriteDistributionWriteStride[i] != 0 ||
         P::systemWriteDistributionWriteXlineStride[i] > 0 ||
         P::systemWriteDistributionWriteYlineStride[i] > 0 ||
         P::systemWriteDistributionWriteZlineStride[i] > 0) {
         includefSaved = true;
      }
   }
   for(uint i=0; i<P::systemWriteDistributionWriteShellRadius.size(); i++) {
      if(P::systemWriteDistributionWriteShellRadius[i] > 0) {
         includefSaved = true;
      }
   }


   vector<string> mpiioKeys, mpiioValues;
   RP::get("io.system_write_mpiio_hint_key", mpiioKeys);
   RP::get("io.system_write_mpiio_hint_value", mpiioValues);

   if (mpiioKeys.size() != mpiioValues.size()) {
      if (myRank == MASTER_RANK) {
         cerr << "WARNING the number of io.system_write_mpiio_hint_key and io.system_write_mpiio_hint_value do not "
                 "match. Disregarding these options."
              << endl;
      }
   } else {
      for (uint i = 0; i < mpiioKeys.size(); i++) {
         P::systemWriteHints.push_back({mpiioKeys[i], mpiioValues[i]});
      }
   }

   mpiioKeys.clear();
   mpiioValues.clear();
   RP::get("io.restart_write_mpiio_hint_key", mpiioKeys);
   RP::get("io.restart_write_mpiio_hint_value", mpiioValues);

   if (mpiioKeys.size() != mpiioValues.size()) {
      if (myRank == MASTER_RANK) {
         cerr << "WARNING the number of io.restart_write_mpiio_hint_key and io.restart_write_mpiio_hint_value do not "
                 "match. Disregarding these options."
              << endl;
      }
   } else {
      for (uint i = 0; i < mpiioKeys.size(); i++) {
         P::restartWriteHints.push_back({mpiioKeys[i], mpiioValues[i]});
      }
   }

   mpiioKeys.clear();
   mpiioValues.clear();
   RP::get("io.restart_read_mpiio_hint_key", mpiioKeys);
   RP::get("io.restart_read_mpiio_hint_value", mpiioValues);

   if (mpiioKeys.size() != mpiioValues.size()) {
      if (myRank == MASTER_RANK) {
         cerr << "WARNING the number of io.restart_read_mpiio_hint_key and io.restart_read_mpiio_hint_value do not "
                 "match. Disregarding these options."
              << endl;
      }
   } else {
      for (uint i = 0; i < mpiioKeys.size(); i++) {
         P::restartReadHints.push_back({mpiioKeys[i], mpiioValues[i]});
      }
   }

   RP::get("propagate_field", P::propagateField);
   RP::get("propagate_vlasov_acceleration", P::propagateVlasovAcceleration);
   RP::get("propagate_vlasov_translation", P::propagateVlasovTranslation);
   RP::get("dynamic_timestep", P::dynamicTimestep);
   Real hallRho;
   RP::get("hallMinimumRho", hallRho);
   P::hallMinimumRhom = hallRho * physicalconstants::MASS_PROTON;
   P::hallMinimumRhoq = hallRho * physicalconstants::CHARGE;
   RP::get("restart.write_as_float", P::writeRestartAsFloat);
   RP::get("restart.filename", P::restartFileName);
   P::isRestart = (P::restartFileName != string(""));

   // manual FsGrid decomposition should be complete with three values. If at least one is set but all are not set, abort
   if ((RP::isSet("restart.overrideReadFsGridDecompositionX")||RP::isSet("restart.overrideReadFsGridDecompositionY")||RP::isSet("restart.overrideReadFsGridDecompositionZ")) &&
        !(RP::isSet("restart.overrideReadFsGridDecompositionX")&&RP::isSet("restart.overrideReadFsGridDecompositionY")&&RP::isSet("restart.overrideReadFsGridDecompositionZ")) ) {
      cerr << "ERROR all of restart.overrideReadFsGridDecompositionX,Y,Z should be defined." << endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
   }
   FsGridTools::Task_t temp_task_t;
   RP::get("restart.overrideReadFsGridDecompositionX", temp_task_t);
   P::overrideReadFsGridDecomposition[0] = temp_task_t;
   RP::get("restart.overrideReadFsGridDecompositionY", temp_task_t);
   P::overrideReadFsGridDecomposition[1] = temp_task_t;
   RP::get("restart.overrideReadFsGridDecompositionZ", temp_task_t);
   P::overrideReadFsGridDecomposition[2] = temp_task_t;

   RP::get("project", P::projectName);
   if (RP::helpRequested) {
      P::projectName = string("Magnetosphere");
   }

   /*get numerical values, let Readparameters handle the conversions*/
   string geometryString;

   RP::get("gridbuilder.x_min", P::xmin);
   RP::get("gridbuilder.x_max", P::xmax);
   RP::get("gridbuilder.y_min", P::ymin);
   RP::get("gridbuilder.y_max", P::ymax);
   RP::get("gridbuilder.z_min", P::zmin);
   RP::get("gridbuilder.z_max", P::zmax);
   RP::get("gridbuilder.x_length", P::xcells_ini);
   RP::get("gridbuilder.y_length", P::ycells_ini);
   RP::get("gridbuilder.z_length", P::zcells_ini);

   RP::get("AMR.max_spatial_level", P::amrMaxSpatialRefLevel);
   RP::get("AMR.max_allowed_spatial_level", P::amrMaxAllowedSpatialRefLevel);
   if(P::amrMaxAllowedSpatialRefLevel < 0) { // negative (default is -1) just goes to max
      P::amrMaxAllowedSpatialRefLevel = P::amrMaxSpatialRefLevel; // set max allowed to the same as the absolute max
   }
   if(P::amrMaxSpatialRefLevel < P::amrMaxAllowedSpatialRefLevel) {
      if(myRank == MASTER_RANK) {
         cerr << "AMR.max_allowed_spatial_level cannot be greater than AMR.max_spatial_level!\n";
      }
      MPI_Abort(MPI_COMM_WORLD, 1);
   }
   RP::get("AMR.adapt_refinement",P::adaptRefinement);
   RP::get("AMR.refine_on_restart",P::refineOnRestart);
   RP::get("AMR.force_refinement",P::forceRefinement);
   RP::get("AMR.should_filter",P::shouldFilter);
   RP::get("AMR.use_alpha1",P::useAlpha1);
   RP::get("AMR.alpha1_refine_threshold",P::alpha1RefineThreshold);
   RP::get("AMR.alpha1_coarsen_threshold",P::alpha1CoarsenThreshold);
   if (P::useAlpha1 && P::alpha1CoarsenThreshold < 0) {
      P::alpha1CoarsenThreshold = P::alpha1RefineThreshold / 2.0;
   }
   if (P::useAlpha1 && P::alpha1RefineThreshold < 0) {
      if (myRank == MASTER_RANK) {
         cerr << "ERROR invalid alpha_1 refine threshold" << endl;
      }
      MPI_Abort(MPI_COMM_WORLD, 1);
   }

   RP::get("AMR.use_alpha2",P::useAlpha2);
   RP::get("AMR.alpha2_refine_threshold",P::alpha2RefineThreshold);
   RP::get("AMR.alpha2_coarsen_threshold",P::alpha2CoarsenThreshold);
   if (P::useAlpha2 && P::alpha2CoarsenThreshold < 0) {
      P::alpha2CoarsenThreshold = P::alpha2RefineThreshold / 2.0;
   }
   if (P::useAlpha2 && P::alpha2RefineThreshold < 0) {
      if (myRank == MASTER_RANK) {
         cerr << "ERROR invalid alpha_2 refine threshold" << endl;
      }
      MPI_Abort(MPI_COMM_WORLD, 1);
   }

   RP::get("AMR.use_vorticity",P::useVorticity);
   RP::get("AMR.vorticity_refine_threshold",P::vorticityRefineThreshold);
   RP::get("AMR.vorticity_coarsen_threshold",P::vorticityCoarsenThreshold);
   if (P::useVorticity && P::vorticityCoarsenThreshold < 0) {
      P::vorticityCoarsenThreshold = P::vorticityRefineThreshold / 2.0;
   }
   if (P::useVorticity && P::vorticityRefineThreshold < 0) {
      if (myRank == MASTER_RANK) {
         cerr << "ERROR invalid vorticity refine threshold" << endl;
      }
      MPI_Abort(MPI_COMM_WORLD, 1);
   }

   RP::get("AMR.use_anisotropy",P::useAnisotropy);
   RP::get("AMR.anisotropy_refine_threshold", P::anisotropyRefineThreshold);
   RP::get("AMR.anisotropy_coarsen_threshold", P::anisotropyCoarsenThreshold);
   RP::get("AMR.anisotropy_max_reflevel", P::anisotropyMaxReflevel);
   if (P::useAnisotropy && P::anisotropyCoarsenThreshold < 0) {
      P::anisotropyCoarsenThreshold = P::anisotropyRefineThreshold * 2.0;
   }
   if (P::useAnisotropy && P::anisotropyRefineThreshold < 0) {
      if (myRank == MASTER_RANK) {
         cerr << "ERROR invalid anisotropy refine threshold" << endl;
      }
      MPI_Abort(MPI_COMM_WORLD, 1);
   }

   RP::get("AMR.refine_cadence",P::refineCadence);
   RP::get("AMR.refine_after",P::refineAfter);
   RP::get("AMR.refine_radius",P::refineRadius);
   RP::get("AMR.number_of_refine_boxes", P::refineBoxNumber);
   RP::get("AMR.refinement_min_x", P::refinementMinX);
   RP::get("AMR.refinement_min_y", P::refinementMinY);
   RP::get("AMR.refinement_min_z", P::refinementMinZ);
   RP::get("AMR.refinement_max_x", P::refinementMaxX);
   RP::get("AMR.refinement_max_y", P::refinementMaxY);
   RP::get("AMR.refinement_max_z", P::refinementMaxZ);
   RP::get("AMR.alpha1_drho_weight", P::alphaDRhoWeight);
   RP::get("AMR.alpha1_du_weight", P::alphaDUWeight);
   RP::get("AMR.alpha1_dpsq_weight", P::alphaDPSqWeight);
   RP::get("AMR.alpha1_dbsq_weight", P::alphaDBSqWeight);
   RP::get("AMR.alpha1_db_weight", P::alphaDBWeight);
   RP::get("AMR.number_of_boxes", P::amrBoxNumber);
   RP::get("AMR.box_max_level", P::amrBoxMaxLevel);
   RP::get("AMR.box_half_width_x", P::amrBoxHalfWidthX);
   RP::get("AMR.box_half_width_y", P::amrBoxHalfWidthY);
   RP::get("AMR.box_half_width_z", P::amrBoxHalfWidthZ);
   RP::get("AMR.box_center_x", P::amrBoxCenterX);
   RP::get("AMR.box_center_y", P::amrBoxCenterY);
   RP::get("AMR.box_center_z", P::amrBoxCenterZ);
   RP::get("AMR.transShortPencils", P::amrTransShortPencils);
   RP::get("AMR.filterpasses", P::blurPassString);
   RP::get("adaptGPUWID", P::adaptGPUWID);
   RP::get("GPUallocations", P::GPUallocations);

   // We need the correct number of parameters for the AMR boxes
   if(   P::amrBoxNumber != (int)P::amrBoxHalfWidthX.size()
      || P::amrBoxNumber != (int)P::amrBoxHalfWidthY.size()
      || P::amrBoxNumber != (int)P::amrBoxHalfWidthZ.size()
      || P::amrBoxNumber != (int)P::amrBoxCenterX.size()
      || P::amrBoxNumber != (int)P::amrBoxCenterY.size()
      || P::amrBoxNumber != (int)P::amrBoxCenterZ.size()
      || P::amrBoxNumber != (int)P::amrBoxMaxLevel.size()
   ) {
      cerr << "AMR.number_of_boxes is set to " << P::amrBoxNumber << " so the same number of values is required for AMR.box_half_width_[xyz] and AMR.box_center_[xyz]." << endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
   }

   // We need the correct number of parameters for the allowed refinement boxes
   if(   P::refineBoxNumber != (int)P::refinementMinX.size()
      || P::refineBoxNumber != (int)P::refinementMinY.size()
      || P::refineBoxNumber != (int)P::refinementMinZ.size()
      || P::refineBoxNumber != (int)P::refinementMaxX.size()
      || P::refineBoxNumber != (int)P::refinementMaxY.size()
      || P::refineBoxNumber != (int)P::refinementMaxZ.size()
   ) {
      cerr << "AMR.number_of_refine_boxes is set to " << P::refineBoxNumber << " so the same number of values is required for AMR.refinementMin[XYZ] and AMR.refinementMax[XYZ]." << endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
   }

   // If we are in an AMR run we need to set up the filtering scheme.
   if (P::amrMaxSpatialRefLevel>0){
      bool isEmpty = blurPassString.size() == 0;
      if (!isEmpty){
         //sanity check=> user should define a pass for every level
         if ((int)blurPassString.size() != P::amrMaxSpatialRefLevel + 1) {
            cerr << "Filter Passes=" << blurPassString.size() << "\t" << "AMR Levels=" << P::amrMaxSpatialRefLevel + 1 << endl;
            cerr << "FilterPasses do not match AMR levels. \t" << " in " << __FILE__ << ":" << __LINE__ << endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
         }
         //sort the filtering passes per refLevel
         numPasses.clear();
         //Parse to a vector of ints
         for (auto pass : blurPassString){
            P::numPasses.push_back(stoi(pass));
         }
         sort(numPasses.begin(),numPasses.end(),greater<int>());
      }else{
         //here we will default to manually constructing the number of passes
         numPasses.clear();
         auto g_sequence=[](int size){
            int retval=1;
            while(size!=0){
               retval*=2;
               size-=1;
            }
            return retval;
         };
         int maxPasses=g_sequence(P::amrMaxSpatialRefLevel-1);
         for (int refLevel=0; refLevel<=P::amrMaxSpatialRefLevel; refLevel++){
            numPasses.push_back(maxPasses);
            maxPasses/=2;
         }
         //Overwrite passes for the highest refLevel. We do not want to filter there.
         numPasses.at(P::amrMaxSpatialRefLevel) = 0;
      }
         P::maxFilteringPasses = numPasses[0];
   }

   if (P::xmax < P::xmin || (P::ymax < P::ymin || P::zmax < P::zmin)) {
      cerr << "Box domain error!" << endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
   }

   // Set some parameter values.
   P::dx_ini = (P::xmax - P::xmin) / P::xcells_ini;
   P::dy_ini = (P::ymax - P::ymin) / P::ycells_ini;
   P::dz_ini = (P::zmax - P::zmin) / P::zcells_ini;

   RP::get("gridbuilder.dt", P::dt);

   RP::get("gridbuilder.t_max", P::t_max);
   RP::get("gridbuilder.timestep_max", P::tstep_max);

   RP::get("gridbuilder.dt_ceil", P::dt_ceil);

   if (P::dynamicTimestep)
      P::dt = 0.0; // if dynamic timestep then first dt is always 0

   // if we are restarting, t,t_min, tstep, tstep_min will be overwritten in readGrid
   P::t_min = 0;
   P::t = P::t_min;
   P::tstep_min = 0;
   P::tstep = P::tstep_min;

   // Get field solver parameters
   RP::get("fieldsolver.maxWaveVelocity", P::maxWaveVelocity);
   RP::get("fieldsolver.maxSubcycles", P::maxFieldSolverSubcycles);
   RP::get("fieldsolver.resistivity", P::resistivity);
   RP::get("fieldsolver.diffusiveEterms", P::fieldSolverDiffusiveEterms);
   RP::get("fieldsolver.ohmHallTerm", P::ohmHallTerm);
   RP::get("fieldsolver.ohmGradPeTerm", P::ohmGradPeTerm); // Which order solver to use for fieldsolver eGradPe term (supported: 0 for off, 1 for first-order)
   RP::get("fieldsolver.electronTemperature", P::electronTemperature); // Electron temperature associated with anchor point, e.g. incoming solar wind
   RP::get("fieldsolver.electronDensity", P::electronDensity); // Electron density associated with anchor point, e.g. incoming solar wind
   RP::get("fieldsolver.electronPTindex", P::electronPTindex); // Polytropic index for solving electron equation of state to use in eGradPe term
   RP::get("fieldsolver.maxCFL", P::fieldSolverMaxCFL);
   RP::get("fieldsolver.minCFL", P::fieldSolverMinCFL);

   // manual FsGrid decomposition should be complete with three values. If at least one is set but all are not set, abort
   if ((RP::isSet("fieldsolver.manualFsGridDecompositionX")||RP::isSet("fieldsolver.manualFsGridDecompositionY")||RP::isSet("fieldsolver.manualFsGridDecompositionZ")) &&
        !(RP::isSet("fieldsolver.manualFsGridDecompositionX")&&RP::isSet("fieldsolver.manualFsGridDecompositionY")&&RP::isSet("fieldsolver.manualFsGridDecompositionZ")) ) {
      cerr << "ERROR all of fieldsolver.manualFsGridDecompositionX,Y,Z should be defined." << endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
   }
   RP::get("fieldsolver.manualFsGridDecompositionX", temp_task_t);
   P::manualFsGridDecomposition[0] = temp_task_t;
   RP::get("fieldsolver.manualFsGridDecompositionY", temp_task_t);
   P::manualFsGridDecomposition[1] = temp_task_t;
   RP::get("fieldsolver.manualFsGridDecompositionZ", temp_task_t);
   P::manualFsGridDecomposition[2] = temp_task_t;



   // Get Vlasov solver parameters
   RP::get("vlasovsolver.maxSlAccelerationRotation", P::maxSlAccelerationRotation);
   RP::get("vlasovsolver.maxSlAccelerationSubcycles", P::maxSlAccelerationSubcycles);
   RP::get("vlasovsolver.maxCFL", P::vlasovSolverMaxCFL);
   RP::get("vlasovsolver.minCFL", P::vlasovSolverMinCFL);
   RP::get("vlasovsolver.GhostTranslate",P::vlasovSolverGhostTranslate);
   RP::get("vlasovsolver.GhostTranslateExtent",P::vlasovSolverGhostTranslateExtent);
   RP::get("vlasovsolver.accelerateMaxwellianBoundaries",  P::vlasovAccelerateMaxwellianBoundaries);
   if (P::vlasovSolverGhostTranslate==true) {
      if (myRank == MASTER_RANK) {
         logFile<<"Performing spatial translation using ghost cell information with coalesced MPI updates."<<endl;
      }
      if (P::vlasovSolverGhostTranslateExtent == 0) {
         P::vlasovSolverGhostTranslateExtent = VLASOV_STENCIL_WIDTH+1;
      } else {
         if (P::vlasovSolverGhostTranslateExtent == VLASOV_STENCIL_WIDTH+1) {
            if (myRank == MASTER_RANK) {
               logFile<<"Ghost translating full stencil of size "<<P::vlasovSolverGhostTranslateExtent<<" around local domain."<<endl;
            }
         } else if (P::vlasovSolverGhostTranslateExtent > VLASOV_STENCIL_WIDTH+1) {
            P::vlasovSolverGhostTranslateExtent = VLASOV_STENCIL_WIDTH+1;
            if (myRank == MASTER_RANK) {
               logFile<<"Capping ghost translation stencil size to VLASOV_STENCIL_WIDTH+1 around local domain."<<endl;
            }
         } else {
            if (myRank == MASTER_RANK) {
               logFile<<"Ghost translating reduced stencil of size "<<P::vlasovSolverGhostTranslateExtent<<" around local domain."<<endl;
            }
         }
      }
   }
   // Get load balance parameters
   RP::get("loadBalance.algorithm", P::loadBalanceAlgorithm);
   loadBalanceOptions["IMBALANCE_TOL"] = "";
   RP::get("loadBalance.tolerance", loadBalanceOptions["IMBALANCE_TOL"]);
   RP::get("loadBalance.rebalanceInterval", P::rebalanceInterval);

   std::vector<std::string> loadBalanceKeys;
   std::vector<std::string> loadBalanceValues;
   RP::get("loadBalance.optionKey", loadBalanceKeys);
   RP::get("loadBalance.optionValue", loadBalanceValues);
   if (loadBalanceKeys.size() != loadBalanceValues.size()) {
      if (myRank == MASTER_RANK) {
         cerr << "WARNING the number of load balance keys and values do not match. Disregarding these options." << endl;
      }
   } else {
      for (size_t i = 0; i < loadBalanceKeys.size(); ++i) {
         loadBalanceOptions[loadBalanceKeys[i]] = loadBalanceValues[i];
      }
   }

   // Get output variable parameters
   RP::get("variables.output", P::outputVariableList);
   RP::get("variables.diagnostic", P::diagnosticVariableList);

   // Insert vg_f_saved to the list if necessary
   if(includefSaved) {
      P::outputVariableList.push_back("vg_f_saved");
   }

   // Filter duplicate variable names
   set<string> dummy(P::outputVariableList.begin(), P::outputVariableList.end());
   P::outputVariableList.clear();
   P::outputVariableList.insert(P::outputVariableList.end(), dummy.begin(), dummy.end());
   dummy.clear();

   dummy.insert(P::diagnosticVariableList.begin(), P::diagnosticVariableList.end());
   P::diagnosticVariableList.clear();
   P::diagnosticVariableList.insert(P::diagnosticVariableList.end(), dummy.begin(), dummy.end());

   // Get parameters related to bailout
   RP::get("bailout.write_restart", P::bailout_write_restart);
   RP::get("bailout.min_dt", P::bailout_min_dt);
   RP::get("bailout.max_memory", P::bailout_max_memory);
   RP::get("bailout.velocity_space_wall_block_margin", P::bailout_velocity_space_wall_margin);
   if(P::bailout_velocity_space_wall_margin > MAX_BLOCKS_PER_DIM / 2 && myRank == MASTER_RANK) {
      std::cerr << "bailout.velocity_space_wall_block_margin is larger than 0.5 * MAX_BLOCKS_PER_DIM, aborting." << std::endl;
      abort();
   }

   for (size_t s = 0; s < P::systemWriteName.size(); ++s) {
      P::systemWrites.push_back(0);
   }

   RP::get("PAD.enable", P::artificialPADiff);   
   RP::get("PAD.coefficient", P::PADcoefficient);
   RP::get("PAD.CFL",P::PADCFL);
   RP::get("PAD.vbins",P::PADvbins);
   RP::get("PAD.mubins",P::PADmubins);
   RP::get("PAD.file",P::PADnu0);
   RP::get("PAD.fudge",P::PADfudge);

   RP::get("fieldtracing.fieldLineTracer", tracerString);
   RP::get("fieldtracing.tracer_max_allowed_error", FieldTracing::fieldTracingParameters.max_allowed_error);
   RP::get("fieldtracing.tracer_max_attempts", FieldTracing::fieldTracingParameters.max_field_tracer_attempts);
   RP::get("fieldtracing.tracer_min_dx", FieldTracing::fieldTracingParameters.min_tracer_dx_full_box);
   RP::get("fieldtracing.fullbox_max_incomplete_cells", FieldTracing::fieldTracingParameters.fullbox_max_incomplete_cells);
   RP::get("fieldtracing.fluxrope_max_incomplete_cells", FieldTracing::fieldTracingParameters.fluxrope_max_incomplete_cells);
   RP::get("fieldtracing.fullbox_and_fluxrope_max_absolute_distance_to_trace", FieldTracing::fieldTracingParameters.fullbox_and_fluxrope_max_distance);
   RP::get("fieldtracing.use_reconstruction_cache", FieldTracing::fieldTracingParameters.useCache);
   RP::get("fieldtracing.fluxrope_max_curvature_radii_to_trace", FieldTracing::fieldTracingParameters.fluxrope_max_curvature_radii_to_trace);
   RP::get("fieldtracing.fluxrope_max_curvature_radii_extent", FieldTracing::fieldTracingParameters.fluxrope_max_curvature_radii_extent);
   RP::get("fieldtracing.min_allowed_x", FieldTracing::fieldTracingParameters.x_min);
   RP::get("fieldtracing.min_allowed_y", FieldTracing::fieldTracingParameters.y_min);
   RP::get("fieldtracing.min_allowed_z", FieldTracing::fieldTracingParameters.z_min);
   RP::get("fieldtracing.max_allowed_x", FieldTracing::fieldTracingParameters.x_max);
   RP::get("fieldtracing.max_allowed_y", FieldTracing::fieldTracingParameters.y_max);
   RP::get("fieldtracing.max_allowed_z", FieldTracing::fieldTracingParameters.z_max);
   FieldTracing::fieldTracingParameters.x_min = max((double)FieldTracing::fieldTracingParameters.x_min, P::xmin + 4.0*P::dx_ini);
   FieldTracing::fieldTracingParameters.y_min = max((double)FieldTracing::fieldTracingParameters.y_min, P::ymin + 4.0*P::dy_ini);
   FieldTracing::fieldTracingParameters.z_min = max((double)FieldTracing::fieldTracingParameters.z_min, P::zmin + 4.0*P::dz_ini);
   FieldTracing::fieldTracingParameters.x_max = min((double)FieldTracing::fieldTracingParameters.x_max, P::xmax - 4.0*P::dx_ini);
   FieldTracing::fieldTracingParameters.y_max = min((double)FieldTracing::fieldTracingParameters.y_max, P::ymax - 4.0*P::dy_ini);
   FieldTracing::fieldTracingParameters.z_max = min((double)FieldTracing::fieldTracingParameters.z_max, P::zmax - 4.0*P::dz_ini);

   if(tracerString == "Euler") {
      FieldTracing::fieldTracingParameters.tracingMethod = FieldTracing::Euler;
   } else if (tracerString == "ADPT_Euler") {
      FieldTracing::fieldTracingParameters.tracingMethod = FieldTracing::ADPT_Euler;
   } else if (tracerString == "BS") {
      FieldTracing::fieldTracingParameters.tracingMethod = FieldTracing::BS;
   } else if (tracerString == "DP") {
      FieldTracing::fieldTracingParameters.tracingMethod = FieldTracing::DPrince;
   } else {
      cerr << __FILE__ << ":" << __LINE__ << " ERROR: Unknown value for fieldtracing.fieldLineTracer: " << tracerString << endl;
      abort();
   }
}


