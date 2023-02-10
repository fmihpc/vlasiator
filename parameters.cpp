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

#ifndef NAN
#define NAN 0
#endif

using namespace std;

typedef Parameters P;

// Using numeric_limits<Real>::max() leads to FP exceptions inside boost programoptions, use a slightly smaller value to
// avoid...

const Real LARGE_REAL = 1e20;
// Define static members:
int P::geometry = geometry::XYZ6D;
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
Real P::dt = NAN;
Real P::vlasovSolverMaxCFL = NAN;
Real P::vlasovSolverMinCFL = NAN;
Real P::fieldSolverMaxCFL = NAN;
Real P::fieldSolverMinCFL = NAN;
uint P::fieldSolverSubcycles = 1;

bool P::amrTransShortPencils = false;

uint P::tstep = 0;
uint P::tstep_min = 0;
uint P::tstep_max = 0;
uint P::diagnosticInterval = numeric_limits<uint>::max();
bool P::writeInitialState = true;

bool P::meshRepartitioned = true;
bool P::prepareForRebalance = false;
vector<CellID> P::localCells;

vector<string> P::systemWriteName;
vector<string> P::systemWritePath;
vector<Real> P::systemWriteTimeInterval;
vector<int> P::systemWriteDistributionWriteStride;
vector<int> P::systemWriteDistributionWriteXlineStride;
vector<int> P::systemWriteDistributionWriteYlineStride;
vector<int> P::systemWriteDistributionWriteZlineStride;
vector<Real> P::systemWriteDistributionWriteShellRadius;
vector<int> P::systemWriteDistributionWriteShellStride;
vector<int> P::systemWrites;
vector<pair<string, string>> P::systemWriteHints;
vector<pair<string, string>> P::restartWriteHints;

Real P::saveRestartWalltimeInterval = -1.0;
uint P::exitAfterRestarts = numeric_limits<uint>::max();
uint64_t P::vlsvBufferSize = 0;
int P::restartStripeFactor = 0;
int P::systemStripeFactor = 0;
string P::restartWritePath = string("");

uint P::transmit = 0;

bool P::recalculateStencils = true;
bool P::propagateVlasovAcceleration = true;
bool P::propagateVlasovTranslation = true;
bool P::propagateField = true;

bool P::dynamicTimestep = true;

Real P::maxWaveVelocity = 0.0;
uint P::maxFieldSolverSubcycles = 0.0;
int P::maxSlAccelerationSubcycles = 0.0;
bool P::ResolvePlasmaPeriod = false;
int P::maxResonantSubcycleFactor = 100;
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
string P::loadBalanceTolerance = string("");
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

uint P::amrMaxVelocityRefLevel = 0;
Realf P::amrRefineLimit = 1.0;
Realf P::amrCoarsenLimit = 0.5;
string P::amrVelRefCriterion = string("");
uint P::amrMaxSpatialRefLevel = 0;
int P::maxFilteringPasses = 0;
uint P::amrBoxHalfWidthX = 1;
uint P::amrBoxHalfWidthY = 1;
uint P::amrBoxHalfWidthZ = 1;
Realf P::amrBoxCenterX = 0.0;
Realf P::amrBoxCenterY = 0.0;
Realf P::amrBoxCenterZ = 0.0;
vector<string> P::blurPassString;
std::vector<int> P::numPasses; //numpasses

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

   RP::add("io.write_initial_state",
           "Write initial state, not even the 0.5 dt propagation is done. Do not use for restarting. ", false);

   RP::add("io.restart_walltime_interval",
           "Save the complete simulation in given walltime intervals. Negative values disable writes.", -1.0);
   RP::add("io.number_of_restarts", "Exit the simulation after certain number of walltime-based restarts.",
           numeric_limits<uint>::max());
   RP::add("io.vlsv_buffer_size",
           "Buffer size passed to VLSV writer (bytes, up to uint64_t), default 0 as this is sensible on sisu", 0);
   RP::add("io.write_restart_stripe_factor", "Stripe factor for restart and initial grid writing. Default 0 to inherit.", 0);
   RP::add("io.write_system_stripe_factor", "Stripe factor for bulk file writing. Default 0 to inherit.", 0);
   RP::add("io.write_as_float", "If true, write in floats instead of doubles", false);
   RP::add("io.restart_write_path",
           "Path to the location where restart files should be written. Defaults to the local directory, also if the "
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
           "Shocktest Template test_fp testHall test_trans VelocityBox verificationLarmor",
           string(""));

   RP::add("restart.write_as_float", "If true, write restart fields in floats instead of doubles", false);
   RP::add("restart.filename", "Restart from this vlsv file. No restart if empty file.", string(""));

   RP::add("gridbuilder.geometry", "Simulation geometry XY4D,XZ4D,XY5D,XZ5D,XYZ6D", string("XYZ6D"));
   RP::add("gridbuilder.x_min", "Minimum value of the x-coordinate.", NAN);
   RP::add("gridbuilder.x_max", "Minimum value of the x-coordinate.", NAN);
   RP::add("gridbuilder.y_min", "Minimum value of the y-coordinate.", NAN);
   RP::add("gridbuilder.y_max", "Minimum value of the y-coordinate.", NAN);
   RP::add("gridbuilder.z_min", "Minimum value of the z-coordinate.", NAN);
   RP::add("gridbuilder.z_max", "Minimum value of the z-coordinate.", NAN);
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
           "Upstream electron temperature to be used for the electron pressure gradient term (K).", 0.0);
   RP::add("fieldsolver.electronDensity",
           "Upstream electron density to be used for the electron pressure gradient term (m^-3).", 0.0);
   RP::add("fieldsolver.electronPTindex",
           "Polytropic index for electron pressure gradient term. 0 is isobaric, 1 is isothermal, 1.667 is adiabatic "
           "electrons, ",
           0.0);
   RP::add("fieldsolver.maxCFL",
           "The maximum CFL limit for field propagation. Used to set timestep if dynamic_timestep is true.", 0.5);
   RP::add("fieldsolver.minCFL",
           "The minimum CFL limit for field propagation. Used to set timestep if dynamic_timestep is true.", 0.4);

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
   RP::add("vlasovsolver.ResolvePlasmaPeriod","Require semi-lagrangian solver to resolve plasma oscillation (default false)",false);
   RP::add("vlasovsolver.maxResonantSubcycleFactor","Resolve (electron) plasma oscillation-cyclotron resonance with additional subcycles, max multiplicative substep number (default 100)",100);


   // Load balancing parameters
   RP::add("loadBalance.algorithm", "Load balancing algorithm to be used", string("RCB"));
   RP::add("loadBalance.tolerance", "Load imbalance tolerance", string("1.05"));
   RP::add("loadBalance.rebalanceInterval", "Load rebalance interval (steps)", 10);

   // Output variable parameters
   // NOTE Do not remove the : before the list of variable names as this is parsed by tools/check_vlasiator_cfg.sh

   RP::addComposing("variables.output",
                    string() +
                        "List of data reduction operators (DROs) to add to the grid file output.  Each variable to be "
                        "added has to be on a new line output = XXX. Names are case insensitive.  " +
                        "Available (20210125): " + "fg_b fg_b_background fg_b_perturbed fg_e " +
                        "vg_rhom vg_rhoq populations_vg_rho " + "fg_rhom fg_rhoq " + "vg_v fg_v populations_vg_v " +
                        "populations_vg_moments_thermal populations_vg_moments_nonthermal " +
                        "populations_vg_effectivesparsitythreshold populations_vg_rho_loss_adjust " +
                        "populations_vg_energydensity populations_vg_precipitationdifferentialflux " +
                        "populations_vg_heatflux" +  
                        "vg_maxdt_acceleration vg_maxdt_translation populations_vg_maxdt_acceleration "
                        "populations_vg_maxdt_translation " +
                        "fg_maxdt_fieldsolver " + "vg_rank fg_rank fg_amr_level vg_loadbalance_weight " +
                        "vg_eje vg_rhoqe fg_rhoqe "+
                        "vg_boundarytype fg_boundarytype vg_boundarylayer fg_boundarylayer " +
                        "populations_vg_blocks vg_f_saved " + "populations_vg_acceleration_subcycles " +
                        "vg_e_vol fg_e_vol " +
                        "fg_e_hall fg_e_gradpe vg_e_gradpe fg_b_vol vg_b_vol vg_b_background_vol vg_b_perturbed_vol " +
                        "vg_pressure fg_pressure populations_vg_ptensor " + "b_vol_derivatives " +
                        "ig_fac ig_latitude ig_cellarea ig_upmappedarea ig_sigmap ig_sigmah ig_rhom " +
                        "ig_electronTemp ig_potential ig_solverinternals ig_upmappedcodecoords ig_upmappedb ig_potential"+
                        "ig_inplanecurrent ig_e"+
                        "vg_gridcoordinates fg_gridcoordinates meshdata");

   RP::addComposing(
       "variables_deprecated.output",
       string() + "List of deprecated names for data reduction operators (DROs). Names are case insensitive. " +
           "Available (20201211): " + "B BackgroundB fg_BackgroundB PerturbedB fg_PerturbedB " + "E " +
           "Rhom Rhoq populations_Rho " + "V populations_V " +
				"eje rhoqe "+
           "populations_moments_Backstream populations_moments_NonBackstream " +
           "populations_moments_thermal populations_moments_nonthermal " +
           "populations_minvalue populations_EffectiveSparsityThreshold populations_RhoLossAdjust "
           "populations_rho_loss_adjust" +
           "populations_EnergyDensity populations_PrecipitationFlux populations_precipitationdifferentialflux" +
           "LBweight vg_lbweight vg_loadbalanceweight MaxVdt MaxRdt populations_MaxVdt populations_MaxRdt " +
           "populations_maxdt_acceleration populations_maxdt_translation MaxFieldsdt fg_maxfieldsdt" +
           "MPIrank FsGridRank " + "FsGridBoundaryType BoundaryType FsGridBoundaryLayer BoundaryLayer " +
           "populations_Blocks fSaved vg_fsaved" + "populations_accSubcycles populations_acceleration_subcycles" +
           "VolE vg_VolE Evol E_vol fg_VolE fg_Evol " +
           "HallE fg_HallE GradPeE e_gradpe VolB vg_VolB fg_VolB B_vol Bvol vg_Bvol fg_volB fg_Bvol" +
           "BackgroundVolB PerturbedVolB " + "Pressure vg_Pressure fg_Pressure populations_PTensor " +
           "BVOLderivs b_vol_derivs");

   // NOTE Do not remove the : before the list of variable names as this is parsed by tools/check_vlasiator_cfg.sh
   RP::addComposing("variables.diagnostic",
                    string() +
                        "List of data reduction operators (DROs) to add to the diagnostic runtime output. Each "
                        "variable to be added has to be on a new line diagnostic = XXX. Names are case insensitive. " +
                        "Available (20201111): " + "populations_vg_blocks " +
                        "vg_rhom populations_vg_rho_loss_adjust " + "vg_loadbalance_weight vg_eje " +
                        "vg_maxdt_acceleration vg_maxdt_translation " + "fg_maxdt_fieldsolver " +
                        "populations_vg_maxdt_acceleration populations_vg_maxdt_translation " +
                        "populations_vg_maxdistributionfunction populations_vg_mindistributionfunction");

   RP::addComposing("variables_deprecated.diagnostic",
                    string() +
                        "List of deprecated data reduction operators (DROs) to add to the diagnostic runtime output. "
                        "Names are case insensitive. " +
                        "Available (20201111): " + "rhom populations_rholossadjust populations_rho_loss_adjust " +
                        "populations_blocks lbweight loadbalance_weight " + "vg_lbweight vg_loadbalanceweight eje " +
                        "maxvdt maxdt_acceleration " + "maxrdt maxdt_translation " +
                        "populations_maxvdt populations_maxrdt " +
                        "populations_maxdt_acceleration populations_maxdt_translation " +
                        "populations_maxdistributionfunction populations_mindistributionfunction " +
                        "maxfieldsdt maxdt_fieldsolver fg_maxfieldsdt");

   // bailout parameters
   RP::add("bailout.write_restart",
           "If 1, write a restart file on bailout. Gets reset when sending a STOP (1) or a KILL (0).", true);
   RP::add("bailout.min_dt", "Minimum time step below which bailout occurs (s).", 1e-6);
   RP::add("bailout.max_memory", "Maximum amount of memory used per node (in GiB) over which bailout occurs.",
           1073741824.);
   RP::add("bailout.velocity_space_wall_block_margin", "Distance from the velocity space limits in blocks, if the distribution function reaches that distance from the wall we bail out to avoid hitting the wall.", 1);

   // Refinement parameters
   RP::add("AMR.vel_refinement_criterion", "Name of the velocity refinement criterion", string(""));
   RP::add("AMR.max_velocity_level", "Maximum velocity mesh refinement level", (uint)0);
   RP::add("AMR.refine_limit",
           "If the refinement criterion function returns a larger value than this, block is refined", (Realf)1.0);
   RP::add("AMR.coarsen_limit",
           "If the refinement criterion function returns a smaller value than this, block can be coarsened",
           (Realf)0.5);
   RP::add("AMR.max_spatial_level", "Maximum spatial mesh refinement level", (uint)0);
   RP::add("AMR.box_half_width_x", "Half width of the box that is refined (for testing)", (uint)1);
   RP::add("AMR.box_half_width_y", "Half width of the box that is refined (for testing)", (uint)1);
   RP::add("AMR.box_half_width_z", "Half width of the box that is refined (for testing)", (uint)1);
   RP::add("AMR.box_center_x", "x coordinate of the center of the box that is refined (for testing)", 0.0);
   RP::add("AMR.box_center_y", "y coordinate of the center of the box that is refined (for testing)", 0.0);
   RP::add("AMR.box_center_z", "z coordinate of the center of the box that is refined (for testing)", 0.0);
   RP::add("AMR.transShortPencils", "if true, use one-cell pencils", false);
   RP::addComposing("AMR.filterpasses", string("AMR filter passes for each individual refinement level"));

   return true;
}

void Parameters::getParameters() {
   typedef Readparameters RP;
   // get numerical values of the parameters
   RP::get("io.diagnostic_write_interval", P::diagnosticInterval);
   RP::get("io.system_write_t_interval", P::systemWriteTimeInterval);
   RP::get("io.system_write_file_name", P::systemWriteName);
   RP::get("io.system_write_path", P::systemWritePath);
   RP::get("io.system_write_distribution_stride", P::systemWriteDistributionWriteStride);
   RP::get("io.system_write_distribution_xline_stride", P::systemWriteDistributionWriteXlineStride);
   RP::get("io.system_write_distribution_yline_stride", P::systemWriteDistributionWriteYlineStride);
   RP::get("io.system_write_distribution_zline_stride", P::systemWriteDistributionWriteZlineStride);
   RP::get("io.system_write_distribution_shell_radius", P::systemWriteDistributionWriteShellRadius);
   RP::get("io.system_write_distribution_shell_stride", P::systemWriteDistributionWriteShellStride);
   RP::get("io.write_initial_state", P::writeInitialState);
   RP::get("io.restart_walltime_interval", P::saveRestartWalltimeInterval);
   RP::get("io.number_of_restarts", P::exitAfterRestarts);
   RP::get("io.vlsv_buffer_size", P::vlsvBufferSize);
   RP::get("io.write_restart_stripe_factor", P::restartStripeFactor);
   RP::get("io.write_system_stripe_factor", P::systemStripeFactor);
   RP::get("io.restart_write_path", P::restartWritePath);
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

   RP::get("project", P::projectName);
   if (RP::helpRequested) {
      P::projectName = string("Magnetosphere");
   }

   /*get numerical values, let Readparameters handle the conversions*/
   string geometryString;

   RP::get("gridbuilder.geometry", geometryString);
   RP::get("gridbuilder.x_min", P::xmin);
   RP::get("gridbuilder.x_max", P::xmax);
   RP::get("gridbuilder.y_min", P::ymin);
   RP::get("gridbuilder.y_max", P::ymax);
   RP::get("gridbuilder.z_min", P::zmin);
   RP::get("gridbuilder.z_max", P::zmax);
   RP::get("gridbuilder.x_length", P::xcells_ini);
   RP::get("gridbuilder.y_length", P::ycells_ini);
   RP::get("gridbuilder.z_length", P::zcells_ini);

   RP::get("AMR.max_velocity_level", P::amrMaxVelocityRefLevel);
   RP::get("AMR.max_spatial_level", P::amrMaxSpatialRefLevel);
   RP::get("AMR.box_half_width_x", P::amrBoxHalfWidthX);
   RP::get("AMR.box_half_width_y", P::amrBoxHalfWidthY);
   RP::get("AMR.box_half_width_z", P::amrBoxHalfWidthZ);
   RP::get("AMR.box_center_x", P::amrBoxCenterX);
   RP::get("AMR.box_center_y", P::amrBoxCenterY);
   RP::get("AMR.box_center_z", P::amrBoxCenterZ);
   RP::get("AMR.vel_refinement_criterion", P::amrVelRefCriterion);
   RP::get("AMR.refine_limit", P::amrRefineLimit);
   RP::get("AMR.coarsen_limit", P::amrCoarsenLimit);
   RP::get("AMR.transShortPencils", P::amrTransShortPencils);
   RP::get("AMR.filterpasses", P::blurPassString);

   // If we are in an AMR run we need to set up the filtering scheme.
   if (P::amrMaxSpatialRefLevel>0){
      bool isEmpty = blurPassString.size() == 0;
      if (!isEmpty){
         //sanity check=> user should define a pass for every level
         if (blurPassString.size() != P::amrMaxSpatialRefLevel + 1) {
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
         for (uint refLevel=0; refLevel<=P::amrMaxSpatialRefLevel; refLevel++){
            numPasses.push_back(maxPasses);
            maxPasses/=2; 
         }
         //Overwrite passes for the highest refLevel. We do not want to filter there.
         numPasses[P::amrMaxSpatialRefLevel] = 0;
      }
         P::maxFilteringPasses = numPasses[0];
   }

   if (geometryString == "XY4D") {
      P::geometry = geometry::XY4D;
   } else if (geometryString == "XZ4D") {
      P::geometry = geometry::XZ4D;
   } else if (geometryString == "XY5D") {
      P::geometry = geometry::XY5D;
   } else if (geometryString == "XZ5D") {
      P::geometry = geometry::XZ5D;
   } else if (geometryString == "XYZ6D") {
      P::geometry = geometry::XYZ6D;
   } else {
      cerr << "Unknown simulation geometry " << geometryString << " in " << __FILE__ << ":" << __LINE__ << endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
   }

   if (P::amrCoarsenLimit >= P::amrRefineLimit) {
      cerr << "amrRefineLimit must be smaller than amrCoarsenLimit!" << endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
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
   RP::get("fieldsolver.ohmGradPeTerm", P::ohmGradPeTerm);
   RP::get("fieldsolver.electronTemperature", P::electronTemperature);
   RP::get("fieldsolver.electronDensity", P::electronDensity);
   RP::get("fieldsolver.electronPTindex", P::electronPTindex);
   RP::get("fieldsolver.maxCFL", P::fieldSolverMaxCFL);
   RP::get("fieldsolver.minCFL", P::fieldSolverMinCFL);
   // Get Vlasov solver parameters
   RP::get("vlasovsolver.maxSlAccelerationRotation", P::maxSlAccelerationRotation);
   RP::get("vlasovsolver.maxSlAccelerationSubcycles", P::maxSlAccelerationSubcycles);
   RP::get("vlasovsolver.maxCFL", P::vlasovSolverMaxCFL);
   RP::get("vlasovsolver.minCFL", P::vlasovSolverMinCFL);
   RP::get("vlasovsolver.accelerateMaxwellianBoundaries",  P::vlasovAccelerateMaxwellianBoundaries);
   RP::get("vlasovsolver.ResolvePlasmaPeriod",P::ResolvePlasmaPeriod);
   RP::get("vlasovsolver.maxResonantSubcycleFactor",P::maxResonantSubcycleFactor);


   // Get load balance parameters
   RP::get("loadBalance.algorithm", P::loadBalanceAlgorithm);
   RP::get("loadBalance.tolerance", P::loadBalanceTolerance);
   RP::get("loadBalance.rebalanceInterval", P::rebalanceInterval);

   // Get output variable parameters
   RP::get("variables.output", P::outputVariableList);
   RP::get("variables.diagnostic", P::diagnosticVariableList);

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

   for (size_t s = 0; s < P::systemWriteName.size(); ++s)
      P::systemWrites.push_back(0);
}
