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
#include <iostream>
#include "parameters.h"
#include "readparameters.h"
#include <limits>
#include <set>
#include <unistd.h>

#ifndef NAN
   #define NAN 0
#endif

using namespace std;

typedef Parameters P;

//Using numeric_limits<Real>::max() leads to FP exceptions inside boost programoptions, use a slightly smaller value to avoid...

const Real LARGE_REAL=1e20;
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

Real P::backstreamradius = NAN;
Real P::backstreamvx = NAN;
Real P::backstreamvy = NAN;
Real P::backstreamvz = NAN;

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
int P::fieldSolverSubcycles = 1;

uint P::tstep = 0;
uint P::tstep_min = 0;
uint P::tstep_max = 0;
uint P::diagnosticInterval = numeric_limits<uint>::max();
bool P::writeInitialState = true;

bool P::meshRepartitioned = true;
bool P::prepareForRebalance = false;
std::vector<CellID> P::localCells;

vector<string> P::systemWriteName;
vector<string> P::systemWritePath;
vector<Real> P::systemWriteTimeInterval;
vector<int> P::systemWriteDistributionWriteStride;
vector<int> P::systemWriteDistributionWriteXlineStride;
vector<int> P::systemWriteDistributionWriteYlineStride;
vector<int> P::systemWriteDistributionWriteZlineStride;
vector<int> P::systemWrites;

Real P::saveRestartWalltimeInterval = -1.0;
uint P::exitAfterRestarts = numeric_limits<uint>::max();
int P::restartStripeFactor = -1;
string P::restartWritePath = string("");

uint P::transmit = 0;

bool P::recalculateStencils = true;
bool P::propagateVlasovAcceleration = true;
bool P::propagateVlasovTranslation = true;
bool P::propagateField = true;
bool P::propagatePotential = false;

bool P::dynamicTimestep = true;

Real P::maxWaveVelocity = 0.0;
int P::maxFieldSolverSubcycles = 0.0;
int P::maxSlAccelerationSubcycles = 0.0;
Real P::resistivity = NAN;
bool P::fieldSolverDiffusiveEterms = true;
uint P::ohmHallTerm = 0;
uint P::ohmGradPeTerm = 0;
Real P::electronTemperature = 0.0;

Real P::sparseMinValue = NAN;
int  P::sparseDynamicAlgorithm = 0;

#warning sparseDynamicBulkValue1 needs to be defined per species
Real P::sparseDynamicBulkValue1 = 1;
#warning sparseDynamicBulkValue2 needs to be defined per species
Real P::sparseDynamicBulkValue2 = 1;
#warning sparseDynamicMinValue1 needs to be defined per species
Real P::sparseDynamicMinValue1 = 1;
#warning sparseDynamicMinValue2 needs to be defined per species
Real P::sparseDynamicMinValue2 = 1;

int P::sparseBlockAddWidthV = 1;
bool P::sparse_conserve_mass = false;

string P::restartFileName = string("");
bool P::isRestart=false;
int P::writeAsFloat = false;
string P::loadBalanceAlgorithm = string("");
string P::loadBalanceTolerance = string("");
uint P::rebalanceInterval = numeric_limits<uint>::max();

vector<string> P::outputVariableList;
vector<string> P::diagnosticVariableList;

string P::projectName = string("");

Real P::maxSlAccelerationRotation=10.0;
Real P::hallMinimumRho=1.0;

bool P::bailout_write_restart = false;
Real P::bailout_min_dt = NAN;
Real P::bailout_max_memory = 1073741824.;

uint P::amrMaxVelocityRefLevel = 0;
Realf P::amrRefineLimit = 1.0;
Realf P::amrCoarsenLimit = 0.5;
string P::amrVelRefCriterion = "";

bool Parameters::addParameters(){
   //the other default parameters we read through the add/get interface
   Readparameters::add("io.diagnostic_write_interval", "Write diagnostic output every arg time steps",numeric_limits<uint>::max());
   

   Readparameters::addComposing("io.system_write_t_interval", "Save the simulation every arg simulated seconds. Negative values disable writes. [Define for all groups.]");
   Readparameters::addComposing("io.system_write_file_name", "Save the simulation to this file name series. [Define for all groups.]");
   Readparameters::addComposing("io.system_write_path", "Save this series in this location. Default is ./ [Define for all groups or none.]");
   Readparameters::addComposing("io.system_write_distribution_stride", "Every this many cells write out their velocity space. 0 is none. [Define for all groups.]");
   Readparameters::addComposing("io.system_write_distribution_xline_stride", "Every this many lines of cells along the x direction write out their velocity space. 0 is none. [Define for all groups.]");
   Readparameters::addComposing("io.system_write_distribution_yline_stride", "Every this many lines of cells along the y direction write out their velocity space. 0 is none. [Define for all groups.]");
   Readparameters::addComposing("io.system_write_distribution_zline_stride", "Every this many lines of cells along the z direction write out their velocity space. 0 is none. [Define for all groups.]");

   Readparameters::add("io.write_initial_state","Write initial state, not even the 0.5 dt propagation is done. Do not use for restarting. ",false);

   Readparameters::add("io.restart_walltime_interval","Save the complete simulation in given walltime intervals. Negative values disable writes.",-1.0);
   Readparameters::add("io.number_of_restarts","Exit the simulation after certain number of walltime-based restarts.",numeric_limits<uint>::max());
   Readparameters::add("io.write_restart_stripe_factor","Stripe factor for restart writing.", -1);
   Readparameters::add("io.write_as_float","If true, write in floats instead of doubles", false);
   Readparameters::add("io.restart_write_path", "Path to the location where restart files should be written. Defaults to the local directory, also if the specified destination is not writeable.", string("./"));
   
   Readparameters::add("propagate_potential","Propagate electrostatic potential during the simulation",false);
   Readparameters::add("propagate_field","Propagate magnetic field during the simulation",true);
   Readparameters::add("propagate_vlasov_acceleration","Propagate distribution functions during the simulation in velocity space. If false, it is propagated with zero length timesteps.",true);
   Readparameters::add("propagate_vlasov_translation","Propagate distribution functions during the simulation in ordinary space. If false, it is propagated with zero length timesteps.",true);
   Readparameters::add("dynamic_timestep","If true,  timestep is set based on  CFL limits (default on)",true);
   Readparameters::add("hallMinimumRho", "Minimum rho value used for the Hall and electron pressure gradient terms in the Lorentz force and in the field solver. Default is very low and has no effect in practice.", 1.0);
   Readparameters::add("project", "Specify the name of the project to use. Supported to date (20150610): Alfven Diffusion Dispersion Distributions Firehose Flowthrough Fluctuations Harris KHB Larmor Magnetosphere Multipeak PoissonTest Riemann1 Shock Shocktest Template test_fp testHall test_trans VelocityBox verificationLarmor", string(""));

   Readparameters::add("restart.filename","Restart from this vlsv file. No restart if empty file.",string(""));
   
   Readparameters::add("gridbuilder.geometry","Simulation geometry XY4D,XZ4D,XY5D,XZ5D,XYZ6D",string("XYZ6D"));
   Readparameters::add("gridbuilder.x_min","Minimum value of the x-coordinate.","");
   Readparameters::add("gridbuilder.x_max","Minimum value of the x-coordinate.","");
   Readparameters::add("gridbuilder.y_min","Minimum value of the y-coordinate.","");
   Readparameters::add("gridbuilder.y_max","Minimum value of the y-coordinate.","");
   Readparameters::add("gridbuilder.z_min","Minimum value of the z-coordinate.","");
   Readparameters::add("gridbuilder.z_max","Minimum value of the z-coordinate.","");
   Readparameters::add("gridbuilder.x_length","Number of cells in x-direction in initial grid.","");
   Readparameters::add("gridbuilder.y_length","Number of cells in y-direction in initial grid.","");
   Readparameters::add("gridbuilder.z_length","Number of cells in z-direction in initial grid.","");
   
   Readparameters::add("gridbuilder.dt","Initial timestep in seconds.",0.0);

   Readparameters::add("gridbuilder.t_max","Maximum simulation time, in seconds. If timestep_max limit is hit first this time will never be reached",LARGE_REAL);
   Readparameters::add("gridbuilder.timestep_max","Max. value for timesteps. If t_max limit is hit first, this step will never be reached",numeric_limits<uint>::max());
   
   // Field solver parameters
   Readparameters::add("fieldsolver.maxWaveVelocity", "Maximum wave velocity allowed in the fastest velocity determination in m/s, default unlimited", LARGE_REAL);
   Readparameters::add("fieldsolver.maxSubcycles", "Maximum allowed field solver subcycles", 1);
   Readparameters::add("fieldsolver.resistivity", "Resistivity for the eta*J term in Ohm's law.", 0.0);
   Readparameters::add("fieldsolver.diffusiveEterms", "Enable diffusive terms in the computation of E",true);
   Readparameters::add("fieldsolver.ohmHallTerm", "Enable/choose spatial order of the Hall term in Ohm's law. 0: off, 1: 1st spatial order, 2: 2nd spatial order", 0);
   Readparameters::add("fieldsolver.ohmGradPeTerm", "Enable/choose spatial order of the electron pressure gradient term in Ohm's law. 0: off, 1: 1st spatial order.", 0);
   Readparameters::add("fieldsolver.electronTemperature", "Constant electron temperature to be used for the electron pressure gradient term (K).", 0.0);
   Readparameters::add("fieldsolver.maxCFL","The maximum CFL limit for field propagation. Used to set timestep if dynamic_timestep is true.",0.5);
   Readparameters::add("fieldsolver.minCFL","The minimum CFL limit for field propagation. Used to set timestep if dynamic_timestep is true.",0.4);

   // Vlasov solver parameters
   Readparameters::add("vlasovsolver.maxSlAccelerationRotation","Maximum rotation angle (degrees) allowed by the Semi-Lagrangian solver (Use >25 values with care)",25.0);
   Readparameters::add("vlasovsolver.maxSlAccelerationSubcycles","Maximum number of subcycles for acceleration",1);
   Readparameters::add("vlasovsolver.maxCFL","The maximum CFL limit for vlasov propagation in ordinary space. Used to set timestep if dynamic_timestep is true.",0.99);
   Readparameters::add("vlasovsolver.minCFL","The minimum CFL limit for vlasov propagation in ordinary space. Used to set timestep if dynamic_timestep is true.",0.8);
   
   // Grid sparsity parameters
   Readparameters::add("sparse.minValue", "Minimum value of distribution function in any cell of a velocity block for the block to be considered to have contents", 1);
   Readparameters::add("sparse.blockAddWidthV", "Number of layers of blocks that are kept in velocity space around the blocks with content",1);
   Readparameters::add("sparse.conserve_mass", "If true, then mass is conserved by scaling the dist. func. in the remaining blocks", false);
   Readparameters::add("sparse.dynamicAlgorithm", "Type of algorithm used for calculating the dynamic minValue; 0 = none, 1 = linear algorithm based on rho, 2 = linear algorithm based on Blocks, (Example linear algorithm: y = kx+b, where dynamicMinValue1=k*dynamicBulkValue1 + b, and dynamicMinValue2 = k*dynamicBulkValue2 + b", 0);
   Readparameters::add("sparse.dynamicMinValue1", "The minimum value for the dynamic minValue", 1);
   Readparameters::add("sparse.dynamicMinValue2", "The maximum value (value 2) for the dynamic minValue", 1);
   Readparameters::add("sparse.dynamicBulkValue1", "Minimum value for the dynamic algorithm range, so for example if dynamicAlgorithm=1 then for sparse.dynamicBulkValue1 = 1e3, sparse.dynamicBulkValue2=1e5, we apply the algorithm to cells for which 1e3<cell.rho<1e5", 0);
   Readparameters::add("sparse.dynamicBulkValue2", "Maximum value for the dynamic algorithm range, so for example if dynamicAlgorithm=1 then for sparse.dynamicBulkValue1 = 1e3, sparse.dynamicBulkValue2=1e5, we apply the algorithm to cells for which 1e3<cell.rho<1e5", 0);

   
   

   // Load balancing parameters
   Readparameters::add("loadBalance.algorithm", "Load balancing algorithm to be used", string("RCB"));
   Readparameters::add("loadBalance.tolerance", "Load imbalance tolerance", string("1.05"));
   Readparameters::add("loadBalance.rebalanceInterval", "Load rebalance interval (steps)", 10);
   
// Output variable parameters
   Readparameters::addComposing("variables.output", "List of data reduction operators (DROs) to add to the grid file output. Each variable to be added has to be on a new line output = XXX. Available are (20141218) B BackgroundB PerturbedB E Rho RhoBackstream RhoV RhoVBackstream RhoVNonBackstream PressureBackstream  PTensorBackstreamDiagonal PTensorNonBackstreamDiagonal PTensorBackstreamOffDiagonal PTensorNonBackstreamOffDiagonal PTensorBackstream PTensorNonBackstream RhoNonBackstream RhoLossAdjust RhoLossVelBoundary LBweight MaxVdt MaxRdt MaxFieldsdt accSubcycles MPIrank FsGridRank FsGridBoundaryType BoundaryType BoundaryLayer Blocks fSaved VolE HallE BackgroundBedge VolB BackgroundVolB PerturbedVolB Pressure PTensor derivs BVOLderivs.");
   Readparameters::addComposing("variables.diagnostic", "List of data reduction operators (DROs) to add to the diagnostic runtime output. Each variable to be added has to be on a new line diagnostic = XXX. Available (20141218) are FluxB FluxE Blocks Pressure Rho RhoLossAdjust RhoLossVelBoundary LBweight MaxVdt MaxRdt MaxFieldsdt MaxDistributionFunction MinDistributionFunction BoundaryType BoundaryLayer.");
   Readparameters::add("variables.dr_backstream_vx", "Center coordinate for the maxwellian distribution. Used for calculating the backstream contriution for rho.", -500000.0);
   Readparameters::add("variables.dr_backstream_vy", "Center coordinate for the maxwellian distribution. Used for calculating the backstream contriution for rho.", 0.0);
   Readparameters::add("variables.dr_backstream_vz", "Center coordinate for the maxwellian distribution. Used for calculating the backstream contriution for rho.", 0.0);
   Readparameters::add("variables.dr_backstream_radius", "Radius of the maxwellian distribution. Used for calculating the backstream contribution for rho", 468621.0);

   // bailout parameters
   Readparameters::add("bailout.write_restart", "If 1, write a restart file on bailout. Gets reset when sending a STOP (1) or a KILL (0).", true);
   Readparameters::add("bailout.min_dt", "Minimum time step below which bailout occurs (s).", 1e-6);
   Readparameters::add("bailout.max_memory", "Maximum amount of memory used per node (in GiB) over which bailout occurs.", 1073741824.);

   // Refinement parameters
   Readparameters::add("AMR.vel_refinement_criterion","Name of the velocity refinement criterion",string(""));
   Readparameters::add("AMR.max_velocity_level","Maximum velocity mesh refinement level",(uint)0);
   Readparameters::add("AMR.refine_limit","If the refinement criterion function returns a larger value than this, block is refined",(Realf)1.0);
   Readparameters::add("AMR.coarsen_limit","If the refinement criterion function returns a smaller value than this, block can be coarsened",(Realf)0.5);
   return true;
}

bool Parameters::getParameters(){
   //get numerical values of the parameters
   Readparameters::get("io.diagnostic_write_interval", P::diagnosticInterval);
   Readparameters::get("io.system_write_t_interval", P::systemWriteTimeInterval);
   Readparameters::get("io.system_write_file_name", P::systemWriteName);
   Readparameters::get("io.system_write_path", P::systemWritePath);
   Readparameters::get("io.system_write_distribution_stride", P::systemWriteDistributionWriteStride);
   Readparameters::get("io.system_write_distribution_xline_stride", P::systemWriteDistributionWriteXlineStride);
   Readparameters::get("io.system_write_distribution_yline_stride", P::systemWriteDistributionWriteYlineStride);
   Readparameters::get("io.system_write_distribution_zline_stride", P::systemWriteDistributionWriteZlineStride);
   Readparameters::get("io.write_initial_state", P::writeInitialState);
   Readparameters::get("io.restart_walltime_interval", P::saveRestartWalltimeInterval);
   Readparameters::get("io.number_of_restarts", P::exitAfterRestarts);
   Readparameters::get("io.write_restart_stripe_factor", P::restartStripeFactor);
   Readparameters::get("io.restart_write_path", P::restartWritePath);
   Readparameters::get("io.write_as_float", P::writeAsFloat);
   
   // Checks for validity of io and restart parameters
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
   const string prefix = string("./");
   if (access(&(P::restartWritePath[0]), W_OK) != 0) {
      if(myRank == MASTER_RANK) {
         cerr << "ERROR restart write path " << P::restartWritePath << " not writeable, defaulting to local directory." << endl;
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
   if ( P::systemWriteTimeInterval.size() != maxSize) {
      if(myRank == MASTER_RANK) {
         cerr << "ERROR io.system_write_t_interval should be defined for all file types." << endl;
      }
      return false;
   }
   if ( P::systemWriteName.size() != maxSize) {
      if(myRank == MASTER_RANK) {
      cerr << "ERROR io.system_write_file_name should be defined for all file types." << endl;
      }
      return false;
   }
   if ( P::systemWritePath.size() != maxSize && P::systemWritePath.size() != 0) {
      if(myRank == MASTER_RANK) {
         cerr << "ERROR io.system_write_path should be defined for all file types or none at all." << endl;
      }
      return false;
   }
   if ( P::systemWriteDistributionWriteStride.size() != maxSize) {
      if(myRank == MASTER_RANK) {
         cerr << "ERROR io.system_write_distribution_stride should be defined for all file types." << endl;
      }
      return false;
   }
   if ( P::systemWriteDistributionWriteXlineStride.size() != maxSize) {
      if(myRank == MASTER_RANK) {
         cerr << "ERROR io.system_write_distribution_xline_stride should be defined for all file types." << endl;
      }
      return false;
   }
   if ( P::systemWriteDistributionWriteYlineStride.size() != maxSize) {
      if(myRank == MASTER_RANK) {
         cerr << "ERROR io.system_write_distribution_yline_stride should be defined for all file types." << endl;
      }
      return false;
   }
   if ( P::systemWriteDistributionWriteZlineStride.size() != maxSize) {
      if(myRank == MASTER_RANK) {
         cerr << "ERROR io.system_write_distribution_zline_stride should be defined for all file types." << endl;
      }
      return false;
   }
   if ( P::systemWritePath.size() == 0 ) {
      for (uint i = 0; i < P::systemWriteName.size(); i++) {
         P::systemWritePath.push_back(string("./"));
      }
   } else {
      for (uint i = 0; i < P::systemWritePath.size(); i++) {
         if (access(&(P::systemWritePath.at(i)[0]), W_OK) != 0) {
            if(myRank == MASTER_RANK) {
               cerr << "ERROR " << P::systemWriteName.at(i) << " write path " << P::systemWritePath.at(i) << " not writeable, defaulting to local directory." << endl;
            }
            P::systemWritePath.at(i) = prefix;
         }
      }
   }
   
   Readparameters::get("propagate_field",P::propagateField);
   Readparameters::get("propagate_potential",P::propagatePotential);
   Readparameters::get("propagate_vlasov_acceleration",P::propagateVlasovAcceleration);
   Readparameters::get("propagate_vlasov_translation",P::propagateVlasovTranslation);
   Readparameters::get("dynamic_timestep",P::dynamicTimestep);
   Readparameters::get("hallMinimumRho",P::hallMinimumRho);
   Readparameters::get("restart.filename",P::restartFileName);
   P::isRestart=(P::restartFileName!=string(""));

   Readparameters::get("project", P::projectName);
 
   /*get numerical values, let Readparameters handle the conversions*/
   string geometryString;
   Readparameters::get("gridbuilder.geometry",geometryString);
   Readparameters::get("gridbuilder.x_min",P::xmin);
   Readparameters::get("gridbuilder.x_max",P::xmax);
   Readparameters::get("gridbuilder.y_min",P::ymin);
   Readparameters::get("gridbuilder.y_max",P::ymax);
   Readparameters::get("gridbuilder.z_min",P::zmin);
   Readparameters::get("gridbuilder.z_max",P::zmax);
   Readparameters::get("gridbuilder.x_length",P::xcells_ini);
   Readparameters::get("gridbuilder.y_length",P::ycells_ini);
   Readparameters::get("gridbuilder.z_length",P::zcells_ini);
   Readparameters::get("AMR.max_velocity_level",P::amrMaxVelocityRefLevel);
   Readparameters::get("AMR.vel_refinement_criterion",P::amrVelRefCriterion);
   Readparameters::get("AMR.refine_limit",P::amrRefineLimit);
   Readparameters::get("AMR.coarsen_limit",P::amrCoarsenLimit);
   
   if (geometryString == "XY4D") P::geometry = geometry::XY4D;
   else if (geometryString == "XZ4D") P::geometry = geometry::XZ4D;
   else if (geometryString == "XY5D") P::geometry = geometry::XY5D;
   else if (geometryString == "XZ5D") P::geometry = geometry::XZ5D;
   else if (geometryString == "XYZ6D") P::geometry = geometry::XYZ6D;
   else {
      cerr << "Unknown simulation geometry in " << __FILE__ << ":" << __LINE__ << endl;
      return false;
   }
   
   if (P::amrCoarsenLimit >= P::amrRefineLimit) return false;
   if (P::xmax < P::xmin || (P::ymax < P::ymin || P::zmax < P::zmin)) return false;
   
   // Set some parameter values. 
   P::dx_ini = (P::xmax-P::xmin)/P::xcells_ini;
   P::dy_ini = (P::ymax-P::ymin)/P::ycells_ini;
   P::dz_ini = (P::zmax-P::zmin)/P::zcells_ini;
   
   Readparameters::get("gridbuilder.dt",P::dt);
   
   Readparameters::get("gridbuilder.t_max",P::t_max);
   Readparameters::get("gridbuilder.timestep_max",P::tstep_max);
   
   if(P::dynamicTimestep)
      P::dt=0.0; //if dynamic timestep then first dt is always 0 
   
   //if we are restarting, t,t_min, tstep, tstep_min will be overwritten in readGrid
   P::t_min=0;
   P::t = P::t_min;
   P::tstep_min=0;
   P::tstep = P::tstep_min;
   
   // Get field solver parameters
   Readparameters::get("fieldsolver.maxWaveVelocity", P::maxWaveVelocity);
   Readparameters::get("fieldsolver.maxSubcycles", P::maxFieldSolverSubcycles);
   Readparameters::get("fieldsolver.resistivity", P::resistivity);
   Readparameters::get("fieldsolver.diffusiveEterms", P::fieldSolverDiffusiveEterms);
   Readparameters::get("fieldsolver.ohmHallTerm", P::ohmHallTerm);
   Readparameters::get("fieldsolver.ohmGradPeTerm", P::ohmGradPeTerm);
   Readparameters::get("fieldsolver.electronTemperature", P::electronTemperature);
   Readparameters::get("fieldsolver.maxCFL",P::fieldSolverMaxCFL);
   Readparameters::get("fieldsolver.minCFL",P::fieldSolverMinCFL);
   // Get Vlasov solver parameters
   Readparameters::get("vlasovsolver.maxSlAccelerationRotation",P::maxSlAccelerationRotation);
   Readparameters::get("vlasovsolver.maxSlAccelerationSubcycles",P::maxSlAccelerationSubcycles);
   Readparameters::get("vlasovsolver.maxCFL",P::vlasovSolverMaxCFL);
   Readparameters::get("vlasovsolver.minCFL",P::vlasovSolverMinCFL);
   
   // Get sparsity parameters
   Readparameters::get("sparse.minValue", P::sparseMinValue);
   Readparameters::get("sparse.blockAddWidthV", P::sparseBlockAddWidthV); 
   Readparameters::get("sparse.conserve_mass", P::sparse_conserve_mass);
   Readparameters::get("sparse.dynamicAlgorithm", P::sparseDynamicAlgorithm);
   Readparameters::get("sparse.dynamicBulkValue1", P::sparseDynamicBulkValue1);
   Readparameters::get("sparse.dynamicBulkValue2", P::sparseDynamicBulkValue2);
   Readparameters::get("sparse.dynamicMinValue1", P::sparseDynamicMinValue1);
   Readparameters::get("sparse.dynamicMinValue2", P::sparseDynamicMinValue2);

   
   // Get load balance parameters
   Readparameters::get("loadBalance.algorithm", P::loadBalanceAlgorithm);
   Readparameters::get("loadBalance.tolerance", P::loadBalanceTolerance);
   Readparameters::get("loadBalance.rebalanceInterval", P::rebalanceInterval);
   
   // Get output variable parameters
   Readparameters::get("variables.output", P::outputVariableList);
   Readparameters::get("variables.diagnostic", P::diagnosticVariableList);

   // Filter duplicate variable names
   set<string> dummy(P::outputVariableList.begin(),P::outputVariableList.end());
   P::outputVariableList.clear();
   P::outputVariableList.insert(P::outputVariableList.end(),dummy.begin(),dummy.end());
   dummy.clear();
   
   dummy.insert(P::diagnosticVariableList.begin(),P::diagnosticVariableList.end());
   P::diagnosticVariableList.clear();
   P::diagnosticVariableList.insert(P::diagnosticVariableList.end(),dummy.begin(),dummy.end());

   //Get parameters related to calculating backstream contributions
   Readparameters::get("variables.dr_backstream_radius", P::backstreamradius);
   Readparameters::get("variables.dr_backstream_vx", P::backstreamvx);
   Readparameters::get("variables.dr_backstream_vy", P::backstreamvy);
   Readparameters::get("variables.dr_backstream_vz", P::backstreamvz);
   
   // Get parameters related to bailout
   Readparameters::get("bailout.write_restart", P::bailout_write_restart);
   Readparameters::get("bailout.min_dt", P::bailout_min_dt);
   Readparameters::get("bailout.max_memory", P::bailout_max_memory);

   for (size_t s=0; s<P::systemWriteName.size(); ++s) P::systemWrites.push_back(0);
   
   return true;
}
