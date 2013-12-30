/*
This file is part of Vlasiator.

Copyright 2010, 2011, 2012, 2013 Finnish Meteorological Institute
*/

#include "parameters.h"
#include "readparameters.h"
#include <limits>

#ifndef NAN
#define NAN 0
#endif


using namespace std;

typedef Parameters P;

//Using numeric_limits<Real>::max() leads to FP exceptions inside boost programoptions, use a slightly smaller value to avoid...

const Real LARGE_REAL=1e20;
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

Real P::vxmin = NAN;
Real P::vxmax = NAN;
Real P::vymin = NAN;
Real P::vymax = NAN;
Real P::vzmin = NAN;
Real P::vzmax = NAN;

Real P::backstreamradius = NAN;
Real P::backstreamvx = NAN;
Real P::backstreamvy = NAN;
Real P::backstreamvz = NAN;

uint P::xcells_ini = numeric_limits<uint>::max();
uint P::ycells_ini = numeric_limits<uint>::max();
uint P::zcells_ini = numeric_limits<uint>::max();
uint P::vxblocks_ini = numeric_limits<uint>::max();
uint P::vyblocks_ini = numeric_limits<uint>::max();
uint P::vzblocks_ini = numeric_limits<uint>::max();

Real P::q = NAN;
Real P::m = NAN;
Real P::q_per_m = NAN;
Real P::t = 0;
Real P::t_min = 0;
Real P::t_max = LARGE_REAL;
Real P::dt = NAN;
Real P::vlasovSolverMaxCFL = NAN;
Real P::vlasovSolverMinCFL = NAN;
Real P::fieldSolverMaxCFL = NAN;
Real P::fieldSolverMinCFL = NAN;


luint P::tstep = 0;
luint P::tstep_min = 0;
luint P::tstep_max = 0;
luint P::diagnosticInterval = numeric_limits<uint>::max();
bool P::writeInitialState = true;

std::vector<std::string> P::systemWriteName; 
std::vector<Real> P::systemWriteTimeInterval;
std::vector<int> P::systemWriteDistributionWriteStride;
std::vector<int> P::systemWriteDistributionWriteXlineStride;
std::vector<int> P::systemWriteDistributionWriteYlineStride;
std::vector<int> P::systemWriteDistributionWriteZlineStride;
std::vector<int> P::systemWrites;
   

Real P::saveRestartWalltimeInterval = -1.0;
uint P::exitAfterRestarts = numeric_limits<uint>::max();
int P::restartStripeFactor = -1;

uint P::transmit = 0;

bool P::recalculateStencils = true;
bool P::propagateVlasov = true;
bool P::propagateField = true;

uint P::maxAccelerationSubsteps=1;
bool P::dynamicTimestep = true;

Real P::maxAlfvenVelocity = 0.0;
Real P::resistivity = NAN;
bool P::fieldSolverDiffusiveEterms = true;
bool P::ohmHallTerm = false;

Real P::sparseMinValue = NAN;
int P::sparseBlockAddWidthV = 1;

string P::restartFileName = string("");                
bool P::isRestart=false;
int P::writeAsFloat = false;
string P::loadBalanceAlgorithm = string("");
string P::loadBalanceTolerance = string("");
uint P::rebalanceInterval = numeric_limits<uint>::max();

Real P::loadBalanceAlpha = 1.0;
Real P::loadBalanceGamma = 0.0;

vector<string> P::outputVariableList;
vector<string> P::diagnosticVariableList;

string P::projectName = string("");

Real P::maxSlAccelerationRotation=10.0;
bool P::lorentzHallTerm=false;
Real P::lorentzHallMinimumRho=1.0;

bool Parameters::addParameters(){
   //the other default parameters we read through the add/get interface
   Readparameters::add("io.diagnostic_write_interval", "Write diagnostic output every arg time steps",numeric_limits<uint>::max());
   

   Readparameters::addComposing("io.system_write_t_interval", "Save the simulation every arg simulated seconds. Negative values disable writes.");
   Readparameters::addComposing("io.system_write_file_name", "Save the simulation to this filename series");
   Readparameters::addComposing("io.system_write_distribution_stride", "Every this many cells write out their velocity space. 0 is none.");
   Readparameters::addComposing("io.system_write_distribution_xline_stride", "Every this many lines of cells along the x direction write out their velocity space. 0 is none.");
   Readparameters::addComposing("io.system_write_distribution_yline_stride", "Every this many lines of cells along the y direction write out their velocity space. 0 is none.");
   Readparameters::addComposing("io.system_write_distribution_zline_stride", "Every this many lines of cells along the z direction write out their velocity space. 0 is none.");

   Readparameters::add("io.write_initial_state","Write initial state, not even the 0.5 dt propagation is done. Do not use for restarting. ",false);

   Readparameters::add("io.restart_walltime_interval","Save the complete simulation in given walltime intervals. Negative values disable writes.",-1.0);
   Readparameters::add("io.number_of_restarts","Exit the simulation after certain number of walltime-based restarts.",numeric_limits<uint>::max());
   Readparameters::add("io.write_restart_stripe_factor","Stripe factor for restar writing.", -1);
   Readparameters::add("io.write_as_float","If true, write in floats instead of doubles", false);
   
   Readparameters::add("propagate_field","Propagate magnetic field during the simulation",true);
   Readparameters::add("propagate_vlasov","Propagate distribution functions during the simulation",true);
   Readparameters::add("max_acceleration_substeps","Maximum number of  acceleration substeps that are allowed to be taken in acceleration. The default number of 1 disables substepping and the acceleration is always done in one step. A value of 0 has a special meaning, it activates unlimited substepping",1);
   Readparameters::add("dynamic_timestep","If true,  timestep is set based on  CFL limits (default)",true);
   Readparameters::add("project", "Specify the name of the project to use. Supported to date (20121112): Alfven Diffusion Dispersion Firehose Flowthrough Fluctuations harm1D KelvinHelmholtz Magnetosphere", "Fluctuations");

   Readparameters::add("restart.filename","Restart from this vlsv file. No restart if empty file.",string(""));     
   
   Readparameters::add("gridbuilder.x_min","Minimum value of the x-coordinate.","");
   Readparameters::add("gridbuilder.x_max","Minimum value of the x-coordinate.","");
   Readparameters::add("gridbuilder.y_min","Minimum value of the y-coordinate.","");
   Readparameters::add("gridbuilder.y_max","Minimum value of the y-coordinate.","");
   Readparameters::add("gridbuilder.z_min","Minimum value of the z-coordinate.","");
   Readparameters::add("gridbuilder.z_max","Minimum value of the z-coordinate.","");
   Readparameters::add("gridbuilder.x_length","Number of cells in x-direction in initial grid.","");
   Readparameters::add("gridbuilder.y_length","Number of cells in y-direction in initial grid.","");
   Readparameters::add("gridbuilder.z_length","Number of cells in z-direction in initial grid.","");
   Readparameters::add("gridbuilder.vx_min","Minimum value for velocity block vx-coordinates.","");
   Readparameters::add("gridbuilder.vx_max","Maximum value for velocity block vx-coordinates.","");
   Readparameters::add("gridbuilder.vy_min","Minimum value for velocity block vy-coordinates.","");
   Readparameters::add("gridbuilder.vy_max","Maximum value for velocity block vy-coordinates.","");
   Readparameters::add("gridbuilder.vz_min","Minimum value for velocity block vz-coordinates.","");
   Readparameters::add("gridbuilder.vz_max","Maximum value for velocity block vz-coordinates.","");
   Readparameters::add("gridbuilder.vx_length","Initial number of velocity blocks in vx-direction.","");
   Readparameters::add("gridbuilder.vy_length","Initial number of velocity blocks in vy-direction.","");
   Readparameters::add("gridbuilder.vz_length","Initial number of velocity blocks in vz-direction.","");
   
   Readparameters::add("gridbuilder.q","Charge of simulated particle species, in Coulombs.",1.60217653e-19);
   Readparameters::add("gridbuilder.m","Mass of simulated particle species, in kilograms.",1.67262171e-27);
   Readparameters::add("gridbuilder.dt","Initial timestep in seconds.",0.0);

   Readparameters::add("gridbuilder.t_max","Maximum simulation time, in seconds. If timestep_max limit is hit first this time will never be reached",LARGE_REAL);
   Readparameters::add("gridbuilder.timestep_max","Max. value for timesteps. If t_max limit is hit first, this step will never be reached",numeric_limits<uint>::max());
   
   // Field solver parameters
   Readparameters::add("fieldsolver.maxAlfvenVelocity", "Maximum Alfven velocity allowed in the fast MS velocity determination in m/s, default unlimited", LARGE_REAL);
   Readparameters::add("fieldsolver.resistivity", "Resistivity for the eta*J term in Ohm's law.", 0.0);
   Readparameters::add("fieldsolver.diffusiveEterms", "Enable diffusive terms in the computation of E",true);
   Readparameters::add("fieldsolver.ohmHallTerm", "Enable the Hall term in Ohm's law", false);
   Readparameters::add("fieldsolver.maxCFL","The maximum CFL limit for field propagation. Used to set timestep if dynamic_timestep is true.",0.5);
   Readparameters::add("fieldsolver.minCFL","The minimum CFL limit for field propagation. Used to set timestep if dynamic_timestep is true.",0.4);

   // Vlasov solver parameters
   Readparameters::add("vlasovsolver.maxSlAccelerationRotation","Maximum rotation angle allowed by the Semi-Lagrangian solver (Do not use too large values, test properly)",10.0);
   Readparameters::add("vlasovsolver.lorentzHallTerm", "Add JxB term to Lorentz force",true);
   Readparameters::add("vlasovsolver.lorentzHallMinimumRho", "Minimum rho value used for Hall term in Lorentz force. Default is very low and has no effect in practice.",1.0);
   Readparameters::add("vlasovsolver.maxCFL","The maximum CFL limit for vlasov-leveque propagation. Used to set timestep if dynamic_timestep is true.",0.9);
   Readparameters::add("vlasovsolver.minCFL","The minimum CFL limit for vlasov-leveque propagation. Used to set timestep if dynamic_timestep is true.",0.7);



   
   // Grid sparsity parameters
   Readparameters::add("sparse.minValue", "Minimum value of distribution function in any cell of a velocity block for the block to be considered to have contents", 0);
   Readparameters::add("sparse.blockAddWidthV", "Number of layers of blocks that are kept in velocity space around the blocks with content",1);
   // Load balancing parameters
   Readparameters::add("loadBalance.algorithm", "Load balancing algorithm to be used", std::string("RCB"));
   Readparameters::add("loadBalance.tolerance", "Load imbalance tolerance", std::string("1.05"));
   Readparameters::add("loadBalance.rebalanceInterval", "Load rebalance interval (steps)", 10);
   Readparameters::add("loadBalance.alpha", "alpha in LB weight = gamma + blocks * alpha",1.0);
   Readparameters::add("loadBalance.gamma", "gamma in LB weight = gamma + blocks * alpha",0);
   
   
// Output variable parameters
   Readparameters::addComposing("variables.output", "List of data reduction operators (DROs) to add to the grid file output. Each variable to be added has to be on a new line output = XXX. Available are B BackgroundB PerturbedB E Rho RhoV RhoLossAdjust RhoLossVelBoundary MPIrank Blocks BoundaryType BoundaryLayer VolE VolB Pressure PTensor derivs BVOLderivs MaxVdt MaxRdt MaxFieldsdt LBweight VelocitySubSteps.");
   Readparameters::addComposing("variables.diagnostic", "List of data reduction operators (DROs) to add to the diagnostic runtime output. Each variable to be added has to be on a new line diagnostic = XXX. Available (20121005) are Blocks FluxB FluxE Rho RhoLossAdjust RhoLossVelBoundary  MaxDistributionFunction MinDistributionFunction  BoundaryType BoundaryLayer  MaxVdt MaxRdt MaxFieldsdt LBweight.");
   Readparameters::add("variables.dr_backstream_vx", "Center coordinate for the maxwellian distribution. Used for calculating the backstream contriution for rho.", -500000.0);
   Readparameters::add("variables.dr_backstream_vy", "Center coordinate for the maxwellian distribution. Used for calculating the backstream contriution for rho.", 0.0);
   Readparameters::add("variables.dr_backstream_vz", "Center coordinate for the maxwellian distribution. Used for calculating the backstream contriution for rho.", 0.0);
   Readparameters::add("variables.dr_backstream_radius", "Radius of the maxwellian distribution. Used for calculating the backstream contriution for rho", 468621.0);

   return true;
}

bool Parameters::getParameters(){
   //get numerical values of the parameters

   Readparameters::get("io.diagnostic_write_interval", P::diagnosticInterval);
   Readparameters::get("io.system_write_t_interval", P::systemWriteTimeInterval);
   Readparameters::get("io.system_write_file_name", P::systemWriteName);
   Readparameters::get("io.system_write_distribution_stride", P::systemWriteDistributionWriteStride);
   Readparameters::get("io.system_write_distribution_xline_stride", P::systemWriteDistributionWriteXlineStride);
   Readparameters::get("io.system_write_distribution_xline_stride", P::systemWriteDistributionWriteYlineStride);
   Readparameters::get("io.system_write_distribution_xline_stride", P::systemWriteDistributionWriteZlineStride);
   //TODO, check that the systemWrite vectors are of equal length
   Readparameters::get("io.write_initial_state", P::writeInitialState);
   Readparameters::get("io.restart_walltime_interval", P::saveRestartWalltimeInterval);
   Readparameters::get("io.number_of_restarts", P::exitAfterRestarts);
   Readparameters::get("io.write_restart_stripe_factor", P::restartStripeFactor);
   Readparameters::get("io.write_as_float", P::writeAsFloat);
   
   Readparameters::get("propagate_field",P::propagateField);
   Readparameters::get("propagate_vlasov",P::propagateVlasov);
   Readparameters::get("max_acceleration_substeps",P::maxAccelerationSubsteps);
   Readparameters::get("dynamic_timestep",P::dynamicTimestep);
   Readparameters::get("restart.filename",P::restartFileName);
   P::isRestart=(P::restartFileName!=string(""));
   
   Readparameters::get("project", projectName);
   
   /*get numerical values, let Readparameters handle the conversions*/
   Readparameters::get("gridbuilder.x_min",P::xmin);
   Readparameters::get("gridbuilder.x_max",P::xmax);
   Readparameters::get("gridbuilder.y_min",P::ymin);
   Readparameters::get("gridbuilder.y_max",P::ymax);
   Readparameters::get("gridbuilder.z_min",P::zmin);
   Readparameters::get("gridbuilder.z_max",P::zmax);
   Readparameters::get("gridbuilder.x_length",P::xcells_ini);
   Readparameters::get("gridbuilder.y_length",P::ycells_ini);
   Readparameters::get("gridbuilder.z_length",P::zcells_ini);
   Readparameters::get("gridbuilder.vx_min",P::vxmin);
   Readparameters::get("gridbuilder.vx_max",P::vxmax);
   Readparameters::get("gridbuilder.vy_min",P::vymin);
   Readparameters::get("gridbuilder.vy_max",P::vymax);
   Readparameters::get("gridbuilder.vz_min",P::vzmin);
   Readparameters::get("gridbuilder.vz_max",P::vzmax);
   Readparameters::get("gridbuilder.vx_length",P::vxblocks_ini);
   Readparameters::get("gridbuilder.vy_length",P::vyblocks_ini);
   Readparameters::get("gridbuilder.vz_length",P::vzblocks_ini);
   
   if (P::xmax < P::xmin || (P::ymax < P::ymin || P::zmax < P::zmin)) return false;
   if (P::vxmax < P::vxmin || (P::vymax < P::vymin || P::vzmax < P::vzmin)) return false;
   
   // Set some parameter values. 
   P::dx_ini = (P::xmax-P::xmin)/P::xcells_ini;
   P::dy_ini = (P::ymax-P::ymin)/P::ycells_ini;
   P::dz_ini = (P::zmax-P::zmin)/P::zcells_ini;
   
   Readparameters::get("gridbuilder.q",P::q);
   Readparameters::get("gridbuilder.m",P::m);
   Readparameters::get("gridbuilder.dt",P::dt);

   Readparameters::get("gridbuilder.t_max",P::t_max);
   Readparameters::get("gridbuilder.timestep_max",P::tstep_max);

   if(P::dynamicTimestep)
      P::dt=0.0; //if dynamic timestep then first dt is always 0 
   P::q_per_m = P::q/P::m;
   //if we are restarting, t,t_min, tstep, tstep_min will be overwritten in readGrid
   P::t_min=0;          
   P::t = P::t_min;
   P::tstep_min=0;
   P::tstep = P::tstep_min;
   
   // Get field solver parameters
   Readparameters::get("fieldsolver.maxAlfvenVelocity", P::maxAlfvenVelocity);
   Readparameters::get("fieldsolver.resistivity", P::resistivity);
   Readparameters::get("fieldsolver.diffusiveEterms", P::fieldSolverDiffusiveEterms);
   Readparameters::get("fieldsolver.ohmHallTerm", P::ohmHallTerm);
   Readparameters::get("fieldsolver.maxCFL",P::fieldSolverMaxCFL);
   Readparameters::get("fieldsolver.minCFL",P::fieldSolverMinCFL);
   // Get Vlasov solver parameters
   Readparameters::get("vlasovsolver.maxSlAccelerationRotation",P::maxSlAccelerationRotation);
   Readparameters::get("vlasovsolver.lorentzHallTerm", P::lorentzHallTerm);
   Readparameters::get("vlasovsolver.lorentzHallMinimumRho",P::lorentzHallMinimumRho);
   Readparameters::get("vlasovsolver.maxCFL",P::vlasovSolverMaxCFL);
   Readparameters::get("vlasovsolver.minCFL",P::vlasovSolverMinCFL);

   // Get sparsity parameters
   Readparameters::get("sparse.minValue", P::sparseMinValue);
   Readparameters::get("sparse.blockAddWidthV", P::sparseBlockAddWidthV); 
   // Get load balance parameters
   Readparameters::get("loadBalance.algorithm", P::loadBalanceAlgorithm);
   Readparameters::get("loadBalance.tolerance", P::loadBalanceTolerance);
   Readparameters::get("loadBalance.rebalanceInterval", P::rebalanceInterval);
   Readparameters::get("loadBalance.alpha", P::loadBalanceAlpha);
   Readparameters::get("loadBalance.gamma", P::loadBalanceGamma);

   // Get output variable parameters
   Readparameters::get("variables.output", P::outputVariableList);
   Readparameters::get("variables.diagnostic", P::diagnosticVariableList);

   //Get parameters related to calculating backstream contributions
   Readparameters::get("variables.dr_backstream_radius", P::backstreamradius);
   Readparameters::get("variables.dr_backstream_vx", P::backstreamvx);
   Readparameters::get("variables.dr_backstream_vy", P::backstreamvy);
   Readparameters::get("variables.dr_backstream_vz", P::backstreamvz);

   
   return true;
}
