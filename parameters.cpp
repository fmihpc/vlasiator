/*
This file is part of Vlasiator.

Copyright 2010, 2011, 2012 Finnish Meteorological Institute

Vlasiator is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3
as published by the Free Software Foundation.

Vlasiator is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include "parameters.h"
#include "readparameters.h"
#include <limits>
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
Real P::CFL_max = NAN;
Real P::CFL_min = NAN;

luint P::tstep = 0;
luint P::tstep_min = 0;
luint P::tstep_max = 0;
luint P::diagnosticInterval = numeric_limits<uint>::max();
Real P::saveRestartTimeInterval = -1.0;
Real P::saveSystemTimeInterval = -1.0;

uint P::transmit = 0;

bool P::recalculateStencils = true;
bool P::propagateVlasov = true;
bool P::propagateField = true;

uint P::maxAccelerationSubsteps=1;
bool P::dynamicTimestep = true;


Real P::sparseMinValue = NAN;
uint P::blockAdjustmentInterval = numeric_limits<uint>::max();

string P::restartFileName = string("");                
bool P::isRestart=false;
string P::loadBalanceAlgorithm = string("");
string P::loadBalanceTolerance = string("");
uint P::rebalanceInterval = numeric_limits<uint>::max();

vector<string> P::outputVariableList;
vector<string> P::diagnosticVariableList;

string P::projectName = string("");


bool Parameters::addParameters(){
   //the other default parameters we read through the add/get interface
   Readparameters::add("diagnostic_write_interval", "Write diagnostic output every arg time steps",numeric_limits<uint>::max());
   Readparameters::add("system_write_t_interval", "Save the simulation every arg simulated seconds. Negative values disable writes.",-1.0);
   Readparameters::add("restart_write_t_interval","Save the complete simulation every arg simulated seconds. Negative values disable writes.",-1.0);
   //TODO Readparameters::add("output.restart_walltime_interval","Save the complete simulation every arg wall-time seconds",numeric_limits<uint>::max());

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
   Readparameters::add("gridbuilder.CFL_max","The maximum CFL limit for propagation. Used to set timestep if dynamic_timestep is true. Also used to compute substeps in acceleration",0.9);
   Readparameters::add("gridbuilder.CFL_min","The minimum CFL limit for propagation. Used to set timestep if dynamic_timestep is true. Also used to compute substeps in acceleration",0.7);
   Readparameters::add("gridbuilder.t_max","Maximum simulation time, in seconds. If timestep_max limit is hit first this time will never be reached",LARGE_REAL);
   Readparameters::add("gridbuilder.timestep_max","Max. value for timesteps. If t_max limit is hit first, this step will never be reached",numeric_limits<uint>::max());
   
   // Field solver parameters
   
   // Grid sparsity parameters
   Readparameters::add("sparse.minValue", "Minimum value of distribution function in any cell of a velocity block for the block to be considered to have contents", 0);
   Readparameters::add("sparse.blockAdjustmentInterval", "Block adjustment interval (steps)", 1);
   
   // Load balancing parameters
   Readparameters::add("loadBalance.algorithm", "Load balancing algorithm to be used", std::string("RCB"));
   Readparameters::add("loadBalance.tolerance", "Load imbalance tolerance", std::string("1.05"));
   Readparameters::add("loadBalance.rebalanceInterval", "Load rebalance interval (steps)", 10);

// Output variable parameters
   Readparameters::addComposing("variables.output", "List of data reduction operators (DROs) to add to the grid file output. Each variable to be added has to be on a new line output = XXX. Available are B BackgroundB PerturbedB E Rho RhoV RhoLossAdjust RhoLossVelBoundary MPIrank Blocks BoundaryType VolE VolB Pressure PTensor Bderivs MaxVdt MaxRdt MaxFieldsdt LBweight.");
   Readparameters::addComposing("variables.diagnostic", "List of data reduction operators (DROs) to add to the diagnostic runtime output. Each variable to be added has to be on a new line diagnostic = XXX. Available (20121005) are Blocks FluxB FluxE Rho RhoLossAdjust RhoLossVelBoundary  MaxDistributionFunction MinDistributionFunction MaxVdt MaxRdt MaxFieldsdt LBweight.");

   
   return true;
}

bool Parameters::getParameters(){
   //get numerical values of the parameters

   Readparameters::get("diagnostic_write_interval", P::diagnosticInterval);
   Readparameters::get("system_write_t_interval", P::saveSystemTimeInterval);
   Readparameters::get("restart_write_t_interval", P::saveRestartTimeInterval);
   
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
   Readparameters::get("gridbuilder.CFL_max",P::CFL_max);
   Readparameters::get("gridbuilder.CFL_min",P::CFL_min);

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
   
   // Get sparsity parameters
   Readparameters::get("sparse.minValue", P::sparseMinValue);
   Readparameters::get("sparse.blockAdjustmentInterval", P::blockAdjustmentInterval);
   
   // Get load balance parameters
   Readparameters::get("loadBalance.algorithm", P::loadBalanceAlgorithm);
   Readparameters::get("loadBalance.tolerance", P::loadBalanceTolerance);
   Readparameters::get("loadBalance.rebalanceInterval", P::rebalanceInterval);
   
   // Get output variable parameters
   Readparameters::get("variables.output", P::outputVariableList);
   Readparameters::get("variables.diagnostic", P::diagnosticVariableList);
   
   return true;
}
