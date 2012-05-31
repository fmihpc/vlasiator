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
Real P::dt = NAN;
luint P::tstep = 0;
luint P::tstep_min = 0;
luint P::tsteps = 0;
luint P::saveRestartInterval = numeric_limits<uint>::max();
luint P::diagnosticInterval = numeric_limits<uint>::max();
luint P::saveInterval = numeric_limits<uint>::max();

bool P::save_spatial_grid;
bool P::save_velocity_grid;

uint P::transmit = 0;

bool P::recalculateStencils = true;
bool P::propagateVlasov = true;
bool P::propagateField = true;
uint P::splitMethod=1;

bool P::periodic_x = false;
bool P::periodic_y = false;
bool P::periodic_z = false;

Real P::sparseMinValue = NAN;
Real P::sparseMinAvgValue = NAN;
uint P::blockAdjustmentInterval = numeric_limits<uint>::max();

string P::loadBalanceAlgorithm = string("");
string P::loadBalanceTolerance = string("");
uint P::rebalanceInterval = numeric_limits<uint>::max();

vector<string> P::outputVariableList;
vector<string> P::diagnosticVariableList;

bool Parameters::addParameters(){
        //the other default parameters we read through the add/get interface
        Readparameters::add("save_interval", "Save the simulation every arg time steps",1);
	Readparameters::add("diagnostic_interval", "Write diagnostic output every arg time steps", 1);
	Readparameters::add("restart_interval","Save the complete simulation every arg time steps",numeric_limits<uint>::max());
        Readparameters::add("save_spatial_grid", "Save spatial cell averages for the whole simulation",true);
        Readparameters::add("save_velocity_grid","Save velocity grid from every spatial cell in the simulation",false);
        
        Readparameters::add("propagate_field","Propagate magnetic field during the simulation",true);
        Readparameters::add("propagate_vlasov","Propagate distribution functions during the simulation",true);
        Readparameters::add("split_method","Split method for splitting spatial/velocity space solvers. 0: first order, 1: strang splitting with half-steps for spatial space, 2: strang splitting with half-steps for velocity space",1);
        

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
        Readparameters::add("gridbuilder.periodic_x","If 'yes' the grid is periodic in x-direction. Defaults to 'no'.","no");
        Readparameters::add("gridbuilder.periodic_y","If 'yes' the grid is periodic in y-direction. Defaults to 'no'.","no");
        Readparameters::add("gridbuilder.periodic_z","If 'yes' the grid is periodic in z-direction. Defaults to 'no'.","no");
   
        Readparameters::add("gridbuilder.q","Charge of simulated particle species, in Coulombs.",numeric_limits<Real>::max());
        Readparameters::add("gridbuilder.m","Mass of simulated particle species, in kilograms.",numeric_limits<Real>::max());
        Readparameters::add("gridbuilder.dt","Timestep in seconds.",numeric_limits<Real>::max());
        Readparameters::add("gridbuilder.t_min","Simulation time at timestep 0, in seconds.",numeric_limits<Real>::max());
        Readparameters::add("gridbuilder.timestep","Timestep when grid is loaded. Defaults to value zero.",0);
        Readparameters::add("gridbuilder.max_timesteps","Max. value for timesteps. Defaults to value zero.",0);
   
   // Grid sparsity parameters
        Readparameters::add("sparse.minValue", "Minimum value of distribution function in any cell of a velocity block for the block to be considered to have contents", 1e-5);
        Readparameters::add("sparse.minAvgValue", "Minimum value of the average of distribution function within a velocity block for the block to be considered to have contents", 0.5e-5);
        Readparameters::add("sparse.blockAdjustmentInterval", "Block adjustment interval (steps)", 1);
   
        // Load balancing parameters
        Readparameters::add("loadBalance.algorithm", "Load balancing algorithm to be used", std::string("HYPERGRAPH"));
        Readparameters::add("loadBalance.tolerance", "Load imbalance tolerance", std::string("1.05"));
        Readparameters::add("loadBalance.rebalanceInterval", "Load rebalance interval (steps)", 10);
	
	// Output variable parameters
	Readparameters::addComposing("variables.output", "List of data reduction operators (DROs) to add to the grid file output. Each variable to be added has to be on a new line output = XXX. Available (20120531) are B E Rho RhoV RhoLossAdjust RhoLossVelBoundary MPIrank Blocks VolE VolB Pressure dBxdz.");
	Readparameters::addComposing("variables.diagnostic", "List of data reduction operators (DROs) to add to the diagnostic runtime output. Each variable to be added has to be on a new line diagnostic = XXX. Available (20120531) are Blocks FluxB Rho RhoLossAdjust RhoLossVelBoundary MaxVi.");
        return true;
}

bool Parameters::getParameters(){
   //get numerical values of the parameters
   
   Readparameters::get("save_interval", P::saveInterval);
   Readparameters::get("diagnostic_interval", P::diagnosticInterval);
   Readparameters::get("restart_interval", P::saveRestartInterval);
   Readparameters::get("save_spatial_grid", P::save_spatial_grid);
   Readparameters::get("save_velocity_grid", P::save_velocity_grid);
   Readparameters::get("propagate_field",P::propagateField);
   Readparameters::get("propagate_vlasov",P::propagateVlasov);
   Readparameters::get("split_method",P::splitMethod);
    
   
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
    
    std::string periodic_x,periodic_y,periodic_z;
    Readparameters::get("gridbuilder.periodic_x",periodic_x);
    Readparameters::get("gridbuilder.periodic_y",periodic_y);
    Readparameters::get("gridbuilder.periodic_z",periodic_z);
    P::periodic_x = false;
    P::periodic_y = false;
    P::periodic_z = false;
    if (periodic_x == "yes") P::periodic_x = true;
    if (periodic_y == "yes") P::periodic_y = true;
    if (periodic_z == "yes") P::periodic_z = true;
    
    // Set some parameter values. 
    P::dx_ini = (P::xmax-P::xmin)/P::xcells_ini;
    P::dy_ini = (P::ymax-P::ymin)/P::ycells_ini;
    P::dz_ini = (P::zmax-P::zmin)/P::zcells_ini;
    
   Real t_min;
   Readparameters::get("gridbuilder.q",P::q);
   Readparameters::get("gridbuilder.m",P::m);
   Readparameters::get("gridbuilder.dt",P::dt);
   Readparameters::get("gridbuilder.t_min",t_min);
   Readparameters::get("gridbuilder.timestep",P::tstep);
   Readparameters::get("gridbuilder.max_timesteps",P::tsteps);
   
   P::q_per_m = P::q/P::m;
   P::t = t_min + P::tstep*P::dt;
   P::tstep_min = P::tstep;
   
   // Get sparsity parameters
   Readparameters::get("sparse.minValue", P::sparseMinValue);
   Readparameters::get("sparse.minAvgValue", P::sparseMinAvgValue);
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


