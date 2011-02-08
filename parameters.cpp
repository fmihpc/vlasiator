#include <boost/program_options.hpp>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>

#include "parameters.h"

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
uint P::tstep = 0;
uint P::tsteps = 0;
uint P::saveRestartInterval = numeric_limits<uint>::max();
uint P::diagnInterval = numeric_limits<uint>::max();

bool P::save_spatial_grid;
bool P::save_velocity_grid;

std::string P::solar_wind_file;

std::vector<Real> P::save_spatial_cells_x;
std::vector<Real> P::save_spatial_cells_y;
std::vector<Real> P::save_spatial_cells_z;

uint P::transmit = 0;

// Handles parameter processing from the user
namespace PO = boost::program_options;
//static PO::options_description options 
//  = PO::options_description("Usage: main [options (options given on the command line override options given everywhere else)], where options are:");

static bool initialized = false;
static PO::options_description* descriptions = NULL;
static PO::variables_map* variables = NULL;

static string global_config_file_name = "";
static string user_config_file_name = "";
static string run_config_file_name = "";

int Parameters::argc;
char** Parameters::argv;

/** Constructor for class Parameters. The constructor defines some 
 * default parameters and parses the input files.
 * @param argc Command line argc.
 * @param argv Command line argv.
 */
Parameters::Parameters(int argc, char* argv[]) {
   Parameters::argc = argc;
   Parameters::argv = argv;
   if (initialized == false) {
      descriptions = new PO::options_description("Usage: main [options (options given on the command line override options given everywhere else)], where options are:");
      variables = new PO::variables_map;
      initialized = true;

      addDefaultParameters();
   }
}

/** Add a new input parameter to Parameters. Note that Parameters::parse must be called 
 * in order for the input file(s) to be re-read.
 * @param name The name of the parameter, as given in the input file(s).
 * @param desc Description for the parameter.
 * @param var A variable where the value of the parameter is to be written.
 * @param defValue Default value for variable var.
 * @return If true, the new parameter was added successfully.
 */
bool Parameters::add(const string& name,const string& desc,char& var,const char& defValue) {
   if (initialized == false) return false;
   descriptions->add_options()(name.c_str(), PO::value<char>(&var)->default_value(defValue), desc.c_str());
   return true;
}

/** Add a new input parameter to Parameters. Note that Parameters::parse must be called
 * in order for the input file(s) to be re-read.
 * @param name The name of the parameter, as given in the input file(s).
 * @param desc Description for the parameter.
 * @param var A variable where the value of the parameter is to be written.
 * @param defValue Default value for variable var.
 * @return If true, the new parameter was added successfully.
 */
bool Parameters::add(const string& name,const string& desc,int& var,const int& defValue) {
   if (initialized == false) return false;
   descriptions->add_options()(name.c_str(), PO::value<int>(&var)->default_value(defValue), desc.c_str());
   return true;
}

/** Add a new input parameter to Parameters. Note that Parameters::parse must be called
 * in order for the input file(s) to be re-read.
 * @param name The name of the parameter, as given in the input file(s).
 * @param desc Description for the parameter.
 * @param var A variable where the value of the parameter is to be written.
 * @param defValue Default value for variable var.
 * @return If true, the new parameter was added successfully.
 */
bool Parameters::add(const string& name,const string& desc,unsigned int& var,const unsigned int& defValue) {
   if (initialized == false) return false;
   descriptions->add_options()(name.c_str(), PO::value<unsigned int>(&var)->default_value(defValue), desc.c_str());
   return true;
}

/** Add a new input parameter to Parameters. Note that Parameters::parse must be called
 * in order for the input file(s) to be re-read.
 * @param name The name of the parameter, as given in the input file(s).
 * @param desc Description for the parameter.
 * @param var A variable where the value of the parameter is to be written.
 * @param defValue Default value for variable var.
 * @return If true, the new parameter was added successfully.
 */
bool Parameters::add(const string& name,const string& desc,float& var,const float& defValue) {
   if (initialized == false) return false;
   descriptions->add_options()(name.c_str(), PO::value<float>(&var)->default_value(defValue), desc.c_str());
   return true;
}

/** Add a new input parameter to Parameters. Note that Parameters::parse must be called
 * in order for the input file(s) to be re-read.
 * @param name The name of the parameter, as given in the input file(s).
 * @param desc Description for the parameter.
 * @param var A variable where the value of the parameter is to be written.
 * @param defValue Default value for variable var.
 * @return If true, the new parameter was added successfully.
 */
bool Parameters::add(const string& name,const string& desc,double& var,const double& defValue) {
   if (initialized == false) return false;
   descriptions->add_options()(name.c_str(), PO::value<double>(&var)->default_value(defValue), desc.c_str());
   return true;
}

/** Add a new input parameter to Parameters. Note that Parameters::parse must be called
 * in order for the input file(s) to be re-read.
 * @param name The name of the parameter, as given in the input file(s).
 * @param desc Description for the parameter.
 * @param var A variable where the value of the parameter is to be written.
 * @param defValue Default value for variable var.
 * @param multitoken If true, this parameter is allowed to have multiple values. For example, 
 * command line parameter --email address1@mail.com address2@mail.com would be written to 
 * var as "address1@mail.com address2@mail.com".
 * @return If true, the new parameter was added successfully.
 */
bool Parameters::add(const string& name,const string& desc,std::string& var,const std::string& defValue,const bool& multitoken) {
   if (initialized == false) return false;
   if (multitoken == true)
     descriptions->add_options()(name.c_str(), PO::value<string>(&var)->default_value(defValue)->multitoken(), desc.c_str());
   else
     descriptions->add_options()(name.c_str(), PO::value<string>(&var)->default_value(defValue), desc.c_str());
   return true;
}

bool Parameters::addDefaultParameters() {
   if (initialized == false) return false;
   descriptions->add_options()
   ("help", "print this help message");
   
   // Parameters which set the names of the configuration file(s):
   descriptions->add_options()
   ("global_config", PO::value<string>(&global_config_file_name)->default_value(""),"read options from the global configuration file arg (relative to the current working directory). Options given in this file are overridden by options given in the user's and run's configuration files and by options given in environment variables (prefixed with MAIN_) and the command line")
   ("user_config", PO::value<string>(&user_config_file_name)->default_value(""), "read options from the user's configuration file arg (relative to the current working directory). Options given in this file override options given in the global configuration file. Options given in this file are overridden by options given in the run's configuration file and by options given in environment variables (prefixed with MAIN_) and the command line")
   ("run_config", PO::value<string>(&run_config_file_name)->default_value(""), "read options from the run's configuration file arg (relative to the current working directory). Options given in this file override options given in the user's and global configuration files. Options given in this override options given in the user's and global configuration files. Options given in this file are overridden by options given in environment variables (prefixed with MAIN_) and the command line");
   
   // Parameters related to the spatial grid 
   // DEPRECATED: These are removed in the future. The correct place 
   // to define (and cache) these variables is GridBuilder.
   descriptions->add_options()
   ("x_min", PO::value<Real>(&P::xmin)->default_value(0), "Starting corner of spatial grid in x direction in meters")
   ("x_max", PO::value<Real>(&P::xmax)->default_value(1), "End corner of spatial grid in x direction in meters")
   ("y_min", PO::value<Real>(&P::ymin)->default_value(0), "Starting corner of spatial grid in y direction in meters")
   ("y_max", PO::value<Real>(&P::ymax)->default_value(1), "End corner of spatial grid in y direction in meters")
   ("z_min", PO::value<Real>(&P::zmin)->default_value(0), "Starting corner of spatial grid in z direction in meters")
   ("z_max", PO::value<Real>(&P::zmax)->default_value(1), "End corner of spatial grid in z direction in meters")
   ("x_length", PO::value<uint>(&P::xcells_ini)->default_value(1), "Length of the spatial grid in unrefined cells in the x direction")
   ("y_length", PO::value<uint>(&P::ycells_ini)->default_value(1), "Length of the spatial grid in unrefined cells in the y direction")
   ("z_length", PO::value<uint>(&P::zcells_ini)->default_value(1), "Length of the spatial grid in unrefined cells in the z direction");
	     
   // Parameters related to velocity block grid.
   // DEPRECATED: These will be moved to somewhere else in the future.
   descriptions->add_options()
   ("vx_min", PO::value<Real>(&P::vxmin)->default_value(-1), "Starting corner of velocity grids in x direction in m/s")
   ("vx_max", PO::value<Real>(&P::vxmax)->default_value(1), "End corner of velocity grids in x direction in m/s")
   ("vy_min", PO::value<Real>(&P::vymin)->default_value(-1), "Starting corner of velocity grids in y direction in m/s")
   ("vy_max", PO::value<Real>(&P::vymax)->default_value(1), "End corner of velocity grids in y direction in m/s")
   ("vz_min", PO::value<Real>(&P::vzmin)->default_value(-1), "Starting corner of velocity grids in z direction in m/s")
   ("vz_max", PO::value<Real>(&P::vzmax)->default_value(1), "End corner of velocity grids in z direction in m/s")
   ("vx_length", PO::value<uint>(&P::vxblocks_ini)->default_value(1), "Length of the velocity grid in unrefined blocks in the x direction")
   ("vy_length", PO::value<uint>(&P::vyblocks_ini)->default_value(1), "Length of the velocity grid in unrefined blocks in the y direction")
   ("vz_length", PO::value<uint>(&P::vzblocks_ini)->default_value(1), "Length of the velocity grid in unrefined blocks in the z direction");
  
   // Parameters related to time stepping:
   descriptions->add_options()
   ("q", PO::value<Real>(&P::q)->default_value(1.0), "Charge (C) of simulated particle species")
   ("m", PO::value<Real>(&P::m)->default_value(1.0), "Mass (kg) of simulated particle species")
   ("dt", PO::value<Real>(&P::dt)->default_value(1), "Length of one time step in seconds")
   ("time_steps", PO::value<uint>(&P::tsteps)->default_value(1), "Number of time steps to take");
   
   // Parameters related to solar wind simulations:
   // DEPRECATED: These will be moved to somewhere else in the future.
   descriptions->add_options()
   ("solar_wind_file", PO::value<std::string>(&P::solar_wind_file)->default_value(""), "Read solar wind data from the file arg");
   
   // Parameters related to saving data:
   // WARNING: Some of these parameters may become deprecated in the future.
   descriptions->add_options()
   ("save_interval", PO::value<uint>(&P::diagnInterval)->default_value(1), "Save the simulation every arg time steps")
   ("restart_interval", PO::value<uint>(&P::saveRestartInterval)->default_value(numeric_limits<uint>::max()), "Save the complete simulation every arg time steps")
   ("save_spatial_grid", PO::value<bool>(&P::save_spatial_grid)->default_value(true), "Save spatial cell averages for the whole simulation")
   ("save_velocity_grid", PO::value<bool>(&P::save_velocity_grid)->default_value(false), "Save velocity grid from every spatial cell in the simulation")
   ("save_spatial_cells_at_x,X", PO::value<std::vector<Real> >(&P::save_spatial_cells_x)->composing(), "Save the velocity grid in spatial cells at these coordinates (x components, also give as many y and z components, values from command line, configuration files and environment variables are added together [short version only works on command line])")
   ("save_spatial_cells_at_y,Y", PO::value<std::vector<Real> >(&P::save_spatial_cells_y)->composing(), "Save the velocity grid in spatial cells at these (y components, also give as many x and z components, values from command line, configuration files and environment variables are added together [short version only works on command line])")
   ("save_spatial_cells_at_z,Z", PO::value<std::vector<Real> >(&P::save_spatial_cells_z)->composing(), "Save the velocity grid in spatial cells at these coordinates (z components, also give as many x and y components, values from command line, configuration files and environment variables are added together [short version only works on command line])");
   
   return true;
}

/** Deallocate memory reserved by Parameters.
 * @return If true, class Parameters finalized successfully.
 */
bool Parameters::finalize() {
   delete descriptions;
   delete variables;
   descriptions = NULL;
   variables = NULL;
   initialized = false;
   return true;
}

/** Get the value of the given parameter.
 * @param name The name of the parameter.
 * @param value A variable where the value of the parameter is written.
 * @return If true, the given parameter was found and its value was written to value.
 */
bool Parameters::get(const std::string& name,char& value) {
   if (variables->count(name) > 0) {value = (*variables)[name].as<char>(); return true;}
   return false;
}

/** Get the value of the given parameter.
 * @param name The name of the parameter.
 * @param value A variable where the value of the parameter is written.
 * @return If true, the given parameter was found and its value was written to value.
 */
bool Parameters::get(const std::string& name,int& value) {
   if (variables->count(name) > 0) {value = (*variables)[name].as<int>(); return true;}
   return false;
}

/** Get the value of the given parameter.
 * @param name The name of the parameter.
 * @param value A variable where the value of the parameter is written.
 * @return If true, the given parameter was found and its value was written to value.
 */
bool Parameters::get(const std::string& name,unsigned int& value) {
   if (variables->count(name) > 0) {value = (*variables)[name].as<unsigned int>(); return true;}
   return false;
}

/** Get the value of the given parameter.
 * @param name The name of the parameter.
 * @param value A variable where the value of the parameter is written.
 * @return If true, the given parameter was found and its value was written to value.
 */
bool Parameters::get(const std::string& name,float& value) {
   if (variables->count(name) > 0) {value = (*variables)[name].as<float>(); return true;}
   return false;
}

/** Get the value of the given parameter.
 * @param name The name of the parameter.
 * @param value A variable where the value of the parameter is written.
 * @return If true, the given parameter was found and its value was written to value.
 */
bool Parameters::get(const std::string& name,double& value) {
   if (variables->count(name) > 0) {value = (*variables)[name].as<double>(); return true;}
   return false;
}

/** Get the value of the given parameter.
 * @param name The name of the parameter.
 * @param value A variable where the value of the parameter is written.
 * @return If true, the given parameter was found and its value was written to value.
 */
bool Parameters::get(const std::string& name,std::string& value) {
   if (variables->count(name) > 0) {value = (*variables)[name].as<string>(); return true;}
   return false;
}

/** Write the descriptions of known input options to standard output if 
 * an option called "help" has been read.
 * @return If true, option called "help" was found and descriptions were written 
 * to standard output.
 */
bool Parameters::helpMessage() {
   if (variables->count("help") > 0) {
      cout << *descriptions << endl;
      return true;
   }
   return false;
}

/** Query is Parameters has initialized successfully.
 * @return If true, Parameters is ready for use.
 */
bool Parameters::isInitialized() {return initialized;}

/** Request Parameters to reparse input file(s). This function needs 
 * to be called after new options have been added via Parameters:add functions.
 * Otherwise the values of the new options are not read.
 * @return If true, input file(s) were parsed successfully.
 */
bool Parameters::parse() {
   if (initialized == false) return false;
   // Tell Boost to allow undescribed options (throws exception otherwise)
   const bool ALLOW_UNKNOWN = true;
   
   // Read options from command line:
   PO::store(PO::parse_command_line(argc, argv, *descriptions), *variables);
   PO::notify(*variables);
   // Read options from environment variables:
   PO::store(PO::parse_environment(*descriptions, "MAIN_"), *variables);
   PO::notify(*variables);
   // Read options from run config file:
   if (run_config_file_name.size() > 0) {
      ifstream run_config_file(run_config_file_name.c_str(), fstream::in);
      if (run_config_file.good() == true) {
	 PO::store(PO::parse_config_file(run_config_file, *descriptions, ALLOW_UNKNOWN), *variables);
	 PO::notify(*variables);
      }
      run_config_file.close();
   }
   // Read options from user config file:
   if (user_config_file_name.size() > 0) {
      ifstream user_config_file(user_config_file_name.c_str(), fstream::in);
      if (user_config_file.good() == true) {
	 PO::store(PO::parse_config_file(user_config_file, *descriptions, ALLOW_UNKNOWN), *variables);
	 PO::notify(*variables);
      }
      user_config_file.close();
   }
   // Read options from global config file:
   if (global_config_file_name.size() > 0) {
      ifstream global_config_file(global_config_file_name.c_str(), fstream::in);
      if (global_config_file.good() == true) {
	 PO::store(PO::parse_config_file(global_config_file, *descriptions, ALLOW_UNKNOWN), *variables);
	 PO::notify(*variables);
      }
      global_config_file.close();
   }
   
   // Sanity checks (DEPRECATED):
   if (P::save_spatial_cells_x.size() != P::save_spatial_cells_y.size()
       || P::save_spatial_cells_x.size() != P::save_spatial_cells_z.size()) {
      cerr << "Must have equal number of values for x, y and z components, but given: x " << P::save_spatial_cells_x.size();
      cerr << ", y " << P::save_spatial_cells_y.size() << ", z " << P::save_spatial_cells_z.size() << endl;
      exit(1);
   }
   
   dx_ini = (xmax-xmin)/xcells_ini;
   dy_ini = (ymax-ymin)/ycells_ini;
   dz_ini = (zmax-zmin)/zcells_ini;
   
   q_per_m = q/m;
   return true;
}





