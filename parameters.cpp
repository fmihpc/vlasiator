#include "boost/program_options.hpp"
#include <cmath>
#include "fstream"
#include "iostream"
#include <limits>
#include "string"

#include "parameters.h"

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

Real P::dt = NAN;
uint P::tstep = 0;
uint P::tsteps = 0;
uint P::saveInterval = numeric_limits<uint>::max();
uint P::diagnInterval = numeric_limits<uint>::max();

bool P::save_spatial_grid;
bool P::save_velocity_grid;

std::vector<Real> P::save_spatial_cells_x;
std::vector<Real> P::save_spatial_cells_y;
std::vector<Real> P::save_spatial_cells_z;

uint P::transmit = 0;

boost::program_options::options_description P::options = boost::program_options::options_description("Usage: main [options (options given on the command line override options given everywhere else)], where options are:");

Parameters::Parameters(int argc, char* argv[]) {

	string global_config_file_name, user_config_file_name, run_config_file_name;

	options.add_options()
		("help", "print this help message")

		// configuration file options
		("global_config", boost::program_options::value<std::string>(&global_config_file_name)->default_value(""), "read options from the global configuration file arg (relative to the current working directory). Options given in this file are overridden by options given in the user's and run's configuration files and by options given in environment variables (prefixed with MAIN_) and the command line")
		("user_config", boost::program_options::value<std::string>(&user_config_file_name)->default_value(""), "read options from the user's configuration file arg (relative to the current working directory). Options given in this file override options given in the global configuration file. Options given in this file are overridden by options given in the run's configuration file and by options given in environment variables (prefixed with MAIN_) and the command line")
		("run_config", boost::program_options::value<std::string>(&run_config_file_name)->default_value(""), "read options from the run's configuration file arg (relative to the current working directory). Options given in this file override options given in the user's and global configuration files. Options given in this file are overridden by options given in environment variables (prefixed with MAIN_) and the command line")

		// spatial grid options
		("x_min", boost::program_options::value<Real>(&P::xmin)->default_value(0), "Starting corner of spatial grid in x direction in meters")
		("x_max", boost::program_options::value<Real>(&P::xmax)->default_value(1), "End corner of spatial grid in x direction in meters")
		("y_min", boost::program_options::value<Real>(&P::ymin)->default_value(0), "Starting corner of spatial grid in y direction in meters")
		("y_max", boost::program_options::value<Real>(&P::ymax)->default_value(1), "End corner of spatial grid in y direction in meters")
		("z_min", boost::program_options::value<Real>(&P::zmin)->default_value(0), "Starting corner of spatial grid in z direction in meters")
		("z_max", boost::program_options::value<Real>(&P::zmax)->default_value(1), "End corner of spatial grid in z direction in meters")
		("x_length", boost::program_options::value<uint>(&P::xcells_ini)->default_value(1), "Length of the spatial grid in unrefined cells in the x direction")
		("y_length", boost::program_options::value<uint>(&P::ycells_ini)->default_value(1), "Length of the spatial grid in unrefined cells in the y direction")
		("z_length", boost::program_options::value<uint>(&P::zcells_ini)->default_value(1), "Length of the spatial grid in unrefined cells in the z direction")

		// velocity grid options
		("vx_min", boost::program_options::value<Real>(&P::vxmin)->default_value(-1), "Starting corner of velocity grids in x direction in m/s")
		("vx_max", boost::program_options::value<Real>(&P::vxmax)->default_value(1), "End corner of velocity grids in x direction in m/s")
		("vy_min", boost::program_options::value<Real>(&P::vymin)->default_value(-1), "Starting corner of velocity grids in y direction in m/s")
		("vy_max", boost::program_options::value<Real>(&P::vymax)->default_value(1), "End corner of velocity grids in y direction in m/s")
		("vz_min", boost::program_options::value<Real>(&P::vzmin)->default_value(-1), "Starting corner of velocity grids in z direction in m/s")
		("vz_max", boost::program_options::value<Real>(&P::vzmax)->default_value(1), "End corner of velocity grids in z direction in m/s")
		("vx_length", boost::program_options::value<uint>(&P::vxblocks_ini)->default_value(1), "Length of the velocity grid in unrefined blocks in the x direction")
		("vy_length", boost::program_options::value<uint>(&P::vyblocks_ini)->default_value(1), "Length of the velocity grid in unrefined blocks in the y direction")
		("vz_length", boost::program_options::value<uint>(&P::vzblocks_ini)->default_value(1), "Length of the velocity grid in unrefined blocks in the z direction")

		// time stepping parameters
		("dt", boost::program_options::value<Real>(&P::dt)->default_value(1), "Length of one time step in seconds")
		("time_steps", boost::program_options::value<uint>(&P::tsteps)->default_value(1), "Number of time steps to take")

		// saving options
		("save_interval", boost::program_options::value<uint>(&P::diagnInterval)->default_value(1), "Save the simulation every arg time steps")
		("save_spatial_grid", boost::program_options::value<bool>(&P::save_spatial_grid)->default_value(true), "Save spatial cell averages for the whole simulation")
		("save_velocity_grid", boost::program_options::value<bool>(&P::save_velocity_grid)->default_value(false), "Save velocity grid from every spatial cell in the simulation")
		("save_spatial_cells_at_x,X", boost::program_options::value<std::vector<Real> >(&P::save_spatial_cells_x)->composing(), "Save the velocity grid in spatial cells at these coordinates (x components, also give as many y and z components, values from command line, configuration files and environment variables are added together [short version only works on command line])")
		("save_spatial_cells_at_y,Y", boost::program_options::value<std::vector<Real> >(&P::save_spatial_cells_y)->composing(), "Save the velocity grid in spatial cells at these coordinates (y components, also give as many x and z components, values from command line, configuration files and environment variables are added together [short version only works on command line])")
		("save_spatial_cells_at_z,Z", boost::program_options::value<std::vector<Real> >(&P::save_spatial_cells_z)->composing(), "Save the velocity grid in spatial cells at these coordinates (z components, also give as many x and y components, values from command line, configuration files and environment variables are added together [short version only works on command line])");


	boost::program_options::variables_map variables;

	// read options from command line
	boost::program_options::store(boost::program_options::parse_command_line(argc, argv, options), variables);
	boost::program_options::notify(variables);

	// read options from environment
	boost::program_options::store(boost::program_options::parse_environment(options, "MAIN_"), variables);
	boost::program_options::notify(variables);

	// read options from run config file
	if (run_config_file_name.size() > 0) {
		std::ifstream run_config_file(run_config_file_name.c_str(), std::fstream::in);
		if (run_config_file) {
			boost::program_options::store(boost::program_options::parse_config_file(run_config_file, options), variables);
			boost::program_options::notify(variables);
			run_config_file.close();
		}
	}

	// read options from user config file
	if (user_config_file_name.size() > 0) {
		std::ifstream user_config_file(user_config_file_name.c_str(), std::fstream::in);
		if (user_config_file) {
			boost::program_options::store(boost::program_options::parse_config_file(user_config_file, options), variables);
			boost::program_options::notify(variables);
			user_config_file.close();
		}
	}

	// read options from global config file
	if (global_config_file_name.size() > 0) {
		std::ifstream global_config_file(global_config_file_name.c_str(), std::fstream::in);
		if (global_config_file) {
			boost::program_options::store(boost::program_options::parse_config_file(global_config_file, options), variables);
			boost::program_options::notify(variables);
			global_config_file.close();
		}
	}

	if (P::save_spatial_cells_x.size() != P::save_spatial_cells_y.size()
	|| P::save_spatial_cells_x.size() != P::save_spatial_cells_z.size()) {
		cerr << "Must have equal number of values for x, y and z components, but given: x " << P::save_spatial_cells_x.size() << ", y " << P::save_spatial_cells_y.size() << ", z " << P::save_spatial_cells_z.size() << endl;
		exit(1);
	}

	dx_ini = (xmax-xmin)/xcells_ini;
	dy_ini = (ymax-ymin)/ycells_ini;
	dz_ini = (zmax-zmin)/zcells_ini;

	// help the user if requested and exit successfully
	if (variables.count("help")) {
		cout << options << endl;
		exit(0);
	}
}







