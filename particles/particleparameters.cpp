#include "particleparameters.h"
#include "../readparameters.h"

typedef ParticleParameters P;

std::string P::input_filename_pattern;
std::string P::output_filename_pattern;
std::string P::mode = std::string("distribution");
Real P::init_x = 0;
Real P::init_y = 0;
Real P::init_z = 0;

Real P::start_time = 0;
Real P::end_time = 0;
uint64_t P::num_particles = 0;

std::string P::distribution = std::string("maxwell");
Real P::temperature = 1e6;

bool ParticleParameters::addParameters() {
	Readparameters::add("particles.input_filename_pattern","Printf() like pattern giving the field input filenames.", std::string("bulk.%07i.vlsv"));
	Readparameters::add("particles.output_filename_pattern","Printf() like pattern giving the particle output filenames.", std::string("particles.%07i.vlsv"));
	Readparameters::add("particles.mode","Mode to run the particle pusher in.",std::string("distribution"));

	Readparameters::add("particles.init_x", "Particle starting point, x-coordinate (meters).", 0);
	Readparameters::add("particles.init_y", "Particle starting point, y-coordinate (meters).", 0);
	Readparameters::add("particles.init_z", "Particle starting point, z-coordinate (meters).", 0);

	Readparameters::add("particles.start_time", "Simulation time (seconds) for particle start.",0);
	Readparameters::add("particles.end_time", "Simulation time (seconds) at which particle simulation stops.",0);
	Readparameters::add("particles.num_particles", "Number of particles to simulate.",10000);
	Readparameters::add("particles.distribution", "Type of distribution function to sample particles from.",std::string("maxwell"));
	Readparameters::add("particles.temperature", "Temperature of the particle distribution",1e6);
	return true;
}

bool ParticleParameters::getParameters() {
	Readparameters::get("particles.input_filename_pattern",P::input_filename_pattern);
	Readparameters::get("particles.output_filename_pattern",P::output_filename_pattern);

	Readparameters::get("particles.mode",P::mode);

	Readparameters::get("particles.init_x",P::init_x);
	Readparameters::get("particles.init_y",P::init_y);
	Readparameters::get("particles.init_z",P::init_z);

	Readparameters::get("particles.start_time",P::start_time);
	Readparameters::get("particles.end_time",P::end_time);
	Readparameters::get("particles.num_particles",P::num_particles);

	Readparameters::get("particles.distribution",P::distribution);
	Readparameters::get("particles.temperature",P::temperature);
	return true;
}
