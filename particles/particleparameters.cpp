#include <map>
#include <iostream>
#include "particleparameters.h"
#include "distribution.h"
#include "../readparameters.h"
#include "physconst.h"

typedef ParticleParameters P;

std::string P::input_filename_pattern;
std::string P::output_filename_pattern;
std::string P::mode = std::string("distribution");
Real P::init_x = 0;
Real P::init_y = 0;
Real P::init_z = 0;

Real P::dt = 0;
Real P::start_time = 0;
Real P::end_time = 0;
uint64_t P::num_particles = 0;

std::default_random_engine::result_type P::random_seed = 1;
Distribution* (*P::distribution)(std::default_random_engine&) = NULL;
Real P::temperature = 1e6;
Real P::particle_vel = 0;
Real P::mass = PhysicalConstantsSI::mp;
Real P::charge = PhysicalConstantsSI::e;

bool ParticleParameters::addParameters() {
   Readparameters::add("particles.input_filename_pattern","Printf() like pattern giving the field input filenames.", std::string("bulk.%07i.vlsv"));
   Readparameters::add("particles.output_filename_pattern","Printf() like pattern giving the particle output filenames.", std::string("particles.%07i.vlsv"));
   Readparameters::add("particles.mode","Mode to run the particle pusher in.",std::string("distribution"));

   Readparameters::add("particles.init_x", "Particle starting point, x-coordinate (meters).", 0);
   Readparameters::add("particles.init_y", "Particle starting point, y-coordinate (meters).", 0);
   Readparameters::add("particles.init_z", "Particle starting point, z-coordinate (meters).", 0);

   Readparameters::add("particles.dt", "Particle pusher timestep",0);
   Readparameters::add("particles.start_time", "Simulation time (seconds) for particle start.",0);
   Readparameters::add("particles.end_time", "Simulation time (seconds) at which particle simulation stops.",0);
   Readparameters::add("particles.num_particles", "Number of particles to simulate.",10000);
   Readparameters::add("particles.random_seed", "Random seed for particle creation.",1);
   Readparameters::add("particles.distribution", "Type of distribution function to sample particles from.",std::string("maxwell"));
   Readparameters::add("particles.temperature", "Temperature of the particle distribution",1e6);
   Readparameters::add("particles.particle_vel", "Initial velocity of the particles (in the plasma rest frame)",0);
   Readparameters::add("particles.mass", "Mass of the test particles",PhysicalConstantsSI::mp);
   Readparameters::add("particles.charge", "Charge of the test particles",PhysicalConstantsSI::e);
   return true;
}

bool ParticleParameters::getParameters() {
   Readparameters::get("particles.input_filename_pattern",P::input_filename_pattern);
   Readparameters::get("particles.output_filename_pattern",P::output_filename_pattern);

   Readparameters::get("particles.mode",P::mode);

   Readparameters::get("particles.init_x",P::init_x);
   Readparameters::get("particles.init_y",P::init_y);
   Readparameters::get("particles.init_z",P::init_z);

   Readparameters::get("particles.dt",P::dt);
   Readparameters::get("particles.start_time",P::start_time);
   Readparameters::get("particles.end_time",P::end_time);
   Readparameters::get("particles.num_particles",P::num_particles);
   if(P::dt == 0 || P::end_time <= P::start_time) {
      return false;
   }

   Readparameters::get("particles.random_seed",P::random_seed);

   /* Look up particle distribution generator */
   std::string distribution_name;
   Readparameters::get("particles.distribution",distribution_name);

   std::map<std::string, Distribution*(*)(std::default_random_engine&)> distribution_lookup;
   distribution_lookup["maxwell"]=&create_distribution<Maxwell_Boltzmann>;
   distribution_lookup["monoenergetic"]=&create_distribution<Monoenergetic>;
   distribution_lookup["kappa2"]=&create_distribution<Kappa2>;
   distribution_lookup["kappa6"]=&create_distribution<Kappa6>;

   if(distribution_lookup.find(distribution_name) == distribution_lookup.end()) {
      std::cerr << "Error: particles.distribution value \"" << distribution_name << "\" does not specify a valid distribution!" << std::endl;
      return false;
   } else {
      P::distribution = distribution_lookup[distribution_name];
   }

   Readparameters::get("particles.temperature",P::temperature);
   Readparameters::get("particles.particle_vel",P::particle_vel);
   return true;
}
