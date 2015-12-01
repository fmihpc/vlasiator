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
Real P::input_dt = 1;
Real P::start_time = 0;
Real P::end_time = 0;
uint64_t P::num_particles = 0;

std::default_random_engine::result_type P::random_seed = 1;
Distribution* (*P::distribution)(std::default_random_engine&) = NULL;
Real P::temperature = 1e6;
Real P::particle_vel = 0;
Real P::mass = PhysicalConstantsSI::mp;
Real P::charge = PhysicalConstantsSI::e;

Boundary* P::boundary_behaviour_x = NULL;
Boundary* P::boundary_behaviour_y = NULL;
Boundary* P::boundary_behaviour_z = NULL;

Real P::precip_inner_boundary;
Real P::precip_start_x;
Real P::precip_stop_x;

Real P::reflect_start_y;
Real P::reflect_stop_y;
Real P::reflect_y_scale;
Real P::reflect_x_offset;
Real P::reflect_upstream_boundary;
Real P::reflect_downstream_boundary;

bool ParticleParameters::addParameters() {
   Readparameters::add("particles.input_filename_pattern","Printf() like pattern giving the field input filenames.", std::string("bulk.%07i.vlsv"));
   Readparameters::add("particles.output_filename_pattern","Printf() like pattern giving the particle output filenames.", std::string("particles.%07i.vlsv"));
   Readparameters::add("particles.mode","Mode to run the particle pusher in.",std::string("distribution"));

   Readparameters::add("particles.init_x", "Particle starting point, x-coordinate (meters).", 0);
   Readparameters::add("particles.init_y", "Particle starting point, y-coordinate (meters).", 0);
   Readparameters::add("particles.init_z", "Particle starting point, z-coordinate (meters).", 0);

   Readparameters::add("particles.dt", "Particle pusher timestep",0);
   Readparameters::add("particles.input_dt", "Time spacing (seconds) of input files",1.);
   Readparameters::add("particles.start_time", "Simulation time (seconds) for particle start.",0);
   Readparameters::add("particles.end_time", "Simulation time (seconds) at which particle simulation stops.",0);
   Readparameters::add("particles.num_particles", "Number of particles to simulate.",10000);
   Readparameters::add("particles.random_seed", "Random seed for particle creation.",1);
   Readparameters::add("particles.distribution", "Type of distribution function to sample particles from.",std::string("maxwell"));
   Readparameters::add("particles.temperature", "Temperature of the particle distribution",1e6);
   Readparameters::add("particles.particle_vel", "Initial velocity of the particles (in the plasma rest frame)",0);
   Readparameters::add("particles.mass", "Mass of the test particles",PhysicalConstantsSI::mp);
   Readparameters::add("particles.charge", "Charge of the test particles",PhysicalConstantsSI::e);

   Readparameters::add("particles.boundary_behaviour_x", "What to do with particles that reach the x boundaries (DELETE/REFLECT/PERIODIC)",std::string("DELETE"));
   Readparameters::add("particles.boundary_behaviour_y", "What to do with particles that reach the y boundaries (DELETE/REFLECT/PERIODIC)",std::string("PERIODIC"));
   Readparameters::add("particles.boundary_behaviour_z", "What to do with particles that reach the z boundaries (DELETE/REFLECT/PERIODIC)",std::string("PERIODIC"));

   // Parameters for the precipitation mode
   Readparameters::add("particles.inner_boundary", "Distance of the inner boundary from the coordinate centre (meters)", 30e6);
   Readparameters::add("particles.precipitation_start_x", "X-Coordinate at which precipitation injection starts (meters)", -200e6);
   Readparameters::add("particles.precipitation_stop_x", "X-Coordinate at which precipitation injection stops (meters)", -50e6);

   // Parameters for shock reflection mode
   Readparameters::add("particles.reflect_start_y", "Y-Coordinate of the bottom end of the parabola, at which shock reflection scenario particles are injected", 60e6);
   Readparameters::add("particles.reflect_stop_y", "Y-Coordinate of the bottom end of the parabola, at which shock reflection scenario particles are injected", 60e6);
   Readparameters::add("particles.reflect_y_scale", "Curvature scale of the injection parabola for the reflection scenario", 60e6);
   Readparameters::add("particles.reflect_x_offset", "X-Coordinate of the tip of the injection parabola for the reflection scenario", 40e6);
   Readparameters::add("particles.reflect_upstream_boundary", "Distance from particle injection point at which particles are to be counted as 'reflected'", 10e6);
   Readparameters::add("particles.reflect_downstream_boundary", "Distance from particle injection point at which particles are to be counted as 'transmitted'", 10e6);

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
   Readparameters::get("particles.input_dt", P::input_dt);
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

   // Boundaries
   std::map<std::string, Boundary*(*)(int)> boundary_lookup;
   boundary_lookup["DELETE"] = &create_boundary<OpenBoundary>;
   boundary_lookup["REFLECT"] = &create_boundary<ReflectBoundary>;
   boundary_lookup["PERIODIC"] = &create_boundary<PeriodicBoundary>;
   std::string tempstring;
   Readparameters::get("particles.boundary_behaviour_x",tempstring);
   if(boundary_lookup.find(tempstring) == boundary_lookup.end()) {
      std::cerr << "Error: invalid boundary condition \"" << tempstring << "\" in x-direction" << std::endl;
      exit(0);
   } else {
      P::boundary_behaviour_x = boundary_lookup[tempstring](0);
   }
   Readparameters::get("particles.boundary_behaviour_y",tempstring);
   if(boundary_lookup.find(tempstring) == boundary_lookup.end()) {
      std::cerr << "Error: invalid boundary condition \"" << tempstring << "\" in y-direction" << std::endl;
      exit(0);
   } else {
      P::boundary_behaviour_y = boundary_lookup[tempstring](1);
   }
   Readparameters::get("particles.boundary_behaviour_z",tempstring);
   if(boundary_lookup.find(tempstring) == boundary_lookup.end()) {
      std::cerr << "Error: invalid boundary condition \"" << tempstring << "\" in z-direction" << std::endl;
      exit(0);
   } else {
      P::boundary_behaviour_z = boundary_lookup[tempstring](2);
   }

   Readparameters::get("particles.inner_boundary", P::precip_inner_boundary);
   Readparameters::get("particles.precipitation_start_x", P::precip_start_x);
   Readparameters::get("particles.precipitation_stop_x", P::precip_stop_x);

   Readparameters::get("particles.reflect_start_y", P::reflect_start_y);
   Readparameters::get("particles.reflect_stop_y", P::reflect_stop_y);
   Readparameters::get("particles.reflect_y_scale", P::reflect_y_scale);
   Readparameters::get("particles.reflect_x_offset", P::reflect_x_offset);
   Readparameters::get("particles.reflect_upstream_boundary", P::reflect_upstream_boundary);
   Readparameters::get("particles.reflect_downstream_boundary", P::reflect_downstream_boundary);

   return true;
}
