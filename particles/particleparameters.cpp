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
std::string P::V_field_name = "V";
std::string P::rho_field_name = "rho";
bool P::divide_rhov_by_rho = false;

std::default_random_engine::result_type P::random_seed = 1;
Distribution* (*P::distribution)(std::default_random_engine&) = NULL;
std::string distribution_name;
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

Real P::ipshock_inject_x0;
Real P::ipshock_inject_x1;
Real P::ipshock_inject_y0;
Real P::ipshock_inject_y1;
Real P::ipshock_inject_z0;
Real P::ipshock_inject_z1;
Real P::ipshock_transmit;
Real P::ipshock_reflect;

bool ParticleParameters::addParameters() {
   Readparameters::add<std::string>("particles.input_filename_pattern","Printf() like pattern giving the field input filenames.", P::input_filename_pattern,
         std::string("bulk.%07i.vlsv"));
   Readparameters::add<std::string>("particles.output_filename_pattern","Printf() like pattern giving the particle output filenames.", P::output_filename_pattern,
         std::string("particles.%07i.vlsv"));
   Readparameters::add<std::string>("particles.mode","Mode to run the particle pusher in.", P::mode,std::string("distribution"));

   Readparameters::add<Real>("particles.init_x", "Particle starting point, x-coordinate (meters).", P::init_x, 0);
   Readparameters::add<Real>("particles.init_y", "Particle starting point, y-coordinate (meters).", P::init_y, 0);
   Readparameters::add<Real>("particles.init_z", "Particle starting point, z-coordinate (meters).", P::init_z, 0);

   Readparameters::add<Real>("particles.dt", "Particle pusher timestep", P::dt,0);
   Readparameters::add<Real>("particles.input_dt", "Time spacing (seconds) of input files", P::input_dt,1.);
   Readparameters::add<Real>("particles.start_time", "Simulation time (seconds) for particle start.", P::start_time,0);
   Readparameters::add<Real>("particles.end_time", "Simulation time (seconds) at which particle simulation stops.", P::end_time,0);
   Readparameters::add<uint64_t>("particles.num_particles", "Number of particles to simulate.", P::num_particles,10000);
   Readparameters::add<std::string>("particles.V_field_name", "Name of the Velocity data set in the input files", P::V_field_name, "V");
   Readparameters::add<std::string>("particles.rho_field_name", "Name of the Density data set in the input files", P::rho_field_name, "rho");
   Readparameters::add<bool>("particles.divide_rhov_by_rho", "Do the input file store rho_v and rho separately?", P::divide_rhov_by_rho, false);
   Readparameters::add<int>("particles.random_seed", "Random seed for particle creation.", P::random_seed,1);
   Readparameters::add<std::string>("particles.distribution", "Type of distribution function to sample particles from.", P::distribution_name,
         std::string("maxwell"));
   Readparameters::add<Real>("particles.temperature", "Temperature of the particle distribution", P::temperature,1e6);
   Readparameters::add<Real>("particles.particle_vel", "Initial velocity of the particles (in the plasma rest frame)", P::particle_vel,0);
   Readparameters::add<Real>("particles.mass", "Mass of the test particles", P::mass,PhysicalConstantsSI::mp);
   Readparameters::add<Real>("particles.charge", "Charge of the test particles", P::charge,PhysicalConstantsSI::e);

   Readparameters::add<std::string>("particles.boundary_behaviour_x",
         "What to do with particles that reach the x boundaries (DELETE/REFLECT/PERIODIC)",P::boundary_behaviour_x_string,std::string("DELETE"));
   Readparameters::add<std::string>("particles.boundary_behaviour_y",
         "What to do with particles that reach the y boundaries (DELETE/REFLECT/PERIODIC)",P::boundary_behaviour_y,std::string("PERIODIC"));
   Readparameters::add<std::string>("particles.boundary_behaviour_z",
         "What to do with particles that reach the z boundaries (DELETE/REFLECT/PERIODIC)",P::boundary_behaviour_z,std::string("PERIODIC"));

   // Parameters for the precipitation mode
   Readparameters::add<Real>("particles.inner_boundary", "Distance of the inner boundary from the coordinate centre (meters)", P::inner_boundary,
         30e6);
   Readparameters::add<Real>("particles.precipitation_start_x", "X-Coordinate at which precipitation injection starts (meters)", P::precipitation_start_x,
         -200e6);
   Readparameters::add<Real>("particles.precipitation_stop_x", "X-Coordinate at which precipitation injection stops (meters)", P::precipitation_stop_x,
         -50e6);

   // Parameters for shock reflection mode
   Readparameters::add<Real>("particles.reflect_start_y",
         "Y-Coordinate of the bottom end of the parabola, at which shock reflection scenario particles are injected",P::reflect_start_y, 60e6);
   Readparameters::add<Real>("particles.reflect_stop_y",
         "Y-Coordinate of the bottom end of the parabola, at which shock reflection scenario particles are injected",P::reflect_stop_y, 60e6);
   Readparameters::add<Real>("particles.reflect_y_scale",
         "Curvature scale of the injection parabola for the reflection scenario",P::reflect_y_scale, 60e6);
   Readparameters::add<Real>("particles.reflect_x_offset",
         "X-Coordinate of the tip of the injection parabola for the reflection scenario",P::reflect_x_offset, 40e6);
   Readparameters::add<Real>("particles.reflect_upstream_boundary",
         "Distance from particle injection point at which particles are to be counted as 'reflected'",P::reflect_upstream_boundary, 10e6);
   Readparameters::add<Real>("particles.reflect_downstream_boundary",
         "Distance from particle injection point at which particles are to be counted as 'transmitted'",P::reflect_downstream_boundary, 10e6);

   // Parameters for ip shock injection test mode
   Readparameters::add<Real>("particles.ipshock_inject_x0",
         "X-Coordinate of the lower edge of particle injection region for the ipShock scenario",P::ipshock_inject_x0, -1.e6);
   Readparameters::add<Real>("particles.ipshock_inject_x1",
         "X-Coordinate of the upper edge of particle injection region for the ipShock scenario",P::ipshock_inject_x1, 1.e6);
   Readparameters::add<Real>("particles.ipshock_inject_y0",
         "Y-Coordinate of the lower edge of particle injection region for the ipShock scenario",P::ipshock_inject_y0, -1.e6);
   Readparameters::add<Real>("particles.ipshock_inject_y1",
         "Y-Coordinate of the upper edge of particle injection region for the ipShock scenario",P::ipshock_inject_y1, 1.e6);
   Readparameters::add<Real>("particles.ipshock_inject_z0",
         "Z-Coordinate of the lower edge of particle injection region for the ipShock scenario",P::ipshock_inject_z0, -1.e6);
   Readparameters::add<Real>("particles.ipshock_inject_z1",
         "Z-Coordinate of the upper edge of particle injection region for the ipShock scenario",P::ipshock_inject_z1, 1.e6);
   Readparameters::add<Real>("particles.ipshock_transmit",
         "X-Coordinate of threshold for where particles are counted  as transmitted for the ipShock scenario",P::ipshock_transmit, -10.e6);
   Readparameters::add<Real>("particles.ipshock_reflect",
         "X-Coordinate of threshold for where particles are counted  as reflected for the ipShock scenario",P::ipshock_reflect, 10.e6);

   return true;
}

bool ParticleParameters::getParameters() {

   if(P::dt == 0 || P::end_time == P::start_time) {
      std::cerr << "Error end_time == start_time! Won't do anything (and will probably crash now)." << std::endl;
      return false;
   }


   /* Look up particle distribution generator */

   std::map<std::string, Distribution*(*)(std::default_random_engine&)> distribution_lookup;
   distribution_lookup["maxwell"]=&createDistribution<Maxwell_Boltzmann>;
   distribution_lookup["monoenergetic"]=&createDistribution<Monoenergetic>;
   distribution_lookup["kappa2"]=&createDistribution<Kappa2>;
   distribution_lookup["kappa6"]=&createDistribution<Kappa6>;

   if(distribution_lookup.find(distribution_name) == distribution_lookup.end()) {
      std::cerr << "Error: particles.distribution value \"" << distribution_name
         << "\" does not specify a valid distribution!" << std::endl;
      return false;
   } else {
      P::distribution = distribution_lookup[distribution_name];
   }

   // Boundaries
   std::map<std::string, Boundary*(*)(int)> boundaryLookup;
   boundaryLookup["DELETE"] = &createBoundary<OpenBoundary>;
   boundaryLookup["REFLECT"] = &createBoundary<ReflectBoundary>;
   boundaryLookup["PERIODIC"] = &createBoundary<PeriodicBoundary>;

   if(boundaryLookup.find(boundary_behaviour_x_string) == boundaryLookup.end()) {
      std::cerr << "Error: invalid boundary condition \"" << boundary_behaviour_x_string << "\" in x-direction" << std::endl;
      exit(0);
   } else {
      P::boundary_behaviour_x = boundaryLookup[boundary_behaviour_x_string](0);
   }

   if(boundaryLookup.find(boundary_behaviour_y_string) == boundaryLookup.end()) {
      std::cerr << "Error: invalid boundary condition \"" << boundary_behaviour_y_string << "\" in y-direction" << std::endl;
      exit(0);
   } else {
      P::boundary_behaviour_y = boundaryLookup[boundary_behaviour_y_string](1);
   }

   if(boundaryLookup.find(boundary_behaviour_z_string) == boundaryLookup.end()) {
      std::cerr << "Error: invalid boundary condition \"" << boundary_behaviour_z_string << "\" in z-direction" << std::endl;
      exit(0);
   } else {
      P::boundary_behaviour_z = boundaryLookup[boundary_behaviour_z_string](2);
   }

   return true;
}
