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
   Readparameters::add("particles.input_filename_pattern","Printf() like pattern giving the field input filenames.",
         std::string("bulk.%07i.vlsv"));
   Readparameters::add("particles.output_filename_pattern","Printf() like pattern giving the particle output filenames.",
         std::string("particles.%07i.vlsv"));
   Readparameters::add("particles.mode","Mode to run the particle pusher in.",std::string("distribution"));

   Readparameters::add("particles.init_x", "Particle starting point, x-coordinate (meters).", 0);
   Readparameters::add("particles.init_y", "Particle starting point, y-coordinate (meters).", 0);
   Readparameters::add("particles.init_z", "Particle starting point, z-coordinate (meters).", 0);

   Readparameters::add("particles.dt", "Particle pusher timestep",0);
   Readparameters::add("particles.input_dt", "Time spacing (seconds) of input files",1.);
   Readparameters::add("particles.start_time", "Simulation time (seconds) for particle start.",0);
   Readparameters::add("particles.end_time", "Simulation time (seconds) at which particle simulation stops.",0);
   Readparameters::add("particles.num_particles", "Number of particles to simulate.",10000);
   Readparameters::add("particles.V_field_name", "Name of the Velocity data set in the input files", "V");
   Readparameters::add("particles.rho_field_name", "Name of the Density data set in the input files", "rho");
   Readparameters::add("particles.divide_rhov_by_rho", "Do the input file store rho_v and rho separately?", false);
   Readparameters::add("particles.random_seed", "Random seed for particle creation.",1);
   Readparameters::add("particles.distribution", "Type of distribution function to sample particles from.",
         std::string("maxwell"));
   Readparameters::add("particles.temperature", "Temperature of the particle distribution",1e6);
   Readparameters::add("particles.particle_vel", "Initial velocity of the particles (in the plasma rest frame)",0);
   Readparameters::add("particles.mass", "Mass of the test particles",PhysicalConstantsSI::mp);
   Readparameters::add("particles.charge", "Charge of the test particles",PhysicalConstantsSI::e);

   Readparameters::add("particles.boundary_behaviour_x",
         "What to do with particles that reach the x boundaries (DELETE/REFLECT/PERIODIC)",std::string("DELETE"));
   Readparameters::add("particles.boundary_behaviour_y",
         "What to do with particles that reach the y boundaries (DELETE/REFLECT/PERIODIC)",std::string("PERIODIC"));
   Readparameters::add("particles.boundary_behaviour_z",
         "What to do with particles that reach the z boundaries (DELETE/REFLECT/PERIODIC)",std::string("PERIODIC"));

   // Parameters for the precipitation mode
   Readparameters::add("particles.inner_boundary", "Distance of the inner boundary from the coordinate centre (meters)",
         30e6);
   Readparameters::add("particles.precipitation_start_x", "X-Coordinate at which precipitation injection starts (meters)",
         -200e6);
   Readparameters::add("particles.precipitation_stop_x", "X-Coordinate at which precipitation injection stops (meters)",
         -50e6);

   // Parameters for shock reflection mode
   Readparameters::add("particles.reflect_start_y",
         "Y-Coordinate of the bottom end of the parabola, at which shock reflection scenario particles are injected", 60e6);
   Readparameters::add("particles.reflect_stop_y",
         "Y-Coordinate of the bottom end of the parabola, at which shock reflection scenario particles are injected", 60e6);
   Readparameters::add("particles.reflect_y_scale",
         "Curvature scale of the injection parabola for the reflection scenario", 60e6);
   Readparameters::add("particles.reflect_x_offset",
         "X-Coordinate of the tip of the injection parabola for the reflection scenario", 40e6);
   Readparameters::add("particles.reflect_upstream_boundary",
         "Distance from particle injection point at which particles are to be counted as 'reflected'", 10e6);
   Readparameters::add("particles.reflect_downstream_boundary",
         "Distance from particle injection point at which particles are to be counted as 'transmitted'", 10e6);

   // Parameters for ip shock injection test mode
   Readparameters::add("particles.ipshock_inject_x0",
         "X-Coordinate of the lower edge of particle injection region for the ipShock scenario", -1.e6);
   Readparameters::add("particles.ipshock_inject_x1",
         "X-Coordinate of the upper edge of particle injection region for the ipShock scenario", 1.e6);
   Readparameters::add("particles.ipshock_inject_y0",
         "Y-Coordinate of the lower edge of particle injection region for the ipShock scenario", -1.e6);
   Readparameters::add("particles.ipshock_inject_y1",
         "Y-Coordinate of the upper edge of particle injection region for the ipShock scenario", 1.e6);
   Readparameters::add("particles.ipshock_inject_z0",
         "Z-Coordinate of the lower edge of particle injection region for the ipShock scenario", -1.e6);
   Readparameters::add("particles.ipshock_inject_z1",
         "Z-Coordinate of the upper edge of particle injection region for the ipShock scenario", 1.e6);
   Readparameters::add("particles.ipshock_transmit",
         "X-Coordinate of threshold for where particles are counted  as transmitted for the ipShock scenario", -10.e6);
   Readparameters::add("particles.ipshock_reflect",
         "X-Coordinate of threshold for where particles are counted  as reflected for the ipShock scenario", 10.e6);

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
   if(P::dt == 0 || P::end_time == P::start_time) {
      std::cerr << "Error end_time == start_time! Won't do anything (and will probably crash now)." << std::endl;
      return false;
   }
   Readparameters::get("particles.V_field_name",P::V_field_name);
   Readparameters::get("particles.rho_field_name",P::rho_field_name);
   Readparameters::get("particles.divide_rhov_by_rho",P::divide_rhov_by_rho);

   Readparameters::get("particles.random_seed",P::random_seed);

   /* Look up particle distribution generator */
   std::string distribution_name;
   Readparameters::get("particles.distribution",distribution_name);

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

   Readparameters::get("particles.temperature",P::temperature);
   Readparameters::get("particles.particle_vel",P::particle_vel);

   // Boundaries
   std::map<std::string, Boundary*(*)(int)> boundaryLookup;
   boundaryLookup["DELETE"] = &createBoundary<OpenBoundary>;
   boundaryLookup["REFLECT"] = &createBoundary<ReflectBoundary>;
   boundaryLookup["PERIODIC"] = &createBoundary<PeriodicBoundary>;
   std::string tempstring;
   Readparameters::get("particles.boundary_behaviour_x",tempstring);
   if(boundaryLookup.find(tempstring) == boundaryLookup.end()) {
      std::cerr << "Error: invalid boundary condition \"" << tempstring << "\" in x-direction" << std::endl;
      exit(0);
   } else {
      P::boundary_behaviour_x = boundaryLookup[tempstring](0);
   }
   Readparameters::get("particles.boundary_behaviour_y",tempstring);
   if(boundaryLookup.find(tempstring) == boundaryLookup.end()) {
      std::cerr << "Error: invalid boundary condition \"" << tempstring << "\" in y-direction" << std::endl;
      exit(0);
   } else {
      P::boundary_behaviour_y = boundaryLookup[tempstring](1);
   }
   Readparameters::get("particles.boundary_behaviour_z",tempstring);
   if(boundaryLookup.find(tempstring) == boundaryLookup.end()) {
      std::cerr << "Error: invalid boundary condition \"" << tempstring << "\" in z-direction" << std::endl;
      exit(0);
   } else {
      P::boundary_behaviour_z = boundaryLookup[tempstring](2);
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

   Readparameters::get("particles.ipshock_inject_x0", P::ipshock_inject_x0);
   Readparameters::get("particles.ipshock_inject_x1", P::ipshock_inject_x1);
   Readparameters::get("particles.ipshock_inject_y0", P::ipshock_inject_y0);
   Readparameters::get("particles.ipshock_inject_y1", P::ipshock_inject_y1);
   Readparameters::get("particles.ipshock_inject_z0", P::ipshock_inject_z0);
   Readparameters::get("particles.ipshock_inject_z1", P::ipshock_inject_z1);
   Readparameters::get("particles.ipshock_transmit", P::ipshock_transmit);
   Readparameters::get("particles.ipshock_reflect", P::ipshock_reflect);

   return true;
}
