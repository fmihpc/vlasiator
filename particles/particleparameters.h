#pragma once

#include <random>
#include <string>

#include "distribution.h"
#include "../definitions.h"
#include "boundaries.h"

class Distribution; //forward declaration

struct ParticleParameters {
   static std::string input_filename_pattern; /*!< printf() - like pattern giving the input filenames */
   static std::string output_filename_pattern; /*!< printf() - like pattern giving the output filenames */
   static std::string mode; /*!< Particle tracing mode  */

   static Real init_x; /*!< Particle starting point, x-coordinate */
   static Real init_y; /*!< Particle starting point, y-coordinate */
   static Real init_z; /*!< Particle starting point, z-coordinate */

   static Real dt; /*!< Particle pusher timestep */
   static Real start_time; /*!< Simulation time at which the particles are injected */
   static Real end_time;  /*!< Simulation time at which the particle-simulation should be stopped */
   static Real input_dt; /*!< Time interval between input files */

   static uint64_t num_particles; /*!< Number of particles to generate */

   static Boundary* boundary_behaviour_x; /*!< What to do with particles that reach the x boundary */
   static Boundary* boundary_behaviour_y; /*!< What to do with particles that reach the y boundary */
   static Boundary* boundary_behaviour_z; /*!< What to do with particles that reach the z boundary */

   static Real precip_inner_boundary; /*!< Distance of the inner boundary from the coordinate centre (meters) */
   static Real precip_start_x; /*!< X-Coordinate at which precipitation injection starts (meters) */
   static Real precip_stop_x; /*!< X-Coordinate at which precipitation injection stops (meters) */

   static Real reflect_start_y; /*!< Y-Coordinate of the bottom end of the parabola, at which shock reflection scenario particles are injected */
   static Real reflect_stop_y; /*!< Y-Coordinate of the top end of the parabola, at which shock reflection scenario particles are injected */
   static Real reflect_y_scale; /*!< Curvature scale of the injection parabola for the reflection scenario */
   static Real reflect_x_offset; /*!< X-Coordinate of the tip of the injection parabola for the reflection scenario */
   static Real reflect_upstream_boundary; /*!< Distance from particle injection point at which particles are to be counted as 'reflected' */
   static Real reflect_downstream_boundary; /*!< Distance from particle injection point at which particles are to be counted as 'transmitted' */

   static std::default_random_engine::result_type random_seed; /*!< Random seed for particle creation */
   static Distribution* (*distribution)(std::default_random_engine&); /*!< Type of distribution from which to sample the particles */
   static Real temperature; /*!< Initial particle temperature (for distributions where a temperature is meaningful) */
   static Real particle_vel; /*!< Initial particle velocity (for distributions with a single initial velocity) */

   static Real mass; /*!< Mass of the test particles */
   static Real charge; /*!< Charge of the test particles */

   static bool addParameters();
   static bool getParameters();
};
