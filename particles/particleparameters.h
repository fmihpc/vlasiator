#pragma once

#include <random>
#include <string>

#include "distribution.h"
#include "../definitions.h"

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

   static uint64_t num_particles; /*!< Number of particles to generate */

   static std::default_random_engine::result_type random_seed; /*!< Random seed for particle creation */
   static Distribution* (*distribution)(std::default_random_engine&); /*!< Type of distribution from which to sample the particles */
   static Real temperature; /*!< Initial particle temperature (for distributions where a temperature is meaningful) */
   static Real particle_vel; /*!< Initial particle velocity (for distributions with a single initial velocity) */

   static Real mass; /*!< Mass of the test particles */
   static Real charge; /*!< Charge of the test particles */

   static bool addParameters();
   static bool getParameters();
};
