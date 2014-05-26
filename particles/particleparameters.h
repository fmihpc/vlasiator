#pragma once

#include <string>

#include "../definitions.h"

struct ParticleParameters {
	static std::string input_filename_pattern; /*!< printf() - like pattern giving the input filenames */
	static std::string output_filename_pattern; /*!< printf() - like pattern giving the output filenames */
	static std::string mode; /*!< Particle tracing mode  */

	static Real init_x; /*!< Particle starting point, x-coordinate */
	static Real init_y; /*!< Particle starting point, y-coordinate */
	static Real init_z; /*!< Particle starting point, z-coordinate */

	static Real start_time; /*!< Simulation time at which the particles are injected */
	static Real end_time;  /*!< Simulation time at which the particle-simulation should be stopped */

	static uint64_t num_particles; /*!< Number of particles to generate */

	static std::string distribution; /*!< Type of distribution from which to sample the particles */
	static Real temperature; /*!< Initial particle temperature */

	static bool addParameters();
	static bool getParameters();
};
