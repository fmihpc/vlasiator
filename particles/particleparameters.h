#pragma once
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
   static bool staticfields; /*!< Flag for using static fields from first time step */

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

   static Real ipshock_inject_x0; /*!< Smaller x-coordinate for particle injection (meters) */
   static Real ipshock_inject_x1; /*!< Larger x-coordinate for particle injection (meters) */
   static Real ipshock_inject_y0; /*!< Smaller y-coordinate for particle injection (meters) */
   static Real ipshock_inject_y1; /*!< Larger y-coordinate for particle injection (meters) */
   static Real ipshock_inject_z0; /*!< Smaller z-coordinate for particle injection (meters) */
   static Real ipshock_inject_z1; /*!< Larger z-coordinate for particle injection (meters) */
   static Real ipshock_transmit;  /*!< X-Coordinate for particle transmission (downstream) (meters) */
   static Real ipshock_reflect;   /*!< X-Coordinate for particle reflection (upstream) (meters) */

/*    static Real injection_r0;           /\*!< injection scenario: Bow shock fit standoff distance in metres *\/ */
/*    static Real injection_alpha;        /\*!< injection scenario: Bow shock fit tail flare parameter alpha  *\/ */
/*    static Real injection_ecc;          /\*!< injection scenario: Bow shock fit eccentricity parameter ecc*\/ */
   static Real injection_bs_p0;        /*!< injection scenario: Bow shock fit p0 (at start) */
   static Real injection_bs_p1;        /*!< injection scenario: Bow shock fit p1 (at start) */
   static Real injection_bs_p2;        /*!< injection scenario: Bow shock fit p2 (at start) */
   static Real injection_bs_p3;        /*!< injection scenario: Bow shock fit p3 (at start) */
   static Real injection_bs_p4;        /*!< injection scenario: Bow shock fit p4 (at start) */
   static Real injection_bs2_p0;        /*!< injection scenario: Bow shock fit p0 (at end) */
   static Real injection_bs2_p1;        /*!< injection scenario: Bow shock fit p1 (at end) */
   static Real injection_bs2_p2;        /*!< injection scenario: Bow shock fit p2 (at end) */
   static Real injection_bs2_p3;        /*!< injection scenario: Bow shock fit p3 (at end) */
   static Real injection_bs2_p4;        /*!< injection scenario: Bow shock fit p4 (at end) */

   static Real injection_rho_meet;     /*!< injection scenario: Number density value to assume particle meets shock */
   static Real injection_TNBS_meet;    /*!< injection scenario: Non-backstreaming temperature value to assume particle meets shock */
   static Real injection_Mms_meet;     /*!< injection scenario: Magnetosonic Mach number to assume particle meets shock */
   static Real injection_tbn_meet;     /*!< injection scenario: Shock-normal angle value to assume particle meets shock */
   // No parameter required for flipmu calculation

   static Real injection_r_bound_ds;   /*!< injection scenario: Downstream transmission boundary radial distance in metres */
   static Real injection_r_bound_us;   /*!< injection scenario: Upstream reflection boundary radial distance in metres */
   static Real injection_x_bound_ds;   /*!< injection scenario: discard boundary X-coordinate in metres */
   static Real injection_start_deg0;   /*!< injection scenario: initialisation arc start angle in degrees */
   static Real injection_start_deg1;   /*!< injection scenario: initialisation arc finish angle in degrees */
   static Real injection_start_rplus;  /*!< injection scenario: initialisation arc distance from shock in metres */
   static Real injection_init_vx;        /*!< injection scenario: V_x for initialisation solar wind */
   static Real injection_init_vy;        /*!< injection scenario: V_y for initialisation solar wind */
   static Real injection_init_vz;        /*!< injection scenario: V_z for initialisation solar wind */

   static std::default_random_engine::result_type random_seed; /*!< Random seed for particle creation */
   static Distribution* (*distribution)(std::default_random_engine&); /*!< Type of distribution from which to sample the particles */
   static Real temperature; /*!< Initial particle temperature (for distributions where a temperature is meaningful) */
   static Real particle_vel; /*!< Initial particle velocity (for distributions with a single initial velocity) */

   static Real mass; /*!< Mass of the test particles */
   static Real charge; /*!< Charge of the test particles */

   static bool addParameters();
   static bool getParameters();
};
