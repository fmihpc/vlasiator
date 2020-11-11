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
#include <vector>
#include <map>
#include "particles.h"
#include "particleparameters.h"
#include "field.h"
#include "histogram.h"

// A "Scenario" describes one possible class of simulation setup, and
// Output quantities we are interested in.
// There are multiple different ways in which particles can be used on
// vlasiator data, this one abstracts these ways.

struct Scenario {

	// Fill in the initial particles, given electromagnetic- and velocity fields
  // Further parameters, depending on the scenario, are given by the parameter object
	virtual ParticleContainer initialParticles(Field& E, Field& B, Field& V) {return ParticleContainer();};

  // Do something when a new input timestep has been opened
  virtual void newTimestep(int input_file_counter, int step, double time, ParticleContainer& particles, Field& E,
        Field& B, Field& V) {};

  // Modify or analyze particle behaviour before they are moved by the particle pusher.
  virtual void beforePush(ParticleContainer& particles, Field& E, Field& B, Field& V) {};

  // Modify or analzye particle behaviour after they are moved by the particle pusher
  virtual void afterPush(int step, double time, ParticleContainer& particles, Field& E, Field& B, Field& V) {};

  // Analyze the final state and / or write output
  virtual void finalize(ParticleContainer& particles, Field& E, Field& B, Field& V) {};

  // Flags specifying which fields are required for this scenario
  bool needV; // Velocity field (it's always available for initialization)

  Scenario() : needV(false) {};
};





// Specific scenarios below here

// Trace a single particle, it's initial position and velocity given in the parameter file
struct singleParticleScenario : Scenario {
  ParticleContainer initialParticles(Field& E, Field& B, Field& V);
  virtual void afterPush(int step, double time, ParticleContainer& particles, Field& E, Field& B, Field& V);

  singleParticleScenario() {needV = false;};
};

// Sample a bunch of particles from a distribution, create them at a given point, then trace them
struct distributionScenario : Scenario {
  ParticleContainer initialParticles(Field& E, Field& B, Field& V);
  void newTimestep(int input_file_counter, int step, double time, ParticleContainer& particles, Field& E, Field& B,
        Field& V);
  void finalize(ParticleContainer& particles, Field& E, Field& B, Field& V);

  distributionScenario() {needV = false;};
};

// Inject particles in the tail continuously, see where they precipitate
struct precipitationScenario : Scenario {
  void newTimestep(int input_file_counter, int step, double time, ParticleContainer& particles, Field& E, Field& B,
        Field& V);
  void afterPush(int step, double time, ParticleContainer& particles, Field& E, Field& B, Field& V);

  precipitationScenario() {needV = true;};
};

// For interactive usage from analysator: read positions and velocities from stdin, push those particles.
struct analysatorScenario : Scenario {
  ParticleContainer initialParticles(Field& E, Field& B, Field& V);
  void newTimestep(int input_file_counter, int step, double time, ParticleContainer& particles, Field& E, Field& B,
        Field& V);
  analysatorScenario() {needV = false;};
};

// Analyzation of shock reflectivity in radial IMF runs
struct shockReflectivityScenario : Scenario {

  LinearHistogram2D transmitted;
  LinearHistogram2D reflected;

  void newTimestep(int input_file_counter, int step, double time, ParticleContainer& particles, Field& E, Field& B,
        Field& V);
  void afterPush(int step, double time, ParticleContainer& particles, Field& E, Field& B, Field& V);
  void finalize(ParticleContainer& particles, Field& E, Field& B, Field& V);

  shockReflectivityScenario() :
      transmitted(200,300, Vec2d(ParticleParameters::reflect_start_y,ParticleParameters::start_time),
            Vec2d(ParticleParameters::reflect_stop_y,ParticleParameters::end_time)),
      reflected(200,300, Vec2d(ParticleParameters::reflect_start_y,ParticleParameters::start_time),
            Vec2d(ParticleParameters::reflect_stop_y,ParticleParameters::end_time)) {
        needV= true;
      }
};

// Initialize particles on a plane in front of the shock, track their precipitation upstream or downstream
struct ipShockScenario : Scenario {
  FILE * traFile;
  FILE * refFile;

  ParticleContainer initialParticles(Field& E, Field& B, Field& V);
  void newTimestep(int input_file_counter, int step, double time, ParticleContainer& particles, Field& E, Field& B,
        Field& V);
  void afterPush(int step, double time, ParticleContainer& particles, Field& E, Field& B, Field& V);
  void finalize(ParticleContainer& particles, Field& E, Field& B, Field& V);

  ipShockScenario() {needV = true;};
};

template<typename T> Scenario* createScenario() {
  return new T;
}

Scenario* createScenario(std::string name);
