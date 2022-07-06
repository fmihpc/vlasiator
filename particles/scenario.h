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
	virtual ParticleContainer initialParticles(Field& E, Field& B, Field& V, Field& R) {return ParticleContainer();};

  // Do something when a new input timestep has been opened
  virtual void newTimestep(int input_file_counter, int step, double time, ParticleContainer& particles, Interpolated_Field& E,
        Interpolated_Field& B, Interpolated_Field& V, Interpolated_Field& R) {};

  // Modify or analyze particle behaviour before they are moved by the particle pusher.
  virtual void beforePush(ParticleContainer& particles, Interpolated_Field& E, Interpolated_Field& B, Interpolated_Field& V, Interpolated_Field& R) {};

  // Modify or analzye particle behaviour after they are moved by the particle pusher
  virtual void afterPush(int step, double time, ParticleContainer& particles, Interpolated_Field& E, Interpolated_Field& B, Interpolated_Field& V, Interpolated_Field& R) {};

  // Analyze the final state and / or write output
  virtual void finalize(ParticleContainer& particles, Field& E, Field& B, Field& V, Field& R) {};

  // Flags specifying which fields are required for this scenario
  bool needV; // Velocity field (it's always available for initialization)
  bool needRho; // mass density, for e.g. finding a shock position

Scenario() : needV(false), needRho(false) {};
};





// Specific scenarios below here

// Trace a single particle, it's initial position and velocity given in the parameter file
struct singleParticleScenario : Scenario {
  ParticleContainer initialParticles(Field& E, Field& B, Field& V, Field& R);
  virtual void afterPush(int step, double time, ParticleContainer& particles, Interpolated_Field& E, Interpolated_Field& B, Interpolated_Field& V, Interpolated_Field& R);

  singleParticleScenario() {needV = false; needRho = false;};
};

// Sample a bunch of particles from a distribution, create them at a given point, then trace them
struct distributionScenario : Scenario {
  ParticleContainer initialParticles(Field& E, Field& B, Field& V, Field& R);
  void newTimestep(int input_file_counter, int step, double time, ParticleContainer& particles, Interpolated_Field& E, Interpolated_Field& B,
        Interpolated_Field& V, Interpolated_Field& R);
  void finalize(ParticleContainer& particles, Field& E, Field& B, Field& V, Field& R);

  distributionScenario() {needV = false; needRho = false;};
};

// Inject particles in the tail continuously, see where they precipitate
struct precipitationScenario : Scenario {
  void newTimestep(int input_file_counter, int step, double time, ParticleContainer& particles, Interpolated_Field& E, Interpolated_Field& B,
        Interpolated_Field& V, Interpolated_Field& R);
  void afterPush(int step, double time, ParticleContainer& particles, Interpolated_Field& E, Interpolated_Field& B, Interpolated_Field& V, Interpolated_Field& R);

  precipitationScenario() {needV = true; needRho = false;};
};

// For interactive usage from analysator: read positions and velocities from stdin, push those particles.
struct analysatorScenario : Scenario {
  ParticleContainer initialParticles(Field& E, Field& B, Field& V, Field& R);
  void newTimestep(int input_file_counter, int step, double time, ParticleContainer& particles, Interpolated_Field& E, Interpolated_Field& B,
        Interpolated_Field& V, Interpolated_Field& R);
  analysatorScenario() {needV = false; needRho = false;};
};

// Analyzation of shock reflectivity in radial IMF runs
struct shockReflectivityScenario : Scenario {

  LinearHistogram2D transmitted;
  LinearHistogram2D reflected;

  void newTimestep(int input_file_counter, int step, double time, ParticleContainer& particles, Interpolated_Field& E, Interpolated_Field& B,
        Interpolated_Field& V, Interpolated_Field& R);
  void afterPush(int step, double time, ParticleContainer& particles, Interpolated_Field& E, Interpolated_Field& B, Interpolated_Field& V, Interpolated_Field& R);
  void finalize(ParticleContainer& particles, Field& E, Field& B, Field& V, Field& R);

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

  ParticleContainer initialParticles(Field& E, Field& B, Field& V, Field& R);
  void newTimestep(int input_file_counter, int step, double time, ParticleContainer& particles, Interpolated_Field& E, Interpolated_Field& B,
        Interpolated_Field& V, Interpolated_Field& R);
  void afterPush(int step, double time, ParticleContainer& particles, Interpolated_Field& E, Interpolated_Field& B, Interpolated_Field& V, Interpolated_Field& R);
  void finalize(ParticleContainer& particles, Field& E, Field& B, Field& V, Field& R);

  ipShockScenario() {needV = true; needRho = false;};
};

// Initialize particles at a given distance away from the fitted bow shock shape, track their precipitation upstream or downstream
struct InjectionScenario : Scenario {
  FILE * initFile;
  FILE * traFile;
  FILE * refFile;
  FILE * lostFile;

  FILE * meet_rho_File;
  FILE * meet_TNBS_File;
  FILE * meet_Mms_File;
  FILE * meet_tbn_File;
  FILE * meet_flipmu_File;
  FILE * minxFile;
  FILE * boostFile;

  bool * fs_hasmet_rho;
  bool * fs_hasmet_TNBS;
  bool * fs_hasmet_Mms;
  bool * fs_hasmet_tbn;
  bool * fs_hasmet_flipmu;

  int particlespertimestep;

  ParticleContainer initialParticles(Field& E, Field& B, Field& V, Field& R);
  void beforePush(ParticleContainer& particles, Interpolated_Field& E, Interpolated_Field& B, Interpolated_Field& V, Interpolated_Field& R);
  void newTimestep(int input_file_counter, int step, double time, ParticleContainer& particles, Interpolated_Field& E, Interpolated_Field& B,
        Interpolated_Field& V, Interpolated_Field& R);
  void afterPush(int step, double time, ParticleContainer& particles, Interpolated_Field& E, Interpolated_Field& B, Interpolated_Field& V, Interpolated_Field& R);
  void finalize(ParticleContainer& particles, Field& E, Field& B, Field& V, Field& R);

  InjectionScenario() {needV = true; needRho = true;};
};

template<typename T> Scenario* createScenario() {
  return new T;
}

Scenario* createScenario(std::string name);
