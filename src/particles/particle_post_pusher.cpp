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
 *
 * Postprocessing particle trajectory analyzer
 * for Vlasiator
 *
 */
#include <mpi.h>
#include <iostream>
#include <random>
#include <string.h>
#include "particles.h"
#include "field.h"
#include "physconst.h"
#include "readfields.h"
#include "distribution.h"
#include "particleparameters.h"
#include "../readparameters.h"
#include "scenario.h"
#include "boundaries.h"

int main(int argc, char** argv) {

   MPI_Init(&argc, &argv);

   /* Parse commandline and config*/
   Readparameters parameters(argc, argv, MPI_COMM_WORLD);
   ParticleParameters::addParameters();
   parameters.parse(false);  // Parse parameters and don't require run_config
   if(!ParticleParameters::getParameters()) {
      std::cerr << "Parsing parameters failed, aborting." << std::endl;
      std::cerr << "Did you add a --run_config=file.cfg parameter?" << std::endl;
      return 1;
   }

   /* Read starting fields from specified input file */
   std::string filename_pattern = ParticleParameters::input_filename_pattern;
   char filename_buffer[256];

   int input_file_counter=floor(ParticleParameters::start_time / ParticleParameters::input_dt);
   Field E[2],B[2],V;
   std::cerr << "Loading first file with index " << ParticleParameters::start_time / ParticleParameters::input_dt
      << std::endl;
   snprintf(filename_buffer,256,filename_pattern.c_str(),input_file_counter-1);
   E[0].dimension[0] = E[1].dimension[0] = B[0].dimension[0] = B[1].dimension[0] = V.dimension[0] = ParticleParameters::boundary_behaviour_x;
   E[0].dimension[1] = E[1].dimension[1] = B[0].dimension[1] = B[1].dimension[1] = V.dimension[1] = ParticleParameters::boundary_behaviour_y;
   E[0].dimension[2] = E[1].dimension[2] = B[0].dimension[2] = B[1].dimension[2] = V.dimension[2] = ParticleParameters::boundary_behaviour_z;
   readfields(filename_buffer,E[1],B[1],V);
   E[0]=E[1]; B[0]=B[1];

   // Set boundary conditions based on sizes
   if(B[0].dimension[0]->cells <= 1) {
      delete ParticleParameters::boundary_behaviour_x;
      ParticleParameters::boundary_behaviour_x = createBoundary<CompactSpatialDimension>(0);
   }
   if(B[0].dimension[1]->cells <= 1) {
      delete ParticleParameters::boundary_behaviour_y;
      ParticleParameters::boundary_behaviour_y = createBoundary<CompactSpatialDimension>(1);
   }
   if(B[0].dimension[2]->cells <= 1) {
      delete ParticleParameters::boundary_behaviour_z;
      ParticleParameters::boundary_behaviour_z = createBoundary<CompactSpatialDimension>(2);
   }

   // Make sure updated boundary conditions are also correctly known to the fields
   E[0].dimension[0] = E[1].dimension[0] = B[0].dimension[0] = B[1].dimension[0] = V.dimension[0] = ParticleParameters::boundary_behaviour_x;
   E[0].dimension[1] = E[1].dimension[1] = B[0].dimension[1] = B[1].dimension[1] = V.dimension[1] = ParticleParameters::boundary_behaviour_y;
   E[0].dimension[2] = E[1].dimension[2] = B[0].dimension[2] = B[1].dimension[2] = V.dimension[2] = ParticleParameters::boundary_behaviour_z;

   ParticleParameters::boundary_behaviour_x->setExtent(B[0].dimension[0]->min, B[0].dimension[0]->max, B[0].dimension[0]->cells);
   ParticleParameters::boundary_behaviour_y->setExtent(B[0].dimension[1]->min, B[0].dimension[1]->max, B[0].dimension[1]->cells);
   ParticleParameters::boundary_behaviour_z->setExtent(B[0].dimension[2]->min, B[0].dimension[2]->max, B[0].dimension[2]->cells);

   /* Init particles */
   double dt=ParticleParameters::dt;
   double maxtime=ParticleParameters::end_time - ParticleParameters::start_time;
   int maxsteps = maxtime/dt;

   Scenario* scenario = createScenario(ParticleParameters::mode);
   ParticleContainer particles = scenario->initialParticles(E[0],B[0],V);

   std::cerr << "Pushing " << particles.size() << " particles for " << maxsteps << " steps..." << std::endl;
   std::cerr << "[                                                                        ]\x0d[";

   /* Push them around */
   for(int step=0; step<maxsteps; step++) {

      bool newfile;
      /* Load newer fields, if neccessary */
      if(step >= 0) {
         newfile = readNextTimestep(filename_pattern, ParticleParameters::start_time + step*dt, 1,E[0], E[1],
               B[0], B[1], V, scenario->needV, input_file_counter);
      } else {
         newfile = readNextTimestep(filename_pattern, ParticleParameters::start_time + step*dt, -1,E[1], E[0],
               B[1], B[0], V, scenario->needV, input_file_counter);
      }

      Interpolated_Field cur_E(E[0],E[1],ParticleParameters::start_time + step*dt);
      Interpolated_Field cur_B(B[0],B[1],ParticleParameters::start_time + step*dt);

      // If a new timestep has been opened, add a new bunch of particles
      if(newfile) {
         scenario->newTimestep(input_file_counter, step, step*dt, particles, cur_E, cur_B, V);
      }

      scenario->beforePush(particles,cur_E,cur_B,V);

#pragma omp parallel for
      for(unsigned int i=0; i< particles.size(); i++) {

         if(isnan(vector_length(particles[i].x))) {
            // Skip disabled particles.
            continue;
         }

         /* Get E- and B-Field at their position */
         Vec3d Eval,Bval;

         Eval = cur_E(particles[i].x);
         Bval = cur_B(particles[i].x);

         if(dt < 0) {
           // If propagating backwards in time, flip B-field pseudovector
           Bval *= -1;
         }

         /* Push them around */
         particles[i].push(Bval,Eval,dt);

      }

      // Remove all particles that have left the simulation box after this step
      // (unfortunately, this can not be done in parallel, so it better be fast!)
      for(auto i = particles.begin(); i != particles.end(); ) {

         // Boundaries are allowed to mangle the particles here.
         // If they return false, particles are deleted.
         bool do_erase = false;
         if(!ParticleParameters::boundary_behaviour_x->handleParticle(*i)) {
            do_erase = true;
         }
         if(!ParticleParameters::boundary_behaviour_y->handleParticle(*i)) {
            do_erase = true;
         }
         if(!ParticleParameters::boundary_behaviour_z->handleParticle(*i)) {
            do_erase = true;
         }
         if(do_erase) {
            particles.erase(i);
         } else {
            i++;
         }
      }

      scenario->afterPush(step, step*dt, particles, cur_E, cur_B, V);

      /* Draw progress bar */
      if((step % (maxsteps/71))==0) {
         std::cerr << "=";
      }
   }

   scenario->finalize(particles,E[1],B[1],V);

   std::cerr << std::endl;

   MPI_Finalize();
   return 0;
}
