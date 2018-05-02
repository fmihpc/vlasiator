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

   MPI::Init(argc, argv);

   /* Parse commandline and config*/
   Readparameters parameters(argc, argv, MPI_COMM_WORLD);
   ParticleParameters::addParameters();
   parameters.parse(false);  // Parse parameters and don't require run_config
   if(!ParticleParameters::getParameters()) {
      std::cerr << "Parsing parameters failed, aborting." << std::endl;
      return 1;
   }

   /* Read starting fields from specified input file */
   std::string filename_pattern = ParticleParameters::input_filename_pattern;
   char filename_buffer[256];

   int input_file_counter=floor(ParticleParameters::start_time / ParticleParameters::input_dt);
   Field E[2],B[2],V[2],R[2];
   std::cerr << "Loading first file with index " << ParticleParameters::start_time / ParticleParameters::input_dt
      << std::endl;
   snprintf(filename_buffer,256,filename_pattern.c_str(),input_file_counter-1);
   E[0].dimension[0] = E[1].dimension[0] = B[0].dimension[0] = B[1].dimension[0] = V[0].dimension[0] = V[1].dimension[0] = R[0].dimension[0] = R[1].dimension[0] = ParticleParameters::boundary_behaviour_x;
   E[0].dimension[1] = E[1].dimension[1] = B[0].dimension[1] = B[1].dimension[1] = V[0].dimension[1] = V[1].dimension[1] = R[0].dimension[1] = R[1].dimension[1] = ParticleParameters::boundary_behaviour_y;
   E[0].dimension[2] = E[1].dimension[2] = B[0].dimension[2] = B[1].dimension[2] = V[0].dimension[2] = V[1].dimension[2] = R[0].dimension[2] = R[1].dimension[2] = ParticleParameters::boundary_behaviour_z;
   readfields(filename_buffer,E[1],B[1],V[1], R[1]);
   E[0]=E[1]; B[0]=B[1]; V[0]=V[1];  R[0]=R[1];

   // Set boundary conditions based on sizes
   if(B[0].dimension[0]->cells <= 1) {
      delete ParticleParameters::boundary_behaviour_x;
      ParticleParameters::boundary_behaviour_x = createBoundary<CompactSpatialDimension>(0);      
      //ParticleParameters::boundary_behaviour_x->cells=1;
   }
   if(B[0].dimension[1]->cells <= 1) {
      delete ParticleParameters::boundary_behaviour_y;
      ParticleParameters::boundary_behaviour_y = createBoundary<CompactSpatialDimension>(1);
      //ParticleParameters::boundary_behaviour_y->cells=1;
   }
   if(B[0].dimension[2]->cells <= 1) {
      delete ParticleParameters::boundary_behaviour_z;
      ParticleParameters::boundary_behaviour_z = createBoundary<CompactSpatialDimension>(2);
      //ParticleParameters::boundary_behaviour_z->cells=1;
   }

   // Make sure updated boundary conditions are also correctly known to the fields
   E[0].dimension[0] = E[1].dimension[0] = B[0].dimension[0] = B[1].dimension[0] = V[0].dimension[0] = V[1].dimension[0] = R[0].dimension[0] = R[1].dimension[0] = ParticleParameters::boundary_behaviour_x;
   E[0].dimension[1] = E[1].dimension[1] = B[0].dimension[1] = B[1].dimension[1] = V[0].dimension[1] = V[1].dimension[1] = R[0].dimension[1] = R[1].dimension[1] = ParticleParameters::boundary_behaviour_y;
   E[0].dimension[2] = E[1].dimension[2] = B[0].dimension[2] = B[1].dimension[2] = V[0].dimension[2] = V[1].dimension[2] = R[0].dimension[2] = R[1].dimension[2] = ParticleParameters::boundary_behaviour_z;

   ParticleParameters::boundary_behaviour_x->setExtent(B[0].dimension[0]->min, B[0].dimension[0]->max, B[0].dimension[0]->cells);
   ParticleParameters::boundary_behaviour_y->setExtent(B[0].dimension[1]->min, B[0].dimension[1]->max, B[0].dimension[1]->cells);
   ParticleParameters::boundary_behaviour_z->setExtent(B[0].dimension[2]->min, B[0].dimension[2]->max, B[0].dimension[2]->cells);

   /* Init particles */
   double dt=ParticleParameters::dt;
   double maxtime=ParticleParameters::end_time - ParticleParameters::start_time;
   int maxsteps = maxtime/dt;
   Real filetime = ParticleParameters::start_time;    /* for use in static fields */

   Scenario* scenario = createScenario(ParticleParameters::mode);
   ParticleContainer particles = scenario->initialParticles(E[0],B[0],V[0],R[0]);

   std::cerr << "Pushing " << particles.size() << " particles for " << maxsteps << " steps..." << std::endl;
   std::cerr << "[                                                                        ]\x0d[";

   // Initial fields
   Interpolated_Field cur_E(E[0],E[1],ParticleParameters::start_time);
   Interpolated_Field cur_B(B[0],B[1],ParticleParameters::start_time);
   Interpolated_Field cur_V(V[0],V[1],ParticleParameters::start_time);
   Interpolated_Field cur_R(R[0],R[1],ParticleParameters::start_time);
   if (ParticleParameters::staticfields) {
     E[1].time = ParticleParameters::start_time + ParticleParameters::input_dt;
     B[1].time = ParticleParameters::start_time + ParticleParameters::input_dt;
     V[1].time = ParticleParameters::start_time + ParticleParameters::input_dt;
     R[1].time = ParticleParameters::start_time + ParticleParameters::input_dt;
     bool initfields = readNextTimestep(filename_pattern, ParticleParameters::start_time, 1,E[0], E[1],
					B[0], B[1], V[0], V[1], R[0], R[1], scenario->needV, scenario->needRho, input_file_counter);
   }

   /* Push them around */
   for(int step=0; step<maxsteps; step++) {
     
     bool newfile = false;
     if (!ParticleParameters::staticfields) {
       /* Load newer fields, if neccessary */
       if(step >= 0) {
	 newfile = readNextTimestep(filename_pattern, ParticleParameters::start_time + step*dt, 1,E[0], E[1],
				    B[0], B[1], V[0], V[1], R[0], R[1], scenario->needV, scenario->needRho, input_file_counter);
       } else {
	 newfile = readNextTimestep(filename_pattern, ParticleParameters::start_time + step*dt, -1,E[1], E[0],
				    B[1], B[0], V[0], V[1], R[0], R[1], scenario->needV, scenario->needRho, input_file_counter);
       }

       cur_E.setfields(E[0],E[1],ParticleParameters::start_time + step*dt);
       cur_B.setfields(B[0],B[1],ParticleParameters::start_time + step*dt);       
       cur_V.setfields(V[0],V[1],ParticleParameters::start_time + step*dt);       
       cur_R.setfields(R[0],R[1],ParticleParameters::start_time + step*dt);       
     } else {
       if (ParticleParameters::start_time + step*dt > filetime) {
	 if (step >= 0) {
	   input_file_counter += 1;
	 } else {
	   input_file_counter += -1;
	 }
	 filetime = input_file_counter * ParticleParameters::input_dt;
	 newfile = true;
       }
     }
     // If a new timestep has been opened, add a new bunch of particles
     if(newfile) {
       std::cerr << "New file " << input_file_counter<<" step "<< step<<" time "<< step*dt<<std::endl;
       //std::cerr << "   fields "<<ParticleParameters::staticfields<<std::endl;

       scenario->newTimestep(input_file_counter, step, step*dt, particles, cur_E, cur_B, cur_V, cur_R);
     }

     scenario->beforePush(particles,cur_E,cur_B,cur_V, cur_R);

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
	 //std::cerr<<"Eval "<<Eval[0]<<" "<<Eval[1]<<" "<<Eval[2]<<"Bval "<<Bval[0]<<" "<<Bval[1]<<" "<<Bval[2]<<std::endl;

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

      scenario->afterPush(step, step*dt, particles, cur_E, cur_B, cur_V, cur_R);

      /* Draw progress bar */
      if((step % (maxsteps/71))==0) {
         std::cerr << "=";
      }
   }

   scenario->finalize(particles,E[1],B[1],V[1],R[1]);

   std::cerr << std::endl;

   MPI::Finalize();
   return 0;
}
