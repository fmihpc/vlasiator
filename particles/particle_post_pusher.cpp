/*
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

int main(int argc, char** argv) {

   MPI::Init(argc, argv);

   /* Parse commandline and config*/
   Readparameters parameters(argc, argv, MPI_COMM_WORLD);
   ParticleParameters::addParameters();
   parameters.parse(false);  // Parse parameters and don't require run_config
   ParticleParameters::getParameters();

   /* Read starting fields from specified input file */
   std::string filename_pattern = ParticleParameters::input_filename_pattern;
   char filename_buffer[256];

   int input_file_counter=floor(ParticleParameters::start_time / ParticleParameters::input_dt);
   Field E[2],B[2],V;
   std::cerr << "Loading first file with index " << ParticleParameters::start_time / ParticleParameters::input_dt << std::endl;
   snprintf(filename_buffer,256,filename_pattern.c_str(),input_file_counter-1);
   readfields(filename_buffer,E[1],B[1],V);
   E[0]=E[1]; B[0]=B[1];

   /* Init particles */

   double dt=ParticleParameters::dt;
   double maxtime=ParticleParameters::end_time - ParticleParameters::start_time;
   int maxsteps = maxtime/dt;

   Scenario* scenario = create_scenario(ParticleParameters::mode);
   std::vector<Particle> particles = scenario->initial_particles(E[0],B[0],V);

   std::cerr << "Pushing " << particles.size() << " particles for " << maxsteps << " steps..." << std::endl;
        std::cerr << "[                                                                        ]\x0d[";

   /* Push them around */
   for(int step=0; step<maxsteps; step++) {

      bool newfile;
      /* Load newer fields, if neccessary */
      if(step >= 0) {
         newfile = read_next_timestep(filename_pattern, ParticleParameters::start_time + step*dt, 1,E[0], E[1], B[0], B[1], V, scenario->needV, input_file_counter);
      } else {
         newfile = read_next_timestep(filename_pattern, ParticleParameters::start_time + step*dt, -1,E[1], E[0], B[1], B[0], V, scenario->needV, input_file_counter);
      }

      Interpolated_Field cur_E(E[0],E[1],step*dt);
      Interpolated_Field cur_B(B[0],B[1],step*dt);

      // If a new timestep has been opened, add a new bunch of particles
      if(newfile) {
        scenario->new_timestep(input_file_counter, step, step*dt, particles, cur_E, cur_B, V);
      }

      scenario->before_push(particles,cur_E,cur_B,V);

      #pragma omp parallel for
      for(unsigned int i=0; i< particles.size(); i++) {

          if(vector_length(particles[i].x) == 0) {
            // Skip disabled particles.
            continue;
          }

         /* Get E- and B-Field at their position */
         Vec3d Eval,Bval;

         Eval = cur_E(particles[i].x);
         Bval = cur_B(particles[i].x);

         /* Push them around */
         particles[i].push(Bval,Eval,dt);

      }

			// Remove all particles that have left the simulation box after this step
			// (unfortunately, this can not be done in parallel, so it better be fast!)
			double inner_boundary_min[3]; // Inner boundary: minimum + 2*dx
			inner_boundary_min[0] = B[0].min[0] + 2.*B[0].dx[0];
			inner_boundary_min[1] = B[0].min[1] + 2.*B[0].dx[1];
			inner_boundary_min[2] = B[0].min[2] + 2.*B[0].dx[2];
			double inner_boundary_max[3]; // Inner boundary maximum - 2*dx
			inner_boundary_max[0] = B[0].max[0] - 2.*B[0].dx[0];
			inner_boundary_max[1] = B[0].max[1] - 2.*B[0].dx[1];
			inner_boundary_max[2] = B[0].max[2] - 2.*B[0].dx[2];
			for(std::vector<Particle>::iterator i = particles.begin(); i != particles.end(); ) {
				bool do_delete=false;
				Vec3d x =i->x;
				if(B[0].cells[0] > 1 && ( x[0] <= inner_boundary_min[0] || x[0] >= inner_boundary_max[0])) {
					switch(ParticleParameters::boundary_behaviour_x) {
						case DELETE:
							// mark for deletion
							do_delete=true;
							break;
						case REFLECT:
							i->v *= Vec3d(-1.,1.,1.);
							break;
						case PERIODIC:
							if(i->x[0] <= inner_boundary_min[0]) {
								i->x += Vec3d(inner_boundary_max[0]-inner_boundary_min[0],0,0);
							} else if(i->x[0] >= inner_boundary_max[0]) {
								i->x -= Vec3d(inner_boundary_max[0]-inner_boundary_min[0],0,0);
							}
							break;
					}
				}
				if(B[0].cells[1] > 1 && (x[1] <= inner_boundary_min[1] || x[1] >= inner_boundary_max[1])) {
					switch(ParticleParameters::boundary_behaviour_y) {
						case DELETE:
							// mark for deletion
							do_delete=true;
							break;
						case REFLECT:
							i->v *= Vec3d(1.,-1.,1.);
							break;
						case PERIODIC:
							if(i->x[1] <= inner_boundary_min[1]) {
								i->x += Vec3d(0,inner_boundary_max[1]-inner_boundary_min[1],0);
							} else if(i->x[1] >= inner_boundary_max[1]) {
								i->x -= Vec3d(0,inner_boundary_max[1]-inner_boundary_min[1],0);
							}
							break;
					}
				}
			  if(B[0].cells[2] > 1 && (x[2] <= inner_boundary_min[2] || x[2] >= inner_boundary_max[2])) {
					switch(ParticleParameters::boundary_behaviour_y) {
						case DELETE:
							// mark for deletion
							do_delete=true;
							break;
						case REFLECT:
							i->v *= Vec3d(1.,1.,-1.);
							break;
						case PERIODIC:
							if(i->x[2] <= inner_boundary_min[2]) {
								i->x += Vec3d(0,0,inner_boundary_max[2]-inner_boundary_min[2]);
							} else if(i->x[2] >= inner_boundary_max[2]) {
								i->x -= Vec3d(0,0,inner_boundary_max[2]-inner_boundary_min[2]);
							}
							break;
					}
				}

				if(do_delete) {
					particles.erase(i);
				} else {
					i++;
				}
			}

      scenario->after_push(step, step*dt, particles, cur_E, cur_B, V);

      /* Draw progress bar */
      if((step % (maxsteps/71))==0) {
         std::cerr << "=";
      }
   }

   scenario->finalize(particles,E[1],B[1],V);

   std::cerr << std::endl;

   MPI::Finalize();
   return 0;
}
