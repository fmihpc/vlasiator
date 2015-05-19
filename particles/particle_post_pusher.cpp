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
#include "physconst.h"
#include "readfields.h"
#include "distribution.h"
#include "particleparameters.h"
#include "../readparameters.h"

// "mode" that the particle pusher operates in
enum pusher_mode {
   single = 0,    // Single particle trajectory tracking
   distribution,  // Particle distribution created at a given point
   analysator     // Particle distribution input as ascii, positions written to stdout as ascii.
};

int main(int argc, char** argv) {

   MPI::Init(argc, argv);

   /* Parse commandline and config*/
   Readparameters parameters(argc, argv, MPI_COMM_WORLD);
   ParticleParameters::addParameters();
   parameters.parse();
   ParticleParameters::getParameters();

   /* Read starting fields from specified input file */
   std::string filename_pattern = ParticleParameters::input_filename_pattern;
   char filename_buffer[256];

   int input_file_counter=floor(ParticleParameters::start_time / ParticleParameters::input_dt);
   Field E[2],B[2],V;
   snprintf(filename_buffer,256,filename_pattern.c_str(),input_file_counter-1);
   readfields(filename_buffer,E[1],B[1],V);
   E[0]=E[1]; B[0]=B[1];

   /* Init particles */
   std::vector<Particle> particles;
   pusher_mode mode;

   double dt=ParticleParameters::dt;
   double maxtime=ParticleParameters::end_time - ParticleParameters::start_time;
   int maxsteps = maxtime/dt;

   if(ParticleParameters::mode == "single") {

      Vec3d vpos(ParticleParameters::init_x, ParticleParameters::init_y, ParticleParameters::init_z);
      /* Look up builk velocity in the V-field */
      Vec3d bulk_vel = V(vpos);

      mode = single;
      particles.push_back(Particle(PhysicalConstantsSI::mp, PhysicalConstantsSI::e, vpos, bulk_vel));
   } else if(ParticleParameters::mode == "distribution") {

      mode = distribution;

      std::default_random_engine generator(ParticleParameters::random_seed);
      Distribution* velocity_distribution=ParticleParameters::distribution(generator);

      Vec3d vpos(ParticleParameters::init_x, ParticleParameters::init_y, ParticleParameters::init_z);

      /* Look up builk velocity in the V-field */
      Vec3d bulk_vel = V(vpos);

      for(unsigned int i=0; i< ParticleParameters::num_particles; i++) {
         /* Create a particle with velocity drawn from the given distribution ... */
         Particle p = velocity_distribution->next_particle();
         /* Shift it by the bulk velocity ... */
         p.v += bulk_vel;
         /* And put it in place. */
         p.x=vpos;
         particles.push_back(p);
      }

      delete velocity_distribution;
   } else if(ParticleParameters::mode == "analysator") {

     mode = analysator;

     std::cerr << "Reading initial particle data from stdin" << std::endl
        << "(format: x y z vx vy vz)" << std::endl;

     while(std::cin) {
        double x0,x1,x2,v0,v1,v2;
        std::cin >> x0 >> x1 >> x2 >> v0 >> v1 >> v2;
        if(std::cin) {
           particles.push_back(Particle(PhysicalConstantsSI::mp, PhysicalConstantsSI::e, Vec3d(x0,x1,x2), Vec3d(v0,v1,v2)));
        }
     }

   } else {
      std::cerr << "Config option particles.mode not set correctly!" << std::endl;
      return 1;
   }

   /* Dump the initial state */
   if(mode == distribution) {
      write_particles(particles, "particles_initial.vlsv");
   }

   std::cerr << "Pushing " << particles.size() << " particles for " << maxsteps << " steps..." << std::endl;
        std::cerr << "[                                                                        ]\x0d[";

   /* Push them around */
   for(int step=0; step<maxsteps; step++) {

      bool newfile;
      /* Load newer fields, if neccessary */
      if(step >= 0) {
         newfile = read_next_timestep(filename_pattern, ParticleParameters::start_time + step*dt, 1,E[0], E[1], B[0], B[1], input_file_counter);
      } else {
         newfile = read_next_timestep(filename_pattern, ParticleParameters::start_time + step*dt, -1,E[1], E[0], B[1], B[0], input_file_counter);
      }

      Interpolated_Field cur_E(E[0],E[1],step*dt);
      Interpolated_Field cur_B(B[0],B[1],step*dt);

      //#pragma omp parallel for
      for(unsigned int i=0; i< particles.size(); i++) {
         /* Get E- and B-Field at their position */
         Vec3d Eval,Bval;

         Eval = cur_E(particles[i].x);
         Bval = cur_B(particles[i].x);

         /* Push them around */
         particles[i].push(Bval,Eval,dt);
      }

      /* Write output */
      if(mode == distribution && newfile) {
         snprintf(filename_buffer,256, ParticleParameters::output_filename_pattern.c_str(),input_file_counter-1);
         write_particles(particles, filename_buffer);
      } else if(mode == analysator && newfile) {
        for(unsigned int i=0; i< particles.size(); i++) {
           Vec3d& x = particles[i].x;
           Vec3d& v = particles[i].v;
          std::cout << i << " " << step*dt << "\t" <<  x[0] << " " << x[1] << " " << x[2] << "\t" << v[0] << " " << v[1] << " " << v[2] << std::endl;
        }
      }

      /* Draw progress bar */
      if((step % (maxsteps/71))==0) {
         std::cerr << "=";
      }
   }

   /* Dump the final state */
   if(mode == distribution) {
      write_particles(particles, "particles_final.vlsv");
   }

   std::cerr << std::endl;

   MPI::Finalize();
   return 0;
}
