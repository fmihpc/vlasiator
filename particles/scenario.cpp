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
#include <iostream>
#include "scenario.h"

ParticleContainer singleParticleScenario::initialParticles(Field& E, Field& B, Field& V, Field& R) {

   ParticleContainer particles;

   Vec3d vpos(ParticleParameters::init_x, ParticleParameters::init_y, ParticleParameters::init_z);
   /* Look up builk velocity in the V-field */
   Vec3d bulk_vel = V(vpos);

   particles.push_back(Particle(PhysicalConstantsSI::mp, PhysicalConstantsSI::e, vpos, bulk_vel));

   return particles;
}

void singleParticleScenario::afterPush(int step, double time, ParticleContainer& particles, 
      Interpolated_Field& E, Interpolated_Field& B, Interpolated_Field& V, Interpolated_Field& R) {

   Vec3d& x = particles[0].x;
   Vec3d& v = particles[0].v;

   std::cout << 0 << " " << time << "\t" <<  x[0] << " " << x[1] << " " << x[2] << "\t"
      << v[0] << " " << v[1] << " " << v[2] << std::endl;
}



ParticleContainer distributionScenario::initialParticles(Field& E, Field& B, Field& V, Field& R) {

   ParticleContainer particles;

   std::default_random_engine generator(ParticleParameters::random_seed);
   Distribution* velocity_distribution=ParticleParameters::distribution(generator);

   Vec3d vpos(ParticleParameters::init_x, ParticleParameters::init_y, ParticleParameters::init_z);

   // Look up bulk velocity in the V-field 
   Vec3d bulk_vel = V(vpos);
   // Lookup B vector
   Vec3d Bval = B(vpos);

   // Build velocity space coordinate axes
   Vec3d velspace_x, velspace_y, velspace_z;
   if(ParticleParameters::vel_BcrossVframe) {
     // BxV frame
     velspace_x = normalize_vector(Bval);
     velspace_y = normalize_vector(cross_product(Bval,bulk_vel));
     velspace_z = normalize_vector(cross_product(Bval,velspace_y));
   } else {
     // Cartesian simulation frame
     velspace_x = Vec3d(1,0,0);
     velspace_y = Vec3d(0,1,0);
     velspace_z = Vec3d(0,0,1);
   }

   for(unsigned int i=0; i< ParticleParameters::num_particles; i++) {
      // Create a particle with velocity drawn from the given distribution ...
      Particle p = velocity_distribution->next_particle();

      // Potentially give it a drift velocity
      p.v += Vec3d(ParticleParameters::parallelDriftVel,ParticleParameters::perpDriftVel1,ParticleParameters::perpDriftVel2);

      // Rotate it into the chosen coordinate frame
      p.v = p.v[0] * velspace_x + p.v[1] * velspace_y + p.v[2] * velspace_z;
      // Shift it by the bulk velocity ...
//      p.v += bulk_vel;
      // And put it in place.
      p.x=vpos;
      particles.push_back(p);
   }

   delete velocity_distribution;

   return particles;
}

void distributionScenario::newTimestep(int input_file_counter, int step, double time, ParticleContainer& particles,
      Interpolated_Field& E, Interpolated_Field& B, Interpolated_Field& V, Interpolated_Field& R) {

   char filename_buffer[256];

   // /!\ Go one step opposite to the propagation time direction
   snprintf(filename_buffer,256, ParticleParameters::output_filename_pattern.c_str(),input_file_counter - ParticleParameters::propagation_direction);
   writeParticles(particles, filename_buffer);
}

void distributionScenario::finalize(ParticleContainer& particles, Field& E, Field& B, Field& V, Field& R) {
   writeParticles(particles, "particles_final.vlsv");
}

void precipitationScenario::afterPush(int step, double time, ParticleContainer& particles,
      Interpolated_Field& E, Interpolated_Field& B, Interpolated_Field& V, Interpolated_Field& R) {

   for(unsigned int i=0; i<particles.size(); i++) {

      if(isnan(vector_length(particles[i].x))) {
         // skip disabled particles
         continue;
      }

      // Check if the particle hit a boundary. If yes, mark it as disabled.
      // Original starting x of this particle
      double start_pos = ParticleParameters::precip_start_x +
         ((double)(i%ParticleParameters::num_particles))/ParticleParameters::num_particles *
          (ParticleParameters::precip_stop_x - ParticleParameters::precip_start_x);
      int start_timestep = i / ParticleParameters::num_particles;
      if(vector_length(particles[i].x) <= ParticleParameters::precip_inner_boundary) {

         // Record latitude and energy
         double latitude = atan2(particles[i].x[2],particles[i].x[0]);
         printf("%u %i %lf %lf %lf\n",i, start_timestep, start_pos, latitude, .5*particles[i].m *
               dot_product(particles[i].v,particles[i].v)/PhysicalConstantsSI::e);

         // Disable by setting position to NaN and velocity to 0
         particles[i].x = Vec3d(std::numeric_limits<double>::quiet_NaN(),0.,0.);
         particles[i].v = Vec3d(0,0,0);
      } else if (particles[i].x[0] <= ParticleParameters::precip_start_x) {

         // Record marker value for lost particle
         printf("%u %i %lf -5. -1.\n", i, start_timestep, start_pos);
         // Disable by setting position to NaN and velocity to 0
         particles[i].x = Vec3d(std::numeric_limits<double>::quiet_NaN(),0.,0.);
         particles[i].v = Vec3d(0,0,0);
      }
   }
}

void precipitationScenario::newTimestep(int input_file_counter, int step, double time, ParticleContainer& particles,
      Interpolated_Field& E, Interpolated_Field& B, Interpolated_Field& V, Interpolated_Field& R) {

   // Create particles along the negative x-axis, from inner boundary
   // up to outer one
   for(unsigned int i=0; i< ParticleParameters::num_particles; i++) {

      // Choose x coordinate
      double start_x = ParticleParameters::precip_start_x +
         ((double)i)/ParticleParameters::num_particles *
          (ParticleParameters::precip_stop_x - ParticleParameters::precip_start_x);
      Vec3d pos(start_x,0,0);

      // Find cell with minimum B value in this plane
      double min_B = 99999999999.;
      for(double z=-1e7; z<1e7; z+=1e5) {
         Vec3d candidate_pos(start_x,0,z);
         double B_here = vector_length(B(candidate_pos));
         if(B_here < min_B) {
            pos = candidate_pos;
            min_B = B_here;
         }
      }

      // Add a particle at this location, with bulk velocity at its starting point
      particles.push_back(Particle(PhysicalConstantsSI::mp, PhysicalConstantsSI::e, pos, V(pos)));
   }

   // Write out the state
   char filename_buffer[256];
   
   // /!\ Go one step opposite to the propagation time direction
   snprintf(filename_buffer,256, ParticleParameters::output_filename_pattern.c_str(),input_file_counter - ParticleParameters::propagation_direction);
   writeParticles(particles, filename_buffer);
}


ParticleContainer analysatorScenario::initialParticles(Field& E, Field& B, Field& V, Field& R) {

   ParticleContainer particles;

   std::cerr << "Reading initial particle data from stdin" << std::endl
      << "(format: x y z vx vy vz)" << std::endl;

   while(std::cin) {
      double x0,x1,x2,v0,v1,v2;
      std::cin >> x0 >> x1 >> x2 >> v0 >> v1 >> v2;
      if(std::cin) {
         particles.push_back(Particle(PhysicalConstantsSI::mp, PhysicalConstantsSI::e,
                  Vec3d(x0,x1,x2), Vec3d(v0,v1,v2)));
      }
   }

   return particles;
}

void analysatorScenario::newTimestep(int input_file_counter, int step, double time, ParticleContainer& particles,
      Interpolated_Field& E, Interpolated_Field& B, Interpolated_Field& V, Interpolated_Field& R) {

   for(unsigned int i=0; i< particles.size(); i++) {
      Vec3d& x = particles[i].x;
      Vec3d& v = particles[i].v;
      std::cout << i << " " << time << "\t" <<  x[0] << " " << x[1] << " " << x[2] << "\t"
         << v[0] << " " << v[1] << " " << v[2] << std::endl;
   }
}

void shockReflectivityScenario::newTimestep(int input_file_counter, int step, double time,
      ParticleContainer& particles, Interpolated_Field& E, Interpolated_Field& B, Interpolated_Field& V, Interpolated_Field& R) {

   const int num_points = 200;

   std::default_random_engine generator(ParticleParameters::random_seed+step);
   Distribution* velocity_distribution=ParticleParameters::distribution(generator);

   // Create particles along a parabola, in front of the shock
   for(unsigned int i=0; i< num_points; i++) {

      // Choose y coordinate
      double start_y = ParticleParameters::reflect_start_y +
         ((double)i)/num_points *
          (ParticleParameters::reflect_stop_y - ParticleParameters::reflect_start_y);

      // Calc x-coordinate from it
      double x = start_y / ParticleParameters::reflect_start_y;
      x*=-x;
      x *= ParticleParameters::reflect_y_scale - 10e6*(time-250.)/435.;
      x += ParticleParameters::reflect_x_offset + 10e6*(time-250.)/435.;

      Vec3d pos(x,start_y,0);
      // Add a particle at this location, with bulk velocity at its starting point
      // TODO: Multiple
      //particles.push_back(Particle(PhysicalConstantsSI::mp, PhysicalConstantsSI::e, pos, V(pos)));

      /* Look up builk velocity in the V-field */
      Vec3d bulk_vel = V(pos);

      for(unsigned int i=0; i< ParticleParameters::num_particles; i++) {
         /* Create a particle with velocity drawn from the given distribution ... */
         Particle p = velocity_distribution->next_particle();
         /* Shift it by the bulk velocity ... */
         p.v += bulk_vel;
         /* And put it in place. */
         p.x=pos;
         particles.push_back(p);
      }

   }

   delete velocity_distribution;

   // Write out the state
   char filename_buffer[256];
   
   // /!\ Go one step opposite to the propagation time direction
   snprintf(filename_buffer,256, ParticleParameters::output_filename_pattern.c_str(),input_file_counter - ParticleParameters::propagation_direction);
   writeParticles(particles, filename_buffer);
}

void shockReflectivityScenario::afterPush(int step, double time, ParticleContainer& particles,
      Interpolated_Field& E, Interpolated_Field& B, Interpolated_Field& V, Interpolated_Field& R) {

   for(unsigned int i=0; i<particles.size(); i++) {

      if(isnan(vector_length(particles[i].x))) {
         // skip disabled particles
         continue;
      }

      //Get particle's y-coordinate
      double y = particles[i].x[1];

      // Get x for it's shock boundary (approx)
      double x = y / ParticleParameters::reflect_start_y;
      x*=-x;
      x *= ParticleParameters::reflect_y_scale - 10e6*(time-250.)/435.;
      x += ParticleParameters::reflect_x_offset + 10e6*(time-250.)/435.;

      // Boundaries are somewhat left or right of it
      double boundary_left = x - ParticleParameters::reflect_downstream_boundary;
      double boundary_right = x + ParticleParameters::reflect_upstream_boundary;

      // Check if the particle hit a boundary. If yes, mark it as disabled.
      // Original starting x of this particle
      int start_timestep = i / 200 / ParticleParameters::num_particles;
      double start_time = ParticleParameters::start_time + start_timestep * ParticleParameters::input_dt;
      if(particles[i].x[0] < boundary_left) {
         // Record it is transmitted.
         transmitted.addValue(Vec2d(y,start_time));

         // Disable by setting position to NaN and velocity to 0
         particles[i].x = Vec3d(std::numeric_limits<double>::quiet_NaN(),0.,0.);
         particles[i].v = Vec3d(0,0,0);
      } else if (particles[i].x[0] > boundary_right) {

         //Record it as reflected
         reflected.addValue(Vec2d(y,start_time));

         // Disable by setting position to NaN and velocity to 0
         particles[i].x = Vec3d(std::numeric_limits<double>::quiet_NaN(),0.,0.);
         particles[i].v = Vec3d(0.,0.,0.);
      }
   }
}

void shockReflectivityScenario::finalize(ParticleContainer& particles, Field& E, Field& B, Field& V, Field& R) {
   transmitted.save("transmitted.dat");
   transmitted.writeBovAscii("transmitted.dat.bov",0,"transmitted.dat");
   reflected.save("reflected.dat");
   reflected.writeBovAscii("reflected.dat.bov",0,"reflected.dat");
}


ParticleContainer ipShockScenario::initialParticles(Field& E, Field& B, Field& V, Field& R) {

  // Open output files for transmission and reflection
  traFile = fopen("transmitted.dat","w"); 
  refFile = fopen("reflected.dat","w"); 
  
   ParticleContainer particles;

   /* Prepare randomization engines */
   std::default_random_engine generator(ParticleParameters::random_seed);
   Distribution* velocity_distribution=ParticleParameters::distribution(generator);

   std::random_device rd;
   std::mt19937 gen(rd());
   std::uniform_real_distribution<> disx(ParticleParameters::ipshock_inject_x0, ParticleParameters::ipshock_inject_x1);
   std::uniform_real_distribution<> disy(ParticleParameters::ipshock_inject_y0, ParticleParameters::ipshock_inject_y1);
   std::uniform_real_distribution<> disz(ParticleParameters::ipshock_inject_z0, ParticleParameters::ipshock_inject_z1);

   /* Loop over particles to generate */
   for(unsigned int i=0; i< ParticleParameters::num_particles; i++) {
     /* Create a particle at a random position within the initialisation box */
     Real posx = disx(gen);
     Real posy = disy(gen);
     Real posz = disz(gen);
     Vec3d vpos(posx, posy, posz);

     /* Look up bulk velocity in the V-field */
     Vec3d bulk_vel = V(vpos); 
     
     /* Create a particle with velocity drawn from the given distribution ... */
     Particle p = velocity_distribution->next_particle();
     /* Shift it by the bulk velocity ... */
     p.v += bulk_vel;
     /* And put it in place. */
     p.x=vpos;
     particles.push_back(p);
   }

   delete velocity_distribution;
   return particles;
}

void ipShockScenario::newTimestep(int input_file_counter, int step, double time, ParticleContainer& particles,
      Interpolated_Field& E, Interpolated_Field& B, Interpolated_Field& V, Interpolated_Field& R) {

   char filename_buffer[256];

   // /!\ Go one step opposite to the propagation time direction
   snprintf(filename_buffer,256, ParticleParameters::output_filename_pattern.c_str(),input_file_counter - ParticleParameters::propagation_direction);
   writeParticles(particles, filename_buffer); //Generates VLSV file
}

void ipShockScenario::afterPush(int step, double time, ParticleContainer& particles,
      Interpolated_Field& E, Interpolated_Field& B, Interpolated_Field& V, Interpolated_Field& R) {
  
  /* Perform transmission / reflection check for each particle */
   for(unsigned int i=0; i<particles.size(); i++) {

      if(isnan(vector_length(particles[i].x))) {
         // skip disabled particles
         continue;
      }

      //Get particle's x-coordinate
      double x = particles[i].x[0];

      // Check if the particle hit a boundary. 
      // If yes, print it and mark it as disabled.
      if(particles[i].x[0] < ParticleParameters::ipshock_transmit) {
	// Record it as transmitted.
	//transmitted.addValue(Vec2d(y,start_time));
	
	// Write particle information to a file
	fprintf(traFile,"%12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e\n", time, 
		particles[i].x[0], particles[i].x[1], particles[i].x[2],
		particles[i].v[0], particles[i].v[1], particles[i].v[2],
		.5 * particles[i].m * dot_product(particles[i].v, particles[i].v) / PhysicalConstantsSI::e,
		dot_product(normalize_vector(particles[i].v), normalize_vector(B(particles[i].x))) );

	// Disable by setting position to NaN and velocity to 0
	particles[i].x = Vec3d(std::numeric_limits<double>::quiet_NaN(),0.,0.);
	particles[i].v = Vec3d(0,0,0);
      } else if (particles[i].x[0] > ParticleParameters::ipshock_reflect) {
	// Record it as reflected
	//reflected.addValue(Vec2d(y,start_time));
	
	// Write particle information to a file
	// Write particle information to a file
	fprintf(refFile,"%12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e\n", time, 
		particles[i].x[0], particles[i].x[1], particles[i].x[2],
		particles[i].v[0], particles[i].v[1], particles[i].v[2],
		.5 * particles[i].m * dot_product(particles[i].v, particles[i].v) / PhysicalConstantsSI::e,
		dot_product(normalize_vector(particles[i].v), normalize_vector(B(particles[i].x))) );

	// Disable by setting position to NaN and velocity to 0
	particles[i].x = Vec3d(std::numeric_limits<double>::quiet_NaN(),0.,0.);
	particles[i].v = Vec3d(0.,0.,0.);
      }
   }
   fflush(traFile);
   fflush(refFile);
}

void ipShockScenario::finalize(ParticleContainer& particles, Field& E, Field& B, Field& V, Field& R) {
   writeParticles(particles, "particles_final.vlsv");
   /* histograms */
   //transmitted.save("transmitted.dat");
   //transmitted.writeBovAscii("transmitted.dat.bov",0,"transmitted.dat");
   //reflected.save("reflected.dat");
   //reflected.writeBovAscii("reflected.dat.bov",0,"reflected.dat");

   fclose(traFile);
   fclose(refFile);
}







ParticleContainer InjectionScenario::initialParticles(Field& E, Field& B, Field& V, Field& R) {

   // Open output files for transmission and reflection
   initFile = fopen("inj_init.dat","w"); 
   traFile = fopen("inj_trans.dat","w"); 
   refFile = fopen("inj_refl.dat","w"); 
   lostFile = fopen("inj_lost.dat","w"); 

   meet_rho_File = fopen("inj_meet_rho.dat","w"); 
   meet_TNBS_File = fopen("inj_meet_TNBS.dat","w"); 
   meet_Mms_File = fopen("inj_meet_Mms.dat","w"); 
   meet_tbn_File = fopen("inj_meet_tbn.dat","w"); 
   meet_flipmu_File = fopen("inj_meet_flipmu.dat","w"); 

   minxFile = fopen("inj_minx.dat","w"); 
   boostFile = fopen("inj_boost.dat","w"); 

   ParticleContainer particles;

   /* Prepare randomization engines */
   std::default_random_engine generator(ParticleParameters::random_seed);
   Distribution* velocity_distribution=ParticleParameters::distribution(generator);

   std::random_device rd;
   std::mt19937 gen(rd());
   std::uniform_real_distribution<> disrad((M_PI/180.)*ParticleParameters::injection_start_deg0, (M_PI/180.)*ParticleParameters::injection_start_deg1);

   fs_hasmet_rho = new bool [ParticleParameters::num_particles];
   fs_hasmet_TNBS = new bool [ParticleParameters::num_particles];
   fs_hasmet_Mms = new bool [ParticleParameters::num_particles];
   fs_hasmet_tbn = new bool [ParticleParameters::num_particles];
   fs_hasmet_flipmu = new bool [ParticleParameters::num_particles];

   // Use solar wind frame from simulation data, or defined from config file?
   bool useVlsvV = true;
   Vec3d bulk_vel(ParticleParameters::injection_init_vx, ParticleParameters::injection_init_vy, ParticleParameters::injection_init_vz);
   if (sqrt( std::pow(ParticleParameters::injection_init_vx,2)+std::pow(ParticleParameters::injection_init_vy,2)
	     +std::pow(ParticleParameters::injection_init_vz,2) ) < 1e10) useVlsvV = false;

   /* Check for continuous initialisation */
   if (abs(ParticleParameters::injection_end_time) > 1e-10) {
     // Initialisation end time is set
     if (abs(ParticleParameters::injection_end_time) < abs(ParticleParameters::start_time)) {
       // If end time is invalid
       particlespertimestep = ParticleParameters::num_particles;
       std::cout << " Invalid init end time " << ParticleParameters::injection_end_time << std::endl;
       ParticleParameters::injection_end_time = 0;
       std::cout << " Initialising " << particlespertimestep << " all at time " << ParticleParameters::start_time << std::endl;
     } else {
       // Split particle count across time steps
       particlespertimestep = ParticleParameters::num_particles * ParticleParameters::input_dt 
	 / (abs(ParticleParameters::injection_end_time) - abs(ParticleParameters::start_time) );
       std::cout << " Initialising " << particlespertimestep << " at each of "
		 << (abs(ParticleParameters::injection_end_time) - abs(ParticleParameters::start_time))/ParticleParameters::input_dt
		 << " init steps" << std::endl;
       std::cout << "Between times " << ParticleParameters::start_time << " and " << ParticleParameters::injection_end_time << std::endl;
     }
   } else {
     // End time is not set, initialise all at once
     particlespertimestep = ParticleParameters::num_particles;
     std::cout << " Initialising " << particlespertimestep << " all at time " << ParticleParameters::start_time << std::endl;
   }

   /* Set memory addresses for all particles as inactive */
   for(unsigned int i=0; i< ParticleParameters::num_particles; i++) {
     fs_hasmet_rho[i]=false;
     fs_hasmet_TNBS[i]=false;
     fs_hasmet_Mms[i]=false;
     fs_hasmet_tbn[i]=false;
     fs_hasmet_flipmu[i]=false;

   }

//    /* Loop over particles to generate */
//    for(unsigned int i=0; i< particlespertimestep; i++) {
//      /* Create a particle with a random nose angle.
// 	Nose angle is defined as going from -pi to +pi, with negative values placed in the
// 	-y or -z hemisphere. */
//      Real rad = disrad(gen);
//      Real sinalpha = sin(rad);
//      Real cosalpha = cos(rad);

//      // At initialisation, use the first values for the bow shock fit
//      Real initr = ParticleParameters::injection_bs_p0
//        + ParticleParameters::injection_bs_p1 * std::pow(rad,1)
//        + ParticleParameters::injection_bs_p2 * std::pow(rad,2)
//        + ParticleParameters::injection_bs_p3 * std::pow(rad,3)
//        + ParticleParameters::injection_bs_p4 * std::pow(rad,4)
//        + ParticleParameters::injection_start_rplus;

//      Real posx;
//      Real posy;
//      Real posz;
//      if (B.dimension[1]->cells <= 1) {
//        /* polar x-z simulation */
//        posx = cosalpha * initr;
//        posy = 0.0;
//        posz = sinalpha * initr;
//      }
//      if (B.dimension[2]->cells <= 1) {
//        /* equatorial x-y simulation */
//        posx = cosalpha * initr;
//        posy = sinalpha * initr;
//        posz = 0.0;
//      }
//      Vec3d vpos(posx, posy, posz);

//      if (useVlsvV==true) {
//        /* Look up bulk velocity in the V-field */     
//        bulk_vel = V(vpos); 
//      }

//      /* Create a particle with velocity drawn from the given distribution ... */
//      Particle p = velocity_distribution->next_particle();
//      /* Shift it by the bulk velocity ... */
//      p.v += bulk_vel;
//      /* And put it in place. */
//      p.x = vpos;
//      Vec3d minx_init(1.e20, 0.0, 0.0);
//      p.minx_x = minx_init;

//      Real mu = dot_product(normalize_vector(p.v), normalize_vector(B(p.x)));
//      Real vsq = dot_product(p.v, p.v);

//      p.boost_vsq = vsq;
//      p.org_mu  = mu;

//      particles.push_back(p);

//      fprintf(initFile,"%d %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e\n", i, ParticleParameters::start_time, 
// 	     p.x[0], p.x[1], p.x[2], p.v[0], p.v[1], p.v[2], .5 * p.m * vsq / PhysicalConstantsSI::e, mu );

//    }

//    fflush(initFile);
//    fclose(initFile);

   delete velocity_distribution;
   return particles;
}

void InjectionScenario::newTimestep(int input_file_counter, int step, double time, ParticleContainer& particles,
      Interpolated_Field& E, Interpolated_Field& B, Interpolated_Field& V, Interpolated_Field& R) {


  if ((abs(ParticleParameters::injection_end_time) > 1.e-10) && (ParticleParameters::start_time+time<ParticleParameters::injection_end_time)) {
    std::cout << " Initialising " << particlespertimestep << " at time " << time << std::endl;
    std::cout << "       " << ParticleParameters::start_time << ParticleParameters::start_time +time << ParticleParameters::injection_end_time << std::endl;
    
    initFile = fopen("inj_init.dat","a"); 
    /* Prepare randomization engines */
    std::default_random_engine generator(ParticleParameters::random_seed);
    Distribution* velocity_distribution=ParticleParameters::distribution(generator);
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> disrad((M_PI/180.)*ParticleParameters::injection_start_deg0, (M_PI/180.)*ParticleParameters::injection_start_deg1);
    // Use solar wind frame from simulation data, or defined from config file?
    bool useVlsvV = true;
    Vec3d bulk_vel(ParticleParameters::injection_init_vx, ParticleParameters::injection_init_vy, ParticleParameters::injection_init_vz);
    if (sqrt( std::pow(ParticleParameters::injection_init_vx,2)+std::pow(ParticleParameters::injection_init_vy,2)
	      +std::pow(ParticleParameters::injection_init_vz,2) ) < 1e10) useVlsvV = false;

    /* Loop over particles to generate */
    for(unsigned int i=0; i< particlespertimestep; i++) {
      // Sanity check
      if (particles.size() >= ParticleParameters::num_particles) break;

      /* Create a particle with a random nose angle.
	 Nose angle is defined as going from -pi to +pi, with negative values placed in the
	 -y or -z hemisphere. */
      Real rad = disrad(gen);
      Real sinalpha = sin(rad);
      Real cosalpha = cos(rad);

      // Interpolate fit parameters between _bs_ and _bs2_ values
      Real interp_dist = time/(ParticleParameters::end_time - ParticleParameters::start_time);
      Real interp_bs_p0 = (1.-interp_dist)*ParticleParameters::injection_bs_p0 + interp_dist*ParticleParameters::injection_bs2_p0;
      Real interp_bs_p1 = (1.-interp_dist)*ParticleParameters::injection_bs_p1 + interp_dist*ParticleParameters::injection_bs2_p1;
      Real interp_bs_p2 = (1.-interp_dist)*ParticleParameters::injection_bs_p2 + interp_dist*ParticleParameters::injection_bs2_p2;
      Real interp_bs_p3 = (1.-interp_dist)*ParticleParameters::injection_bs_p3 + interp_dist*ParticleParameters::injection_bs2_p3;
      Real interp_bs_p4 = (1.-interp_dist)*ParticleParameters::injection_bs_p4 + interp_dist*ParticleParameters::injection_bs2_p4;

      Real initr = interp_bs_p0 
	+ interp_bs_p1 * std::pow(rad,1)
	+ interp_bs_p2 * std::pow(rad,2)
	+ interp_bs_p3 * std::pow(rad,3)
	+ interp_bs_p4 * std::pow(rad,4)
	+ ParticleParameters::injection_start_rplus;

      Real posx;
      Real posy;
      Real posz;
      if (B.a.dimension[1]->cells <= 1) {
	/* polar x-z simulation */
	posx = cosalpha * initr;
	posy = 0.0;
	posz = sinalpha * initr;
      }
      if (B.a.dimension[2]->cells <= 1) {
	/* equatorial x-y simulation */
	posx = cosalpha * initr;
	posy = sinalpha * initr;
	posz = 0.0;
      }
      Vec3d vpos(posx, posy, posz);

      if (useVlsvV==true) {
	/* Look up bulk velocity in the V-field */     
	bulk_vel = V(vpos); 
      }

      /* Create a particle with velocity drawn from the given distribution ... */
      Particle p = velocity_distribution->next_particle();
      /* Shift it by the bulk velocity ... */
      p.v += bulk_vel;
      /* And put it in place. */
      p.x = vpos;
      Vec3d minx_init(1.e20, 0.0, 0.0);
      p.minx_x = minx_init;

      Real mu = dot_product(normalize_vector(p.v), normalize_vector(B(p.x)));
      Real vsq = dot_product(p.v, p.v);

      p.boost_vsq = vsq;
      p.org_mu  = mu;
      particles.push_back(p);


      fprintf(initFile,"%d %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e\n", i, time, 
	      p.x[0], p.x[1], p.x[2], p.v[0], p.v[1], p.v[2], .5 * p.m * vsq / PhysicalConstantsSI::e, mu );

    }
    fflush(initFile);
    fclose(initFile);

    delete velocity_distribution;
  }
  
  char filename_buffer[256];
  snprintf(filename_buffer,256, ParticleParameters::output_filename_pattern.c_str(),input_file_counter-1);
  writeParticles(particles, filename_buffer); //Generates VLSV file
}

void InjectionScenario::beforePush(ParticleContainer& particles,
				   Interpolated_Field& E, Interpolated_Field& B, Interpolated_Field& V, Interpolated_Field& R) {
  /* Save pitch-angles for tracking purposes */
  for(unsigned int i=0; i<particles.size(); i++) {
    if(isnan(vector_length(particles[i].x))) {
      // skip disabled particles
      continue;
    }
    //particles[i].mu = dot_product(particles[i].v,Bval);
    //current_mu[i] = dot_product(normalize_vector(particles[i].v), normalize_vector(B(particles[i].x)));
  }
}

void InjectionScenario::afterPush(int step, double time, ParticleContainer& particles,
      Interpolated_Field& E, Interpolated_Field& B, Interpolated_Field& V, Interpolated_Field& R) {
  
  /* Perform transmission / reflection check for each particle */
   for(unsigned int i=0; i<particles.size(); i++) {

      if(isnan(vector_length(particles[i].x))) {
         // skip disabled particles
         continue;
      }

      //Get particle coordinates
      double x = particles[i].x[0];      
      double y = particles[i].x[1];
      double z = particles[i].x[2];
      
      // calculate particle costheta and bow shock position
      Real r = 0.0;
      Real rad = 0.0;
      if (B.a.dimension[1]->cells <= 1) {
	/* polar x-z simulation */
	rad = atan2(z,x);
	r = sqrt(x*x+z*z);
	//std::cerr<<"polar"<<std::endl;
      }
      if (B.a.dimension[2]->cells <= 1) {
	/* equatorial x-y simulation */
	rad = atan2(y,x);
	r = sqrt(x*x+y*y);
	//std::cerr<<"ecliptic"<<std::endl;
      }

      //Evaluate exact values at particle position via interpolation
      Vec3d Bval;
      Vec3d Eval;
      Vec3d Rval;
      Vec3d Vval;
      Bval = B(particles[i].x);
      Eval = E(particles[i].x);
      Rval = R(particles[i].x);
      Vval = V(particles[i].x);

      // Calculate current pitch-cosine and square of velocity
      Real curr_mu = dot_product(normalize_vector(particles[i].v), normalize_vector(Bval));
      Real curr_vsq = dot_product(particles[i].v, particles[i].v);
      
      // Interpolate fit parameters between _bs_ and _bs2_ values
      Real interp_dist = time/(ParticleParameters::end_time - ParticleParameters::start_time);
      Real interp_bs_p0 = (1.-interp_dist)*ParticleParameters::injection_bs_p0 + interp_dist*ParticleParameters::injection_bs2_p0;
      Real interp_bs_p1 = (1.-interp_dist)*ParticleParameters::injection_bs_p1 + interp_dist*ParticleParameters::injection_bs2_p1;
      Real interp_bs_p2 = (1.-interp_dist)*ParticleParameters::injection_bs_p2 + interp_dist*ParticleParameters::injection_bs2_p2;
      Real interp_bs_p3 = (1.-interp_dist)*ParticleParameters::injection_bs_p3 + interp_dist*ParticleParameters::injection_bs2_p3;
      Real interp_bs_p4 = (1.-interp_dist)*ParticleParameters::injection_bs_p4 + interp_dist*ParticleParameters::injection_bs2_p4;

      Real r0 = interp_bs_p0 
	+ interp_bs_p1 * std::pow(rad,1)
	+ interp_bs_p2 * std::pow(rad,2)
	+ interp_bs_p3 * std::pow(rad,3)
	+ interp_bs_p4 * std::pow(rad,4);

      // Check if particle is past rear boundary
      if (particles[i].x[0] < ParticleParameters::injection_x_bound_ds) {
	// Record it as lost
	fprintf(lostFile,"%d %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e\n", i, ParticleParameters::start_time+time, 
		particles[i].x[0], particles[i].x[1], particles[i].x[2],
		particles[i].v[0], particles[i].v[1], particles[i].v[2],
		.5 * particles[i].m * curr_vsq / PhysicalConstantsSI::e,
		curr_mu );

	// Disable by setting position to NaN and velocity to 0
	particles[i].x = Vec3d(std::numeric_limits<double>::quiet_NaN(),0.,0.);
	particles[i].v = Vec3d(0,0,0);
	continue;
      }

      // Check if particle escaped upstream
      if (r > r0 + ParticleParameters::injection_r_bound_us) {
	// Record it as reflected
	//reflected.addValue(Vec2d(y,start_time));
	
	// Write particle information to a file
	fprintf(refFile,"%d %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e\n",
		i, ParticleParameters::start_time+time, 
		particles[i].x[0], particles[i].x[1], particles[i].x[2],
		particles[i].v[0], particles[i].v[1], particles[i].v[2],
		.5 * particles[i].m * curr_vsq / PhysicalConstantsSI::e, curr_mu,
		Bval[0], Bval[1], Bval[2], Eval[0], Eval[1], Eval[2]);


	// Write information of particle at minimum x-position to file
	fprintf(minxFile,"%d %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e\n",
		i, ParticleParameters::start_time+particles[i].minx_t, 
		particles[i].minx_x[0], particles[i].minx_x[1], particles[i].minx_x[2],
		particles[i].minx_v[0], particles[i].minx_v[1], particles[i].minx_v[2],
		.5 * particles[i].m * dot_product(particles[i].minx_v, particles[i].minx_v) / PhysicalConstantsSI::e,
		particles[i].minx_mu,
		Bval[0], Bval[1], Bval[2], Eval[0], Eval[1], Eval[2]);

	// Disable by setting position to NaN and velocity to 0
	particles[i].x = Vec3d(std::numeric_limits<double>::quiet_NaN(),0.,0.);
	particles[i].v = Vec3d(0.,0.,0.);
	continue;
      }

      // Check if particle escaped downstream
      if (r < r0 - ParticleParameters::injection_r_bound_ds) {
	// Record it as transmitted.
	//transmitted.addValue(Vec2d(y,start_time));
	
	// Write particle information to a file
	fprintf(traFile,"%d %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e\n",
		i, ParticleParameters::start_time+time, 
		particles[i].x[0], particles[i].x[1], particles[i].x[2],
		particles[i].v[0], particles[i].v[1], particles[i].v[2],
		.5 * particles[i].m * curr_vsq / PhysicalConstantsSI::e, curr_mu,
		Bval[0], Bval[1], Bval[2], Eval[0], Eval[1], Eval[2]);

	// Disable by setting position to NaN and velocity to 0
	particles[i].x = Vec3d(std::numeric_limits<double>::quiet_NaN(),0.,0.);
	particles[i].v = Vec3d(0,0,0);
	continue;
      }

      // Store position and velocity at position with smallest x-coordinate
      if (x < particles[i].minx_x[0]) {
	particles[i].minx_t = time;
	particles[i].minx_mu = curr_mu;
	particles[i].minx_x = particles[i].x;
	particles[i].minx_v = particles[i].v;
      }

      // Check if particle energy has been doubled
      if (curr_vsq > 2. * particles[i].boost_vsq) {
	particles[i].boost_vsq = 2. * particles[i].boost_vsq;
	// Write particle information to a file
	fprintf(boostFile,"%d %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e\n",
		i, ParticleParameters::start_time+time, 
		particles[i].x[0], particles[i].x[1], particles[i].x[2],
		particles[i].v[0], particles[i].v[1], particles[i].v[2],
		.5 * particles[i].m * curr_vsq / PhysicalConstantsSI::e, curr_mu,
		Bval[0], Bval[1], Bval[2], Eval[0], Eval[1], Eval[2]);
      }

      // Check if particle has just now met shock for the first time (based on flip of mu)
      if ( (curr_mu * particles[i].org_mu < 0)  && (fs_hasmet_flipmu[i]==false) ) {
	fs_hasmet_flipmu[i]=true;
	fprintf(meet_flipmu_File,"%d %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e\n",
		i, ParticleParameters::start_time+time, 
		particles[i].x[0], particles[i].x[1], particles[i].x[2],
		particles[i].v[0], particles[i].v[1], particles[i].v[2],
		.5 * particles[i].m * curr_vsq / PhysicalConstantsSI::e, curr_mu, 
		Bval[0], Bval[1], Bval[2], Eval[0], Eval[1], Eval[2]);
      }

      // Check if particle has just now met shock for the first time (based on density)
      // Rval[0] is total density for this position
      if ( (Rval[0] > ParticleParameters::injection_rho_meet) && (fs_hasmet_rho[i]==false) ) {
	fs_hasmet_rho[i]=true;
	fprintf(meet_rho_File,"%d %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e\n",
		i, ParticleParameters::start_time+time, 
		particles[i].x[0], particles[i].x[1], particles[i].x[2],
		particles[i].v[0], particles[i].v[1], particles[i].v[2],
		.5 * particles[i].m * curr_vsq / PhysicalConstantsSI::e, curr_mu, 
		Bval[0], Bval[1], Bval[2], Eval[0], Eval[1], Eval[2]);
      }
      
      // Check if particle has just now met shock for the first time (based on non-backstreaming temperature)
      // Rval[1] is non-backstreaming temperature for this position
      if ( (Rval[1] > ParticleParameters::injection_TNBS_meet) && (fs_hasmet_TNBS[i]==false) ) {
	fs_hasmet_TNBS[i]=true;
	fprintf(meet_TNBS_File,"%d %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e\n",
		i, ParticleParameters::start_time+time, 
		particles[i].x[0], particles[i].x[1], particles[i].x[2],
		particles[i].v[0], particles[i].v[1], particles[i].v[2],
		.5 * particles[i].m * curr_vsq / PhysicalConstantsSI::e, curr_mu, 
		Bval[0], Bval[1], Bval[2], Eval[0], Eval[1], Eval[2]);
      }

      /* As per https://www.math24.net/tangent-normal-lines/
	 calculate the direction of the shock-normal based on the shock shape fit
	 dy/dx = -(r' cos theta - r sin theta)/( r' sin theta + r cos theta) */
      Real rder = interp_bs_p1
	+ 2. * interp_bs_p2 * std::pow(rad,1)
	+ 3. * interp_bs_p3 * std::pow(rad,2)
	+ 4. * interp_bs_p4 * std::pow(rad,3);
      Real costheta = cos(rad);
      Real sintheta = sin(rad);
      Vec3d normalvect((rder*sintheta + r0*costheta), -(rder*costheta - r0*sintheta), 0);

      // Calculate shock-normal bulk velocity from shock shape fit
      // Rval [2] is magnetosonic speed v_ms for this position. Calculate Mms assuming shock normal direction from fit.
      Real bulkV_normal = abs(dot_product(normalize_vector(normalvect), Vval));
      Real Mms_local = bulkV_normal / Rval[2];

      // Check if particle has just now met shock for the first time (based on Mms)
      if ( (Mms_local < ParticleParameters::injection_Mms_meet) && (fs_hasmet_Mms[i]==false) ) {
	fs_hasmet_Mms[i]=true;
	fprintf(meet_Mms_File,"%d %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e\n",
		i, ParticleParameters::start_time+time, 
		particles[i].x[0], particles[i].x[1], particles[i].x[2],
		particles[i].v[0], particles[i].v[1], particles[i].v[2],
		.5 * particles[i].m * curr_vsq / PhysicalConstantsSI::e, curr_mu, 
		Bval[0], Bval[1], Bval[2], Eval[0], Eval[1], Eval[2]);
      }

      // Check if particle has just now met shock for the first time (based on thetaBn)
      // Calculate theta_Bn for current position from magnetic field and shock normal direction from fit
      Real thetabn_local = acos(abs(dot_product(normalize_vector(normalvect), normalize_vector(Bval))))*(180./3.14159);
      if ( (thetabn_local > ParticleParameters::injection_tbn_meet) && (fs_hasmet_tbn[i]==false) ) {
	fs_hasmet_tbn[i]=true;
	fprintf(meet_tbn_File,"%d %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e\n",
		i, ParticleParameters::start_time+time, 
		particles[i].x[0], particles[i].x[1], particles[i].x[2],
		particles[i].v[0], particles[i].v[1], particles[i].v[2],
		.5 * particles[i].m * curr_vsq / PhysicalConstantsSI::e, curr_mu, 
		Bval[0], Bval[1], Bval[2], Eval[0], Eval[1], Eval[2]);
      }

   }
   fflush(traFile);
   fflush(refFile);
   fflush(lostFile);

   fflush(meet_rho_File);
   fflush(meet_TNBS_File);
   fflush(meet_Mms_File);
   fflush(meet_tbn_File);
   fflush(meet_flipmu_File);

   fflush(minxFile);
   fflush(boostFile);
}

void InjectionScenario::finalize(ParticleContainer& particles, Field& E, Field& B, Field& V, Field& R) {
   writeParticles(particles, "particles_final.vlsv");
   /* histograms */
   //transmitted.save("transmitted.dat");
   //transmitted.writeBovAscii("transmitted.dat.bov",0,"transmitted.dat");
   //reflected.save("reflected.dat");
   //reflected.writeBovAscii("reflected.dat.bov",0,"reflected.dat");

   fclose(traFile);
   fclose(refFile);
   fclose(lostFile);

   fclose(meet_rho_File);
   fclose(meet_TNBS_File);
   fclose(meet_Mms_File);
   fclose(meet_tbn_File);
   fclose(meet_flipmu_File);

   fclose(minxFile);
   fclose(boostFile);
}



Scenario* createScenario(std::string name) {
   std::map<std::string, Scenario*(*)()> scenario_lookup;
   scenario_lookup["single"]=&createScenario<singleParticleScenario>;
   scenario_lookup["distribution"]=&createScenario<distributionScenario>;
   scenario_lookup["precipitation"]=&createScenario<precipitationScenario>;
   scenario_lookup["analysator"]=&createScenario<analysatorScenario>;
   scenario_lookup["reflectivity"]=&createScenario<shockReflectivityScenario>;
   scenario_lookup["ipshock"]=&createScenario<ipShockScenario>;
   scenario_lookup["injection"]=&createScenario<InjectionScenario>;

   if(scenario_lookup.find(name) == scenario_lookup.end()) {
      std::cerr << "Error: can't find particle pusher mode \"" << name << "\". Aborting." << std::endl;
      exit(0);
   }

   return scenario_lookup[name]();
}




