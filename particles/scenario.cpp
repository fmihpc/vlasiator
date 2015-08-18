#include <iostream>
#include "scenario.h"

std::vector<Particle> singleParticleScenario::initial_particles(Field& E, Field& B, Field& V) {

  std::vector<Particle> particles;

  Vec3d vpos(ParticleParameters::init_x, ParticleParameters::init_y, ParticleParameters::init_z);
  /* Look up builk velocity in the V-field */
  Vec3d bulk_vel = V(vpos);

  particles.push_back(Particle(PhysicalConstantsSI::mp, PhysicalConstantsSI::e, vpos, bulk_vel));

  return particles;
}

void singleParticleScenario::after_push(int step, double time, std::vector<Particle>& particles, Field& E, Field& B, Field& V) {

  Vec3d& x = particles[0].x;
  Vec3d& v = particles[0].v;

  std::cout << 0 << " " << time << "\t" <<  x[0] << " " << x[1] << " " << x[2] << "\t" << v[0] << " " << v[1] << " " << v[2] << std::endl;
}



std::vector<Particle> distributionScenario::initial_particles(Field& E, Field& B, Field& V) {

  std::vector<Particle> particles;

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

  return particles;
}

void distributionScenario::new_timestep(int input_file_counter, int step, double time, std::vector<Particle>& particles, Field& E, Field& B, Field& V) {

  char filename_buffer[256];

  snprintf(filename_buffer,256, ParticleParameters::output_filename_pattern.c_str(),input_file_counter-1);
  write_particles(particles, filename_buffer);
}

void distributionScenario::finalize(std::vector<Particle>& particles, Field& E, Field& B, Field& V) {
  write_particles(particles, "particles_final.vlsv");
}

void precipitationScenario::after_push(int step, double time, std::vector<Particle>& particles, Field& E, Field& B, Field& V) {

  for(unsigned int i=0; i<particles.size(); i++) {

    if(vector_length(particles[i].x) == 0) {
      // skip disabled particles
      continue;
    }

    // Check if the particle hit a boundary. If yes, mark it as disabled.
    // Original starting x of this particle
    double start_pos = ParticleParameters::precip_start_x + ((double)(i%ParticleParameters::num_particles))/ParticleParameters::num_particles * (ParticleParameters::precip_stop_x - ParticleParameters::precip_start_x);
    int start_timestep = i / ParticleParameters::num_particles;
    if(vector_length(particles[i].x) <= ParticleParameters::precip_inner_boundary) {

      // Record latitude and energy
      double latitude = atan2(particles[i].x[2],particles[i].x[0]);
      printf("%u %i %lf %lf %lf\n",i, start_timestep, start_pos, latitude, .5*particles[i].m * dot_product(particles[i].v,particles[i].v)/PhysicalConstantsSI::e);

      // Disable by setting position and velocity to 0
      particles[i].x = Vec3d(0,0,0);
      particles[i].v = Vec3d(0,0,0);
    } else if (particles[i].x[0] <= ParticleParameters::precip_start_x) {

      // Record marker value for lost particle
      printf("%u %i %lf -5. -1.\n", i, start_timestep, start_pos);
      // Disable by setting position and velocity to 0
      particles[i].x = Vec3d(0,0,0);
      particles[i].v = Vec3d(0,0,0);
    }
  }
}

void precipitationScenario::new_timestep(int input_file_counter, int step, double time, std::vector<Particle>& particles, Field& E, Field& B, Field& V) {

  // Create particles along the negative x-axis, from inner boundary
  // up to outer one
  for(unsigned int i=0; i< ParticleParameters::num_particles; i++) {

    // Choose x coordinate
    double start_x = ParticleParameters::precip_start_x + ((double)i)/ParticleParameters::num_particles * (ParticleParameters::precip_stop_x - ParticleParameters::precip_start_x);
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

  snprintf(filename_buffer,256, ParticleParameters::output_filename_pattern.c_str(),input_file_counter-1);
  write_particles(particles, filename_buffer);
}


std::vector<Particle> analysatorScenario::initial_particles(Field& E, Field& B, Field& V) {

  std::vector<Particle> particles;

  std::cerr << "Reading initial particle data from stdin" << std::endl
    << "(format: x y z vx vy vz)" << std::endl;

  while(std::cin) {
    double x0,x1,x2,v0,v1,v2;
    std::cin >> x0 >> x1 >> x2 >> v0 >> v1 >> v2;
    if(std::cin) {
      particles.push_back(Particle(PhysicalConstantsSI::mp, PhysicalConstantsSI::e, Vec3d(x0,x1,x2), Vec3d(v0,v1,v2)));
    }
  }

  return particles;
}

void analysatorScenario::new_timestep(int input_file_counter, int step, double time, std::vector<Particle>& particles, Field& E, Field& B, Field& V) {

  for(unsigned int i=0; i< particles.size(); i++) {
    Vec3d& x = particles[i].x;
    Vec3d& v = particles[i].v;
    std::cout << i << " " << time << "\t" <<  x[0] << " " << x[1] << " " << x[2] << "\t" << v[0] << " " << v[1] << " " << v[2] << std::endl;
  }
}

Scenario* create_scenario(std::string name) {
  std::map<std::string, Scenario*(*)()> scenario_lookup;
  scenario_lookup["single"]=&create_scenario<singleParticleScenario>;
  scenario_lookup["distribution"]=&create_scenario<distributionScenario>;
  scenario_lookup["precipitation"]=&create_scenario<precipitationScenario>;
  scenario_lookup["analysator"]=&create_scenario<analysatorScenario>;

  if(scenario_lookup.find(name) == scenario_lookup.end()) {
    std::cerr << "Error: can't find particle pusher mode \"" << name << "\". Aborting." << std::endl;
    exit(0);
  }

  return scenario_lookup[name]();
}
