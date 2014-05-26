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
#include "particleparameters.h"
#include "../readparameters.h"

// "mode" that the particle pusher operates in
enum pusher_mode {
	single = 0,    // Single particle trajectory tracking
	distribution  // Particle distribution created at a given point
};

void help_message() {
	std::cerr << "Usage: particle_post_pusher <filename> <mode> [args]" << std::endl;
	std::cerr << "Where <filename> is a .vlsv file providing field input and <mode> is one of:" << std::endl;
	std::cerr << std::endl;
	std::cerr  << "  single <x> <y> <z> <vx> <vy> <vz>" << std::endl;
	std::cerr  << "  distribution <x> <y> <z> <temperature> <num_particles>" << std::endl;
}

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

	/* TODO: This assumes that output is written every second. Beware. */
	int input_file_counter=floor(ParticleParameters::start_time)+1;
	Field E[2],B[2],V;
	snprintf(filename_buffer,256,filename_pattern.c_str(),floor(ParticleParameters::start_time));
	if(checkVersion(filename_buffer)) {
		readfields<newVlsv::Reader>(filename_buffer,E[1],B[1],V);
	} else {
		readfields<oldVlsv::Reader>(filename_buffer,E[1],B[1],V);
	}

	/* Init particles */
	std::vector<Particle> particles;
	pusher_mode mode;

	double dt=0.0004012841091492777/10;
	double maxtime=ParticleParameters::end_time - ParticleParameters::start_time;
	int maxsteps = maxtime/dt;

	if(ParticleParameters::mode == "single") {

		std::cerr << "Sorry, single-particle mode is currently broken, because urs was to lazy to implement it yet!" << std::endl;
		return 1;

		double pos[3], vel[3];
		pos[0] = strtod(argv[3],NULL);
		pos[1] = strtod(argv[4],NULL);
		pos[2] = strtod(argv[5],NULL);
		vel[0] = strtod(argv[6],NULL);
		vel[1] = strtod(argv[7],NULL);
		vel[2] = strtod(argv[8],NULL);

		Vec3d Vpos, Vvel;
		Vpos.load(pos);
		Vvel.load(vel);

		mode = single;
		particles.push_back(Particle(PhysicalConstantsSI::mp, PhysicalConstantsSI::e, Vpos, Vvel));
	} else if(ParticleParameters::mode == "distribution") {

		double temperature;
		uint64_t num_particles;

		temperature = ParticleParameters::temperature;
		num_particles = ParticleParameters::num_particles;

		Vec3d vpos(ParticleParameters::init_x, ParticleParameters::init_y, ParticleParameters::init_z);
		mode = distribution;

		/* Look up builk velocity in the V-field */
		Vec3d bulk_vel = V(vpos);

		std::normal_distribution<Real> velocity_distribution(0,sqrt(temperature*PhysicalConstantsSI::k/PhysicalConstantsSI::mp));
		std::default_random_engine generator;

		for(unsigned int i=0; i< num_particles; i++) {
			// Generate normal-distribution...
			Vec3d vel(velocity_distribution(generator),velocity_distribution(generator),velocity_distribution(generator));
			// .. and shift it by the bulk velocity.
			vel += bulk_vel;
			particles.push_back(Particle(PhysicalConstantsSI::mp, PhysicalConstantsSI::e, vpos, vel));
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

		/* Load newer fields, if neccessary */
		while(ParticleParameters::start_time + step*dt >= E[1].time) {
			E[0]=E[1];
			B[0]=B[1];

			/* TODO: Don't reload V here! */
			snprintf(filename_buffer,256,filename_pattern.c_str(),input_file_counter++);
			if(checkVersion(filename_buffer)) {
				readfields<newVlsv::Reader>(filename_buffer,E[1],B[1],V);
			} else {
				readfields<oldVlsv::Reader>(filename_buffer,E[1],B[1],V);
			}

			/* This is also a good opportunity to write some output.
			 * Note the -2 here, since we're writing this as we have just
			 * passed that given timestep, and incremented it by 1. */
			snprintf(filename_buffer,256, ParticleParameters::output_filename_pattern.c_str(),input_file_counter-2);
			write_particles(particles, filename_buffer);
		}

		Interpolated_Field cur_E(E[0],E[1],step*dt);
		Interpolated_Field cur_B(B[0],B[1],step*dt);

		#pragma omp parallel for
		for(unsigned int i=0; i< particles.size(); i++) {
			/* Get E- and B-Field at their position */
			Vec3d Eval,Bval;

			Eval = cur_E(particles[i].x);
			Bval = cur_B(particles[i].x);

			/* Push them around */
			particles[i].push(Bval,Eval,dt);
		}

		/* Draw progress bar */
		if((step % (maxsteps/72))==0) {
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
