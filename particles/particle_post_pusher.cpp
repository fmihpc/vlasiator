/*
 * Postprocessing particle trajectory analyzer
 * for Vlasiator
 *
 */
#include <mpi.h>
#include <iostream>
#include <random>
#include "particles.h"
#include "physconst.h"
#include "readfields.h"

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

	/* Parse commandline */
	if(argc < 3) {
		help_message();
		return 1;
	}

	/* Read fields from specified input file */
	std::string filename = argv[1];
	Field E,B,V;
	if(checkVersion(filename)) {
		readfields<newVlsv::Reader>(filename,E,B,V);
	} else {
		readfields<oldVlsv::Reader>(filename,E,B,V);
	}

	//debug_output(B, "B.png");


	/* Init particles */
	std::vector<Particle> particles;
	pusher_mode mode;

	double dt=0.0004012841091492777/10;
	double maxtime=100;
	int maxsteps = maxtime/dt;

	if(!strcmp(argv[2],"single")) {

		if(argc != 9) {
			std::cerr<< "Mode 'single' requires 6 additional arguments!" << std::endl;
			help_message();
			return 1;
		}

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
	} else if(!strcmp(argv[2], "distribution")) {

		if(argc != 8) {
			std::cerr<< "Mode 'distribution' requires 5 additional arguments!" << std::endl;
			help_message();
			return 1;
		}

		double pos[3];
		double temperature;
		unsigned int num_particles;

		pos[0] = strtod(argv[3],NULL);
		pos[1] = strtod(argv[4],NULL);
		pos[2] = strtod(argv[5],NULL);
		temperature = strtod(argv[6],NULL);
		num_particles = strtoul(argv[7],NULL,0);

		Vec3d vpos;
		vpos.load(pos);

		mode = distribution;

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
		std::cerr << "Unknown operation mode '" << argv[2] << "'." << std::endl;
		help_message();
		return 1;
	}

	/* Dump the initial state */
	if(mode == distribution) {
		write_particles(particles, "particles_initial.vlsv");
	}

	std::cerr << "Pushing " << particles.size() << " particles for " << maxsteps << " steps..." << std::endl;
        std::cerr << "[                                                                        ]\x0d[";

	char output_filename[256];

	/* Push them around */
	for(int step=0; step<maxsteps; step++) {

		if(step%(10000) == 0) {
			snprintf(output_filename,256, "particles_%06i.vlsv",step);
			write_particles(particles, output_filename);
		}

		#pragma omp parallel for
		for(unsigned int i=0; i< particles.size(); i++) {
			/* Get E- and B-Field at their position */
			Vec3d Eval,Bval;

			Eval = E(particles[i].x);
			Bval = B(particles[i].x);

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
