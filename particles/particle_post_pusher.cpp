/*
 * Postprocessing particle trajectory analyzer
 * for Vlasiator
 *
 */
#include <iostream>
#include "particles.h"
#include "physconst.h"
#include "readfields.h"

int main(int argc, char** argv) {

	/* Parse commandline */
	if(argc < 2) {
		std::cerr << "Usage: particle_post_pusher <filename>" << std::endl;
		return 1;
	}

	std::string filename = argv[1];

	Field E,B;

	if(checkVersion(filename)) {
			readfields<newVlsv::Reader>(filename,E,B);
	} else {
			readfields<oldVlsv::Reader>(filename,E,B);
	}

	debug_output(B,"B.png");

	double dt=0.0293778 / 10.;
	int num_particles = 1;
	int maxsteps = 10000;

	for(int i=0; i< num_particles; i++) {
		/* Initialize particles */
		Particle p(PhysicalConstantsSI::mp, PhysicalConstantsSI::e, glm::dvec3(1e6,-1e6,0), glm::dvec3(-1e5,0,0));

		for(int step=0; step<maxsteps; step++) {
			/* Get E- and B-Field at their position */
			glm::dvec3 Eval,Bval;

			Eval = E(p.x);
			Bval = B(p.x);

			/* Push them around */
			p.push(Bval,Eval,dt);

			/* Output position and velocities */
			//std::cout << "step " << step << ": B=[" << Bval.x << ", " << Bval.y << ", " << Bval.z << "], E=["
			//	<< Eval.x << ", " << Eval.y << ", " << Eval.z << "]" << std::endl;
			std::cout << i << " " << step << " " << p.x.x << " " << p.x.y << " " << p.x.z << std::endl;
			//std::cerr << ".";
		}
	}

	std::cerr << std::endl;
	return 0;
}
