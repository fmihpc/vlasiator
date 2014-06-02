#include <vector>
#include "particles.h"
#include "physconst.h"
#include "relativistic_math.h"
#include "vectorclass.h"
#include "vector3d.h"
#include "vlsv_writer.h"


/* Particle propagation given E- and B-Field at the particle location
 * with the Boris-Method */
void Particle::push(Vec3d& B, Vec3d& E, double dt) {

	Vec3d uminus = v + (q * E * dt)/(2. * m);
	Vec3d h = (q * B * dt)/(2. * m * gamma(uminus));
	Vec3d uprime = uminus + cross_product(uminus, h);
	h = (2.* h)/(1. + dot_product(h,h));
	Vec3d uplus = uminus + cross_product(uprime, h);

	v = uplus + (q * E * dt)/(2. * m);
	x += dt * v;
}

void write_particles(std::vector<Particle>& p,const char* filename) {

	vlsv::Writer vlsvWriter;
	vlsvWriter.open(filename,MPI_COMM_WORLD,0);

	std::vector<double> writebuf(p.size() * 3);

	/* First, store particle positions */
	for(unsigned int i=0; i < p.size(); i++) {
		p[i].x.store(&(writebuf[3*i]));
	}

	std::map<std::string,std::string> attribs;
	attribs["name"] = "proton_position";
	attribs["type"] = vlsv::mesh::STRING_POINT;
	if (vlsvWriter.writeArray("MESH",attribs,p.size(),3,writebuf.data()) == false) {
		std::cerr << "\t ERROR failed to write particle positions!" << std::endl;
	}

	/* Then, velocities */
	for(unsigned int i=0; i < p.size(); i++) {
		p[i].v.store(&(writebuf[3*i]));
	}

	attribs["name"] = "proton_velocity";
	if (vlsvWriter.writeArray("MESH",attribs,p.size(),3,writebuf.data()) == false) {
		std::cerr << "\t ERROR failed to write particle velocities!" << std::endl;
	}
	vlsvWriter.close();
}
