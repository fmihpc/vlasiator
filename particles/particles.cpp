#include "particles.h"
#include "physconst.h"
#include "relativistic_math.h"
#include <glm/glm.hpp>

using namespace glm;

/* Particle propagation given E- and B-Field at the particle location
 * with the Boris-Method */
void Particle::push(dvec3& B, dvec3& E, double dt) {

	dvec3 uminus = v + (q * E * dt)/(2. * m);
	dvec3 h = (q * B * dt)/(2. * m * gamma(uminus));
	dvec3 uprime = uminus + cross(uminus, h);
	h = (2.* h)/(1. + dot(h,h));
	dvec3 uplus = uminus + cross(uprime, h);

	v = uplus + (q * E * dt)/(2. * m);
	x += dt * v;
}
