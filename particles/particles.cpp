#include "particles.h"
#include "physconst.h"
#include "relativistic_math.h"
#include "vectorclass.h"
#include "vector3d.h"

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
