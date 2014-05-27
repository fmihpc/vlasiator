#include <random>
#include "vector3d.h"
#include "distribution.h"
#include "particles.h"

Distribution::Distribution(std::default_random_engine& _rand) : rand(_rand) {
	mass = ParticleParameters::mass;
	charge = ParticleParameters::charge;
}
Maxwell_Boltzmann::Maxwell_Boltzmann(std::default_random_engine& _rand) : Distribution(_rand) {
	Real kT = ParticleParameters::temperature * PhysicalConstantsSI::k;
	velocity_distribution=std::normal_distribution<Real>(0.,sqrt(kT/mass));
}
Monoenergetic::Monoenergetic(std::default_random_engine& _rand) : Distribution(_rand) {
	vel = ParticleParameters::particle_vel;
}

Particle Maxwell_Boltzmann::next_particle() {
	Vec3d v(velocity_distribution(rand), velocity_distribution(rand), velocity_distribution(rand));

	return Particle(mass, charge, Vec3d(0.), v);
}
