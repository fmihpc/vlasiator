#pragma once

#include <random>
#include "particles.h"
#include "particleparameters.h"
#include "physconst.h"

/* virtual base of particle distributions */
class Distribution {
	protected:
		Real mass,charge;
		std::default_random_engine& rand;
	public:
		Distribution(std::default_random_engine& _rand, Real _mass, Real _charge) : mass(_mass), charge(_charge),rand(_rand) {};
		Distribution(std::default_random_engine& _rand);
		virtual Particle next_particle() = 0;
};


/* Thermal maxwell-boltzmann distribution */
class Maxwell_Boltzmann : public Distribution {
	private:
		std::normal_distribution<Real> velocity_distribution;

	public:
		/* Constructor where all parameters are explicitly given */
		Maxwell_Boltzmann(std::default_random_engine& _rand,Real _mass, Real _charge, Real kT) : 
			Distribution(_rand,_mass,_charge),
			velocity_distribution(0,sqrt(kT/mass)) {};

		/* Constructor that fishes parameters from the parameter-class by itself */
		Maxwell_Boltzmann(std::default_random_engine& _rand);
		virtual Particle next_particle();
};

/* Monoenergetic isotropic particles */
class Monoenergetic : public Distribution {
	private:

		std::uniform_real_distribution<Real> direction_random;
		Real vel;

	public:
		Monoenergetic(std::default_random_engine& _rand,Real _mass, Real _charge, Real _vel) :
			Distribution(_rand, _mass, _charge),
			vel(_vel) {
		}
		Monoenergetic(std::default_random_engine& _rand);

		virtual Particle next_particle() {
			/* Sphere point-picking to get isotropic direction (from wolfram Mathworld) */
			Real u,v;
			u = 2.*direction_random(rand)-1.;
			v = direction_random(rand) * 2. * M_PI;

			Vec3d dir(sqrt(1-u*u) * cos(v),
					sqrt(1-u*u) * sin(v),
					u);
			return Particle(mass, charge, Vec3d(0.), vel*dir);
		}
};

template<typename T> Distribution* create_distribution(std::default_random_engine& rand) {
	return new T(rand);
}
