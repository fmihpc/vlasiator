#pragma once
#include <glm/glm.hpp>
#include "../definitions.h"

struct Particle {
		Real m;
		Real q;
		glm::dvec3 x;
		glm::dvec3 v;

		Particle(Real mass, Real charge, glm::dvec3 _x, glm::dvec3 _v) :
			m(mass),q(charge),x(_x),v(_v) {}

		/* Particle propagation given E- and B-Field at the particle location
		 * with the Boris-Method */
		void push(glm::dvec3& B, glm::dvec3& E, double dt);
};


