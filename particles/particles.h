#pragma once
#include <vector>
#include "vectorclass.h"
#include "vector3d.h"
#include "../definitions.h"

struct Particle {
      Real m;
      Real q;
      Vec3d x;
      Vec3d v;

      Particle(Real mass, Real charge, Vec3d _x, Vec3d _v) :
         m(mass),q(charge),x(_x),v(_v) {}

      /* Particle propagation given E- and B-Field at the particle location
       * with the Boris-Method */
      void push(Vec3d& B, Vec3d& E, double dt);
};


void write_particles(std::vector<Particle>& p, const char* filename);
