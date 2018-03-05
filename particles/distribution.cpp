/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://www.physics.helsinki.fi/vlasiator/
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */
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

double Kappa::find_v_for_r(double rand) {

   /* Use binary search to find this value */
   int step = lookup_size/2;
   int a=lookup_size/2;
   double diff;
   do {
      step = ceil(step*.5);

      diff = rand-lookup[a];
      if(diff < 0) {
         a-=step;
      } else {
         a+=step;
      }
      a = (a<lookup_size-1)?(a):(lookup_size-1);
   } while(step > 1);

   /* Final step adjustment */
   if(rand < lookup[a]) {
      a--;
   }

   /* Special case for the upper bound: Extrapolate the last bin */
   if(a+1 > lookup_size-1) {
      return maxw0/lookup_size*
         (lookup_size-1 + (rand-lookup[lookup_size-2])/(lookup[lookup_size-1]-lookup[lookup_size-1]));
   }
   /* Linear interpolate between the two neighbouring entries */
   return maxw0/lookup_size*(a + (rand-lookup[a])/(lookup[a+1]-lookup[a]));
}

Kappa6::Kappa6(std::default_random_engine& _rand) : Kappa(_rand) {
   Real kT = ParticleParameters::temperature * PhysicalConstantsSI::k;
   w0 = sqrt(2.* kT *(6. - 1.5)/(6. * mass));
   generate_lookup();
}
Kappa2::Kappa2(std::default_random_engine& _rand) : Kappa(_rand) {
   Real kT = ParticleParameters::temperature * PhysicalConstantsSI::k;
   w0 = sqrt(2.* kT *(2. - 1.5)/(2. * mass));
   generate_lookup();
}


Particle Maxwell_Boltzmann::next_particle() {
   Vec3d v(velocity_distribution(rand), velocity_distribution(rand), velocity_distribution(rand));

   return Particle(mass, charge, Vec3d(0.), v);
}
