#pragma once
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
#include "particles.h"
#include "particleparameters.h"
#include "physconst.h"

/* virtual base of particle distributions */
class Distribution
{
   public:
      Distribution(std::default_random_engine& _rand, Real _mass, Real _charge)
         : mass(_mass), charge(_charge),rand(_rand) {};
      Distribution(std::default_random_engine& _rand);
      virtual Particle next_particle() = 0;
   protected:
      Real mass,charge;
      std::default_random_engine& rand;
};


/* Thermal maxwell-boltzmann distribution */
class Maxwell_Boltzmann : public Distribution
{
   public:
      /* Constructor where all parameters are explicitly given */
      Maxwell_Boltzmann(std::default_random_engine& _rand,Real _mass, Real _charge, Real kT) :
         Distribution(_rand,_mass,_charge),
         velocity_distribution(0,sqrt(kT/mass)) {};

      /* Constructor that fishes parameters from the parameter-class by itself */
      Maxwell_Boltzmann(std::default_random_engine& _rand);
      virtual Particle next_particle();
   private:
      std::normal_distribution<Real> velocity_distribution;
};

/* Monoenergetic isotropic particles */
class Monoenergetic : public Distribution
{
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
   private:

      std::uniform_real_distribution<Real> direction_random;
      Real vel;

};

/* Kappa distributions as taken from the CSA code */
class Kappa : public Distribution
{
   public:
      Kappa(std::default_random_engine& _rand, Real _mass, Real _charge, Real _w0, Real _maxw0)
         : Distribution(_rand,_mass,_charge),w0(_w0),maxw0(_maxw0) {}

      Kappa(std::default_random_engine& _rand) : Distribution(_rand) {
         maxw0=50;
      }

      virtual Particle next_particle() {

         Real vel = w0 * find_v_for_r(r(rand));

         /* Sphere point-picking to get isotropic direction (from wolfram Mathworld) */
         Real u,v;
         u = 2.*r(rand)-1.;
         v = r(rand) * 2. * M_PI;

         Vec3d dir(sqrt(1-u*u) * cos(v),
               sqrt(1-u*u) * sin(v),
               u);
         return Particle(mass, charge, Vec3d(0.), vel*dir);
      }

   protected:
      std::uniform_real_distribution<Real> r;
      std::vector<Real> lookup;
      const static int lookup_size = 65536;
      Real w0;
      Real maxw0;

      virtual void generate_lookup() = 0;

      Real find_v_for_r(Real rand);

};

/* Kappa = 6 version */
class Kappa6 : public Kappa
{

   public:
      Kappa6(std::default_random_engine& _rand, double _mass, double _charge, double _w0, double _maxw0)
         : Kappa(_rand,_mass,_charge,_w0,_maxw0){
            generate_lookup();
      }
      Kappa6(std::default_random_engine& _rand);

   private:
      virtual void generate_lookup() {
         double prefix = (512. * sqrt(2./3.))/(63.*M_PI);

         for(int i=0; i<lookup_size; i++) {
            double vkappa = (maxw0/lookup_size)*i;
            double under = vkappa*vkappa + 6.;

            lookup.push_back( prefix*(-23328.*vkappa/(pow(under,6.))
                     + 1944. * vkappa/(5. *pow(under,5.))
                     + 729.  * vkappa/(10.*pow(under,4.))
                     + 567.  * vkappa/(40.*pow(under,3.))
                     + 189.  * vkappa/(64.*under*under)
                     + 189.  * vkappa/(256.*under)
                     + (63. * sqrt(3./2.)) * atan2(vkappa,sqrt(6.))/256.));
         }

      }
};

/* Kappa = 2 version */
class Kappa2 : public Kappa
{
   public:
      Kappa2(std::default_random_engine& _rand, double _mass, double _charge, double _w0, double _maxw0)
         : Kappa(_rand,_mass,_charge,_w0,_maxw0){
            generate_lookup();
      }
      Kappa2(std::default_random_engine& _rand);

   private:
      virtual void generate_lookup() {
         double prefix = 4./M_PI * sqrt(2.);

         for(int i=0; i<lookup_size; i++) {
            double vkappa = (maxw0/lookup_size)*i;
            double under = vkappa*vkappa + 2.;

            lookup.push_back( prefix*(
                     - 2. * vkappa/(under*under)
                     +      vkappa/(2.*under)
                     + atan2(vkappa,sqrt(2.))/(2.*sqrt(2.))
                     )
                  );
         }

      }
};

template<typename T> Distribution* createDistribution(std::default_random_engine& rand) {
   return new T(rand);
}
