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

/* Kappa distributions as taken from the CSA code */
class Kappa : public Distribution {

   protected:
      std::uniform_real_distribution<Real> r;
      std::vector<Real> lookup;
      const static int lookup_size = 65536;
      Real w0;
      Real maxw0;

      virtual void generate_lookup() = 0;

      Real find_v_for_r(Real rand);

   public:
      Kappa(std::default_random_engine& _rand, Real _mass, Real _charge, Real _w0, Real _maxw0) : Distribution(_rand,_mass,_charge),w0(_w0),maxw0(_maxw0) {
      }
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
};

/* Kappa = 6 version */
class Kappa6 : public Kappa {

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
   public:
      Kappa6(std::default_random_engine& _rand, double _mass, double _charge, double _w0, double _maxw0) : Kappa(_rand,_mass,_charge,_w0,_maxw0){
         generate_lookup();
      }
      Kappa6(std::default_random_engine& _rand);
};

/* Kappa = 2 version */
class Kappa2 : public Kappa {

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
   public:
      Kappa2(std::default_random_engine& _rand, double _mass, double _charge, double _w0, double _maxw0) : Kappa(_rand,_mass,_charge,_w0,_maxw0){
         generate_lookup();
      }
      Kappa2(std::default_random_engine& _rand);
};
template<typename T> Distribution* create_distribution(std::default_random_engine& rand) {
   return new T(rand);
}
