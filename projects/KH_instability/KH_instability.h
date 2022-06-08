
#ifndef KH_INSTABILITY_H
#define KH_INSTABILITY_H

#include <array>
#include <vector>

#include "../../definitions.h"
#include "../project.h"

namespace projects {

   enum {
      LEFT,
      RIGHT
   };

   struct KH_instabilitySpeciesParameters {
      Real rho[2];
      Real P;
      uint nSpaceSamples;
      uint nVelocitySamples;
   };

   class KH_instability: public Project {
    public:
      KH_instability();
      virtual ~KH_instability();

      virtual bool initialize(void);
      static void addParameters(void);
      virtual void getParameters(void);
      virtual void setProjectBField(
         FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2>& perBGrid,
         FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, 2>& BgBGrid,
         FsGrid< fsgrids::technical, 2>& technicalGrid 
      );

    protected:

      Real profile(creal left, creal right, creal x, creal y, creal z) const;
      Real getDistribValue(
                           creal& x, creal& y, creal& z,
                           creal& vx, creal& vy, creal& vz,
                           creal& dvx, creal& dvy, creal& dvz,
                           const uint popID
                          ) const;
      virtual void calcCellParameters(spatial_cell::SpatialCell* cell,creal& t);
      virtual Real calcPhaseSpaceDensity(
                                         creal& x, creal& y, creal& z,
                                         creal& dx, creal& dy, creal& dz,
                                         creal& vx, creal& vy, creal& vz,
                                         creal& rand_x, creal& rand_y, creal& rand_z,const uint popID
                                        ) const;

      Real Vx[2];
      Real Vy[2];
      Real Vz[2];
      Real Bx[2];
      Real By[2];
      Real Bz[2];
      Real transitionWidth;
      Real perturbationAmplitude;
      Real perturbationWaveNumber;
      Real noiseAmplitude;
      std::vector<Real> random_x;
      std::vector<Real> random_y;
      std::vector<Real> random_z;
      std::vector<KH_instabilitySpeciesParameters> speciesParams; 
   }; // class KH_instability 
} // namespace projects

#endif
