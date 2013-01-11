#ifndef ISOTROPICMAXWELLIAN_H
#define ISOTROPICMAXWELLIAN_H

#include "project.h"

using namespace std;

namespace projects {
   class IsotropicMaxwellian: public Project {
      public:
         
      protected:
         /*! \brief Find blocks above the threshold centred isotropically around a bulk velocity.
          * 
          * Instead of looping through the whole velocity space this function starts from the project's bulk velocity V0[3].
          * It then proceeds along VX at fixed V0[1], V0[2] to determine at what radius a block falls below (0.1 times) the threshold.
          * 
          * This radius is used to determine all blocks within that radius of V0, create them and return their list for initialisation.
          */
         virtual vector<uint> findBlocksToInitialize(SpatialCell* cell);
         virtual Real getV0(
            creal x,
            creal y,
            creal z,
            cuint component
         );
   }; // class IsotropicMaxwellian
} // namespace


#endif
