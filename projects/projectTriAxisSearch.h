#ifndef TRIAXISSEARCH_H
#define TRIAXISSEARCH_H

#include "project.h"

using namespace std;

namespace projects {
   class TriAxisSearch: public Project {
      public:
         
      protected:
         /*! \brief Find blocks above the threshold centred isotropically around a bulk velocity.
          * 
          * Instead of looping through the whole velocity space this function starts from the project's bulk velocity V0[3].
          * It then proceeds along V[XYZ] successively to determine at what maximum radius a block falls below (0.1 times) the threshold.
          * 
          * This radius is used to determine all blocks within that radius of V0, create them and return their list for initialisation.
          */
         virtual vector<uint> findBlocksToInitialize(SpatialCell* cell,const int& popID);
         
         /*! \brief Return a vector containing the velocity coordinate of the centre of each ion population in the distribution.
          * 
          * This function is used by findBlocksToInitialize() to start the search of the extent of the distribution along each axis.
          * The extension to a vector allows to have more than one population in each spatial cell.
          * 
          * \sa findBlocksToInitialize
          */
         virtual vector<std::array<Real, 3>> getV0(
            creal x,
            creal y,
            creal z
         );
   }; // class TriAxisSearch
} // namespace


#endif
