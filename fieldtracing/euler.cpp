#include "euler.h"

namespace FieldTracing {
   
   /*! Take an Euler step */
   void eulerStep(
      std::array<Real, 3>& x,
      std::array<Real, 3>& v,
      Real& stepsize,
      TracingFieldFunction& BFieldFunction,
      const bool outwards
   ) {
      // Get field direction
      BFieldFunction(x,outwards,v);
      
      for(int c=0; c<3; c++) {
         x[c] += stepsize * v[c];
      }
      
   } // Euler Step

} // namespace FieldTracing
