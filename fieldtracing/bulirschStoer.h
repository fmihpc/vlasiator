
#ifndef BULIRSCHSTOER_H
#define BULIRSCHSTOER_H

#include "fieldtracing.h"

namespace FieldTracing {
   
   bool bulirschStoerStep(
      std::array<Real, 3>& r,
      std::array<Real, 3>& b,
      Real& stepSize,
      creal minStepSize,
      creal maxStepSize,
      TracingFieldFunction& BFieldFunction,
      const bool outwards=true
   ); //Bulirsch Stoer step

   void modifiedMidpointMethod(
      std::array<Real,3> r,
      std::array<Real,3>& r1,
      int n,
      Real stepSize,
      TracingFieldFunction& BFieldFunction,
      const bool outwards=true
   ); // Modified Midpoint Method used by BS step

   void richardsonExtrapolation(
      int i,
      std::vector<Real>& table,
      Real& maxError,
      std::array<int,3>dims
   ); //Richardson extrapolation method used by BS step
   
} // namespace FieldTracing

#endif
