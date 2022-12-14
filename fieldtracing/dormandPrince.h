


#ifndef DORMANDPRINCE_H
#define DORMANDPRINCE_H

#include "fieldtracing.h"

namespace FieldTracing {
   
   bool dormandPrinceStep(
      std::array<Real, 3>& r,
      std::array<Real, 3>& b,
      Real& stepSize,
      creal minStepSize,
      creal maxStepSize,
      TracingFieldFunction& BFieldFunction,
      const bool outwards=true
   ); //Dormand Prince step

} // namespace FieldTracing

#endif
