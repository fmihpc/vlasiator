

#ifndef EULERADAPTIVE_H
#define EULERADAPTIVE_H

#include "fieldtracing.h"

namespace FieldTracing {
   
   bool adaptiveEulerStep(
      std::array<Real, 3>& r,
      std::array<Real, 3>& b,
      Real& stepSize,
      creal minStepSize,
      creal maxStepSize,
      TracingFieldFunction& BFieldFunction,
      const bool outwards=true
   ); //Adaptive Euler step
   
} // namespace FieldTracing

#endif
