

#ifndef EULER_H
#define EULER_H

#include "fieldtracing.h"

namespace FieldTracing {
   
   void eulerStep(
      std::array<Real, 3>& x,
      std::array<Real, 3>& v,
      Real& stepSize,
      TracingFieldFunction& BFieldFunction,
      const bool outwards=true
   ); //Euler step
   
} // namespace FieldTracing

#endif
