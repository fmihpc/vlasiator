// SPDX-License-Identifier: GPL-2.0-or-later
#ifndef FIELDTRACING_STEP_SIZE_H
#define FIELDTRACING_STEP_SIZE_H

#include <algorithm>
#include <cmath>
#include <sstream>
#include <string>

[[noreturn]] void abort_mpi(const std::string str, const int err_type);

namespace FieldTracing {
   template<typename REAL> REAL checkedStepSize(
      const REAL stepSize,
      const REAL minStepSize,
      const REAL maxStepSize,
      const bool adaptive = true
   ) {
      if (!std::isfinite(stepSize) || !std::isfinite(minStepSize) || !std::isfinite(maxStepSize)
          || stepSize <= 0 || minStepSize <= 0 || maxStepSize < minStepSize) {
         std::ostringstream message;
         message << "(fieldtracing) Error: step lengths must be finite and positive, with minimum <= maximum."
                 << " Initial=" << stepSize << ", minimum=" << minStepSize << ", maximum=" << maxStepSize
                 << ". Check tracer limits and grid spacing.";
         abort_mpi(message.str(), 0);
      }
      // The configured minimum applies to adaptive tracers, not fixed-step Euler.
      return adaptive ? std::clamp(stepSize, minStepSize, maxStepSize) : std::min(stepSize, maxStepSize);
   }
}

#endif
