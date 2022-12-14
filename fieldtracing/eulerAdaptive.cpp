#include "eulerAdaptive.h"

namespace FieldTracing {
//    extern FieldTracingParameters fieldTracingParameters;
   
   bool adaptiveEulerStep( 
      std::array<Real, 3>& r,
      std::array<Real, 3>& b,
      Real& stepSize,
      creal minStepSize,
      creal maxStepSize,
      TracingFieldFunction& BFieldFunction,
      const bool outwards
   ) {
      //First evaluation
      std::array<Real, 3> r1;
      BFieldFunction(r,outwards,b);
      
      for(int c=0; c<3; c++) {
         r1[c] = r[c]+ stepSize * b[c];
      }
      
      //Second more accurate evaluation
      std::array<Real, 3> r2,b2;
      for(int c=0; c<3; c++) {
         r2[c] = r[c]+ 0.5*stepSize * b[c];
      }
      
      BFieldFunction(r2,outwards,b2);
      for(int c=0; c<3; c++) {
         r2[c] = r2[c]+ 0.5*stepSize * b2[c];
      }
      
      //Local error estimate
      std::array<Real,3> error_xyz{
         fabs(r2[0]-r1[0]),
         fabs(r2[1]-r1[1]),
         fabs(r2[2]-r1[2])
      };
      
      //Max Error and step adjustment
      creal err=std::max( std::max(error_xyz[0], error_xyz[1]), error_xyz[2]);
      stepSize=stepSize*sqrt(fieldTracingParameters.max_allowed_error/err);
      stepSize = stepSize > maxStepSize ? maxStepSize : stepSize;
      stepSize = stepSize < minStepSize ? minStepSize : stepSize;
      
      if (err<=fieldTracingParameters.max_allowed_error || stepSize == minStepSize) {
         // Note: B at r has been evaluated above so no need to do it here and we're not returning it anyway as it's not needed.
         r=r2;
         return true;
      } else {
         stepSize *= 0.9; // This is to avoid asymptotic convergence when the ratio above is very close to 1.
         return false;
      }
   } //Adaptive Euler Step

} // namespace FieldTracing
