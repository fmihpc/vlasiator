
#include "dormandPrince.h"

namespace FieldTracing {
//    extern FieldTracingParameters fieldTracingParameters;
   
   bool dormandPrinceStep( 
   std::array<Real, 3>& r,
   std::array<Real, 3>& b,
   Real& stepSize,
   creal minStepSize,
   creal maxStepSize,
   TracingFieldFunction& BFieldFunction,
   const bool outwards
   ){
      // BFieldFunction can return false if it steps out of the "comfort zone"
      bool proceed = true;
      std::array<Real,7> kx,ky,kz;
      std::array<Real,3> b_unit;
      std::array<Real,3> _r{0,0,0};
      
      //K1 slope
      proceed = BFieldFunction(r,outwards,b_unit);
      kx[0]=stepSize*b_unit[0];
      ky[0]=stepSize*b_unit[1];
      kz[0]=stepSize*b_unit[2];
      
      //K2 slope
      _r[0]=r[0]+(1./5.)*kx[0];
      _r[1]=r[1]+(1./5.)*ky[0];
      _r[2]=r[2]+(1./5.)*kz[0];
      if (proceed) {
         proceed = BFieldFunction(_r,outwards,b_unit);
      }
      kx[1]=stepSize*b_unit[0];
      ky[1]=stepSize*b_unit[1];
      kz[1]=stepSize*b_unit[2];
      
      //K3 slope  
      _r[0]=r[0]+ (3./10.)*kx[1]; 
      _r[1]=r[1]+ (3./10.)*ky[1]; 
      _r[2]=r[2]+ (3./10.)*kz[1]; 
      if (proceed) {
         proceed = BFieldFunction(_r,outwards,b_unit);
      }
      kx[2]=stepSize*b_unit[0];
      ky[2]=stepSize*b_unit[1];
      kz[2]=stepSize*b_unit[2];
      
      //K4 slope
      _r[0]=r[0]+(4./5.)*kx[2];  
      _r[1]=r[1]+(4./5.)*ky[2];  
      _r[2]=r[2]+(4./5.)*kz[2];  
      if (proceed) {
         proceed = BFieldFunction(_r,outwards,b_unit);
      }
      kx[3]=stepSize*b_unit[0];
      ky[3]=stepSize*b_unit[1];
      kz[3]=stepSize*b_unit[2];
      
      //K5 slope  
      _r[0]=r[0]+(8./9.)*kx[3];   
      _r[1]=r[1]+(8./9.)*ky[3];   
      _r[2]=r[2]+(8./9.)*kz[3];   
      if (proceed) {
         proceed = BFieldFunction(_r,outwards,b_unit);
      }
      kx[4]=stepSize*b_unit[0];
      ky[4]=stepSize*b_unit[1];
      kz[4]=stepSize*b_unit[2];
      
      //K6 slope
      _r[0]=r[0]+kx[4];  
      _r[1]=r[1]+ky[4];  
      _r[2]=r[2]+kz[4];  
      if (proceed) {
         proceed = BFieldFunction(_r,outwards,b_unit);
      }
      kx[5]=stepSize*b_unit[0];
      ky[5]=stepSize*b_unit[1];
      kz[5]=stepSize*b_unit[2];
      
      //K7 slope 
      _r[0]=r[0]+ kx[5];  
      _r[1]=r[1]+ ky[5];  
      _r[2]=r[2]+ kz[5];  
      if (proceed) {
         proceed = BFieldFunction(_r,outwards,b_unit);
      }
      kx[6]=stepSize*b_unit[0];
      ky[6]=stepSize*b_unit[1];
      kz[6]=stepSize*b_unit[2];
      
      Real err=0;
      std::array<Real,3>rf;
      if (proceed) {
         //Error calculation
         std::array<Real,3>error_xyz;
         rf[0]=r[0] +(35./384.)*kx[0] + (500./1113.)*kx[2] + (125./192.)*kx[3] - (2187./6784.)*kx[4] +(11./84.)*kx[5];
         rf[1]=r[1] +(35./384.)*ky[0] + (500./1113.)*ky[2] + (125./192.)*ky[3] - (2187./6784.)*ky[4] +(11./84.)*ky[5];
         rf[2]=r[2] +(35./384.)*kz[0] + (500./1113.)*kz[2] + (125./192.)*kz[3] - (2187./6784.)*kz[4] +(11./84.)*kz[5];
         
         error_xyz[0]=abs((71./57600.)*kx[0] -(71./16695.)*kx[2] + (71./1920.)*kx[3] -(17253./339200.)*kx[4]+(22./525.)*kx[5] -(1./40.)*kx[6] );
         error_xyz[1]=abs((71./57600.)*ky[0] -(71./16695.)*ky[2] + (71./1920.)*ky[3] -(17253./339200.)*ky[4]+(22./525.)*ky[5] -(1./40.)*ky[6] );
         error_xyz[2]=abs((71./57600.)*kz[0] -(71./16695.)*kz[2] + (71./1920.)*kz[3] -(17253./339200.)*kz[4]+(22./525.)*kz[5] -(1./40.)*kz[6] );
         
         //Estimate proper stepsize
         err=std::max( std::max(error_xyz[0], error_xyz[1]), error_xyz[2]);
         Real s=pow((fieldTracingParameters.max_allowed_error/(2*err)),1./5.);
         stepSize=stepSize*s;
      } else { // proceed is false, we probably stepped too far
         stepSize /= 2;
      }
      
      stepSize = stepSize > maxStepSize ? maxStepSize : stepSize;
      stepSize = stepSize < minStepSize ? minStepSize : stepSize;
      if ((err>fieldTracingParameters.max_allowed_error && stepSize > minStepSize) || !proceed) {
         return false;
      }
      // No explicit else so the compiler sees a return at the end of the non-void function.
      // With the return in if it's the same.
      // Evaluate B at the final stepping point to get the b vector of the fieldline
      BFieldFunction(r,outwards,b);
      r=rf;
      return true;
   }//Dormand Prince Step
   
} // namespace FieldTracing
