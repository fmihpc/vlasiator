
#include "bulirschStoer.h"

namespace FieldTracing {
//    extern FieldTracingParameters fieldTracingParameters;
   
   /*! Richardson extrapolation using polynomial fitting used by the Bulirsch-Stoer Method */
   void richardsonExtrapolation(int i, std::vector<Real>&table , Real& maxError, std::array<int,3>dims ){
      int k;
      maxError = 0;
      for (int dim=0; dim<3; dim++){
         for(k =1; k<i+1; k++){
            
            table.at(ijk2Index(i,k,dim,dims)) = table.at(ijk2Index(i,k-1,dim,dims))  +(table.at(ijk2Index(i,k-1,dim,dims)) -table.at(ijk2Index(i-1,k-1,dim,dims)))/(std::pow(4,i) -1);

         }
         
         Real thisError = fabs(table.at(ijk2Index(k-1,k-1,dim,dims))   -  table.at(ijk2Index(k-2,k-2,dim,dims)));
         if(thisError > maxError) {
            maxError = thisError;
         }
      }
   }

   /*Modified Midpoint Method used by the Bulirsch Stoer integrations
   * stepsize: initial step  size
   * r: initial position 
   * r1: new position
   * n: number of substeps
   * stepsize: big stepsize to use
   * z0,zmid,z2: intermediate approximations
   * */
   void modifiedMidpointMethod(
      std::array<Real,3> r,
      std::array<Real,3>& r1,
      int n,
      Real stepsize,
      TracingFieldFunction& BFieldFunction,
      const bool outwards
   ) {
      //Allocate some memory.
      std::array<Real,3> bunit,crd,z0,zmid,z1;
      //Divide by number of sub steps      
      Real h= stepsize/(Real)n;
      //Real norm;
      
      //First step 
      BFieldFunction(r,outwards,bunit);
      z0=r;
      z1={ r[0]+h*bunit[0], r[1]+h*bunit[1], r[2]+h*bunit[2]  };
      BFieldFunction(z1,outwards,bunit);
      
      crd = { r[0]+h*bunit[0], r[1]+h*bunit[1], r[2]+h*bunit[2] };
      
      for (int m =0; m<=n; m++){
         zmid= { z0[0]+2*h*bunit[0], z0[1]+2*h*bunit[1], z0[2]+2*h*bunit[2] };
         z0=z1;
         z1=zmid;
         crd = { crd[0]+h*bunit[0], crd[1]+h*bunit[1], crd[2]+h*bunit[2] };
         BFieldFunction(crd,outwards,bunit);
      }
      
      //These are now are new position
      for (int c=0; c<3; c++){
         r1[c] = 0.5*(z0[c]+z1[c]+h*bunit[c]);
      }
      
   } //modifiedMidpoint Method

   /*! Bulirsch-Stoer Method to trace field line to next point along it */
   bool bulirschStoerStep(
      std::array<Real, 3>& r,
      std::array<Real, 3>& b,
      Real& stepSize,
      creal minStepSize,
      creal maxStepSize,
      TracingFieldFunction& BFieldFunction,
      const bool outwards
   ) {
      //Factors by which the stepsize is multiplied 
      Real shrink = 0.95;
      Real grow = 1.2; 
      //Max substeps for midpoint method
      int kMax = 8;
      //Optimal row to converge at 
      int kOpt = 6;
      
      const int ndim = kMax*kMax*3;
      std::array<int,3>  dims={kMax,kMax,3};
      std::vector<Real>table(ndim);
      std::array<Real,3> rnew,r1;
      Real error;
      
      //Get B field unit vector in case we don't converge yet
      BFieldFunction(r,outwards,b);
      
      //Let's start things up with 2 substeps
      int n =2;
      //Take a first Step
      modifiedMidpointMethod(r,r1,n,stepSize,BFieldFunction,outwards);
      
      //Save values in table
      for (int c =0; c<3; ++c) {
         table[ijk2Index(0,0,c,dims)] = r1[c];
      }
      
      for (int i=1; i<kMax; ++i) {
         
         //Increment n by 2 at every iteration.
         n+=2;
         modifiedMidpointMethod(r,rnew,n,stepSize,BFieldFunction,outwards);
         
         //Save values in table
         for (int c =0; c<3; ++c) {
            table[ijk2Index(i,0,c,dims)] = rnew[c];
         }
         
         //Now let's perform a Richardson extrapolatation
         richardsonExtrapolation(i,table,error,dims);
         
         //Normalize error
         error/=fieldTracingParameters.max_allowed_error;
         
         //If we are below eps good, let's return but also let's modify the stepSize accordingly 
         if (error<1. || stepSize == minStepSize) {
            if (i> kOpt) {
               stepSize*=shrink;
            }
            if (i< kOpt) {
               stepSize*=grow;
            }
            //Make sure stepsize does not exceed maxStepsize
            stepSize = stepSize > maxStepSize ? maxStepSize : stepSize;
            stepSize = stepSize < minStepSize ? minStepSize : stepSize;
            //Keep our new position
            r=rnew;
            //Evaluate B here
            BFieldFunction(r,outwards,b);
            return true;
         }
      }
      
      //If we end up  here it means our tracer did not converge so we need to reduce the stepSize all along and try again
      stepSize*=shrink;
      stepSize = stepSize < minStepSize ? minStepSize : stepSize;
      return false;
      
   } //Bulirsch-Stoer Step
   
} // namespace FieldTracing
