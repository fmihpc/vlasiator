/*
 * This file is part of Vlasiator.
 * 
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://www.physics.helsinki.fi/vlasiator/
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */


#ifndef FIELDTRACING_H
#define FIELDTRACING_H

#include <cstdlib>
#include <array>
#include "../common.h"
#include "../fieldsolver/fs_common.h"
#include "../fieldsolver/derivatives.hpp"
#include "../sysboundary/ionosphere.h"

// Used in full box + flux rope tracing, the others used in coupling should use Real as double probably.
typedef float TReal;

// Get the (integer valued) global fsgrid cell index (i,j,k) for the magnetic-field traced mapping point that node n is
// associated with
template<class T> std::array<FsGridTools::FsSize_t, 3> getGlobalFsGridCellIndexForCoord(T& grid,const std::array<Real, 3>& x) {
   std::array<FsGridTools::FsSize_t, 3> retval;
   retval[0] = floor((x[0] - grid.physicalGlobalStart[0]) / grid.DX);
   retval[1] = floor((x[1] - grid.physicalGlobalStart[1]) / grid.DY);
   retval[2] = floor((x[2] - grid.physicalGlobalStart[2]) / grid.DZ);
   return retval;
}
// Get the (integer valued) local fsgrid cell index (i,j,k) for the magnetic-field traced mapping point that node n is
// associated with If the cell is not in our local domain, will return {-1,-1,-1}
template<class T> std::array<FsGridTools::FsIndex_t, 3> getLocalFsGridCellIndexForCoord(T& grid, const std::array<Real, 3>& x) {
   std::array<FsGridTools::FsSize_t, 3> globalInd = getGlobalFsGridCellIndexForCoord(grid,x);
   std::array<FsGridTools::FsIndex_t, 3> retval = grid.globalToLocal(globalInd[0], globalInd[1], globalInd[2]);
   return retval;
}
// Get the (integer valued) local fsgrid cell index (i,j,k) for the magnetic-field traced mapping point that node n is associated with
// This includes indices beyond local size (positive and negative) as we need to access ghost cells
template<class T> std::array<FsGridTools::FsIndex_t, 3> getLocalFsGridCellIndexWithGhostsForCoord(T& grid, const std::array<Real, 3>& x) {
   std::array<FsGridTools::FsSize_t, 3> globalInd = getGlobalFsGridCellIndexForCoord(grid,x);
   const std::array<FsGridTools::FsIndex_t, 3> localStart = grid.getLocalStart();
   std::array<FsGridTools::FsIndex_t,3> retval = {(FsGridTools::FsIndex_t)globalInd[0]-localStart[0], (FsGridTools::FsIndex_t)globalInd[1]-localStart[1], (FsGridTools::FsIndex_t)globalInd[2]-localStart[2]};
   return retval;
}
// Get the fraction fsgrid cell index for the magnetic-field traced mapping point that node n is associated with.
// Note that these are floating point values between 0 and 1
template<class T> std::array<Real, 3> getFractionalFsGridCellForCoord(T& grid, const std::array<Real, 3>& x) {
   std::array<Real, 3> retval;
   std::array<FsGridTools::FsSize_t, 3> fsgridCell = getGlobalFsGridCellIndexForCoord(grid,x);
   retval[0] = (x[0] - grid.physicalGlobalStart[0]) / grid.DX - fsgridCell[0];
   retval[1] = (x[1] - grid.physicalGlobalStart[1]) / grid.DY - fsgridCell[1];
   retval[2] = (x[2] - grid.physicalGlobalStart[2]) / grid.DZ - fsgridCell[2];
   return retval;
}

namespace FieldTracing {
   
   enum Direction {
      FORWARD,
      BACKWARD
   };
   
   /*! Field line integrator for Magnetosphere<->Ionosphere coupling */
   enum TracingMethod { 
      Euler,        // Euler stepping (constant stepsize)
      ADPT_Euler,   // Adaptive Euler stepping (adaptive stepsize)
      BS,           // Bulirsch-Stoer Stepping (adaptive stepsize)
      DPrince       // Dormand-Prince Stepping (adaptive stepsize) 
   };
   
   struct FieldTracingParameters {
      bool doTraceOpenClosed=false;
      bool doTraceFullBox=false;
      bool useCache=false;
      TracingMethod tracingMethod;
      Real max_allowed_error; /*!< Maximum alowed error for the adaptive field line tracing methods */
      uint32_t max_field_tracer_attempts; /*!< Max allowed attempts for the iterative field tracers */
      Real min_tracer_dx_full_box; /*!< Min allowed tracer dx for tracing in the full domain to avoid getting bogged down in the archipelago */
      const Real min_tracer_dx_ionospere_coupling=50e3; /*!< Min allowed tracer dx for tracing between the Vlasov domain and the ionosphere */
      const Real max_tracer_dx_ionospere_coupling=100e3; /*!< Max allowed tracer dx for tracing between the Vlasov domain and the ionosphere */
      Real fullbox_max_incomplete_cells; /*!< Max allowed fraction of cells left unfinished before exiting tracing loop, fullbox */
      Real fluxrope_max_incomplete_cells; /*!< Max allowed fraction of cells left unfinished before exiting tracing loop, fluxrope */
      Real fullbox_and_fluxrope_max_distance; /*!< Max allowed tracing distance before ending tracing, fullbox and fluxrope tracing */
      std::map< std::array<int, 3>, std::array<Real, Rec::N_REC_COEFFICIENTS> > reconstructionCoefficientsCache; /*!< cache for Balsara reconstruction coefficients */
      Real fluxrope_max_curvature_radii_to_trace;
      Real fluxrope_max_curvature_radii_extent;
      Real innerBoundaryRadius=0; /*!< If non-zero this will be used to determine CLOSED field lines. */
   };
   
   extern FieldTracingParameters fieldTracingParameters;
   
   /*! Type of field line ending, used to classify the ionospheric nodes and the forward and backward field lines in full.box tracing.
   * CLOSED: ends in the ionosphere
   * OPEN: exits the simulation domain
   * DANGLING: has not exited, might keep looping or would exit given enough time/steps. Not called LOOP to avoid confusion with fluxrope tracing.
   * UNPROCESSED: cells inside the ionosphere or outside the outer limits that weren't even processed in the first place
   */
   enum TracingLineEndType {
      UNPROCESSED, // Keep first for the reductions to work!
      CLOSED,
      OPEN,
      DANGLING,
      OUTSIDE,
      N_TYPES
   };

   /*! Type of connection for a point traced forward and backward.
   * The first is for the forward end, the second for the backward end.
   * Used for full-box tracing.
   * See TracingLineEndTypes for explanation of types.
   * INVALID used for UNPROCESSED or other unparsed values ==> bug?
   * \sa TracingLineEndTypes
   */
   enum TracingPointConnectionType {
      CLOSED_CLOSED,
      CLOSED_OPEN,
      OPEN_CLOSED,
      OPEN_OPEN,
      CLOSED_DANGLING,
      DANGLING_CLOSED,
      OPEN_DANGLING,
      DANGLING_OPEN,
      DANGLING_DANGLING,
      INVALID // Keep last for the reductions to work!
   };

   /*! Simple method to translate 3D to 1D indices */
   inline int ijk2Index(
      int i,
      int j,
      int k,
      std::array<int,3>dims
   ) {
      return i + j*dims[0] +k*dims[0]*dims[1];
   }

   /*! Handler function for field line tracing */
   template<typename REAL> using TracingFieldFunction = std::function<bool(std::array<REAL,3>&, const bool, std::array<REAL, 3>&)>;
   
   template<typename REAL> bool traceFullFieldFunction(
      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
      std::array<REAL,3>& r,
      const bool alongB,
      std::array<REAL,3>& b
   ) {
      if(   r[0] > P::xmax - 2*P::dx_ini
         || r[0] < P::xmin + 2*P::dx_ini
         || r[1] > P::ymax - 2*P::dy_ini
         || r[1] < P::ymin + 2*P::dy_ini
         || r[2] > P::zmax - 2*P::dz_ini
         || r[2] < P::zmin + 2*P::dz_ini
      ) {
         cerr << (string)("(fieldtracing) Error: fsgrid coupling trying to step outside of the global domain?\n");
         return false;
      }
      
      // Get field direction
      b[0] = SBC::ionosphereGrid.dipoleField(r[0],r[1],r[2],X,0,X) + SBC::ionosphereGrid.BGB[0];
      b[1] = SBC::ionosphereGrid.dipoleField(r[0],r[1],r[2],Y,0,Y) + SBC::ionosphereGrid.BGB[1];
      b[2] = SBC::ionosphereGrid.dipoleField(r[0],r[1],r[2],Z,0,Z) + SBC::ionosphereGrid.BGB[2];
      
      std::array<FsGridTools::FsSize_t, 3> fsgridCellu = getGlobalFsGridCellIndexForCoord(technicalGrid,{(TReal)r[0], (TReal)r[1], (TReal)r[2]});
      std::array<FsGridTools::FsIndex_t,3> fsgridCell = {(FsGridTools::FsIndex_t)fsgridCellu[0],(FsGridTools::FsIndex_t)fsgridCellu[1],(FsGridTools::FsIndex_t)fsgridCellu[2]};
      const std::array<FsGridTools::FsIndex_t, 3> localStart = technicalGrid.getLocalStart();
      const std::array<FsGridTools::FsIndex_t, 3> localSize = technicalGrid.getLocalSize();
      // Make the global index a local one, bypass the fsgrid function that yields (-1,-1,-1) also for ghost cells.
      fsgridCell[0] -= localStart[0];
      fsgridCell[1] -= localStart[1];
      fsgridCell[2] -= localStart[2];
      
      if(fsgridCell[0] > localSize[0] || fsgridCell[1] > localSize[1] || fsgridCell[2] > localSize[2]
         || fsgridCell[0] < -1 || fsgridCell[1] < -1 || fsgridCell[2] < -1) {
         cerr << (string)("(fieldtracing) Error: fsgrid coupling trying to access local ID " + to_string(fsgridCell[0]) + " " + to_string(fsgridCell[1]) + " " + to_string(fsgridCell[2])
         + " for local domain size " + to_string(localSize[0]) + " " + to_string(localSize[1]) + " " + to_string(localSize[2])
         + " at position " + to_string(r[0]) + " " + to_string(r[1]) + " " + to_string(r[2]) + " radius " + to_string(sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]))
         + "\n");
         abort();
         return false;
      } else {
         if(technicalGrid.get(fsgridCell[0],fsgridCell[1],fsgridCell[2])->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
            const std::array<Real, 3> perB = interpolatePerturbedB(
               perBGrid,
               dPerBGrid,
               technicalGrid,
               fieldTracingParameters.reconstructionCoefficientsCache,
               fsgridCell[0],fsgridCell[1],fsgridCell[2],
               {(Real)r[0], (Real)r[1], (Real)r[2]}
            );
            b[0] += perB[0];
            b[1] += perB[1];
            b[2] += perB[2];
         }
      }
      
      // Normalize
      REAL  norm = 1. / sqrt(b[0]*b[0] + b[1]*b[1] + b[2]*b[2]);
      for(int c=0; c<3; c++) {
         b[c] = b[c] * norm;
      }
      
      // Make sure motion is outwards. Flip b if dot(r,b) < 0
      if(!(std::isfinite(b[0]) && std::isfinite(b[1]) && std::isfinite(b[2]))) {
         cerr << "(fieldtracing) Error: magnetic field is nan or inf in getRadialBfieldDirection at location "
         << r[0] << ", " << r[1] << ", " << r[2] << ", with B = " << b[0] << ", " << b[1] << ", " << b[2] << endl;
         b[0] = 0;
         b[1] = 0;
         b[2] = 0;
      }
      if(!alongB) { // In this function, outwards indicates whether we trace along (true) or against (false) the field direction
         b[0] *= -1;
         b[1] *= -1;
         b[2] *= -1;
      }
      return true;
   }


   /*Modified Midpoint Method used by the Bulirsch Stoer integrations
   * stepsize: initial step  size
   * r: initial position 
   * r1: new position
   * n: number of substeps
   * stepsize: big stepsize to use
   * z0,zmid,z2: intermediate approximations
   * */
   template <typename REAL> void modifiedMidpointMethod(
      std::array<REAL,3> r,
      std::array<REAL,3>& r1,
      int n,
      REAL stepSize,
      TracingFieldFunction<REAL>& BFieldFunction,
      const bool outwards=true
   ){
      //Allocate some memory.
      std::array<REAL,3> bunit,crd,z0,zmid,z1;
      //Divide by number of sub steps      
      REAL h= stepSize/(REAL)n;
      
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
   }; // Modified Midpoint Method used by BS step

   /*! Richardson extrapolation using polynomial fitting used by the Bulirsch-Stoer Method */
   template <typename REAL> void richardsonExtrapolation(
      int i,
      std::vector<REAL>& table,
      REAL& maxError,
      std::array<int,3>dims
   ) {
      int k;
      maxError = 0;
      for (int dim=0; dim<3; dim++){
         for(k =1; k<i+1; k++){
            
            table.at(ijk2Index(i,k,dim,dims)) = table.at(ijk2Index(i,k-1,dim,dims))  +(table.at(ijk2Index(i,k-1,dim,dims)) -table.at(ijk2Index(i-1,k-1,dim,dims)))/(std::pow(4,i) -1);
         }
         
         REAL thisError = fabs(table.at(ijk2Index(k-1,k-1,dim,dims))   -  table.at(ijk2Index(k-2,k-2,dim,dims)));
         if(thisError > maxError) {
            maxError = thisError;
         }
      }

   }; //Richardson extrapolation method used by BS step
   
   template <typename REAL> bool bulirschStoerStep(
      std::array<REAL, 3>& r,
      std::array<REAL, 3>& b,
      REAL& stepSize,
      const REAL minStepSize,
      const REAL maxStepSize,
      TracingFieldFunction<REAL>& BFieldFunction,
      const bool outwards=true
   ) {
      //Factors by which the stepsize is multiplied 
      REAL shrink = 0.95;
      REAL grow = 1.2; 
      //Max substeps for midpoint method
      int kMax = 8;
      //Optimal row to converge at 
      int kOpt = 6;
      
      const int ndim = kMax*kMax*3;
      std::array<int,3>  dims={kMax,kMax,3};
      std::vector<REAL>table(ndim);
      std::array<REAL,3> rold,rnew,r1;
      REAL error;
      
      //Get B field unit vector in case we don't converge yet
      BFieldFunction(r,outwards,b);
      
      //Let's start things up with 2 substeps
      int n =2;
      //Save old state
      rold = r;
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
   }; //Bulirsch Stoer step


   template<typename REAL> bool dormandPrinceStep(
      std::array<REAL, 3>& r,
      std::array<REAL, 3>& b,
      REAL& stepSize,
      const REAL minStepSize,
      const REAL maxStepSize,
      TracingFieldFunction<REAL>& BFieldFunction,
      const bool outwards=true
   ) {
      // BFieldFunction can return false if it steps out of the "comfort zone"
      bool proceed = true;
      std::array<REAL,7> kx,ky,kz;
      std::array<REAL,3> b_unit;
      std::array<REAL,3> _r{0,0,0};
      
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
      
      REAL err=0;
      std::array<REAL,3>rf;
      if (proceed) {
         //Error calculation
         std::array<REAL,3>error_xyz;
         rf[0]=r[0] +(35./384.)*kx[0] + (500./1113.)*kx[2] + (125./192.)*kx[3] - (2187./6784.)*kx[4] +(11./84.)*kx[5];
         rf[1]=r[1] +(35./384.)*ky[0] + (500./1113.)*ky[2] + (125./192.)*ky[3] - (2187./6784.)*ky[4] +(11./84.)*ky[5];
         rf[2]=r[2] +(35./384.)*kz[0] + (500./1113.)*kz[2] + (125./192.)*kz[3] - (2187./6784.)*kz[4] +(11./84.)*kz[5];
         
         error_xyz[0]=abs((71./57600.)*kx[0] -(71./16695.)*kx[2] + (71./1920.)*kx[3] -(17253./339200.)*kx[4]+(22./525.)*kx[5] -(1./40.)*kx[6] );
         error_xyz[1]=abs((71./57600.)*ky[0] -(71./16695.)*ky[2] + (71./1920.)*ky[3] -(17253./339200.)*ky[4]+(22./525.)*ky[5] -(1./40.)*ky[6] );
         error_xyz[2]=abs((71./57600.)*kz[0] -(71./16695.)*kz[2] + (71./1920.)*kz[3] -(17253./339200.)*kz[4]+(22./525.)*kz[5] -(1./40.)*kz[6] );
         
         //Estimate proper stepsize
         err=std::max( std::max(error_xyz[0], error_xyz[1]), error_xyz[2]);
         REAL s=pow((fieldTracingParameters.max_allowed_error/(2*err)),1./5.);
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

   }; //Dormand Prince step


   template<typename REAL> bool adaptiveEulerStep(
      std::array<REAL, 3>& r,
      std::array<REAL, 3>& b,
      REAL& stepSize,
      const REAL minStepSize,
      const REAL maxStepSize,
      TracingFieldFunction<REAL>& BFieldFunction,
      const bool outwards=true
   ) {
      //First evaluation
      std::array<REAL, 3> r1;
      BFieldFunction(r,outwards,b);
      
      for(int c=0; c<3; c++) {
         r1[c] = r[c]+ stepSize * b[c];
      }
      
      //Second more accurate evaluation
      std::array<REAL, 3> r2,b2;
      for(int c=0; c<3; c++) {
         r2[c] = r[c]+ 0.5*stepSize * b[c];
      }
      
      BFieldFunction(r2,outwards,b2);
      for(int c=0; c<3; c++) {
         r2[c] = r2[c]+ 0.5*stepSize * b2[c];
      }
      
      //Local error estimate
      std::array<REAL,3> error_xyz{
         fabs(r2[0]-r1[0]),
         fabs(r2[1]-r1[1]),
         fabs(r2[2]-r1[2])
      };
      
      //Max Error and step adjustment
      const REAL err=std::max( std::max(error_xyz[0], error_xyz[1]), error_xyz[2]);
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
   } ; //Adaptive Euler step


   template<typename REAL> void eulerStep(
      std::array<REAL, 3>& x,
      std::array<REAL, 3>& v,
      REAL& stepSize,
      TracingFieldFunction<REAL>& BFieldFunction,
      const bool outwards=true
   ) {
      // Get field direction
      BFieldFunction(x,outwards,v);
      
      for(int c=0; c<3; c++) {
         x[c] += stepSize * v[c];
      }
   }; //Euler step


   /*! Take a step along the field line*/
   template<typename REAL> void stepFieldLine(
      std::array<REAL, 3>& x,
      std::array<REAL, 3>& v,
      REAL& stepsize,
      const REAL minStepSize,
      const REAL maxStepSize,
      TracingMethod method,
      TracingFieldFunction<REAL>& BFieldFunction,
      const bool outwards
   ) {
      bool reTrace;
      uint32_t attempts=0;
      switch(method) {
         case Euler:
            eulerStep(x, v,stepsize, BFieldFunction, outwards);
            break;
         case ADPT_Euler:
            do {
               reTrace=!adaptiveEulerStep(x, v, stepsize, minStepSize, maxStepSize, BFieldFunction, outwards);
               attempts++;
            } while (reTrace && attempts<= fieldTracingParameters.max_field_tracer_attempts);
            if (reTrace) {
               logFile << "(fieldtracing) Warning: Adaptive Euler field line tracer exhausted all available attempts and still did not converge." << std::endl;
            }
            break;
         case BS:
            do{
               reTrace=!bulirschStoerStep(x, v, stepsize, minStepSize, maxStepSize, BFieldFunction, outwards);
               attempts++;
            } while (reTrace && attempts<= fieldTracingParameters.max_field_tracer_attempts);
            break;
         case DPrince:
            do {
               reTrace=!dormandPrinceStep(x, v, stepsize, minStepSize, maxStepSize, BFieldFunction, outwards);
               attempts++;
            } while (reTrace && attempts<= fieldTracingParameters.max_field_tracer_attempts);
            if (reTrace) {
               logFile << "(fieldtracing) Warning: Dormand Prince field line tracer exhausted all available attempts and still did not converge..." << std::endl;
            }
            break;
         default:
            std::cerr << "(fieldtracing) Error: No field line tracing method defined."<<std::endl;
            abort();
            break;
      }
   }//stepFieldLine
   
   /*! function to empty the Balsara reconstruction coefficient cache at a new time step */
   inline void resetReconstructionCoefficientsCache() {
      fieldTracingParameters.reconstructionCoefficientsCache.clear();
   }
   
   /*! Link each ionospheric node to fsgrid cells for coupling */
   void calculateIonosphereFsgridCoupling(
      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
      std::vector<SBC::SphericalTriGrid::Node> & nodes,
      creal radius
   );

   /*! Find coupled ionosphere mesh node for given location */
   std::array<std::pair<int, Real>, 3> calculateIonosphereVlasovGridCoupling(
      std::array<Real,3> x,
      std::vector<SBC::SphericalTriGrid::Node> & nodes,
      creal couplingRadius
   );

   /*! Compute whether a node is connected to the ionosphere or the IMF. */
   void traceOpenClosedConnection(
      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
      std::vector<SBC::SphericalTriGrid::Node> & nodes
   );

   /*! Trace magnetic field lines forward and backward from each DCCRG cell to record the connectivity and detect flux ropes. */
   void traceFullBoxConnectionAndFluxRopes(
      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
      dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid
   );

   void reduceData(
      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
      dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry> & mpiGrid,
      std::vector<SBC::SphericalTriGrid::Node> & nodes
   );

} // namespace FieldTracing




#endif
