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

/*!\file fieldtracing.cpp
 * \brief Implementation of the field tracing algorithms used in Magnetosphere runs and in magnetosphere-ionosphere coupling.
 */

#include "../fieldsolver/fs_common.h"
#include "fieldtracing.h"
#include "bulirschStoer.h"
#include "dormandPrince.h"
#include "euler.h"
#include "eulerAdaptive.h"

#include <Eigen/Dense>

#define Vec3d Eigen::Vector3d
#define cross_product(av,bv) (av).cross(bv)
#define dot_product(av,bv) (av).dot(bv)
#define vector_length(v) (v).norm()
#define normalize_vector(v) (v).normalized()

#include "../logger.h"
extern Logger logFile;

namespace FieldTracing {
   FieldTracingParameters fieldTracingParameters;

   /* Call the heavier operations for DROs to be called only if needed, before an IO.
    */
   void reduceData(
      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
      dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry> & mpiGrid,
      std::vector<SBC::SphericalTriGrid::Node> & nodes
   ) {
      if(fieldTracingParameters.doTraceOpenClosed) {
         traceOpenClosedConnection(technicalGrid, perBGrid, dPerBGrid, nodes);
      }
      if(fieldTracingParameters.doTraceFullBox) {
         traceFullBoxConnection(technicalGrid, perBGrid, dPerBGrid, mpiGrid);
      }
      if(fieldTracingParameters.doTraceFluxRopes) {
         traceFluxRopes(technicalGrid, perBGrid, dPerBGrid, mpiGrid);
      }
   }
   
   /*! Take a step along the field line*/
   void stepFieldLine(
      std::array<Real, 3>& x,
      std::array<Real, 3>& v,
      Real& stepsize,
      creal minStepSize,
      creal maxStepSize,
      TracingMethod method,
      TracingFieldFunction& BFieldFunction,
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
   
   bool traceFullFieldFunction(
      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
      std::array<Real,3>& r,
      const bool alongB,
      std::array<Real,3>& b
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
      
      std::array<int32_t, 3> fsgridCell = getGlobalFsGridCellIndexForCoord(technicalGrid,r);
      const std::array<int32_t, 3> localStart = technicalGrid.getLocalStart();
      const std::array<int32_t, 3> localSize = technicalGrid.getLocalSize();
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
         return false;
      } else {
         if(technicalGrid.get(fsgridCell[0],fsgridCell[1],fsgridCell[2])->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
            const std::array<Real, 3> perB = interpolatePerturbedB(
               perBGrid,
               dPerBGrid,
               technicalGrid,
               fieldTracingParameters.reconstructionCoefficientsCache,
               fsgridCell[0],fsgridCell[1],fsgridCell[2],
               r
            );
            b[0] += perB[0];
            b[1] += perB[1];
            b[2] += perB[2];
         }
      }
      
      // Normalize
      Real  norm = 1. / sqrt(b[0]*b[0] + b[1]*b[1] + b[2]*b[2]);
      for(int c=0; c<3; c++) {
         b[c] = b[c] * norm;
      }
      
      // Make sure motion is outwards. Flip b if dot(r,b) < 0
      if(std::isnan(b[0]) || std::isnan(b[1]) || std::isnan(b[2])) {
         cerr << "(fieldtracing) Error: magnetic field is nan in getRadialBfieldDirection at location "
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
   
   /*! Calculate mapping between ionospheric nodes and fsGrid cells.
   * To do so, the magnetic field lines are traced from all mesh nodes
   * outwards until a non-boundary cell is encountered. Their proportional
   * coupling values are recorded in the grid nodes.
   */
   void calculateIonosphereFsgridCoupling(
      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
      std::vector<SBC::SphericalTriGrid::Node> & nodes,
      creal couplingRadius
   ) {
      
      // we don't need to do anything if we have no nodes
      if(nodes.size() == 0) {
         return;
      }
      
      phiprof::start("fieldtracing-ionosphere-fsgridCoupling");
      // Pick an initial stepsize
      creal stepSize = min(100e3, technicalGrid.DX / 2.);
      std::vector<Real> nodeTracingStepSize(nodes.size(), stepSize); // In-flight storage of step size, needed when crossing into next MPI domain
      std::vector<Real> reducedNodeTracingStepSize(nodes.size());
      
      std::vector<Real> nodeDistance(nodes.size(), std::numeric_limits<Real>::max()); // For reduction of node coordinate in case of multiple hits
      std::vector<int> nodeNeedsContinuedTracing(nodes.size(), 1);                    // Flag, whether tracing needs to continue on another task
      std::vector<std::array<Real, 3>> nodeTracingCoordinates(nodes.size());          // In-flight node upmapping coordinates (for global reduction)
      for(uint n=0; n<nodes.size(); n++) {
         nodeTracingCoordinates.at(n) = nodes.at(n).x;
         nodes.at(n).haveCouplingData = 0;
         for (uint c=0; c<3; c++) {
            nodes.at(n).xMapped.at(c) = 0;
            nodes.at(n).parameters.at(ionosphereParameters::UPMAPPED_BX+c) = 0;
         }
      }
      bool anyNodeNeedsTracing;
      
      TracingFieldFunction tracingFullField = [&perBGrid, &dPerBGrid, &technicalGrid](std::array<Real,3>& r, const bool alongB, std::array<Real,3>& b)->bool{
         return traceFullFieldFunction(perBGrid, dPerBGrid, technicalGrid, r, alongB, b);
      };
      
      int itCount = 0;
      do {
         itCount++;
         anyNodeNeedsTracing = false;
         
         #pragma omp parallel
         {
            // Trace node coordinates outwards until a non-sysboundary cell is encountered or the local fsgrid domain has been left.
            #pragma omp for schedule(dynamic)
            for(uint n=0; n<nodes.size(); n++) {
               
               if(!nodeNeedsContinuedTracing[n]) {
                  // This node has already found its target, no need for us to do anything about it.
                  continue;
               }
               SBC::SphericalTriGrid::Node& no = nodes[n];
               
               std::array<Real, 3> x = nodeTracingCoordinates[n];
               std::array<Real, 3> v({0,0,0});
               
               while( true ) {
                  
                  // Check if the current coordinates (pre-step) are in our own domain.
                  std::array<int, 3> fsgridCell = getLocalFsGridCellIndexForCoord(technicalGrid,x);
                  // If it is not in our domain, somebody else takes care of it.
                  if(fsgridCell[0] == -1) {
                     nodeNeedsContinuedTracing[n] = 0;
                     nodeTracingCoordinates[n] = {0,0,0};
                     nodeTracingStepSize[n]=0;
                     break;
                  }
                  
                  
                  // Make one step along the fieldline
                  stepFieldLine(x,v, nodeTracingStepSize[n],fieldTracingParameters.min_tracer_dx,technicalGrid.DX/2,fieldTracingParameters.tracingMethod,tracingFullField,(no.x[2] < 0));
                  
                  // Look up the fsgrid cell belonging to these coordinates
                  fsgridCell = getLocalFsGridCellIndexForCoord(technicalGrid,x);
                  std::array<Real, 3> interpolationFactor=getFractionalFsGridCellForCoord(technicalGrid,x);
                  
                  creal distance = sqrt((x[0]-no.x[0])*(x[0]-no.x[0])+(x[1]-no.x[1])*(x[1]-no.x[1])+(x[2]-no.x[2])*(x[2]-no.x[2]));
                  
                  // TODO I simplified by just looking when we change hemispheres now.
                  // This WILL fail as soon as there is a dipole tilt.
                  // But do we need it beyond debugging? Tracing back for closed/non-mapping lines is perfectly legit (once the tracer is debugged).
                  if(sign(x[2]) != sign(no.x[2])) {
                     nodeNeedsContinuedTracing.at(n) = 0;
                     nodeTracingCoordinates.at(n) = {0,0,0};
                     break;
                  }
                  
                  // If we somehow still map into the ionosphere, we missed the 88 degree criterion but shouldn't couple there.
                  if(sqrt(x.at(0)*x.at(0) + x.at(1)*x.at(1) + x.at(2)*x.at(2)) < SBC::Ionosphere::innerRadius) {
                     // TODO drop this warning if it never occurs? To be followed.
                     cerr << (string)("(fieldtracing) Warning: Triggered mapping back into Earth from node " + to_string(n) + " at z " + to_string(no.x[2]) + "\n");
                     nodeNeedsContinuedTracing.at(n) = 0;
                     nodeTracingCoordinates.at(n) = {0,0,0};
                     break;
                  }
                  
                  // Now, after stepping, if it is no longer in our domain, another MPI rank will pick up later.
                  if(fsgridCell[0] == -1) {
                     nodeNeedsContinuedTracing[n] = 1;
                     nodeTracingCoordinates[n] = x;
                     break;
                  }
                  
                  if(
                     technicalGrid.get( fsgridCell[0], fsgridCell[1], fsgridCell[2])->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY
                     && x[0]*x[0]+x[1]*x[1]+x[2]*x[2] > SBC::Ionosphere::downmapRadius*SBC::Ionosphere::downmapRadius
                  ) {
                     
                     // Store the cells mapped coordinates and upmapped magnetic field
                     no.xMapped = x;
                     no.haveCouplingData = 1;
                     nodeDistance[n] = distance;
                     const std::array<Real, 3> perB = interpolatePerturbedB(
                        perBGrid,
                        dPerBGrid,
                        technicalGrid,
                        fieldTracingParameters.reconstructionCoefficientsCache,
                        fsgridCell[0], fsgridCell[1], fsgridCell[2],
                        x
                     );
                     no.parameters[ionosphereParameters::UPMAPPED_BX] = SBC::ionosphereGrid.dipoleField(x[0],x[1],x[2],X,0,X) + perB[0];
                     no.parameters[ionosphereParameters::UPMAPPED_BY] = SBC::ionosphereGrid.dipoleField(x[0],x[1],x[2],Y,0,Y) + perB[1];
                     no.parameters[ionosphereParameters::UPMAPPED_BZ] = SBC::ionosphereGrid.dipoleField(x[0],x[1],x[2],Z,0,Z) + perB[2];
                     
                     nodeNeedsContinuedTracing[n] = 0;
                     nodeTracingCoordinates[n] = {0,0,0};
                     break;
                  }
               } // while(true)
            } // for
         } // pragma omp parallel
         
         // Globally reduce whether any node still needs to be picked up and traced onwards
         std::vector<int> sumNodeNeedsContinuedTracing(nodes.size(), 0);
         std::vector<std::array<Real, 3>> sumNodeTracingCoordinates(nodes.size());
         MPI_Allreduce(nodeNeedsContinuedTracing.data(), sumNodeNeedsContinuedTracing.data(), nodes.size(), MPI_INT, MPI_SUM, MPI_COMM_WORLD);
         if(sizeof(Real) == sizeof(double)) {
            MPI_Allreduce(nodeTracingCoordinates.data(), sumNodeTracingCoordinates.data(), 3*nodes.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(nodeTracingStepSize.data(), reducedNodeTracingStepSize.data(), nodes.size(), MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
         } else {
            MPI_Allreduce(nodeTracingCoordinates.data(), sumNodeTracingCoordinates.data(), 3*nodes.size(), MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(nodeTracingStepSize.data(), reducedNodeTracingStepSize.data(), nodes.size(), MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
         }
         for(uint n=0; n<nodes.size(); n++) {
            if(sumNodeNeedsContinuedTracing[n] > 0) {
               anyNodeNeedsTracing=true;
               nodeNeedsContinuedTracing[n] = 1;
               
               // Update that nodes' tracing coordinates
               nodeTracingCoordinates[n][0] = sumNodeTracingCoordinates[n][0] / sumNodeNeedsContinuedTracing[n];
               nodeTracingCoordinates[n][1] = sumNodeTracingCoordinates[n][1] / sumNodeNeedsContinuedTracing[n];
               nodeTracingCoordinates[n][2] = sumNodeTracingCoordinates[n][2] / sumNodeNeedsContinuedTracing[n];
            }
            nodeTracingStepSize[n] = reducedNodeTracingStepSize[n];
         }
         
      } while(anyNodeNeedsTracing);
      
      logFile << "(fieldtracing) fsgrid coupling traced in " << itCount << " iterations of the tracing loop." << endl;
      
      std::vector<Real> reducedNodeDistance(nodes.size());
      if(sizeof(Real) == sizeof(double)) {
         MPI_Allreduce(nodeDistance.data(), reducedNodeDistance.data(), nodes.size(), MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      } else {
         MPI_Allreduce(nodeDistance.data(), reducedNodeDistance.data(), nodes.size(), MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
      }
      
      // Reduce upmapped magnetic field to be consistent on all nodes
      std::vector<Real> sendUpmappedB(3 * nodes.size());
      std::vector<Real> reducedUpmappedB(3 * nodes.size());
      // Likewise, reduce upmapped coordinates
      std::vector<Real> sendxMapped(3 * nodes.size());
      std::vector<Real> reducedxMapped(3 * nodes.size());
      // And coupling rank number
      std::vector<int> sendCouplingNum(nodes.size());
      std::vector<int> reducedCouplingNum(nodes.size());
      
      for(uint n=0; n<nodes.size(); n++) {
         SBC::SphericalTriGrid::Node& no = nodes[n];
         // Discard false hits from cells that are further out from the node
         if(nodeDistance[n] > reducedNodeDistance[n]) {
            no.haveCouplingData = 0;
            for(int c=0; c<3; c++) {
               no.parameters[ionosphereParameters::UPMAPPED_BX+c] = 0;
               no.xMapped[c] = 0;
            }
         } else {
            // Cell found, add association.
            SBC::ionosphereGrid.isCouplingInwards = true;
         }
         
         
         sendUpmappedB[3*n] = no.parameters[ionosphereParameters::UPMAPPED_BX];
         sendUpmappedB[3*n+1] = no.parameters[ionosphereParameters::UPMAPPED_BY];
         sendUpmappedB[3*n+2] = no.parameters[ionosphereParameters::UPMAPPED_BZ];
         sendxMapped[3*n] = no.xMapped[0];
         sendxMapped[3*n+1] = no.xMapped[1];
         sendxMapped[3*n+2] = no.xMapped[2];
         sendCouplingNum[n] = no.haveCouplingData;
      }
      if(sizeof(Real) == sizeof(double)) { 
         MPI_Allreduce(sendUpmappedB.data(), reducedUpmappedB.data(), 3*nodes.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(sendxMapped.data(), reducedxMapped.data(), 3*nodes.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      } else {
         MPI_Allreduce(sendUpmappedB.data(), reducedUpmappedB.data(), 3*nodes.size(), MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(sendxMapped.data(), reducedxMapped.data(), 3*nodes.size(), MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
      }
      MPI_Allreduce(sendCouplingNum.data(), reducedCouplingNum.data(), nodes.size(), MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      
      for(uint n=0; n<nodes.size(); n++) {
         SBC::SphericalTriGrid::Node& no = nodes[n];
         
         // We don't even care about nodes that couple nowhere.
         if(reducedCouplingNum[n] == 0) {
            continue;
         }
         no.parameters[ionosphereParameters::UPMAPPED_BX] = reducedUpmappedB[3*n] / reducedCouplingNum[n];
         no.parameters[ionosphereParameters::UPMAPPED_BY] = reducedUpmappedB[3*n+1] / reducedCouplingNum[n];
         no.parameters[ionosphereParameters::UPMAPPED_BZ] = reducedUpmappedB[3*n+2] / reducedCouplingNum[n];
         no.xMapped[0] = reducedxMapped[3*n] / reducedCouplingNum[n];
         no.xMapped[1] = reducedxMapped[3*n+1] / reducedCouplingNum[n];
         no.xMapped[2] = reducedxMapped[3*n+2] / reducedCouplingNum[n];
      }
      
      phiprof::stop("fieldtracing-ionosphere-fsgridCoupling");
   }

   /*! Calculate mapping between ionospheric nodes and Vlasov grid cells.
   * Input is the cell coordinate of the vlasov grid cell.
   * To do so, magnetic field lines are traced inwords from the Vlasov grid
   * IONOSPHERE boundary cells to the ionosphere shell.
   *
   * The return value is a pair of nodeID and coupling factor for the three
   * corners of the containing element.
   */
   std::array<std::pair<int, Real>, 3> calculateIonosphereVlasovGridCoupling(
      std::array<Real,3> x,
      std::vector<SBC::SphericalTriGrid::Node> & nodes,
      creal couplingRadius
   ) {
      
      std::array<std::pair<int, Real>, 3> coupling;
      
      Real stepSize = 100e3;
      std::array<Real,3> v;
      phiprof::start("fieldtracing-ionosphere-VlasovGridCoupling");
      
      // For tracing towards the vlasov boundary, we only require the dipole field.
      TracingFieldFunction dipoleFieldOnly = [](std::array<Real,3>& r, const bool outwards, std::array<Real,3>& b)->bool {
         
         // Get field direction
         b[0] = SBC::ionosphereGrid.dipoleField(r[0],r[1],r[2],X,0,X) + SBC::ionosphereGrid.BGB[0];
         b[1] = SBC::ionosphereGrid.dipoleField(r[0],r[1],r[2],Y,0,Y) + SBC::ionosphereGrid.BGB[1];
         b[2] = SBC::ionosphereGrid.dipoleField(r[0],r[1],r[2],Z,0,Z) + SBC::ionosphereGrid.BGB[2];
         
         // Normalize
         Real  norm = 1. / sqrt(b[0]*b[0] + b[1]*b[1] + b[2]*b[2]);
         for(int c=0; c<3; c++) {
            b[c] = b[c] * norm;
         }
         
         // Make sure motion is outwards. Flip b if dot(r,b) < 0
         if(outwards) {
            if(b[0]*r[0] + b[1]*r[1] + b[2]*r[2] < 0) {
               b[0]*=-1;
               b[1]*=-1;
               b[2]*=-1;
            }
         } else {
            if(b[0]*r[0] + b[1]*r[1] + b[2]*r[2] > 0) {
               b[0]*=-1;
               b[1]*=-1;
               b[2]*=-1;
            }
         }
         return true;
      };
      
      while(sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]) > SBC::Ionosphere::innerRadius) {
         
         // Make one step along the fieldline
         stepFieldLine(x,v, stepSize,50e3,100e3,fieldTracingParameters.tracingMethod,dipoleFieldOnly,false);
         
         // If the field lines is moving even further outwards, abort.
         // (this shouldn't happen under normal magnetospheric conditions, but who
         // knows what crazy driving this will be run with)
         if(sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]) > 1.5*couplingRadius) {
            cerr << "(fieldtracing) Warning: coupling of Vlasov grid cell failed due to weird magnetic field topology." << endl;
            
            // Return a coupling that has 0 value and results in zero potential
            phiprof::stop("fieldtracing-ionosphere-VlasovGridCoupling");
            return coupling;
         }
      }
      
      // Determine the nearest ionosphere node to this point.
      uint32_t nearestNode = SBC::ionosphereGrid.findNodeAtCoordinates(x);
      int32_t elementIndex = nodes[nearestNode].touchingElements[0];
      int32_t oldElementIndex;
      
      std::unordered_set<int32_t> elementHistory;
      bool override=false;
      
      for (uint toto=0; toto<15; toto++) {
         const SBC::SphericalTriGrid::Element& el = SBC::ionosphereGrid.elements[elementIndex];
         oldElementIndex = elementIndex;
         
         if(elementHistory.find(elementIndex) == elementHistory.end()) {
            elementHistory.insert(elementIndex);
         } else {
            // This element was already seen, entering a loop, let's get out
            // It happens when the projection rx is left seen from the right triangle and right seen from the left triangle.
            cerr << "Entered a loop, taking the current element " << elementIndex << "." << endl;
            override=true;
         }
         
         // Calculate barycentric coordinates for x in this element.
         Vec3d r1(nodes[el.corners[0]].x.data());
         Vec3d r2(nodes[el.corners[1]].x.data());
         Vec3d r3(nodes[el.corners[2]].x.data());
         Vec3d rx(x[0],x[1],x[2]);
         
         cint handedness = sign(dot_product(cross_product(r2-r1, r3-r1), r1));
         
         creal kappa1 = handedness*sign(dot_product(cross_product(r1, r2-r1), rx-r1));
         creal kappa2 = handedness*sign(dot_product(cross_product(r2, r3-r2), rx-r2));
         creal kappa3 = handedness*sign(dot_product(cross_product(r3, r1-r3), rx-r3));
         
         if(override || (kappa1 > 0 && kappa2 > 0 && kappa3 > 0)) {
            // Total area
            Real A = vector_length(cross_product(r2-r1,r3-r1));
            
            // Project x into the plane of this triangle
            Vec3d normal = normalize_vector(cross_product(r2-r1, r3-r1));
            rx -= normal*dot_product(rx-r1, normal);
            
            // Area of the sub-triangles
            Real lambda1 = vector_length(cross_product(r2-rx, r3-rx)) / A;
            Real lambda2 = vector_length(cross_product(r1-rx, r3-rx)) / A;
            Real lambda3 = vector_length(cross_product(r1-rx, r2-rx)) / A;
            
            coupling[0] = {el.corners[0], lambda1};
            coupling[1] = {el.corners[1], lambda2};
            coupling[2] = {el.corners[2], lambda3};
            phiprof::stop("fieldtracing-ionosphere-VlasovGridCoupling");
            return coupling;
         } else if (kappa1 > 0 && kappa2 > 0 && kappa3 < 0) {
            elementIndex = SBC::ionosphereGrid.findElementNeighbour(elementIndex,0,2);
         } else if (kappa2 > 0 && kappa3 > 0 && kappa1 < 0) {
            elementIndex = SBC::ionosphereGrid.findElementNeighbour(elementIndex,0,1);
         } else if (kappa3 > 0 && kappa1 > 0 && kappa2 < 0) {
            elementIndex = SBC::ionosphereGrid.findElementNeighbour(elementIndex,1,2);
         } else if (kappa1 < 0 && kappa2 < 0 && kappa3 > 0) {
            if (handedness > 0) {
               elementIndex = SBC::ionosphereGrid.findElementNeighbour(elementIndex,0,1);
            } else {
               elementIndex = SBC::ionosphereGrid.findElementNeighbour(elementIndex,1,2);
            }
         } else if (kappa1 < 0 && kappa2 > 0 && kappa3 < 0) {
            if (handedness > 0) {
               elementIndex = SBC::ionosphereGrid.findElementNeighbour(elementIndex,0,2);
            } else {
               elementIndex = SBC::ionosphereGrid.findElementNeighbour(elementIndex,0,1);
            }
         } else if (kappa1 > 0 && kappa2 < 0 && kappa3 < 0) {
            if (handedness > 0) {
               elementIndex = SBC::ionosphereGrid.findElementNeighbour(elementIndex,1,2);
            } else {
               elementIndex = SBC::ionosphereGrid.findElementNeighbour(elementIndex,0,2);
            }
         }  else {
            cerr << "This fell through, strange."
            << " kappas " << kappa1 << " " << kappa2 << " " << kappa3
            << " handedness " << handedness
            << " r1 " << r1[0] << " " << r1[1] << " " << r1[2]
            << " r2 " << r2[0] << " " << r2[1] << " " << r2[2]
            << " r3 " << r3[0] << " " << r3[1] << " " << r3[2]
            << " rx " << rx[0] << " " << rx[1] << " " << rx[2]
            << endl;
         }
         if(elementIndex == -1) {
            cerr << __FILE__ << ":" << __LINE__ << ": invalid elementIndex returned for coordinate "
            << x[0] << " " << x[1] << " " << x[2] << " projected to rx " << rx[0] << " " << rx[1] << " " << rx[2]
            << ". Last valid elementIndex: " << oldElementIndex << "." << endl;
            phiprof::stop("ionosphere-VlasovGridCoupling");
            return coupling;
         }
      }
      
      // If we arrived here, we did not find an element to couple to (why?)
      // Return an empty coupling instead
      cerr << "(fieldtracing) Failed to find an ionosphere element to couple to for coordinate " <<
      x[0] << " " << x[1] << " " << x[2] << endl;
      phiprof::stop("fieldtracing-ionosphere-VlasovGridCoupling");
      return coupling;
   }

   /*! Trace magnetic field lines out from ionospheric nodes to record whether they are on an open or closed field line.
   */
   void traceOpenClosedConnection(
      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
      std::vector<SBC::SphericalTriGrid::Node> & nodes
   ) {
      
      // we don't need to do anything if we have no nodes
      if(nodes.size() == 0) {
         return;
      }
      
      phiprof::start("fieldtracing-ionosphere-openclosedTracing");
      // Pick an initial stepsize
      creal stepSize = min(1000e3, technicalGrid.DX / 2.);
      std::vector<Real> nodeTracingStepSize(nodes.size(), stepSize); // In-flight storage of step size, needed when crossing into next MPI domain
      std::vector<Real> reducedNodeTracingStepSize(nodes.size());
      std::array<int, 3> gridSize = technicalGrid.getGlobalSize();
      uint64_t maxTracingSteps = 8 * (gridSize[0] * technicalGrid.DX + gridSize[1] * technicalGrid.DY + gridSize[2] * technicalGrid.DZ) / stepSize;
      
      std::vector<int> nodeMapping(nodes.size(), TracingLineEndType::UNPROCESSED);                                 /*!< For reduction of node coupling */
      std::vector<uint64_t> nodeStepCounter(nodes.size(), 0);                                 /*!< Count number of field line tracing steps */
      std::vector<int> nodeNeedsContinuedTracing(nodes.size(), 1);                    /*!< Flag, whether tracing needs to continue on another task */
      std::vector<std::array<Real, 3>> nodeTracingCoordinates(nodes.size());          /*!< In-flight node upmapping coordinates (for global reduction) */
      
      std::vector<int> nodeTracingStepCount(nodes.size());
      
      for(uint n=0; n<nodes.size(); n++) {
         nodeTracingCoordinates.at(n) = nodes.at(n).x;
      }
      bool anyNodeNeedsTracing;
      
      TracingFieldFunction tracingFullField = [&perBGrid, &dPerBGrid, &technicalGrid](std::array<Real,3>& r, const bool alongB, std::array<Real,3>& b)->bool{
         return traceFullFieldFunction(perBGrid, dPerBGrid, technicalGrid, r, alongB, b);
      };
      
      int itCount=0;
      bool warnMaxStepsExceeded = false;
      do {
         anyNodeNeedsTracing = false;
         itCount++;
         
         #pragma omp parallel
         {
            // Trace node coordinates outwards until a non-sysboundary cell is encountered or the local fsgrid domain has been left.
            #pragma omp for schedule(dynamic)
            for(uint n=0; n<nodes.size(); n++) {
               
               if(!nodeNeedsContinuedTracing[n]) {
                  // This node has already found its target, no need for us to do anything about it.
                  continue;
               }
               SBC::SphericalTriGrid::Node& no = nodes[n];
               
               std::array<Real, 3> x = nodeTracingCoordinates[n];
               std::array<Real, 3> v({0,0,0});
               
               while( true ) {
                  nodeStepCounter[n]++;
                  
                  // Check if the current coordinates (pre-step) are in our own domain.
                  std::array<int, 3> fsgridCell = getLocalFsGridCellIndexForCoord(technicalGrid,x);
                  // If it is not in our domain, somebody else takes care of it.
                  if(fsgridCell[0] == -1) {
                     nodeNeedsContinuedTracing[n] = 0;
                     nodeTracingCoordinates[n] = {0,0,0};
                     nodeTracingStepSize[n]=0;
                     break;
                  }
                  
                  if(nodeStepCounter[n] > maxTracingSteps) {
                     nodeNeedsContinuedTracing[n] = 0;
                     nodeTracingCoordinates[n] = {0,0,0};
                     #pragma omp critical
                     {
                        warnMaxStepsExceeded = true;
                     }
                     break;
                  }
                  
                  // Make one step along the fieldline
                  // If the node is in the North, trace along -B (false for last argument), in the South, trace along B
                  stepFieldLine(x,v, nodeTracingStepSize[n],fieldTracingParameters.min_tracer_dx,technicalGrid.DX/2,fieldTracingParameters.tracingMethod,tracingFullField,(no.x[2] < 0));
                  nodeTracingStepCount[n]++;
                  
                  // Look up the fsgrid cell belonging to these coordinates
                  fsgridCell = getLocalFsGridCellIndexForCoord(technicalGrid,x);
                  std::array<Real, 3> interpolationFactor=getFractionalFsGridCellForCoord(technicalGrid,x);
                  
                  // If we map into the ionosphere, this node is on a closed field line.
                  if(sqrt(x.at(0)*x.at(0) + x.at(1)*x.at(1) + x.at(2)*x.at(2)) < SBC::Ionosphere::innerRadius) {
                     nodeNeedsContinuedTracing[n] = 0;
                     nodeTracingCoordinates[n] = {0,0,0};
                     nodeMapping[n] = TracingLineEndType::CLOSED;
                     break;
                  }
                  
                  // If we map out of the box, this node is on an open field line.
                  if(   x[0] > P::xmax - 4*P::dx_ini
                     || x[0] < P::xmin + 4*P::dx_ini
                     || x[1] > P::ymax - 4*P::dy_ini
                     || x[1] < P::ymin + 4*P::dy_ini
                     || x[2] > P::zmax - 4*P::dz_ini
                     || x[2] < P::zmin + 4*P::dz_ini
                  ) {
                     nodeNeedsContinuedTracing[n] = 0;
                     nodeTracingCoordinates[n] = {0,0,0};
                     nodeMapping[n] = TracingLineEndType::OPEN;
                     break;
                  }
                  
                  // Now, after stepping, if it is no longer in our domain, another MPI rank will pick up later.
                  if(fsgridCell[0] == -1) {
                     nodeNeedsContinuedTracing[n] = 1;
                     nodeTracingCoordinates[n] = x;
                     break;
                  }
               }
            } // pragma omp parallel
         }
         
         // Globally reduce whether any node still needs to be picked up and traced onwards
         std::vector<int> sumNodeNeedsContinuedTracing(nodes.size(), 0);
         std::vector<std::array<Real, 3>> sumNodeTracingCoordinates(nodes.size());
         std::vector<uint64_t> maxNodeStepCounter(nodes.size(), 0);
         MPI_Allreduce(nodeNeedsContinuedTracing.data(), sumNodeNeedsContinuedTracing.data(), nodes.size(), MPI_INT, MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(nodeStepCounter.data(), maxNodeStepCounter.data(), nodes.size(), MPI_UINT64_T, MPI_MAX, MPI_COMM_WORLD);
         if(sizeof(Real) == sizeof(double)) {
            MPI_Allreduce(nodeTracingCoordinates.data(), sumNodeTracingCoordinates.data(), 3*nodes.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(nodeTracingStepSize.data(), reducedNodeTracingStepSize.data(), nodes.size(), MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
         } else {
            MPI_Allreduce(nodeTracingCoordinates.data(), sumNodeTracingCoordinates.data(), 3*nodes.size(), MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(nodeTracingStepSize.data(), reducedNodeTracingStepSize.data(), nodes.size(), MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
         }
         for(uint n=0; n<nodes.size(); n++) {
            if(sumNodeNeedsContinuedTracing[n] > 0) {
               anyNodeNeedsTracing=true;
               nodeNeedsContinuedTracing[n] = 1;
               
               // Update that nodes' tracing coordinates
               nodeTracingCoordinates[n][0] = sumNodeTracingCoordinates[n][0] / sumNodeNeedsContinuedTracing[n];
               nodeTracingCoordinates[n][1] = sumNodeTracingCoordinates[n][1] / sumNodeNeedsContinuedTracing[n];
               nodeTracingCoordinates[n][2] = sumNodeTracingCoordinates[n][2] / sumNodeNeedsContinuedTracing[n];
               
               nodeStepCounter[n] = maxNodeStepCounter[n];
            }
            nodeTracingStepSize[n] = reducedNodeTracingStepSize[n];
         }
      } while(anyNodeNeedsTracing);
      
      logFile << "(fieldtracing) open-closed tracing traced in " << itCount << " iterations of the tracing loop." << endl;
      
      bool redWarning = false;
      MPI_Allreduce(&warnMaxStepsExceeded, &redWarning, 1, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD);
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      if(redWarning && rank == MASTER_RANK) {
         logFile << "(fieldtracing) Warning: reached the maximum number of tracing steps " << maxTracingSteps << " allowed for open-closed ionosphere tracing." << endl;
      }
      
      std::vector<int> reducedNodeMapping(nodes.size());
      std::vector<int> reducedNodeTracingStepCount(nodes.size());
      MPI_Allreduce(nodeMapping.data(), reducedNodeMapping.data(), nodes.size(), MPI_INT, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(nodeTracingStepCount.data(), reducedNodeTracingStepCount.data(), nodes.size(), MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      
      for(uint n=0; n<nodes.size(); n++) {
         nodes[n].openFieldLine = reducedNodeMapping.at(n);
      }
      
      phiprof::stop("fieldtracing-ionosphere-openclosedTracing");
   }

   /*! Trace magnetic field lines forward and backward from each DCCRG cell to record the connectivity.
   */
   void traceFullBoxConnection(
      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
      dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid
   ) {
      phiprof::start("fieldtracing-fullTracing");
      
      std::vector<CellID> localDccrgCells = getLocalCells();
      int localDccrgSize = localDccrgCells.size();
      int globalDccrgSize;
      MPI_Allreduce(&localDccrgSize, &globalDccrgSize, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      int commSize;
      MPI_Comm_size(MPI_COMM_WORLD, &commSize);
      std::vector<int> amounts(commSize);
      std::vector<int> displacements(commSize);
      std::vector<CellID> allDccrgCells(globalDccrgSize);
      MPI_Allgather(&localDccrgSize, 1, MPI_INT, amounts.data(), 1, MPI_INT, MPI_COMM_WORLD);
      for(int i=1; i<commSize; i++) {
         displacements[i] = displacements[i-1] + amounts[i-1];
      }
      MPI_Allgatherv(localDccrgCells.data(), localDccrgSize, MPI_UINT64_T, allDccrgCells.data(), amounts.data(), displacements.data(), MPI_UINT64_T, MPI_COMM_WORLD);
      
      // Pick an initial stepsize
      creal stepSize = min(1000e3, technicalGrid.DX / 2.);
      std::vector<Real> cellFWTracingStepSize(globalDccrgSize, stepSize); // In-flight storage of step size, needed when crossing into next MPI domain
      std::vector<Real> cellBWTracingStepSize(globalDccrgSize, stepSize); // In-flight storage of step size, needed when crossing into next MPI domain
      std::array<int, 3> gridSize = technicalGrid.getGlobalSize();
      uint64_t maxTracingSteps = 4 * (gridSize[0] * technicalGrid.DX + gridSize[1] * technicalGrid.DY + gridSize[2] * technicalGrid.DZ) / stepSize;
      
      std::vector<int> cellFWConnection(globalDccrgSize, TracingLineEndType::UNPROCESSED);                                 /*!< For reduction of node coupling */
      std::vector<int> cellBWConnection(globalDccrgSize, TracingLineEndType::UNPROCESSED);                                 /*!< For reduction of node coupling */
      std::vector<uint64_t> cellFWStepCounter(globalDccrgSize, 0);                                 /*!< Count number of field line tracing steps */
      std::vector<uint64_t> cellBWStepCounter(globalDccrgSize, 0);                                 /*!< Count number of field line tracing steps */
      std::vector<int> cellNeedsContinuedFWTracing(globalDccrgSize, 1);                    /*!< Flag, whether tracing needs to continue on another task */
      std::vector<int> cellNeedsContinuedBWTracing(globalDccrgSize, 1);                    /*!< Flag, whether tracing needs to continue on another task */
      std::vector<std::array<Real, 3>> cellFWTracingCoordinates(globalDccrgSize);          /*!< In-flight node upmapping coordinates (for global reduction) */
      std::vector<std::array<Real, 3>> cellBWTracingCoordinates(globalDccrgSize);          /*!< In-flight node upmapping coordinates (for global reduction) */
      
      // These guys are needed in the reductions at the bottom of the tracing loop.
      std::vector<int> reducedCellNeedsContinuedFWTracing(globalDccrgSize, 0);
      std::vector<int> reducedCellNeedsContinuedBWTracing(globalDccrgSize, 0);
      std::vector<std::array<Real, 3>> sumCellFWTracingCoordinates(globalDccrgSize);
      std::vector<std::array<Real, 3>> sumCellBWTracingCoordinates(globalDccrgSize);
      std::vector<uint64_t> maxCellFWStepCounter(globalDccrgSize, 0);
      std::vector<uint64_t> maxCellBWStepCounter(globalDccrgSize, 0);
      std::vector<Real> reducedCellFWTracingStepSize(globalDccrgSize);
      std::vector<Real> reducedCellBWTracingStepSize(globalDccrgSize);
      
      phiprof::start("first-loop");
      for(int n=0; n<globalDccrgSize; n++) {
         const CellID id = allDccrgCells[n];
         cellFWTracingCoordinates.at(n) = mpiGrid.get_center(id);
         cellBWTracingCoordinates.at(n) = cellFWTracingCoordinates.at(n);
         
         if(mpiGrid.is_local(id)) {
            if((mpiGrid[id]->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY)
               || cellFWTracingCoordinates[n][0] > P::xmax - 4*P::dx_ini
               || cellFWTracingCoordinates[n][0] < P::xmin + 4*P::dx_ini
               || cellFWTracingCoordinates[n][1] > P::ymax - 4*P::dy_ini
               || cellFWTracingCoordinates[n][1] < P::ymin + 4*P::dy_ini
               || cellFWTracingCoordinates[n][2] > P::zmax - 4*P::dz_ini
               || cellFWTracingCoordinates[n][2] < P::zmin + 4*P::dz_ini
            ) {
               cellNeedsContinuedFWTracing[n] = 0;
               cellNeedsContinuedBWTracing[n] = 0;
               cellFWTracingCoordinates[n] = {0,0,0};
               cellBWTracingCoordinates[n] = {0,0,0};
               cellFWTracingStepSize[n] = 0;
               cellBWTracingStepSize[n] = 0;
            }
         }
      }
      phiprof::stop("first-loop");
      MPI_Allreduce(cellNeedsContinuedFWTracing.data(), reducedCellNeedsContinuedFWTracing.data(), globalDccrgSize, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(cellNeedsContinuedBWTracing.data(), reducedCellNeedsContinuedBWTracing.data(), globalDccrgSize, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
      if(sizeof(Real) == sizeof(double)) {
         MPI_Allreduce(cellFWTracingStepSize.data(), reducedCellFWTracingStepSize.data(), globalDccrgSize, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
         MPI_Allreduce(cellBWTracingStepSize.data(), reducedCellBWTracingStepSize.data(), globalDccrgSize, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      } else {
         MPI_Allreduce(cellFWTracingStepSize.data(), reducedCellFWTracingStepSize.data(), globalDccrgSize, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
         MPI_Allreduce(cellBWTracingStepSize.data(), reducedCellBWTracingStepSize.data(), globalDccrgSize, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
      }
      cellNeedsContinuedFWTracing = reducedCellNeedsContinuedFWTracing;
      cellNeedsContinuedBWTracing = reducedCellNeedsContinuedBWTracing;
      cellFWTracingStepSize = reducedCellFWTracingStepSize;
      cellBWTracingStepSize = reducedCellBWTracingStepSize;
      bool anyCellNeedsTracing;
      
      TracingFieldFunction tracingFullField = [&perBGrid, &dPerBGrid, &technicalGrid](std::array<Real,3>& r, const bool alongB, std::array<Real,3>& b)->bool{
         return traceFullFieldFunction(perBGrid, dPerBGrid, technicalGrid, r, alongB, b);
      };
      
      int itCount = 0;
      bool warnMaxStepsExceeded = false;
      int cellsToDo;
      phiprof::start("loop");
      #pragma omp parallel shared(cellsToDo)
      {
         do { // while(anyCellNeedsTracing)
            #pragma omp single
            {
               itCount++;
            }
            // Trace node coordinates forward and backwards until a non-sysboundary cell is encountered or the local fsgrid domain has been left.
            #pragma omp for schedule(dynamic)
            for(int n=0; n<globalDccrgSize; n++) {
               
               if(cellNeedsContinuedFWTracing[n]) {
                  
                  std::array<Real, 3> x = cellFWTracingCoordinates[n];
                  std::array<Real, 3> v({0,0,0});
                  
                  while( true ) {
                     // Check if the current coordinates (pre-step) are in our own domain.
                     std::array<int, 3> fsgridCell = getLocalFsGridCellIndexForCoord(technicalGrid,x);
                     // If it is not in our domain, somebody else takes care of it.
                     if(fsgridCell[0] == -1) {
                        cellNeedsContinuedFWTracing[n] = 0;
                        cellFWTracingCoordinates[n] = {0,0,0};
                        cellFWTracingStepSize[n]=0;
                        break;
                     }
                     
                     if(cellFWStepCounter[n] > maxTracingSteps) {
                        cellNeedsContinuedFWTracing[n] = 0;
                        cellFWTracingCoordinates[n] = {0,0,0};
                        cellFWConnection[n] = TracingLineEndType::LOOP;
                        #pragma omp critical
                        {
                           warnMaxStepsExceeded = true;
                        }
                        break;
                     }
                     
                     cellFWStepCounter[n]++;
                     
                     // Make one step along the fieldline
                     // Forward tracing means true for last argument
                     stepFieldLine(x,v, cellFWTracingStepSize[n],100e3,technicalGrid.DX/2,fieldTracingParameters.tracingMethod,tracingFullField,true);
                     
                     // Look up the fsgrid cell belonging to these coordinates
                     fsgridCell = getLocalFsGridCellIndexForCoord(technicalGrid,x);
                     
                     // If we map into the ionosphere, this node is on a closed field line.
                     if(sqrt(x.at(0)*x.at(0) + x.at(1)*x.at(1) + x.at(2)*x.at(2)) < SBC::Ionosphere::innerRadius) {
                        cellNeedsContinuedFWTracing[n] = 0;
                        cellFWTracingCoordinates[n] = {0,0,0};
                        cellFWConnection[n] = TracingLineEndType::CLOSED;
                        break;
                     }
                     
                     // If we map out of the box, this node is on an open field line.
                     if(
                        x[0] > P::xmax - 4*P::dx_ini
                        || x[0] < P::xmin + 4*P::dx_ini
                        || x[1] > P::ymax - 4*P::dy_ini
                        || x[1] < P::ymin + 4*P::dy_ini
                        || x[2] > P::zmax - 4*P::dz_ini
                        || x[2] < P::zmin + 4*P::dz_ini
                     ) {
                        cellNeedsContinuedFWTracing[n] = 0;
                        cellFWTracingCoordinates[n] = {0,0,0};
                        cellFWConnection[n] = TracingLineEndType::OPEN;
                        break;
                     }
                     
                     // Now, after stepping, if it is no longer in our domain, another MPI rank will pick up later.
                     if(fsgridCell[0] == -1) {
                        cellNeedsContinuedFWTracing[n] = 1;
                        cellFWTracingCoordinates[n] = x;
                        break;
                     }
                  }
               } // if FW
               if(cellNeedsContinuedBWTracing[n]) {
                  
                  std::array<Real, 3> x = cellBWTracingCoordinates[n];
                  std::array<Real, 3> v({0,0,0});
                  
                  while( true ) {
                     // Check if the current coordinates (pre-step) are in our own domain.
                     std::array<int, 3> fsgridCell = getLocalFsGridCellIndexForCoord(technicalGrid,x);
                     // If it is not in our domain, somebody else takes care of it.
                     if(fsgridCell[0] == -1) {
                        cellNeedsContinuedBWTracing[n] = 0;
                        cellBWTracingCoordinates[n] = {0,0,0};
                        cellBWTracingStepSize[n]=0;
                        break;
                     }
                     
                     if(cellBWStepCounter[n] > maxTracingSteps) {
                        cellNeedsContinuedBWTracing[n] = 0;
                        cellBWTracingCoordinates[n] = {0,0,0};
                        cellBWConnection[n] = TracingLineEndType::LOOP;
                        #pragma omp critical
                        {
                           warnMaxStepsExceeded = true;
                        }
                        break;
                     }
                     
                     cellBWStepCounter[n]++;
                     
                     // Make one step along the fieldline
                     // Backward tracing means false for last argument
                     stepFieldLine(x,v, cellBWTracingStepSize[n],fieldTracingParameters.min_tracer_dx,technicalGrid.DX/2,fieldTracingParameters.tracingMethod,tracingFullField,false);
                     
                     // Look up the fsgrid cell belonging to these coordinates
                     fsgridCell = getLocalFsGridCellIndexForCoord(technicalGrid,x);
                     
                     // If we map into the ionosphere, this node is on a closed field line.
                     if(sqrt(x.at(0)*x.at(0) + x.at(1)*x.at(1) + x.at(2)*x.at(2)) < SBC::Ionosphere::innerRadius) {
                        cellNeedsContinuedBWTracing[n] = 0;
                        cellBWTracingCoordinates[n] = {0,0,0};
                        cellBWConnection[n] = TracingLineEndType::CLOSED;
                        break;
                     }
                     
                     // If we map out of the box, this node is on an open field line.
                     if(
                        x[0] > P::xmax - 4*P::dx_ini
                        || x[0] < P::xmin + 4*P::dx_ini
                        || x[1] > P::ymax - 4*P::dy_ini
                        || x[1] < P::ymin + 4*P::dy_ini
                        || x[2] > P::zmax - 4*P::dz_ini
                        || x[2] < P::zmin + 4*P::dz_ini
                     ) {
                        cellNeedsContinuedBWTracing[n] = 0;
                        cellBWTracingCoordinates[n] = {0,0,0};
                        cellBWConnection[n] = TracingLineEndType::OPEN;
                        break;
                     }
                     
                     // Now, after stepping, if it is no longer in our domain, another MPI rank will pick up later.
                     if(fsgridCell[0] == -1) {
                        cellNeedsContinuedBWTracing[n] = 1;
                        cellBWTracingCoordinates[n] = x;
                        break;
                     }
                  }
               } // if BW
            } // for
            
            //          string stringi = to_string(rank) + " arrived at " + (string)(__FILE__) + ":" + to_string(__LINE__) + "\n";
            //          cerr << stringi;
            
            // Globally reduce whether any node still needs to be picked up and traced onwards
            #pragma omp barrier
            phiprof::start("MPI-loop");
            #pragma omp master
            {
               MPI_Allreduce(cellNeedsContinuedFWTracing.data(), reducedCellNeedsContinuedFWTracing.data(), globalDccrgSize, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
               MPI_Allreduce(cellNeedsContinuedBWTracing.data(), reducedCellNeedsContinuedBWTracing.data(), globalDccrgSize, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
               MPI_Allreduce(cellFWStepCounter.data(), maxCellFWStepCounter.data(), globalDccrgSize, MPI_UINT64_T, MPI_MAX, MPI_COMM_WORLD);
               MPI_Allreduce(cellBWStepCounter.data(), maxCellBWStepCounter.data(), globalDccrgSize, MPI_UINT64_T, MPI_MAX, MPI_COMM_WORLD);
               if(sizeof(Real) == sizeof(double)) {
                  MPI_Allreduce(cellFWTracingCoordinates.data(), sumCellFWTracingCoordinates.data(), 3*globalDccrgSize, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                  MPI_Allreduce(cellBWTracingCoordinates.data(), sumCellBWTracingCoordinates.data(), 3*globalDccrgSize, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                  MPI_Allreduce(cellFWTracingStepSize.data(), reducedCellFWTracingStepSize.data(), globalDccrgSize, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
                  MPI_Allreduce(cellBWTracingStepSize.data(), reducedCellBWTracingStepSize.data(), globalDccrgSize, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
               } else {
                  MPI_Allreduce(cellFWTracingCoordinates.data(), sumCellFWTracingCoordinates.data(), 3*globalDccrgSize, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
                  MPI_Allreduce(cellBWTracingCoordinates.data(), sumCellBWTracingCoordinates.data(), 3*globalDccrgSize, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
                  MPI_Allreduce(cellFWTracingStepSize.data(), reducedCellFWTracingStepSize.data(), globalDccrgSize, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
                  MPI_Allreduce(cellBWTracingStepSize.data(), reducedCellBWTracingStepSize.data(), globalDccrgSize, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
               }
               anyCellNeedsTracing = false;
            }
            #pragma omp barrier
            phiprof::stop("MPI-loop");
            #pragma omp single
            {
               cellsToDo = 0;
            }
            #pragma omp for schedule(dynamic) reduction(||:anyCellNeedsTracing) reduction(+:cellsToDo)
            for(int n=0; n<globalDccrgSize; n++) {
               if(reducedCellNeedsContinuedFWTracing[n] > 0) {
                  anyCellNeedsTracing=true;
                  cellNeedsContinuedFWTracing[n] = 1;
                  cellsToDo++;
                  
                  // Update that nodes' tracing coordinates
                  cellFWTracingCoordinates[n][0] = sumCellFWTracingCoordinates[n][0] / reducedCellNeedsContinuedFWTracing[n];
                  cellFWTracingCoordinates[n][1] = sumCellFWTracingCoordinates[n][1] / reducedCellNeedsContinuedFWTracing[n];
                  cellFWTracingCoordinates[n][2] = sumCellFWTracingCoordinates[n][2] / reducedCellNeedsContinuedFWTracing[n];
                  
                  cellFWStepCounter[n] = maxCellFWStepCounter[n];
               }
               if(reducedCellNeedsContinuedBWTracing[n] > 0) {
                  anyCellNeedsTracing=true;
                  cellNeedsContinuedBWTracing[n] = 1;
                  cellsToDo++;
                  
                  // Update that nodes' tracing coordinates
                  cellBWTracingCoordinates[n][0] = sumCellBWTracingCoordinates[n][0] / reducedCellNeedsContinuedBWTracing[n];
                  cellBWTracingCoordinates[n][1] = sumCellBWTracingCoordinates[n][1] / reducedCellNeedsContinuedBWTracing[n];
                  cellBWTracingCoordinates[n][2] = sumCellBWTracingCoordinates[n][2] / reducedCellNeedsContinuedBWTracing[n];
                  
                  cellBWStepCounter[n] = maxCellBWStepCounter[n];
               }
               cellFWTracingStepSize[n] = reducedCellFWTracingStepSize[n];
               cellBWTracingStepSize[n] = reducedCellBWTracingStepSize[n];
            }
            #pragma omp barrier
         } while(anyCellNeedsTracing && (cellsToDo >= fieldTracingParameters.max_incomplete_lines_fullbox * 2 * globalDccrgSize));
         
         // Last pass to sort cells that would still continue but won't as we exited as TracingLineEndType::LOOP instead of UNPROCESSED.
         #pragma omp for schedule(dynamic)
         for(int n=0; n<globalDccrgSize; n++) {
            if(cellNeedsContinuedFWTracing[n] == 1) {
               cellFWConnection[n] = TracingLineEndType::LOOP;
            }
            if(cellNeedsContinuedBWTracing[n] == 1) {
               cellBWConnection[n] = TracingLineEndType::LOOP;
            }
         }
      } // pragma omp parallel
      phiprof::stop("loop");
      
      logFile << "(fieldtracing) full box tracing traced in " << itCount << " iterations of the tracing loop with " << cellsToDo << " remaining incomplete field lines (total spatial cells " << globalDccrgSize <<  ")." << endl;
      
      bool redWarning = false;
      MPI_Allreduce(&warnMaxStepsExceeded, &redWarning, 1, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD);
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      if(redWarning && rank == MASTER_RANK) {
         logFile << "(fieldtracing) Warning: reached the maximum number of tracing steps " << maxTracingSteps << " allowed for full-box ionosphere tracing." << endl;
      }
      
      std::vector<int> reducedCellFWConnection(globalDccrgSize);
      std::vector<int> reducedCellBWConnection(globalDccrgSize);
      MPI_Allreduce(cellFWConnection.data(), reducedCellFWConnection.data(), globalDccrgSize, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(cellBWConnection.data(), reducedCellBWConnection.data(), globalDccrgSize, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
      phiprof::start("final-loop");
      for(int n=0; n<globalDccrgSize; n++) {
         const CellID id = allDccrgCells.at(n);
         if(mpiGrid.is_local(id)) {
            mpiGrid[id]->parameters[CellParams::CONNECTION] = TracingPointConnectionType::INVALID;
            if (reducedCellFWConnection[n] == TracingLineEndType::CLOSED && reducedCellBWConnection[n] == TracingLineEndType::CLOSED) {
               mpiGrid[id]->parameters[CellParams::CONNECTION] = TracingPointConnectionType::CLOSED_CLOSED;
            }
            if (reducedCellFWConnection[n] == TracingLineEndType::CLOSED && reducedCellBWConnection[n] == TracingLineEndType::OPEN) {
               mpiGrid[id]->parameters[CellParams::CONNECTION] = TracingPointConnectionType::CLOSED_OPEN;
            }
            if (reducedCellFWConnection[n] == TracingLineEndType::OPEN && reducedCellBWConnection[n] == TracingLineEndType::CLOSED) {
               mpiGrid[id]->parameters[CellParams::CONNECTION] = TracingPointConnectionType::OPEN_CLOSED;
            }
            if (reducedCellFWConnection[n] == TracingLineEndType::OPEN && reducedCellBWConnection[n] == TracingLineEndType::OPEN) {
               mpiGrid[id]->parameters[CellParams::CONNECTION] = TracingPointConnectionType::OPEN_OPEN;
            }
            if (reducedCellFWConnection[n] == TracingLineEndType::CLOSED && reducedCellBWConnection[n] == TracingLineEndType::LOOP) {
               mpiGrid[id]->parameters[CellParams::CONNECTION] = TracingPointConnectionType::CLOSED_LOOP;
            }
            if (reducedCellFWConnection[n] == TracingLineEndType::LOOP && reducedCellBWConnection[n] == TracingLineEndType::CLOSED) {
               mpiGrid[id]->parameters[CellParams::CONNECTION] = TracingPointConnectionType::LOOP_CLOSED;
            }
            if (reducedCellFWConnection[n] == TracingLineEndType::OPEN && reducedCellBWConnection[n] == TracingLineEndType::LOOP) {
               mpiGrid[id]->parameters[CellParams::CONNECTION] = TracingPointConnectionType::OPEN_LOOP;
            }
            if (reducedCellFWConnection[n] == TracingLineEndType::LOOP && reducedCellBWConnection[n] == TracingLineEndType::OPEN) {
               mpiGrid[id]->parameters[CellParams::CONNECTION] = TracingPointConnectionType::LOOP_OPEN;
            }
            if (reducedCellFWConnection[n] == TracingLineEndType::LOOP && reducedCellBWConnection[n] == TracingLineEndType::LOOP) {
               mpiGrid[id]->parameters[CellParams::CONNECTION] = TracingPointConnectionType::LOOP_LOOP;
            }
         }
      }
      phiprof::stop("final-loop");
      
      phiprof::stop("fieldtracing-fullTracing");
   }
   
   
   /*! Trace magnetic field lines forward and backward from each DCCRG cell to record the connectivity.
    */
   void traceFluxRopes(
      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
      dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid
   ) {
      phiprof::start("fieldtracing-fluxropeTracing");
      
      std::vector<CellID> localDccrgCells = getLocalCells();
      int localDccrgSize = localDccrgCells.size();
      int globalDccrgSize;
      MPI_Allreduce(&localDccrgSize, &globalDccrgSize, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      int commSize;
      MPI_Comm_size(MPI_COMM_WORLD, &commSize);
      std::vector<int> amounts(commSize);
      std::vector<int> displacements(commSize);
      std::vector<CellID> allDccrgCells(globalDccrgSize);
      MPI_Allgather(&localDccrgSize, 1, MPI_INT, amounts.data(), 1, MPI_INT, MPI_COMM_WORLD);
      for(int i=1; i<commSize; i++) {
         displacements[i] = displacements[i-1] + amounts[i-1];
      }
      MPI_Allgatherv(localDccrgCells.data(), localDccrgSize, MPI_UINT64_T, allDccrgCells.data(), amounts.data(), displacements.data(), MPI_UINT64_T, MPI_COMM_WORLD);
      
      // Pick an initial stepsize
      creal stepSize = min(1000e3, technicalGrid.DX / 2.);
      std::vector<Real> cellFWTracingStepSize(globalDccrgSize, stepSize); // In-flight storage of step size, needed when crossing into next MPI domain
      std::vector<Real> cellBWTracingStepSize(globalDccrgSize, stepSize); // In-flight storage of step size, needed when crossing into next MPI domain
      std::array<int, 3> gridSize = technicalGrid.getGlobalSize();
      
      std::vector<Real> cellCurvatureRadius(globalDccrgSize);
      std::vector<Real> reducedCellCurvatureRadius(globalDccrgSize);
      std::vector<int> cellNeedsContinuedFWTracing(globalDccrgSize, 1);                    /*!< Flag, whether tracing needs to continue on another task */
      std::vector<int> cellNeedsContinuedBWTracing(globalDccrgSize, 1);                    /*!< Flag, whether tracing needs to continue on another task */
      std::vector<std::array<Real, 3>> cellFWTracingCoordinates(globalDccrgSize);          /*!< In-flight node upmapping coordinates (for global reduction) */
      std::vector<std::array<Real, 3>> cellBWTracingCoordinates(globalDccrgSize);          /*!< In-flight node upmapping coordinates (for global reduction) */
      std::vector<Real> cellFWRunningDistance(globalDccrgSize, 0);
      std::vector<Real> cellBWRunningDistance(globalDccrgSize, 0);
      std::vector<Real> cellFWMaxDistance(globalDccrgSize, 0);
      std::vector<Real> cellBWMaxDistance(globalDccrgSize, 0);
      std::vector<std::array<Real, 3>> cellFWMaxCoordinates(globalDccrgSize);
      std::vector<std::array<Real, 3>> cellBWMaxCoordinates(globalDccrgSize);
      
      // These guys are needed in the reductions at the bottom of the tracing loop.
      std::vector<int> reducedCellNeedsContinuedFWTracing(globalDccrgSize, 0);
      std::vector<int> reducedCellNeedsContinuedBWTracing(globalDccrgSize, 0);
      std::vector<std::array<Real, 3>> sumCellFWTracingCoordinates(globalDccrgSize);
      std::vector<std::array<Real, 3>> sumCellBWTracingCoordinates(globalDccrgSize);
      std::vector<Real> reducedCellFWRunningDistance(globalDccrgSize, 0);
      std::vector<Real> reducedCellBWRunningDistance(globalDccrgSize, 0);
      std::vector<Real> reducedCellFWMaxDistance(globalDccrgSize, 0);
      std::vector<Real> reducedCellBWMaxDistance(globalDccrgSize, 0);
      std::vector<Real> reducedCellFWTracingStepSize(globalDccrgSize);
      std::vector<Real> reducedCellBWTracingStepSize(globalDccrgSize);
      std::vector<std::array<Real, 3>> reducedCellFWMaxCoordinates(globalDccrgSize);
      std::vector<std::array<Real, 3>> reducedCellBWMaxCoordinates(globalDccrgSize);
      
      phiprof::start("first-loop");
      for(int n=0; n<globalDccrgSize; n++) {
         const CellID id = allDccrgCells[n];
         cellFWTracingCoordinates.at(n) = mpiGrid.get_center(id);
         cellBWTracingCoordinates.at(n) = cellFWTracingCoordinates.at(n);
         
         if(mpiGrid.is_local(id)) {
            if((mpiGrid[id]->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY)
               || cellFWTracingCoordinates[n][0] > P::xmax - 4*P::dx_ini
               || cellFWTracingCoordinates[n][0] < P::xmin + 4*P::dx_ini
               || cellFWTracingCoordinates[n][1] > P::ymax - 4*P::dy_ini
               || cellFWTracingCoordinates[n][1] < P::ymin + 4*P::dy_ini
               || cellFWTracingCoordinates[n][2] > P::zmax - 4*P::dz_ini
               || cellFWTracingCoordinates[n][2] < P::zmin + 4*P::dz_ini
            ) {
               cellNeedsContinuedFWTracing[n] = 0;
               cellNeedsContinuedBWTracing[n] = 0;
               cellFWTracingCoordinates[n] = {0,0,0};
               cellBWTracingCoordinates[n] = {0,0,0};
               cellFWTracingStepSize[n] = 0;
               cellBWTracingStepSize[n] = 0;
            } else {
               cellCurvatureRadius[n] = 1 / sqrt(mpiGrid[id]->parameters[CellParams::CURVATUREX]*mpiGrid[id]->parameters[CellParams::CURVATUREX] + mpiGrid[id]->parameters[CellParams::CURVATUREY]*mpiGrid[id]->parameters[CellParams::CURVATUREY] + mpiGrid[id]->parameters[CellParams::CURVATUREZ]*mpiGrid[id]->parameters[CellParams::CURVATUREZ]);
            }
         }
      }
      phiprof::stop("first-loop");
      
      std::vector<std::array<Real,3>> cellInitialCoordinates = cellFWTracingCoordinates;
      
      MPI_Allreduce(cellNeedsContinuedFWTracing.data(), reducedCellNeedsContinuedFWTracing.data(), globalDccrgSize, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(cellNeedsContinuedBWTracing.data(), reducedCellNeedsContinuedBWTracing.data(), globalDccrgSize, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
      if(sizeof(Real) == sizeof(double)) {
         MPI_Allreduce(cellCurvatureRadius.data(), reducedCellCurvatureRadius.data(), globalDccrgSize, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
         MPI_Allreduce(cellFWTracingStepSize.data(), reducedCellFWTracingStepSize.data(), globalDccrgSize, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
         MPI_Allreduce(cellBWTracingStepSize.data(), reducedCellBWTracingStepSize.data(), globalDccrgSize, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      } else {
         MPI_Allreduce(cellCurvatureRadius.data(), reducedCellCurvatureRadius.data(), globalDccrgSize, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
         MPI_Allreduce(cellFWTracingStepSize.data(), reducedCellFWTracingStepSize.data(), globalDccrgSize, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
         MPI_Allreduce(cellBWTracingStepSize.data(), reducedCellBWTracingStepSize.data(), globalDccrgSize, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
      }
      cellCurvatureRadius = reducedCellCurvatureRadius;
      cellNeedsContinuedFWTracing = reducedCellNeedsContinuedFWTracing;
      cellNeedsContinuedBWTracing = reducedCellNeedsContinuedBWTracing;
      cellFWTracingStepSize = reducedCellFWTracingStepSize;
      cellBWTracingStepSize = reducedCellBWTracingStepSize;
      bool anyCellNeedsTracing;
      
      TracingFieldFunction tracingFullField = [&perBGrid, &dPerBGrid, &technicalGrid](std::array<Real,3>& r, const bool alongB, std::array<Real,3>& b)->bool{
         return traceFullFieldFunction(perBGrid, dPerBGrid, technicalGrid, r, alongB, b);
      };
      
      int itCount = 0;
      phiprof::start("loop");
      #pragma omp parallel
      {
         do { // while(anyCellNeedsTracing)
            #pragma omp single
            {
               itCount++;
            }
            // Trace node coordinates forward and backwards until a non-sysboundary cell is encountered or the local fsgrid domain has been left.
            #pragma omp for schedule(dynamic)
            for(int n=0; n<globalDccrgSize; n++) {
               
               if(cellNeedsContinuedFWTracing[n]) {
                  
                  std::array<Real, 3> x = cellFWTracingCoordinates[n];
                  std::array<Real, 3> v({0,0,0});
                  
                  while( true ) {
                     // Check if the current coordinates (pre-step) are in our own domain.
                     std::array<int, 3> fsgridCell = getLocalFsGridCellIndexForCoord(technicalGrid,x);
                     // If it is not in our domain, somebody else takes care of it.
                     if(fsgridCell[0] == -1) {
                        cellNeedsContinuedFWTracing[n] = 0;
                        cellFWTracingCoordinates[n] = {0,0,0};
                        cellFWTracingStepSize[n]=0;
                        break;
                     }
                     
                     if(cellFWRunningDistance[n] > min(fieldTracingParameters.fte_max_curvature_radii_to_trace*cellCurvatureRadius[n],fieldTracingParameters.fte_max_m_to_trace)) {
                        cellNeedsContinuedFWTracing[n] = 0;
                        cellFWTracingCoordinates[n] = {0,0,0};
                        break;
                     }
                     
                     // Make one step along the fieldline
                     // Forward tracing means true for last argument
                     stepFieldLine(x,v, cellFWTracingStepSize[n],100e3,technicalGrid.DX/2,fieldTracingParameters.tracingMethod,tracingFullField,true);
                     
                     cellFWRunningDistance[n] += cellFWTracingStepSize[n];
                     creal distance = sqrt(
                          (x[0]-(cellInitialCoordinates[n])[0])*(x[0]-(cellInitialCoordinates[n])[0])
                        + (x[1]-(cellInitialCoordinates[n])[1])*(x[1]-(cellInitialCoordinates[n])[1])
                        + (x[2]-(cellInitialCoordinates[n])[2])*(x[2]-(cellInitialCoordinates[n])[2])
                     );
                     if(distance > cellFWMaxDistance[n]) {
                        cellFWMaxDistance[n] = distance;
                        cellFWMaxCoordinates[n] = x;
                     }
                     
                     // Look up the fsgrid cell belonging to these coordinates
                     fsgridCell = getLocalFsGridCellIndexForCoord(technicalGrid,x);
                     
                     // If we map into the ionosphere, discard this field line.
                     if(sqrt(x.at(0)*x.at(0) + x.at(1)*x.at(1) + x.at(2)*x.at(2)) < SBC::Ionosphere::innerRadius) {
                        cellNeedsContinuedFWTracing[n] = 0;
                        cellFWTracingCoordinates[n] = {0,0,0};
                        cellFWMaxDistance[n] = fieldTracingParameters.fte_max_m_to_trace;
                        cellFWMaxCoordinates[n] = cellInitialCoordinates[n];
                        break;
                     }
                     
                     // If we map out of the box, discard this field line.
                     if(
                        x[0] > P::xmax - 4*P::dx_ini
                        || x[0] < P::xmin + 4*P::dx_ini
                        || x[1] > P::ymax - 4*P::dy_ini
                        || x[1] < P::ymin + 4*P::dy_ini
                        || x[2] > P::zmax - 4*P::dz_ini
                        || x[2] < P::zmin + 4*P::dz_ini
                     ) {
                        cellNeedsContinuedFWTracing[n] = 0;
                        cellFWTracingCoordinates[n] = {0,0,0};
                        cellFWMaxDistance[n] = fieldTracingParameters.fte_max_m_to_trace;
                        cellFWMaxCoordinates[n] = cellInitialCoordinates[n];
                        break;
                     }
                     
                     // Now, after stepping, if it is no longer in our domain, another MPI rank will pick up later.
                     if(fsgridCell[0] == -1) {
                        cellNeedsContinuedFWTracing[n] = 1;
                        cellFWTracingCoordinates[n] = x;
                        break;
                     }
                  }
               } // if FW
               if(cellNeedsContinuedBWTracing[n]) {
                  
                  std::array<Real, 3> x = cellBWTracingCoordinates[n];
                  std::array<Real, 3> v({0,0,0});
                  
                  while( true ) {
                     // Check if the current coordinates (pre-step) are in our own domain.
                     std::array<int, 3> fsgridCell = getLocalFsGridCellIndexForCoord(technicalGrid,x);
                     // If it is not in our domain, somebody else takes care of it.
                     if(fsgridCell[0] == -1) {
                        cellNeedsContinuedBWTracing[n] = 0;
                        cellBWTracingCoordinates[n] = {0,0,0};
                        cellBWTracingStepSize[n]=0;
                        break;
                     }
                     
                     if(cellBWRunningDistance[n] > min(fieldTracingParameters.fte_max_curvature_radii_to_trace*cellCurvatureRadius[n], fieldTracingParameters.fte_max_m_to_trace)) {
                        cellNeedsContinuedBWTracing[n] = 0;
                        cellBWTracingCoordinates[n] = {0,0,0};
                        break;
                     }
                     
                     // Make one step along the fieldline
                     // Forward tracing means true for last argument
                     stepFieldLine(x,v, cellBWTracingStepSize[n],100e3,technicalGrid.DX/2,fieldTracingParameters.tracingMethod,tracingFullField,false);
                     
                     cellBWRunningDistance[n] += cellBWTracingStepSize[n];
                     creal distance = sqrt(
                          (x[0]-cellInitialCoordinates[n][0])*(x[0]-cellInitialCoordinates[n][0])
                        + (x[1]-cellInitialCoordinates[n][1])*(x[1]-cellInitialCoordinates[n][1])
                        + (x[2]-cellInitialCoordinates[n][2])*(x[2]-cellInitialCoordinates[n][2])
                     );
                     if(distance > cellBWMaxDistance[n]) {
                        cellBWMaxDistance[n] = distance;
                        cellBWMaxCoordinates[n] = x;
                     }
                     
                     // Look up the fsgrid cell belonging to these coordinates
                     fsgridCell = getLocalFsGridCellIndexForCoord(technicalGrid,x);
                     
                     // If we map into the ionosphere, discard this field line.
                     if(sqrt(x.at(0)*x.at(0) + x.at(1)*x.at(1) + x.at(2)*x.at(2)) < SBC::Ionosphere::innerRadius) {
                        cellNeedsContinuedBWTracing[n] = 0;
                        cellBWTracingCoordinates[n] = {0,0,0};
                        cellBWMaxDistance[n] = fieldTracingParameters.fte_max_m_to_trace;
                        cellBWMaxCoordinates[n] = cellInitialCoordinates[n];
                        break;
                     }
                     
                     // If we map out of the box, discard this field line.
                     if(
                        x[0] > P::xmax - 4*P::dx_ini
                        || x[0] < P::xmin + 4*P::dx_ini
                        || x[1] > P::ymax - 4*P::dy_ini
                        || x[1] < P::ymin + 4*P::dy_ini
                        || x[2] > P::zmax - 4*P::dz_ini
                        || x[2] < P::zmin + 4*P::dz_ini
                     ) {
                        cellNeedsContinuedBWTracing[n] = 0;
                        cellBWTracingCoordinates[n] = {0,0,0};
                        cellBWMaxDistance[n] = fieldTracingParameters.fte_max_m_to_trace;
                        cellBWMaxCoordinates[n] = cellInitialCoordinates[n];
                        break;
                     }
                     
                     // Now, after stepping, if it is no longer in our domain, another MPI rank will pick up later.
                     if(fsgridCell[0] == -1) {
                        cellNeedsContinuedBWTracing[n] = 1;
                        cellBWTracingCoordinates[n] = x;
                        break;
                     }
                  }
               } // if BW
            } // for
            
            //          string stringi = to_string(rank) + " arrived at " + (string)(__FILE__) + ":" + to_string(__LINE__) + "\n";
            //          cerr << stringi;
            
            // Globally reduce whether any node still needs to be picked up and traced onwards
            #pragma omp barrier
            phiprof::start("MPI-loop");
            #pragma omp master
            {
               MPI_Allreduce(cellNeedsContinuedFWTracing.data(), reducedCellNeedsContinuedFWTracing.data(), globalDccrgSize, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
               MPI_Allreduce(cellNeedsContinuedBWTracing.data(), reducedCellNeedsContinuedBWTracing.data(), globalDccrgSize, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
               if(sizeof(Real) == sizeof(double)) {
                  MPI_Allreduce(cellFWTracingCoordinates.data(), sumCellFWTracingCoordinates.data(), 3*globalDccrgSize, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                  MPI_Allreduce(cellBWTracingCoordinates.data(), sumCellBWTracingCoordinates.data(), 3*globalDccrgSize, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                  MPI_Allreduce(cellFWMaxCoordinates.data(), reducedCellFWMaxCoordinates.data(), 3*globalDccrgSize, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                  MPI_Allreduce(cellBWMaxCoordinates.data(), reducedCellBWMaxCoordinates.data(), 3*globalDccrgSize, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                  MPI_Allreduce(cellFWTracingStepSize.data(), reducedCellFWTracingStepSize.data(), globalDccrgSize, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
                  MPI_Allreduce(cellBWTracingStepSize.data(), reducedCellBWTracingStepSize.data(), globalDccrgSize, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
                  MPI_Allreduce(cellFWRunningDistance.data(), reducedCellFWRunningDistance.data(), globalDccrgSize, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
                  MPI_Allreduce(cellBWRunningDistance.data(), reducedCellBWRunningDistance.data(), globalDccrgSize, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
                  MPI_Allreduce(cellFWMaxDistance.data(), reducedCellFWMaxDistance.data(), globalDccrgSize, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
                  MPI_Allreduce(cellBWMaxDistance.data(), reducedCellBWMaxDistance.data(), globalDccrgSize, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
               } else {
                  MPI_Allreduce(cellFWTracingCoordinates.data(), sumCellFWTracingCoordinates.data(), 3*globalDccrgSize, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
                  MPI_Allreduce(cellBWTracingCoordinates.data(), sumCellBWTracingCoordinates.data(), 3*globalDccrgSize, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
                  MPI_Allreduce(cellFWMaxCoordinates.data(), reducedCellFWMaxCoordinates.data(), 3*globalDccrgSize, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
                  MPI_Allreduce(cellBWMaxCoordinates.data(), reducedCellBWMaxCoordinates.data(), 3*globalDccrgSize, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
                  MPI_Allreduce(cellFWTracingStepSize.data(), reducedCellFWTracingStepSize.data(), globalDccrgSize, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
                  MPI_Allreduce(cellBWTracingStepSize.data(), reducedCellBWTracingStepSize.data(), globalDccrgSize, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
                  MPI_Allreduce(cellFWRunningDistance.data(), reducedCellFWRunningDistance.data(), globalDccrgSize, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
                  MPI_Allreduce(cellBWRunningDistance.data(), reducedCellBWRunningDistance.data(), globalDccrgSize, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
                  MPI_Allreduce(cellFWMaxDistance.data(), reducedCellFWMaxDistance.data(), globalDccrgSize, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
                  MPI_Allreduce(cellBWMaxDistance.data(), reducedCellBWMaxDistance.data(), globalDccrgSize, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
               }
               anyCellNeedsTracing = false;
            }
            #pragma omp barrier
            phiprof::stop("MPI-loop");
            #pragma omp for schedule(dynamic) reduction(||:anyCellNeedsTracing)
            for(int n=0; n<globalDccrgSize; n++) {
               if(reducedCellNeedsContinuedFWTracing[n] > 0) {
                  anyCellNeedsTracing=true;
                  cellNeedsContinuedFWTracing[n] = 1;
                  
                  // Update that nodes' tracing coordinates
                  cellFWTracingCoordinates[n][0] = sumCellFWTracingCoordinates[n][0] / reducedCellNeedsContinuedFWTracing[n];
                  cellFWTracingCoordinates[n][1] = sumCellFWTracingCoordinates[n][1] / reducedCellNeedsContinuedFWTracing[n];
                  cellFWTracingCoordinates[n][2] = sumCellFWTracingCoordinates[n][2] / reducedCellNeedsContinuedFWTracing[n];
                  
                  cellFWRunningDistance[n] = reducedCellFWRunningDistance[n];
               }
               if(reducedCellNeedsContinuedBWTracing[n] > 0) {
                  anyCellNeedsTracing=true;
                  cellNeedsContinuedBWTracing[n] = 1;
                  
                  // Update that nodes' tracing coordinates
                  cellBWTracingCoordinates[n][0] = sumCellBWTracingCoordinates[n][0] / reducedCellNeedsContinuedBWTracing[n];
                  cellBWTracingCoordinates[n][1] = sumCellBWTracingCoordinates[n][1] / reducedCellNeedsContinuedBWTracing[n];
                  cellBWTracingCoordinates[n][2] = sumCellBWTracingCoordinates[n][2] / reducedCellNeedsContinuedBWTracing[n];
                  
                  cellBWRunningDistance[n] = reducedCellBWRunningDistance[n];
               }
               cellFWTracingStepSize[n] = reducedCellFWTracingStepSize[n];
               cellBWTracingStepSize[n] = reducedCellBWTracingStepSize[n];
               cellFWMaxDistance[n] = reducedCellFWMaxDistance[n];
               cellBWMaxDistance[n] = reducedCellBWMaxDistance[n];
               cellFWMaxCoordinates[n] = reducedCellFWMaxCoordinates[n];
               cellBWMaxCoordinates[n] = reducedCellBWMaxCoordinates[n];
            }
            #pragma omp barrier
         } while(anyCellNeedsTracing);
         
      } // pragma omp parallel
      phiprof::stop("loop");
      
      logFile << "(fieldtracing) flux rope tracing traced in " << itCount << " iterations of the tracing loop." << endl;
      
      phiprof::start("final-loop");
      for(int n=0; n<globalDccrgSize; n++) {
         const CellID id = allDccrgCells.at(n);
         if(mpiGrid.is_local(id)) {
            mpiGrid[id]->parameters[CellParams::FLUXROPE] = 0;
            if(   reducedCellFWMaxDistance[n] < min(fieldTracingParameters.fte_max_curvature_radii_extent*reducedCellCurvatureRadius[n], fieldTracingParameters.fte_max_m_to_trace)
               && reducedCellBWMaxDistance[n] < min(fieldTracingParameters.fte_max_curvature_radii_extent*reducedCellCurvatureRadius[n], fieldTracingParameters.fte_max_m_to_trace)
            ) {
               mpiGrid[id]->parameters[CellParams::FLUXROPE] = 1;
            }
         }
      }
      phiprof::stop("final-loop");
      
      phiprof::stop("fieldtracing-fluxropeTracing");
   }
   
} // namespace FieldTracing
