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
#include "../sysboundary/ionosphere.h"

// Get the (integer valued) global fsgrid cell index (i,j,k) for the magnetic-field traced mapping point that node n is
// associated with
template<class T> std::array<int32_t, 3> getGlobalFsGridCellIndexForCoord(T& grid,const std::array<Real, 3>& x) {
   std::array<int32_t, 3> retval;
   retval[0] = floor((x[0] - grid.physicalGlobalStart[0]) / grid.DX);
   retval[1] = floor((x[1] - grid.physicalGlobalStart[1]) / grid.DY);
   retval[2] = floor((x[2] - grid.physicalGlobalStart[2]) / grid.DZ);
   return retval;
}
// Get the (integer valued) local fsgrid cell index (i,j,k) for the magnetic-field traced mapping point that node n is
// associated with If the cell is not in our local domain, will return {-1,-1,-1}
template<class T> std::array<int32_t, 3> getLocalFsGridCellIndexForCoord(T& grid, const std::array<Real, 3>& x) {
   std::array<int32_t, 3> retval = getGlobalFsGridCellIndexForCoord(grid,x);
   retval = grid.globalToLocal(retval[0], retval[1], retval[2]);
   return retval;
}
// Get the fraction fsgrid cell index for the magnetic-field traced mapping point that node n is associated with.
// Note that these are floating point values between 0 and 1
template<class T> std::array<Real, 3> getFractionalFsGridCellForCoord(T& grid, const std::array<Real, 3>& x) {
   std::array<Real, 3> retval;
   std::array<int, 3> fsgridCell = getGlobalFsGridCellIndexForCoord(grid,x);
   retval[0] = (x[0] - grid.physicalGlobalStart[0]) / grid.DX - fsgridCell[0];
   retval[1] = (x[1] - grid.physicalGlobalStart[1]) / grid.DY - fsgridCell[1];
   retval[2] = (x[2] - grid.physicalGlobalStart[2]) / grid.DZ - fsgridCell[2];
   return retval;
}

namespace FieldTracing {
      
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
      TracingMethod tracingMethod;
      Real max_allowed_error; /*!< Maximum alowed error for the adaptive field line tracing methods */
      uint32_t max_field_tracer_attempts; /*!< Max allowed attempts for the iterative field tracers */
      Real min_tracer_dx; /*!< Min allowed tracer dx to avoid getting bogged down in the archipelago */
      Real max_incomplete_lines_fullbox; /*!< Max allowed fraction of field lines left unfinished before exiting tracing loop */
      std::map< std::array<int, 3>, std::array<Real, Rec::N_REC_COEFFICIENTS> > reconstructionCoefficientsCache; /*!< cache for Balsara reconstruction coefficients */
   };
   
   extern FieldTracingParameters fieldTracingParameters;
   
   /*! Type of field line ending, used to classify the ionospheric nodes and the forward and backward field lines in full.box tracing.
   * CLOSED: ends in the ionosphere
   * OPEN: exits the simulation domain
   * LOOP: has not exited, might keep looping or would exit given enough time/steps
   * UNPROCESSED: cells inside the ionosphere or outside the outer limits that weren't even processed in the first place
   */
   enum TracingLineEndType {
      CLOSED,
      OPEN,
      LOOP,
      UNPROCESSED // Keep last for the reductions to work!
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
      CLOSED_LOOP,
      LOOP_CLOSED,
      OPEN_LOOP,
      LOOP_OPEN,
      LOOP_LOOP,
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
   typedef std::function<bool(std::array<Real,3>&, const bool, std::array<Real, 3>&)> TracingFieldFunction;
   
   void stepFieldLine(
      std::array<Real, 3>& x,
      std::array<Real, 3>& v,
      Real& stepSize,
      creal minStepSize,
      creal maxStepSize,
      TracingMethod method,
      TracingFieldFunction& BFieldFunction,
      const bool outwards=true
   );
   
   /*! function to empty the Balsara reconstruction coefficient cache at a new time step */
   inline void resetReconstructionCoefficientsCache() {
      fieldTracingParameters.reconstructionCoefficientsCache.clear();
   }
   
   /*! Link each ionospheric node to fsgrid cells for coupling */
   void calculateFsgridCoupling(
      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
      std::vector<SBC::SphericalTriGrid::Node> & nodes,
      creal radius
   );

   /*! Find coupled ionosphere mesh node for given location */
   std::array<std::pair<int, Real>, 3> calculateVlasovGridCoupling(
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

   /*! Compute the forward and backward connection of all DCCRG cells, tracing done on fsgrid. */
   void traceFullBoxConnection(
      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
      dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid
   );

   void reduceData(
      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
      dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      std::vector<SBC::SphericalTriGrid::Node> & nodes
   );

} // namespace FieldTracing




#endif
