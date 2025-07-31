/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
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

// clang-format off
#include "fs_common.h"
#include <limits>
#include "ldz_gradpe.hpp"
// clang-format on

#ifdef DEBUG_VLASIATOR
#define DEBUG_FSOLVER
#endif

using namespace std;

/**
   Vlasiator considers electrons an inertia-free charge-neutralizing fluid. To model the effects of electron pressure on
   plasma dynamics (such as the cross-shock potential, an electric field induced by electron pressure changes), we consider the
   electron pressure gradient term in the generalized Ohm's law:

   - nabla dot P_e /(n_e e)

   To model this, we can assume some equation of state for electrons. We now assume that the electron fluid has an isotropic
   pressure tensor (so pressure is a scalar), and that it is governed by a polytropic process, where P * V^k is constant.
   Here P is the pressure, V is the volume of an unit of fluid, and k is a polytropic index describing the process.

   Using this equation, for any position in the simulation, the electron temperature can be evaluated based on the local
   electron number density (assumed to be equal to the ion charge density fsgrids::moments::RHOQ divided by the elementary
   charge) and some normalization or anchoring values.

   The polytropic index can be, e.g.:
   k = 0     : Isobaric process (pressure is constant)
   k = 1     : Isothermal process (temperature is constant)
   k = gamma : Adiabatic process (no energy transfer)

   Here gamma is the ratio of specific heats, that is the heat capacity at constant pressure (C_P) divided by the
   heat capacity at constant volume (C_V). For an ideal monoatomic gas, gamma = 5/3.

   Now because P * V^k is constant, we can evaluate that constant for the whole simulation based on the values at
   some hypothetical anchor point, with values usually declared based on the upstream inflow conditions. These values are
   only used for normalizing the relationship between electron density and electron temperature for the given polytropic
   relation. These anchor point values are provided by the user as config parameters - the anchor point is not an actual
   physical point in the simulation, but just imagined for normalization purposes..

   Example in a cfg-script:
   [fieldsolver]
   ohmGradPeTerm = 1 # active
   electronPTindex = 1.666667 # adiabatic, 3/2
   electronDensity = 1.0e6 # anchor point (here: inflow solar wind) electron density value
   electronTemperature = 0.125e6 # electron temperature associated with the anchor point density value. Here
   # the value is set to 1/4 of the solar wind proton temperature.

   The user may, for example, use inflow solar wind electron values as the anchor point values. A good rule of thumb is
   that electron number density times elementary charge should match the ion charge density in the solar wind. A zeroth-order
   approximation is to set electron temperature to equal that of solar wind ions, but alternatively an observation-based factor
   (e.g. 0.25) may be applied to decrease electron temperatures at the anchor point.

   In derivatives.cpp we first calculate the electron pressure at this anchor point:
   Real Pe_anchor = Parameters::electronTemperature * Parameters::electronDensity * physicalconstants::K_B;
   And then calculate the constant used further on:
   Real Pe_const = Pe_anchor * pow(Parameters::electronDensity, -Parameters::electronPTindex);

   Remembering that P * V^k = P * (n_e)^-k is constant, we can calculate that constant from the anchor point values,
   and use that for solving the electron pressure at any given position as
   P(r) = Pe_const * (n_e(r))^k

   Thus, the gradient of electron pressure is calculated in derivatives.cpp using a standard slope limiter, for each Cartesian direction.

   Note:
   The value stored in e.g. dMomentsGrid.get(i,j,k)->at(fsgrids::dmoments::drhoqdx) is a difference, not a derivative. Thus, in the next step,
   evaluating the electron pressure gradient term in the general Ohm's law as
   E_gradPe(x,y,z) = - dPe/d(x,y,z) / (n_e * e * Delta(x,y,z))
   requires adding EGradPeGrid.DX/DY/DZ in the denominator.

   Similar to the Hall term, we use the Parameters::hallMinimumRhoq value as a lower bound for electron number density in order to
   prevent runaway electric field terms at depletion zones.

 */


/** Calculate the electron pressure gradient term on all given cells.
 * @param sysBoundaries System boundary condition functions.
 */
void calculateGradPeTerm(std::span<std::array<Real, fsgrids::egradpe::N_EGRADPE>> egradpes,
                         fsgrids::constmomentsspan moments,
                         std::span<const std::array<Real, fsgrids::dmoments::N_DMOMENTS>> dmoments,
                         fsgrids::consttechnicalspan technical, const fsgrid::FsStencil& stencil,
                         const auto& gridSpacing, SysBoundary& sysBoundaries) {
#ifdef DEBUG_FSOLVER
   if (stencil.ooo() >= moments.size()) {
      cerr << "Out-of-bounds access in " << __FILE__ << ":" << __LINE__ << endl;
      exit(1);
   }
#endif
   auto& egradpe = egradpes[stencil.ooo()];
   const auto& moment = moments[stencil.ooo()];
   const auto& dmoment = dmoments[stencil.ooo()];
   const auto& tech = technical[stencil.ooo()];

   cuint cellSysBoundaryFlag = tech.sysBoundaryFlag;
   cuint cellSysBoundaryLayer = tech.sysBoundaryLayer;

   if (cellSysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE ||
       cellSysBoundaryFlag == sysboundarytype::OUTER_BOUNDARY_PADDING) {
      return;
   }

   if ((cellSysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY) && (cellSysBoundaryLayer != 1)) {
      sysBoundaries.getSysBoundary(cellSysBoundaryFlag)
          ->fieldSolverBoundaryCondGradPeElectricField(egradpes, stencil, 0);
      sysBoundaries.getSysBoundary(cellSysBoundaryFlag)
          ->fieldSolverBoundaryCondGradPeElectricField(egradpes, stencil, 1);
      sysBoundaries.getSysBoundary(cellSysBoundaryFlag)
          ->fieldSolverBoundaryCondGradPeElectricField(egradpes, stencil, 2);
   } else {
      if (Parameters::ohmGradPeTerm == 0) {
         cerr << __FILE__ << __LINE__
              << "You shouldn't be in a electron pressure gradient term function if Parameters::ohmGradPeTerm == 0."
              << endl;
      } else if (Parameters::ohmGradPeTerm == 1) {
         auto calculateEdgeGradPeTermComponent = [&egradpe, &moment, &dmoment](size_t i, size_t j, Real spacing) {
            const Real rhoq = moment[fsgrids::moments::RHOQ];
            const Real min = Parameters::hallMinimumRhoq;
            const Real max = std::numeric_limits<Real>::max();
            const Real hallRhoq = std::clamp(rhoq, min, max);
            egradpe[i] = -dmoment[j] / (hallRhoq * spacing);
         };
         calculateEdgeGradPeTermComponent(fsgrids::egradpe::EXGRADPE, fsgrids::dmoments::dPedx, gridSpacing[0]);
         calculateEdgeGradPeTermComponent(fsgrids::egradpe::EYGRADPE, fsgrids::dmoments::dPedy, gridSpacing[1]);
         calculateEdgeGradPeTermComponent(fsgrids::egradpe::EZGRADPE, fsgrids::dmoments::dPedz, gridSpacing[2]);
      } else {
         cerr << __FILE__ << ":" << __LINE__ << "You are welcome to code higher-order Hall term correction terms."
              << endl;
      }
   }
}

void calculateGradPeTermSimple(std::span<std::array<Real, fsgrids::egradpe::N_EGRADPE>> egradpe,
                               std::span<std::array<Real, fsgrids::egradpe::N_EGRADPE>> egradpedt2,
                               fsgrids::momentsspan moments,
                               fsgrids::momentsspan momentsdt2,
                               fsgrids::dmomentsspan dmoments,
                               fsgrids::dmomentsspan dmomentsdt2,
                               fsgrids::technicalspan technical, FieldSolverGrid &fsgrid,
                               SysBoundary& sysBoundaries, cint& RKCase) {
   phiprof::Timer gradPeTimer{"Calculate GradPe term"};

   const auto& gridSpacing = fsgrid.getGridSpacing();
   const size_t numCells = fsgrid.getNumCells();

   if (not(RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2)) {
      egradpe = egradpedt2;
      moments = momentsdt2;
      dmoments = dmomentsdt2;
   }

   phiprof::Timer mpiTimer{"EgradPe field update ghosts MPI", {"MPI"}};
   fsgrid.updateGhostCells(dmoments);
   mpiTimer.stop();

   // Calculate GradPe term
   fsgrid.parallel_for([](int timerId) -> phiprof::Timer { return phiprof::Timer{timerId}; },
                       phiprof::initializeTimer("EgradPe compute cells"), technical,
                       [=, &sysBoundaries](const fsgrid::FsStencil& stencil, cuint sysBoundaryFlag, cuint sysBoundaryLayer) {
                          calculateGradPeTerm(egradpe, moments, dmoments, technical, stencil, gridSpacing, sysBoundaries);
                       });

   gradPeTimer.stop(numCells, "Spatial Cells");
}
