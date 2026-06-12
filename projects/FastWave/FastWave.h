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

#ifndef LINEARWAVE_H
#define LINEARWAVE_H

#include "../../definitions.h"
#include "../projectTriAxisSearch.h"

namespace projects {

class FastWave : public TriAxisSearch {
public:
   FastWave();
   virtual ~FastWave();

   virtual bool initialize(void);
   static void addParameters(void);
   virtual void getParameters(void);
   virtual void setProjectBField(FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH>& perBGrid,
                                 FsGrid<std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH>& BgBGrid,
                                 FsGrid<fsgrids::technical, FS_STENCIL_WIDTH>& technicalGrid);
   virtual Real calcPhaseSpaceDensity(creal& x, creal& y, creal& z, creal& dx, creal& dy, creal& dz, creal& vx,
                                      creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz, const uint popID) const;

protected:
   Real getMaxwellian(creal& x, creal& y, creal& z, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz,
                      const uint popID) const;
   Real getBiMaxwellian(creal& x, creal& y, creal& z, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy,
                        creal& dvz, const uint popID) const;
   virtual void calcCellParameters(spatial_cell::SpatialCell* cell, creal& t);
   virtual std::vector<std::array<Real, 3>> getV0(creal x, creal y, creal z, const uint popID) const;
   virtual Realf fillPhaseSpace(spatial_cell::SpatialCell *cell,
                                  const uint popID,
                                  const uint nRequested) const;
   virtual Realf probePhaseSpace(spatial_cell::SpatialCell *cell,
                                 const uint popID,
                                 Real vx_in, Real vy_in, Real vz_in) const;
   Real rho0, p0, B0, rho1, u1, p1, ppar1, pperp1;
   Real A, nRef, BRef, lRef, tRef, uRef, TRef, pRef;
   Real lambda, kwave, alpha, cosalpha, sinalpha;
   bool verbose;
   bool isBiMaxwellian;
}; // class FastWave
} // namespace projects

#endif