/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://vlasiator.fmi.fi/
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

//#pragma once

#ifndef IPSHOCK_H
#define IPSHOCK_H

#include "../../definitions.h"
#include "../projectTriAxisSearch.h"

namespace projects {
   class IPShock: public TriAxisSearch {
      public:
         IPShock();
         virtual ~IPShock();
      
         virtual bool initialize(void);
         static void addParameters(void);
         virtual void getParameters(void);

	 virtual void setCellBackgroundField(spatial_cell::SpatialCell* cell) const;
	 virtual Real calcPhaseSpaceDensity(
					    creal& x, creal& y, creal& z,
					    creal& dx, creal& dy, creal& dz,
					    creal& vx, creal& vy, creal& vz,
					    creal& dvx, creal& dvy, creal& dvz,
					    const int& popID
					    ) const;
      


      protected:
         Real getDistribValue(
            creal& x,creal& y, creal& z,
            creal& vx, creal& vy, creal& vz,
	    creal& dvx, creal& dvy, creal& dvz
         ) const;
	 virtual std::vector<std::array<Real, 3>> getV0(creal x, creal y, creal z) const;
	 //virtual void calcCellParameters(Real* cellParams,creal& t);
	 virtual void calcCellParameters(spatial_cell::SpatialCell* cell, creal& t);
         // Interpolate between up- and downstream quantities
         // based on position
         Real interpolate(Real u, Real d, Real x) const;

         // Upstream bulk values
         Real B0u[3];
         Real V0u[3];
         Real DENSITYu;
         Real TEMPERATUREu;
         Real B0utangential;
         Real V0utangential;

         // Downstream bulk values
         Real B0d[3];
         Real V0d[3];
         Real DENSITYd;
         Real TEMPERATUREd;
         Real B0dtangential;
         Real V0dtangential;

         // Flow direction definitions
	 Real Bucosphi;
	 Real Bdcosphi;
	 Real Vucosphi;
	 Real Vdcosphi;
	 int Byusign;
	 int Bydsign;
	 int Bzusign;
	 int Bzdsign;
	 int Vyusign;
	 int Vydsign;
	 int Vzusign;
	 int Vzdsign;

         Real maxwCutoff;
         uint nSpaceSamples;
         uint nVelocitySamples;
         Real Shockwidth;

  } ; //class IPShock
} // namespace projects

#endif
