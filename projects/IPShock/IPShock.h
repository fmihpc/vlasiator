/*
This file is part of Vlasiator.

Copyright 2016 Finnish Meteorological Institute

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

         // Downstream bulk values
         Real B0d[3];
         Real V0d[3];
         Real DENSITYd;
         Real TEMPERATUREd;

         Real maxwCutoff;
         uint nSpaceSamples;
         uint nVelocitySamples;
         Real Shockwidth;

         Real BPerturbationAmp;
         Real BPerturbationScale;
         Real BPerturbationOctaves;

  } ; //class IPShock
} // namespace projects

#endif
