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

#ifndef ALFVENCASCADE_H
#define ALFVENCASCADE_H

#include "../../definitions.h"
#include "../project.h"

namespace projects {

   struct WaveParameters {
      Real wavelength;
      Real amplitude;
      Real phase;
   };

   class AlfvenCascade : public Project {
   public:
      AlfvenCascade();
      virtual ~AlfvenCascade();

      virtual bool initialize(void);
      static void addParameters(void);
      virtual void getParameters(void);
      virtual void setProjectBField(FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH>& perBGrid,
                                    FsGrid<std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH>& BgBGrid,
                                    FsGrid<fsgrids::technical, FS_STENCIL_WIDTH>& technicalGrid);
      virtual Realf fillPhaseSpace(spatial_cell::SpatialCell *cell,
                                    const uint popID,
                                    const uint nRequested) const override;
      virtual void calcCellParameters(spatial_cell::SpatialCell* cell, creal& t);
      // virtual std::vector<std::array<Real, 3>> getV0(creal x, creal y, creal z, const uint popID) const;
      

      // Basic plasma parameters
      Real rho0;    // Background mass density
      Real T;       // Temperature
      Real B;       // Background magnetic field strength
      Real p0;      // Thermal pressure
      Real n0;       // Background number density
      Real VA;      // Alfvén speed
      Real angle;   // Wave vector angle
      Real m;       // Particle mass

      // Turbulence parameters
      int nWaves; // Number of waves in simulation
      std::vector<Real> wavelength;  // Vector of wavelengths
      std::vector<Real> amplitude;   // Vector of velocity amplitudes
      std::vector<Real> phase;       // Vector of initial phases
      
      Real spectralIndex;            // Power law index for initial spectrum
      int randomSeed;               // Seed for random phase generation

      bool verbose;
      std::vector<WaveParameters> waves;
   }; // class AlfvenCascade
} // namespace projects

#endif