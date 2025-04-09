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

#ifndef SHOCK_H
#define SHOCK_H

#include "../../definitions.h"
#include "../project.h"

namespace projects {
   class Shock: public Project {
      public:
         Shock();
         virtual ~Shock();
      
         virtual bool initialize(void) override;
         static void addParameters(void);
         virtual void getParameters(void) override;
         virtual void setProjectBField(
            FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
            FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
            FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid
         ) override;
      
         virtual Realf fillPhaseSpace(spatial_cell::SpatialCell *cell,
                                  const uint popID,
                                  const uint nRequested) const override;
         virtual void calcCellParameters(spatial_cell::SpatialCell* cell,creal& t) override;


         Real BX0;
         Real BY0;
         Real BZ0;
         Real EX0;
         Real VX0;
         Real VY0;
         Real VZ0;
         Real DENSITY;
         Real TEMPERATURE;
         Real magPertAmp;
         Real densityPertAmp;
         Real velocityPertAmp;
         Real maxwCutoff;
         Real SCA_X;
         Real SCA_Y;
         Real Sharp_Y;
   } ; //class Shock
} // namespace projects
#endif
