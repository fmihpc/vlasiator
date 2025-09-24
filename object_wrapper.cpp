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

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#ifndef OBJECT_WRAPPER_H
#define OBJECT_WRAPPER_H

// Forward declarations (only possible for classes for which we don't need to know the size)
namespace spatial_cell {
   class SpatialCell;
}
namespace projects {
   class Project;
}

#include <vector>

#include "definitions.h"
#ifndef USE_GPU
#include "object_factory.h"
#endif

#include "particle_species.h"
#include "projects/project.h"
#include "sysboundary/sysboundary.h"

struct ObjectWrapper {
   ObjectWrapper() {}

   std::vector<species::Species> particleSpecies; /**< Parameters for all particle species.*/
   projects::Project* project;                    /**< Simulated project.*/
   SysBoundary sysBoundaryContainer;              /**< Container for sysboundaries.*/

   bool addParameters();           /**< Add config file parameters for objects held in this wrapper */
   bool addPopulationParameters(); /**< After parsing the names of populations, create parameters for each of them */
   bool getPopulationParameters(); /**< Use parsed config file parameters for objects held in this wrapper */

private:
   ObjectWrapper(const ObjectWrapper& ow);
   ObjectWrapper& operator=(const ObjectWrapper& ow);
};

// Currently defined in vlasiator.cpp
ObjectWrapper& getObjectWrapper();

#endif
