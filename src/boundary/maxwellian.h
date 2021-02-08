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

#ifndef MAXWELLIAN_H
#define MAXWELLIAN_H

#include "../definitions.h"
#include "../readparameters.h"
#include "../spatial_cell.hpp"
#include "boundarycondition.h"
#include "inflow.h"
#include <vector>

using namespace std;

namespace BC {
/*!\brief Maxwellian is a class applying fixed Maxwellian conditions according to parameters read from an input file.
 *
 * Maxwellian is a class handling cells tagged as boundarytype::MAXWELLIAN by this
 * boundary condition.
 *
 * It applies fixed Maxwellian settings to the boundary cells, the parameters of
 * which are being read from an input file.
 *
 * The class inherits most of its machinery from BoundaryCondition::Inflow.
 * The parameters are more general than for Maxwellian and could be put in
 * BoundaryCondition::Inflow but this way they can have a specific prefix which
 * is needed if several inheriting classes are needed.
 */
class Maxwellian : public Inflow {
public:
   Maxwellian();
   ~Maxwellian() override;

   static void addParameters();
   void getParameters() override;

   string getName() const override;
   uint getIndex() const override;

protected:
   void generateTemplateCell(spatial_cell::SpatialCell &templateCell, Real (&B)[3], int inputDataIndex, creal t);

   Real maxwellianDistribution(const uint popID, creal &rho, creal &T, creal &vx, creal &vy, creal &vz);

   vector<vmesh::GlobalID> findBlocksToInitialize(const uint popID, SpatialCell &cell, creal &rho, creal &T, creal &VX,
                                                  creal &VY, creal &VZ);
};
} // namespace BC

#endif
