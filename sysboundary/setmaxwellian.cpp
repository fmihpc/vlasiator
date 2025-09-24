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

#ifndef SETMAXWELLIAN_H
#define SETMAXWELLIAN_H

#include "../definitions.h"
#include "../readparameters.h"
#include "../spatial_cells/spatial_cell_wrapper.hpp"
#include "inflow.h"
#include "sysboundarycondition.h"
#include <vector>

using namespace std;

namespace SBC {
   /*!\brief Maxwellian is a class applying fixed Maxwellian conditions according to parameters read from an input file.
    *
    * Maxwellian is a class handling cells tagged as sysboundarytype::MAXWELLIAN by this boundary condition.
    *
    * It applies fixed Maxwellian settings to the inflow boundary cells, the parameters of
    * which are being read from an input file.
    *
    */
   class Maxwellian : public Inflow {
   public:
      Maxwellian();
      virtual ~Maxwellian();

      static void addParameters();
      virtual void getParameters() override;

      virtual string getName() const override;
      virtual uint getIndex() const override;

      void generateTemplateCell(spatial_cell::SpatialCell& templateCell, Real (&B)[3], int inputDataIndex, creal t) override;
   };
} // namespace SBC

#endif
