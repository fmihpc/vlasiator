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
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * File:   particle_species.h
 * Author: sandroos
 *
 * Created on January 27, 2015.
 */

#ifndef PARTICLE_SPECIES_H
#define	PARTICLE_SPECIES_H

#include <omp.h>
#include "definitions.h"

namespace species {
    
   // Aliases for the maximum dt values stored in spatial_cell::Population::max_dt array.
   enum Dt_Elements {
      MAXRDT,                           /**< Maximum spatial translation dt.*/
      MAXVDT,                           /**< Maximum acceleration dt.*/
      SIZE_DT_ELEMENTS                  /**< Number of elements in array.*/
   };
   
    /** Variables common to a particle species.*/
    struct Species {
       std::string name;               /**< Name of the species.*/
       Real charge;                    /**< Particle species charge, in simulation units.*/
       Real mass;                      /**< Particle species mass, in simulation units.*/
       Real sparseMinValue;            /**< Sparse mesh threshold value for the population.*/
       size_t velocityMesh;            /**< ID of the velocity mesh (parameters) this species uses.*/
       
       Species();
       Species(const Species& other);
       ~Species();
    };

} // namespace species

#endif	// PARTICLE_SPECIES_H

