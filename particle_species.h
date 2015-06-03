/* This file is part of Vlasiator.
 * Copyright 2015 Finnish Meteorological Institute
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
       
       #warning These are specific to Esail simulations
       Real* inflowCounters;           /**< Inflow counters for this species, three counters per OpenMP thread.*/
       Real* outflowCounters;          /**< Outflow counters for this species, three counters per OpenMP thread.*/
       
       Species();
       Species(const Species& other);
       ~Species();
    };

} // namespace species

#endif	// PARTICLE_SPECIES_H

