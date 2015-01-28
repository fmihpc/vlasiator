/* This file is part of Vlasiator.
 * Copyright 2015 Finnish Meteorological Institute
 * File:   particle_species.h
 * Author: sandroos
 *
 * Created on January 27, 2015.
 */

#ifndef PARTICLE_SPECIES_H
#define	PARTICLE_SPECIES_H

#include "definitions.h"

namespace species {
    
    /** Variables common to a particle species.*/
    struct Species {
        std::string name;               /**< Name of the species.*/
        Real charge;                    /**< Particle species charge, in simulation units.*/
        Real mass;                      /**< Particle species mass, in simulation units.*/
        Real sparseMinValue;            /**< Sparse mesh threshold value for the population.*/
    };

} // namespace species

#endif	// PARTICLE_SPECIES_H

