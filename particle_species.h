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
#include <array>

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

      int sparseBlockAddWidthV;        /*!< Number of layers of blocks that are kept in velocity space around the blocks with content */
      bool sparse_conserve_mass;       /*!< If true, density is scaled to conserve mass when removing blocks*/
      int  sparseDynamicAlgorithm;     /*!< Type of algorithm used for calculating the dynamic minValue; 0 = none, 1 = linear algorithm based on minValue and rho, 2 = linear algorithm based on minValue and Blocks, (Example linear algorithm: minValue = rho / sparse.dynamicValue * sparse.minValue)*/
      Real sparseDynamicBulkValue1;    /*!< Minimum value for the dynamic algorithm range, so for example if dynamicAlgorithm=1 then for sparse.dynamicMinValue = 1e3, sparse.dynamicMaxValue=1e5, we apply the algorithm to cells for which 1e3<cell.rho<1e5*/
      Real sparseDynamicBulkValue2;    /*!< Maximum value for the dynamic algorithm range, so for example if dynamicAlgorithm=1 then for sparse.dynamicMinValue = 1e3, sparse.dynamicMaxValue=1e5, we apply the algorithm to cells for which 1e3<cell.rho<1e5*/
      Real sparseDynamicMinValue1;     /*!< The minimum value for the minValue*/
      Real sparseDynamicMinValue2;     /*!< The maximum value for the minValue*/

      Real thermalRadius;           /*!< Radius of sphere to split the distribution into thermal and suprathermal. 0 (default in cfg) disables the DRO. */
      std::array<Real, 3> thermalV; /*!< Centre of sphere to split the distribution into thermal and suprathermal. 0 (default in cfg) disables the DRO. */

      Real EnergyDensityLimit1;   /*!< Lower bound for second Energy density bin in units of solar wind ram energy. Default 5. */
      Real EnergyDensityLimit2;   /*!< Lower bound forthird Energy density bin in units of solar wind ram energy. Default 10. */
      Real SolarWindEnergy;       /*!< Solar wind ram energy, used for calculating energy density bins. Default value of 0 attempts to use SolarWindSpeed instead.  */
      Real SolarWindSpeed;        /*!< Solar wind speed, used for calculating energy density bins if solar wind ram energy wasn't given. Default 0. */

      int precipitationNChannels;              /*!< Number of energy channels for precipitation differential flux evaluation. Default 16. */
      Real precipitationEmin;                  /*!< Lowest energy channel (in keV) for precipitation differential flux evaluation. Default 0.1. */
      Real precipitationEmax;                  /*!< Highest energy channel (in keV) for precipitation differential flux evaluation. Default 100. */
      Real precipitationLossConeAngle;         /*!< Fixed loss cone opening angle (in deg) for precipitation differential flux evaluation. Default 10. */

       Species();
       Species(const Species& other);
       ~Species();
    };

} // namespace species

#endif	// PARTICLE_SPECIES_H

