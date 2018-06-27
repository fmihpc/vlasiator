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
 * File:   read_gaussian_population.h
 * Author: sandroos
 *
 * Created on March 19, 2015.
 */

#ifndef READ_GAUSSIAN_POPULATION_H
#define	READ_GAUSSIAN_POPULATION_H

#include <cstdlib>
#include <vector>

#include "../definitions.h"

namespace projects {
   
   struct GaussianPopulation {
      uint numberOfPopulations;
      std::vector<Real> rho;
      std::vector<Real> rhoPertAbsAmp;
      std::vector<Real> Tx;
      std::vector<Real> Ty;
      std::vector<Real> Tz;
      std::vector<Real> Vx;
      std::vector<Real> Vy;
      std::vector<Real> Vz;
   };

   class ReadGaussianPopulation {
   public:
      bool addParameters(const std::string& prefix);
      bool getParameters(const std::string& prefix,projects::GaussianPopulation& populations);
   };

} // namespace projects

#endif	/* READ_GAUSSIAN_H */

