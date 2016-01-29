/* This file is part of Vlasiator.
 * 
 * File:   read_gaussian_population.h
 * Author: sandroos
 *
 * Created on March 19, 2015.
 * 
 * Copyright 2015 Finnish Meteorological Institute
 */

#ifndef READ_GAUSSIAN_POPULATION_H
#define	READ_GAUSSIAN_POPULATION_H

#include <cstdlib>
#include <vector>

#include "../definitions.h"

namespace projects {
   
   struct GaussianPopulation {
      int numberOfPopulations;
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

