/* This file is part of Vlasiator.
 * 
 * File:   read_gaussian_population.cpp
 * Author: sandroos
 *
 * Created on March 19, 2015.
 * 
 * Copyright 2015 Finnish Meteorological Institute
 */

#include <cstdlib>
#include <iostream>
#include "../readparameters.h"

#include "read_gaussian_population.h"

using namespace std;

namespace projects {
   
   bool ReadGaussianPopulation::addParameters(const std::string& prefix) {
      // Add input variables to config file reader
      typedef Readparameters RP;
      RP::add(prefix+".n", "Number of populations to use", 0);
      RP::addComposing(prefix+".rho", "Number density (m^-3)");
      RP::addComposing(prefix+".rhoPertAbsAmp", "Absolute amplitude of the density perturbation");
      RP::addComposing(prefix+".Tx", "Temperature (K)");
      RP::addComposing(prefix+".Ty", "Temperature");
      RP::addComposing(prefix+".Tz", "Temperature");
      RP::addComposing(prefix+".Vx", "Bulk velocity x component (m/s)");
      RP::addComposing(prefix+".Vy", "Bulk velocity y component (m/s)");
      RP::addComposing(prefix+".Vz", "Bulk velocity z component (m/s)");
   }

   bool ReadGaussianPopulation::getParameters(const std::string& prefix,projects::GaussianPopulation& populations) {
      // Read values of input variables
      typedef Readparameters RP;
      RP::get(prefix+".n", populations.numberOfPopulations);
      RP::get(prefix+".rho", populations.rho);
      RP::get(prefix+".rhoPertAbsAmp", populations.rhoPertAbsAmp);
      RP::get(prefix+".Tx", populations.Tx);
      RP::get(prefix+".Ty", populations.Ty);
      RP::get(prefix+".Tz", populations.Tz);
      RP::get(prefix+".Vx", populations.Vx);
      RP::get(prefix+".Vy", populations.Vy);
      RP::get(prefix+".Vz", populations.Vz);
      
      // Do some sanity check on input variables
      bool success = true;
      if (populations.numberOfPopulations < 1) {
         cerr << "ERROR, you need to define at least one particle population " << __FILE__ << ":" << __LINE__ << endl;
         success = false;
      }
      if (populations.numberOfPopulations != populations.rho.size()) {
         cerr << "ERROR, number of populations does not match the size of rho in " << __FILE__ << ":" << __LINE__ << endl;
         success = false;
      }
      if (populations.numberOfPopulations != populations.rhoPertAbsAmp.size()) {
         cerr << "ERROR, number of populations does not match the size of rhoPertAbsAmp in " << __FILE__ << ":" << __LINE__ << endl;
         success = false;
      }
      if (populations.numberOfPopulations != populations.Tx.size()) {
         cerr << "ERROR, number of populations does not match the size of Tx in " << __FILE__ << ":" << __LINE__ << endl;
         success = false;
      }
      if (populations.numberOfPopulations != populations.Ty.size()) {
         cerr << "ERROR, number of populations does not match the size of Ty in " << __FILE__ << ":" << __LINE__ << endl;
         success = false;
      }
      if (populations.numberOfPopulations != populations.Tz.size()) {
         cerr << "ERROR, number of populations does not match the size of Tz in " << __FILE__ << ":" << __LINE__ << endl;
         success = false;
      }
      if (populations.numberOfPopulations != populations.Vx.size()) {
         cerr << "ERROR, number of populations does not match the size of Vx in " << __FILE__ << ":" << __LINE__ << endl;
         success = false;
      }
      if (populations.numberOfPopulations != populations.Vy.size()) {
         cerr << "ERROR, number of populations does not match the size of Vy in " << __FILE__ << ":" << __LINE__ << endl;
         success = false;
      }
      if (populations.numberOfPopulations != populations.Vz.size()) {
         cerr << "ERROR, number of populations does not match the size of Vz in " << __FILE__ << ":" << __LINE__ << endl;
         success = false;
      }
      return success;
   }

} // namespace projects