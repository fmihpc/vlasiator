/*
 This file is part of Vlasiator.
 Copyright 2014, 2015 Finnish Meteorological Institute
 */

#ifndef OBJECT_WRAPPER_H
#define OBJECT_WRAPPER_H

#include <vector>

#include "definitions.h"
#include "item_storage.h"
#include "object_factory.h"
#include "amr_refinement_criteria.h"
#include "particle_species.h"

struct ObjectWrapper {
   ObjectWrapper() { }

   ObjectFactory<amr_ref_criteria::Base> amrVelRefCriteria; /**< Factory for all known AMR refinement criteria.*/
   std::vector<species::Species> particleSpecies;           /**< Parameters for all particle species.*/

 private:
   ObjectWrapper(const ObjectWrapper& ow);
   ObjectWrapper& operator=(const ObjectWrapper& ow);
};

#endif
