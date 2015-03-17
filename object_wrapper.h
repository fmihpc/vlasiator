/*
 This file is part of Vlasiator.
 Copyright 2013, 2014 Finnish Meteorological Institute
 */

#ifndef OBJECT_WRAPPER_H
#define OBJECT_WRAPPER_H

#include "definitions.h"
#include "item_storage.h"
#include "object_factory.h"
#include "amr_refinement_criteria.h"
#include "projects/project.h"

struct ObjectWrapper {
   ObjectWrapper() { }
   ObjectFactory<amr_ref_criteria::Base> amrVelRefCriteria;
   ObjectFactory<projects::Project>      projects;
   
 private:
   ObjectWrapper(const ObjectWrapper& ow);
   ObjectWrapper& operator=(const ObjectWrapper& ow);
};

// Currently defined in vlasiator.cpp
ObjectWrapper& getObjectWrapper();

#endif
