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
#include "mesh_data_container.h"
#include "particle_species.h"
#include "projects/project.h"
#include "velocity_mesh_parameters.h"

struct ObjectWrapper {
   ObjectWrapper() { }

   ObjectFactory<amr_ref_criteria::Base> amrVelRefCriteria; /**< Factory for all known AMR refinement criteria.*/
   mesh::MeshDataContainer meshData;                        /**< Container for user-defined mesh data.*/
   std::vector<species::Species> particleSpecies;           /**< Parameters for all particle species.*/
   projects::Project*                    project;           /**< Simulated project.*/
   std::vector<vmesh::MeshParameters> velocityMeshes;       /**< Parameters for velocity mesh(es).*/

 private:
   ObjectWrapper(const ObjectWrapper& ow);
   ObjectWrapper& operator=(const ObjectWrapper& ow);
};

// Currently defined in vlasiator.cpp
ObjectWrapper& getObjectWrapper();

#endif
