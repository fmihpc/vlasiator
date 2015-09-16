#include <cstdlib>
#include <iostream>

#include "particle_species.h"

using namespace std;

species::Species::Species() {
   inflowCounters  = new Real[omp_get_num_threads()*3];
   outflowCounters = new Real[omp_get_num_threads()*3];
}

species::Species::Species(const Species& other) {
   name = other.name;
   charge = other.charge;
   mass = other.mass;   
   sparseMinValue = other.sparseMinValue;
   velocityMesh = other.velocityMesh;
   inflowCounters  = new Real[omp_get_num_threads()*3];
   outflowCounters = new Real[omp_get_num_threads()*3];
   for (int i=0; i<omp_get_num_threads()*3; ++i) inflowCounters[i] = other.inflowCounters[i];
   for (int i=0; i<omp_get_num_threads()*3; ++i) outflowCounters[i] = other.outflowCounters[i];
}

species::Species::~Species() {
   delete [] inflowCounters;
   delete [] outflowCounters;
}

