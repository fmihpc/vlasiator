#include <cstdlib>
#include <iostream>

#include "particle_species.h"

using namespace std;

species::Species::Species() { }

species::Species::Species(const Species& other) {
   name = other.name;
   charge = other.charge;
   mass = other.mass;   
   sparseMinValue = other.sparseMinValue;
   velocityMesh = other.velocityMesh;
}

species::Species::~Species() { }
