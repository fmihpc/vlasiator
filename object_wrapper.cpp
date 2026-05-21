#include "object_wrapper.h"
#include <CLI11.hpp>
#include "particle_species.h"
#include "readparameters.h"
#include "velocity_mesh_parameters.h"
#include <array>
#include <iostream>
#include <string>

bool ObjectWrapper::addParameters() {
   typedef Readparameters RP;
   RP::add("ParticlePopulations", "Name of the simulated particle populations (string)", RP::populations);
   return true;
}
bool ObjectWrapper::addHelp() {

   typedef Readparameters RP;
   if (RP::helpRequested) { // dummy name for the help message
      RP::populations.clear();
      RP::populations.push_back("POPULATION");
   }
   return true;
}
void ObjectWrapper::initpop(std::string pop) {

   typedef Readparameters RP;
   vmesh::MeshParameters* newVMeshinit = new vmesh::MeshParameters();
   // Originally, there was support for species and velocity meshes to be separate.
   // This was abandoned, since there wasn't really any use for it.
   newVMeshinit->name = pop;
   size_t meshsize = vmesh::getMeshWrapper()->velocityMeshesCreation->size();
   std::array<Real, 3> thermv = {-500000.0, 0, 0};
   species::Species* newSpecies = new species::Species();
   newSpecies->name=pop;
   newSpecies->velocityMesh = meshsize;
   this->particleSpeciesRead.push_back(newSpecies);

   auto species_i = this->particleSpeciesRead.size();
   vmesh::getMeshWrapper()->velocityMeshesCreation->push_back(newVMeshinit);

   auto newVMesh = vmesh::getMeshWrapper()->velocityMeshesCreation->at(species_i - 1);

   RP::add<Real>(pop + "_properties.charge", "Particle charge, in units of elementary charges (int)", newSpecies->charge,1.0);
   RP::add<string>(pop + "_properties.mass_units",
           "Units in which particle mass is given, either 'PROTON' or 'ELECTRON' (string)", newSpecies->mass_units,"PROTON");
   // currently seems like these params are not read from config, debug try adding lambda?
   RP::add<Real>(pop + "_properties.mass", "Particle mass in given units (float)", newSpecies->mass,1.0);
   // Grid sparsity parameters
   RP::add<Real>(pop + "_sparse.minValue",
           "Minimum value of distribution function in any cell of a velocity block for the block to be considered to "
           "have contents",
           newSpecies->sparseMinValue,1.0e-15);
   RP::add<int>(pop + "_sparse.blockAddWidthV",
           "Number of layers of blocks that are kept in velocity space around the blocks with content",
           newSpecies->sparseBlockAddWidthV,1);
   RP::add<bool>(pop + "_sparse.conserve_mass",
           "If true, then mass is conserved by scaling the dist. func. in the remaining blocks",
           newSpecies->sparse_conserve_mass,false);
   RP::add<int>(pop + "_sparse.dynamicAlgorithm",
           "Type of algorithm used for calculating the dynamic minValue; 0 = none, 1 = linear algorithm based on rho, "
           "2 = linear algorithm based on Blocks, (Example linear algorithm: y = kx+b, where "
           "dynamicMinValue1=k*dynamicBulkValue1 + b, and dynamicMinValue2 = k*dynamicBulkValue2 + b",
           newSpecies->sparseDynamicAlgorithm,0);
   RP::add<Real>(pop + "_sparse.dynamicMinValue1", "The minimum value for the dynamic minValue",
           newSpecies->sparseDynamicMinValue1,1.0);
   RP::add<Real>(pop + "_sparse.dynamicMinValue2", "The maximum value (value 2) for the dynamic minValue",
           newSpecies->sparseDynamicMinValue2,1.0);
   RP::add<Real>(pop + "_sparse.dynamicBulkValue1",
           "Minimum value for the dynamic algorithm range, so for example if dynamicAlgorithm=1 then for "
           "sparse.dynamicBulkValue1 = 1e3, sparse.dynamicBulkValue2=1e5, we apply the algorithm to cells for which "
           "1e3<cell.rho<1e5",
           newSpecies->sparseDynamicBulkValue1,0.0);
   RP::add<Real>(pop + "_sparse.dynamicBulkValue2",
           "Maximum value for the dynamic algorithm range, so for example if dynamicAlgorithm=1 then for "
           "sparse.dynamicBulkValue1 = 1e3, sparse.dynamicBulkValue2=1e5, we apply the algorithm to cells for which "
           "1e3<cell.rho<1e5",
           newSpecies->sparseDynamicBulkValue2,0.0);

   // Grid parameters
   RP::add(pop + "_vspace.vx_min", "Minimum value for velocity mesh vx-coordinates.", newVMesh->meshLimits[0]);
   RP::add(pop + "_vspace.vx_max", "Maximum value for velocity mesh vx-coordinates.", newVMesh->meshLimits[1]);
   RP::add(pop + "_vspace.vy_min", "Minimum value for velocity mesh vy-coordinates.", newVMesh->meshLimits[2]);
   RP::add(pop + "_vspace.vy_max", "Maximum value for velocity mesh vx-coordinates.", newVMesh->meshLimits[3]);
   RP::add(pop + "_vspace.vz_min", "Minimum value for velocity mesh vz-coordinates.", newVMesh->meshLimits[4]);
   RP::add(pop + "_vspace.vz_max", "Maximum value for velocity mesh vx-coordinates.", newVMesh->meshLimits[5]);
   RP::add(pop + "_vspace.vx_length", "Initial number of velocity blocks in vx-direction.", newVMesh->gridLength[0]);
   RP::add(pop + "_vspace.vy_length", "Initial number of velocity blocks in vy-direction.", newVMesh->gridLength[1]);
   RP::add(pop + "_vspace.vz_length", "Initial number of velocity blocks in vz-direction.", newVMesh->gridLength[2]);
   // RP::add(pop + "_vspace.max_refinement_level","Maximum allowed mesh refinement level.", (*newVMesh).meshMinLimits);
   // //Was not even used?

   // Thermal / suprathermal parameters
   RP::add<Real>(pop + "_thermal.vx",
           "Center coordinate for the maxwellian distribution. Used for calculating the suprathermal moments.",
           newSpecies->thermalV[0],thermv.at(0));
   RP::add<Real>(pop + "_thermal.vy",
           "Center coordinate for the maxwellian distribution. Used for calculating the suprathermal moments.",
           newSpecies->thermalV[1],thermv.at(1));
   RP::add<Real>(pop + "_thermal.vz",
           "Center coordinate for the maxwellian distribution. Used for calculating the suprathermal moments.",
           newSpecies->thermalV[2],thermv.at(2));

   RP::add<Real>(pop + "_thermal.radius",
           "Radius of the maxwellian distribution. Used for calculating the suprathermal moments. If set to 0 "
           "(default), the thermal/suprathermal DROs are skipped.",
           newSpecies->thermalRadius,0.0);

   // Precipitation parameters
   RP::add<int>(pop + "_precipitation.nChannels", "Number of energy channels for precipitation differential flux evaluation",
           newSpecies->precipitationNChannels,16);
   RP::add<Real>(pop + "_precipitation.emin", "Lowest energy channel (in eV) for precipitation differential flux evaluation",
           newSpecies->precipitationEmin,0.1);
   RP::add<Real>(pop + "_precipitation.emax", "Highest energy channel (in eV) for precipitation differential flux evaluation",
           newSpecies->precipitationEmax,100.0);
   RP::add<Real>(pop + "_precipitation.lossConeAngle",
           "Fixed loss cone opening angle (in deg) for precipitation differential flux evaluation",
           newSpecies->precipitationLossConeAngle,10.0);

   // Energy density parameters
   RP::add<Real>(pop + "_energydensity.limit1",
           "Lower limit of second bin for energy density, given in units of solar wind ram energy.",
           newSpecies->EnergyDensityLimit1,5.0);
   RP::add<Real>(pop + "_energydensity.limit2",
           "Lower limit of third bin for energy density, given in units of solar wind ram energy.",
           newSpecies->EnergyDensityLimit2,10.0);
   RP::add<Real>(pop + "_energydensity.solarwindspeed",
           "Incoming solar wind velocity magnitude in m/s. Used for calculating energy densities.",
           newSpecies->SolarWindSpeed,0.0);
   RP::add<Real>(pop + "_energydensity.solarwindenergy",
           "Incoming solar wind ram energy in eV. Used for calculating energy densities.",
           newSpecies->SolarWindEnergy,0.0);
}

bool ObjectWrapper::addPopulationParameters() {
  typedef Readparameters RP;
  for (auto pop : RP::populations){
    this->initpop(pop);
  }
  return true;
}

bool ObjectWrapper::getPopulationParameters() {
   typedef Readparameters RP;

   // Particle population parameters
   for (unsigned int i = 0; i < getObjectWrapper().particleSpeciesRead.size(); i++) {

      species::Species& species =*getObjectWrapper().particleSpeciesRead[i];
    
      vmesh::MeshParameters& vMesh = *vmesh::getMeshWrapper()->velocityMeshesCreation->at(i);
      const std::string& pop = species.name;
      // Sanity check name
      if (species.name != vMesh.name) {
         std::cerr << "ParticlePopulation parse error: Name " << species.name << " != " << vMesh.name << std::endl;
         return false;
      }

      // Elementary particle parameters
      species.charge *= physicalconstants::CHARGE;
      if (species.mass_units == "PROTON") {
         species.mass *= physicalconstants::MASS_PROTON;
      } else if (species.mass_units == "ELECTRON") {
         species.mass *= physicalconstants::MASS_ELECTRON;
      } else {
         std::cerr << "Invalid mass unit for species " << pop << ": '" << species.mass_units << "'" << std::endl;
         return false;
      }

      // Particle velocity space properties
      if (vMesh.gridLength[0] > MAX_BLOCKS_PER_DIM || vMesh.gridLength[1] > MAX_BLOCKS_PER_DIM ||
          vMesh.gridLength[2] > MAX_BLOCKS_PER_DIM) {

         // Build error message as a string first, so that the cerr output hapens atomically and we don't spam
         // thousands of unreadable lines through each other
         species.mass = 1.5;
         std::string errormsg = "(VSPACE) ERROR: Velocity mesh for population " + species.name +
                                " has too many blocks per dimension. Maximum defined in MAX_BLOCKS_PER_DIM is " +
                                std::to_string(MAX_BLOCKS_PER_DIM) + " " + std::string(__FILE__) + ":" +
                                std::to_string(__LINE__) + "\n";
         if (!Readparameters::helpRequested) {
            std::cerr << errormsg;
         }
      }

      /* Special handling of WID=8; halve the number of blocks */
      int myRank;
      MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
      if ((WID == 8 && P::adaptGPUWID)) {
         // First verify that we can halve the value2
         if ((vMesh.gridLength[0] % 2 == 0) && (vMesh.gridLength[1] % 2 == 0) && (vMesh.gridLength[2] % 2 == 0)) {
            vMesh.gridLength[0] /= 2;
            vMesh.gridLength[1] /= 2;
            vMesh.gridLength[2] /= 2;
            if (myRank == MASTER_RANK) {
               std::cerr << " Note: Using WID=8; Halving velocity block counts per dimension. Deactivate with "
                            "parameter adaptGPUWID=false."
                         << std::endl;
            }
         } else {
            if (myRank == MASTER_RANK) {
               std::cerr << " Warning: Using WID=8 but odd number of velocity blocks! Cannot halve the blocks count."
                         << std::endl;
            }
         }
      }

      vMesh.blockLength[0] = vMesh.blockLength[1] = vMesh.blockLength[2] = WID;

      const Real EPSILON = 1.e-25;
      if (species.SolarWindEnergy < EPSILON) {
         // Energy stored internally in SI units
         species.SolarWindEnergy = 0.5 * species.mass * species.SolarWindSpeed * species.SolarWindSpeed;
      } else {
         species.SolarWindEnergy = species.SolarWindEnergy*physicalconstants::CHARGE;
      }
      // Convert from eV to SI units
      species.precipitationEmin = species.precipitationEmin*physicalconstants::CHARGE;
      species.precipitationEmax = species.precipitationEmax*physicalconstants::CHARGE;

      getObjectWrapper().particleSpecies.resize(getObjectWrapper().particleSpeciesRead.size());
      getObjectWrapper().particleSpecies.at(i)=species;


   }

   return true;
}
