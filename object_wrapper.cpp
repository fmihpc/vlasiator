#include "object_wrapper.h"
#include "CLI11.hpp"
#include "particle_species.h"
#include "projects/project.h"
#include "readparameters.h"
#include "velocity_mesh_parameters.h"
#include <array>
#include <functional>
#include <iostream>
#include <string>

bool ObjectWrapper::addParameters() {
   typedef Readparameters RP;
   // Parameters needed to create particle populations
   //NOTE was this still necessary?
   std::function<void(const std::string)> lambda_fun = [this](std::string s) {
      this->initpop(s);
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD,&rank);

   };
   RP::add_each_lambda("ParticlePopulations", "Name of the simulated particle populations (string)", RP::populations,
                       lambda_fun);

   return true;
}
bool ObjectWrapper::addHelp() {

   typedef Readparameters RP;
   if (RP::helpRequested) { // dummy name for the help message
      RP::populations.push_back("POPULATION");
      this->initpop("POPULATION");
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
   // species::Species initnewSpecies(pop, "PROTON", 1.0, 1.0, 1e-15, meshsize, 1, false, 0, 0, 0, 1, 1, 0.0, thermv, 5.0,
                                  // 10.0, 0.0, 0.0, 16, 0.1, 100.0, 10.0);
   species::Species* newSpecies = new species::Species();
   newSpecies->name=pop;
   newSpecies->velocityMesh = meshsize;
   this->particleSpecies.push_back(newSpecies);

   auto species_i = this->particleSpecies.size();
   // species::Species* newSpecies = &(this->particleSpecies.at(species_i - 1));
   vmesh::getMeshWrapper()->velocityMeshesCreation->push_back(newVMeshinit);
// std::cout << "TEST2" << std::endl;

   auto newVMesh = vmesh::getMeshWrapper()->velocityMeshesCreation->at(species_i - 1);
   // std::cout << "VELMESH SIZE!!=" << newSpecies->velocityMesh << std::endl;
   // std::cout << (*newSpecies).name << " and " << (*newSpecies).mass_units << " " << species_i - 1 << std::endl;


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

   // (*newVMesh).meshLimits[0]=0;
   // (*newVMesh).meshLimits[1]=0;
   // (*newVMesh).meshLimits[2]=0;
   // (*newVMesh).meshLimits[3]=0;
   // (*newVMesh).meshLimits[4]=0;
   // (*newVMesh).meshLimits[5]=0;
   // (*newVMesh).gridLength[0]=1;
   // (*newVMesh).gridLength[1]=1;
   // newVMesh->gridLength[2]=1;
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
   RP::add<std::array<Real,3>>(pop + "_thermal.vx",
           "Center coordinate for the maxwellian distribution. Used for calculating the suprathermal moments.",
           newSpecies->thermalV,thermv);
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
   //
   //   std::vector<std::string> popNames=RP::populations;
   //   if (RP::helpRequested) {
   //      popNames.push_back(std::string("<population>"));
   //   } else {
   //      RP::add("", const std::string &desc, const T &defValue)
   //   }
   //
   //  // Create appropriate subparameters for each population
   //  for(auto& pop : popNames) {
   // }
   //
   return true;
}

bool ObjectWrapper::getPopulationParameters() {
   typedef Readparameters RP;

   // Particle population parameters
   for (unsigned int i = 0; i < getObjectWrapper().particleSpecies.size(); i++) {

      species::Species& species =*getObjectWrapper().particleSpecies[i];
      vmesh::MeshParameters& vMesh = *vmesh::getMeshWrapper()->velocityMeshesCreation->at(i);
      const std::string& pop = species.name;

      // Sanity check name
      if (species.name != vMesh.name) {
         std::cerr << "ParticlePopulation parse error: Name " << species.name << " != " << vMesh.name << std::endl;
         return false;
      }

      // Elementary particle parameters
      // RP::get(pop + "_properties.charge", species.charge);

      species.charge *= physicalconstants::CHARGE;

      // RP::get(pop + "_properties.mass", species.mass);
      // RP::get(pop + "_properties.mass_units", massUnit);

      if (species.mass_units == "PROTON") {
         species.mass *= physicalconstants::MASS_PROTON;
      } else if (species.mass_units == "ELECTRON") {
         species.mass *= physicalconstants::MASS_ELECTRON;
      } else {
         std::cerr << "Invalid mass unit for species " << pop << ": '" << species.mass_units << "'" << std::endl;
         return false;
      }

      // sparsity parameters
      // RP::get(pop + "_sparse.minValue", species.sparseMinValue);
      // RP::get(pop + "_sparse.blockAddWidthV", species.sparseBlockAddWidthV);
      // RP::get(pop + "_sparse.conserve_mass", species.sparse_conserve_mass);
      // RP::get(pop + "_sparse.dynamicAlgorithm", species.sparseDynamicAlgorithm);
      // RP::get(pop + "_sparse.dynamicBulkValue1", species.sparseDynamicBulkValue1);
      // RP::get(pop + "_sparse.dynamicBulkValue2", species.sparseDynamicBulkValue2);
      // RP::get(pop + "_sparse.dynamicMinValue1", species.sparseDynamicMinValue1);
      // RP::get(pop + "_sparse.dynamicMinValue2", species.sparseDynamicMinValue2);

      // Particle velocity space properties
      // RP::get(pop + "_vspace.vx_min",vMesh.meshLimits[0]);
      // RP::get(pop + "_vspace.vx_max",vMesh.meshLimits[1]);
      // RP::get(pop + "_vspace.vy_min",vMesh.meshLimits[2]);
      // RP::get(pop + "_vspace.vy_max",vMesh.meshLimits[3]);
      // RP::get(pop + "_vspace.vz_min",vMesh.meshLimits[4]);
      // RP::get(pop + "_vspace.vz_max",vMesh.meshLimits[5]);
      // RP::get(pop + "_vspace.vx_length",vMesh.gridLength[0]);
      // RP::get(pop + "_vspace.vy_length",vMesh.gridLength[1]);
      // RP::get(pop + "_vspace.vz_length",vMesh.gridLength[2]);
      std::cout << vMesh.gridLength[0] << std::endl;
      if (vMesh.gridLength[0] > MAX_BLOCKS_PER_DIM || vMesh.gridLength[1] > MAX_BLOCKS_PER_DIM ||
          vMesh.gridLength[2] > MAX_BLOCKS_PER_DIM) {

         // Build error message as a string first, so that the cerr output hapens atomically and we don't spam
         // thousands of unreadable lines through each other
         species.mass = 1.5;
         std::string errormsg = "(VSPACE) ERROR: Velocity mesh for population " + species.name +
                                " has too many blocks per dimension. Maximum defined in MAX_BLOCKS_PER_DIM is " +
                                std::to_string(MAX_BLOCKS_PER_DIM) + " " + std::string(__FILE__) + ":" +
                                std::to_string(__LINE__) + "\n";
         std::cerr << errormsg;
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

      ////Get thermal / suprathermal moments parameters
      ////Readparameters::get(pop + "_thermal.radius", species.thermalRadius);
      ////Readparameters::get(pop + "_thermal.vx", species.thermalV[0]);
      ////Readparameters::get(pop + "_thermal.vy", species.thermalV[1]);
      ////Readparameters::get(pop + "_thermal.vz", species.thermalV[2]);
      //
      ////Get energy density parameters
      ////Readparameters::get(pop + "_energydensity.limit1", species.EnergyDensityLimit1);
      ////Readparameters::get(pop + "_energydensity.limit2", species.EnergyDensityLimit2);
      ////Readparameters::get(pop + "_energydensity.solarwindenergy", species.SolarWindEnergy);
      ////Readparameters::get(pop + "_energydensity.solarwindspeed", species.SolarWindSpeed);
      //
      const Real EPSILON = 1.e-25;
      if (species.SolarWindEnergy < EPSILON) {
         // Energy stored internally in SI units
         species.SolarWindEnergy = 0.5 * species.mass * species.SolarWindSpeed * species.SolarWindSpeed;
      } else {
         species.SolarWindEnergy = species.SolarWindEnergy * physicalconstants::CHARGE;
      }

      // // Get precipitation parameters
      // //Readparameters::get(pop + "_precipitation.nChannels", species.precipitationNChannels);
      // //Readparameters::get(pop + "_precipitation.emin", species.precipitationEmin);
      // //Readparameters::get(pop + "_precipitation.emax", species.precipitationEmax);
      // //Readparameters::get(pop + "_precipitation.lossConeAngle", species.precipitationLossConeAngle);
      // Convert from eV to SI units
      species.precipitationEmin = species.precipitationEmin * physicalconstants::CHARGE;
      species.precipitationEmax = species.precipitationEmax * physicalconstants::CHARGE;


   }

   return true;
}
