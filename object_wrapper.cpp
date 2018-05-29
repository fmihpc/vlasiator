#include "object_wrapper.h"
#include "readparameters.h"
#include <iostream>

bool ObjectWrapper::addParameters() {
   typedef Readparameters RP;

   // Parameters needed to create particle populations
   
   if (RP::helpRequested) { // dummy name for the help message
      RP::add("ParticlePopulations","Name of the simulated particle populations (string)","<population>");
   } else {
      RP::addComposing("ParticlePopulations","Name of the simulated particle populations (string)");
   }

   return true;
}

bool ObjectWrapper::addPopulationParameters() {
   typedef Readparameters RP;

   std::vector<std::string> popNames;
   if (RP::helpRequested) {
      popNames.push_back(std::string("<population>"));
   } else {
      RP::get("ParticlePopulations", popNames);
   }

  // Create appropriate subparameters for each population
  for(auto& pop : popNames) {
     species::Species newSpecies;
     vmesh::MeshParameters newVMesh;

     // Originally, there was support for species and velocity meshes to be separate.
     // This was abandoned, since there wasn't really any use for it.
     newSpecies.name = newVMesh.name = pop;
     newSpecies.velocityMesh = getObjectWrapper().velocityMeshes.size();

     getObjectWrapper().particleSpecies.push_back(newSpecies);
     getObjectWrapper().velocityMeshes.push_back(newVMesh);

     RP::add(pop + "_properties.charge", "Particle charge, in units of elementary charges (int)", 1);
     RP::add(pop + "_properties.mass_units", "Units in which particle mass is given, either 'PROTON' or 'ELECTRON' (string)", std::string("PROTON"));
     RP::add(pop + "_properties.mass","Particle mass in given units (float)", 1);
     RP::add(pop + "_properties.test_population","Considered a test particle population? (bool)", 0);
     RP::add(pop + "_properties.propagate_population","Perform spatial translation and acceleration to this species? (bool)", 1);

     // Grid sparsity parameters
     RP::add(pop + "_sparse.minValue", "Minimum value of distribution function in any cell of a velocity block for the block to be considered to have contents", 1e-15);
     RP::add(pop + "_sparse.blockAddWidthV", "Number of layers of blocks that are kept in velocity space around the blocks with content",1);
     RP::add(pop + "_sparse.conserve_mass", "If true, then mass is conserved by scaling the dist. func. in the remaining blocks", false);
     RP::add(pop + "_sparse.dynamicAlgorithm", "Type of algorithm used for calculating the dynamic minValue; 0 = none, 1 = linear algorithm based on rho, 2 = linear algorithm based on Blocks, (Example linear algorithm: y = kx+b, where dynamicMinValue1=k*dynamicBulkValue1 + b, and dynamicMinValue2 = k*dynamicBulkValue2 + b", 0);
     RP::add(pop + "_sparse.dynamicMinValue1", "The minimum value for the dynamic minValue", 1);
     RP::add(pop + "_sparse.dynamicMinValue2", "The maximum value (value 2) for the dynamic minValue", 1);
     RP::add(pop + "_sparse.dynamicBulkValue1", "Minimum value for the dynamic algorithm range, so for example if dynamicAlgorithm=1 then for sparse.dynamicBulkValue1 = 1e3, sparse.dynamicBulkValue2=1e5, we apply the algorithm to cells for which 1e3<cell.rho<1e5", 0);
     RP::add(pop + "_sparse.dynamicBulkValue2", "Maximum value for the dynamic algorithm range, so for example if dynamicAlgorithm=1 then for sparse.dynamicBulkValue1 = 1e3, sparse.dynamicBulkValue2=1e5, we apply the algorithm to cells for which 1e3<cell.rho<1e5", 0);

     // Grid parameters
     RP::add(pop + "_vspace.vx_min","Minimum value for velocity mesh vx-coordinates.",0);
     RP::add(pop + "_vspace.vx_max","Maximum value for velocity mesh vx-coordinates.",0);
     RP::add(pop + "_vspace.vy_min","Minimum value for velocity mesh vy-coordinates.",0);
     RP::add(pop + "_vspace.vy_max","Maximum value for velocity mesh vx-coordinates.",0);
     RP::add(pop + "_vspace.vz_min","Minimum value for velocity mesh vz-coordinates.",0);
     RP::add(pop + "_vspace.vz_max","Maximum value for velocity mesh vx-coordinates.",0);
     RP::add(pop + "_vspace.vx_length","Initial number of velocity blocks in vx-direction.",1);
     RP::add(pop + "_vspace.vy_length","Initial number of velocity blocks in vy-direction.",1);
     RP::add(pop + "_vspace.vz_length","Initial number of velocity blocks in vz-direction.",1);
     RP::add(pop + "_vspace.max_refinement_level","Maximum allowed mesh refinement level.", 1);
     
     // Backstreaming parameters
     Readparameters::add(pop + "_backstream.vx", "Center coordinate for the maxwellian distribution. Used for calculating the backstream moments.", -500000.0);
     Readparameters::add(pop + "_backstream.vy", "Center coordinate for the maxwellian distribution. Used for calculating the backstream moments.", 0.0);
     Readparameters::add(pop + "_backstream.vz", "Center coordinate for the maxwellian distribution. Used for calculating the backstream moments.", 0.0);
     Readparameters::add(pop + "_backstream.radius", "Radius of the maxwellian distribution. Used for calculating the backstream moments. If set to 0 (default), the backstream/non-backstream DROs are skipped.", 0.0);
  }

  return true;
}


bool ObjectWrapper::getParameters() {
   typedef Readparameters RP;
   
   // Particle population parameters
   unsigned int popCounter=0;

   for(unsigned int i =0; i < getObjectWrapper().particleSpecies.size(); i++) {

      species::Species& species=getObjectWrapper().particleSpecies[i];
      vmesh::MeshParameters& vMesh=getObjectWrapper().velocityMeshes[i];

      const std::string& pop = species.name;

      // Sanity check name
      if(species.name != vMesh.name) {
         std::cerr << "ParticlePopulation parse error: Name " << species.name << " != " << vMesh.name << std::endl;
         return false;
      }

      // Elementary particle parameters
      RP::get(pop + "_properties.charge", species.charge);
      species.charge *= physicalconstants::CHARGE;

      RP::get(pop + "_properties.mass", species.mass);
      std::string massUnit;
      RP::get(pop + "_properties.mass_units", massUnit);
      bool testPop;
      RP::get(pop + "_properties.test_population", testPop);
      species.isTestSpecies = testPop;
      bool propagatePop;
      RP::get(pop + "_properties.propagate_population", propagatePop);
      species.propagateSpecies = propagatePop;
      if(massUnit == "PROTON") {
         species.mass *= physicalconstants::MASS_PROTON;
      } else if(massUnit == "ELECTRON") {
         species.mass *= physicalconstants::MASS_ELECTRON;
      } else {
         std::cerr << "Invalid mass unit for species " << pop << ": '" << massUnit << "'" << std::endl;
         return false;
      }

      // sparsity parameters
      RP::get(pop + "_sparse.minValue", species.sparseMinValue);
      RP::get(pop + "_sparse.blockAddWidthV", species.sparseBlockAddWidthV);
      RP::get(pop + "_sparse.conserve_mass", species.sparse_conserve_mass);
      RP::get(pop + "_sparse.dynamicAlgorithm", species.sparseDynamicAlgorithm);
      RP::get(pop + "_sparse.dynamicBulkValue1", species.sparseDynamicBulkValue1);
      RP::get(pop + "_sparse.dynamicBulkValue2", species.sparseDynamicBulkValue2);
      RP::get(pop + "_sparse.dynamicMinValue1", species.sparseDynamicMinValue1);
      RP::get(pop + "_sparse.dynamicMinValue2", species.sparseDynamicMinValue2);


      // Particle velocity space properties
      RP::get(pop + "_vspace.vx_min",vMesh.meshLimits[0]);
      RP::get(pop + "_vspace.vx_max",vMesh.meshLimits[1]);
      RP::get(pop + "_vspace.vy_min",vMesh.meshLimits[2]);
      RP::get(pop + "_vspace.vy_max",vMesh.meshLimits[3]);
      RP::get(pop + "_vspace.vz_min",vMesh.meshLimits[4]);
      RP::get(pop + "_vspace.vz_max",vMesh.meshLimits[5]);
      RP::get(pop + "_vspace.vx_length",vMesh.gridLength[0]);
      RP::get(pop + "_vspace.vy_length",vMesh.gridLength[1]);
      RP::get(pop + "_vspace.vz_length",vMesh.gridLength[2]);
      if(vMesh.gridLength[0] > MAX_BLOCKS_PER_DIM  ||
            vMesh.gridLength[1] > MAX_BLOCKS_PER_DIM  ||
            vMesh.gridLength[2] > MAX_BLOCKS_PER_DIM ) {

         // Build error message as a string first, so that the cerr output hapens atomically and we don't spam
         // thousands of unreadable lines through each other
         std::string errormsg = "(VSPACE) ERROR: Velocity mesh for population " + species.name + " has too many blocks per dimension. Maximum defined in MAX_BLOCKS_PER_DIM is " + std::to_string(MAX_BLOCKS_PER_DIM) + " "
            + std::string(__FILE__) + ":" + std::to_string(__LINE__) + "\n";
         std::cerr << errormsg;
      }

      vMesh.blockLength[0] = vMesh.blockLength[1] = vMesh.blockLength[2] = WID;
      int maxRefLevel; // Temporary variable, since target value is a uint8_t
      RP::get(pop + "_vspace.max_refinement_level",maxRefLevel);
      vMesh.refLevelMaxAllowed = maxRefLevel;

      
      //Get backstream/non-backstream moments parameters
      Readparameters::get(pop + "_backstream.radius", species.backstreamRadius);
      Readparameters::get(pop + "_backstream.vx", species.backstreamV[0]);
      Readparameters::get(pop + "_backstream.vy", species.backstreamV[1]);
      Readparameters::get(pop + "_backstream.vz", species.backstreamV[2]);
   }

   return true;
}
