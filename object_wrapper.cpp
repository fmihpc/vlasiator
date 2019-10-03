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
     
     // Thermal / suprathermal parameters
     Readparameters::add(pop + "_thermal.vx", "Center coordinate for the maxwellian distribution. Used for calculating the suprathermal moments.", -500000.0);
     Readparameters::add(pop + "_thermal.vy", "Center coordinate for the maxwellian distribution. Used for calculating the suprathermal moments.", 0.0);
     Readparameters::add(pop + "_thermal.vz", "Center coordinate for the maxwellian distribution. Used for calculating the suprathermal moments.", 0.0);
     Readparameters::add(pop + "_thermal.radius", "Radius of the maxwellian distribution. Used for calculating the suprathermal moments. If set to 0 (default), the thermal/suprathermal DROs are skipped.", 0.0);

     // Precipitation parameters
     Readparameters::add(pop + "_precipitation.nChannels", "Number of energy channels for precipitation differential flux evaluation", 16);
     Readparameters::add(pop + "_precipitation.emin", "Lowest energy channel (in eV) for precipitation differential flux evaluation", 0.1);
     Readparameters::add(pop + "_precipitation.emax", "Highest energy channel (in eV) for precipitation differential flux evaluation", 100.0);
     Readparameters::add(pop + "_precipitation.lossConeAngle", "Fixed loss cone opening angle (in deg) for precipitation differential flux evaluation", 10.0);

     // Energy density parameters
     Readparameters::add(pop + "_energydensity.limit1", "Lower limit of second bin for energy density, given in units of solar wind ram energy.", 5.0);
     Readparameters::add(pop + "_energydensity.limit2", "Lower limit of third bin for energy density, given in units of solar wind ram energy.", 10.0);
     Readparameters::add(pop + "_energydensity.solarwindspeed", "Incoming solar wind velocity magnitude in m/s. Used for calculating energy densities.", 0.0);
     Readparameters::add(pop + "_energydensity.solarwindenergy", "Incoming solar wind ram energy in eV. Used for calculating energy densities.", 0.0);
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

      
      //Get thermal / suprathermal moments parameters
      Readparameters::get(pop + "_thermal.radius", species.thermalRadius);
      Readparameters::get(pop + "_thermal.vx", species.thermalV[0]);
      Readparameters::get(pop + "_thermal.vy", species.thermalV[1]);
      Readparameters::get(pop + "_thermal.vz", species.thermalV[2]);

      //Get energy density parameters
      Readparameters::get(pop + "_energydensity.limit1", species.EnergyDensityLimit1);
      Readparameters::get(pop + "_energydensity.limit2", species.EnergyDensityLimit2);
      Readparameters::get(pop + "_energydensity.solarwindenergy", species.SolarWindEnergy);
      Readparameters::get(pop + "_energydensity.solarwindspeed", species.SolarWindSpeed);
      
      const Real EPSILON = 1.e-25;
      if (species.SolarWindEnergy < EPSILON) {
	 // Energy stored internally in SI units
	 species.SolarWindEnergy = 0.5 * species.mass * species.SolarWindSpeed * species.SolarWindSpeed;
      } else {
	 species.SolarWindEnergy = species.SolarWindEnergy*physicalconstants::CHARGE;
      }

      // Get precipitation parameters
      Readparameters::get(pop + "_precipitation.nChannels", species.precipitationNChannels);
      Readparameters::get(pop + "_precipitation.emin", species.precipitationEmin);
      Readparameters::get(pop + "_precipitation.emax", species.precipitationEmax);
      Readparameters::get(pop + "_precipitation.lossConeAngle", species.precipitationLossConeAngle);
      // Convert from eV to SI units
      species.precipitationEmin = species.precipitationEmin*physicalconstants::CHARGE;
      species.precipitationEmax = species.precipitationEmax*physicalconstants::CHARGE;
   }

   return true;
}
