/*
 * This file is part of Vlasiator.
 * 
 * Copyright 2011-2015 Finnish Meteorological Institute
 * 
 */

#include "project.h"
#include <cstdlib>
#include "../common.h"
#include "../parameters.h"
#include "../readparameters.h"
#include "../vlasovmover.h"
#include "../particle_species.h"
#include "../logger.h"
#include "../object_wrapper.h"

#include "Alfven/Alfven.h"
#include "Diffusion/Diffusion.h"
#include "Dispersion/Dispersion.h"
#include "Distributions/Distributions.h"
#include "ElectricSail/electric_sail.h"
#include "Firehose/Firehose.h"
#include "Flowthrough/Flowthrough.h"
#include "Fluctuations/Fluctuations.h"
#include "Harris/Harris.h"
#include "KHB/KHB.h"
#include "Larmor/Larmor.h"
#include "Magnetosphere/Magnetosphere.h"
#include "MultiPeak/MultiPeak.h"
#include "VelocityBox/VelocityBox.h"
#include "Riemann1/Riemann1.h"
#include "Shock/Shock.h"
#include "Template/Template.h"
#include "test_fp/test_fp.h"
#include "testHall/testHall.h"
#include "test_trans/test_trans.h"
#include "verificationLarmor/verificationLarmor.h"
#include "../backgroundfield/backgroundfield.h"
#include "../backgroundfield/constantfield.hpp"
#include "Shocktest/Shocktest.h"
#include "Poisson/poisson_test.h"

using namespace std;

extern Logger logFile;

char projects::Project::rngStateBuffer[256];
random_data projects::Project::rngDataBuffer;

/** Struct for creating a new velocity mesh.
 * The values are read from the configuration file and 
 * copied to ObjectWrapper::velocityMeshes.*/
struct VelocityMeshParams {
   vector<string> name;
   vector<double> vx_min;
   vector<double> vy_min;
   vector<double> vz_min;
   vector<double> vx_max;
   vector<double> vy_max;
   vector<double> vz_max;
   vector<double> vx_length;
   vector<double> vy_length;
   vector<double> vz_length;
   vector<unsigned int> maxRefLevels;
   
   void resize(const size_t& size) {
      name.resize(1);
      vx_min.resize(1);
      vy_min.resize(1);
      vz_min.resize(1);
      vx_max.resize(1);
      vy_max.resize(1);
      vz_max.resize(1);
      vx_length.resize(1);
      vy_length.resize(1);
      vz_length.resize(1);
      maxRefLevels.resize(1);
   }
};

static VelocityMeshParams* velMeshParams = NULL;

namespace projects {
   Project::Project() { 
      baseClassInitialized = false;
   }
   
   Project::~Project() { }
   
   void Project::addParameters() {
      typedef Readparameters RP;
      // TODO add all projects' static addParameters() functions here.
      projects::Alfven::addParameters();
      projects::Diffusion::addParameters();
      projects::Dispersion::addParameters();
      projects::Distributions::addParameters();
      projects::ElectricSail::addParameters();
      projects::Firehose::addParameters();
      projects::Flowthrough::addParameters();
      projects::Fluctuations::addParameters();
      projects::Harris::addParameters();
      projects::KHB::addParameters();
      projects::Larmor::addParameters();
      projects::Magnetosphere::addParameters();
      projects::MultiPeak::addParameters();
      projects::VelocityBox::addParameters();
      projects::Riemann1::addParameters();
      projects::Shock::addParameters();
      projects::Template::addParameters();
      projects::test_fp::addParameters();
      projects::TestHall::addParameters();
      projects::test_trans::addParameters();
      projects::verificationLarmor::addParameters();
      projects::Shocktest::addParameters();
      projects::PoissonTest::addParameters();
      RP::add("Project_common.seed", "Seed for the RNG", 42);
      
      // Add parameters needed to create particle populations
      RP::addComposing("ParticlePopulation.name","Name of the simulated particle population (string)");
      RP::addComposing("ParticlePopulation.charge","Particle charge, in units of elementary charges (int)");
      RP::addComposing("ParticlePopulation.mass_units","Units in which particle mass is given, either 'PROTON' or 'ELECTRON' (string)");
      RP::addComposing("ParticlePopulation.mass","Particle mass in given units (float)");
      RP::addComposing("ParticlePopulation.sparse_min_value","Minimum value of distribution function in any cell of a velocity block for the block to be considered to have content");
      RP::addComposing("ParticlePopulation.mesh","Name of the velocity mesh the species should use (string)");
      
      // Add parameters needed to create velocity meshes
      RP::addComposing("velocitymesh.name","Name of the mesh (unique,string)");
      RP::addComposing("velocitymesh.vx_min","Minimum value for velocity mesh vx-coordinates.");
      RP::addComposing("velocitymesh.vx_max","Maximum value for velocity mesh vx-coordinates.");
      RP::addComposing("velocitymesh.vy_min","Minimum value for velocity mesh vy-coordinates.");
      RP::addComposing("velocitymesh.vy_max","Maximum value for velocity mesh vx-coordinates.");
      RP::addComposing("velocitymesh.vz_min","Minimum value for velocity mesh vz-coordinates.");
      RP::addComposing("velocitymesh.vz_max","Maximum value for velocity mesh vx-coordinates.");
      RP::addComposing("velocitymesh.vx_length","Initial number of velocity blocks in vx-direction.");
      RP::addComposing("velocitymesh.vy_length","Initial number of velocity blocks in vy-direction.");
      RP::addComposing("velocitymesh.vz_length","Initial number of velocity blocks in vz-direction.");
      RP::addComposing("velocitymesh.max_refinement_level","Maximum allowed mesh refinement level.");

      // These parameters are only read if the 'velocitymesh.' parameters are not defined 
      // in order to support older configuration files.
      Real defValue = numeric_limits<Real>::infinity();
      Readparameters::add("gridbuilder.vx_min","Minimum value for velocity mesh vx-coordinates.",defValue);
      Readparameters::add("gridbuilder.vx_max","Maximum value for velocity mesh vx-coordinates.",defValue);
      Readparameters::add("gridbuilder.vy_min","Minimum value for velocity mesh vy-coordinates.",defValue);
      Readparameters::add("gridbuilder.vy_max","Maximum value for velocity mesh vy-coordinates.",defValue);
      Readparameters::add("gridbuilder.vz_min","Minimum value for velocity mesh vz-coordinates.",defValue);
      Readparameters::add("gridbuilder.vz_max","Maximum value for velocity mesh vz-coordinates.",defValue);
      Readparameters::add("gridbuilder.vx_length","Initial number of velocity blocks in vx-direction.",(vmesh::LocalID)0);
      Readparameters::add("gridbuilder.vy_length","Initial number of velocity blocks in vy-direction.",(vmesh::LocalID)0);
      Readparameters::add("gridbuilder.vz_length","Initial number of velocity blocks in vz-direction.",(vmesh::LocalID)0);
   }

   void Project::getParameters() {
      typedef Readparameters RP;
      RP::get("Project_common.seed", this->seed);
      RP::get("ParticlePopulation.name",popNames);
      RP::get("ParticlePopulation.charge",popCharges);
      RP::get("ParticlePopulation.mass_units",popMassUnits);
      RP::get("ParticlePopulation.mass",popMasses);
      RP::get("ParticlePopulation.sparse_min_value",popSparseMinValue);
      RP::get("ParticlePopulation.mesh",popMeshNames);

      if (velMeshParams == NULL) velMeshParams = new VelocityMeshParams();
      RP::get("velocitymesh.name",velMeshParams->name);
      RP::get("velocitymesh.vx_min",velMeshParams->vx_min);
      RP::get("velocitymesh.vy_min",velMeshParams->vy_min);
      RP::get("velocitymesh.vz_min",velMeshParams->vz_min);
      RP::get("velocitymesh.vx_max",velMeshParams->vx_max);
      RP::get("velocitymesh.vy_max",velMeshParams->vy_max);
      RP::get("velocitymesh.vz_max",velMeshParams->vz_max);
      RP::get("velocitymesh.vx_length",velMeshParams->vx_length);
      RP::get("velocitymesh.vy_length",velMeshParams->vy_length);
      RP::get("velocitymesh.vz_length",velMeshParams->vz_length);
      RP::get("velocitymesh.max_refinement_level",velMeshParams->maxRefLevels);
   }

   /** Initialize the Project. Velocity mesh and particle population 
    * parameters are read from the configuration file, and corresponding internal 
    * variables are created here.
    * NOTE: Each project must call this function!
    * @return If true, particle species and velocity meshes were created successfully.*/
   bool Project::initialize() {
      typedef Readparameters RP;
      
      // Basic error checking
      bool success = true;
      if (popNames.size() != popCharges.size()) success = false;
      if (popNames.size() != popMassUnits.size()) success = false;
      if (popNames.size() != popMasses.size()) success = false;
      if (popNames.size() != popSparseMinValue.size()) success = false;
      if (popNames.size() != popMeshNames.size()) success = false;
      if (success == false) {
         stringstream ss;
         ss << "(PROJECT) ERROR in configuration file particle population definitions at ";
         ss << __FILE__ << ":" << __LINE__ << endl;
         ss << "\t vector sizes are: " << popNames.size() << ' ' << popMassUnits.size();
         ss << ' ' << popMasses.size() << ' ' << popSparseMinValue.size() << ' ';
         ss << popMeshNames.size() << endl;
         cerr << ss.str(); return success;
      }

      if (velMeshParams->vx_min.size() != velMeshParams->name.size()) success = false;
      if (velMeshParams->vy_min.size() != velMeshParams->name.size()) success = false;
      if (velMeshParams->vz_min.size() != velMeshParams->name.size()) success = false;
      if (velMeshParams->vx_max.size() != velMeshParams->name.size()) success = false;
      if (velMeshParams->vy_max.size() != velMeshParams->name.size()) success = false;
      if (velMeshParams->vz_max.size() != velMeshParams->name.size()) success = false;
      if (velMeshParams->vx_length.size() != velMeshParams->name.size()) success = false;
      if (velMeshParams->vy_length.size() != velMeshParams->name.size()) success = false;
      if (velMeshParams->vz_length.size() != velMeshParams->name.size()) success = false;
      if (velMeshParams->maxRefLevels.size() != velMeshParams->name.size()) success = false;
      if (success == false) {
         stringstream ss;
         ss << "(PROJECT) ERROR in configuration file velocity mesh definitions at ";
         ss << __FILE__ << ":" << __LINE__ << endl;
         cerr << ss.str(); return success;
      }

      ObjectWrapper& owrapper = getObjectWrapper();
      
      // ********** VELOCITY MESHES  ********** //
      
      // If velocity meshes were not defined under 'velocitymesh' config file region, 
      // read the parameters from 'gridbuilder'
      if (velMeshParams->name.size() == 0) {
         velMeshParams->resize(1);
         velMeshParams->name[0] = "gridbuilder";
         RP::get("gridbuilder.vx_min",velMeshParams->vx_min[0]);
         RP::get("gridbuilder.vy_min",velMeshParams->vy_min[0]);
         RP::get("gridbuilder.vz_min",velMeshParams->vz_min[0]);
         RP::get("gridbuilder.vx_max",velMeshParams->vx_max[0]);
         RP::get("gridbuilder.vy_max",velMeshParams->vy_max[0]);
         RP::get("gridbuilder.vz_max",velMeshParams->vz_max[0]);
         RP::get("gridbuilder.vx_length",velMeshParams->vx_length[0]);
         RP::get("gridbuilder.vy_length",velMeshParams->vy_length[0]);
         RP::get("gridbuilder.vz_length",velMeshParams->vz_length[0]);
         velMeshParams->maxRefLevels[0] = 0;
      }

      // Store velocity mesh parameters
      for (size_t m=0; m<velMeshParams->name.size(); ++m) {         
         // Check that a mesh with the same name doesn't already exists:
         bool addMesh = true;         
         for (size_t i=0; i<owrapper.velocityMeshes.size(); ++i) {
            if (velMeshParams->name[m] == owrapper.velocityMeshes[i].name) {
               addMesh = false;
               break;
            }
         }
         if (addMesh == false) {
            stringstream ss;
            ss << "(PROJECT) ERROR: Velocity mesh called '" << velMeshParams->name[m] << "' already exists in ";
            ss << __FILE__ << ":" << __LINE__ << endl;
            cerr << ss.str(); success = false;
            continue;
         }

         vmesh::MeshParameters meshParams;
         meshParams.name = velMeshParams->name[m];
         meshParams.meshLimits[0] = velMeshParams->vx_min[m];
         meshParams.meshLimits[1] = velMeshParams->vx_max[m];
         meshParams.meshLimits[2] = velMeshParams->vy_min[m];
         meshParams.meshLimits[3] = velMeshParams->vy_max[m];
         meshParams.meshLimits[4] = velMeshParams->vz_min[m];
         meshParams.meshLimits[5] = velMeshParams->vz_max[m];
         meshParams.gridLength[0] = velMeshParams->vx_length[m];
         meshParams.gridLength[1] = velMeshParams->vy_length[m];
         meshParams.gridLength[2] = velMeshParams->vz_length[m];
         meshParams.blockLength[0] = WID;
         meshParams.blockLength[1] = WID;
         meshParams.blockLength[2] = WID;
         meshParams.refLevelMaxAllowed = velMeshParams->maxRefLevels[m];
         owrapper.velocityMeshes.push_back(meshParams);
      }

      delete velMeshParams; velMeshParams = NULL;
      
      // ********** PARTICLE SPECIES ********** //

      // If particle population(s) have not been defined, add protons as a default population
      // and assume that the velocity mesh is called 'gridbuilder'
      if (popNames.size() == 0) {
         species::Species population;
         population.name   = "proton";
         population.charge = physicalconstants::CHARGE;
         population.mass   = physicalconstants::MASS_PROTON;
         population.sparseMinValue = Parameters::sparseMinValue;
         
         size_t index=owrapper.velocityMeshes.size();
         for (size_t m=0; m<owrapper.velocityMeshes.size(); ++m) {
            if (owrapper.velocityMeshes[m].name == "gridbuilder") {
               index = m; break;
            }
         }
         if (index >= owrapper.velocityMeshes.size()) {
            stringstream ss;
            ss << "(PROJECT) ERROR: Could not associate default particle population with a velocity ";
            ss << "mesh in " << __FILE__ << ":" << __LINE__ << endl;
            cerr << ss.str(); success = false;
            return success;
         }
         population.velocityMesh = index;

         owrapper.particleSpecies.push_back(population);
         printPopulations();
         baseClassInitialized = success;
         return success;
      }

      // Parse populations from configuration file parameters:
      for (size_t p=0; p<popNames.size(); ++p) {
         species::Species population;
         population.name = popNames[p];
         population.charge = popCharges[p]*physicalconstants::CHARGE;
         double massUnits = 0;
         if (popMassUnits[p] == "PROTON") massUnits = physicalconstants::MASS_PROTON;
         else if (popMassUnits[p] == "ELECTRON") massUnits = physicalconstants::MASS_ELECTRON;
         else {
            stringstream ss;
            ss << "(PROJECT) ERROR: Could not determine species '" << popNames[p] << "' mass units in ";
            ss << __FILE__ << ":" << __LINE__ << endl;
            cerr << ss.str(); success = false;
         }
         population.mass = massUnits*popMasses[p];
         population.sparseMinValue = popSparseMinValue[p];
         
         bool meshFound = false;
         for (size_t m=0; m<owrapper.velocityMeshes.size(); ++m) {
            if (owrapper.velocityMeshes[m].name == popMeshNames[p]) {
               population.velocityMesh = m;
               meshFound = true; break;
            }
         }
         if (meshFound == false) {
            stringstream ss;
            ss << "(PROJECT) ERROR: Could not associate population '" << popNames[p] << "' with a velocity mesh in ";
            ss << __FILE__ << ":" << __LINE__ << endl; 
            cerr << ss.str(); success = false;
         }

         if (success == false) {
            stringstream ss;
            ss << "ERROR in population '" << popNames[p] << "' parameters" << endl;
            cerr << ss.str(); continue;
         }
         
         owrapper.particleSpecies.push_back(population);
      }

      if (success == false) {
         logFile << "ERROR in configuration file particle population definitions!" << endl << writeVerbose;
      } else {
         printPopulations();
      }

      baseClassInitialized = success;
      return success;
   }
   
   /** Check if base class has been initialized.
    * @return If true, base class was successfully initialized.*/
   bool Project::initialized() {return baseClassInitialized;}

   /*! Base class sets zero background field */
   void Project::setCellBackgroundField(SpatialCell* cell) const {
      ConstantField bgField;
      bgField.initialize(0,0,0); //bg bx, by,bz
      setBackgroundField(bgField,cell->parameters, cell->derivatives,cell->derivativesBVOL);
      
      // Print a warning message to stderr so that the user knows to check 
      // that it is still correct to call the base class setCellBackgroundField function.
      static int printed = false;
      if (printed == false) {
         int rank;
         MPI_Comm_rank(MPI_COMM_WORLD,&rank);
         if (rank == 0) {
            cerr << "(Project.cpp) WARNING: Base class 'setCellBackgroundField' in " << __FILE__ << ":" << __LINE__;
            cerr << " called, make sure this is intentional" << endl;
            printed = true;
         }
      }
   }

   void Project::setCell(SpatialCell* cell) {
      // Set up cell parameters:
      calcCellParameters(cell,0.0);

      cell->parameters[CellParams::RHOLOSSADJUST] = 0.0;
      cell->parameters[CellParams::RHOLOSSVELBOUNDARY] = 0.0;

      for (size_t p=0; p<getObjectWrapper().particleSpecies.size(); ++p) {
         this->setVelocitySpace(p,cell);
      }

      //let's get rid of blocks not fulfilling the criteria here to save memory.
      //cell->adjustSingleCellVelocityBlocks();

      // Passing true for the doNotSkip argument as we want to calculate 
      // the moment no matter what when this function is called.
      calculateCellMoments(cell,true,true);
   }

   std::vector<vmesh::GlobalID> Project::findBlocksToInitialize(spatial_cell::SpatialCell* cell,const int& popID) const {
      vector<vmesh::GlobalID> blocksToInitialize;
      const uint8_t refLevel = 0;

      const vmesh::LocalID* vblocks_ini = cell->get_velocity_grid_length(popID,refLevel);
      
      for (uint kv=0; kv<vblocks_ini[2]; ++kv) 
         for (uint jv=0; jv<vblocks_ini[1]; ++jv)
            for (uint iv=0; iv<vblocks_ini[0]; ++iv) {
               vmesh::LocalID blockIndices[3];
               blockIndices[0] = iv;
               blockIndices[1] = jv;
               blockIndices[2] = kv;
               const vmesh::GlobalID blockGID = cell->get_velocity_block(popID,blockIndices,refLevel);

               cell->add_velocity_block(blockGID,popID);
               blocksToInitialize.push_back(blockGID);
      }

      return blocksToInitialize;
   }
   
   /** Write simulated particle populations to logfile.*/
   void Project::printPopulations() {
      logFile << "(PROJECT): Loaded particle populations are:" << endl;
      
      for (size_t p=0; p<getObjectWrapper().particleSpecies.size(); ++p) {
         const species::Species& spec = getObjectWrapper().particleSpecies[p];
         logFile << "Population #" << p << endl;
         logFile << "\t name             : '" << spec.name << "'" << endl;
         logFile << "\t charge           : '" << spec.charge << "'" << endl;
         logFile << "\t mass             : '" << spec.mass << "'" << endl;
         logFile << "\t sparse threshold : '" << spec.sparseMinValue << "'" << endl;
         logFile << "\t velocity mesh    : '" << getObjectWrapper().velocityMeshes[spec.velocityMesh].name << "'" << endl;
         logFile << endl;
      }
      logFile << write;
   }
   
   /** Calculate the volume averages of distribution function for the 
    * given particle population in the given spatial cell. The velocity block 
    * is defined by its local ID. The function returns the maximum value of the 
    * distribution function within the velocity block. If it is below the sparse 
    * min value for the population, this block should be removed or marked 
    * as a no-content block.
    * @param cell Spatial cell.
    * @param blockLID Velocity block local ID within the spatial cell.
    * @param popID Population ID.
    * @return Maximum value of the calculated distribution function.*/
   Real Project::setVelocityBlock(spatial_cell::SpatialCell* cell,const vmesh::LocalID& blockLID,const int& popID) const {
      // If simulation doesn't use one or more velocity coordinates, 
      // only calculate the distribution function for one layer of cells.
      uint WID_VX = WID;
      uint WID_VY = WID;
      uint WID_VZ = WID;
      switch (Parameters::geometry) {         
         case geometry::XY4D:
            WID_VZ=1;
            break;
         case geometry::XZ4D:
            WID_VY=1;
            break;
         default:
            break;
      }

      // Fetch spatial cell coordinates and size
      creal x  = cell->parameters[CellParams::XCRD];
      creal y  = cell->parameters[CellParams::YCRD];
      creal z  = cell->parameters[CellParams::ZCRD];
      creal dx = cell->parameters[CellParams::DX];
      creal dy = cell->parameters[CellParams::DY];
      creal dz = cell->parameters[CellParams::DZ];

      const Real* parameters = cell->get_block_parameters(popID);
      Realf* data = cell->get_data(popID);
      
      creal vxBlock = parameters[blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VXCRD];
      creal vyBlock = parameters[blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VYCRD];
      creal vzBlock = parameters[blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VZCRD];
      creal dvxCell = parameters[blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVX];
      creal dvyCell = parameters[blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVY];
      creal dvzCell = parameters[blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVZ];
      
      // Calculate volume average of distribution function for each phase-space cell in the block.
      Real maxValue = 0.0;
      for (uint kc=0; kc<WID_VZ; ++kc) for (uint jc=0; jc<WID_VY; ++jc) for (uint ic=0; ic<WID_VX; ++ic) {
         creal vxCell = vxBlock + ic*dvxCell;
         creal vyCell = vyBlock + jc*dvyCell;
         creal vzCell = vzBlock + kc*dvzCell;
         creal average =
            calcPhaseSpaceDensity(
               x, y, z, dx, dy, dz,
               vxCell,vyCell,vzCell,
               dvxCell,dvyCell,dvzCell,popID);
         if (average != 0.0) {
            data[blockLID*SIZE_VELBLOCK+cellIndex(ic,jc,kc)] = average;
            maxValue = max(maxValue,average);
         }
      }
      
      return maxValue;
   }
   
   void Project::setVelocitySpace(const int& popID,SpatialCell* cell) const {
      vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& vmesh = cell->get_velocity_mesh(popID);

      vector<vmesh::GlobalID> blocksToInitialize = this->findBlocksToInitialize(cell,popID);
      vector<vmesh::GlobalID> removeList;
      for (uint i=0; i<blocksToInitialize.size(); ++i) {
         const vmesh::GlobalID blockGID = blocksToInitialize[i];
         const vmesh::LocalID blockLID = vmesh.getLocalID(blockGID);
         if (blockLID == vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>::invalidLocalID()) {
            cerr << "ERROR, invalid local ID in " << __FILE__ << ":" << __LINE__ << endl;
            exit(1);
         }

         const Real maxValue = setVelocityBlock(cell,blockLID,popID);
         if (maxValue < getObjectWrapper().particleSpecies[popID].sparseMinValue) removeList.push_back(blockGID);
      }

      // Get AMR refinement criterion and use it to test which blocks should be refined
      amr_ref_criteria::Base* refCriterion = getObjectWrapper().amrVelRefCriteria.create(Parameters::amrVelRefCriterion);
      if (refCriterion == NULL) {
         if (rescalesDensity(popID) == true) rescaleDensity(cell,popID);
         return;
      }
      refCriterion->initialize("");

      // Remove blocks with f below sparse min value
      for (size_t b=0; b<removeList.size(); ++b) cell->remove_velocity_block(removeList[b],popID);

      // Loop over blocks in the spatial cell until we reach the maximum
      // refinement level, or until there are no more blocks left to refine
      bool refine = true;
      uint currentLevel = 0;
      if (currentLevel == Parameters::amrMaxVelocityRefLevel) refine = false;
      while (refine == true) {
         removeList.clear();
         
         // Loop over blocks and add blocks to be refined to vector refineList
         vector<vmesh::GlobalID> refineList;
         const vmesh::LocalID startIndex = 0;
         const vmesh::LocalID endIndex   = cell->get_number_of_velocity_blocks(popID);
         for (vmesh::LocalID blockLID=startIndex; blockLID<endIndex; ++blockLID) {
            vector<vmesh::GlobalID> nbrs;
            int32_t refLevelDifference;
            const vmesh::GlobalID blockGID = vmesh.getGlobalID(blockLID);

            // Fetch block data and nearest neighbors
            Realf array[(WID+2)*(WID+2)*(WID+2)];
            cell->fetch_data<1>(blockGID,vmesh,cell->get_data(0,popID),array);

            // If block should be refined, add it to refine list
            if (refCriterion->evaluate(array,popID) > Parameters::amrRefineLimit) {
               refineList.push_back(blockGID);
            }
         }

         // Refine blocks in vector refineList. All blocks that were created 
         // as a result of the refine, including blocks created because of induced 
         // refinement, are added to map insertedBlocks
         map<vmesh::GlobalID,vmesh::LocalID> insertedBlocks;
         for (size_t b=0; b<refineList.size(); ++b) {
            cell->refine_block(refineList[b],insertedBlocks,popID);
         }

         // Loop over blocks in map insertedBlocks and recalculate 
         // values of distribution functions
         for (map<vmesh::GlobalID,vmesh::LocalID>::const_iterator it=insertedBlocks.begin(); it!=insertedBlocks.end(); ++it) {
            const vmesh::GlobalID blockGID = it->first;
            const vmesh::LocalID blockLID = it->second;
            const Real maxValue = setVelocityBlock(cell,blockLID,popID);
            if (maxValue <= getObjectWrapper().particleSpecies[popID].sparseMinValue) 
              removeList.push_back(it->first);
         }

         // Remove blocks with f below sparse min value
         for (size_t b=0; b<removeList.size(); ++b) cell->remove_velocity_block(removeList[b],popID);

         if (refineList.size() == 0) refine = false;
         ++currentLevel;
         if (currentLevel == Parameters::amrMaxVelocityRefLevel) refine = false;
      }

      delete refCriterion;

      if (rescalesDensity(popID) == true) rescaleDensity(cell,popID);
   }

   /** Check if the project wants to rescale densities.
    * @param popID ID of the particle species.
    * @return If true, rescaleDensity is called for this species.*/
   bool Project::rescalesDensity(const int& popID) const {
      return false;
   }

   /** Rescale the distribution function of the given particle species so that 
    * the number density corresponds to the value returned by getCorrectNumberDensity.
    * @param cell Spatial cell.
    * @param popID ID of the particle species.*/
   void Project::rescaleDensity(spatial_cell::SpatialCell* cell,const int& popID) const {      
      // Re-scale densities
      Real sum = 0.0;
      Realf* data = cell->get_data(popID);
      const Real* blockParams = cell->get_block_parameters(popID);
      for (vmesh::LocalID blockLID=0; blockLID<cell->get_number_of_velocity_blocks(popID); ++blockLID) {
         Real tmp = 0.0;
         for (int i=0; i<WID3; ++i) tmp += data[blockLID*WID3+i];
         const Real DV3 = blockParams[BlockParams::DVX]*blockParams[BlockParams::DVY]*blockParams[BlockParams::DVZ];
         sum += tmp*DV3;
         blockParams += BlockParams::N_VELOCITY_BLOCK_PARAMS;
      }
      
      const Real correctSum = getCorrectNumberDensity(cell,popID);
      const Real ratio = correctSum / sum;
      
      for (size_t i=0; i<cell->get_number_of_velocity_blocks(popID)*WID3; ++i) {
         data[i] *= ratio;
      }
   }
   
   /*default one does not compute any parameters*/
   void Project::calcCellParameters(spatial_cell::SpatialCell* cell,creal& t) {
      static int printed = false;
      if (printed == false) {
         int rank;
         MPI_Comm_rank(MPI_COMM_WORLD,&rank);
         if (rank == 0) {
            cerr << "(Project.cpp) WARNING: Base class 'calcCellParameters' in " << __FILE__ << ":" << __LINE__;
            cerr << " called, make sure this is intentional" << endl;
            printed = true;
         }
      }
   }

   Real Project::calcPhaseSpaceDensity(
      creal& x, creal& y, creal& z,
      creal& dx, creal& dy, creal& dz,
      creal& vx, creal& vy, creal& vz,
      creal& dvx, creal& dvy, creal& dvz,
      const int& popID) const {
      cerr << "ERROR: Project::calcPhaseSpaceDensity called instead of derived class function!" << endl;
      exit(1);
      return -1.0;
   }
   /*!
     Get random number between 0 and 1.0. One should always first initialize the rng.
   */
   
   Real Project::getCorrectNumberDensity(spatial_cell::SpatialCell* cell,const int& popID) const {
      cerr << "ERROR: Project::getCorrectNumberDensity called instead of derived class function!" << endl;
      exit(1);
      return 0.0;
   }

   /** Get random number between 0 and 1.0. One should always first initialize the rng.
    * @param cell Spatial cell.
    * @return Uniformly distributed random number between 0 and 1.*/
   Real Project::getRandomNumber(spatial_cell::SpatialCell* cell) const {
#ifdef _AIX
      int64_t rndInt;
      random_r(&rndInt, &rngDataBuffer);
#else
      int32_t rndInt;
      random_r(&rngDataBuffer, &rndInt);
#endif
      Real rnd = (Real) rndInt / RAND_MAX;
      return rnd;
   }

   /*!  Set random seed (thread-safe). Seed is based on the seed read
     in from cfg + the seedModifier parameter

     \param seedModifier d. Seed is based on the seed read in from cfg + the seedModifier parameter
   */

   void Project::setRandomSeed(spatial_cell::SpatialCell* cell,CellID seedModifier) const {
      memset(&(this->rngDataBuffer), 0, sizeof(this->rngDataBuffer));
#ifdef _AIX
      initstate_r(this->seed+seedModifier, &(this->rngStateBuffer[0]), 256, NULL, &(this->rngDataBuffer));
#else
      initstate_r(this->seed+seedModifier, &(this->rngStateBuffer[0]), 256, &(this->rngDataBuffer));
#endif
   }

   /*!
     Set random seed (thread-safe) that is always the same for
     this particular cellID. Can be used to make reproducible
     simulations that do not depend on number of processes or threads.

     \param  cellParams The cell parameters list in each spatial cell
   */
   void Project::setRandomCellSeed(spatial_cell::SpatialCell* cell,const Real* const cellParams) const {
      const creal x = cellParams[CellParams::XCRD];
      const creal y = cellParams[CellParams::YCRD];
      const creal z = cellParams[CellParams::ZCRD];
      const creal dx = cellParams[CellParams::DX];
      const creal dy = cellParams[CellParams::DY];
      const creal dz = cellParams[CellParams::DZ];
      
      const CellID cellID = (int) ((x - Parameters::xmin) / dx) +
         (int) ((y - Parameters::ymin) / dy) * Parameters::xcells_ini +
         (int) ((z - Parameters::zmin) / dz) * Parameters::xcells_ini * Parameters::ycells_ini;
      setRandomSeed(cell,cellID);
   }
   
Project* createProject() {
   Project* rvalue = NULL;
   if(Parameters::projectName == "") {
      cerr << "No project specified! Please set 'project' parameter!" << endl;
      abort();
   }
   if(Parameters::projectName == "Alfven") {
      rvalue = new projects::Alfven;
   }
   if(Parameters::projectName == "Diffusion") {
      rvalue = new projects::Diffusion;
   }
   if(Parameters::projectName == "Dispersion") {
      rvalue = new projects::Dispersion;
   }
   if(Parameters::projectName == "Distributions") {
      rvalue = new projects::Distributions;
   }
   if (Parameters::projectName == "ElectricSail") {
      return new projects::ElectricSail;
   }
   if(Parameters::projectName == "Firehose") {
      rvalue = new projects::Firehose;
   }
   if(Parameters::projectName == "Flowthrough") {
      rvalue = new projects::Flowthrough;
   }
   if(Parameters::projectName == "Fluctuations") {
      rvalue = new projects::Fluctuations;
   }
   if(Parameters::projectName == "Harris") {
      rvalue = new projects::Harris;
   }
   if(Parameters::projectName == "KHB") {
      rvalue = new projects::KHB;
   }
   if(Parameters::projectName == "Larmor") {
      rvalue = new projects::Larmor;
   }
   if(Parameters::projectName == "Magnetosphere") {
      rvalue = new projects::Magnetosphere;
   }
   if(Parameters::projectName == "MultiPeak") {
      rvalue = new projects::MultiPeak;
   } 
   if(Parameters::projectName == "VelocityBox") {
      rvalue = new projects::VelocityBox;
   } 
   if(Parameters::projectName == "Riemann1") {
      rvalue = new projects::Riemann1;
   }
   if(Parameters::projectName == "Shock") {
      rvalue = new projects::Shock;
   }
   if(Parameters::projectName == "Template") {
      rvalue = new projects::Template;
   }
   if(Parameters::projectName == "test_fp") {
      rvalue = new projects::test_fp;
   }
   if(Parameters::projectName == "testHall") {
      rvalue = new projects::TestHall;
   }
   if(Parameters::projectName == "test_trans") {
      rvalue = new projects::test_trans;
   }
   if(Parameters::projectName == "verificationLarmor") {
      rvalue = new projects::verificationLarmor;
   }
   if(Parameters::projectName == "Shocktest") {
      rvalue = new projects::Shocktest;
   }
   if (Parameters::projectName == "PoissonTest") {
      rvalue = new projects::PoissonTest;
   }
   if (rvalue == NULL) {
      cerr << "Unknown project name!" << endl;
      abort();
   }

   getObjectWrapper().project = rvalue;
   return rvalue;
}

} // namespace projects
