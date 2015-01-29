#include "project.h"
#include <cstdlib>
#include "../common.h"
#include "../parameters.h"
#include "../readparameters.h"
#include "../vlasovmover.h"
#include "../particle_species.h"
#include "../logger.h"

#include "Alfven/Alfven.h"
#include "Diffusion/Diffusion.h"
#include "Dispersion/Dispersion.h"
#include "Distributions/Distributions.h"
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

using namespace std;

extern Logger logFile;

char projects::Project::rngStateBuffer[256];
random_data projects::Project::rngDataBuffer;

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
      RP::add("Project_common.seed", "Seed for the RNG", 42);
      
      // Add parameters needed to create particle populations
      RP::addComposing("ParticlePopulation.name","Name of the simulated particle population (string)");
      RP::addComposing("ParticlePopulation.charge","Particle charge, in units of elementary charges (int)");
      RP::addComposing("ParticlePopulation.mass_units","Units in which particle mass is given, either 'PROTON' or 'ELECTRON' (string)");
      RP::addComposing("ParticlePopulation.mass","Particle mass in given units (float)");
      RP::addComposing("ParticlePopulation.sparse_min_value","Minimum value of distribution function in any cell of a velocity block for the block to be considered to have content");
   }
   
   void Project::getParameters() {
      typedef Readparameters RP;
      RP::get("Project_common.seed", this->seed);
      RP::get("ParticlePopulation.name",popNames);
      RP::get("ParticlePopulation.charge",popCharges);
      RP::get("ParticlePopulation.mass_units",popMassUnits);
      RP::get("ParticlePopulation.mass",popMasses);
      RP::get("ParticlePopulation.sparse_min_value",popSparseMinValue);
   }

   bool Project::initialize() {
      // Basic error checking
      bool success = true;
      if (popNames.size() != popCharges.size()) success = false;
      if (popNames.size() != popMassUnits.size()) success = false;
      if (popNames.size() != popMasses.size()) success = false;
      if (popNames.size() != popSparseMinValue.size()) success = false;
      if (success == false) {
         cerr << "ERROR in configuration file particle population definitions!" << endl;
         cerr << "\t vector sizes are: " << popNames.size() << ' ' << popMassUnits.size();
         cerr << ' ' << popMasses.size() << ' ' << popSparseMinValue.size() << endl;
         return success;
      }
      
      // If particle population(s) have not been defined, add protons as a default population
      ObjectWrapper& owrapper = getObjectWrapper();
      if (popNames.size() == 0) {
         species::Species population;
         population.name   = "proton";
         population.charge = physicalconstants::CHARGE;
         population.mass   = physicalconstants::MASS_PROTON;
         population.sparseMinValue = Parameters::sparseMinValue;
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
         else success = false;
         population.mass = massUnits*popMasses[p];
         population.sparseMinValue = popSparseMinValue[p];
         
         if (success == false) {
            cerr << "ERROR in population '" << popNames[p] << "' parameters" << endl;
         }
         
         owrapper.particleSpecies.push_back(population);
      }

      if (success == false) {
         cerr << "ERROR in configuration file particle population definitions!" << endl;
      } else {
         printPopulations();
      }

      baseClassInitialized = success;
      return success;
   }
   
   /** Check if base class has been initialized.
    * @return If true, base class was successfully initialized.*/
   bool Project::initialized() {return baseClassInitialized;}

   /** Set active particle population. Successive calls to member functions 
    * should then return values for the selected population.
    * @param popID Population ID.*/
   void Project::setActivePopulation(const int& popID) { }
   
   /*! Base class sets zero background field */
   void Project::setCellBackgroundField(SpatialCell* cell) {
      ConstantField bgField;
      bgField.initialize(0,0,0); //bg bx, by,bz
      setBackgroundField(bgField,cell->parameters, cell->derivatives,cell->derivativesBVOL);
   }
   
   void Project::setCell(SpatialCell* cell) {
      // Set up cell parameters:
      this->calcCellParameters(&((*cell).parameters[0]), 0.0);

      cell->parameters[CellParams::RHOLOSSADJUST] = 0.0;
      cell->parameters[CellParams::RHOLOSSVELBOUNDARY] = 0.0;

      for (size_t p=0; p<getObjectWrapper().particleSpecies.size(); ++p) {
         this->setVelocitySpace(p,cell);
         //cell->adjustSingleCellVelocityBlocks();
         //cell->printMeshSizes();
      }

      //let's get rid of blocks not fulfilling the criteria here to save memory.
      //cell->adjustSingleCellVelocityBlocks();

      // Passing true for the doNotSkip argument as we want to calculate the moment no matter what when this function is called.
      calculateCellVelocityMoments(cell, true);
   }

   vector<vmesh::GlobalID> Project::findBlocksToInitialize(SpatialCell* cell) {
      vector<vmesh::GlobalID> blocksToInitialize;

      for (uint kv=0; kv<P::vzblocks_ini; ++kv) 
         for (uint jv=0; jv<P::vyblocks_ini; ++jv)
            for (uint iv=0; iv<P::vxblocks_ini; ++iv) {
               creal vx = P::vxmin + (iv+0.5) * SpatialCell::get_velocity_grid_block_size()[0]; // vx-coordinate of the centre
               creal vy = P::vymin + (jv+0.5) * SpatialCell::get_velocity_grid_block_size()[1]; // vy-
               creal vz = P::vzmin + (kv+0.5) * SpatialCell::get_velocity_grid_block_size()[2];

               //FIXME, add_velocity_blocks should  not be needed as set_value handles it!!
               //FIXME,  We should get_velocity_block based on indices, not v
               cell->add_velocity_block(cell->get_velocity_block(vx, vy, vz));
               blocksToInitialize.push_back(cell->get_velocity_block(vx, vy, vz));
      }

      return blocksToInitialize;
   }
   
   /** Write simulated particle populations to logfile.*/
   void Project::printPopulations() {
      logFile << "(PROJECT): Loaded particle populations are:" << endl;
      
      for (size_t p=0; p<getObjectWrapper().particleSpecies.size(); ++p) {
         logFile << "Population #" << p << endl;
         logFile << "\t name             : '" << getObjectWrapper().particleSpecies[p].name << "'" << endl;
         logFile << "\t charge           : '" << getObjectWrapper().particleSpecies[p].charge << "'" << endl;
         logFile << "\t mass             : '" << getObjectWrapper().particleSpecies[p].mass << "'" << endl;
         logFile << "\t sparse threshold : " << getObjectWrapper().particleSpecies[p].sparseMinValue << "'" << endl;
         logFile << endl;
      }
      logFile << write;
   }
   
   void Project::setVelocitySpace(const int& popID,SpatialCell* cell) {
      setActivePopulation(popID);
      cell->setActivePopulation(popID);

      vector<vmesh::GlobalID> blocksToInitialize = this->findBlocksToInitialize(cell);
      Real* parameters = cell->get_block_parameters();

      creal x = cell->parameters[CellParams::XCRD];
      creal y = cell->parameters[CellParams::YCRD];
      creal z = cell->parameters[CellParams::ZCRD];
      creal dx = cell->parameters[CellParams::DX];
      creal dy = cell->parameters[CellParams::DY];
      creal dz = cell->parameters[CellParams::DZ];
      
      for (uint i=0; i<blocksToInitialize.size(); ++i) {
         const vmesh::GlobalID blockGID = blocksToInitialize.at(i);
         const vmesh::LocalID blockLID = cell->get_velocity_block_local_id(blockGID);
         creal vxBlock = parameters[blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VXCRD];
         creal vyBlock = parameters[blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VYCRD];
         creal vzBlock = parameters[blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VZCRD];
         creal dvxCell = parameters[blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVX];
         creal dvyCell = parameters[blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVY];
         creal dvzCell = parameters[blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVZ];

         // Calculate volume average of distrib. function for each cell in the block.
         for (uint kc=0; kc<WID; ++kc) 
            for (uint jc=0; jc<WID; ++jc) 
               for (uint ic=0; ic<WID; ++ic) {
                  //FIXME, block/cell index should be handled by spatial cell function (create if it does not exist)
                  creal vxCell = vxBlock + ic*dvxCell;
                  creal vyCell = vyBlock + jc*dvyCell;
                  creal vzCell = vzBlock + kc*dvzCell;
                  creal average =
                  this->calcPhaseSpaceDensity(
                     x, y, z, dx, dy, dz,
                     vxCell,vyCell,vzCell,
                     dvxCell,dvyCell,dvzCell);

                  if (average != 0.0){
                     //FIXME!!! set_value is slow as we again have to convert v -> index
                     // We should set_value to a specific block index (as we already have it!)
                     creal vxCellCenter = vxBlock + (ic+convert<Real>(0.5))*dvxCell;
                     creal vyCellCenter = vyBlock + (jc+convert<Real>(0.5))*dvyCell;
                     creal vzCellCenter = vzBlock + (kc+convert<Real>(0.5))*dvzCell;
                     cell->set_value(vxCellCenter,vyCellCenter,vzCellCenter,average);
                  }
               }
      }

      // Get AMR refinement criterion and use it to test which blocks should be refined
      amr_ref_criteria::Base* refCriterion = getObjectWrapper().amrVelRefCriteria.create(Parameters::amrVelRefCriterion);
      if (refCriterion == NULL) return;
      refCriterion->initialize("");

      // Loop over blocks in the spatial cell until we reach the maximum
      // refinement level, or until there are no more blocks left to refine
      bool refine = true;
      //refine = false;
      uint currentLevel = 0;
      if (currentLevel == Parameters::amrMaxVelocityRefLevel) refine = false;
      while (refine == true) {
         // Loop over blocks and add blocks to be refined to vector refineList
         vector<vmesh::GlobalID> refineList;
         const vmesh::LocalID startIndex = 0;
         const vmesh::LocalID endIndex   = cell->get_number_of_velocity_blocks();
         for (vmesh::LocalID blockLID=startIndex; blockLID<endIndex; ++blockLID) {
            vector<vmesh::GlobalID> nbrs;
            int32_t refLevelDifference;
            const vmesh::GlobalID blockGID = cell->get_velocity_block_global_id(blockLID);

            // Fetch block data and nearest neighbors
            Realf array[(WID+2)*(WID+2)*(WID+2)];
            cell->fetch_data<1>(blockGID,cell->get_velocity_mesh(popID),cell->get_data(),array);

            // If block should be refined, add it to refine list
            if (refCriterion->evaluate(array) > Parameters::amrRefineLimit) {
               refineList.push_back(blockGID);
            }
         }

         // Refine blocks in vector refineList. All blocks that were created 
         // as a result of the refine, including blocks created because of induced 
         // refinement, are added to map insertedBlocks
         map<vmesh::GlobalID,vmesh::LocalID> insertedBlocks;
         for (size_t b=0; b<refineList.size(); ++b) {
            cell->refine_block(refineList[b],insertedBlocks);
         }

         // Loop over blocks in map insertedBlocks and recalculate 
         // values of distribution functions
         for (map<vmesh::GlobalID,vmesh::LocalID>::const_iterator it=insertedBlocks.begin(); it!=insertedBlocks.end(); ++it) {
            const vmesh::GlobalID blockGID = it->first;
            const vmesh::LocalID blockLID = it->second;
            parameters = cell->get_block_parameters();
            creal vxBlock = parameters[blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VXCRD];
            creal vyBlock = parameters[blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VYCRD];
            creal vzBlock = parameters[blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VZCRD];
            creal dvxCell = parameters[blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVX];
            creal dvyCell = parameters[blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVY];
            creal dvzCell = parameters[blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVZ];

            for (uint kc=0; kc<WID; ++kc) {
               for (uint jc=0; jc<WID; ++jc) {
                  for (uint ic=0; ic<WID; ++ic) {
                     creal vxCell = vxBlock + ic*dvxCell;
                     creal vyCell = vyBlock + jc*dvyCell;
                     creal vzCell = vzBlock + kc*dvzCell;
                     creal average =
                       calcPhaseSpaceDensity(
                                             x, y, z, dx, dy, dz,
                                             vxCell,vyCell,vzCell,
                                             dvxCell,dvyCell,dvzCell);
                     cell->get_data(blockLID)[kc*WID2+jc*WID+ic] = average;
                  }
               }
            }
         }

         //refine = false;
         if (refineList.size() == 0) refine = false;
         ++currentLevel;
         if (currentLevel == Parameters::amrMaxVelocityRefLevel) refine = false;
      }
      delete refCriterion;
   }

   /*default one does not compute any parameters*/
   void Project::calcCellParameters(Real* cellParams,creal& t) { }
   
   Real Project::calcPhaseSpaceDensity(
      creal& x, creal& y, creal& z,
      creal& dx, creal& dy, creal& dz,
      creal& vx, creal& vy, creal& vz,
      creal& dvx, creal& dvy, creal& dvz) {
      cerr << "ERROR: Project::calcPhaseSpaceDensity called instead of derived class function!" << endl;
      return -1.0;
   }
   /*!
     Get random number between 0 and 1.0. One should always first initialize the rng.
   */
   
   Real Project::getRandomNumber() {
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
   
   void Project::setRandomSeed(CellID seedModifier) {
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
   void Project::setRandomCellSeed(const Real* const cellParams) {
      const creal x = cellParams[CellParams::XCRD];
      const creal y = cellParams[CellParams::YCRD];
      const creal z = cellParams[CellParams::ZCRD];
      const creal dx = cellParams[CellParams::DX];
      const creal dy = cellParams[CellParams::DY];
      const creal dz = cellParams[CellParams::DZ];
      
      const CellID cellID = (int) ((x - Parameters::xmin) / dx) +
         (int) ((y - Parameters::ymin) / dy) * Parameters::xcells_ini +
         (int) ((z - Parameters::zmin) / dz) * Parameters::xcells_ini * Parameters::ycells_ini;
      setRandomSeed(cellID);
   }


   
Project* createProject() {
   if(Parameters::projectName == "Alfven") {
      return new projects::Alfven;
   }
   if(Parameters::projectName == "Diffusion") {
      return new projects::Diffusion;
   }
   if(Parameters::projectName == "Dispersion") {
      return new projects::Dispersion;
   }
   if(Parameters::projectName == "Distributions") {
      return new projects::Distributions;
   }
   if(Parameters::projectName == "Firehose") {
      return new projects::Firehose;
   }
   if(Parameters::projectName == "Flowthrough") {
      return new projects::Flowthrough;
   }
   if(Parameters::projectName == "Fluctuations") {
         return new projects::Fluctuations;
   }

   if(Parameters::projectName == "Harris") {
      return new projects::Harris;
   }

   if(Parameters::projectName == "KHB") {
      return new projects::KHB;
   }
   if(Parameters::projectName == "Larmor") {
      return new projects::Larmor;
   }
   if(Parameters::projectName == "Magnetosphere") {
      return new projects::Magnetosphere;
   }
   if(Parameters::projectName == "MultiPeak") {
      return new projects::MultiPeak;
   } 
   if(Parameters::projectName == "VelocityBox") {
      return new projects::VelocityBox;
   } 
   if(Parameters::projectName == "Riemann1") {
      return new projects::Riemann1;
   }
   if(Parameters::projectName == "Shock") {
      return new projects::Shock;
   }
   if(Parameters::projectName == "Template") {
      return new projects::Template;
   }
   if(Parameters::projectName == "test_fp") {
      return new projects::test_fp;
   }
   if(Parameters::projectName == "testHall") {
      return new projects::TestHall;
   }
   if(Parameters::projectName == "test_trans") {
      return new projects::test_trans;
   }
   if(Parameters::projectName == "verificationLarmor") {
      return new projects::verificationLarmor;
   }
   if(Parameters::projectName == "Shocktest") {
      return new projects::Shocktest;
   }
   cerr << "Unknown project name!" << endl;
   abort();
}

} // namespace projects
