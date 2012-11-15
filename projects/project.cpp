#include "project.h"
#include <cstdlib>
#include "../vlasovmover.h"

#include "Alfven/Alfven.h"
#include "Diffusion/Diffusion.h"
#include "Dispersion/Dispersion.h"
#include "Firehose/Firehose.h"
#include "Flowthrough/Flowthrough.h"
#include "Fluctuations/Fluctuations.h"
#include "harm1D/harm1D.h"
#include "Harris/Harris.h"
#include "KelvinHelmholtz/KelvinHelmholtz.h"
#include "Larmor/Larmor.h"
#include "Magnetosphere/Magnetosphere.h"
#include "MultiPeak/MultiPeak.h"
#include "Riemann1/Riemann1.h"
#include "Shock/Shock.h"


using namespace std;

namespace projects {
   Project::Project() { }
   
   Project::~Project() { }
   
   void Project::addParameters() {
      typedef Readparameters RP;
      // TODO add all projects' static addParameters() functions here.
      projects::Alfven::addParameters();
      projects::Diffusion::addParameters();
      projects::Dispersion::addParameters();
      projects::Firehose::addParameters();
      projects::Flowthrough::addParameters();
      projects::Fluctuations::addParameters();
      projects::harm1D::addParameters();
      projects::Harris::addParameters();
      projects::KelvinHelmholtz::addParameters();
      projects::Larmor::addParameters();
      projects::Magnetosphere::addParameters();
      projects::MultiPeak::addParameters();
      projects::Riemann1::addParameters();
      projects::Shock::addParameters();
      RP::add("Project_common.seed", "Seed for the RNG", 42);
   }
   
   void Project::getParameters() {
      typedef Readparameters RP;
      RP::get("Fluctuations.seed", this->seed);
   }
   
   bool Project::initialize() {
      cerr << "ERROR: Project::initialize called instead of derived class function!" << endl;
      return false;
   }
   
   void Project::setCell(SpatialCell* cell) {
      // Set up cell parameters:
      this->calcCellParameters(&((*cell).parameters[0]), 0.0);
      
      cell->parameters[CellParams::RHOLOSSADJUST] = 0.0;
      cell->parameters[CellParams::RHOLOSSVELBOUNDARY] = 0.0;
      
      this->setVelocitySpace(cell);
      
      calculateCellVelocityMoments(cell);
      
      //let's get rid of blocks not fulfilling the criteria here to save memory.
      cell->adjustSingleCellVelocityBlocks();
   }
   
   vector<uint> Project::findBlocksToInitialize(SpatialCell* cell) {
      vector<uint> blocksToInitialize;
      
      for (uint kv=0; kv<P::vzblocks_ini; ++kv) 
         for (uint jv=0; jv<P::vyblocks_ini; ++jv)
            for (uint iv=0; iv<P::vxblocks_ini; ++iv) {
               creal vx = P::vxmin + (iv+0.5) * SpatialCell::block_dvx; // vx-coordinate of the centre
               creal vy = P::vymin + (jv+0.5) * SpatialCell::block_dvy; // vy-
               creal vz = P::vzmin + (kv+0.5) * SpatialCell::block_dvz; // vz-
               
               cell->add_velocity_block(cell->get_velocity_block(vx, vy, vz));
               blocksToInitialize.push_back(cell->get_velocity_block(vx, vy, vz));
      }
      
      return blocksToInitialize;
   }
   
   void Project::setVelocitySpace(SpatialCell* cell) {
      vector<uint> blocksToInitialize = this->findBlocksToInitialize(cell);
      
      for(uint i = 0; i < blocksToInitialize.size(); i++) {
         Velocity_Block* blockPtr = cell->at(blocksToInitialize.at(i));
         creal vxBlock = blockPtr->parameters[BlockParams::VXCRD];
         creal vyBlock = blockPtr->parameters[BlockParams::VYCRD];
         creal vzBlock = blockPtr->parameters[BlockParams::VZCRD];
         creal dvxCell = SpatialCell::cell_dvx; // Size of one cell in a block in vx-direction
         creal dvyCell = SpatialCell::cell_dvy; //                                vy
         creal dvzCell = SpatialCell::cell_dvz; //                                vz
         
         creal x = cell->parameters[CellParams::XCRD];
         creal y = cell->parameters[CellParams::YCRD];
         creal z = cell->parameters[CellParams::ZCRD];
         creal dx = cell->parameters[CellParams::DX];
         creal dy = cell->parameters[CellParams::DY];
         creal dz = cell->parameters[CellParams::DZ];
         
         // Calculate volume average of distrib. function for each cell in the block.
         for (uint kc=0; kc<WID; ++kc) 
            for (uint jc=0; jc<WID; ++jc) 
               for (uint ic=0; ic<WID; ++ic) {
                  creal vxCell = vxBlock + ic*dvxCell;
                  creal vyCell = vyBlock + jc*dvyCell;
                  creal vzCell = vzBlock + kc*dvzCell;
                  creal average =
                  this->calcPhaseSpaceDensity(
                     x, y, z, dx, dy, dz,
                     vxCell,vyCell,vzCell,
                     dvxCell,dvyCell,dvzCell);
                  
                  if(average!=0.0){
                     creal vxCellCenter = vxBlock + (ic+convert<Real>(0.5))*dvxCell;
                     creal vyCellCenter = vyBlock + (jc+convert<Real>(0.5))*dvyCell;
                     creal vzCellCenter = vzBlock + (kc+convert<Real>(0.5))*dvzCell;
                     cell->set_value(vxCellCenter,vyCellCenter,vzCellCenter,average);
                  }
         }
      }
   }
   
   void Project::calcCellParameters(Real* cellParams,creal& t) {
      cerr << "ERROR: Project::calcCellParameters called instead of derived class function!" << endl;
   }

   
   
   Real Project::calcPhaseSpaceDensity(
      creal& x, creal& y, creal& z,
      creal& dx, creal& dy, creal& dz,
      creal& vx, creal& vy, creal& vz,
      creal& dvx, creal& dvy, creal& dvz) {
      cerr << "ERROR: Project::calcPhaseSpaceDensity called instead of derived class function!" << endl;
      return -1.0;
   }
   
   Real Project::getRandomNumber() {
#ifdef _AIX
      int64_t rndInt;
#else
      int32_t rndInt;
#endif
      Real rnd;
      random_r(&rngDataBuffer,&rndInt);
      rnd=(Real)rndInt / RAND_MAX;
      return rnd;
   }
   
   void Project::setRandomSeed(uint64_t seedModifier) {
      memset(&(this->rngDataBuffer), 0, sizeof(this->rngDataBuffer));
#ifdef _AIX
      initstate_r(this->seed+seedModifier, &(this->rngStateBuffer[0]), 256, NULL, &(this->rngDataBuffer));
#else
      initstate_r(this->seed+seedModifier, &(this->rngStateBuffer[0]), 256, &(this->rngDataBuffer));   
#endif
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
   if(Parameters::projectName == "Firehose") {
      return new projects::Firehose;
   }
   if(Parameters::projectName == "Flowthrough") {
      return new projects::Flowthrough;
   }
   if(Parameters::projectName == "Fluctuations") {
         return new projects::Fluctuations;
   }
   if(Parameters::projectName == "harm1D") {
      return new projects::harm1D;
   }
   if(Parameters::projectName == "Harris") {
      return new projects::Harris;
   }
   if(Parameters::projectName == "KelvinHelmholtz") {
      return new projects::KelvinHelmholtz;
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
   if(Parameters::projectName == "Riemann1") {
      return new projects::Riemann1;
   }
   if(Parameters::projectName == "Shock") {
      return new projects::Shock;
   }
   cerr << "Unknown project name!" << endl;
   abort();
}

} // namespace projects
