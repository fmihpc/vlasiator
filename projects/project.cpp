#include "project.h"
#include "../vlasovmover.h"

#include "Alfven/Alfven.h"
#include "Diffusion/Diffusion.h"
#include "Dispersion/Dispersion.h"
#include "Magnetosphere/Magnetosphere.h"

using namespace std;

namespace projects {
   Project::Project() { }
   
   Project::~Project() { }
   
   void Project::addParameters() {
      // TODO add all projects' static addParameters() functions here.
      projects::Alfven::addParameters();
      projects::Diffusion::addParameters();
      projects::Dispersion::addParameters();
      projects::Magnetosphere::addParameters();
   }
   
   void Project::getParameters() {
      cerr << "ERROR: Project::getParameters called instead of derived class function!" << endl;
   }
   
   bool Project::initialize() {
      cerr << "ERROR: Project::initialize called instead of derived class function!" << endl;
      return false;
   }
   
   void Project::setCell(SpatialCell* cell) {
      // Set up cell parameters:
      calcCellParameters(&((*cell).parameters[0]), 0.0);
      
      cell->parameters[CellParams::RHOLOSSADJUST] = 0.0;
      cell->parameters[CellParams::RHOLOSSVELBOUNDARY] = 0.0;
      
      loopThroughFullVelocitySpace(cell);
      
      calculateCellVelocityMoments(cell);
      
      //let's get rid of blocks not fulfilling the criteria here to save memory.
      cell->adjustSingleCellVelocityBlocks();
   }
   
   void Project::loopThroughFullVelocitySpace(SpatialCell* cell) {
      // Go through each velocity block in the velocity phase space grid.
      // Set the initial state and block parameters:
      creal dvx_block = SpatialCell::block_dvx; // Size of a block in vx-direction
      creal dvy_block = SpatialCell::block_dvy; //                    vy
      creal dvz_block = SpatialCell::block_dvz; //                    vz
      creal dvx_blockCell = SpatialCell::cell_dvx; // Size of one cell in a block in vx-direction
      creal dvy_blockCell = SpatialCell::cell_dvy; //                                vy
      creal dvz_blockCell = SpatialCell::cell_dvz; //                                vz
      
      for (uint kv=0; kv<P::vzblocks_ini; ++kv) 
         for (uint jv=0; jv<P::vyblocks_ini; ++jv)
            for (uint iv=0; iv<P::vxblocks_ini; ++iv) {
               creal vx_block = P::vxmin + iv*dvx_block; // vx-coordinate of the lower left corner
               creal vy_block = P::vymin + jv*dvy_block; // vy-
               creal vz_block = P::vzmin + kv*dvz_block; // vz-
               
               // Calculate volume average of distrib. function for each cell in the block.
               for (uint kc=0; kc<WID; ++kc) 
                  for (uint jc=0; jc<WID; ++jc) 
                     for (uint ic=0; ic<WID; ++ic) {
                        creal vx_cell = vx_block + ic*dvx_blockCell;
                        creal vy_cell = vy_block + jc*dvy_blockCell;
                        creal vz_cell = vz_block + kc*dvz_blockCell;
                        Real average = 
                        calcPhaseSpaceDensity(cell->parameters[CellParams::XCRD],
                                              cell->parameters[CellParams::YCRD],
                                              cell->parameters[CellParams::ZCRD],
                                              cell->parameters[CellParams::DX],
                                              cell->parameters[CellParams::DY],
                                              cell->parameters[CellParams::DZ],
                                              vx_cell,vy_cell,vz_cell,
                                              dvx_blockCell,dvy_blockCell,dvz_blockCell);
                        
                        if(average!=0.0){
                           creal vx_cell_center = vx_block + (ic+convert<Real>(0.5))*dvx_blockCell;
                           creal vy_cell_center = vy_block + (jc+convert<Real>(0.5))*dvy_blockCell;
                           creal vz_cell_center = vz_block + (kc+convert<Real>(0.5))*dvz_blockCell;
                           cell->set_value(vx_cell_center,vy_cell_center,vz_cell_center,average);
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
    if(Parameters::projectName == "Magnetosphere") {
       return new projects::Magnetosphere;
    }
    cerr << "Unknown project name!" << endl;
    abort();
 }
   
} // namespace projects