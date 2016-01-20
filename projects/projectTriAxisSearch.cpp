#include "projectTriAxisSearch.h"
#include "../object_wrapper.h"

using namespace std;
using namespace spatial_cell;

using namespace std;

namespace projects {
   /*!
    * WARNING This assumes that the velocity space is isotropic (same resolution in vx, vy, vz).
    */
   std::vector<vmesh::GlobalID> TriAxisSearch::findBlocksToInitialize(SpatialCell* cell,const int& popID) const {
      set<vmesh::GlobalID> blocksToInitialize;
      bool search;
      int counter;
      
      creal x = cell->parameters[CellParams::XCRD];
      creal y = cell->parameters[CellParams::YCRD];
      creal z = cell->parameters[CellParams::ZCRD];
      creal dx = cell->parameters[CellParams::DX];
      creal dy = cell->parameters[CellParams::DY];
      creal dz = cell->parameters[CellParams::DZ];
      
      const uint8_t refLevel = 0;
      creal dvxCell = cell->get_velocity_grid_cell_size(popID,refLevel)[0];
      creal dvyCell = cell->get_velocity_grid_cell_size(popID,refLevel)[1];
      creal dvzCell = cell->get_velocity_grid_cell_size(popID,refLevel)[2];
      creal dvxBlock = cell->get_velocity_grid_block_size(popID,refLevel)[0];
      creal dvyBlock = cell->get_velocity_grid_block_size(popID,refLevel)[1];
      creal dvzBlock = cell->get_velocity_grid_block_size(popID,refLevel)[2];
      
      const size_t vxblocks_ini = cell->get_velocity_grid_length(popID,refLevel)[0];
      const size_t vyblocks_ini = cell->get_velocity_grid_length(popID,refLevel)[1];
      const size_t vzblocks_ini = cell->get_velocity_grid_length(popID,refLevel)[2];

      const vector<std::array<Real, 3>> V0 = this->getV0(x+0.5*dx, y+0.5*dy, z+0.5*dz);
      for (vector<std::array<Real, 3>>::const_iterator it = V0.begin(); it != V0.end(); it++) {
         // VX search
         search = true;
         counter = 0;
         #warning TODO: add SpatialCell::getVelocityBlockMinValue() in place of sparseMinValue
         while (search) {
            if (0.1 * getObjectWrapper().particleSpecies[popID].sparseMinValue >
                calcPhaseSpaceDensity(x,
                                      y,
                                      z,
                                      dx,
                                      dy,
                                      dz,
                                      it->at(0) + counter*dvxBlock, it->at(1), it->at(2),
                                      dvxCell, dvyCell, dvzCell, popID
                                     )
               ) {
               search = false;
            }
            ++counter;
            if (counter >= cell->get_velocity_grid_length(popID,refLevel)[0]) search = false;
         }
         counter+=2;
         Real vRadiusSquared = (Real)counter*(Real)counter*dvxBlock*dvxBlock;

         // VY search
         search = true;
         counter = 0;
         while(search) {
            if (0.1 * getObjectWrapper().particleSpecies[popID].sparseMinValue >
               calcPhaseSpaceDensity(
                                     x,
                                     y,
                                     z,
                                     dx,
                                     dy,
                                     dz,
                                     it->at(0), it->at(1) + counter*dvyBlock, it->at(2),
                                     dvxCell, dvyCell, dvzCell, popID
                                    )
               ||
               counter > vxblocks_ini
              ) {
               search = false;
            }
            ++counter;
            if (counter >= cell->get_velocity_grid_length(popID,refLevel)[1]) search = false;
         }
         counter+=2;
         vRadiusSquared = max(vRadiusSquared, (Real)counter*(Real)counter*dvyBlock*dvyBlock);

         // VZ search
         search = true;
         counter = 0;
         while(search) {
            if (0.1 * getObjectWrapper().particleSpecies[popID].sparseMinValue >
               calcPhaseSpaceDensity(
                                     x,
                                     y,
                                     z,
                                     dx,
                                     dy,
                                     dz,
                                     it->at(0), it->at(1), it->at(2) + counter*dvzBlock,
                                     dvxCell, dvyCell, dvzCell, popID
                                    )
               ||
               counter > vxblocks_ini
              ) {
               search = false;
            }
            ++counter;
            if (counter >= cell->get_velocity_grid_length(popID,refLevel)[2]) search = false;
         }
         counter+=2;
         vRadiusSquared = max(vRadiusSquared, (Real)counter*(Real)counter*dvzBlock*dvzBlock);

         // Block listing
         for (uint kv=0; kv<vzblocks_ini; ++kv) 
            for (uint jv=0; jv<vyblocks_ini; ++jv)
               for (uint iv=0; iv<vxblocks_ini; ++iv) {
                  vmesh::GlobalID blockIndices[3];
                  blockIndices[0] = iv;
                  blockIndices[1] = jv;
                  blockIndices[2] = kv;
                  const vmesh::GlobalID blockGID = cell->get_velocity_block(popID,blockIndices,refLevel);
                  
                  Real V_crds[3];
                  cell->get_velocity_block_coordinates(popID,blockGID,V_crds);
                  Real dV[3];
                  cell->get_velocity_block_size(popID,blockGID,dV);
                  V_crds[0] += 0.5*dV[0];
                  V_crds[1] += 0.5*dV[1];
                  V_crds[2] += 0.5*dV[2];
                  Real R2 = ((V_crds[0]-it->at(0))*(V_crds[0]-it->at(0))
                          + (V_crds[1]-it->at(1))*(V_crds[1]-it->at(1))
                          + (V_crds[2]-it->at(2))*(V_crds[2]-it->at(2)));
                  
                  if (R2 < vRadiusSquared) {
                     cell->add_velocity_block(blockGID,popID);
                     blocksToInitialize.insert(blockGID);
                  }
               }
      }

      vector<vmesh::GlobalID> returnVector;
      for (set<vmesh::GlobalID>::const_iterator it=blocksToInitialize.begin(); it!=blocksToInitialize.end(); ++it) {
         returnVector.push_back(*it);
      }

      return returnVector;
   }
   
   vector<std::array<Real, 3>> TriAxisSearch::getV0(
      creal x,
      creal y,
      creal z
   ) const {
      cerr << "ERROR: TriAxisSearch::getV0 called instead of derived class function!" << endl;
      abort();
      vector<std::array<Real, 3>> dummy;
      std::array<Real, 3> dudummy  {{0.0, 0.0, 0.0}};
      dummy.push_back(dudummy);
      return dummy;
   }
   
} // namespace projects
