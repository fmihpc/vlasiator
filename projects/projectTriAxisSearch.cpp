#include "projectTriAxisSearch.h"

namespace projects {
   /*!
    * WARNING This assumes that the velocity space is isotropic (same resolution in vx, vy, vz).
    */
   vector<uint> TriAxisSearch::findBlocksToInitialize(SpatialCell* cell) {
      set<uint> blocksToInitialize;
      bool search;
      int counter;
      
      creal x = cell->parameters[CellParams::XCRD];
      creal y = cell->parameters[CellParams::YCRD];
      creal z = cell->parameters[CellParams::ZCRD];
      creal dx = cell->parameters[CellParams::DX];
      creal dy = cell->parameters[CellParams::DY];
      creal dz = cell->parameters[CellParams::DZ];
      
      creal dvxCell = SpatialCell::get_velocity_grid_cell_size()[0];
      creal dvyCell = SpatialCell::get_velocity_grid_cell_size()[1];
      creal dvzCell = SpatialCell::get_velocity_grid_cell_size()[2];
      creal dvxBlock = SpatialCell::get_velocity_grid_block_size()[0];
      creal dvyBlock = SpatialCell::get_velocity_grid_block_size()[1];
      creal dvzBlock = SpatialCell::get_velocity_grid_block_size()[2];
      
      const vector<std::array<Real, 3>> V0 = this->getV0(x+0.5*dx, y+0.5*dy, z+0.5*dz);
      for (vector<std::array<Real, 3>>::const_iterator it = V0.begin(); it != V0.end(); it++) {
         // VX search
         search = true;
         counter = 0;
         while (search) {
            if (0.1 * P::sparseMinValue >
                calcPhaseSpaceDensity(x,
                                      y,
                                      z,
                                      dx,
                                      dy,
                                      dz,
                                      it->at(0) + counter*dvxBlock, it->at(1), it->at(2),
                                      dvxCell, dvyCell, dvzCell
                                     )
               ) {
               search = false;
            }
            ++counter;
            if (counter >= cell->get_velocity_grid_length(0)[0]) search = false;
         }
         counter+=2;
         Real vRadiusSquared = (Real)counter*(Real)counter*dvxBlock*dvxBlock;
         
         // VY search
         search = true;
         counter = 0;
         while(search) {
            if(0.1 * P::sparseMinValue >
               calcPhaseSpaceDensity(
                                     x,
                                     y,
                                     z,
                                     dx,
                                     dy,
                                     dz,
                                     it->at(0), it->at(1) + counter*dvyBlock, it->at(2),
                                     dvxCell, dvyCell, dvzCell
                                    )
               ||
               counter > P::vxblocks_ini
              ) {
               search = false;
            }
            ++counter;
            if (counter >= cell->get_velocity_grid_length(0)[1]) search = false;
         }
         counter+=2;
         vRadiusSquared = max(vRadiusSquared, (Real)counter*(Real)counter*dvyBlock*dvyBlock);

         // VZ search
         search = true;
         counter = 0;
         while(search) {
            if(0.1 * P::sparseMinValue >
               calcPhaseSpaceDensity(
                                     x,
                                     y,
                                     z,
                                     dx,
                                     dy,
                                     dz,
                                     it->at(0), it->at(1), it->at(2) + counter*dvzBlock,
                                     dvxCell, dvyCell, dvzCell
                                    )
               ||
               counter > P::vxblocks_ini
              ) {
               search = false;
            }
            ++counter;
            if (counter >= cell->get_velocity_grid_length(0)[2]) search = false;
         }
         counter+=2;
         vRadiusSquared = max(vRadiusSquared, (Real)counter*(Real)counter*dvzBlock*dvzBlock);
         
         // Block listing
         for (uint kv=0; kv<P::vzblocks_ini; ++kv) 
            for (uint jv=0; jv<P::vyblocks_ini; ++jv)
               for (uint iv=0; iv<P::vxblocks_ini; ++iv) {
                  creal vx = P::vxmin + (iv+0.5) * dvxBlock; // vx-coordinate of the centre
                  creal vy = P::vymin + (jv+0.5) * dvyBlock; // vy-
                  creal vz = P::vzmin + (kv+0.5) * dvzBlock; // vz-
                  
                  if ((vx-it->at(0))*(vx-it->at(0)) + (vy-it->at(1))*(vy-it->at(1)) + (vz-it->at(2))*(vz-it->at(2)) < vRadiusSquared) {
                     cell->add_velocity_block(cell->get_velocity_block(vx, vy, vz));
                     blocksToInitialize.insert(cell->get_velocity_block(vx, vy, vz));
                  }
               }
      }

      vector<uint> returnVector;
      for (set<uint>::const_iterator it = blocksToInitialize.begin(); it != blocksToInitialize.end(); it++) {
         returnVector.push_back(*it);
      }

      return returnVector;
   }
   
   vector<std::array<Real, 3>> TriAxisSearch::getV0(
      creal x,
      creal y,
      creal z
   ) {
      cerr << "ERROR: TriAxisSearch::getV0 called instead of derived class function!" << endl;
      abort();
      vector<std::array<Real, 3>> dummy;
      std::array<Real, 3> dudummy  {{0.0, 0.0, 0.0}};
      dummy.push_back(dudummy);
      return dummy;
   }
   
} // namespace projects