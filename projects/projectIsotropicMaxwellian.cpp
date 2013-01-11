#include "projectIsotropicMaxwellian.h"

namespace projects {
   /*!
    * WARNING This assumes that the velocity space is isotropic (same resolution in vx, vy, vz).
    */
   vector<uint> IsotropicMaxwellian::findBlocksToInitialize(SpatialCell* cell) {
      vector<uint> blocksToInitialize;
      bool search = true;
      int counter = 0;
      
      creal x = cell->parameters[CellParams::XCRD];
      creal y = cell->parameters[CellParams::YCRD];
      creal z = cell->parameters[CellParams::ZCRD];
      creal dx = cell->parameters[CellParams::DX];
      creal dy = cell->parameters[CellParams::DY];
      creal dz = cell->parameters[CellParams::DZ];
      
      creal dvxCell = SpatialCell::cell_dvx; // Size of one cell in a block in vx-direction
      creal dvyCell = SpatialCell::cell_dvy; //                                vy
      creal dvzCell = SpatialCell::cell_dvz; //                                vz
      
      creal V0[3] = {
         this->getV0(x+0.5*dx, y+0.5*dy, z+0.5*dz, 0),
         this->getV0(x+0.5*dx, y+0.5*dy, z+0.5*dz, 1),
         this->getV0(x+0.5*dx, y+0.5*dy, z+0.5*dz, 2) 
      };
      
      while(search) {
         if(0.1 * P::sparseMinValue >
            calcPhaseSpaceDensity(
               x,
               y,
               z,
               dx,
               dy,
               dz,
               V0[0] + counter*SpatialCell::block_dvx, V0[1], V0[2],
               dvxCell, dvyCell, dvzCell
            )
         ) {
            search = false;
         }
         counter++;
      }
      counter+=2;
      Real vRadiusSquared = (Real)counter*(Real)counter*SpatialCell::block_dvx*SpatialCell::block_dvx;
      
      for (uint kv=0; kv<P::vzblocks_ini; ++kv) 
         for (uint jv=0; jv<P::vyblocks_ini; ++jv)
            for (uint iv=0; iv<P::vxblocks_ini; ++iv) {
               creal vx = P::vxmin + (iv+0.5) * SpatialCell::block_dvx; // vx-coordinate of the centre
               creal vy = P::vymin + (jv+0.5) * SpatialCell::block_dvy; // vy-
               creal vz = P::vzmin + (kv+0.5) * SpatialCell::block_dvz; // vz-
               
               if((vx-V0[0])*(vx-V0[0]) + (vy-V0[1])*(vy-V0[1]) + (vz-V0[2])*(vz-V0[2]) < vRadiusSquared) {
                  cell->add_velocity_block(cell->get_velocity_block(vx, vy, vz));
                  blocksToInitialize.push_back(cell->get_velocity_block(vx, vy, vz));
               }
      }
      
      return blocksToInitialize;
   }
   
   Real IsotropicMaxwellian::getV0(
      creal x,
      creal y,
      creal z,
      cuint component
   ) {
      cerr << "ERROR: projectIsotropicMaxwellian::getV0 called instead of derived class function!" << endl;
      abort();
      return 0.0;
   }
   
} // namespace projects