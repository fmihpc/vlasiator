/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://www.physics.helsinki.fi/vlasiator/
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#include <cstdlib>
#include <iostream>
#include <cmath>

#include "../../common.h"
#include "../../readparameters.h"
#include "../../backgroundfield/backgroundfield.h"
#include "../../object_wrapper.h"

#include "test_fp.h"

enum cases {BXCASE,BYCASE,BZCASE,BALLCASE};

using namespace std;

namespace projects {
   test_fp::test_fp(): TriAxisSearch() { }
   test_fp::~test_fp() { }


   /*typedef test_fpParameters tfP;
   Real this->B0 = NAN;
   Real this->DENSITY = NAN;
   Real this->TEMPERATURE = NAN;
   Real this->ALPHA = NAN;
   int  this->CASE = 5;
   bool this->shear = false;
   */

   bool test_fp::initialize(void) {
      Project::initialize();
      this->ALPHA *= M_PI / 4.0;
      return true;
   }

   void test_fp::addParameters(void){
      typedef Readparameters RP;
      RP::add("test_fp.V0", "Velocity magnitude (m/s)", 1.0e6);
      RP::add("test_fp.B0", "Magnetic field value in the non-zero patch (T)", 1.0e-9);
      RP::add("test_fp.rho", "Number density (m^-3)", 1.0e7);
      RP::add("test_fp.Temperature", "Temperature (K)", 1.0e-6);
      RP::add("test_fp.angle", "Orientation of the propagation expressed in pi/4", 0.0);
      RP::add("test_fp.Bdirection", "Direction of the magnetic field (0:x, 1:y, 2:z, 3:all)", 0);
      RP::add("test_fp.shear", "Add a shear (if false, V=0.5 everywhere).", true);
   }

   void test_fp::getParameters(void){
      Project::getParameters();
      typedef Readparameters RP;

      if(getObjectWrapper().particleSpecies.size() > 1) {
         std::cerr << "The selected project does not support multiple particle populations! Aborting in " << __FILE__ << " line " << __LINE__ << std::endl;
         abort();
      }
      RP::get("test_fp.B0", this->B0);
      RP::get("test_fp.V0", this->V0);
      RP::get("test_fp.rho", this->DENSITY);
      RP::get("test_fp.Temperature", this->TEMPERATURE);
      RP::get("test_fp.angle", this->ALPHA);
      RP::get("test_fp.Bdirection", this->CASE);
      RP::get("test_fp.shear", this->shear);
   }

   Real test_fp::sign(creal value) const {
      if (abs(value) < 1e-5) return 0.0;
      else return value / abs(value);
   }

   Real test_fp::calcPhaseSpaceDensity(creal& x,creal& y,creal& z,creal& dx,creal& dy,creal& dz,
                                       creal& vx,creal& vy,creal& vz,creal& dvx,creal& dvy,creal& dvz,
                                       const uint popID) const {      
      vector<std::array<Real, 3> > V = this->getV0(x,y,z,dx,dy,dz, popID);
      
      creal VX2 = (vx+0.5*dvx-V[0][0])*(vx+0.5*dvx-V[0][0]);
      creal VY2 = (vy+0.5*dvy-V[0][1])*(vy+0.5*dvy-V[0][1]);
      creal VZ2 = (vz+0.5*dvz-V[0][2])*(vz+0.5*dvz-V[0][2]);
      
      creal CONST = physicalconstants::MASS_PROTON / 2.0 / physicalconstants::K_B / this->TEMPERATURE;
      Real NORM = (physicalconstants::MASS_PROTON / 2.0 / M_PI / physicalconstants::K_B / this->TEMPERATURE);
      NORM = this->DENSITY * pow(NORM,1.5);
      
      creal result = NORM*exp(-CONST*(VX2+VY2+VZ2));
      return result;
   }
   
   void test_fp::setProjectBField(
      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2> & perBGrid,
      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, 2>& BgBGrid,
      FsGrid< fsgrids::technical, 2>& technicalGrid
   ) {
      setBackgroundFieldToZero(BgBGrid);
      
      if(!P::isRestart) {
         auto localSize = perBGrid.getLocalSize().data();
         
         creal dx = perBGrid.DX * 3.5;
         creal dy = perBGrid.DY * 3.5;
         creal dz = perBGrid.DZ * 3.5;
         
         Real areaFactor = 1.0;
         
         #pragma omp parallel for collapse(3)
         for (int i = 0; i < localSize[0]; ++i) {
            for (int j = 0; j < localSize[1]; ++j) {
               for (int k = 0; k < localSize[2]; ++k) {
                  const std::array<Real, 3> xyz = perBGrid.getPhysicalCoords(i, j, k);
                  std::array<Real, fsgrids::bfield::N_BFIELD>* cell = perBGrid.get(i, j, k);
                  
                  creal x = xyz[0] + 0.5 * perBGrid.DX;
                  creal y = xyz[1] + 0.5 * perBGrid.DY;
                  creal z = xyz[2] + 0.5 * perBGrid.DZ;
                  
                  switch (this->CASE) {
                     case BXCASE:         
                        cell->at(fsgrids::bfield::PERBX) = 0.1 * this->B0 * areaFactor;
                        //areaFactor = (CellParams::DY * CellParams::DZ) / (dy * dz);
                        if (y >= -dy && y <= dy)
                           if (z >= -dz && z <= dz)
                              cell->at(fsgrids::bfield::PERBX) = this->B0 * areaFactor;
                        break;
                     case BYCASE:
                        cell->at(fsgrids::bfield::PERBY) = 0.1 * this->B0 * areaFactor;
                        //areaFactor = (CellParams::DX * CellParams::DZ) / (dx * dz);
                        if (x >= -dx && x <= dx)
                           if (z >= -dz && z <= dz)
                              cell->at(fsgrids::bfield::PERBY) = this->B0 * areaFactor;
                        break;
                     case BZCASE:
                        cell->at(fsgrids::bfield::PERBZ) = 0.1 * this->B0 * areaFactor;
                        //areaFactor = (CellParams::DX * CellParams::DY) / (dx * dy);
                        if (x >= -dx && x <= dx)
                           if (y >= -dy && y <= dy)
                              cell->at(fsgrids::bfield::PERBZ) = this->B0 * areaFactor;
                        break;
                     case BALLCASE:
                        cell->at(fsgrids::bfield::PERBX) = 0.1 * this->B0 * areaFactor;
                        cell->at(fsgrids::bfield::PERBY) = 0.1 * this->B0 * areaFactor;
                        cell->at(fsgrids::bfield::PERBZ) = 0.1 * this->B0 * areaFactor;
                        
                        //areaFactor = (CellParams::DX * CellParams::DY) / (dx * dy);
                        
                        if (y >= -dy && y <= dy)
                           if (z >= -dz && z <= dz)
                              cell->at(fsgrids::bfield::PERBX) = this->B0 * areaFactor;
                        if (x >= -dx && x <= dx)
                           if (z >= -dz && z <= dz)
                              cell->at(fsgrids::bfield::PERBY) = this->B0 * areaFactor;
                        if (x >= -dx && x <= dx)
                           if (y >= -dy && y <= dy)
                              cell->at(fsgrids::bfield::PERBZ) = this->B0 * areaFactor;
                        break;
                  }
               }
            }
         }
      }
   }
   
   void test_fp::calcCellParameters(spatial_cell::SpatialCell* cell,creal& t) {
      
      typedef Parameters P;
      
   }
   
   vector<std::array<Real, 3>> test_fp::getV0(
      creal x,
      creal y,
      creal z,
      creal dx,
      creal dy,
      creal dz,
      const uint popID
   ) const {
      vector<std::array<Real, 3>> centerPoints;
      
      Real VX=0.0,VY=0.0,VZ=0.0;
      if (this->shear == true)
      {
         Real ksi,eta;
         switch (this->CASE) {
            case BXCASE:
               ksi = ((y + 0.5 * dy)  * cos(this->ALPHA) + (z + 0.5 * dz) * sin(this->ALPHA)) / (2.0 * sqrt(2.0));
               eta = (-(y + 0.5 * dy)  * sin(this->ALPHA) + (z + 0.5 * dz) * cos(this->ALPHA)) / (2.0 * sqrt(2.0));
               VX = 0.0;
               VY = sign(cos(this->ALPHA)) * 0.5 + 0.1*cos(this->ALPHA) * sin(2.0 * M_PI * eta);
               VZ = sign(sin(this->ALPHA)) * 0.5 + 0.1*sin(this->ALPHA) * sin(2.0 * M_PI * eta);
               break;
            case BYCASE:
               ksi = ((z + 0.5 * dz)  * cos(this->ALPHA) + (x + 0.5 * dx) * sin(this->ALPHA)) / (2.0 * sqrt(2.0));
               eta = (-(z + 0.5 * dz)  * sin(this->ALPHA) + (x + 0.5 * dx) * cos(this->ALPHA)) / (2.0 * sqrt(2.0));
               VX = sign(sin(this->ALPHA)) * 0.5 + 0.1*sin(this->ALPHA) * sin(2.0 * M_PI * eta);
               VY = 0.0;
               VZ = sign(cos(this->ALPHA)) * 0.5 + 0.1*cos(this->ALPHA) * sin(2.0 * M_PI * eta);
               break;
            case BZCASE:
               ksi = ((x + 0.5 * dx)  * cos(this->ALPHA) + (y + 0.5 * dy) * sin(this->ALPHA)) / (2.0 * sqrt(2.0));
               eta = (-(x + 0.5 * dx)  * sin(this->ALPHA) + (y + 0.5 * dy) * cos(this->ALPHA)) / (2.0 * sqrt(2.0));
               VX = sign(cos(this->ALPHA)) * 0.5 + 0.1*cos(this->ALPHA) * sin(2.0 * M_PI * eta);
               VY = sign(sin(this->ALPHA)) * 0.5 + 0.1*sin(this->ALPHA) * sin(2.0 * M_PI * eta);
               VZ = 0.0;
            break;
          case BALLCASE:
            std::cerr << "not implemented in " << __FILE__ << ":" << __LINE__ << std::endl;
            exit(1);
            break;
         }
      } else {
         switch (this->CASE) {
          case BXCASE:
            VX = 0.0;
            VY = cos(this->ALPHA) * 0.5;
            VZ = sin(this->ALPHA) * 0.5; 
            break;
          case BYCASE:
            VX = sin(this->ALPHA) * 0.5;
            VY = 0.0;
            VZ = cos(this->ALPHA) * 0.5;
            break;
          case BZCASE:
            VX = cos(this->ALPHA) * 0.5;
            VY = sin(this->ALPHA) * 0.5;
            VZ = 0.0;
            break;
          case BALLCASE:
            VX = 0.5 / sqrt(3.0);
            VY = 0.5 / sqrt(3.0);
            VZ = 0.5 / sqrt(3.0);
            break;
         }
      }
      
      VX *= this->V0 * 2.0;
      VY *= this->V0 * 2.0;
      VZ *= this->V0 * 2.0;
      
      std::array<Real, 3> point {{VX, VY, VZ}};
      centerPoints.push_back(point);
      return centerPoints;
   }
   
   vector<std::array<Real, 3>> test_fp::getV0(
      creal x,
      creal y,
      creal z,
      const uint popID
   ) const {
      vector<std::array<Real, 3>> centerPoints;
      
      creal dx = 0.0;
      creal dy = 0.0;
      creal dz = 0.0;
      
      return this->getV0(x,y,z,dx,dy,dz,popID);
   }

   bool test_fp::refineSpatialCells( dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid ) const {
 
     int myRank;       
     MPI_Comm_rank(MPI_COMM_WORLD,&myRank);

     // mpiGrid.set_maximum_refinement_level(std::min(this->maxSpatialRefinementLevel, mpiGrid.mapping.get_maximum_refinement_level()));

      // cout << "I am at line " << __LINE__ << " of " << __FILE__ <<  endl;
     if(myRank == MASTER_RANK) std::cout << "Maximum refinement level is " << mpiGrid.mapping.get_maximum_refinement_level() << std::endl;


      for (double x = P::amrBoxCenterX - P::amrBoxHalfWidthX * P::dx_ini; x <= P::amrBoxCenterX + P::amrBoxHalfWidthX * P::dx_ini; x += 0.99 * P::dx_ini) {
         for (double y = P::amrBoxCenterY - P::amrBoxHalfWidthY * P::dy_ini; y <= P::amrBoxCenterY + P::amrBoxHalfWidthY * P::dy_ini; y += 0.99 * P::dy_ini) {
            for (double z = P::amrBoxCenterZ - P::amrBoxHalfWidthZ * P::dz_ini; z <= P::amrBoxCenterZ + P::amrBoxHalfWidthZ * P::dz_ini; z += 0.99 * P::dz_ini) {
     
               std::array<double,3> xyz;
               xyz[0] = x;
               xyz[1] = y;
               xyz[2] = z;
               CellID myCell = mpiGrid.get_existing_cell(xyz);
               if (mpiGrid.refine_completely_at(xyz)) {
                  std::cout << "Rank " << myRank << " is refining cell " << myCell << std::endl;
               }
            }
         }
      }

      std::vector<CellID> refinedCells = mpiGrid.stop_refining(true);      
      if(myRank == MASTER_RANK) std::cout << "Finished first level of refinement" << endl;
      if(refinedCells.size() > 0) {
	std::cout << "Refined cells produced by rank " << myRank << " are: ";
	for (auto cellid : refinedCells) {
	  std::cout << cellid << " ";
	}
	std::cout << endl;
      }      
                  
      mpiGrid.balance_load();

//       auto cells = mpiGrid.get_cells();           
//       if(cells.empty()) {
//          std::cout << "Rank " << myRank << " has no cells!" << std::endl;
//       } else {
//          std::cout << "Cells on rank " << myRank << ": ";
//          for (auto c : cells) {
//             std::cout << c << " ";
//          }
//          std::cout << std::endl;
//       }

      return true;
   }

}// namespace projects
