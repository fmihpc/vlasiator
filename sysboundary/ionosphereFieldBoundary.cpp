#include "ionosphereFieldBoundary.h"
#include "sysboundarycondition.h"
#include "../fieldsolver/fs_common.h"

using namespace std;

extern ARCH_MANAGED GridParameters meshParams;

namespace SBC {

   /*! We want here to
    *
    * -- Retain only the boundary-normal projection of perturbed face B
    */
   void IonosphereFieldBoundary::fieldSolverBoundaryCondMagneticFieldProjection(
      FsGrid<Real, fsgrids::bfield::N_BFIELD, FS_STENCIL_WIDTH> & bGrid,
      FsGrid< fsgrids::technical, 1, FS_STENCIL_WIDTH> & technicalGrid,
      cint i,
      cint j,
      cint k
   ) {
      // Projection of B-field to normal direction
      Real BdotN = 0;
      std::array<Real, 3> normalDirection = fieldSolverGetNormalDirection(technicalGrid, i, j, k);
      for(uint component=0; component<3; component++) {
         BdotN += bGrid.get(i,j,k)[fsgrids::bfield::PERBX+component] * normalDirection[component];
      }
      // Apply to any components that were not solved
      if ((technicalGrid.get(i,j,k,0).sysBoundaryLayer == 2) ||
          ((technicalGrid.get(i,j,k,0).sysBoundaryLayer == 1) && ((technicalGrid.get(i,j,k,0).SOLVE & compute::BX) != compute::BX))
         ) {
         bGrid.get(i,j,k)[fsgrids::bfield::PERBX] = BdotN*normalDirection[0];
      }
      if ((technicalGrid.get(i,j,k,0).sysBoundaryLayer == 2) ||
          ((technicalGrid.get(i,j,k,0).sysBoundaryLayer == 1) && ((technicalGrid.get(i,j,k,0).SOLVE & compute::BY) != compute::BY))
         ) {
         bGrid.get(i,j,k)[fsgrids::bfield::PERBY] = BdotN*normalDirection[1];
      }
      if ((technicalGrid.get(i,j,k,0).sysBoundaryLayer == 2) ||
          ((technicalGrid.get(i,j,k,0).sysBoundaryLayer == 1) && ((technicalGrid.get(i,j,k,0).SOLVE & compute::BZ) != compute::BZ))
         ) {
         bGrid.get(i,j,k)[fsgrids::bfield::PERBZ] = BdotN*normalDirection[2];
      }
   }

   void IonosphereFieldBoundary::fieldSolverBoundaryCondElectricField(
      FsGrid<Real, fsgrids::efield::N_EFIELD, FS_STENCIL_WIDTH> & EGrid,
      cint i,
      cint j,
      cint k,
      cuint component
   ) {
      EGrid.get(i,j,k)[fsgrids::efield::EX+component] = 0.0;
   }
   
   void IonosphereFieldBoundary::fieldSolverBoundaryCondHallElectricField(
      FsGrid<Real, fsgrids::ehall::N_EHALL, FS_STENCIL_WIDTH> & EHallGrid,
      cint i,
      cint j,
      cint k,
      cuint component
   ) {
      auto cp = EHallGrid.get(i,j,k);
      switch (component) {
         case 0:
            cp[fsgrids::ehall::EXHALL_000_100] = 0.0;
            cp[fsgrids::ehall::EXHALL_010_110] = 0.0;
            cp[fsgrids::ehall::EXHALL_001_101] = 0.0;
            cp[fsgrids::ehall::EXHALL_011_111] = 0.0;
            break;
         case 1:
            cp[fsgrids::ehall::EYHALL_000_010] = 0.0;
            cp[fsgrids::ehall::EYHALL_100_110] = 0.0;
            cp[fsgrids::ehall::EYHALL_001_011] = 0.0;
            cp[fsgrids::ehall::EYHALL_101_111] = 0.0;
            break;
         case 2:
            cp[fsgrids::ehall::EZHALL_000_001] = 0.0;
            cp[fsgrids::ehall::EZHALL_100_101] = 0.0;
            cp[fsgrids::ehall::EZHALL_010_011] = 0.0;
            cp[fsgrids::ehall::EZHALL_110_111] = 0.0;
            break;
         default:
            cerr << __FILE__ << ":" << __LINE__ << ":" << " Invalid component" << endl;
      }
   }
   
   void IonosphereFieldBoundary::fieldSolverBoundaryCondGradPeElectricField(
      FsGrid<Real, fsgrids::egradpe::N_EGRADPE, FS_STENCIL_WIDTH> & EGradPeGrid,
      cint i,
      cint j,
      cint k,
      cuint component
   ) {
      EGradPeGrid.get(i,j,k)[fsgrids::egradpe::EXGRADPE+component] = 0.0;
   }
   
   void IonosphereFieldBoundary::fieldSolverBoundaryCondDerivatives(
      FsGrid<Real, fsgrids::dperb::N_DPERB, FS_STENCIL_WIDTH> & dPerBGrid,
      FsGrid<Real, fsgrids::dmoments::N_DMOMENTS, FS_STENCIL_WIDTH> & dMomentsGrid,
      cint i,
      cint j,
      cint k,
      cuint& RKCase,
      cuint& component
   ) {
      SysBoundaryCondition::setCellDerivativesToZero(dPerBGrid, dMomentsGrid, i, j, k, component);
      return;
   }
   
   void IonosphereFieldBoundary::fieldSolverBoundaryCondBVOLDerivatives(
      FsGrid<Real, fsgrids::volfields::N_VOL, FS_STENCIL_WIDTH> & volGrid,
      cint i,
      cint j,
      cint k,
      cuint& component
   ) {
      // FIXME This should be OK as the BVOL derivatives are only used for Lorentz force JXB, which is not applied on the ionosphere cells.
      SysBoundaryCondition::setCellBVOLDerivativesToZero(volGrid, i, j, k, component);
   }

    std::array<Real, 3> IonosphereFieldBoundary::fieldSolverGetNormalDirection(
      FsGrid< fsgrids::technical, 1, FS_STENCIL_WIDTH> & technicalGrid,
      cint i,
      cint j,
      cint k
   ) {
      phiprof::start("Ionosphere::fieldSolverGetNormalDirection");
      std::array<Real, 3> normalDirection{{ 0.0, 0.0, 0.0 }};
      
      static creal DIAG2 = 1.0 / sqrt(2.0);
      static creal DIAG3 = 1.0 / sqrt(3.0);
      
      creal dx = technicalGrid.DX;
      creal dy = technicalGrid.DY;
      creal dz = technicalGrid.DZ;
      int globalIndices[3];
      technicalGrid.getGlobalIndices(i,j,k, globalIndices);
      creal x = meshParams.xmin + (convert<Real>(globalIndices[0])+0.5)*dx;
      creal y = meshParams.ymin + (convert<Real>(globalIndices[1])+0.5)*dy;
      creal z = meshParams.zmin + (convert<Real>(globalIndices[2])+0.5)*dz;
      creal xsign = divideIfNonZero(x, fabs(x));
      creal ysign = divideIfNonZero(y, fabs(y));
      creal zsign = divideIfNonZero(z, fabs(z));
      
      Real length = 0.0;
      
      if (meshParams.xcells_ini == 1) {
         if (meshParams.ycells_ini == 1) {
            if (meshParams.zcells_ini == 1) {
               // X,Y,Z
               std::cerr << __FILE__ << ":" << __LINE__ << ":" << "What do you expect to do with a single-cell simulation of ionosphere boundary type? Stop kidding." << std::endl;
               abort();
               // end of X,Y,Z
            } else {
               // X,Y
               normalDirection[2] = zsign;
               // end of X,Y
            }
         } else if (meshParams.zcells_ini == 1) {
            // X,Z
            normalDirection[1] = ysign;
            // end of X,Z
         } else {
            // X
            switch(this->geometry) {
               case 0:
                  normalDirection[1] = DIAG2*ysign;
                  normalDirection[2] = DIAG2*zsign;
                  break;
               case 1:
                  if(fabs(y) == fabs(z)) {
                     normalDirection[1] = ysign*DIAG2;
                     normalDirection[2] = zsign*DIAG2;
                     break;
                  }
                  if(fabs(y) > (this->radius - dy)) {
                     normalDirection[1] = ysign;
                     break;
                  }
                  if(fabs(z) > (this->radius - dz)) {
                     normalDirection[2] = zsign;
                     break;
                  }
                  if(fabs(y) > (this->radius - 2.0*dy)) {
                     normalDirection[1] = ysign;
                     break;
                  }
                  if(fabs(z) > (this->radius - 2.0*dz)) {
                     normalDirection[2] = zsign;
                     break;
                  }
                  break;
               case 2:
                  length = sqrt(y*y + z*z);
                  normalDirection[1] = y / length;
                  normalDirection[2] = z / length;
                  break;
               default:
                  std::cerr << __FILE__ << ":" << __LINE__ << ":" << "ionosphere.geometry has to be 0, 1 or 2 with this grid shape." << std::endl;
                  abort();
            }
            // end of X
         }
      } else if (meshParams.ycells_ini == 1) {
         if (meshParams.zcells_ini == 1) {
            // Y,Z
            normalDirection[0] = xsign;
            // end of Y,Z
         } else {
            // Y
            switch(this->geometry) {
               case 0:
                  normalDirection[0] = DIAG2*xsign;
                  normalDirection[2] = DIAG2*zsign;
                  break;
               case 1:
                  if(fabs(x) == fabs(z)) {
                     normalDirection[0] = xsign*DIAG2;
                     normalDirection[2] = zsign*DIAG2;
                     break;
                  }
                  if(fabs(x) > (this->radius - dx)) {
                     normalDirection[0] = xsign;
                     break;
                  }
                  if(fabs(z) > (this->radius - dz)) {
                     normalDirection[2] = zsign;
                     break;
                  }
                  if(fabs(x) > (this->radius - 2.0*dx)) {
                     normalDirection[0] = xsign;
                     break;
                  }
                  if(fabs(z) > (this->radius - 2.0*dz)) {
                     normalDirection[2] = zsign;
                     break;
                  }
                  break;
               case 2:
               case 3:
                  length = sqrt(x*x + z*z);
                  normalDirection[0] = x / length;
                  normalDirection[2] = z / length;
                  break;
               default:
                  std::cerr << __FILE__ << ":" << __LINE__ << ":" << "ionosphere.geometry has to be 0, 1, 2 or 3 with this grid shape." << std::endl;
                  abort();
            }
            // end of Y
         }
      } else if (meshParams.zcells_ini == 1) {
         // Z
         switch(this->geometry) {
            case 0:
               normalDirection[0] = DIAG2*xsign;
               normalDirection[1] = DIAG2*ysign;
               break;
            case 1:
               if(fabs(x) == fabs(y)) {
                  normalDirection[0] = xsign*DIAG2;
                  normalDirection[1] = ysign*DIAG2;
                  break;
               }
               if(fabs(x) > (this->radius - dx)) {
                  normalDirection[0] = xsign;
                  break;
               }
               if(fabs(y) > (this->radius - dy)) {
                  normalDirection[1] = ysign;
                  break;
               }
               if(fabs(x) > (this->radius - 2.0*dx)) {
                  normalDirection[0] = xsign;
                  break;
               }
               if(fabs(y) > (this->radius - 2.0*dy)) {
                  normalDirection[1] = ysign;
                  break;
               }
               break;
            case 2:
               length = sqrt(x*x + y*y);
               normalDirection[0] = x / length;
               normalDirection[1] = y / length;
               break;
            default:
               std::cerr << __FILE__ << ":" << __LINE__ << ":" << "ionosphere.geometry has to be 0, 1 or 2 with this grid shape." << std::endl;
               abort();
         }
         // end of Z
      } else {
         // 3D
         switch(this->geometry) {
            case 0:
               normalDirection[0] = DIAG3*xsign;
               normalDirection[1] = DIAG3*ysign;
               normalDirection[2] = DIAG3*zsign;
               break;
            case 1:
               if(fabs(x) == fabs(y) && fabs(x) == fabs(z) && fabs(x) > this->radius - dx) {
                  normalDirection[0] = xsign*DIAG3;
                  normalDirection[1] = ysign*DIAG3;
                  normalDirection[2] = zsign*DIAG3;
                  break;
               }
               if(fabs(x) == fabs(y) && fabs(x) == fabs(z) && fabs(x) > this->radius - 2.0*dx) {
                  normalDirection[0] = xsign*DIAG3;
                  normalDirection[1] = ysign*DIAG3;
                  normalDirection[2] = zsign*DIAG3;
                  break;
               }
               if(fabs(x) == fabs(y) && fabs(x) > this->radius - dx && fabs(z) < this->radius - dz) {
                  normalDirection[0] = xsign*DIAG2;
                  normalDirection[1] = ysign*DIAG2;
                  normalDirection[2] = 0.0;
                  break;
               }
               if(fabs(y) == fabs(z) && fabs(y) > this->radius - dy && fabs(x) < this->radius - dx) {
                  normalDirection[0] = 0.0;
                  normalDirection[1] = ysign*DIAG2;
                  normalDirection[2] = zsign*DIAG2;
                  break;
               }
               if(fabs(x) == fabs(z) && fabs(x) > this->radius - dx && fabs(y) < this->radius - dy) {
                  normalDirection[0] = xsign*DIAG2;
                  normalDirection[1] = 0.0;
                  normalDirection[2] = zsign*DIAG2;
                  break;
               }
               if(fabs(x) == fabs(y) && fabs(x) > this->radius - 2.0*dx && fabs(z) < this->radius - 2.0*dz) {
                  normalDirection[0] = xsign*DIAG2;
                  normalDirection[1] = ysign*DIAG2;
                  normalDirection[2] = 0.0;
                  break;
               }
               if(fabs(y) == fabs(z) && fabs(y) > this->radius - 2.0*dy && fabs(x) < this->radius - 2.0*dx) {
                  normalDirection[0] = 0.0;
                  normalDirection[1] = ysign*DIAG2;
                  normalDirection[2] = zsign*DIAG2;
                  break;
               }
               if(fabs(x) == fabs(z) && fabs(x) > this->radius - 2.0*dx && fabs(y) < this->radius - 2.0*dy) {
                  normalDirection[0] = xsign*DIAG2;
                  normalDirection[1] = 0.0;
                  normalDirection[2] = zsign*DIAG2;
                  break;
               }
               if(fabs(x) > (this->radius - dx)) {
                  normalDirection[0] = xsign;
                  break;
               }
               if(fabs(y) > (this->radius - dy)) {
                  normalDirection[1] = ysign;
                  break;
               }
               if(fabs(z) > (this->radius - dz)) {
                  normalDirection[2] = zsign;
                  break;
               }
               if(fabs(x) > (this->radius - 2.0*dx)) {
                  normalDirection[0] = xsign;
                  break;
               }
               if(fabs(y) > (this->radius - 2.0*dy)) {
                  normalDirection[1] = ysign;
                  break;
               }
               if(fabs(z) > (this->radius - 2.0*dz)) {
                  normalDirection[2] = zsign;
                  break;
               }
               break;
            case 2:
               length = sqrt(x*x + y*y + z*z);
               normalDirection[0] = x / length;
               normalDirection[1] = y / length;
               normalDirection[2] = z / length;
               break;
            case 3:
               length = sqrt(x*x + z*z);
               normalDirection[0] = x / length;
               normalDirection[2] = z / length;
               break;
            default:
               std::cerr << __FILE__ << ":" << __LINE__ << ":" << "ionosphere.geometry has to be 0, 1, 2 or 3 with this grid shape." << std::endl;
               abort();
         }
         // end of 3D
      }
      
      phiprof::stop("Ionosphere::fieldSolverGetNormalDirection");
      return normalDirection;
   }


}