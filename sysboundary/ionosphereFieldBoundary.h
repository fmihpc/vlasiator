#include "../common.h"
#include "../fieldsolver/fs_common.h"
#include "../parameters.h"
#include <iostream>

using namespace std;

namespace SBC {
    class IonosphereFieldBoundary {

    public:
    IonosphereFieldBoundary(Real center[3], Real radius, uint geometry): center{center[0], center[1], center[2]}, radius{radius}, geometry{geometry} {}
 
    /*! We want here to
    * 
    * -- Average perturbed face B from the nearest neighbours
    */
    ARCH_HOSTDEV Real fieldSolverBoundaryCondMagneticField(
        const arch::buf<FsGrid<Real, fsgrids::bfield::N_BFIELD, FS_STENCIL_WIDTH>> & bGrid,
        const arch::buf<FsGrid< fsgrids::technical, 1, FS_STENCIL_WIDTH>> & technicalGrid,
        cint i,
        cint j,
        cint k,
        creal& dt,
        cuint& component
    ) {
        if (technicalGrid.get(i,j,k)->sysBoundaryLayer == 1) {
            switch(component) {
                case 0:
                if (  ((technicalGrid.get(i-1,j,k)->SOLVE & compute::BX) == compute::BX)
                    && ((technicalGrid.get(i+1,j,k)->SOLVE & compute::BX) == compute::BX)
                ) {
                    return 0.5 * (bGrid.get(i-1,j,k)[fsgrids::bfield::PERBX] + bGrid.get(i+1,j,k)[fsgrids::bfield::PERBX]);
                } else if ((technicalGrid.get(i-1,j,k)->SOLVE & compute::BX) == compute::BX) {
                    return bGrid.get(i-1,j,k)[fsgrids::bfield::PERBX];
                } else if ((technicalGrid.get(i+1,j,k)->SOLVE & compute::BX) == compute::BX) {
                    return bGrid.get(i+1,j,k)[fsgrids::bfield::PERBX];
                } else {
                    Real retval = 0.0;
                    uint nCells = 0;
                    if ((technicalGrid.get(i,j-1,k)->SOLVE & compute::BX) == compute::BX) {
                        retval += bGrid.get(i,j-1,k)[fsgrids::bfield::PERBX];
                        nCells++;
                    }
                    if ((technicalGrid.get(i,j+1,k)->SOLVE & compute::BX) == compute::BX) {
                        retval += bGrid.get(i,j+1,k)[fsgrids::bfield::PERBX];
                        nCells++;
                    }
                    if ((technicalGrid.get(i,j,k-1)->SOLVE & compute::BX) == compute::BX) {
                        retval += bGrid.get(i,j,k-1)[fsgrids::bfield::PERBX];
                        nCells++;
                    }
                    if ((technicalGrid.get(i,j,k+1)->SOLVE & compute::BX) == compute::BX) {
                        retval += bGrid.get(i,j,k+1)[fsgrids::bfield::PERBX];
                        nCells++;
                    }
                    if (nCells == 0) {
                        for (int a=i-1; a<i+2; a++) {
                            for (int b=j-1; b<j+2; b++) {
                            for (int c=k-1; c<k+2; c++) {
                                if ((technicalGrid.get(a,b,c)->SOLVE & compute::BX) == compute::BX) {
                                    retval += bGrid.get(a,b,c)[fsgrids::bfield::PERBX];
                                    nCells++;
                                }
                            }
                            }
                        }
                    }
                    if (nCells == 0) {
                        #ifndef __CUDA_ARCH__
                        cerr << __FILE__ << ":" << __LINE__ << ": ERROR: this should not have fallen through." << endl;
                        #endif
                        return 0.0;
                    }
                    return retval / nCells;
                }
                case 1:
                if (  (technicalGrid.get(i,j-1,k)->SOLVE & compute::BY) == compute::BY
                    && (technicalGrid.get(i,j+1,k)->SOLVE & compute::BY) == compute::BY
                ) {
                    return 0.5 * (bGrid.get(i,j-1,k)[fsgrids::bfield::PERBY] + bGrid.get(i,j+1,k)[fsgrids::bfield::PERBY]);
                } else if ((technicalGrid.get(i,j-1,k)->SOLVE & compute::BY) == compute::BY) {
                    return bGrid.get(i,j-1,k)[fsgrids::bfield::PERBY];
                } else if ((technicalGrid.get(i,j+1,k)->SOLVE & compute::BY) == compute::BY) {
                    return bGrid.get(i,j+1,k)[fsgrids::bfield::PERBY];
                } else {
                    Real retval = 0.0;
                    uint nCells = 0;
                    if ((technicalGrid.get(i-1,j,k)->SOLVE & compute::BY) == compute::BY) {
                        retval += bGrid.get(i-1,j,k)[fsgrids::bfield::PERBY];
                        nCells++;
                    }
                    if ((technicalGrid.get(i+1,j,k)->SOLVE & compute::BY) == compute::BY) {
                        retval += bGrid.get(i+1,j,k)[fsgrids::bfield::PERBY];
                        nCells++;
                    }
                    if ((technicalGrid.get(i,j,k-1)->SOLVE & compute::BY) == compute::BY) {
                        retval += bGrid.get(i,j,k-1)[fsgrids::bfield::PERBY];
                        nCells++;
                    }
                    if ((technicalGrid.get(i,j,k+1)->SOLVE & compute::BY) == compute::BY) {
                        retval += bGrid.get(i,j,k+1)[fsgrids::bfield::PERBY];
                        nCells++;
                    }
                    if (nCells == 0) {
                        for (int a=i-1; a<i+2; a++) {
                            for (int b=j-1; b<j+2; b++) {
                            for (int c=k-1; c<k+2; c++) {
                                if ((technicalGrid.get(a,b,c)->SOLVE & compute::BY) == compute::BY) {
                                    retval += bGrid.get(a,b,c)[fsgrids::bfield::PERBY];
                                    nCells++;
                                }
                            }
                            }
                        }
                    }
                    if (nCells == 0) {
                        #ifndef __CUDA_ARCH__
                        cerr << __FILE__ << ":" << __LINE__ << ": ERROR: this should not have fallen through." << endl;
                        #endif
                        return 0.0;
                    }
                    return retval / nCells;
                }
                case 2:
                if (  (technicalGrid.get(i,j,k-1)->SOLVE & compute::BZ) == compute::BZ
                    && (technicalGrid.get(i,j,k+1)->SOLVE & compute::BZ) == compute::BZ
                ) {
                    return 0.5 * (bGrid.get(i,j,k-1)[fsgrids::bfield::PERBZ] + bGrid.get(i,j,k+1)[fsgrids::bfield::PERBZ]);
                } else if ((technicalGrid.get(i,j,k-1)->SOLVE & compute::BZ) == compute::BZ) {
                    return bGrid.get(i,j,k-1)[fsgrids::bfield::PERBZ];
                } else if ((technicalGrid.get(i,j,k+1)->SOLVE & compute::BZ) == compute::BZ) {
                    return bGrid.get(i,j,k+1)[fsgrids::bfield::PERBZ];
                } else {
                    Real retval = 0.0;
                    uint nCells = 0;
                    if ((technicalGrid.get(i-1,j,k)->SOLVE & compute::BZ) == compute::BZ) {
                        retval += bGrid.get(i-1,j,k)[fsgrids::bfield::PERBZ];
                        nCells++;
                    }
                    if ((technicalGrid.get(i+1,j,k)->SOLVE & compute::BZ) == compute::BZ) {
                        retval += bGrid.get(i+1,j,k)[fsgrids::bfield::PERBZ];
                        nCells++;
                    }
                    if ((technicalGrid.get(i,j-1,k)->SOLVE & compute::BZ) == compute::BZ) {
                        retval += bGrid.get(i,j-1,k)[fsgrids::bfield::PERBZ];
                        nCells++;
                    }
                    if ((technicalGrid.get(i,j+1,k)->SOLVE & compute::BZ) == compute::BZ) {
                        retval += bGrid.get(i,j+1,k)[fsgrids::bfield::PERBZ];
                        nCells++;
                    }
                    if (nCells == 0) {
                        for (int a=i-1; a<i+2; a++) {
                            for (int b=j-1; b<j+2; b++) {
                            for (int c=k-1; c<k+2; c++) {
                                if ((technicalGrid.get(a,b,c)->SOLVE & compute::BZ) == compute::BZ) {
                                    retval += bGrid.get(a,b,c)[fsgrids::bfield::PERBZ];
                                    nCells++;
                                }
                            }
                            }
                        }
                    }
                    if (nCells == 0) {
                        #ifndef __CUDA_ARCH__ 
                        cerr << __FILE__ << ":" << __LINE__ << ": ERROR: this should not have fallen through." << endl;
                        #endif 
                        return 0.0;
                    }
                    return retval / nCells;
                }
                default:
                #ifndef __CUDA_ARCH__
                cerr << "ERROR: ionosphere boundary tried to copy nonsensical magnetic field component " << component << endl;
                #endif 
                return 0.0;
            }
        } else { // L2 cells
            Real retval = 0.0;
            uint nCells = 0;
            for (int a=i-1; a<i+2; a++) {
                for (int b=j-1; b<j+2; b++) {
                for (int c=k-1; c<k+2; c++) {
                    if (technicalGrid.get(a,b,c)->sysBoundaryLayer == 1) {
                        retval += bGrid.get(a,b,c)[fsgrids::bfield::PERBX + component];
                        nCells++;
                    }
                }
                }
            }
            if (nCells == 0) {
                #ifndef __CUDA_ARCH__
                cerr << __FILE__ << ":" << __LINE__ << ": ERROR: this should not have fallen through." << endl;
                #endif 
                return 0.0;
            }
            return retval / nCells;
        }
    }

    /*! We want here to
        *
        * -- Retain only the boundary-normal projection of perturbed face B
        */
    ARCH_HOSTDEV void fieldSolverBoundaryCondMagneticFieldProjection(
        const arch::buf<FsGrid<Real, fsgrids::bfield::N_BFIELD, FS_STENCIL_WIDTH>> & bGrid,
        const arch::buf<FsGrid< fsgrids::technical, 1, FS_STENCIL_WIDTH>> & technicalGrid,
        cint i,
        cint j,
        cint k
    ) {
        // Projection of B-field to normal direction
        Real BdotN = 0;
        int32_t normalDirection[3]; 
        fieldSolverGetNormalDirection(technicalGrid, i, j, k, normalDirection);
        for(uint component=0; component<3; component++) {
            BdotN += bGrid.get(i,j,k)[fsgrids::bfield::PERBX+component] * normalDirection[component];
        }
        // Apply to any components that were not solved
        if ((technicalGrid.get(i,j,k)->sysBoundaryLayer == 2) ||
            ((technicalGrid.get(i,j,k)->sysBoundaryLayer == 1) && ((technicalGrid.get(i,j,k)->SOLVE & compute::BX) != compute::BX))
            ) {
            bGrid.get(i,j,k)[fsgrids::bfield::PERBX] = BdotN*normalDirection[0];
        }
        if ((technicalGrid.get(i,j,k)->sysBoundaryLayer == 2) ||
            ((technicalGrid.get(i,j,k)->sysBoundaryLayer == 1) && ((technicalGrid.get(i,j,k)->SOLVE & compute::BY) != compute::BY))
            ) {
            bGrid.get(i,j,k)[fsgrids::bfield::PERBY] = BdotN*normalDirection[1];
        }
        if ((technicalGrid.get(i,j,k)->sysBoundaryLayer == 2) ||
            ((technicalGrid.get(i,j,k)->sysBoundaryLayer == 1) && ((technicalGrid.get(i,j,k)->SOLVE & compute::BZ) != compute::BZ))
            ) {
            bGrid.get(i,j,k)[fsgrids::bfield::PERBZ] = BdotN*normalDirection[2];
        }
    }

    ARCH_HOSTDEV void fieldSolverGetNormalDirection(
        const arch::buf<FsGrid< fsgrids::technical, 1, FS_STENCIL_WIDTH>> & technicalGrid,
        cint i,
        cint j,
        cint k,
        int32_t (&normalDirection)[3]
    ) {
        phiprof::start("Ionosphere::fieldSolverGetNormalDirection");
        
        static creal DIAG2 = 0.7071067811865475; // 1.0 / sqrt(2.0);
        static creal DIAG3 = 0.5773502691896257; // 1.0 / sqrt(3.0);
        
        creal dx = technicalGrid.grid()->DX;
        creal dy = technicalGrid.grid()->DY;
        creal dz = technicalGrid.grid()->DZ;
        int globalIndices[3];
        technicalGrid.grid()->getGlobalIndices(i,j,k, globalIndices);
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
                #ifndef __CUDA_ARCH__
                std::cerr << __FILE__ << ":" << __LINE__ << ":" << "What do you expect to do with a single-cell simulation of ionosphere boundary type? Stop kidding." << std::endl;
                #endif
                assert(0);
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
                    #ifndef __CUDA_ARCH__ 
                    std::cerr << __FILE__ << ":" << __LINE__ << ":" << "ionosphere.geometry has to be 0, 1 or 2 with this grid shape." << std::endl;
                    #endif
                    assert(0); 
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
                    #ifndef __CUDA_ARCH__ 
                    std::cerr << __FILE__ << ":" << __LINE__ << ":" << "ionosphere.geometry has to be 0, 1, 2 or 3 with this grid shape." << std::endl;
                    #endif 
                    assert(0);
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
                #ifndef __CUDA_ARCH__
                std::cerr << __FILE__ << ":" << __LINE__ << ":" << "ionosphere.geometry has to be 0, 1 or 2 with this grid shape." << std::endl;
                #endif 
                assert(0); 
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
                #ifndef __CUDA_ARCH__
                std::cerr << __FILE__ << ":" << __LINE__ << ":" << "ionosphere.geometry has to be 0, 1, 2 or 3 with this grid shape." << std::endl;
                #endif 
                assert(0); 
            }
            // end of 3D
        }
        phiprof::stop("Ionosphere::fieldSolverGetNormalDirection");
    }

    ARCH_HOSTDEV void fieldSolverBoundaryCondElectricField(
      const arch::buf<FsGrid<Real, fsgrids::efield::N_EFIELD, FS_STENCIL_WIDTH>> & EGrid,
      cint i,
      cint j,
      cint k,
      cuint component
   ) {
      EGrid.get(i,j,k)[fsgrids::efield::EX+component] = 0.0;
   }
   
   ARCH_HOSTDEV void fieldSolverBoundaryCondHallElectricField(
      const arch::buf<FsGrid<Real, fsgrids::ehall::N_EHALL, FS_STENCIL_WIDTH>> & EHallGrid,
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
            #ifndef __CUDA_ARCH__
            cerr << __FILE__ << ":" << __LINE__ << ":" << " Invalid component" << endl;
            #endif
            return;
      }
   }
   
   ARCH_HOSTDEV void fieldSolverBoundaryCondGradPeElectricField(
      const arch::buf<FsGrid<Real, fsgrids::egradpe::N_EGRADPE, FS_STENCIL_WIDTH>> & EGradPeGrid,
      cint i,
      cint j,
      cint k,
      cuint component
   ) {
      EGradPeGrid.get(i,j,k)[fsgrids::egradpe::EXGRADPE+component] = 0.0;
   }
   
   ARCH_HOSTDEV void fieldSolverBoundaryCondDerivatives(
      const arch::buf<FsGrid<Real, fsgrids::dperb::N_DPERB, FS_STENCIL_WIDTH>> & dPerBGrid,
      const arch::buf<FsGrid<Real, fsgrids::dmoments::N_DMOMENTS, FS_STENCIL_WIDTH>> & dMomentsGrid,
      cint i,
      cint j,
      cint k,
      cuint& RKCase,
      cuint& component
   ) {
      SysBoundaryCondition::setCellDerivativesToZero(dPerBGrid, dMomentsGrid, i, j, k, component);
      return;
   }
   
   ARCH_HOSTDEV void fieldSolverBoundaryCondBVOLDerivatives(
      const arch::buf<FsGrid<Real, fsgrids::volfields::N_VOL, FS_STENCIL_WIDTH>> & volGrid,
      cint i,
      cint j,
      cint k,
      cuint& component
   ) {
      // FIXME This should be OK as the BVOL derivatives are only used for Lorentz force JXB, which is not applied on the ionosphere cells.
      SysBoundaryCondition::setCellBVOLDerivativesToZero(volGrid, i, j, k, component);
   } 

    private:
    Real center[3]; /*!< Coordinates of the centre of the ionosphere. */
    Real radius; /*!< Radius of the ionosphere. */
    uint geometry; /*!< Geometry of the ionosphere, 0: inf-norm (diamond), 1: 1-norm (square), 2: 2-norm (circle, DEFAULT), 3: polar-plane cylinder with line dipole. */


    }; 
}