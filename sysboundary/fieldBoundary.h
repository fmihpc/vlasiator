#ifndef FIELD_BOUNDARY_H
#define FIELD_BOUNDARY_H

#include "../arch/arch_device_api.h"
#include "../fieldsolver/fs_common.h"
namespace SBC {
    class FieldBoundary {
        public:
        ARCH_HOSTDEV void fieldSolverGetNormalDirection(
            const arch::buf<FsGrid< fsgrids::technical, 1, FS_STENCIL_WIDTH>> & technicalGrid,
            cint i,
            cint j,
            cint k,
            Real (&normalDirection)[3]
        ) {
            phiprof::start("Ionosphere::fieldSolverGetNormalDirection");
            
            static creal DIAG2 = 0.7071067811865475; // 1.0 / sqrt(2.0);
            static creal DIAG3 = 0.5773502691896257; // 1.0 / sqrt(3.0);
            
            creal dx = technicalGrid.grid()->DX;
            creal dy = technicalGrid.grid()->DY;
            creal dz = technicalGrid.grid()->DZ;
            int globalIndices[3];
            technicalGrid.grid()->getGlobalIndices(i,j,k, globalIndices);
            creal x = FSParams.xmin + (convert<Real>(globalIndices[0])+0.5)*dx;
            creal y = FSParams.ymin + (convert<Real>(globalIndices[1])+0.5)*dy;
            creal z = FSParams.zmin + (convert<Real>(globalIndices[2])+0.5)*dz;
            creal xsign = divideIfNonZero(x, fabs(x));
            creal ysign = divideIfNonZero(y, fabs(y));
            creal zsign = divideIfNonZero(z, fabs(z));
            
            Real length = 0.0;
            
            if (FSParams.xcells == 1) {
                if (FSParams.ycells == 1) {
                    if (FSParams.zcells == 1) {
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
                } else if (FSParams.zcells == 1) {
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
            } else if (FSParams.ycells == 1) {
                if (FSParams.zcells == 1) {
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
            } else if (FSParams.zcells == 1) {
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

        protected:
            Real center[3]; /*!< Coordinates of the centre of the (iono)sphere. */
            Real radius; /*!< Radius of the (iono)sphere. */
            uint geometry; /*!< Geometry of the (iono)sphere, 0: inf-norm (diamond), 1: 1-norm (square), 2: 2-norm (circle, DEFAULT), 3: polar-plane cylinder with line dipole. */
    };
}
#endif