#ifndef SETBYUSER_FIELD_BOUNDARY_H
#define SETBYUSER_FIELD_BOUNDARY_H

#include "../definitions.h"
#include "sysboundarycondition.h"
#include "../common.h"
#include "../parameters.h"
#include <vector>
#include <array>
// #include "sysboundaryField.h"

extern ARCH_MANAGED GridParameters meshParams;

// define class
namespace SBC {
    class SetByUserFieldBoundary {

    public:
        SetByUserFieldBoundary(bool isPeriodic[3], Real templateB[6][3]): isPeriodic{isPeriodic[0], isPeriodic[1], isPeriodic[2]} {
            for (uint i = 0; i < 6; i++) {
                for (uint j = 0; j < 3; j++) {
                    this->templateB[i][j] = templateB[i][j];
                }
            }
        }

        ARCH_HOSTDEV Real fieldSolverBoundaryCondMagneticField(
            const arch::buf<FsGrid<Real, fsgrids::bfield::N_BFIELD, FS_STENCIL_WIDTH>> & bGrid,
            const arch::buf<FsGrid< fsgrids::technical, 1, FS_STENCIL_WIDTH>> & technicalGrid,
            cint i,
            cint j,
            cint k,
            creal& dt,
            cuint& component
        ) {
            Real result = 0.0;

            creal dx = meshParams.dx_ini;
            creal dy = meshParams.dy_ini;
            creal dz = meshParams.dz_ini;
            int32_t globalIndices[3];
            technicalGrid.grid()->getGlobalIndices(i,j,k, globalIndices);
           
            creal x = (globalIndices[0] + 0.5) * technicalGrid.grid()->DX + meshParams.xmin;
            creal y = (globalIndices[1] + 0.5) * technicalGrid.grid()->DY + meshParams.ymin;
            creal z = (globalIndices[2] + 0.5) * technicalGrid.grid()->DZ + meshParams.zmin;

            bool isThisCellOnAFace[6];
            determineFace(&isThisCellOnAFace[0], x, y, z, dx, dy, dz, isPeriodic, true);

            for (uint i=0; i<6; i++) {
                if (isThisCellOnAFace[i]) {
                result = templateB[i][component];
                break; // This effectively sets the precedence of faces through the order of faces.
                }
            }
            // return x;
            return result;
        }
        void fieldSolverBoundaryCondMagneticFieldProjection(
            FsGrid<Real, fsgrids::bfield::N_BFIELD, FS_STENCIL_WIDTH> & bGrid,
            FsGrid< fsgrids::technical, 1, FS_STENCIL_WIDTH> & technicalGrid,
            cint i,
            cint j,
            cint k
        );
        void fieldSolverBoundaryCondElectricField(
            FsGrid<Real, fsgrids::efield::N_EFIELD, FS_STENCIL_WIDTH> & EGrid,
            cint i,
            cint j,
            cint k,
            cuint component
        );
        void fieldSolverBoundaryCondHallElectricField(
            FsGrid<Real, fsgrids::ehall::N_EHALL, FS_STENCIL_WIDTH> & EHallGrid,
            cint i,
            cint j,
            cint k,
            cuint component
        );
        void fieldSolverBoundaryCondGradPeElectricField(
            FsGrid<Real, fsgrids::egradpe::N_EGRADPE, FS_STENCIL_WIDTH> & EGradPeGrid,
            cint i,
            cint j,
            cint k,
            cuint component
        );
        void fieldSolverBoundaryCondDerivatives(
            FsGrid<Real, fsgrids::dperb::N_DPERB, FS_STENCIL_WIDTH> & dPerBGrid,
            FsGrid<Real, fsgrids::dmoments::N_DMOMENTS, FS_STENCIL_WIDTH> & dMomentsGrid,
            cint i,
            cint j,
            cint k,
            cuint& RKCase,
            cuint& component
        );
        void fieldSolverBoundaryCondBVOLDerivatives(
            FsGrid<Real, fsgrids::volfields::N_VOL, FS_STENCIL_WIDTH> & volGrid,
            cint i,
            cint j,
            cint k,
            cuint& component
        ); 

        ARCH_HOSTDEV void determineFace(
            bool* isThisCellOnAFace,
            creal x, creal y, creal z,
            creal dx, creal dy, creal dz,
            bool (&isPeriodic)[3],
            const bool excludeSlicesAndPeriodicDimensions //=false (default)
        ) {
            
            for(uint i=0; i<6; i++) {
                isThisCellOnAFace[i] = false;
            }
            if(x > meshParams.xmax - 2.0*dx) {
                isThisCellOnAFace[0] = true;
            }
            if(x < meshParams.xmin + 2.0*dx) {
                isThisCellOnAFace[1] = true;
            }
            if(y > meshParams.ymax - 2.0*dy) {
                isThisCellOnAFace[2] = true;
            }
            if(y < meshParams.ymin + 2.0*dy) {
                isThisCellOnAFace[3] = true;
            }
            if(z > meshParams.zmax - 2.0*dz) {
                isThisCellOnAFace[4] = true;
            }
            if(z < meshParams.zmin + 2.0*dz) {
                isThisCellOnAFace[5] = true;
            }
            if(excludeSlicesAndPeriodicDimensions == true) {
                if(meshParams.xcells_ini == 1 || isPeriodic[0]) {
                    isThisCellOnAFace[0] = false;
                    isThisCellOnAFace[1] = false;
                }
                if(meshParams.ycells_ini == 1 || isPeriodic[1]) {
                    isThisCellOnAFace[2] = false;
                    isThisCellOnAFace[3] = false;
                }
                if(meshParams.zcells_ini == 1 || isPeriodic[2]) {
                    isThisCellOnAFace[4] = false;
                    isThisCellOnAFace[5] = false;
                }
            }
        } 

    private:
        bool isPeriodic[3];
        Real templateB[6][3];
    };
}

#endif