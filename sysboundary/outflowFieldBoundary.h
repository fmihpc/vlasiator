#include "../common.h"

using namespace std;

namespace SBC {
    class OutflowFieldBoundary {

    public:
    ARCH_HOSTDEV Real fieldSolverBoundaryCondMagneticField(
        const arch::buf<FsGrid<Real, fsgrids::bfield::N_BFIELD, FS_STENCIL_WIDTH>> & bGrid,
        const arch::buf<FsGrid< fsgrids::technical, 1, FS_STENCIL_WIDTH>> & technicalGrid,
        cint i,
        cint j,
        cint k,
        creal& dt,
        cuint& component
    ) {
        switch(component) {
        case 0:
            return fieldBoundaryCopyFromSolvingNbrMagneticField(bGrid, technicalGrid, i, j, k, component, compute::BX);
        case 1:
            return fieldBoundaryCopyFromSolvingNbrMagneticField(bGrid, technicalGrid, i, j, k, component, compute::BY);
        case 2:
            return fieldBoundaryCopyFromSolvingNbrMagneticField(bGrid, technicalGrid, i, j, k, component, compute::BZ);
        default:
            return 0.0;
        }
    }
      ARCH_HOSTDEV void fieldSolverBoundaryCondMagneticFieldProjection(
         const arch::buf<FsGrid<Real, fsgrids::bfield::N_BFIELD, FS_STENCIL_WIDTH>> & bGrid,
         const arch::buf<FsGrid< fsgrids::technical, 1, FS_STENCIL_WIDTH>> & technicalGrid,
         cint i,
         cint j,
         cint k
      ) {}

        ARCH_HOSTDEV Real fieldBoundaryCopyFromSolvingNbrMagneticField(
            const arch::buf<FsGrid< Real, fsgrids::bfield::N_BFIELD, FS_STENCIL_WIDTH>> & bGrid,
            const arch::buf<FsGrid< fsgrids::technical, 1, FS_STENCIL_WIDTH>> & technicalGrid,
            cint i,
            cint j,
            cint k,
            cuint component,
            cuint mask
        ) {

            int distance = numeric_limits<int>::max();
            int closestCell[3];
            bool found = false;
 
            // vector< array<int,3> > closestCells;

            for (int kk=-2; kk<3; kk++) {
                for (int jj=-2; jj<3; jj++) {
                    for (int ii=-2; ii<3 ; ii++) {
                    if( technicalGrid.get(i+ii,j+jj,k+kk) // skip invalid cells returning NULL
                        && (technicalGrid.get(i+ii,j+jj,k+kk)->SOLVE & mask) == mask // Did that guy solve this component?
                        && technicalGrid.get(i+ii,j+jj,k+kk)->sysBoundaryFlag != sysboundarytype::DO_NOT_COMPUTE // Do not copy from there
                    ) { 
                        distance = min(distance, ii*ii + jj*jj + kk*kk);
                    }
                    }
                }
            }

            for (int kk=-2; kk<3; kk++) {
                for (int jj=-2; jj<3; jj++) {
                    for (int ii=-2; ii<3 ; ii++) {
                    if( technicalGrid.get(i+ii,j+jj,k+kk) // skip invalid cells returning NULL
                        && (technicalGrid.get(i+ii,j+jj,k+kk)->SOLVE & mask) == mask // Did that guy solve this component?
                        && technicalGrid.get(i+ii,j+jj,k+kk)->sysBoundaryFlag != sysboundarytype::DO_NOT_COMPUTE // Do not copy from there
                    ) { 
                        int d = ii*ii + jj*jj + kk*kk;
                        if( d == distance ) {
                            closestCell[0] = i+ii;
                            closestCell[1] = j+jj;
                            closestCell[2] = k+kk;
                            found = true;
                            break;
                        }
                    }
                    }
                    if (found) break;
                }
                if (found) break; 
            }

            if(!found) {
                #ifndef __CUDA_ARCH__ 
                cerr << __FILE__ << ":" << __LINE__ << ": No closest cell found!" << endl;
                abort();
                #endif
            }
            assert(found);
            return bGrid.get(closestCell[0], closestCell[1], closestCell[2])[fsgrids::bfield::PERBX+component];
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
                    cerr << __FILE__ << ":" << __LINE__ << ":" << " Invalid component" << endl;
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
        }
        
        ARCH_HOSTDEV void fieldSolverBoundaryCondBVOLDerivatives(
            const arch::buf<FsGrid<Real, fsgrids::volfields::N_VOL, FS_STENCIL_WIDTH>> & volGrid,
            cint i,
            cint j,
            cint k,
            cuint& component
        ) {
            SysBoundaryCondition::setCellBVOLDerivativesToZero(volGrid, i, j, k, component);
        }
    }; 
}