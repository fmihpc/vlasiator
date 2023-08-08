#include "setbyuserFieldBoundary.h"

using namespace std;

namespace SBC {
    void SetByUserFieldBoundary::fieldSolverBoundaryCondMagneticFieldProjection(
        FsGrid<Real, fsgrids::bfield::N_BFIELD, FS_STENCIL_WIDTH> & bGrid,
        FsGrid< fsgrids::technical, 1, FS_STENCIL_WIDTH> & technicalGrid,
        cint i,
        cint j,
        cint k
    ) {
    }

    void SetByUserFieldBoundary::fieldSolverBoundaryCondElectricField(
        FsGrid<Real, fsgrids::efield::N_EFIELD, FS_STENCIL_WIDTH> & EGrid,
        cint i,
        cint j,
        cint k,
        cuint component
    ) {
        EGrid.get(i,j,k)[fsgrids::efield::EX+component] = 0.0;
    }

    void SetByUserFieldBoundary::fieldSolverBoundaryCondHallElectricField(
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

    void SetByUserFieldBoundary::fieldSolverBoundaryCondGradPeElectricField(
        FsGrid<Real, fsgrids::egradpe::N_EGRADPE, FS_STENCIL_WIDTH> & EGradPeGrid,
        cint i,
        cint j,
        cint k,
        cuint component
    ) {
            EGradPeGrid.get(i,j,k)[fsgrids::egradpe::EXGRADPE+component] = 0.0;
    }

    void SetByUserFieldBoundary::fieldSolverBoundaryCondDerivatives(
        FsGrid<Real, fsgrids::dperb::N_DPERB, FS_STENCIL_WIDTH> & dPerBGrid,
        FsGrid<Real, fsgrids::dmoments::N_DMOMENTS, FS_STENCIL_WIDTH> & dMomentsGrid,
        cint i,
        cint j,
        cint k,
        cuint& RKCase,
        cuint& component
    ) {
        SysBoundaryCondition::setCellDerivativesToZero(dPerBGrid, dMomentsGrid, i, j, k, component);
    }

    void SetByUserFieldBoundary::fieldSolverBoundaryCondBVOLDerivatives(
        FsGrid<Real, fsgrids::volfields::N_VOL, FS_STENCIL_WIDTH> & volGrid,
        cint i,
        cint j,
        cint k,
        cuint& component
    ) {
        SysBoundaryCondition::setCellBVOLDerivativesToZero(volGrid, i, j, k, component);
    }
}