#include "../common.h"
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

      std::array<Real, 3> fieldSolverGetNormalDirection(
         FsGrid< fsgrids::technical, 1, FS_STENCIL_WIDTH> & technicalGrid,
         cint i,
         cint j,
         cint k
      ); 

    private:
    Real center[3]; /*!< Coordinates of the centre of the ionosphere. */
    Real radius; /*!< Radius of the ionosphere. */
    uint geometry; /*!< Geometry of the ionosphere, 0: inf-norm (diamond), 1: 1-norm (square), 2: 2-norm (circle, DEFAULT), 3: polar-plane cylinder with line dipole. */


    }; 
}