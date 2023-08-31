#ifndef ARCH_SYSBOUNDARY_H
#define ARCH_SYSBOUNDARY_H

#include "../sysboundary/sysboundary.h"
#include "../sysboundary/donotcompute.h"
#include "../sysboundary/ionosphere.h"
#include "../sysboundary/outflow.h"
#include "../sysboundary/setmaxwellian.h"
#include "../sysboundary/conductingsphere.h"

namespace arch {

template <>
class buf<SysBoundary> {
    private:  
        SysBoundary *ptr; 

        // list of possible sysboundaries
        SBC::SetMaxwellian* setmaxwellian;
        SBC::Ionosphere* ionosphere;
        SBC::Outflow* outflow;
        SBC::DoNotCompute* doNotCompute;
        SBC::Conductingsphere* conductingsphere;

        uint is_copy = 0;

    public:   

        void syncDeviceData(void){}

        void syncHostData(void){}

        class Proxy {
            public:
                Proxy(int _sysBoundaryFlag, const buf<SysBoundary>* _ptr) : sysBoundaryFlag(_sysBoundaryFlag), bufPtr(_ptr) {}

                Real fieldSolverBoundaryCondMagneticField(
                    const buf<FsGrid<Real, fsgrids::bfield::N_BFIELD, FS_STENCIL_WIDTH>> & bGrid,
                    const buf<FsGrid< fsgrids::technical, 1, FS_STENCIL_WIDTH>> & technicalGrid,
                    cint i,
                    cint j,
                    cint k,
                    creal& dt,
                    cuint& component
                ) {
                    if (sysBoundaryFlag == sysboundarytype::SET_MAXWELLIAN) {
                        return bufPtr->setmaxwellian->fieldSolverBoundaryCondMagneticField(bGrid, technicalGrid, i, j, k, dt,component);
                    } else if (sysBoundaryFlag == sysboundarytype::IONOSPHERE) {
                        return bufPtr->ionosphere->fieldSolverBoundaryCondMagneticField(bGrid, technicalGrid, i, j, k, dt,
                                                                                    component);
                    } else if (sysBoundaryFlag == sysboundarytype::OUTFLOW) {
                        return bufPtr->outflow->fieldSolverBoundaryCondMagneticField(bGrid, technicalGrid, i, j, k, dt, component);
                    } else if (sysBoundaryFlag == sysboundarytype::CONDUCTINGSPHERE) {
                        return bufPtr->conductingsphere->fieldSolverBoundaryCondMagneticField(bGrid, technicalGrid, i, j, k, dt, component);
                    } else {
                        std::cerr << "ERROR: sysboundarytype not found" << std::endl;
                        exit(1);
                    }
                }

                void fieldSolverBoundaryCondMagneticFieldProjection(
                    const arch::buf<FsGrid<Real, fsgrids::bfield::N_BFIELD, FS_STENCIL_WIDTH>> & bGrid,
                    const arch::buf<FsGrid< fsgrids::technical, 1, FS_STENCIL_WIDTH>> & technicalGrid,
                    cint i,
                    cint j,
                    cint k
                ) {
                    if (sysBoundaryFlag == sysboundarytype::SET_MAXWELLIAN) {
                        bufPtr->setmaxwellian->fieldSolverBoundaryCondMagneticFieldProjection(bGrid, technicalGrid, i, j, k);
                    } else if (sysBoundaryFlag == sysboundarytype::IONOSPHERE) {
                        bufPtr->ionosphere->fieldSolverBoundaryCondMagneticFieldProjection(bGrid, technicalGrid, i, j, k);
                    } else if (sysBoundaryFlag == sysboundarytype::OUTFLOW) {
                        bufPtr->outflow->fieldSolverBoundaryCondMagneticFieldProjection(bGrid, technicalGrid, i, j, k);
                    } else if (sysBoundaryFlag == sysboundarytype::CONDUCTINGSPHERE) {
                        bufPtr->conductingsphere->fieldSolverBoundaryCondMagneticFieldProjection(bGrid, technicalGrid, i, j, k);
                    } else {
                        std::cerr << "ERROR: sysboundarytype not found" << std::endl;
                        exit(1);
                    }
                }

                void fieldSolverBoundaryCondDerivatives(
                    const arch::buf<FsGrid<Real, fsgrids::dperb::N_DPERB, FS_STENCIL_WIDTH>> & dPerBGrid,
                    const arch::buf<FsGrid<Real, fsgrids::dmoments::N_DMOMENTS, FS_STENCIL_WIDTH>> & dMomentsGrid,
                    cint i,
                    cint j,
                    cint k,
                    cuint& RKCase,
                    cuint& component
                ) {
                    if (sysBoundaryFlag == sysboundarytype::SET_MAXWELLIAN) {
                        bufPtr->setmaxwellian->fieldSolverBoundaryCondDerivatives(dPerBGrid, dMomentsGrid, i, j, k, RKCase, component);
                    } else if (sysBoundaryFlag == sysboundarytype::IONOSPHERE) {
                        bufPtr->ionosphere->fieldSolverBoundaryCondDerivatives(dPerBGrid, dMomentsGrid, i, j, k, RKCase, component);
                    } else if (sysBoundaryFlag == sysboundarytype::OUTFLOW) {
                        bufPtr->outflow->fieldSolverBoundaryCondDerivatives(dPerBGrid, dMomentsGrid, i, j, k, RKCase, component);
                    } else if (sysBoundaryFlag == sysboundarytype::CONDUCTINGSPHERE) {
                        bufPtr->conductingsphere->fieldSolverBoundaryCondDerivatives(dPerBGrid, dMomentsGrid, i, j, k, RKCase, component);
                    } else {
                        std::cerr << "ERROR: sysboundarytype not found" << std::endl;
                        exit(1);
                    }
                }

                void fieldSolverBoundaryCondGradPeElectricField(
                    const arch::buf<FsGrid<Real, fsgrids::egradpe::N_EGRADPE, FS_STENCIL_WIDTH>> & EGradPeGrid,
                    cint i,
                    cint j,
                    cint k,
                    cuint component
                ) {
                    if (sysBoundaryFlag == sysboundarytype::SET_MAXWELLIAN) {
                        bufPtr->setmaxwellian->fieldSolverBoundaryCondGradPeElectricField(EGradPeGrid, i, j, k, component);
                    } else if (sysBoundaryFlag == sysboundarytype::IONOSPHERE) {
                        bufPtr->ionosphere->fieldSolverBoundaryCondGradPeElectricField(EGradPeGrid, i, j, k, component);
                    } else if (sysBoundaryFlag == sysboundarytype::OUTFLOW) {
                        bufPtr->outflow->fieldSolverBoundaryCondGradPeElectricField(EGradPeGrid, i, j, k, component);
                    } else if (sysBoundaryFlag == sysboundarytype::CONDUCTINGSPHERE) {
                        bufPtr->conductingsphere->fieldSolverBoundaryCondGradPeElectricField(EGradPeGrid, i, j, k, component);
                    } else {
                        std::cerr << "ERROR: sysboundarytype not found" << std::endl;
                        exit(1);
                    }
                }

                void fieldSolverBoundaryCondHallElectricField(
                    const arch::buf<FsGrid<Real, fsgrids::ehall::N_EHALL, FS_STENCIL_WIDTH>> & EHallGrid,
                    cint i,
                    cint j,
                    cint k,
                    cuint component
                ) { 
                    if (sysBoundaryFlag == sysboundarytype::SET_MAXWELLIAN) {
                        bufPtr->setmaxwellian->fieldSolverBoundaryCondHallElectricField(EHallGrid, i, j, k, component);
                    } else if (sysBoundaryFlag == sysboundarytype::IONOSPHERE) {
                        bufPtr->ionosphere->fieldSolverBoundaryCondHallElectricField(EHallGrid, i, j, k, component);
                    } else if (sysBoundaryFlag == sysboundarytype::OUTFLOW) {
                        bufPtr->outflow->fieldSolverBoundaryCondHallElectricField(EHallGrid, i, j, k, component);
                    } else if (sysBoundaryFlag == sysboundarytype::CONDUCTINGSPHERE) {
                        bufPtr->conductingsphere->fieldSolverBoundaryCondHallElectricField(EHallGrid, i, j, k, component);
                    } else {
                        std::cerr << "ERROR: sysboundarytype not found" << std::endl;
                        exit(1);
                    }
                }


                void fieldSolverBoundaryCondElectricField(
                    const arch::buf<FsGrid<Real, fsgrids::efield::N_EFIELD, FS_STENCIL_WIDTH>> & EGrid,
                    cint i,
                    cint j,
                    cint k,
                    cuint component
                ) {
                    if (sysBoundaryFlag == sysboundarytype::SET_MAXWELLIAN) {
                        bufPtr->setmaxwellian->fieldSolverBoundaryCondElectricField(EGrid, i, j, k, component);
                    } else if (sysBoundaryFlag == sysboundarytype::IONOSPHERE) {
                        bufPtr->ionosphere->fieldSolverBoundaryCondElectricField(EGrid, i, j, k, component);
                    } else if (sysBoundaryFlag == sysboundarytype::OUTFLOW) {
                        bufPtr->outflow->fieldSolverBoundaryCondElectricField(EGrid, i, j, k, component);
                    } else if (sysBoundaryFlag == sysboundarytype::CONDUCTINGSPHERE) {
                        bufPtr->conductingsphere->fieldSolverBoundaryCondElectricField(EGrid, i, j, k, component);
                    } else {
                        std::cerr << "ERROR: sysboundarytype not found" << std::endl;
                        exit(1);
                    }
                }

                void fieldSolverBoundaryCondBVOLDerivatives(
                    const arch::buf<FsGrid<Real, fsgrids::volfields::N_VOL, FS_STENCIL_WIDTH>> & volGrid,
                    cint i,
                    cint j,
                    cint k,
                    cuint& component
                ) {
                    if (sysBoundaryFlag == sysboundarytype::SET_MAXWELLIAN) {
                        bufPtr->setmaxwellian->fieldSolverBoundaryCondBVOLDerivatives(volGrid, i, j, k, component);
                    } else if (sysBoundaryFlag == sysboundarytype::IONOSPHERE) {
                        bufPtr->ionosphere->fieldSolverBoundaryCondBVOLDerivatives(volGrid, i, j, k, component);
                    } else if (sysBoundaryFlag == sysboundarytype::OUTFLOW) {
                        bufPtr->outflow->fieldSolverBoundaryCondBVOLDerivatives(volGrid, i, j, k, component);
                    } else if (sysBoundaryFlag == sysboundarytype::CONDUCTINGSPHERE) {
                        bufPtr->conductingsphere->fieldSolverBoundaryCondBVOLDerivatives(volGrid, i, j, k, component);
                    } else {
                        std::cerr << "ERROR: sysboundarytype not found" << std::endl;
                        exit(1);
                    }
                }



            private:
                int sysBoundaryFlag;
                const buf<SysBoundary>* bufPtr;

        };

        buf(SysBoundary* const _ptr) : ptr(_ptr) {
            // copy each sysboundary to the device
            // loop through all sysboundaries in ptr->sysboundaries
            // and copy them to the device
            for (auto sbc : ptr->getSysBoundaries()) {
                if (sbc->getIndex() == sysboundarytype::SET_MAXWELLIAN) {
                    setmaxwellian = dynamic_cast<SBC::SetMaxwellian*>(sbc);
                } else if (sbc->getIndex() == sysboundarytype::IONOSPHERE) {
                    ionosphere = dynamic_cast<SBC::Ionosphere*>(sbc);
                } else if (sbc->getIndex() == sysboundarytype::OUTFLOW) {
                    outflow = dynamic_cast<SBC::Outflow*>(sbc);
                } else if (sbc->getIndex() == sysboundarytype::CONDUCTINGSPHERE) {
                    conductingsphere = dynamic_cast<SBC::Conductingsphere*>(sbc);
                } else if (sbc->getIndex() == sysboundarytype::DO_NOT_COMPUTE) {
                    doNotCompute = dynamic_cast<SBC::DoNotCompute*>(sbc);
                } else {
                    std::cerr << "ERROR: sysboundarytype not found" << std::endl;
                    exit(1);
                }
            }
        }

        buf(const buf& u) : ptr(u.ptr), is_copy(1), setmaxwellian(u.setmaxwellian), ionosphere(u.ionosphere), outflow(u.outflow), conductingsphere(u.conductingsphere) {}

        Proxy getSysBoundary(int sysBoundaryFlag) const {
            return Proxy(sysBoundaryFlag, this);
        }
};

}

#endif
