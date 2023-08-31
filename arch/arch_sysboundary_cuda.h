#ifndef ARCH_SYSBOUNDARY_CUDA_H
#define ARCH_SYSBOUNDARY_CUDA_H

#include "../sysboundary/sysboundary.h"
#include "../sysboundary/donotcompute.h"
#include "../sysboundary/ionosphere.h"
#include "../sysboundary/outflow.h"
#include "../sysboundary/setmaxwellian.h"
#include "../sysboundary/conductingsphere.h"

#include "arch_device_api.h"

namespace arch {

template <>
class buf<SysBoundary> {
    private:  
        SysBoundary *ptr; 

        // list of possible sysboundaries
        SBC::SetByUserFieldBoundary* setbyUser = NULL;
        SBC::IonosphereFieldBoundary* ionosphere = NULL;
        SBC::OutflowFieldBoundary* outflow = NULL;
        SBC::ConductingSphereFieldBoundary* conductingsphere = NULL; 
        SBC::DoNotCompute* doNotCompute;

        // list of possible sysboundaries on device
        SBC::SetByUserFieldBoundary* setByUser_d;
        SBC::IonosphereFieldBoundary* ionosphere_d;
        SBC::OutflowFieldBoundary* outflow_d;
        SBC::ConductingSphereFieldBoundary* conductingsphere_d;
        SBC::DoNotCompute* doNotCompute_d;

        uint is_copy = 0;

    public:   

        class Proxy {
            public:
                ARCH_HOSTDEV Proxy(int _sysBoundaryFlag, const buf<SysBoundary>* _ptr) : sysBoundaryFlag(_sysBoundaryFlag), bufPtr(_ptr) {}

                ARCH_HOSTDEV Real fieldSolverBoundaryCondMagneticField(
                    const buf<FsGrid<Real, fsgrids::bfield::N_BFIELD, FS_STENCIL_WIDTH>> & bGrid,
                    const buf<FsGrid< fsgrids::technical, 1, FS_STENCIL_WIDTH>> & technicalGrid,
                    cint i,
                    cint j,
                    cint k,
                    creal& dt,
                    cuint& component
                ) {
                    if (sysBoundaryFlag == sysboundarytype::SET_MAXWELLIAN) {
                        #ifdef __CUDA_ARCH__
                            return bufPtr->setByUser_d->fieldSolverBoundaryCondMagneticField(bGrid, technicalGrid, i, j, k, dt,component);
                        #else
                            return bufPtr->setbyUser->fieldSolverBoundaryCondMagneticField(bGrid, technicalGrid, i, j, k, dt,component);
                        #endif
                    } else if (sysBoundaryFlag == sysboundarytype::IONOSPHERE) {
                        #ifdef __CUDA_ARCH__
                            return bufPtr->ionosphere_d->fieldSolverBoundaryCondMagneticField(bGrid, technicalGrid, i, j, k, dt, component);
                        #else
                            return bufPtr->ionosphere->fieldSolverBoundaryCondMagneticField(bGrid, technicalGrid, i, j, k, dt, component);
                        #endif
                    } else if (sysBoundaryFlag == sysboundarytype::CONDUCTINGSPHERE) {
                        #ifdef __CUDA_ARCH__
                            return bufPtr->conductingsphere_d->fieldSolverBoundaryCondMagneticField(bGrid, technicalGrid, i, j, k, dt, component);
                        #else
                            return bufPtr->conductingsphere->fieldSolverBoundaryCondMagneticField(bGrid, technicalGrid, i, j, k, dt, component);
                        #endif
                    } else if (sysBoundaryFlag == sysboundarytype::OUTFLOW) {
                        #ifdef __CUDA_ARCH__
                            return bufPtr->outflow_d->fieldSolverBoundaryCondMagneticField(bGrid, technicalGrid, i, j, k, dt, component);
                        #else
                            return bufPtr->outflow->fieldSolverBoundaryCondMagneticField(bGrid, technicalGrid, i, j, k, dt, component);
                        #endif
                    } else {
                        #ifndef __CUDA_ARCH__
                            std::cerr << "ERROR: sysboundarytype not found bound cuda" << std::endl;
                            exit(1);
                        #endif
                    }
                }

                ARCH_HOSTDEV void fieldSolverBoundaryCondMagneticFieldProjection(
                    const arch::buf<FsGrid<Real, fsgrids::bfield::N_BFIELD, FS_STENCIL_WIDTH>> & bGrid,
                    const arch::buf<FsGrid< fsgrids::technical, 1, FS_STENCIL_WIDTH>> & technicalGrid,
                    cint i,
                    cint j,
                    cint k
                ) {
                    if (sysBoundaryFlag == sysboundarytype::SET_MAXWELLIAN) {
                        #ifdef __CUDA_ARCH__
                            bufPtr->setByUser_d->fieldSolverBoundaryCondMagneticFieldProjection(bGrid, technicalGrid, i, j, k);
                        #else
                            bufPtr->setbyUser->fieldSolverBoundaryCondMagneticFieldProjection(bGrid, technicalGrid, i, j, k);
                        #endif
                    } else if (sysBoundaryFlag == sysboundarytype::IONOSPHERE) {
                        #ifdef __CUDA_ARCH__
                            bufPtr->ionosphere_d->fieldSolverBoundaryCondMagneticFieldProjection(bGrid, technicalGrid, i, j, k);
                        #else
                            bufPtr->ionosphere->fieldSolverBoundaryCondMagneticFieldProjection(bGrid, technicalGrid, i, j, k);
                        #endif
                    } else if (sysBoundaryFlag == sysboundarytype::CONDUCTINGSPHERE) {
                        #ifdef __CUDA_ARCH__
                            bufPtr->conductingsphere_d->fieldSolverBoundaryCondMagneticFieldProjection(bGrid, technicalGrid, i, j, k);
                        #else
                            bufPtr->conductingsphere->fieldSolverBoundaryCondMagneticFieldProjection(bGrid, technicalGrid, i, j, k);
                        #endif
                    } else if (sysBoundaryFlag == sysboundarytype::OUTFLOW) {
                        #ifdef __CUDA_ARCH__
                            bufPtr->outflow_d->fieldSolverBoundaryCondMagneticFieldProjection(bGrid, technicalGrid, i, j, k);
                        #else
                            bufPtr->outflow->fieldSolverBoundaryCondMagneticFieldProjection(bGrid, technicalGrid, i, j, k);
                        #endif
                    } else {
                        #ifndef __CUDA_ARCH__
                            std::cerr << "ERROR: sysboundarytype not found bound cuda" << std::endl;
                            exit(1);
                        #endif
                    }
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
                    if (sysBoundaryFlag == sysboundarytype::SET_MAXWELLIAN) {
                        #ifdef __CUDA_ARCH__
                            bufPtr->setByUser_d->fieldSolverBoundaryCondDerivatives(dPerBGrid, dMomentsGrid, i, j, k, RKCase, component);
                        #else
                            bufPtr->setbyUser->fieldSolverBoundaryCondDerivatives(dPerBGrid, dMomentsGrid, i, j, k, RKCase, component);
                        #endif
                    } else if (sysBoundaryFlag == sysboundarytype::IONOSPHERE) {
                        #ifdef __CUDA_ARCH__
                            bufPtr->ionosphere_d->fieldSolverBoundaryCondDerivatives(dPerBGrid, dMomentsGrid, i, j, k, RKCase, component);
                        #else
                            bufPtr->ionosphere->fieldSolverBoundaryCondDerivatives(dPerBGrid, dMomentsGrid, i, j, k, RKCase, component);
                        #endif
                    } else if (sysBoundaryFlag == sysboundarytype::CONDUCTINGSPHERE) {
                        #ifdef __CUDA_ARCH__
                            bufPtr->conductingsphere_d->fieldSolverBoundaryCondDerivatives(dPerBGrid, dMomentsGrid, i, j, k, RKCase, component);
                        #else
                            bufPtr->conductingsphere->fieldSolverBoundaryCondDerivatives(dPerBGrid, dMomentsGrid, i, j, k, RKCase, component);
                        #endif
                    } else if (sysBoundaryFlag == sysboundarytype::OUTFLOW) {
                        #ifdef __CUDA_ARCH__
                            bufPtr->outflow_d->fieldSolverBoundaryCondDerivatives(dPerBGrid, dMomentsGrid, i, j, k, RKCase, component);
                        #else
                            bufPtr->outflow->fieldSolverBoundaryCondDerivatives(dPerBGrid, dMomentsGrid, i, j, k, RKCase, component);
                        #endif
                    } else {
                        #ifndef __CUDA_ARCH__
                            std::cerr << "ERROR: sysboundarytype not found bound cuda" << std::endl;
                            exit(1);
                        #endif
                    }
                }

                ARCH_HOSTDEV void fieldSolverBoundaryCondGradPeElectricField(
                    const arch::buf<FsGrid<Real, fsgrids::egradpe::N_EGRADPE, FS_STENCIL_WIDTH>> & EGradPeGrid,
                    cint i,
                    cint j,
                    cint k,
                    cuint component
                ) {
                    if (sysBoundaryFlag == sysboundarytype::SET_MAXWELLIAN) {
                        #ifdef __CUDA_ARCH__
                            bufPtr->setByUser_d->fieldSolverBoundaryCondGradPeElectricField(EGradPeGrid, i, j, k, component);
                        #else
                            bufPtr->setbyUser->fieldSolverBoundaryCondGradPeElectricField(EGradPeGrid, i, j, k, component);
                        #endif
                    } else if (sysBoundaryFlag == sysboundarytype::IONOSPHERE) {
                        #ifdef __CUDA_ARCH__
                            bufPtr->ionosphere_d->fieldSolverBoundaryCondGradPeElectricField(EGradPeGrid, i, j, k, component);
                        #else
                            bufPtr->ionosphere->fieldSolverBoundaryCondGradPeElectricField(EGradPeGrid, i, j, k, component);
                        #endif
                    } else if (sysBoundaryFlag == sysboundarytype::CONDUCTINGSPHERE) {
                        #ifdef __CUDA_ARCH__
                            bufPtr->conductingsphere_d->fieldSolverBoundaryCondGradPeElectricField(EGradPeGrid, i, j, k, component);
                        #else
                            bufPtr->conductingsphere->fieldSolverBoundaryCondGradPeElectricField(EGradPeGrid, i, j, k, component);
                        #endif
                    } else if (sysBoundaryFlag == sysboundarytype::OUTFLOW) {
                        #ifdef __CUDA_ARCH__
                            bufPtr->outflow_d->fieldSolverBoundaryCondGradPeElectricField(EGradPeGrid, i, j, k, component);
                        #else
                            bufPtr->outflow->fieldSolverBoundaryCondGradPeElectricField(EGradPeGrid, i, j, k, component);
                        #endif
                    } else {
                        #ifndef __CUDA_ARCH__
                            std::cerr << "ERROR: sysboundarytype not found bound cuda" << std::endl;
                            exit(1);
                        #endif
                    }
                }

                ARCH_HOSTDEV void fieldSolverBoundaryCondHallElectricField(
                    const arch::buf<FsGrid<Real, fsgrids::ehall::N_EHALL, FS_STENCIL_WIDTH>> & EHallGrid,
                    cint i,
                    cint j,
                    cint k,
                    cuint component
                ) { 
                    if (sysBoundaryFlag == sysboundarytype::SET_MAXWELLIAN) {
                        #ifdef __CUDA_ARCH__
                            bufPtr->setByUser_d->fieldSolverBoundaryCondHallElectricField(EHallGrid, i, j, k, component);
                        #else
                            bufPtr->setbyUser->fieldSolverBoundaryCondHallElectricField(EHallGrid, i, j, k, component);
                        #endif
                    } else if (sysBoundaryFlag == sysboundarytype::IONOSPHERE) {
                        #ifdef __CUDA_ARCH__
                            bufPtr->ionosphere_d->fieldSolverBoundaryCondHallElectricField(EHallGrid, i, j, k, component);
                        #else
                            bufPtr->ionosphere->fieldSolverBoundaryCondHallElectricField(EHallGrid, i, j, k, component);
                        #endif
                    } else if (sysBoundaryFlag == sysboundarytype::CONDUCTINGSPHERE) {
                        #ifdef __CUDA_ARCH__
                            bufPtr->conductingsphere_d->fieldSolverBoundaryCondHallElectricField(EHallGrid, i, j, k, component);
                        #else
                            bufPtr->conductingsphere->fieldSolverBoundaryCondHallElectricField(EHallGrid, i, j, k, component);
                        #endif
                    } else if (sysBoundaryFlag == sysboundarytype::OUTFLOW) {
                        #ifdef __CUDA_ARCH__
                            bufPtr->outflow_d->fieldSolverBoundaryCondHallElectricField(EHallGrid, i, j, k, component);
                        #else
                            bufPtr->outflow->fieldSolverBoundaryCondHallElectricField(EHallGrid, i, j, k, component);
                        #endif
                    } else {
                        #ifndef __CUDA_ARCH__
                            std::cerr << "ERROR: sysboundarytype not found bound cuda" << std::endl;
                            exit(1);
                        #endif
                    }
                }


                ARCH_HOSTDEV void fieldSolverBoundaryCondElectricField(
                    const arch::buf<FsGrid<Real, fsgrids::efield::N_EFIELD, FS_STENCIL_WIDTH>> & EGrid,
                    cint i,
                    cint j,
                    cint k,
                    cuint component
                ) {
                    if (sysBoundaryFlag == sysboundarytype::SET_MAXWELLIAN) {
                        #ifdef __CUDA_ARCH__
                            bufPtr->setByUser_d->fieldSolverBoundaryCondElectricField(EGrid, i, j, k, component);
                        #else
                            bufPtr->setbyUser->fieldSolverBoundaryCondElectricField(EGrid, i, j, k, component);
                        #endif
                    } else if (sysBoundaryFlag == sysboundarytype::IONOSPHERE) {
                        #ifdef __CUDA_ARCH__
                            bufPtr->ionosphere_d->fieldSolverBoundaryCondElectricField(EGrid, i, j, k, component);
                        #else
                            bufPtr->ionosphere->fieldSolverBoundaryCondElectricField(EGrid, i, j, k, component);
                        #endif
                    } else if (sysBoundaryFlag == sysboundarytype::CONDUCTINGSPHERE) {
                        #ifdef __CUDA_ARCH__
                            bufPtr->conductingsphere_d->fieldSolverBoundaryCondElectricField(EGrid, i, j, k, component);
                        #else
                            bufPtr->conductingsphere->fieldSolverBoundaryCondElectricField(EGrid, i, j, k, component);
                        #endif
                    } else if (sysBoundaryFlag == sysboundarytype::OUTFLOW) {
                        #ifdef __CUDA_ARCH__
                            bufPtr->outflow_d->fieldSolverBoundaryCondElectricField(EGrid, i, j, k, component);
                        #else
                            bufPtr->outflow->fieldSolverBoundaryCondElectricField(EGrid, i, j, k, component);
                        #endif
                    } else {
                        #ifndef __CUDA_ARCH__
                            std::cerr << "ERROR: sysboundarytype not found bound cuda" << std::endl;
                            exit(1);
                        #endif
                    }
                }

                ARCH_HOSTDEV void fieldSolverBoundaryCondBVOLDerivatives(
                    const arch::buf<FsGrid<Real, fsgrids::volfields::N_VOL, FS_STENCIL_WIDTH>> & volGrid,
                    cint i,
                    cint j,
                    cint k,
                    cuint& component
                ) {
                    if (sysBoundaryFlag == sysboundarytype::SET_MAXWELLIAN) {
                        #ifdef __CUDA_ARCH__
                            bufPtr->setByUser_d->fieldSolverBoundaryCondBVOLDerivatives(volGrid, i, j, k, component);
                        #else
                            bufPtr->setbyUser->fieldSolverBoundaryCondBVOLDerivatives(volGrid, i, j, k, component);
                        #endif
                    } else if (sysBoundaryFlag == sysboundarytype::IONOSPHERE) {
                        #ifdef __CUDA_ARCH__
                            bufPtr->ionosphere_d->fieldSolverBoundaryCondBVOLDerivatives(volGrid, i, j, k, component);
                        #else
                            bufPtr->ionosphere->fieldSolverBoundaryCondBVOLDerivatives(volGrid, i, j, k, component);
                        #endif
                    } else if (sysBoundaryFlag == sysboundarytype::CONDUCTINGSPHERE) {
                        #ifdef __CUDA_ARCH__
                            bufPtr->conductingsphere_d->fieldSolverBoundaryCondBVOLDerivatives(volGrid, i, j, k, component);
                        #else
                            bufPtr->conductingsphere->fieldSolverBoundaryCondBVOLDerivatives(volGrid, i, j, k, component);
                        #endif
                    } else if (sysBoundaryFlag == sysboundarytype::OUTFLOW) {
                        #ifdef __CUDA_ARCH__
                            bufPtr->outflow_d->fieldSolverBoundaryCondBVOLDerivatives(volGrid, i, j, k, component);
                        #else
                            bufPtr->outflow->fieldSolverBoundaryCondBVOLDerivatives(volGrid, i, j, k, component);
                        #endif
                    } else {
                        #ifndef __CUDA_ARCH__
                            std::cerr << "ERROR: sysboundarytype not found bound cuda" << std::endl;
                            exit(1);
                        #endif
                    }
                } 
 
            private:
                int sysBoundaryFlag;
                const buf<SysBoundary>* bufPtr;

        };

        // sync the data from the host to the device
        void syncHostToDevice(void) {
            if (setbyUser != NULL) {
                // copy the data from the host to the device
                CHK_ERR(cudaMemcpy(setByUser_d, setbyUser, sizeof(SBC::SetByUserFieldBoundary), cudaMemcpyHostToDevice));
            }
            if (ionosphere != NULL) {
                // copy the data from the host to the device
                CHK_ERR(cudaMemcpy(ionosphere_d, ionosphere, sizeof(SBC::Ionosphere), cudaMemcpyHostToDevice));
            }
            if (outflow != NULL) {
                // copy the data from the host to the device
                CHK_ERR(cudaMemcpy(outflow_d, outflow, sizeof(SBC::Outflow), cudaMemcpyHostToDevice));
            }
        }

        // sync the data from the device to the host
        void syncDeviceToHost(void) {
            if (setbyUser != NULL) {
                // copy the data from the device to the host
                CHK_ERR(cudaMemcpy(setbyUser, setByUser_d, sizeof(SBC::SetByUserFieldBoundary), cudaMemcpyDeviceToHost));
            }
            if (ionosphere != NULL) {
                // copy the data from the device to the host
                CHK_ERR(cudaMemcpy(ionosphere, ionosphere_d, sizeof(SBC::Ionosphere), cudaMemcpyDeviceToHost));
            }
            if (outflow != NULL) {
                // copy the data from the device to the host
                CHK_ERR(cudaMemcpy(outflow, outflow_d, sizeof(SBC::Outflow), cudaMemcpyDeviceToHost));
            }
        }

        buf(SysBoundary* const _ptr) : ptr(_ptr) {
            // copy each sysboundary to the device
            // loop through all sysboundaries in ptr->sysboundaries
            // and copy them to the device
            for (auto sbc : ptr->getSysBoundaries()) {
                if (sbc->getIndex() == sysboundarytype::SET_MAXWELLIAN) {
                    setbyUser = dynamic_cast<SBC::SetByUser*>(sbc)->getFieldBoundary();
                    CHK_ERR(cudaMalloc(&setByUser_d, sizeof(SBC::SetByUserFieldBoundary)));
                } else if (sbc->getIndex() == sysboundarytype::IONOSPHERE) {
                    ionosphere = dynamic_cast<SBC::Ionosphere*>(sbc)->getFieldBoundary();
                    CHK_ERR(cudaMalloc(&ionosphere_d, sizeof(SBC::Ionosphere)));
                } else if (sbc->getIndex() == sysboundarytype::OUTFLOW) {
                    outflow = dynamic_cast<SBC::Outflow*>(sbc)->getFieldBoundary();
                    CHK_ERR(cudaMalloc(&outflow_d, sizeof(SBC::Outflow)));
                } else if (sbc->getIndex() == sysboundarytype::CONDUCTINGSPHERE) {
                    conductingsphere = dynamic_cast<SBC::Conductingsphere*>(sbc)->getFieldBoundary();
                    CHK_ERR(cudaMalloc(&conductingsphere_d, sizeof(SBC::Conductingsphere)));
                } else if (sbc->getIndex() == sysboundarytype::DO_NOT_COMPUTE) {
                    doNotCompute = dynamic_cast<SBC::DoNotCompute*>(sbc);
                } else {
                    std::cerr << "ERROR: sysboundarytype not found buf cuda" << std::endl;
                    exit(1);
                }

                syncHostToDevice(); 
            }
        }

        ARCH_HOSTDEV buf(const buf& u) : ptr(u.ptr), is_copy(1), setbyUser(u.setbyUser), ionosphere(u.ionosphere), outflow(u.outflow), setByUser_d(u.setByUser_d), ionosphere_d(u.ionosphere_d), outflow_d(u.outflow_d), doNotCompute(u.doNotCompute), doNotCompute_d(u.doNotCompute_d) {}

        ARCH_HOSTDEV Proxy getSysBoundary(int sysBoundaryFlag) const {
            return Proxy(sysBoundaryFlag, this);
        }
};

}

#endif
