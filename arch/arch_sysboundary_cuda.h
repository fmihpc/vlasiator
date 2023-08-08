#ifndef ARCH_SYSBOUNDARY_CUDA_H
#define ARCH_SYSBOUNDARY_CUDA_H

#include "../sysboundary/sysboundary.h"
#include "../sysboundary/donotcompute.h"
#include "../sysboundary/ionosphere.h"
#include "../sysboundary/outflow.h"
#include "../sysboundary/setmaxwellian.h"

#include "arch_device_api.h"

namespace arch {

template <>
class buf<SysBoundary> {
    private:  
        SysBoundary *ptr; 

        // list of possible sysboundaries
        SBC::SetByUserFieldBoundary* setbyUser = NULL;
        SBC::Ionosphere* ionosphere = NULL;
        SBC::Outflow* outflow = NULL;

        // list of possible sysboundaries on device
        SBC::SetByUserFieldBoundary* setByUser_d;
        SBC::Ionosphere* ionosphere_d;
        SBC::Outflow* outflow_d;

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
                    // return 0;
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
                    } else if (sysBoundaryFlag == sysboundarytype::OUTFLOW) {
                        #ifdef __CUDA_ARCH__
                            return bufPtr->outflow->fieldSolverBoundaryCondMagneticField(bGrid, technicalGrid, i, j, k, dt, component);
                        #else
                            return bufPtr->outflow->fieldSolverBoundaryCondMagneticField(bGrid, technicalGrid, i, j, k, dt, component);
                        #endif
                    } else {
                        #ifndef __CUDA_ARCH__
                            std::cerr << "ERROR: sysboundarytype not found" << std::endl;
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
            fflush(stdout);
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
                    ionosphere = dynamic_cast<SBC::Ionosphere*>(sbc);
                    CHK_ERR(cudaMalloc(&ionosphere_d, sizeof(SBC::Ionosphere)));
                } else if (sbc->getIndex() == sysboundarytype::OUTFLOW) {
                    outflow = dynamic_cast<SBC::Outflow*>(sbc);
                    CHK_ERR(cudaMalloc(&outflow_d, sizeof(SBC::Outflow)));
                } else {
                    std::cerr << "ERROR: sysboundarytype not found" << std::endl;
                    exit(1);
                }

                syncHostToDevice(); 
            }
        }

        ARCH_HOSTDEV buf(const buf& u) : ptr(u.ptr), is_copy(1), setbyUser(u.setbyUser), ionosphere(u.ionosphere), outflow(u.outflow), setByUser_d(u.setByUser_d), ionosphere_d(u.ionosphere_d), outflow_d(u.outflow_d) {}

        ARCH_HOSTDEV Proxy getSysBoundary(int sysBoundaryFlag) const {
            return Proxy(sysBoundaryFlag, this);
        }
};

}

#endif