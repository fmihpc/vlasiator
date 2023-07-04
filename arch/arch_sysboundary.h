#ifndef ARCH_SYSBOUNDARY_H
#define ARCH_SYSBOUNDARY_H

#include "arch_device_api.h"
#include "../sysboundary/sysboundary.h"
#include "../sysboundary/donotcompute.h"
#include "../sysboundary/ionosphere.h"
#include "../sysboundary/outflow.h"
#include "../sysboundary/setmaxwellian.h"

#if !defined(USE_CUDA) || !defined(CUDACC)

namespace arch {

template <>
class buf<SysBoundary> {
    private:  
        SysBoundary *ptr; 

        // list of possible sysboundaries
        SBC::SetMaxwellian* setmaxwellian;
        SBC::Ionosphere* ionosphere;
        SBC::Outflow* outflow;

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
                    return bufPtr->fieldSolverBoundaryCondMagneticFieldFactory(sysBoundaryFlag, bGrid, technicalGrid, i, j, k, dt, component);
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
                } else {
                    std::cerr << "ERROR: sysboundarytype not found" << std::endl;
                    exit(1);
                }
            }
        }

        buf(const buf& u) : ptr(u.ptr), is_copy(1), setmaxwellian(u.setmaxwellian), ionosphere(u.ionosphere), outflow(u.outflow) {}

        Proxy getSysBoundary(int sysBoundaryFlag) const {
            return Proxy(sysBoundaryFlag, this);
        }

        Real fieldSolverBoundaryCondMagneticFieldFactory (
            int sysBoundaryFlag, 
            const buf<FsGrid<Real, fsgrids::bfield::N_BFIELD, FS_STENCIL_WIDTH>>& bGrid,
            const buf<FsGrid<fsgrids::technical, 1, FS_STENCIL_WIDTH>>& technicalGrid, cint i, cint j, cint k,
            creal& dt, cuint& component) const {

            if (sysBoundaryFlag == sysboundarytype::SET_MAXWELLIAN) {
                return setmaxwellian->fieldSolverBoundaryCondMagneticField(bGrid, technicalGrid, i, j, k, dt,component);
            } else if (sysBoundaryFlag == sysboundarytype::IONOSPHERE) {
                return ionosphere->fieldSolverBoundaryCondMagneticField(bGrid, technicalGrid, i, j, k, dt,
                                                                             component);
            } else if (sysBoundaryFlag == sysboundarytype::OUTFLOW) {
                return outflow->fieldSolverBoundaryCondMagneticField(bGrid, technicalGrid, i, j, k, dt, component);
            } else {
                std::cerr << "ERROR: sysboundarytype not found" << std::endl;
                exit(1);
            }
        }
};

}

#endif

#endif