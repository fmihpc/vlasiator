/*
This file is part of Vlasiator.

Copyright 2015 Finnish Meteorological Institute

*/

#include <phiprof.hpp>
#include "cpu_moments.h"
#include "../vlasovmover.h"

using namespace std;

void calculateCellMoments(
        spatial_cell::SpatialCell* cell,
        //dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
        //const CellID& cellID,
        const bool& computeSecond,
        const bool& doNotSkip) {

//        const CellID cellID = cells[c];
        //SpatialCell* cell = mpiGrid[cellID];
        
        // if doNotSkip == true then the first clause is false and we will never return,
        // i.e. always compute, otherwise we skip DO_NOT_COMPUTE cells
        // or boundary cells of layer larger than 1.
        bool skipMoments = false;
        if (!doNotSkip &&
            (cell->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE ||
            (cell->sysBoundaryLayer != 1  &&
             cell->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY))
            ) {
            skipMoments = true;
        }

        // Clear old moments to zero value
        if (skipMoments == false) {
            cell->parameters[CellParams::RHO  ] = 0.0;
            cell->parameters[CellParams::RHOVX] = 0.0;
            cell->parameters[CellParams::RHOVY] = 0.0;
            cell->parameters[CellParams::RHOVZ] = 0.0;
            cell->parameters[CellParams::P_11] = 0.0;
            cell->parameters[CellParams::P_22] = 0.0;
            cell->parameters[CellParams::P_33] = 0.0;
        }

        // Loop over all particle species
        if (skipMoments == false) {
        //for (int popID=0; popID<getObjectWrapper.particleSpecies.size(); ++popID) {
            const int popID=0;        
            const vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& vmesh = cell->get_velocity_mesh(popID);
            vmesh::VelocityBlockContainer<vmesh::LocalID>& blockContainer    = cell->get_velocity_blocks(popID);
            const Realf* data       = blockContainer.getData();
            const Real* blockParams = blockContainer.getParameters();
            
            // Calculate species' contribution to first velocity moments
            for (vmesh::LocalID blockLID=0; blockLID<vmesh.size(); ++blockLID) {
                cpu_blockVelocityFirstMoments(
                        data+blockLID*WID3,
                        blockParams+blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS,
                        cell->parameters,
                        CellParams::RHO,
                        CellParams::RHOVX,
                        CellParams::RHOVY,
                        CellParams::RHOVZ);
            }
        //}
        }
        
        // Compute second moments only if requested
        if (computeSecond == false) {
            skipMoments = true;
        }
            
        // Loop over all particle species
        if (skipMoments == false) {
        //for (int popID=0; popID<getObjectWrapper.particleSpecies.size(); ++popID) {
            const int popID=0;
            const vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& vmesh = cell->get_velocity_mesh(popID);
            vmesh::VelocityBlockContainer<vmesh::LocalID>& blockContainer    = cell->get_velocity_blocks(popID);
            const Realf* data       = blockContainer.getData();
            const Real* blockParams = blockContainer.getParameters();
            
            // Calculate species' contribution to second velocity moments
            for (vmesh::LocalID blockLID=0; blockLID<vmesh.size(); ++blockLID) {
                cpu_blockVelocitySecondMoments(
                        data+blockLID*WID3,
                        blockParams+blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS,
                        cell->parameters,
                        CellParams::RHO,
                        CellParams::RHOVX,
                        CellParams::RHOVY,
                        CellParams::RHOVZ,
                        CellParams::P_11,
                        CellParams::P_22,
                        CellParams::P_33);
            }
        //}
        }
    //} // for-loop over spatial cells
}

void calculateMoments_R_maxdt(
        dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
        const std::vector<CellID>& cells,
        const bool& computeSecond) {
 
    phiprof::start("compute-moments-n-maxdt");
    creal HALF = 0.5;
    
    #pragma omp parallel for
    for (size_t c=0; c<cells.size(); ++c) {
        const CellID cellID = cells[c];
        SpatialCell* cell = mpiGrid[cells[c]];
    
        // Clear old moments to zero value
        cell->parameters[CellParams::RHO_R  ] = 0.0;
        cell->parameters[CellParams::RHOVX_R] = 0.0;
        cell->parameters[CellParams::RHOVY_R] = 0.0;
        cell->parameters[CellParams::RHOVZ_R] = 0.0;
        cell->parameters[CellParams::P_11_R] = 0.0;
        cell->parameters[CellParams::P_22_R] = 0.0;
        cell->parameters[CellParams::P_33_R] = 0.0;
        
        const Real dx = cell->parameters[CellParams::DX];
        const Real dy = cell->parameters[CellParams::DY];
        const Real dz = cell->parameters[CellParams::DZ];
        
        //Reset spatial max DT
        cell->parameters[CellParams::MAXRDT] = numeric_limits<Real>::max();
        
        // Loop over all particle species
        //for (int popID=0; popID<getObjectWrapper.particleSpecies.size(); ++popID) {
            const int popID=0;        
            const vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& vmesh = cell->get_velocity_mesh(popID);
            vmesh::VelocityBlockContainer<vmesh::LocalID>& blockContainer    = cell->get_velocity_blocks(popID);
            const Realf* data       = blockContainer.getData();
            const Real* blockParams = blockContainer.getParameters();
            
            // compute maximum dt. Algorithm has a CFL condition, since it
            // is written only for the case where we have a stencil
            // supporting max translation of one cell
            const Real EPS = numeric_limits<Real>::min()*1000;
            for (unsigned int i=0; i<WID;i+=WID-1) {
                const Real Vx = blockParams[BlockParams::VXCRD] + (i+HALF)*blockParams[BlockParams::DVX] + EPS;
                const Real Vy = blockParams[BlockParams::VYCRD] + (i+HALF)*blockParams[BlockParams::DVY] + EPS;
                const Real Vz = blockParams[BlockParams::VZCRD] + (i+HALF)*blockParams[BlockParams::DVZ] + EPS;

                Real dt_max_cell = min(dx/fabs(Vx),min(dy/fabs(Vy),dz/fabs(Vz)));
                cell->parameters[CellParams::MAXRDT] = min(dt_max_cell,cell->parameters[CellParams::MAXRDT]);
            }

            // Calculate species' contribution to first velocity moments
            for (vmesh::LocalID blockLID=0; blockLID<vmesh.size(); ++blockLID) {
                cpu_blockVelocityFirstMoments(
                        data+blockLID*WID3,
                        blockParams+blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS,
                        cell->parameters,
                        CellParams::RHO_R,
                        CellParams::RHOVX_R,
                        CellParams::RHOVY_R,
                        CellParams::RHOVZ_R);
            }
        //}
        
        // Compute second moments only if requested
        if (computeSecond == false) continue;
            
        // Loop over all particle species
        //for (int popID=0; popID<getObjectWrapper.particleSpecies.size(); ++popID) {
            //const vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& vmesh = cell->get_velocity_mesh(popID);
            //vmesh::VelocityBlockContainer<vmesh::LocalID>& blockContainer    = cell->get_velocity_blocks(popID);
            //const Realf* data       = blockContainer.getData();
            //const Real* blockParams = blockContainer.getParameters();
            
            // Calculate species' contribution to second velocity moments
            for (vmesh::LocalID blockLID=0; blockLID<vmesh.size(); ++blockLID) {
                cpu_blockVelocitySecondMoments(
                        data+blockLID*WID3,
                        blockParams+blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS,
                        cell->parameters,
                        CellParams::RHO_R,
                        CellParams::RHOVX_R,
                        CellParams::RHOVY_R,
                        CellParams::RHOVZ_R,
                        CellParams::P_11_R,
                        CellParams::P_22_R,
                        CellParams::P_33_R);
            }
        //}
    }
    phiprof::stop("compute-moments-n-maxdt");
}

void calculateMoments_V(
        dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
        const std::vector<CellID>& cells,
        const bool& computeSecond) {
 
    phiprof::start("Compute moments");

    #pragma omp parallel for
    for (size_t c=0; c<cells.size(); ++c) {
        const CellID cellID = cells[c];
        SpatialCell* cell = mpiGrid[cells[c]];
    
        // Clear old moments to zero value
        cell->parameters[CellParams::RHO_V  ] = 0.0;
        cell->parameters[CellParams::RHOVX_V] = 0.0;
        cell->parameters[CellParams::RHOVY_V] = 0.0;
        cell->parameters[CellParams::RHOVZ_V] = 0.0;
        cell->parameters[CellParams::P_11_V] = 0.0;
        cell->parameters[CellParams::P_22_V] = 0.0;
        cell->parameters[CellParams::P_33_V] = 0.0;
        
        // Loop over all particle species
        //for (int popID=0; popID<getObjectWrapper.particleSpecies.size(); ++popID) {
            const int popID=0;        
            const vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& vmesh = cell->get_velocity_mesh(popID);
            vmesh::VelocityBlockContainer<vmesh::LocalID>& blockContainer    = cell->get_velocity_blocks(popID);
            const Realf* data       = blockContainer.getData();
            const Real* blockParams = blockContainer.getParameters();
            
            // Calculate species' contribution to first velocity moments
            for (vmesh::LocalID blockLID=0; blockLID<vmesh.size(); ++blockLID) {
                cpu_blockVelocityFirstMoments(
                        data+blockLID*WID3,
                        blockParams+blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS,
                        cell->parameters,
                        CellParams::RHO_V,
                        CellParams::RHOVX_V,
                        CellParams::RHOVY_V,
                        CellParams::RHOVZ_V);
            }
        //}
        
        // Compute second moments only if requested
        if (computeSecond == false) continue;
            
        // Loop over all particle species
        //for (int popID=0; popID<getObjectWrapper.particleSpecies.size(); ++popID) {
            //const vmesh::VelocityMesh<vmesh::GlobalID,vmesh::LocalID>& vmesh = cell->get_velocity_mesh(popID);
            //vmesh::VelocityBlockContainer<vmesh::LocalID>& blockContainer    = cell->get_velocity_blocks(popID);
            //const Realf* data       = blockContainer.getData();
            //const Real* blockParams = blockContainer.getParameters();
            
            // Calculate species' contribution to second velocity moments
            for (vmesh::LocalID blockLID=0; blockLID<vmesh.size(); ++blockLID) {
                cpu_blockVelocitySecondMoments(
                        data+blockLID*WID3,
                        blockParams+blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS,
                        cell->parameters,
                        CellParams::RHO_V,
                        CellParams::RHOVX_V,
                        CellParams::RHOVY_V,
                        CellParams::RHOVZ_V,
                        CellParams::P_11_V,
                        CellParams::P_22_V,
                        CellParams::P_33_V);
            }
        //}
    }
    phiprof::stop("Compute moments");
}
