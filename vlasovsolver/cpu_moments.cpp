/*
This file is part of Vlasiator.

Copyright 2015 Finnish Meteorological Institute

*/

#include <phiprof.hpp>
#include "cpu_moments.h"

using namespace std;

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
