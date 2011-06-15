#ifndef VLASOVMOVER_H
#define VLASOVMOVER_H

#include <vector>

#include "definitions.h"
#include "cell_spatial.h"

bool finalizeMover();

#ifdef PARGRID 

#include "pargrid.h"

bool initializeMover(ParGrid<SpatialCell>& mpiGrid);
void calculateSimParameters(ParGrid<SpatialCell>& mpiGrid,creal& t,Real& dt);
void calculateCellParameters(ParGrid<SpatialCell>& mpiGrid,creal& t,ID::type cell);
void calculateAcceleration(ParGrid<SpatialCell>& mpiGrid);
void calculateSpatialDerivatives(ParGrid<SpatialCell>& mpiGrid);
void calculateSpatialFluxes(ParGrid<SpatialCell>& mpiGrid);
void calculateSpatialPropagation(ParGrid<SpatialCell>& mpiGrid,const bool& secondStep,const bool& transferAvgs);
void initialLoadBalance(ParGrid<SpatialCell>& mpiGrid);
void calculateVelocityMoments(ParGrid<SpatialCell>& mpiGrid);

#else // ifdef PARGRID

#include <stdint.h>

#define DCCRG_SEND_SINGLE_CELLS
#define DCCRG_CELL_DATA_SIZE_FROM_USER
#define DCCRG_USER_MPI_DATA_TYPE
#include <dccrg.hpp>

bool initializeMover(dccrg<SpatialCell>& mpiGrid);
void calculateSimParameters(dccrg<SpatialCell>& mpiGrid,creal& t,Real& dt);
void calculateCellParameters(dccrg<SpatialCell>& mpiGrid,creal& t,uint64_t& cell);
void calculateAcceleration(dccrg<SpatialCell>& mpiGrid);
void calculateSpatialDerivatives(dccrg<SpatialCell>& mpiGrid);
void calculateSpatialFluxes(dccrg<SpatialCell>& mpiGrid);
void calculateSpatialPropagation(dccrg<SpatialCell>& mpiGrid,const bool& secondStep,const bool& transferAvgs);
void initialLoadBalance(dccrg<SpatialCell>& mpiGrid);
void calculateVelocityMoments(dccrg<SpatialCell>& mpiGrid);

#endif	// ifdef PARGRID

#endif


