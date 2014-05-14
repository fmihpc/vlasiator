/*
This file is part of Vlasiator.

Copyright 2010, 2011, 2012, 2013 Finnish Meteorological Institute
*/

#ifndef VLASOVMOVER_H
#define VLASOVMOVER_H

#include <vector>

#include "definitions.h"
#include "spatial_cell.hpp"
using namespace spatial_cell;

#include <stdint.h>
#include <dccrg.hpp>
#include <dccrg_cartesian_geometry.hpp>

void calculateAcceleration(
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   Real dt
);

void calculateSpatialTranslation(
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   creal dt);

/*!
  \brief Compute real-time 1st order accurate moments from the moments after propagation in velocity and spatial space
*/
 
void calculateInterpolatedVelocityMoments(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                           const int cp_rho, const int cp_rhovx, const int cp_rhovy, const int cp_rhovz);

/*!
   \brief Compute 0th and 1st velocity moments (RHO,RHOVX,RHOVY,RHOVZ) for a cell directly from distribution function. The simulation should be at a true time-step!
   \param SC pointer to the spatial cell
   \param doNotSkip Used to override the checks about system boundaries which in some cases cause the moments not to be calculated as they have been e.g. copied. Default value: false, in order to preserve all the calls using the checks.
*/
void calculateCellVelocityMoments(
   SpatialCell* SC,
   bool doNotSkip=false
);


/*!
  \brief Compute 0th and 1st velocity moments (RHO,RHOVX,RHOVY,RHOVZ) for all cells in the grid directly from distribution function.The simulation should be at a true time-step!
  \param mpiGrid Grid of spatial cells for which moments are computed 
  
*/
void calculateVelocityMoments(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid);



#endif


