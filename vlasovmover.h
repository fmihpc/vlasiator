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
                                 Real dt);

/** Calculate velocity moments for the given spatial cell.
 * This function is defined in cpu_moments.cpp file.*/
void calculateCellMoments(
        spatial_cell::SpatialCell* cell,
        const bool& computeSecond,const bool& doNotSkip=false);

/*!
  \brief Compute real-time 1st order accurate moments from the moments after propagation in velocity and spatial space
*/
 
void calculateInterpolatedVelocityMoments(
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   const int cp_rho,
   const int cp_rhovx,
   const int cp_rhovy,
   const int cp_rhovz,
   const int cp_p11,
   const int cp_p22,
   const int cp_p33
);

/*!
  \brief Compute 0th, 1st and 2nd velocity moments (RHO,RHOVX,RHOVY,RHOVZ,P_11,P_22,P_33 and *_DT2) for all cells in the grid directly from distribution function. The simulation should be at a true time-step! This is at the moment only called at initialisation.
  \param mpiGrid Grid of spatial cells for which moments are computed 
  
*/
void calculateInitialVelocityMoments(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid);



#endif


