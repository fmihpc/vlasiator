/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://www.physics.helsinki.fi/vlasiator/
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
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
   const int cp_rhom,
   const int cp_vx,
   const int cp_vy,
   const int cp_vz,
   const int cp_rhoq,
   const int cp_p11,
   const int cp_p22,
   const int cp_p33
);

/*!
  \brief Compute 0th, 1st and 2nd velocity moments (RHO,VX,VY,VZ,P_11,P_22,P_33 and *_DT2) for all cells in the grid directly from distribution function. The simulation should be at a true time-step! This is at the moment only called at initialisation.
  \param mpiGrid Grid of spatial cells for which moments are computed 
  
*/
void calculateInitialVelocityMoments(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid);



#endif


