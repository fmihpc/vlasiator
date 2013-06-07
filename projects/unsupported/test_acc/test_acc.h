/*
This file is part of Vlasiator.

Copyright 2011, 2012 Finnish Meteorological Institute












*/

#ifndef TEST_ACC_H
#define TEST_ACC_H

#include "definitions.h"
#include "spatial_cell.hpp"
#include "projects/projects_common.h"
#include "projects/projects_vlasov_acceleration.h"

#include "dccrg.hpp"

struct test_accParameters {
   static Real SPEED;
   static Real v_min;
   static Real v_max;
} ;

/**
 * Initialize project. Can be used, e.g., to read in parameters from the input file
 */
bool initializeProject(void);

/** Register parameters that should be read in
 */
bool addProjectParameters(void);
/** Get the value that was read in
 */
bool getProjectParameters(void);

/*!\brief Set the fields and distribution of a cell according to the default simulation settings.
 * This is used for the NOT_SYSBOUNDARY cells and some other system boundary conditions (e.g. Outflow).
 * \param cell Pointer to the cell to set.
 */
void setProjectCell(SpatialCell* cell);

template<typename UINT,typename REAL> void calcAccFaceX(
   REAL& ax, REAL& ay, REAL& az,
   const UINT& I, const UINT& J, const UINT& K,
   const REAL* const cellParams,
   const REAL* const blockParams
) {
   ax = blockParams[BlockParams::VXCRD] + (I+convert<REAL>(0.5))*blockParams[BlockParams::DVX];
   ay = blockParams[BlockParams::VYCRD] + (J+convert<REAL>(0.5))*blockParams[BlockParams::DVY];
   az = blockParams[BlockParams::VZCRD] + (K+convert<REAL>(0.5))*blockParams[BlockParams::DVZ];
   
//    const REAL SPEED = convert<REAL>(0.5) + convert<REAL>(0.5)*(convert<REAL>(1.0) / convert<REAL>(6.0));
   
   if (ax > convert<REAL>(0.0)) ax = test_accParameters::SPEED;
   else ax = -test_accParameters::SPEED;
   if (ay > convert<REAL>(0.0)) ay = test_accParameters::SPEED;
   else ay = -test_accParameters::SPEED;
   if (az > convert<REAL>(0.0)) az = test_accParameters::SPEED;
   else az = -test_accParameters::SPEED;
}

template<typename UINT,typename REAL> void calcAccFaceY(
   REAL& ax, REAL& ay, REAL& az,
   const UINT& I,const UINT& J,const UINT& K,
   const REAL* const cellParams,
   const REAL* const blockParams
) {
   ax = blockParams[BlockParams::VXCRD] + (I+convert<REAL>(0.5))*blockParams[BlockParams::DVX];
   ay = blockParams[BlockParams::VYCRD] + (J+convert<REAL>(0.5))*blockParams[BlockParams::DVY];
   az = blockParams[BlockParams::VZCRD] + (K+convert<REAL>(0.5))*blockParams[BlockParams::DVZ];
   
//    const REAL SPEED = convert<REAL>(0.5) + convert<REAL>(0.5)*(convert<REAL>(1.0) / convert<REAL>(6.0));
   
   if (ax > convert<REAL>(0.0)) ax = test_accParameters::SPEED;
   else ax = -test_accParameters::SPEED;
   if (ay > convert<REAL>(0.0)) ay = test_accParameters::SPEED;
   else ay = -test_accParameters::SPEED;
   if (az > convert<REAL>(0.0)) az = test_accParameters::SPEED;
   else az = -test_accParameters::SPEED;
}

template<typename UINT,typename REAL> void calcAccFaceZ(
   REAL& ax, REAL& ay, REAL& az,
   const UINT& I,const UINT& J,const UINT& K,
   const REAL* const cellParams,
   const REAL* const blockParams
) {
   ax = blockParams[BlockParams::VXCRD] + (I+convert<REAL>(0.5))*blockParams[BlockParams::DVX];
   ay = blockParams[BlockParams::VYCRD] + (J+convert<REAL>(0.5))*blockParams[BlockParams::DVY];
   az = blockParams[BlockParams::VZCRD] + (K+convert<REAL>(0.5))*blockParams[BlockParams::DVZ];
   
//    const REAL SPEED = convert<REAL>(0.5) + convert<REAL>(0.5)*(convert<REAL>(1.0) / convert<REAL>(6.0));
   
   if (ax > convert<REAL>(0.0)) ax = test_accParameters::SPEED;
   else ax = -test_accParameters::SPEED;
   if (ay > convert<REAL>(0.0)) ay = test_accParameters::SPEED;
   else ay = -test_accParameters::SPEED;
   if (az > convert<REAL>(0.0)) az = test_accParameters::SPEED;
   else az = -test_accParameters::SPEED;
}

#endif

