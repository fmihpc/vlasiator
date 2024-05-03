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

#ifndef DATAREDUCER_H
#define DATAREDUCER_H

#include <vector>
#include "fsgrid.hpp"

#include "../definitions.h"
#include "../spatial_cell_wrapper.hpp"
#include "datareductionoperator.h"
#include "../sysboundary/ionosphere.h"

/** The purpose of DataReducer is to contain DRO::DataReductionOperators, and apply 
 * them to simulation data when writing output files. Files containing full 
 * distribution functions of every spatial cell require so much disk space 
 * that they cannot be written out so often as user would want. Thus, derived 
 * quantities need to be calculated for every spatial cell, which are then 
 * written to data files. This process is here called data reduction.
 */
class DataReducer {
 public:
   DataReducer();
   ~DataReducer();
   
   bool addOperator(DRO::DataReductionOperator* op);
   bool getDataVectorInfo(const unsigned int& operatorID,std::string& dataType,
                          unsigned int& dataSize,unsigned int& vectorSize) const;
   bool addMetadata(const unsigned int operatorID,std::string unit,std::string unitLaTeX,std::string variableLaTeX,std::string unitConversion);
   bool getMetadata(const unsigned int& operatorID,std::string& unit,std::string& unitLaTeX,std::string& variableLaTeX,std::string& unitConversion) const;

   std::string getName(const unsigned int& operatorID) const;
   bool hasParameters(const unsigned int& operatorID) const;
   bool reduceData(const SpatialCell* cell,const unsigned int& operatorID,char* buffer);
   bool reduceDiagnostic(const SpatialCell* cell,const unsigned int& operatorID,Real * result);
   unsigned int size() const;
   bool writeParameters(const unsigned int& operatorID, vlsv::Writer& vlsvWriter);
   bool writeFsGridData(
                      FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
                      FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
                      FsGrid< std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> & EHallGrid,
                      FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
                      FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
                      FsGrid< std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> & dPerBGrid,
                      FsGrid< std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> & dMomentsGrid,
                      FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
                      FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
                      FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
                      const std::string& meshName, const unsigned int operatorID,
                      vlsv::Writer& vlsvWriter,
                      const bool writeAsFloat = false);
   bool writeIonosphereGridData(SBC::SphericalTriGrid& grid, const std::string& meshName,
         const unsigned int operatorID, vlsv::Writer& vlsvWriter);

 private:
   /** Private copy-constructor to prevent copying the class.
    */
   DataReducer(const DataReducer& dr);
   
   std::vector<DRO::DataReductionOperator*> operators;
   /**< A container for all DRO::DataReductionOperators stored in DataReducer.*/
};

void initializeDataReducers(DataReducer * outputReducer, DataReducer * diagnosticReducer);

#endif
