/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://vlasiator.fmi.fi/
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
 * 
 * File:   dro_species_moments.cpp
 * Author: sandroos
 *
 * Created on March 23, 2015
 */

#include <cstdlib>
#include <iostream>

#ifdef _OPENMP
   #include <omp.h>
#endif

#include "../vlasovsolver/cpu_moments.h"
#include "dro_species_moments.h"

using namespace std;

namespace DRO {

   SpeciesMoments::SpeciesMoments(): DataReductionOperator() { }
   
   SpeciesMoments::~SpeciesMoments() { }

   bool SpeciesMoments::getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const {
      dataType = "float";
      dataSize =  sizeof(Real);
      vectorSize = 1;
      return true;
   }
   
   std::string SpeciesMoments::getName() const {
      return "SpeciesMoments";
   }
      
   bool SpeciesMoments::handlesWriting() const {return true;}
   
   bool SpeciesMoments::reduceData(const spatial_cell::SpatialCell* cell,char* buffer) {
      return false;
   }
     
   bool SpeciesMoments::reduceData(const spatial_cell::SpatialCell* cell,Real * result) {
      return false;
   }
      
   bool SpeciesMoments::setSpatialCell(const spatial_cell::SpatialCell* cell) {
      return true;
   }
   
   bool SpeciesMoments::writeData(const dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                  const std::vector<CellID>& cells,
                                  const std::string& meshName,
                                  vlsv::Writer& vlsvWriter) {
      phiprof::start("SpeciesMoments");
      bool success = true;
      
      Real* bufferRho   = new Real[cells.size()];
      Real* bufferRhoV  = new Real[cells.size()*3];
      uint* blocks      = new uint[cells.size()];
      size_t N_cells = 0;
      size_t thread_N_cells = 0;

      for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
         #pragma omp parallel for reduction(+:thread_N_cells)
         for (size_t c=0; c<cells.size(); ++c) {
            Real array[4];
            for (int i=0; i<4; ++i) array[i] = 0;

            spatial_cell::SpatialCell* cell = mpiGrid[cells[c]];
            vmesh::VelocityBlockContainer<vmesh::LocalID>& blockContainer = cell->get_velocity_blocks(popID);
            const Realf* data       = blockContainer.getData();
            const Real* blockParams = blockContainer.getParameters();

            for (vmesh::LocalID blockLID=0; blockLID<blockContainer.size(); ++blockLID) {
               blockVelocityFirstMoments(data+blockLID*WID3,
                                         blockParams+blockLID*BlockParams::N_VELOCITY_BLOCK_PARAMS,
                                         1.0,array);
            } // for-loop over velocity blocks
            
            bufferRho[c]      = array[0];
            bufferRhoV[c*3+0] = array[1];
            bufferRhoV[c*3+1] = array[2];
            bufferRhoV[c*3+2] = array[3];
            thread_N_cells += blockContainer.size()*WID3;
            
            blocks[c] = blockContainer.size();
         } // for-loop over spatial cells

         N_cells += thread_N_cells;
         
         map<string,string> attribs;
         attribs["mesh"] = meshName;
         attribs["name"] = getObjectWrapper().particleSpecies[popID].name + "/" + "n";
         if (vlsvWriter.writeArray("VARIABLE",attribs,cells.size(),1,bufferRho) == false) success = false;
            
         attribs["name"] = getObjectWrapper().particleSpecies[popID].name + "/" + "nV";
         if (vlsvWriter.writeArray("VARIABLE",attribs,cells.size(),3,bufferRhoV) == false) success = false;
         
         attribs["name"] = getObjectWrapper().particleSpecies[popID].name + "/" + "blocks";
         if (vlsvWriter.writeArray("VARIABLE",attribs,cells.size(),1,blocks) == false) success = false;
      } // for-loop over particle species

      delete [] bufferRho ; bufferRho  = NULL;
      delete [] bufferRhoV; bufferRhoV = NULL;
      delete [] blocks; blocks = NULL;

      phiprof::stop("SpeciesMoments",N_cells,"Phase-Space Cells");
      return success;
   }
   
} // namespace DRO
