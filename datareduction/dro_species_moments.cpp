/* This file is part of Vlasiator.
 * 
 * File:   dro_species_moments.cpp
 * Author: sandroos
 *
 * Created on March 23, 2015
 *
 * Copyright 2015 Finnish Meteorological Institute
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

      for (int popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
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
