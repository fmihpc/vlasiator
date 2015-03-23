/* This file is part of Vlasiator.
 * 
 * File:   dro_species.cpp
 * Author: sandroos
 *
 * Created on March 23, 2015
 *
 * Copyright 2015 Finnish Meteorological Institute
 */

#ifndef DRO_SPECIES_H
#define	DRO_SPECIES_H

#include "datareductionoperator.h"
#include "../object_wrapper.h"

namespace DRO {

   class SpeciesMoments: public DataReductionOperator {
   public:
      SpeciesMoments();
      virtual ~SpeciesMoments();
     
      virtual bool getDataVectorInfo(std::string& dataType,unsigned int& dataSize,unsigned int& vectorSize) const;
      virtual std::string getName() const;
      virtual bool handlesWriting() const;
      virtual bool reduceData(const spatial_cell::SpatialCell* cell,char* buffer);
      virtual bool reduceData(const spatial_cell::SpatialCell* cell,Real * result);
      virtual bool setSpatialCell(const spatial_cell::SpatialCell* cell);
      virtual bool writeData(const dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                             const std::vector<CellID>& cells,const std::string& meshName,
                             vlsv::Writer& vlsvWriter);
     
   protected:
      Real averageVX, averageVY, averageVZ;
      Real rho_m;
   };

} // namespace DRO

#endif	/* DRO_SPECIES_H */

