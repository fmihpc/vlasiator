/*
 This file is part of Vlasiator.
 
 Copyright 2010, 2011, 2012, 2013 Finnish Meteorological Institute
 
 
 
 
 
 
 
 
 
 
 
 
 */

#ifndef SETMAXWELLIAN_H
#define SETMAXWELLIAN_H

#include <vector>
#include "../definitions.h"
#include "../readparameters.h"
#include "../spatial_cell.hpp"
#include "sysboundarycondition.h"
#include "setbyuser.h"

using namespace std;

namespace SBC {
   /*!\brief SetMaxwellian is a class applying fixed Maxwellian conditions according to parameters read from an input file.
    * 
    * Maxwellian is a class handling cells tagged as sysboundarytype::MAXWELLIAN by this
    * system boundary condition.
    * 
    * It applies fixed Maxwellian settings to the system boundary cells, the parameters of
    * which are being read from an input file.
    * 
    * The class inherits most of its machinery from
    * SysBoundaryCondition::SetByUser. The parameters are more general than for Maxwellian
    * and could be put in SysBoundaryCondition::SetByUser but this way they can have a
    * specific prefix which is needed if several inheriting classes are needed.
    * 
    */
   class SetMaxwellian: public SetByUser {
   public:
      SetMaxwellian();
      virtual ~SetMaxwellian();
      
      static void addParameters();
      virtual void getParameters();
      
      virtual string getName() const;
      virtual uint getIndex() const;
      
   protected:
      void generateTemplateCell(spatial_cell::SpatialCell& templateCell, int inputDataIndex, creal& t);
      
      Real maxwellianDistribution(const int& popID,
         creal& rho, creal& T, creal& vx, creal& vy, creal& vz
      );
      
      vector<vmesh::GlobalID> findBlocksToInitialize(
         const int& popID,
         SpatialCell& cell,
         creal& rho,
         creal& T,
         creal& VX,
         creal& VY,
         creal& VZ
      );
      
      uint nSpaceSamples;
      uint nVelocitySamples;
   };
}

#endif
