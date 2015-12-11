/*
 * This file is part of Vlasiator.
 * Copyright 2013, 2014 Finnish Meteorological Institute
 */

#ifndef AMR_REFINEMENT_CRITERIA_H
#define AMR_REFINEMENT_CRITERIA_H

#include <iostream>
#include "definitions.h"

namespace amr_ref_criteria {
   
   class Base {
    public:
      Base();
      virtual ~Base();

      virtual Realf evaluate(const Realf* velBlock,const int& popID) = 0;
      virtual void evaluate(const Realf* velBlost,Realf* result,const int& popID);
      virtual bool initialize(const std::string& configRegion) = 0;
      
    protected:

   };

   void addRefinementCriteria();

   class RelativeDifference: public Base {
    public:
      RelativeDifference();
      ~RelativeDifference();
      
      Realf evaluate(const Realf* velBlock,const int& popID);
      void evaluate(const Realf* velBlost,Realf* result,const int& popID);
      bool initialize(const std::string& configRegion);

    protected:
      Realf df_max;
      
      Realf evaluate(const Realf& f_lef,const Realf& f_cen,const Realf& f_rgt);
   };

} // namespace amr_ref_criteria

#endif


