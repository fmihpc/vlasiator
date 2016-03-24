/*
 * This file is part of Vlasiator.
 * Copyright 2013, 2014 Finnish Meteorological Institute
 */

#include <cstdlib>
#include <cmath>

#include "parameters.h"
#include "amr_refinement_criteria.h"
#include "velocity_blocks.h"
#include "object_wrapper.h"

using namespace std;

namespace amr_ref_criteria {
   
   Base::Base() { }   
   
   Base::~Base() { }
   
   Base* relDiffMaker() {return new RelativeDifference;}
   
   void Base::evaluate(const Realf* velBlost,Realf* result,const int& popID) {
      for (unsigned int i=0; i<WID3; ++i) result[i] = 0.0;
   }

   RelativeDifference::RelativeDifference() { }
   
   RelativeDifference::~RelativeDifference() { }

   Realf RelativeDifference::evaluate(const Realf* array,const int& popID) {
      // How many neighbor data points (per coordinate) the given block includes?
      const int PAD=1;
      Realf maxvalue = 0.0;

      for (unsigned int kc=0; kc<WID; ++kc) for (unsigned int jc=0; jc<WID; ++jc) for (unsigned int ic=0; ic<WID; ++ic) {
         Realf f_cen = array[vblock::padIndex<PAD>(ic+1,jc+1,kc+1)];
         
         #warning In here should we use SpatialCell::getVeloctyBlockMinValue()?
         if (fabs(f_cen) < getObjectWrapper().particleSpecies[popID].sparseMinValue) continue;

         Realf f_lft = array[vblock::padIndex<PAD>(ic  ,jc+1,kc+1)];
         Realf f_rgt = array[vblock::padIndex<PAD>(ic+2,jc+1,kc+1)];

         Realf df = evaluate(f_lft,f_cen,f_rgt);
         if (df > maxvalue) maxvalue = df;

         f_lft = array[vblock::padIndex<PAD>(ic+1,jc  ,kc+1)];
         f_rgt = array[vblock::padIndex<PAD>(ic+1,jc+2,kc+1)];
         df = evaluate(f_lft,f_cen,f_rgt);
         if (df > maxvalue) maxvalue = df;

         //f_lft = array[vblock::padIndex<PAD>(ic+1,jc+1,kc  )];
         //f_rgt = array[vblock::padIndex<PAD>(ic+1,jc+1,kc+2)];
         //df = evaluate(f_lft,f_cen,f_rgt);
         //if (df > maxvalue) maxvalue = df;
      }
      
      return maxvalue;
   }

   void RelativeDifference::evaluate(const Realf* array,Realf* result,const int& popID) {
      const int PAD=1;
      for (unsigned int kc=0; kc<WID; ++kc) for (unsigned int jc=0; jc<WID; ++jc) for (unsigned int ic=0; ic<WID; ++ic) {
         Realf f_cen = array[vblock::padIndex<PAD>(ic+1,jc+1,kc+1)];
         #warning In here should we use SpatialCell::getVeloctyBlockMinValue()?
         if (fabs(f_cen) < getObjectWrapper().particleSpecies[popID].sparseMinValue) {
            result[vblock::index(ic,jc,kc)] = 0;
            continue;
         }
         
         Realf f_lft = array[vblock::padIndex<PAD>(ic  ,jc+1,kc+1)];
         Realf f_rgt = array[vblock::padIndex<PAD>(ic+2,jc+1,kc+1)];
         Realf df = evaluate(f_lft,f_cen,f_rgt);
         
         f_lft = array[vblock::padIndex<PAD>(ic+1,jc  ,kc+1)];
         f_rgt = array[vblock::padIndex<PAD>(ic+1,jc+2,kc+1)];
         df = max(df,evaluate(f_lft,f_cen,f_rgt));

         //f_lft = array[vblock::padIndex<PAD>(ic+1,jc+1,kc  )];
         //f_rgt = array[vblock::padIndex<PAD>(ic+1,jc+1,kc+2)];
         //df = max(df,evaluate(f_lft,f_cen,f_rgt));
         
         result[vblock::index(ic,jc,kc)] = df;
      }
   }

   Realf RelativeDifference::evaluate(const Realf& f_lft,const Realf& f_cen,const Realf& f_rgt) {
      Realf df = max(fabs(f_rgt-f_cen),fabs(f_cen-f_lft));
      df = df / ((f_cen + 1e-30)*df_max);

      return df;
   }

   bool RelativeDifference::initialize(const std::string& configRegion) {
      //df_max = 8.0;
      df_max = 1.0;
      return true;
   }

   void addRefinementCriteria() {
      getObjectWrapper().amrVelRefCriteria.add("relative_difference",relDiffMaker);
   }
}

