/*
 * This file is part of Vlasiator.
 * Copyright 2013, 2014 Finnish Meteorological Institute
 */

#include <cstdlib>
#include <cmath>

#include "common.h"
#include "parameters.h"
#include "amr_refinement_criteria.h"
#include "velocity_blocks.h"

using namespace std;

namespace amr_ref_criteria {
   
   Base::Base() { }   
   
   
   Base::~Base() { }
   
   
   Base* relDiffMaker() {return new RelativeDifference;}

   
   void Base::evaluate(const Realf* velBlost,Realf* result) {
      for (int i=0; i<WID3; ++i) result[i] = 0.0;
   }

   RelativeDifference::RelativeDifference() { }

   
   RelativeDifference::~RelativeDifference() { }


   Realf RelativeDifference::evaluate(const Realf* array) {
      // How many neighbor data points (per coordinate) the given block includes?
      const int PAD=1;
      Realf maxvalue = 0.0;

      for (int kc=0; kc<WID; ++kc) for (int jc=0; jc<WID; ++jc) for (int ic=0; ic<WID; ++ic) {
         Realf f_cen = array[vblock::padIndex<PAD>(ic+1,jc+1,kc+1)];
         if (fabs(f_cen) < Parameters::sparseMinValue) continue;

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


   void RelativeDifference::evaluate(const Realf* array,Realf* result) {
      const int PAD=1;
      for (int kc=0; kc<WID; ++kc) for (int jc=0; jc<WID; ++jc) for (int ic=0; ic<WID; ++ic) {
         Realf f_cen = array[vblock::padIndex<PAD>(ic+1,jc+1,kc+1)];
         if (fabs(f_cen) < Parameters::sparseMinValue) {
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

      /*if (df == numeric_limits<Realf>::infinity()) {
         cerr << "df is " << df << " df_max " << df_max << endl;
         cerr << "\t values are " << f_lft << '\t' << f_cen << '\t' << f_rgt << endl;
         cerr << "\t " << fabs(f_rgt-f_cen) << '\t' << fabs(f_cen-f_lft) << '\t' << max(fabs(f_rgt-f_cen),fabs(f_cen-f_lft)) << endl;
         cerr << "\t " << max(fabs(f_rgt-f_cen),fabs(f_cen-f_lft)) / ((f_cen + 1e-30)*df_max) << endl;
         exit(1);
      }*/

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

