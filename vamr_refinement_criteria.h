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

#ifndef VAMR_REFINEMENT_CRITERIA_H
#define VAMR_REFINEMENT_CRITERIA_H

#include <iostream>
#include "definitions.h"

namespace vamr_ref_criteria {
   
   class Base {
    public:
      Base();
      virtual ~Base();

      virtual Realf evaluate(const Realf* velBlock,const uint popID) = 0;
      virtual void evaluate(const Realf* velBlost,Realf* result,const uint popID);
      virtual bool initialize(const std::string& configRegion) = 0;
      
    protected:

   };

   void addRefinementCriteria();

   class RelativeDifference: public Base {
    public:
      RelativeDifference();
      ~RelativeDifference();
      
      Realf evaluate(const Realf* velBlock,const uint popID);
      void evaluate(const Realf* velBlost,Realf* result,const uint popID);
      bool initialize(const std::string& configRegion);

    protected:
      Realf df_max;
      
      Realf evaluate(const Realf& f_lef,const Realf& f_cen,const Realf& f_rgt);
   };

} // namespace vamr_ref_criteria

#endif


