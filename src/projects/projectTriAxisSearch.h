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
#ifndef TRIAXISSEARCH_H
#define TRIAXISSEARCH_H

#include "project.h"

namespace projects {
   class TriAxisSearch: public Project {
    public:
         
    protected:
        /*! \brief Find blocks above the threshold centred isotropically around a bulk velocity.
         * 
         * Instead of looping through the whole velocity space this function starts from the project's bulk velocity V0[3].
         * It then proceeds along V[XYZ] successively to determine at what maximum radius a block falls below (0.1 times) the threshold.
         * 
         * This radius is used to determine all blocks within that radius of V0, create them and return their list for initialisation.
         */
        virtual std::vector<vmesh::GlobalID> findBlocksToInitialize(spatial_cell::SpatialCell* cell,const uint popID) const;
      
      /*! \brief Return a vector containing the velocity coordinate of the centre of each ion population in the distribution.
       * 
       * This function is used by findBlocksToInitialize() to start the search of the extent of the distribution along each axis.
       * The extension to a vector allows to have more than one population in each spatial cell.
       * 
       * \sa findBlocksToInitialize
       */
      virtual std::vector<std::array<Real, 3>> getV0(
                                                     creal x,
                                                     creal y,
                                                     creal z,
                                                     const uint popID
                                                    ) const = 0;
   }; // class TriAxisSearch
} // namespace


#endif
