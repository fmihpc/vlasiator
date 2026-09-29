#pragma once
/*
 * This file is part of Vlasiator.
 * Copyright 2010-2024 Finnish Meteorological Institute
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

#include "../definitions.h"
#include "../logger.h"
#include "../mpiconversion.h"
#include "../object_wrapper.h"
#include "../parameters.h"
#include "../readparameters.h"
#include "compression_tools.h"
#include "../spatial_cells/spatial_cell_wrapper.hpp"
#include <cstdint>

// External API for Vlasiator
namespace ASTERIX {

/*
  Compresses and reconstructs the VDFs using an Asterix Method. The original
  VDFs are overwritten.

   mpiGrid: Grid with all local spatial cells
   number_of_spatial_cells:
      Used to reduce the global comrpession achieved
   update_weights:
      If the flag is set to true the method will create a feedback loop where
  the weights of the MLP are stored and then re used for the next training
  session. This will? lead to faster convergence down the road. If the flag is
  set to false then every training session starts from randomized weights. I
  think this is what ML people call transfer learning (together with freezing
      and adding extra neuron which we do not do here).
*/
void compress_vdfs(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid, const std::vector<CellID>& local_cells,
                   P::ASTERIX_COMPRESSION_METHODS method, bool update_weights,std::vector<std::vector<char>>&mpl_bytes,uint32_t downsampling_factor=1);

/*
  Compresses the VDFs using an Asterix Method but does not overwrite them. This
  method is used to update the weights of the MLP at regular intervals without
  actually modifying the VDFs.

  mpiGrid: Grid with all local spatial cells

*/
void compress_vdfs_transfer_learning(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid);

std::size_t probe_network_size_in_bytes(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid,
                                        size_t number_of_spatial_cells);

// Function to decompress a compressed array of doubles using ZFP
std::vector<double> decompressArrayDouble(char* compressedData, size_t compressedSize, size_t arraySize, double tol);

// Function to decompress a compressed array of floats using ZFP
std::vector<float> decompressArrayFloat(char* compressedData, size_t compressedSize, size_t arraySize, float tol);

} // namespace ASTERIX
