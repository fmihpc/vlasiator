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

#include "compression_tools.h"
#include <concepts>
#include <stdexcept>
#include <sys/types.h>
#include <unordered_set>
#include <vector>

/*
Extracts VDF from spatial cell
 */
ASTERIX::UnorderedVDF ASTERIX::extract_pop_vdf_from_spatial_cell(spatial_cell::SpatialCell* sc, uint popID) {
   assert(sc && "Invalid Pointer to Spatial Cell !");
   auto blockContainer = sc->get_velocity_blocks(popID);
   const size_t total_blocks = blockContainer->size();
   const Real* max_v_lims = sc->get_population(popID).vmesh->getMeshMaxLimits();
   const Real* min_v_lims = sc->get_population(popID).vmesh->getMeshMinLimits();;
   const Real* blockParams = sc->get_block_parameters(popID);
   Realf* data = blockContainer->getData();
   assert(max_v_lims && "Invalid Pointre to max_v_limits");
   assert(min_v_lims && "Invalid Pointre to min_v_limits");
   assert(data && "Invalid Pointre block container data");
   auto vcoords = std::vector<std::array<Real, 3>>(blockContainer->size() * WID3, {Real(0), Real(0), Real(0)});
   auto vspace = std::vector<Realf>(blockContainer->size() * WID3, Realf(0));

   // xmin,ymin,zmin,xmax,ymax,zmax;
   std::array<Real, 6> vlims{std::numeric_limits<Real>::max(),    std::numeric_limits<Real>::max(),
                             std::numeric_limits<Real>::max(),    std::numeric_limits<Real>::lowest(),
                             std::numeric_limits<Real>::lowest(), std::numeric_limits<Real>::lowest()};

   std::size_t cnt = 0;
   for (std::size_t n = 0; n < total_blocks; ++n) {
      auto bp = blockParams + n * BlockParams::N_VELOCITY_BLOCK_PARAMS;
      const Realf* vdf_data = &data[n * WID3];
      for (uint k = 0; k < WID; ++k) {
         for (uint j = 0; j < WID; ++j) {
            for (uint i = 0; i < WID; ++i) {
               const Real vx = bp[BlockParams::VXCRD] + (i + 0.5) * bp[BlockParams::DVX];
               const Real vy = bp[BlockParams::VYCRD] + (j + 0.5) * bp[BlockParams::DVY];
               const Real vz = bp[BlockParams::VZCRD] + (k + 0.5) * bp[BlockParams::DVZ];
               vlims[0] = std::min(vlims[0], vx);
               vlims[1] = std::min(vlims[1], vy);
               vlims[2] = std::min(vlims[2], vz);
               vlims[3] = std::max(vlims[3], vx);
               vlims[4] = std::max(vlims[4], vy);
               vlims[5] = std::max(vlims[5], vz);
               Realf vdf_val = vdf_data[cellIndex(i, j, k)];
               vcoords[cnt] = {vx, vy, vz};
               vspace[cnt] = vdf_val;
               cnt++;
            }
         }
      }
   } // over blocks
   return UnorderedVDF{.vdf_vals = vspace, .vdf_coords = vcoords, .v_limits = vlims};
}

// Simply overwrites the VDF of this population for the give spatial cell with a
// new vspace
void ASTERIX::overwrite_pop_spatial_cell_vdf(spatial_cell::SpatialCell* sc, uint popID, const std::vector<Realf>& new_vspace) {
   assert(sc && "Invalid Pointer to Spatial Cell !");
   auto blockContainer = sc->get_velocity_blocks(popID);
   const size_t total_blocks = blockContainer->size();
   const Real* max_v_lims = sc->get_population(popID).vmesh->getMeshMaxLimits();
   const Real* min_v_lims = sc->get_population(popID).vmesh->getMeshMinLimits();
   const Real* blockParams = sc->get_block_parameters(popID);
   Realf* data = blockContainer->getData();
   assert(max_v_lims && "Invalid Pointre to max_v_limits");
   assert(min_v_lims && "Invalid Pointre to min_v_limits");
   assert(data && "Invalid Pointre block container data");

   std::size_t cnt = 0;
   for (std::size_t n = 0; n < total_blocks; ++n) {
      auto bp = blockParams + n * BlockParams::N_VELOCITY_BLOCK_PARAMS;
      Realf* vdf_data = &data[n * WID3];
      for (uint k = 0; k < WID; ++k) {
         for (uint j = 0; j < WID; ++j) {
            for (uint i = 0; i < WID; ++i) {
               vdf_data[cellIndex(i, j, k)] = new_vspace[cnt];
               cnt++;
            }
         }
      }
   } // over blocks
   return;
}

void ASTERIX::overwrite_pop_spatial_cell_vdf(spatial_cell::SpatialCell* sc, uint popID, const OrderedVDF& vdf) {
   assert(sc && "Invalid Pointer to Spatial Cell !");
   auto blockContainer = sc->get_velocity_blocks(popID);
   const size_t total_blocks = blockContainer->size();
   const Real* blockParams = sc->get_block_parameters(popID);
   Realf* data = blockContainer->getData();

   for (std::size_t n = 0; n < total_blocks; ++n) {
      auto bp = blockParams + n * BlockParams::N_VELOCITY_BLOCK_PARAMS;
      Realf* vdf_data = &data[n * WID3];
      for (uint k = 0; k < WID; ++k) {
         for (uint j = 0; j < WID; ++j) {
            for (uint i = 0; i < WID; ++i) {
               const Real dvx = (blockParams + BlockParams::N_VELOCITY_BLOCK_PARAMS)[BlockParams::DVX];
               const Real dvy = (blockParams + BlockParams::N_VELOCITY_BLOCK_PARAMS)[BlockParams::DVY];
               const Real dvz = (blockParams + BlockParams::N_VELOCITY_BLOCK_PARAMS)[BlockParams::DVZ];
               const std::size_t nx = std::ceil((vdf.v_limits[3] - vdf.v_limits[0]) / dvx);
               const std::size_t ny = std::ceil((vdf.v_limits[4] - vdf.v_limits[1]) / dvy);
               const std::size_t nz = std::ceil((vdf.v_limits[5] - vdf.v_limits[2]) / dvz);
               const Real vx = bp[BlockParams::VXCRD] + (i + 0.5) * bp[BlockParams::DVX];
               const Real vy = bp[BlockParams::VYCRD] + (j + 0.5) * bp[BlockParams::DVY];
               const Real vz = bp[BlockParams::VZCRD] + (k + 0.5) * bp[BlockParams::DVZ];
               const size_t bbox_i = std::min(static_cast<size_t>(std::floor((vx - vdf.v_limits[0]) / dvx)), nx - 1);
               const size_t bbox_j = std::min(static_cast<size_t>(std::floor((vy - vdf.v_limits[1]) / dvy)), ny - 1);
               const size_t bbox_k = std::min(static_cast<size_t>(std::floor((vz - vdf.v_limits[2]) / dvz)), nz - 1);
                    const size_t index = bbox_i * (ny * nz) + bbox_j * nz + bbox_k;
                 // vspace.at(index) += vdf_data[cellIndex(i, j, k)] / ratio;


               vdf_data[cellIndex(i, j, k)] = vdf.vdf_vals.at(index);

            }
         }
      }
   } // over blocks
   return;
}

// Extracts VDF in a cartesian C ordered mesh in a minimum BBOX and with a zoom level used for upsampling/downsampling
ASTERIX::OrderedVDF ASTERIX::extract_pop_vdf_from_spatial_cell_ordered_min_bbox_zoomed(spatial_cell::SpatialCell* sc, uint popID,
                                                                                       int zoom) {
   assert(sc && "Invalid Pointer to Spatial Cell !");
   if (zoom != 1) {
      throw std::runtime_error("Zoom is not supported yet!");
   }
   auto blockContainer = sc->get_velocity_blocks(popID);
   const size_t total_blocks = blockContainer->size();
   const Real* blockParams = sc->get_block_parameters(popID);

   // xmin,ymin,zmin,xmax,ymax,zmax;
   std::array<Real, 6> vlims{std::numeric_limits<Real>::max(),    std::numeric_limits<Real>::max(),
                             std::numeric_limits<Real>::max(),    std::numeric_limits<Real>::lowest(),
                             std::numeric_limits<Real>::lowest(), std::numeric_limits<Real>::lowest()};

   // This pass is computing the active vmesh limits
   // Store dvx,dvy,dvz here
   const Real dvx = (blockParams + BlockParams::N_VELOCITY_BLOCK_PARAMS)[BlockParams::DVX];
   const Real dvy = (blockParams + BlockParams::N_VELOCITY_BLOCK_PARAMS)[BlockParams::DVY];
   const Real dvz = (blockParams + BlockParams::N_VELOCITY_BLOCK_PARAMS)[BlockParams::DVZ];
   for (std::size_t n = 0; n < total_blocks; ++n) {
      const auto bp = blockParams + n * BlockParams::N_VELOCITY_BLOCK_PARAMS;
      for (uint k = 0; k < WID; ++k) {
         for (uint j = 0; j < WID; ++j) {
            for (uint i = 0; i < WID; ++i) {
               const Real vx = bp[BlockParams::VXCRD] + (i + 0.5) * bp[BlockParams::DVX];
               const Real vy = bp[BlockParams::VYCRD] + (j + 0.5) * bp[BlockParams::DVY];
               const Real vz = bp[BlockParams::VZCRD] + (k + 0.5) * bp[BlockParams::DVZ];
               vlims[0] = std::min(vlims[0], vx);
               vlims[1] = std::min(vlims[1], vy);
               vlims[2] = std::min(vlims[2], vz);
               vlims[3] = std::max(vlims[3], vx);
               vlims[4] = std::max(vlims[4], vy);
               vlims[5] = std::max(vlims[5], vz);
            }
         }
      }
   } // over blocks

   assert(isPow2(static_cast<size_t>(std::abs(zoom))));
   float ratio = (zoom > 0) ? static_cast<float>(std::abs(zoom)) : 1.0 / static_cast<float>(std::abs(zoom));
   assert(ratio > 0);

   const Real target_dvx = dvx * ratio;
   const Real target_dvy = dvy * ratio;
   const Real target_dvz = dvz * ratio;
   std::size_t nx = std::ceil((vlims[3] - vlims[0]) / target_dvx);
   std::size_t ny = std::ceil((vlims[4] - vlims[1]) / target_dvy);
   std::size_t nz = std::ceil((vlims[5] - vlims[2]) / target_dvz);
   // printf("VDF min box is %zu , %zu %zu \n ", nx, ny, nz);

   std::unordered_set<vmesh::GlobalID> ignore_list;
   for (std::size_t k = 0; k < nz; ++k) {
      for (std::size_t j = 0; j < ny; ++j) {
         for (std::size_t i = 0; i < nx; ++i) {
               const Real vx = vlims[0] + (i + 0.5) * dvx;
               const Real vy = vlims[1] + (j + 0.5) * dvy;
               const Real vz = vlims[2] + (k + 0.5) * dvz;
               const std::array<Real,3>coords={vx,vy,vz};
               const auto gid=sc->get_velocity_block(popID, &coords[0]);
               ignore_list.insert(gid);
         }
      }
   }

   Realf* data = blockContainer->getData();
   std::vector<Realf> vspace(nx * ny * nz, Realf(0));
   for (std::size_t n = 0; n < total_blocks; ++n) {
      const auto bp = blockParams + n * BlockParams::N_VELOCITY_BLOCK_PARAMS;
      const Realf* vdf_data = &data[n * WID3];
      const vmesh::GlobalID gid = sc->get_velocity_block_global_id(n, popID);
      (void)ignore_list.erase(gid);
      for (uint k = 0; k < WID; ++k) {
         for (uint j = 0; j < WID; ++j) {
            for (uint i = 0; i < WID; ++i) {
               const Real vx = bp[BlockParams::VXCRD] + (i + 0.5) * bp[BlockParams::DVX];
               const Real vy = bp[BlockParams::VYCRD] + (j + 0.5) * bp[BlockParams::DVY];
               const Real vz = bp[BlockParams::VZCRD] + (k + 0.5) * bp[BlockParams::DVZ];
               const size_t bbox_i = std::min(static_cast<size_t>(std::floor((vx - vlims[0]) / target_dvx)), nx - 1);
               const size_t bbox_j = std::min(static_cast<size_t>(std::floor((vy - vlims[1]) / target_dvy)), ny - 1);
               const size_t bbox_k = std::min(static_cast<size_t>(std::floor((vz - vlims[2]) / target_dvz)), nz - 1);

               // Averaging
               if (ratio >= 1.0) {
                  const size_t index = bbox_i * (ny * nz) + bbox_j * nz + bbox_k;
                  vspace.at(index) += vdf_data[cellIndex(i, j, k)] / ratio;
               } else {
                  // Same value in all bins
                  int max_off = 1 / ratio;
                  for (int off_z = 0; off_z <= max_off; off_z++) {
                     for (int off_y = 0; off_y <= max_off; off_y++) {
                        for (int off_x = 0; off_x <= max_off; off_x++) {
                           const size_t index = (bbox_i + off_x) * (ny * nz) + (bbox_j + off_y) * nz + (bbox_k + off_z);
                           if (index < vspace.size()) {
                              vspace.at(index) = vdf_data[cellIndex(i, j, k)];
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   } // over blocks

   std::vector<vmesh::GlobalID >ignored;
   ignored.reserve(ignore_list.size());
   for (auto it = ignore_list.begin(); it != ignore_list.end(); ) {
      ignored.push_back(std::move(ignore_list.extract(it++).value()));
   }

   return ASTERIX::OrderedVDF{.blocks_to_ignore=ignored,.sparse_vdf_bytes=total_blocks*WID*WID*WID*sizeof(Realf),.vdf_vals = vspace, .v_limits = vlims, .shape = {nx, ny, nz}};
}

void ASTERIX::overwrite_cellids_vdfs(const std::span<const CellID> cids, uint popID,
                                     dccrg::Dccrg<spatial_cell::SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid,
                                     const std::vector<std::array<Real, 3>>& vcoords,
                                     const std::vector<Realf>& vspace_union,
                                     const std::unordered_map<vmesh::LocalID, std::size_t>& map_exists_id) {
   const std::size_t nrows = vcoords.size();
   const std::size_t ncols = cids.size();
   // This will be used further down for indexing into the vspace_union
   auto index_2d = [nrows, ncols](std::size_t row, std::size_t col) -> std::size_t { return row * ncols + col; };

   for (std::size_t cc = 0; cc < cids.size(); ++cc) {
      const auto& cid = cids[cc];
      spatial_cell::SpatialCell* sc = mpiGrid[cid];
      auto blockContainer = sc->get_velocity_blocks(popID);
      const size_t total_blocks = blockContainer->size();
      Realf* data = blockContainer->getData();
      const Real* blockParams = sc->get_block_parameters(popID);
      for (std::size_t n = 0; n < total_blocks; ++n) {
         const auto bp = blockParams + n * BlockParams::N_VELOCITY_BLOCK_PARAMS;
         const vmesh::GlobalID gid = sc->get_velocity_block_global_id(n, popID);
         const auto it = map_exists_id.find(gid);
         const bool exists = it != map_exists_id.end();
         if (!exists){
            std::cerr<<"This should not happen!"<<std::endl;
            abort();
         }
         assert(exists && "Someone has a buuuug!");
         const auto index = it->second;
         Realf* vdf_data = &data[n * WID3];
         size_t cnt = 0;
         for (uint k = 0; k < WID; ++k) {
            for (uint j = 0; j < WID; ++j) {
               for (uint i = 0; i < WID; ++i) {
                  const std::size_t index = it->second;
                  vdf_data[cellIndex(i, j, k)] = vspace_union[index_2d(index + cnt, cc)];
                  cnt++;
               }
            }
         }
      }
   }
   return;
}


void ASTERIX::dump_vdf_to_binary_file(const char* filename, CellID cid,
                                      dccrg::Dccrg<spatial_cell::SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid) {

   spatial_cell::SpatialCell* sc = mpiGrid[cid];
   assert(sc && "Invalid Pointer to Spatial Cell !");
   OrderedVDF vdf = extract_pop_vdf_from_spatial_cell_ordered_min_bbox_zoomed(sc, 0, 1);
   vdf.save_to_file(filename);
}

// Taken from DataReducers
Real ASTERIX::get_Non_MaxWellianity(const spatial_cell::SpatialCell* cell, uint popID) {
   Real rho;
   Real V0[3];
   Real b_par[3];
   Real b_perp1[3];
   Real b_perp2[3];
   Real T_par;
   Real T_perp;
   Real epsilon;
   // calculate here rho, v, T
   epsilon = 0.0;

   // get rho and bulk speed
   rho = cell->get_population(popID).RHO;
   V0[0] = cell->get_population(popID).V[0];
   V0[1] = cell->get_population(popID).V[1];
   V0[2] = cell->get_population(popID).V[2];

   // calculate temperature from the pressure tensor
   Real PTensor[3] = {};

   // parallel unit vector (B)
   Real BX = cell->parameters[CellParams::PERBXVOL] + cell->parameters[CellParams::BGBXVOL];
   Real BY = cell->parameters[CellParams::PERBYVOL] + cell->parameters[CellParams::BGBYVOL];
   Real BZ = cell->parameters[CellParams::PERBZVOL] + cell->parameters[CellParams::BGBZVOL];
   Real norm_par = sqrt(BX * BX + BY * BY + BZ * BZ);
   b_par[0] = BX / norm_par;
   b_par[1] = BY / norm_par;
   b_par[2] = BZ / norm_par;

   // perpendicular unit vector 1 (bulk velocity perpendicular to b)
   Real BV0 = sqrt(b_par[0] * V0[0] + b_par[1] * V0[1] + b_par[2] * V0[2]);
   b_perp1[0] = V0[0] - BV0 * b_par[0];
   b_perp1[1] = V0[1] - BV0 * b_par[1];
   b_perp1[2] = V0[2] - BV0 * b_par[2];
   Real norm_perp1 = sqrt(b_perp1[0] * b_perp1[0] + b_perp1[1] * b_perp1[1] + b_perp1[2] * b_perp1[2]);
   if (!(norm_perp1 > 0.0)) {
      // if V0 is aligned with b, take arbitrary perpendicular vector
      b_perp1[0] = +b_par[1] + b_par[2];
      b_perp1[1] = +b_par[2] - b_par[0];
      b_perp1[2] = -b_par[0] - b_par[1];
      norm_perp1 = sqrt(b_perp1[0] * b_perp1[0] + b_perp1[1] * b_perp1[1] + b_perp1[2] * b_perp1[2]);
   }
   b_perp1[0] /= norm_perp1;
   b_perp1[1] /= norm_perp1;
   b_perp1[2] /= norm_perp1;

   // perpendicular unit vector 2 (b_par x b_perp1)
   b_perp2[0] = b_par[1] * b_perp1[2] - b_par[2] * b_perp1[1];
   b_perp2[1] = b_par[2] * b_perp1[0] - b_par[0] * b_perp1[2];
   b_perp2[2] = b_par[0] * b_perp1[1] - b_par[1] * b_perp1[0];
   Real norm_perp2 = sqrt(b_perp2[0] * b_perp2[0] + b_perp2[1] * b_perp2[1] + b_perp2[2] * b_perp2[2]);
   b_perp2[0] /= norm_perp2;
   b_perp2[1] /= norm_perp2;
   b_perp2[2] /= norm_perp2;

   // below calculation is modified from VariablePTensorDiagonal
   constexpr Real HALF = 0.5;
#pragma omp parallel
   {
      Real thread_nvxvx_sum = 0.0;
      Real thread_nvyvy_sum = 0.0;
      Real thread_nvzvz_sum = 0.0;

      const Real* parameters = cell->get_block_parameters(popID);
      const Realf* block_data = cell->get_data(popID);

#pragma omp for
      for (vmesh::LocalID n = 0; n < cell->get_number_of_velocity_blocks(popID); n++) {
         for (uint k = 0; k < WID; ++k)
            for (uint j = 0; j < WID; ++j)
               for (uint i = 0; i < WID; ++i) {
                  const Real VX = parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VXCRD] +
                                  (i + HALF) * parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVX];
                  const Real VY = parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VYCRD] +
                                  (j + HALF) * parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVY];
                  const Real VZ = parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VZCRD] +
                                  (k + HALF) * parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVZ];
                  const Real DV3 = parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVX] *
                                   parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVY] *
                                   parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVZ];

                  const Real V_par = (VX - V0[0]) * b_par[0] + (VY - V0[1]) * b_par[1] + (VZ - V0[2]) * b_par[2];
                  const Real V_perp1 =
                      (VX - V0[0]) * b_perp1[0] + (VY - V0[1]) * b_perp1[1] + (VZ - V0[2]) * b_perp1[2];
                  const Real V_perp2 =
                      (VX - V0[0]) * b_perp2[0] + (VY - V0[1]) * b_perp2[1] + (VZ - V0[2]) * b_perp2[2];

                  thread_nvxvx_sum += block_data[n * SIZE_VELBLOCK + cellIndex(i, j, k)] * V_par * V_par * DV3;
                  thread_nvyvy_sum += block_data[n * SIZE_VELBLOCK + cellIndex(i, j, k)] * V_perp1 * V_perp1 * DV3;
                  thread_nvzvz_sum += block_data[n * SIZE_VELBLOCK + cellIndex(i, j, k)] * V_perp2 * V_perp2 * DV3;
               }
      }
      thread_nvxvx_sum *= getObjectWrapper().particleSpecies[popID].mass;
      thread_nvyvy_sum *= getObjectWrapper().particleSpecies[popID].mass;
      thread_nvzvz_sum *= getObjectWrapper().particleSpecies[popID].mass;

      // Accumulate contributions coming from this velocity block to the
      // spatial cell velocity moments. If multithreading / OpenMP is used,
      // these updates need to be atomic:
#pragma omp critical
      {
         PTensor[0] += thread_nvxvx_sum;
         PTensor[1] += thread_nvyvy_sum;
         PTensor[2] += thread_nvzvz_sum;
      }
   }
   T_par = (PTensor[0]) / (rho * physicalconstants::K_B);
   T_perp = (PTensor[1] + PTensor[2]) / (2.0 * rho * physicalconstants::K_B);

   // thermal speed in parallel direction
   const Real V_par_th_sq = 2.0 * physicalconstants::K_B * T_par / getObjectWrapper().particleSpecies[popID].mass;

#pragma omp parallel
   {
      Real thread_epsilon_sum = 0.0;

      const Real* parameters = cell->get_block_parameters(popID);
      const Realf* block_data = cell->get_data(popID);

#pragma omp for
      for (vmesh::LocalID n = 0; n < cell->get_number_of_velocity_blocks(popID); n++) {
         for (uint k = 0; k < WID; ++k)
            for (uint j = 0; j < WID; ++j)
               for (uint i = 0; i < WID; ++i) {
                  const Real VX = parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VXCRD] +
                                  (i + HALF) * parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVX];
                  const Real VY = parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VYCRD] +
                                  (j + HALF) * parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVY];
                  const Real VZ = parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::VZCRD] +
                                  (k + HALF) * parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVZ];
                  const Real DV3 = parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVX] *
                                   parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVY] *
                                   parameters[n * BlockParams::N_VELOCITY_BLOCK_PARAMS + BlockParams::DVZ];

                  const Real V_par = (VX - V0[0]) * b_par[0] + (VY - V0[1]) * b_par[1] + (VZ - V0[2]) * b_par[2];
                  const Real V_perp1 =
                      (VX - V0[0]) * b_perp1[0] + (VY - V0[1]) * b_perp1[1] + (VZ - V0[2]) * b_perp1[2];
                  const Real V_perp2 =
                      (VX - V0[0]) * b_perp2[0] + (VY - V0[1]) * b_perp2[1] + (VZ - V0[2]) * b_perp2[2];

                  const Real bimaxwellian =
                      rho / sqrt(M_PI * M_PI * M_PI * V_par_th_sq * V_par_th_sq * V_par_th_sq) * (T_par / T_perp) *
                      exp(-(V_par * V_par) / V_par_th_sq -
                          (V_perp1 * V_perp1 + V_perp2 * V_perp2) / (V_par_th_sq * T_perp / T_par));

                  thread_epsilon_sum +=
                      (abs(block_data[n * SIZE_VELBLOCK + cellIndex(i, j, k)] - bimaxwellian) - bimaxwellian) * DV3;
               }
      }

      // Accumulate contributions coming from this velocity block to the
      // spatial cell velocity moments. If multithreading / OpenMP is used,
      // these updates need to be atomic:
#pragma omp critical
      { epsilon += thread_epsilon_sum; }
   }
   epsilon *= HALF / rho;
   epsilon += HALF;
   return epsilon;
}
