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

#include "compression.h"
#include "compression_tools.h"
#include "zfp/array1.hpp"
#include <atomic>
#include <concepts>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <ranges>
#include <stdexcept>
#include <type_traits>
#include <unordered_map>
#include <vector>
#include <zfp.h>
#include <omp.h>

extern Logger logFile;
#define ASTERIX_USE_GPU
#define MEMPOOL_BYTES 60ul*1024ul*1024ul*1024ul 
using namespace ASTERIX;
using namespace spatial_cell;

#ifdef ASTERIX_MLP
#ifdef __cplusplus
extern "C" {
#endif
size_t compress_phasespace6D_f64(GENERIC_TS_POOL::MemPool* p, std::size_t fin,std::size_t fout, double* coords_ptr, double* f_ptr,
                                 std::size_t size, std::size_t max_epochs, std::size_t fourier_order,
                                 size_t* hidden_layers_ptr, size_t n_hidden_layers, double sparsity, double tol,
                                 double* weights_ptr, std::size_t weight_size, bool use_input_weights,
                                 uint32_t downsampling_factor, double& error, uint32_t& epochs, int& status, int rankID);

size_t compress_phasespace6D_f32(GENERIC_TS_POOL::MemPool* p, std::size_t fin,std::size_t fout, float* coords_ptr, float* f_ptr,
                                 std::size_t size, std::size_t max_epochs, std::size_t fourier_order,
                                 size_t* hidden_layers_ptr, size_t n_hidden_layers, float sparsity, float tol,
                                 float* weights_ptr, std::size_t weight_size, bool use_input_weights,
                                 uint32_t downsampling_factor, float& error, uint32_t& epochs, int& status, int rankID);


#ifdef __cplusplus
}
#endif
auto compress_vdfs_fourier_mlp(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid,
                               size_t number_of_spatial_cells, bool update_weights, std::vector<std::vector<char>>&bytes ,uint32_t downsampling_factor)
    -> float;

auto compress_vdfs_fourier_mlp_clustered(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid,
                                         size_t number_of_spatial_cells, bool update_weights, std::vector<std::vector<char>>&bytes,
                                         uint32_t downsampling_factor) -> float;
#endif //ASTERIX_MLP

#ifdef ASTERIX_ZFP


auto compress_vdfs_zfp(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid, const std::vector<CellID>& local_cells)
    -> float;

auto compress(float* array, size_t arraySize, size_t& compressedSize, float tol) -> std::vector<char>;

auto compress(double* array, size_t arraySize, size_t& compressedSize, double tol) -> std::vector<char>;

auto decompressArrayDouble(char* compressedData, size_t compressedSize, size_t arraySize, double tol)
    -> std::vector<double>;

auto decompressArrayFloat(char* compressedData, size_t compressedSize, size_t arraySize, float tol)
    -> std::vector<float>;

#endif //ASTERIX_ZFP

#ifdef ASTERIX_OCTREE
auto compress_vdfs_octree(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid, const std::vector<CellID>& local_cells)
    -> float;
#endif

// Main driver, look at header file  for documentation
void ASTERIX::compress_vdfs(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid, const std::vector<CellID>& cells,
                            P::ASTERIX_COMPRESSION_METHODS method, bool update_weights,
                            std::vector<std::vector<char>>& bytes, uint32_t downsampling_factor /*=1*/) {

   if (downsampling_factor < 1) {
      throw std::runtime_error("Requested downsampling factor in VDF compression makes no sense!");
   }
   float local_compression_ratio = 0.0;
   switch (method) {
#ifdef ASTERIX_MLP
   case P::ASTERIX_COMPRESSION_METHODS::MLP:
      local_compression_ratio =
          compress_vdfs_fourier_mlp(mpiGrid, cells.size(), update_weights, bytes,downsampling_factor);
      break;
   case P::ASTERIX_COMPRESSION_METHODS::MLP_MULTI:
      local_compression_ratio =
          compress_vdfs_fourier_mlp_clustered(mpiGrid, cells.size(), update_weights, bytes, downsampling_factor);
      break;
#endif
#ifdef ASTERIX_ZFP
   case P::ASTERIX_COMPRESSION_METHODS::ZFP:
      local_compression_ratio = compress_vdfs_zfp(mpiGrid, cells);
      break;
#endif
#ifdef ASTERIX_OCTREE
   case P::ASTERIX_COMPRESSION_METHODS::OCTREE:
      local_compression_ratio = compress_vdfs_octree(mpiGrid, cells);
      break;
#endif
   case P::ASTERIX_COMPRESSION_METHODS::NONE:
      break;
   default:
      throw std::runtime_error("This is bad!. Improper Asterix method detected!");
      break;
   };

   // Reduce global compression ratio
   int myRank;
   int mpiProcs;
   float global_compression_ratio = 0.0;
   MPI_Comm_size(MPI_COMM_WORLD, &mpiProcs);
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Reduce(&local_compression_ratio, &global_compression_ratio, 1, MPI_FLOAT, MPI_SUM, MASTER_RANK, MPI_COMM_WORLD);
   MPI_Barrier(MPI_COMM_WORLD);

   if (myRank == MASTER_RANK) {
      logFile << "(VDF COMPRESSION INFO): Compression Ratio = "
              << global_compression_ratio / static_cast<float>(mpiProcs) << std::endl;
   }
}

std::vector<std::pair<std::size_t, std::size_t>> partition(std::size_t array_size, std::size_t chunk_size,
                                                           std::size_t max_chunks) {
   std::vector<std::pair<std::size_t, std::size_t>> result;
   if (array_size == 0 || chunk_size <= 0 || max_chunks <= 0) {
      throw std::runtime_error("ERROR: catastrophic failure in VDF partitioning.");
   }
   std::size_t optimal_chunks = std::max(1.0,std::ceil(array_size / chunk_size));
   std::size_t numChunks = std::min(optimal_chunks, max_chunks);
   std::size_t baseSize = array_size / numChunks;
   std::size_t largerChunks = array_size % numChunks;
   std::size_t start = 0;
   for (std::size_t i = 0; i < numChunks; ++i) {
      std::size_t chunkSize = baseSize + (i < largerChunks ? (1) : (0));
      result.push_back({start, start + chunkSize});
      start += chunkSize;
   }
   return result;
}

std::vector<CellID> sort_cells_based_on_maxwellianity(const std::vector<CellID>& local_cells, uint popID,
                                                      dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid) {
   std::vector<std::pair<CellID, Real>> sorted_vdf;
   sorted_vdf.reserve(local_cells.size());
   for (const auto& cid : local_cells) {
      sorted_vdf.emplace_back(cid, get_Non_MaxWellianity(mpiGrid[cid], popID));
   }
   std::ranges::sort(sorted_vdf, {}, &std::pair<CellID, Real>::second);
   std::vector<CellID> sorted_cells;
   sorted_cells.reserve(sorted_vdf.size());
   for (const auto& [cid, _] : sorted_vdf) {
      sorted_cells.push_back(cid);
   }
   return sorted_cells;
}

#ifdef ASTERIX_MLP
float compress_vdfs_fourier_mlp(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid,
                                size_t number_of_spatial_cells, bool update_weights,
                                std::vector<std::vector<char>>& bytes, uint32_t downsampling_factor) {

   //Grab the rank TODO: remove this later
   GENERIC_TS_POOL::MemPool p{};
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

   if (getObjectWrapper().particleSpecies.size() > 1) {
      throw std::runtime_error("Multi-Pop not implemented yet!");
   }
   std::atomic<float> local_compression_achieved = 0.0;
   std::size_t total_samples = 0;
   for (uint popID = 0; popID < getObjectWrapper().particleSpecies.size(); ++popID) {

      Real sparse = getObjectWrapper().particleSpecies[popID].sparseMinValue;
      const std::vector<CellID>& _local_cells = getLocalCells();
      std::vector<CellID> local_cells;
      for (auto& c : _local_cells) {
         if (mpiGrid[c]->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
            continue;
         }
         auto blockContainer = mpiGrid[c]->get_velocity_blocks(popID);
         const size_t total_blocks = blockContainer->size();
         if (total_blocks==0){
            continue;
         }
         local_cells.push_back(c);
      }
      // local_cells=sort_cells_based_on_maxwellianity(local_cells, popID, mpiGrid);
      const std::size_t num_threads = omp_get_max_threads();
      const auto partitionScheme = partition(local_cells.size(), P::max_vdfs_per_nn, num_threads);

      // printf("Compression on %zu threads\n",num_threads);
      // for (std::size_t ii=0;ii<partitionScheme.size();++ii){
      //    printf("\tChunk(%zu): [%zu,%zu) ->len=(%zu)\n",ii, partitionScheme[ii].first,partitionScheme[ii].second,partitionScheme[ii].second-partitionScheme[ii].first  );
      // }
      const std::size_t threads_needed = partitionScheme.size();
      total_samples += partitionScheme.size();
      omp_set_num_threads(threads_needed);
      std::vector<std::vector<char>> thread_bytes(threads_needed);
#pragma omp parallel
      {
         std::size_t thread_id = omp_get_thread_num();
         const std::pair<std::size_t, std::size_t> index_range = partitionScheme.at(thread_id);
         const std::size_t count = index_range.second - index_range.first;
         printf("Count = %zu \n",count);
         const auto start = local_cells.data() + index_range.first;
         const std::span<const CellID> span(start, count);
         PhaseSpaceUnion<Realf> b(span, popID, mpiGrid, true);
         b.normalize();

         // (2) Do the compression for this VDF
         Realf error = std::numeric_limits<double>::max();
         uint32_t epochs=0;
         int status = 0;
         // Allocate spaced for weights
         auto network_size = calculate_total_size_bytes<Realf>(P::mlp_arch, P::mlp_fourier_order, b._cids.size());
         b._network_weights = (Realf*)malloc(network_size);
         b._n_weights = network_size / sizeof(Realf);

#ifndef DPF
         std::size_t nn_mem_footprint_bytes = compress_phasespace6D_f32(
             &p, 3, span.size(), &b._vcoords[0][0], b._vspace.data(), b._vcoords.size(), P::mlp_max_epochs,
             P::mlp_fourier_order, P::mlp_arch.data(), P::mlp_arch.size(), sparse, P::mlp_tollerance,
             b._network_weights, network_size, false, downsampling_factor, error, epochs ,status,myRank);

#else
         std::size_t nn_mem_footprint_bytes = compress_phasespace6D_f64(
             &p, 3, span.size(), &b._vcoords[0][0], b._vspace.data(), b._vcoords.size(), P::mlp_max_epochs,
             P::mlp_fourier_order, P::mlp_arch.data(), P::mlp_arch.size(), sparse, P::mlp_tollerance,
             b._network_weights, network_size, false, downsampling_factor, error, epochs, status, myRank);
#endif
         for (const auto sc:span){
            mpiGrid[sc]->get_population(popID).mlp_error=static_cast<float>(error);
            mpiGrid[sc]->get_population(popID).mlp_epochs = epochs;
         } 
         assert(network_size == nn_mem_footprint_bytes && "Mismatch betweeen estimated and actual network size!!!");
         thread_bytes.at(thread_id) = std::vector<char>(b.total_serialized_size_bytes());
         b.serialize_into(reinterpret_cast<unsigned char*>(thread_bytes.at(thread_id).data()));

         free(b._network_weights);
         local_compression_achieved +=
             static_cast<float>(b._effective_vdf_size) / static_cast<float>(nn_mem_footprint_bytes);
      }
      bytes = thread_bytes;
   }
   return local_compression_achieved / static_cast<float>(total_samples);
}

std::vector<std::vector<std::pair<CellID, Real>>>
clusterVDFs(const std::vector<CellID>& local_cells, const dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid,
            uint popID) {
   std::vector<Real> non_maxwellianity(local_cells.size(), 0.0);
   std::transform(local_cells.begin(), local_cells.end(), non_maxwellianity.begin(),
                  [&](const auto& cid) { return get_Non_MaxWellianity(mpiGrid[cid], popID); });
   std::vector<std::pair<CellID, Real>> sorted_vdf(local_cells.size());
   for (std::size_t i = 0; i < local_cells.size(); ++i) {
      sorted_vdf[i] = std::pair<CellID, Real>{local_cells[i], non_maxwellianity[i]};
   }
   std::sort(sorted_vdf.begin(), sorted_vdf.end(),
             [=](std::pair<CellID, Real>& a, std::pair<CellID, Real>& b) { return a.second < b.second; });
   std::vector<std::vector<std::pair<CellID, Real>>> clusters;
   std::vector<std::pair<CellID, Real>> current_cluster;
   for (const auto& pair : sorted_vdf) {
      if (current_cluster.empty()) {
         current_cluster.push_back(pair);
      } else {
         Real last_value = current_cluster.back().second;
         Real margin = 0.2f * std::max(last_value, pair.second);
         if (std::fabs(last_value - pair.second) <= margin) {
            current_cluster.push_back(pair);
         } else {
            clusters.push_back(current_cluster);
            current_cluster.clear();
            current_cluster.push_back(pair);
         }
      }
   }
   if (!current_cluster.empty()) {
      clusters.push_back(current_cluster);
   }
   return clusters;
}

float compress_vdfs_fourier_mlp_clustered(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid,
                                          size_t number_of_spatial_cells, bool update_weights, std::vector<std::vector<char>>&bytes,
                                          uint32_t downsampling_factor) {
   
   //Memory allocation
   //TODO remove later
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
   GENERIC_TS_POOL::MemPool p{};
   if (getObjectWrapper().particleSpecies.size() > 1) {
      throw std::runtime_error("Multi-Pop not implemented yet!");
   }
   float local_compression_achieved = 0.0;
   std::size_t total_samples = 0;
   const std::size_t max_span_size = P::max_vdfs_per_nn;
   for (uint popID = 0; popID < getObjectWrapper().particleSpecies.size(); ++popID) {
      Real sparse = getObjectWrapper().particleSpecies[popID].sparseMinValue;
      const std::vector<CellID>& _local_cells = getLocalCells();
      std::vector<CellID> local_cells;
      for (auto& c : _local_cells) {
         if (mpiGrid[c]->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
            continue;
         }
         local_cells.push_back(c);
      }
      const auto clusters = clusterVDFs(local_cells, mpiGrid, popID);
      std::cout << "Generated " << clusters.size() << " clusters" << std::endl;

      bytes.resize(clusters.size());
#pragma omp parallel for reduction(+ : local_compression_achieved)
      for (std::size_t i =0 ;i< clusters.size();++i) {
         auto& cluster = clusters.at(i);
#pragma omp atomic
         total_samples++;

         std::vector<CellID> cids(cluster.size());
         std::transform(cluster.begin(), cluster.end(), cids.begin(), [](const auto& pair) { return pair.first; });

         // Extract this span of VDFs as a union
         const std::span<const CellID> span(cids.data(), cids.size());
         PhaseSpaceUnion<Realf>b(span, popID, mpiGrid, true);
         b.normalize();

         // (2) Do the compression for this VDF
         Realf error = std::numeric_limits<double>::max();
         uint32_t epochs=0;
         int status = 0;
         // Allocate spaced for weights
         auto network_size =
             calculate_total_size_bytes<Realf>(P::mlp_arch, P::mlp_fourier_order, b._cids.size());
         b._network_weights = (Realf*)malloc(network_size);
         b._n_weights = network_size / sizeof(Realf);

         #ifndef DPF
            std::size_t nn_mem_footprint_bytes = compress_phasespace6D_f32(
                &p, 3, span.size(), &b._vcoords[0][0], b._vspace.data(), b._vcoords.size(), P::mlp_max_epochs,
                P::mlp_fourier_order, P::mlp_arch.data(), P::mlp_arch.size(), sparse, P::mlp_tollerance,
                b._network_weights, network_size, false, downsampling_factor, error, epochs, status,myRank);

         #else
            std::size_t nn_mem_footprint_bytes = compress_phasespace6D_f64(
                &p, 3, span.size(), &b._vcoords[0][0], b._vspace.data(), b._vcoords.size(), P::mlp_max_epochs,
                P::mlp_fourier_order, P::mlp_arch.data(), P::mlp_arch.size(), sparse, P::mlp_tollerance,
                b._network_weights, network_size, false, downsampling_factor, error, epochs, status,myRank);
         #endif
         for (const auto sc:span){
            mpiGrid[sc]->get_population(popID).mlp_error=static_cast<float>(error);
            mpiGrid[sc]->get_population(popID).mlp_epochs = epochs;
         } 
         assert(network_size == nn_mem_footprint_bytes && "Mismatch betweeen estimated and actual network size!!!");
         
         bytes.at(i).resize(b.total_serialized_size_bytes());
         b.serialize_into(reinterpret_cast<unsigned char*>(bytes.at(i).data()));
         free(b._network_weights);
         local_compression_achieved += static_cast<float>(b._effective_vdf_size) / static_cast<float>(nn_mem_footprint_bytes);
      }
   } // loop over all populations
   return local_compression_achieved / static_cast<float>(total_samples);
}
#endif //ASTERIX_MLP

#ifdef ASTERIX_ZFP
// Compresses and reconstucts VDFs using ZFP
float compress_vdfs_zfp(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid, const std::vector<CellID>& local_cells) {
   float local_compression_achieved = 0.0;
   std::size_t total_samples = 0;
   for (uint popID = 0; popID < getObjectWrapper().particleSpecies.size(); ++popID) {
      Real sparse = getObjectWrapper().particleSpecies[popID].sparseMinValue;
      // Vlasiator boilerplate
#pragma omp parallel for reduction(+ : local_compression_achieved)
      for (auto& cid : local_cells) { // loop over spatial cells
         SpatialCell* sc = mpiGrid[cid];
         assert(sc && "Invalid Pointer to Spatial Cell !");
         if (sc->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
            continue;
         }
         auto blockContainer = sc->get_velocity_blocks(popID);
         const size_t total_blocks = blockContainer->size();
         if (total_blocks==0){
            sc->get_population(popID).compressed_state_buffer = {};
            continue;
         }

#pragma omp atomic
         total_samples++;
         // (1) Extract and Collect the VDF of this cell
         UnorderedVDF vdf = extract_pop_vdf_from_spatial_cell(sc, popID);
         // (2) Do the compression for this VDF
         // Create spave for the reconstructed VDF
         size_t ss{0};
         sc->get_population(popID).compressed_state_buffer =
             compress(vdf.vdf_vals.data(), vdf.vdf_vals.size(), ss, sparse);
         float ratio = static_cast<float>(vdf.vdf_vals.size() * sizeof(Realf)) / static_cast<float>(ss);
         local_compression_achieved += ratio;
      } // loop over all spatial cells
   }    // loop over all populations
   return local_compression_achieved / static_cast<float>(total_samples);
}

#endif //ASTERIX_ZFP

#ifdef ASTERIX_OCTREE
// Compresses and reconstucts VDFs using ZFP
float compress_vdfs_octree(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid,
                           const std::vector<CellID>& local_cells) {
   int total_bytes = 0;
   int global_total_bytes = 0;
   float local_compression_achieved = 0.0;
   std::size_t total_samples = 0;
   for (uint popID = 0; popID < getObjectWrapper().particleSpecies.size(); ++popID) {
      // Vlasiator boilerplate
#pragma omp parallel for reduction(+ : total_bytes, local_compression_achieved)
      for (auto& cid : local_cells) { // loop over spatial cells
         SpatialCell* sc = mpiGrid[cid];
         if (sc->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
            continue;
         }
         auto blockContainer = sc->get_velocity_blocks(popID);
         const size_t total_blocks = blockContainer->size();
         if (total_blocks==0){
            sc->get_population(popID).compressed_state_buffer = {};
            continue;
         }
         // (1) Extract and Collect the VDF of this cell
         OrderedVDF vdf = extract_pop_vdf_from_spatial_cell_ordered_min_bbox_zoomed(sc, popID, 1);

#pragma omp atomic
         total_samples++;
         // (2) Do the compression for this VDF
         uint8_t* bytes = nullptr;
         std::size_t n_bytes;
         constexpr std::size_t maxiter = 50000;
         constexpr std::size_t skip_levels = 4;
         int status = compress_with_toctree_method(vdf.vdf_vals.data(), vdf.shape[0], vdf.shape[1], vdf.shape[2],
                                                   P::octree_tolerance, &bytes, &n_bytes, maxiter,skip_levels);

         switch(status) {
           case TOCTREE_COMPRESS_STAT_SUCCESS:
             break;
           case TOCTREE_COMPRESS_STAT_FAIL_TOL:
             logFile << "(VDF COMPRESSION INFO): T-Octree failed to reach tolerance " << 
               P::octree_tolerance << " in " << maxiter << " iterations (cid " << cid <<")\n";
             break;
           default:
             throw std::runtime_error("(VDF COMPRESSION ERROR): T-Octree failed.");
             break;
         }

         //Copy compressed state to SC
         sc->get_population(popID).compressed_state_buffer.resize(n_bytes+sizeof(std::size_t) +vdf.blocks_to_ignore.size()*sizeof(vmesh::GlobalID)+3*sizeof(std::size_t)+6*sizeof(Real),0);
         
         std::size_t ignored_blocks=vdf.blocks_to_ignore.size();
         std::size_t write_index=0;
         std::memcpy(sc->get_population(popID).compressed_state_buffer.data()+write_index,&ignored_blocks,sizeof(std::size_t));
         write_index+=sizeof(std::size_t);
         std::memcpy(sc->get_population(popID).compressed_state_buffer.data()+write_index,vdf.blocks_to_ignore.data(),ignored_blocks*sizeof(vmesh::GlobalID));
         write_index+=ignored_blocks*sizeof(vmesh::GlobalID);
         std::memcpy(sc->get_population(popID).compressed_state_buffer.data()+write_index,&vdf.shape[0],3*sizeof(std::size_t));
         write_index+=3*sizeof(std::size_t);
         std::memcpy(&sc->get_population(popID).compressed_state_buffer[write_index],&vdf.v_limits,6*sizeof(Real));
         write_index+=6*sizeof(Real);
         std::memcpy(&sc->get_population(popID).compressed_state_buffer[write_index],bytes,n_bytes);

         if (bytes != NULL) {
            free(bytes);
         }
         total_bytes += n_bytes;
         local_compression_achieved += vdf.sparse_vdf_bytes / static_cast<float>(n_bytes);

      } // loop over all spatial cells
   }    // loop over all populations
   return local_compression_achieved / static_cast<float>(total_samples);
}
#endif //ASTERIX_OCTREE

#ifdef ASTERIX_ZFP
std::vector<char> compress(float* array, size_t arraySize, size_t& compressedSize, float tol) {
   // Allocate memory for compressed data
   zfp_stream* zfp = zfp_stream_open(NULL);
   zfp_field* field = zfp_field_1d(array, zfp_type_float, arraySize);
   size_t maxSize = zfp_stream_maximum_size(zfp, field);
   std::vector<char> compressedData(maxSize);

   // Initialize ZFP compression
   zfp_stream_set_accuracy(zfp, tol);
   bitstream* stream = stream_open(compressedData.data(), compressedSize);
   zfp_stream_set_bit_stream(zfp, stream);
   zfp_stream_rewind(zfp);

   // Compress the array
   compressedSize = zfp_compress(zfp, field);
   compressedData.erase(compressedData.begin() + compressedSize, compressedData.end());
   zfp_field_free(field);
   zfp_stream_close(zfp);
   stream_close(stream);
   return compressedData;
}

// Function to decompress a compressed array of floats using ZFP
std::vector<float> ASTERIX::decompressArrayFloat(char* compressedData, size_t compressedSize, size_t arraySize,
                                                 float tol) {
   // Allocate memory for decompresseFloatd data
   std::vector<float> decompressedArray(arraySize);

   // Initialize ZFP decompression
   zfp_stream* zfp = zfp_stream_open(NULL);
   zfp_stream_set_accuracy(zfp, tol);
   bitstream* stream_decompress = stream_open(compressedData, compressedSize);
   zfp_stream_set_bit_stream(zfp, stream_decompress);
   zfp_stream_rewind(zfp);

   // Decompress the array
   zfp_field* field_decompress = zfp_field_1d(decompressedArray.data(), zfp_type_float, decompressedArray.size());
   size_t retval = zfp_decompress(zfp, field_decompress);
   (void)retval;
   zfp_field_free(field_decompress);
   zfp_stream_close(zfp);
   stream_close(stream_decompress);

   return decompressedArray;
}

// Function to compress a 1D array of doubles using ZFP
std::vector<char> compress(double* array, size_t arraySize, size_t& compressedSize, double tol) {
   zfp_stream* zfp = zfp_stream_open(NULL);
   zfp_field* field = zfp_field_1d(array, zfp_type_double, arraySize);
   size_t maxSize = zfp_stream_maximum_size(zfp, field);
   std::vector<char> compressedData(maxSize);

   zfp_stream_set_accuracy(zfp, tol);
   bitstream* stream = stream_open(compressedData.data(), compressedSize);
   zfp_stream_set_bit_stream(zfp, stream);
   zfp_stream_rewind(zfp);

   compressedSize = zfp_compress(zfp, field);
   compressedData.erase(compressedData.begin() + compressedSize, compressedData.end());
   zfp_field_free(field);
   zfp_stream_close(zfp);
   stream_close(stream);
   return compressedData;
}

// Function to decompress a compressed array of doubles using ZFP
std::vector<double> ASTERIX::decompressArrayDouble(char* compressedData, size_t compressedSize, size_t arraySize,double tol) {
   // Allocate memory for decompressed data
   std::vector<double> decompressedArray(arraySize);

   zfp_stream* zfp = zfp_stream_open(NULL);
   zfp_stream_set_accuracy(zfp, tol);
   bitstream* stream_decompress = stream_open(compressedData, compressedSize);
   zfp_stream_set_bit_stream(zfp, stream_decompress);
   zfp_stream_rewind(zfp);

   zfp_field* field_decompress = zfp_field_1d(decompressedArray.data(), zfp_type_double, decompressedArray.size());
   size_t retval = zfp_decompress(zfp, field_decompress);
   (void)retval;
   zfp_field_free(field_decompress);
   zfp_stream_close(zfp);
   stream_close(stream_decompress);
   return decompressedArray;
}
#endif //ASTERIX_ZFP

