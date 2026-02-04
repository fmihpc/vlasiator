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

// These tools  are fwd declared here and implemented at the end of the file for
// better clarity. They are not for external usage and as such they do not go
// into the header file

#include "../definitions.h"
#include "../object_wrapper.h"
#include <dccrg.hpp>
#include "../logger.h"
#include "../mpiconversion.h"
#include "../spatial_cells/spatial_cell_wrapper.hpp"
#include "../spatial_cells/velocity_block_container.h"
#include "stdlib.h"
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <limits>
#include <span>
#include <stdexcept>
#include <unordered_map>
#include <vector>
#include <cstring>
#ifdef ASTERIX_OCTREE
#include "toctree_compressor.h"
#endif
#ifdef ASTERIX_MLP
#include "lib/usr/local/include/tinyAI/genericTsPool.h"
#endif
#ifdef ASTERIX_ZFP
#include <zfp.h>
#include "zfp/array1.hpp"
#endif

#ifdef ASTERIX_MLP
#ifdef __cplusplus
extern "C" {
#endif
void decompress_phasespace6D_f64(GENERIC_TS_POOL::MemPool* p, std::size_t fin, std::size_t fout, double* vcoords_ptr,
                                 double* vspace_ptr, std::size_t size, std::size_t fourier_order,
                                 size_t* hidden_layers_ptr, size_t n_hidden_layers, double* weights_ptr,
                                 std::size_t weight_size, bool use_input_weights);

void decompress_phasespace6D_f32(GENERIC_TS_POOL::MemPool* p, std::size_t fin, std::size_t fout, float* vcoords_ptr,
                                 float* vspace_ptr, std::size_t size, std::size_t fourier_order,
                                 size_t* hidden_layers_ptr, size_t n_hidden_layers, float* weights_ptr,
                                 std::size_t weight_size, bool use_input_weights);
#ifdef __cplusplus
}
#endif
#endif //ASTERIX_MLP

#define MLP_KEY 42

namespace ASTERIX {
struct VCoords {
   Real vx, vy, vz;
   VCoords operator+(const VCoords& other) { return {vx + other.vx, vy + other.vy, vz + other.vz}; }
   VCoords operator-(const VCoords& other) { return {vx - other.vx, vy - other.vy, vz - other.vz}; }
};

struct OrderedVDF {

   struct OrderedVDFHeader {
      std::size_t ignore_bytes;
      std::size_t octree_bytes;
   };

   std::vector<vmesh::GlobalID> blocks_to_ignore;
   std::size_t sparse_vdf_bytes = {0};
   std::vector<Realf> vdf_vals;
   std::array<Real, 6> v_limits;     // vx_min,vy_min,vz_min,vx_max,vy_max,vz_max
   std::array<std::size_t, 3> shape; // x,y,z
   std::size_t index(std::size_t i, std::size_t j, std::size_t k) const noexcept {
      return i * (shape[1] * shape[2]) + j * shape[2] + k;
   }

   Realf& at(std::size_t i, std::size_t j, std::size_t k) noexcept { return vdf_vals.at(index(i, j, k)); }

   const Realf& at(std::size_t i, std::size_t j, std::size_t k) const noexcept { return vdf_vals.at(index(i, j, k)); }

   bool save_to_file(const char* filename) const noexcept {
      std::ofstream file(filename, std::ios::out | std::ios::binary);
      if (!file) {
         std::cerr << "Could not open file for writting! Exiting!" << std::endl;
         return false;
      }
      file.write((char*)shape.data(), 3 * sizeof(size_t));
      if (!file) {
         std::cerr << "Error writing shape data to file!" << std::endl;
         return false;
      }

      file.write((char*)vdf_vals.data(), vdf_vals.size() * sizeof(Realf));
      if (!file) {
         std::cerr << "Error writing vdf_vals data to file!" << std::endl;
         return false;
      }
      return true;
   }
};

struct UnorderedVDF {
   std::vector<Realf> vdf_vals;
   std::vector<std::array<Real, 3>> vdf_coords;
   std::array<Real, 6> v_limits{std::numeric_limits<Real>::max(),    std::numeric_limits<Real>::max(),
                                std::numeric_limits<Real>::max(),    std::numeric_limits<Real>::lowest(),
                                std::numeric_limits<Real>::lowest(), std::numeric_limits<Real>::lowest()};

   bool save_to_file(const char* filename) const noexcept {
      std::ofstream file(filename, std::ios::out | std::ios::binary);
      if (!file) {
         std::cerr << "Could not open file for writting! Exiting!" << std::endl;
         return false;
      }

      file.write((char*)vdf_vals.size(), sizeof(size_t));
      if (!file) {
         std::cerr << "Error writing size data to file!" << std::endl;
         return false;
      }

      file.write((char*)v_limits.data(), 6 * sizeof(Real));
      if (!file) {
         std::cerr << "Error writing size data to file!" << std::endl;
         return false;
      }

      file.write((char*)vdf_coords.data(), vdf_coords.size() * 3 * sizeof(Real));
      if (!file) {
         std::cerr << "Error writing vdf_coords data to file!" << std::endl;
         return false;
      }

      file.write((char*)vdf_vals.data(), vdf_vals.size() * sizeof(Realf));
      if (!file) {
         std::cerr << "Error writing vdf_vals data to file!" << std::endl;
         return false;
      }
      return true;
   }
};

template <typename T> class PhaseSpaceUnion {
public:
   //---- No need for these------
   PhaseSpaceUnion(const PhaseSpaceUnion& other) = delete;
   PhaseSpaceUnion(PhaseSpaceUnion&& other) = delete;
   PhaseSpaceUnion& operator=(const PhaseSpaceUnion& other) = delete;
   PhaseSpaceUnion& operator=(PhaseSpaceUnion&& other) = delete;
   //----------------------------
   explicit PhaseSpaceUnion(const unsigned char* buffer) { deserialize_from(buffer); }

   explicit PhaseSpaceUnion(const std::span<const CellID> cids, uint popID,
                   const dccrg::Dccrg<spatial_cell::SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid, bool center_vdfs)
       : _center_vdfs(center_vdfs) {

      // Let's find out which of these cellids has the largest VDF
      std::size_t max_cid_block_size = 0;
      std::size_t bytes_of_all_local_vdfs = 0;
      for (const auto& cid : cids) {
         spatial_cell::SpatialCell* sc = mpiGrid[cid];
         const auto blockContainer = sc->get_velocity_blocks(popID);
         const size_t total_size = blockContainer->size();
         max_cid_block_size = std::max(total_size, max_cid_block_size);
         bytes_of_all_local_vdfs += total_size * WID3 * sizeof(Realf);
      }
      _effective_vdf_size = bytes_of_all_local_vdfs;

      std::vector<std::vector<T>> vspaces(cids.size());
      std::vector<double> f_sums(cids.size(), 0);
      const Real sparse = static_cast<double>(getObjectWrapper().particleSpecies[popID].sparseMinValue);
      for (std::size_t cc = 0; cc < cids.size(); ++cc) {
         const auto& cid = cids[cc];
         _cids.push_back(cid);
         spatial_cell::SpatialCell* sc = mpiGrid[cid];
         auto blockContainer = sc->get_velocity_blocks(popID);
         const std::array<T, 3> bulkv{static_cast<T>(sc->get_population(popID).V[0]),
                                      static_cast<T>(sc->get_population(popID).V[1]),
                                      static_cast<T>(sc->get_population(popID).V[2])};

         _vbulks.push_back(bulkv);
         const size_t total_blocks = blockContainer->size();
         Realf* data = blockContainer->getData();
         const Real* blockParams = sc->get_block_parameters(popID);
         for (std::size_t n = 0; n < total_blocks; ++n) {
            const auto bp = blockParams + n * BlockParams::N_VELOCITY_BLOCK_PARAMS;
            const vmesh::GlobalID gid = sc->get_velocity_block_global_id(n, popID);
            const Realf* vdf_data = &data[n * WID3];

            auto [it, block_inserted] = _map.try_emplace(gid, _vcoords.size());
            std::size_t cnt = 0;
            for (uint k = 0; k < WID; ++k) {
               for (uint j = 0; j < WID; ++j) {
                  for (uint i = 0; i < WID; ++i) {

                     std::array<T, 3> coords = {
                         static_cast<T>(bp[BlockParams::VXCRD] + (i + 0.5) * bp[BlockParams::DVX]),
                         static_cast<T>(bp[BlockParams::VYCRD] + (j + 0.5) * bp[BlockParams::DVY]),
                         static_cast<T>(bp[BlockParams::VZCRD] + (k + 0.5) * bp[BlockParams::DVZ])};

                     if (_center_vdfs) {
                        coords[0] = coords[0] - bulkv[0];
                        coords[1] = coords[1] - bulkv[1];
                        coords[2] = coords[2] - bulkv[2];
                     }

                     _v_limits[0] = std::min(_v_limits[0], static_cast<T>(coords[0]));
                     _v_limits[1] = std::min(_v_limits[1], static_cast<T>(coords[1]));
                     _v_limits[2] = std::min(_v_limits[2], static_cast<T>(coords[2]));
                     _v_limits[3] = std::max(_v_limits[3], static_cast<T>(coords[0]));
                     _v_limits[4] = std::max(_v_limits[4], static_cast<T>(coords[1]));
                     _v_limits[5] = std::max(_v_limits[5], static_cast<T>(coords[2]));
                     const double vdf_val = static_cast<double>(vdf_data[cellIndex(i, j, k)]);
                     if (block_inserted) { // which means the block was not there before
                        _vcoords.push_back({coords[0], coords[1], coords[2]});
                        for (std::size_t x = 0; x < cids.size(); ++x) {
                           vspaces[x].push_back((x == cc) ? vdf_val : T(0));
                        }
                     } else { // So the block was there
                        vspaces[cc].at(it->second + cnt) = vdf_val;
                     }
                     f_sums.at(cc) += static_cast<double>(vdf_val);
                     cnt++;
                  }
               }
            }
         }
      }
      _nrows = vspaces.front().size();
      _ncols = cids.size();
      // This will be used further down for indexing into the vspace_union
      auto index_2d = [this](std::size_t row, std::size_t col) -> std::size_t { return row * _ncols + col; };

      // Resize to fit the union of vspace coords and vspace density
      _vspace = std::move(std::vector<T>(_nrows * _ncols, T(0)));
      for (std::size_t i = 0; i < _nrows; ++i) {
         for (std::size_t j = 0; j < _ncols; ++j) {
            _vspace.at(index_2d(i, j)) = vspaces[j][i];
         }
      }
      // Scale now
      scale(sparse);
      // Store norms
      _norms = std::move(std::vector<Norms>(_ncols, Norms{}));
   }

   void scale(T sparse) noexcept {
      std::for_each(_vspace.begin(), _vspace.end(),
                    [sparse](T& value) { value = std::log10(std::max(value, sparse)) - std::log10(sparse); });
   }

   //Here we min max normalize both VDF(per VDF minmax) and VCOORDS
   //TODO remove uneeded stuff from NORMS
   void normalize() noexcept {
      // Vcoords
      std::ranges::for_each(_vcoords, [this](std::array<T, 3>& x) {
         x[0] = 2.0 * ((x[0] - _v_limits[0]) / (_v_limits[3] - _v_limits[0])) - 1.0;
         x[1] = 2.0 * ((x[1] - _v_limits[1]) / (_v_limits[4] - _v_limits[1])) - 1.0;
         x[2] = 2.0 * ((x[2] - _v_limits[2]) / (_v_limits[5] - _v_limits[2])) - 1.0;
      });
      // Per VDF min max scaling
      const std::size_t nVDFS = _ncols;
      for (std::size_t v = 0; v < nVDFS; ++v) {
         T min_val = std::numeric_limits<T>::max();
         T max_val = std::numeric_limits<T>::lowest();
         for (std::size_t i = 0; i < _nrows; ++i) {
            min_val = std::min(min_val, _vspace[index_2d(i, v)]);
            max_val = std::max(max_val, _vspace[index_2d(i, v)]);
         }
         for (std::size_t i = 0; i < _nrows; ++i) {
            _vspace[index_2d(i, v)] /= max_val;
         }
         _norms[v] = Norms{.min = min_val, .max = max_val};
      }
   }

   //We do the reverse of normalize()
   void unormalize_and_unscale(T sparse) noexcept {
      std::ranges::for_each(_vcoords, [this](std::array<T, 3>& x) {
         x[0] = ((x[0] + 1.0) / 2.0) * (_v_limits[3] - _v_limits[0]) + _v_limits[0];
         x[1] = ((x[1] + 1.0) / 2.0) * (_v_limits[4] - _v_limits[1]) + _v_limits[1];
         x[2] = ((x[2] + 1.0) / 2.0) * (_v_limits[5] - _v_limits[2]) + _v_limits[2];
      });

      const std::size_t nVDFS = _ncols;
      for (std::size_t v = 0; v < nVDFS; ++v) {
         const T max_val = _norms[v].max;
         const T min_val = _norms[v].min;
         for (std::size_t i = 0; i < _nrows; ++i) {
            _vspace[index_2d(i, v)] = sparse*std::pow(10.0, _vspace[index_2d(i, v)]*max_val );
         }
      }
   }

   constexpr std::size_t index_2d(std::size_t row, std::size_t col) const noexcept { return row * _ncols + col; };

   void sparsify(T sparse) noexcept {
      std::for_each(_vspace.begin(), _vspace.end(), [sparse](T& x) {
         if (x - sparse < 0.0) {
            x = 0.0;
         }
      });
   }

   std::size_t total_serialized_size_bytes() const {
      return sizeof(Header) + _cids.size() * sizeof(CellID) + _norms.size() * sizeof(Norms) +
             _vbulks.size() * sizeof(std::array<T, 3>) + _vcoords.size() * sizeof(std::array<T, 3>) + 6 * sizeof(T) +
             _n_weights * sizeof(T) + _map.size() * sizeof(std::pair<vmesh::LocalID, std::size_t>);
      ;
   }

   void serialize_into(unsigned char* buffer) const {
      Header header;
      header.key = MLP_KEY;
      header.total_size = total_serialized_size_bytes();
      header.rows = _nrows;
      header.cols = _ncols;
      header.n_weights = _n_weights;
      header.type_size = sizeof(T);
      std::size_t write_index = 0;

      std::memcpy(&buffer[write_index], &header, sizeof(Header));
      write_index += sizeof(Header);

      std::memcpy(&buffer[write_index], &_cids[0], _cids.size() * sizeof(CellID));
      write_index += _cids.size() * sizeof(CellID);

      std::memcpy(&buffer[write_index], &_norms[0], _norms.size() * sizeof(Norms));
      write_index += _norms.size() * sizeof(Norms);

      std::memcpy(&buffer[write_index], &_vbulks[0], _vbulks.size() * sizeof(std::array<T, 3>));
      write_index += _vbulks.size() * sizeof(std::array<T, 3>);

      std::memcpy(&buffer[write_index], &_v_limits[0], 6 * sizeof(T));
      write_index += 6 * sizeof(T);

      std::memcpy(&buffer[write_index], &_vcoords[0], _vcoords.size() * sizeof(std::array<T, 3>));
      write_index += _vcoords.size() * sizeof(std::array<T, 3>);

      std::memcpy(&buffer[write_index], &_network_weights[0], _n_weights * sizeof(T));
      write_index += _n_weights * sizeof(T);

      for (const auto& kval : _map) {
         std::memcpy(&buffer[write_index], &kval, sizeof(std::pair<vmesh::LocalID, std::size_t>));
         write_index += sizeof(std::pair<vmesh::LocalID, std::size_t>);
      }
      assert(header.total_size == write_index);
      if (!(header.total_size == write_index)) {
         throw std::runtime_error("Failed to fully write state");
      }
   }

   void deserialize_from(const unsigned char* buffer) {
      const Header* const header = reinterpret_cast<const Header*>(&buffer[0]);
      assert(header->key == MLP_KEY && "Blame Kostis Papadakis for this!");
      if (!(header->key == MLP_KEY)){
         throw std::runtime_error("Wrong MLP Header KEY");
      }

      // Inflate vspave union
      _vspace.resize(header->cols * header->rows);
      _ncols = header->cols;
      _nrows = header->rows;

      // Recover cids in this union;
      std::size_t read_index = sizeof(Header);
      std::size_t cids_size = header->cols;
      _cids.resize(cids_size);

      std::memcpy(_cids.data(), &buffer[read_index], cids_size * sizeof(CellID));
      read_index += cids_size * sizeof(CellID);

      std::size_t norms_size = cids_size;
      _norms.resize(cids_size);
      std::memcpy(_norms.data(), &buffer[read_index], norms_size * sizeof(Norms));
      read_index += norms_size * sizeof(Norms);

      std::size_t vbulk_size = cids_size;
      _vbulks.resize(vbulk_size);
      std::memcpy(_vbulks.data(), &buffer[read_index], vbulk_size * sizeof(std::array<T, 3>));
      read_index += vbulk_size * sizeof(std::array<T, 3>);

      std::memcpy(&_v_limits[0], &buffer[read_index], 6 * sizeof(T));
      read_index += 6 * sizeof(T);

      std::size_t vcoords_size = header->rows;
      _vcoords.resize(vcoords_size);
      std::memcpy(_vcoords.data(), &buffer[read_index], vcoords_size * sizeof(std::array<T, 3>));
      read_index += vcoords_size * sizeof(std::array<T, 3>);

      if (_network_weights != nullptr) {
         free(_network_weights);
      }

      _network_weights = (T*)malloc(header->n_weights * sizeof(T));
      _n_weights = header->n_weights;
      std::memcpy(_network_weights, &buffer[read_index], header->n_weights * sizeof(T));
      read_index += _n_weights * sizeof(T);

      while (read_index < header->total_size) {
         const std::pair<vmesh::LocalID, std::size_t>* kval =
             reinterpret_cast<const std::pair<vmesh::LocalID, std::size_t>*>(&buffer[read_index]);
         _map[kval->first] = kval->second;
         read_index += sizeof(std::pair<vmesh::LocalID, std::size_t>);
      }
      assert(read_index == header->total_size && "Size mismatch while reading in serialized VDF Union!");
      if (!(read_index == header->total_size)){
         throw std::runtime_error("Failed to fully read state");
      }
   }

   struct Norms {
      double min = std::numeric_limits<double>::max();
      double max = std::numeric_limits<double>::min();
   };

   struct Header {
      std::size_t key;
      std::size_t total_size;
      std::size_t rows;
      std::size_t cols;
      std::size_t n_weights;
      std::size_t type_size;
   };

   std::size_t _nrows = {0};
   std::size_t _ncols = {0};
   bool _center_vdfs;

   std::vector<Norms> _norms;
   std::vector<CellID> _cids;
   std::vector<std::array<T, 3>> _vcoords;
   std::vector<std::array<T, 3>> _vbulks;
   std::vector<T> _vspace;
   std::unordered_map<vmesh::GlobalID, std::size_t> _map;
   T* _network_weights = nullptr;
   std::size_t _effective_vdf_size = {0};
   std::size_t _n_weights = {0};
   std::array<T, 6> _v_limits{std::numeric_limits<T>::max(),    std::numeric_limits<T>::max(),
                              std::numeric_limits<T>::max(),    std::numeric_limits<T>::lowest(),
                              std::numeric_limits<T>::lowest(), std::numeric_limits<T>::lowest()};
};

auto extract_pop_vdf_from_spatial_cell(spatial_cell::SpatialCell* sc, uint popID) -> UnorderedVDF;

auto extract_pop_vdf_from_spatial_cell_ordered_min_bbox_zoomed(spatial_cell::SpatialCell* sc, uint popID, int zoom) -> OrderedVDF;

constexpr auto isPow2(std::unsigned_integral auto val) -> bool { return (val & (val - 1)) == 0; };

auto overwrite_pop_spatial_cell_vdf(spatial_cell::SpatialCell* sc, uint popID, const std::vector<Realf>& new_vspace) -> void;

auto overwrite_pop_spatial_cell_vdf(spatial_cell::SpatialCell* sc, uint popID, const OrderedVDF& vdf) -> void;

auto overwrite_cellids_vdfs(const std::span<const CellID> cids, uint popID,
                            dccrg::Dccrg<spatial_cell::SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid,
                            const std::vector<std::array<Real, 3>>& vcoords, const std::vector<Realf>& vspace_union,
                            const std::unordered_map<vmesh::LocalID, std::size_t>& map_exists_id) -> void;

auto dump_vdf_to_binary_file(const char* filename, CellID cid) -> void;

auto dump_vdf_to_binary_file(const char* filename, CellID cid,
                             dccrg::Dccrg<spatial_cell::SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid) -> void;

// https://en.wikipedia.org/wiki/Entropy_(information_theory)
template <typename T>
requires(std::is_same_v<T, float> || std::is_same_v<T, double>) auto shannon_entropy(const std::vector<T>& data) -> T {
   const std::size_t sz = data.size();
   if (sz == 0) {
      return 0.0;
   }

   using key_t = std::conditional_t<std::is_same_v<T, float>, uint32_t, uint64_t>;
   std::unordered_map<key_t, int> frequency;
   for (std::size_t i = 0; i < sz; ++i) {
      frequency[*(reinterpret_cast<const key_t*>(&data[i]))]++;
   }
   T entropy = 0.0;
   for (const auto& [byte, count] : frequency) {
      T pk = static_cast<T>(count) / sz;
      entropy -= pk * std::log2(pk);
   }
   return entropy;
}

template <typename T>
requires(std::is_same_v<T, float> ||
         std::is_same_v<T, double>) auto theoritical_lossless_compression_ratio(const std::vector<T>& data,
                                                                                std::size_t bits) -> T {
   T entorpy = shannon_entropy(data);
   return static_cast<T>(bits) / entorpy;
}

template <typename NetworkType>
requires(std::is_same_v<NetworkType, float> || std::is_same_v<NetworkType, double>) auto calculate_total_size_bytes(
    const std::vector<std::size_t>& architecture, std::size_t fourier_order, std::size_t output_dim) -> std::size_t {
   if (architecture.empty()) {
      throw std::runtime_error("Architecture cannot be empty.");
   }
   std::size_t input_dim = 2 * fourier_order;
   std::size_t total_size = 0;
   total_size += input_dim * architecture[0];
   total_size += architecture[0];

   for (std::size_t i = 1; i < architecture.size(); ++i) {
      total_size += architecture[i - 1] * architecture[i];
      total_size += architecture[i];
   }

   total_size += architecture.back() * output_dim;
   total_size += output_dim;

   return total_size * sizeof(NetworkType);
}

template <typename NetworkType>
requires(std::is_same_v<NetworkType, float> || std::is_same_v<NetworkType, double>) auto calculate_hidden_neurons(
    std::size_t N_input, std::size_t N_output, std::size_t num_hidden_layers, std::size_t target_size)
    -> std::vector<std::size_t> {
   std::vector<std::size_t> neurons(num_hidden_layers + 2); // 2 input and output
   neurons[0] = N_input;
   neurons[num_hidden_layers + 1] = N_output;

   // We guess this heyuristically
   std::size_t initial_hidden_size = 1;
   for (std::size_t i = 1; i <= num_hidden_layers; ++i) {
      neurons[i] = initial_hidden_size;
   }
   std::size_t current_size = calculate_total_size_bytes<NetworkType>(neurons);

   while (current_size < target_size) {
      for (std::size_t i = 1; i <= num_hidden_layers; ++i) {
         neurons[i]++;
      }
      current_size = calculate_total_size_bytes<NetworkType>(neurons);
   }

   while (current_size > target_size) {
      for (std::size_t i = 1; i <= num_hidden_layers; ++i) {
         if (neurons[i] > 1) {
            neurons[i]--;
         }
      }
      current_size = calculate_total_size_bytes<NetworkType>(neurons);
   }
   return neurons;
}
Real get_Non_MaxWellianity(const spatial_cell::SpatialCell* cell, uint popID);

#ifdef ASTERIX_MLP
template <typename T> void decompressPhaseSpace(PhaseSpaceUnion<T>& rv) {
   // Memory allocation
   GENERIC_TS_POOL::MemPool p{};
   if constexpr (sizeof(T) == sizeof(float)) {
      decompress_phasespace6D_f32(&p, 3, rv._ncols, &rv._vcoords[0][0], rv._vspace.data(), rv._vcoords.size(),
                                  P::mlp_fourier_order, P::mlp_arch.data(), P::mlp_arch.size(), rv._network_weights,
                                  rv._n_weights * sizeof(float), true);
   } else {
      decompress_phasespace6D_f64(&p, 3, rv._ncols, &rv._vcoords[0][0], rv._vspace.data(), rv._vcoords.size(),
                                  P::mlp_fourier_order, P::mlp_arch.data(), P::mlp_arch.size(), rv._network_weights,
                                  rv._n_weights * sizeof(double), true);
   }
}
#endif //ASTERIX_MLP

template <typename T>
void overwrite_cellids_vdf_single_cell(const std::span<const CellID> cids, uint popID, spatial_cell::SpatialCell* sc, size_t cc,
                                       const std::vector<std::array<T, 3>>& vcoords, const std::vector<T>& vspace_union,
                                       const std::unordered_map<vmesh::LocalID, std::size_t>& map_exists_id) {
   const std::size_t nrows = vcoords.size();
   const std::size_t ncols = cids.size();
   // This will be used further down for indexing into the vspace_union
   auto index_2d = [nrows, ncols](std::size_t row, std::size_t col) -> std::size_t { return row * ncols + col; };

   const auto& cid = cids[cc];
   auto blockContainer = sc->get_velocity_blocks(popID);
   const size_t total_blocks = blockContainer->size();
   Realf* data = blockContainer->getData();
   const Real* blockParams = sc->get_block_parameters(popID);
   for (std::size_t n = 0; n < total_blocks; ++n) {
      const auto bp = blockParams + n * BlockParams::N_VELOCITY_BLOCK_PARAMS;
      const vmesh::GlobalID gid = sc->get_velocity_block_global_id(n, popID);
      const auto it = map_exists_id.find(gid);
      const bool exists = it != map_exists_id.end();
      if (!exists) {
         continue;
      }
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
   return;
}

} // namespace ASTERIX
