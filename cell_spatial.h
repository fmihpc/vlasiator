#ifndef CELL_SPATIAL
#define CELL_SPATIAL

#include <cstring>
#include <utility>
#include <vector>
#ifndef PARGRID
   #include <boost/serialization/array.hpp>
   #include <boost/serialization/split_member.hpp>
   #include <boost/mpi/datatype.hpp>
#else
   #include <mpi.h>
#endif
#include "definitions.h"
#include "common.h"
#include "parameters.h"
#include "grid.h"

extern Grid grid;

namespace Cell {
   enum CellType {INNER,BOUNDARY,GHOST,UNINITIALIZED};
   enum Array {Blocks,BlockParams,NbrsVel,CellParams,Fx,Fy,Fz,D1x,D1y,D1z};
   enum Dir {CpuToDev,DevToCpu};
}

struct SpatialCell {
   #ifndef PARGRID
   friend class boost::serialization::access;
   BOOST_SERIALIZATION_SPLIT_MEMBER()
   template<typename Archive> void save(Archive& ar,cuint& version) const;
   template<typename Archive> void load(Archive& ar,cuint& version);
   static size_t size();
   static uint base_address_identifier;
   void* at(void);
   #else
   void* baseAddress;
   void allocate();
   void* getBaseAddress(cuint& identifier);
   static void getMPIdatatype(cuint& identifier,MPI_Datatype& dataType);
   #endif
   uint cpuIndex;         /**< An index to Grid which is used to calculate the data array pointers (cpu_avgs etc.).*/
   uint N_blocks;         /**< Number of velocity blocks in this cell.*/

   // Pointers to arrays containing spatial cell parameters in CPU memory
   Real* cpu_cellParams;  /**< Pointer to physical cell parameters in cpu memory.*/
   uint* cpu_nbrsSpa;     /**< Pointer to spatial  neighbout list in CPU memory.*/
   
   // Pointers to arrays containing velocity block parameters in CPU memory
   uint* cpu_nbrsVel;     /**< Pointer to velocity neighbour list in CPU memory.*/
   Real* cpu_avgs;        /**< Pointer to velocity block array in CPU memory.*/
   Real* cpu_blockParams; /**< Pointer to velocity block parameter array in CPU memory.*/
   Real* cpu_fx;          /**< Pointer to x-flux array in CPU memory.*/
   Real* cpu_fy;          /**< Pointer to y-flux array in CPU memory.*/
   Real* cpu_fz;          /**< Pointer to z-flux array in CPU memory.*/
   Real* cpu_d1x;         /**< Pointer to array in CPU memory that contains 1st derivatives to x-direction.*/
   Real* cpu_d2x;         /**< Pointer to array in CPU memory that contains 2nd derivatives to x-direction.*/
   Real* cpu_d1y;         /**< Pointer to array in CPU memory that contains 1st derivatives to y-direction.*/
   Real* cpu_d2y;         /**< Pointer to array in CPU memory that contains 2nd derivatives to y-direction.*/
   Real* cpu_d1z;         /**< Pointer to array in CPU memory that contains 1st derivatives to z-direction.*/
   Real* cpu_d2z;         /**< Pointer to array in CPU memory that contains 2nd derivatives to z-direction.*/
   
   SpatialCell();
   SpatialCell(const SpatialCell& s);
   ~SpatialCell();
   SpatialCell& operator=(const SpatialCell& s);

   #ifndef NOCUDA
     bool devSync(const Cell::Array& array,const Cell::Dir direction);
   #endif

   bool allocateMemory(const uint& BLOCKS);
   bool freeMemory();
   bool clone(const SpatialCell& s);
   void getMemInfo();
   
 private:
   
};

#ifndef PARGRID
template<typename Archive> void SpatialCell::save(Archive& ar,cuint& version) const {
   typedef Parameters P;
   ar << BOOST_SERIALIZATION_NVP(N_blocks);
   if ((P::transmit & Transmit::CELL_PARAMS))  ar << boost::serialization::make_array(cpu_cellParams,SIZE_CELLPARAMS);
   if ((P::transmit & Transmit::BLOCK_PARAMS)) ar << boost::serialization::make_array(cpu_blockParams,N_blocks*SIZE_BLOCKPARAMS);
   if ((P::transmit & Transmit::AVGS))         ar << boost::serialization::make_array(cpu_avgs,N_blocks*SIZE_VELBLOCK);
   if ((P::transmit & Transmit::FLUXES)) {
      ar << boost::serialization::make_array(cpu_fx,N_blocks*SIZE_FLUXS);
      ar << boost::serialization::make_array(cpu_fy,N_blocks*SIZE_FLUXS);
      ar << boost::serialization::make_array(cpu_fz,N_blocks*SIZE_FLUXS);
   }
   if ((P::transmit & Transmit::DERIV1)) {
      ar << boost::serialization::make_array(cpu_d1x,N_blocks*SIZE_DERIV);
      ar << boost::serialization::make_array(cpu_d1y,N_blocks*SIZE_DERIV);
      ar << boost::serialization::make_array(cpu_d1z,N_blocks*SIZE_DERIV);
   }
   if ((P::transmit & Transmit::DERIV2)) {
      ar << boost::serialization::make_array(cpu_d2x,N_blocks*SIZE_DERIV);
      ar << boost::serialization::make_array(cpu_d2y,N_blocks*SIZE_DERIV);
      ar << boost::serialization::make_array(cpu_d2z,N_blocks*SIZE_DERIV);
   }
   if ((P::transmit & Transmit::NBRSVEL)) ar << boost::serialization::make_array(cpu_nbrsVel,N_blocks*SIZE_NBRS_VEL);
}

template<typename Archive> void SpatialCell::load(Archive& ar,cuint& version) {
   typedef Parameters P;
   uint size;
   ar >> BOOST_SERIALIZATION_NVP(size);
   if (size == 0) { // No incoming data, do nothing
      return;
   }
   
   // Allocate memory & increase reference count:
   N_blocks = size;
   if (cpuIndex == std::numeric_limits<uint>::max()) {
      cpuIndex = grid.getFreeMemory(N_blocks);
      if (cpuIndex == std::numeric_limits<uint>::max()) {
	 std::cerr << "SpatialCell load: failed to allocate memory, aborting" << std::endl; exit(1);
      }
   }
   // Copy incoming data:
   if (P::transmit & Transmit::CELL_PARAMS)    ar >> boost::serialization::make_array(cpu_cellParams,SIZE_CELLPARAMS);
   if ((P::transmit & Transmit::BLOCK_PARAMS)) { // Copy vel. block parameters
      cpu_blockParams = grid.getBlockParams() + cpuIndex*SIZE_BLOCKPARAMS;
      ar >> boost::serialization::make_array(cpu_blockParams,N_blocks*SIZE_BLOCKPARAMS);
   } 
   if ((P::transmit & Transmit::AVGS)) { // Copy vol. averages
      cpu_avgs = grid.getAvgs() + cpuIndex*SIZE_VELBLOCK;
      ar >> boost::serialization::make_array(cpu_avgs,N_blocks*SIZE_VELBLOCK);
   }
   if ((P::transmit & Transmit::FLUXES)) { // Copy flux data
      cpu_fx = grid.getFx() + cpuIndex*SIZE_FLUXS;
      cpu_fy = grid.getFy() + cpuIndex*SIZE_FLUXS;
      cpu_fz = grid.getFz() + cpuIndex*SIZE_FLUXS;
      ar >> boost::serialization::make_array(cpu_fx,N_blocks*SIZE_FLUXS);
      ar >> boost::serialization::make_array(cpu_fy,N_blocks*SIZE_FLUXS);
      ar >> boost::serialization::make_array(cpu_fz,N_blocks*SIZE_FLUXS);
   }
   if ((P::transmit & Transmit::DERIV1)) { // Copy 1st derivativesx
      cpu_d1x = grid.getD1x() + cpuIndex*SIZE_DERIV;
      cpu_d1y = grid.getD1y() + cpuIndex*SIZE_DERIV;
      cpu_d1z = grid.getD1z() + cpuIndex*SIZE_DERIV;
      ar >> boost::serialization::make_array(cpu_d1x,N_blocks*SIZE_DERIV);
      ar >> boost::serialization::make_array(cpu_d1y,N_blocks*SIZE_DERIV);
      ar >> boost::serialization::make_array(cpu_d1z,N_blocks*SIZE_DERIV);
   }
   if ((P::transmit & Transmit::DERIV2)) { // Copy 2nd derivatives
      cpu_d2x = grid.getD2x() + cpuIndex*SIZE_DERIV;
      cpu_d2y = grid.getD2y() + cpuIndex*SIZE_DERIV;
      cpu_d2z = grid.getD2z() + cpuIndex*SIZE_DERIV;
      ar >> boost::serialization::make_array(cpu_d2x,N_blocks*SIZE_DERIV);
      ar >> boost::serialization::make_array(cpu_d2y,N_blocks*SIZE_DERIV);
      ar >> boost::serialization::make_array(cpu_d2z,N_blocks*SIZE_DERIV);
   }
   if ((P::transmit & Transmit::NBRSVEL)) { // Vel. block neighbour lists
      cpu_nbrsVel = grid.getNbrsVel() + cpuIndex*SIZE_NBRS_VEL;
      ar >> boost::serialization::make_array(cpu_nbrsVel,N_blocks*SIZE_NBRS_VEL);
   }
}
#endif // #ifndef PARGRID

#endif

