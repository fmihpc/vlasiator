#ifndef IOWRITE_H
#define IOWRITE_H
#include "mpi.h"
#include <dccrg.hpp>
#include <dccrg_cartesian_geometry.hpp>
#include <string>
#include <vector>
#include <vlsv_writer.h>

#include "spatial_cell.hpp"
#include "datareduction/datareducer.h"

/*!

\brief Write out system into a vlsv file

\param mpiGrid     The DCCRG grid with spatial cells
\param dataReducer Contains datareductionoperators that are used to compute data that is added into file
\param index       Index to call the correct member of the various parameter vectors
\param writeGhosts Write ghost zones
*/
bool writeGrid(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
               DataReducer& dataReducer,
               const uint& index,
               const bool writeGhosts = true
);

/*!

\brief Write out a restart of the simulation into a vlsv file. All block data in remote cells will be reset.

\param mpiGrid   The DCCRG grid with spatial cells
\param dataReducer Contains datareductionoperators that are used to compute data that is added into file
\param name       File name prefix, file will be called "name.index.vlsv"
\param fileIndex  File index, file will be called "name.index.vlsv"
*/
bool writeRestart(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                  DataReducer& dataReducer,
                  const std::string& name,
                  const uint& fileIndex,
                  const int& stripe);

/*!

\brief Write out simulation diagnostics into diagnostic.txt

@param mpiGrid   The DCCRG grid with spatial cells
@param dataReducer Contains datareductionoperators that are used to compute diagnostic data
*/
bool writeDiagnostic(const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,DataReducer& dataReducer);

bool writeVelocitySpace(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                        vlsv::Writer& vlsvWriter,int index,const std::vector<uint64_t>& cells);

bool writeVelocityDistributionData(vlsv::Writer& vlsvWriter,dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                   const std::vector<uint64_t>& cells,MPI_Comm comm);

#endif
