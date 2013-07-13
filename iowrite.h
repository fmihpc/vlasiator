#ifndef IOWRITE_H
#define IOWRITE_H
#include "mpi.h"
#include <dccrg.hpp>
#include <string>

#include "spatial_cell.hpp"
#include "datareduction/datareducer.h"


/*!

\brief Write out system into a vlsv file

\param mpiGrid     The DCCRG grid with spatial cells
\param dataReducer Contains datareductionoperators that are used to compute data that is added into file
\param index       Index to call the correct member of the various parameter vectors
\param newLib      Use the updated version of VLSV library
\param writeGhosts Write ghost zones
*/


bool writeGrid(
   dccrg::Dccrg<spatial_cell::SpatialCell>& mpiGrid,
   DataReducer& dataReducer,
   const uint& index,
   const bool writeSmaller = false,
   const bool writeGhosts = true
);


/*!

\brief Write out a restart of the simulation into a vlsv file

\param mpiGrid   The DCCRG grid with spatial cells
\param dataReducer Contains datareductionoperators that are used to compute data that is added into file
\param name       File name prefix, file will be called "name.index.vlsv"
\param index      File index, file will be called "name.index.vlsv"
*/


bool writeRestart(dccrg::Dccrg<spatial_cell::SpatialCell>& mpiGrid,
               DataReducer& dataReducer,
               const std::string& name,
               const uint& index,
               const int& stripe);





/*!

\brief Write out simulation diagnostics into diagnostic.txt

\param mpiGrid   The DCCRG grid with spatial cells
\param dataReducer Contains datareductionoperators that are used to compute diagnostic data
*/

bool writeDiagnostic(const dccrg::Dccrg<SpatialCell>& mpiGrid, DataReducer& dataReducer);

#endif
