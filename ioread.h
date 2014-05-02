#ifndef IOREAD_H
#define IOREAD_H
#include "mpi.h"
#include <dccrg.hpp>
#include <string>

#include "spatial_cell.hpp"
#include "datareduction/datareducer.h"


/*!

\brief Read in state from a vlsv file in order to restart simulations
\param mpiGrid Vlasiator's grid
\param name Name of the restart file e.g. "restart.00052.vlsv"
*/
bool readGrid(dccrg::Dccrg<spatial_cell::SpatialCell>& mpiGrid,
              const std::string& name);


/*!
 * \brief Check in local directory for external commands passed to the simulation. Only executed by MASTER_RANK
 * \return Returns an int for MPI_Bcast purposes. If the return is 1, a replay is performed.
 */
int checkExternalCommands();


#endif
