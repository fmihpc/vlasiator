#ifndef IOREAD_H
#define IOREAD_H
#include "mpi.h"
#include <dccrg.hpp>
#include <string>

#include "spatial_cell.hpp"
#include "datareduction/datareducer.h"




/*!

\brief Read in state from a vlsv file in order to restart simulations
*/
bool readGrid(dccrg::Dccrg<spatial_cell::SpatialCell>& mpiGrid,
              const std::string& name, bool newLib = false);

#endif
