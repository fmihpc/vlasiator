#pragma once
#include <fsgrid.hpp>
#include <vector>
#include <array>

/*! Take input moments from DCCRG grid and put them into the Fieldsolver grid
 * \param mpiGrid The DCCRG grid carrying rho, rhoV and P
 * \param cells List of local cells
 * \param momentsGrid Fieldsolver grid for these quantities
 *
 * This function assumes that proper grid coupling has been set up.
 */
void feedMomentsIntoFsGrid(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                           const std::vector<CellID>& cells,
                           FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, 2>& momentsGrid);

void getVolumeFieldsFromFsGrid(FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, 2>& momentsGrid,
                           dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                           const std::vector<CellID>& cells);
