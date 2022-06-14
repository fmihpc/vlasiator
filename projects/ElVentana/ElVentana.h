/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
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

#ifndef ELVENTANA_H
#define ELVENTANA_H

#include "vlsv_reader_parallel.h"
#include "vlsv_reader.h"

#include "../../definitions.h"
#include "../projectTriAxisSearch.h"
#include "../../spatial_cell.hpp"

namespace projects {

   struct ElVentanaSpeciesParameters {
      Real rho;
      Real T;
      Real V0[3];
      Real ionosphereV0[3];
      Real ionosphereRho;
      uint nSpaceSamples;
      uint nVelocitySamples;
      Real Temperatureratio;
   };

   class ElVentana: public TriAxisSearch {
   friend class Magnetosphere;
    public:
      ElVentana();
      virtual ~ElVentana();
      
      virtual bool initialize(void);
      static void addParameters(void);
      virtual void getParameters(void);
      //virtual void setCellBackgroundField(spatial_cell::SpatialCell* cell) const;
      virtual void setProjectBField(
         FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2>& perBGrid,
         FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, 2>& BgBGrid,
         FsGrid< fsgrids::technical, 2>& technicalGrid
      );
      virtual Real calcPhaseSpaceDensity(
                                         creal& x, creal& y, creal& z,
                                         creal& dx, creal& dy, creal& dz,
                                         creal& vx, creal& vy, creal& vz,
                                         creal& dvx, creal& dvy, creal& dvz,
                                         const uint popID
                                        ) const;
      virtual void setupBeforeSetCell(const std::vector<CellID>& cells,
				      dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
				      bool& needCurl);

      virtual Real getCorrectNumberDensity(spatial_cell::SpatialCell* cell,const uint popID) const;
      virtual bool rescalesDensity(const uint popID) const;
      //void rescaleDensity(spatial_cell::SpatialCell* cell,const uint popID) const;

    protected:
      Real getDistribValue(
                           creal& x,creal& y, creal& z,
                           creal& vx, creal& vy, creal& vz,
                           creal& dvx, creal& dvy, creal& dvz,
                           const uint popID
                          ) const;
      virtual void calcCellParameters(spatial_cell::SpatialCell* cell,creal& t);
      virtual std::vector<std::array<Real, 3> > getV0(
                                                      creal x,
                                                      creal y,
                                                      creal z,
                                                      const uint popID
                                                     ) const;

      dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry> *newmpiGrid = NULL;
      Real constBgB[3];
      bool noDipoleInSW;
      Real WindowX[2];
      Real WindowY[2];
      Real WindowZ[2];
      std::string StartFile;
      std::vector<ElVentanaSpeciesParameters> speciesParams;
      vlsv::ParallelReader vlsvParaReader;
      vlsv::Reader vlsvSerialReader;
      uint64_t vecsizepressure;
      CellID findCellIDXYZ(creal x, creal y, creal z) const;

      Real* readVar(std::string varname, CellID fileOffset, uint64_t &vecsize);
      std::string pickVarName(const std::string &grid, const std::list<std::string> &varNames);
      bool readGridSize(std::array<double, 3> &fileMin, std::array<double, 3> &fileMax, std::array<uint64_t, 3> &fileCells, std::array<double, 3> &fileDx);
      template<unsigned long int N> bool readFsGridVariable(const std::string& variableName, FsGrid<std::array<Real, N>,2>& targetGrid);
      CellID getOldCellID(CellID newID, dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid, std::array<CellID, 3> fileCells, std::array<double, 3> &fileMin, std::array<double, 3> fileD);
      bool refineSpatialCells( dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid);

      //bool includeIonosphere;
      Real ionosphereRadius;
      uint ionosphereGeometry;
      Real center[3];
      Real dipoleScalingFactor;
      uint dipoleType;
      bool perBSet = false;
      bool totalBRead = false;

   }; // class ElVentana
} // namespace projects

#endif

