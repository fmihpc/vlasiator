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
#include <mpi.h>
#include <iostream>
#include <string>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
#include <string.h>
#include "common.h"
#include "particles/field.h"
#include "particles/readfields.h"

using namespace std;

static bool isInside(Field& B, double R, int x, int y, int z) {
   double R2 = 0;
   std::array<int, 3> xyz{x, y, z};
   for (int i = 0; i < 3; ++i) {
      R2 += pow(B.dimension[i]->min + xyz[i] * B.dx[i], 2);
   }
   return R2 < pow(R, 2);
}

// Calculate fluxfunction by integrating along -y/z boundary first,
// and then going along y/z-direction.
std::vector<double> computeFluxUp(Field& B, int outerBoundary, double innerBoundary) {
   // Create fluxfunction-field to be the same shape as B
   std::vector<double> flux(B.dimension[0]->cells * B.dimension[1]->cells * B.dimension[2]->cells);
   for (int i = 0; i < flux.size(); ++i) {
      flux[i] = NAN;
   }

   bool eqPlane = B.dimension[1]->cells > 1;
   int yCoord = eqPlane ? 1 : 2;

   long double tmp_flux=0.;
   long double bottom_flux=0.;

   // First, fill the y/z=3 cells
   // Then integrate in y/z direction
   for (int x = B.dimension[0]->cells - (outerBoundary+1); x >= outerBoundary; x--) {
      int i = outerBoundary;
      int y = eqPlane ? i : 0;
      int z = eqPlane ? 0 : i;
      if (isInside(B, innerBoundary, x, y, z)) {
         break;
      }

      Vec3d bval = B.getCell(x,y,z);

      bottom_flux -= bval[yCoord] * B.dx[0];
      flux[B.dimension[0]->cells * i + x] = bottom_flux;

      tmp_flux = bottom_flux;
      for(i++; i < B.dimension[yCoord]->cells - outerBoundary; i++) {
         y = eqPlane ? i : 0;
         z = eqPlane ? 0 : i;
         if (isInside(B, innerBoundary, x, y, z)) {
            break;
         }

         bval = B.getCell(x,y,z);

         tmp_flux -= bval[0] * B.dx[yCoord];
         flux[B.dimension[0]->cells * i + x] = tmp_flux;
      }
   }

   return flux;
}

// Calculate fluxfunction by integrating along +y/z boundary first,
// and then going along negative y/z-direction.
std::vector<double> computeFluxDown(Field& B, int outerBoundary, double innerBoundary) {
   // Create fluxfunction-field to be the same shape as B
   std::vector<double> flux(B.dimension[0]->cells * B.dimension[1]->cells * B.dimension[2]->cells);
   for (int i = 0; i < flux.size(); ++i) {
      flux[i] = NAN;
   }

   bool eqPlane = B.dimension[1]->cells > 1;
   int yCoord = eqPlane ? 1 : 2;

   long double tmp_flux=0.;
   long double top_flux=0.;

   // Calculate flux-difference between bottom and top edge
   // of +x boundary (so that values are consistent with computeFluxUp)
   for(int i = outerBoundary; i < B.dimension[yCoord]->cells - outerBoundary; i++) {
      int x = B.dimension[0]->cells - (outerBoundary + 1);
      int y = eqPlane ? i : 0;
      int z = eqPlane ? 0 : i;
      if (isInside(B, innerBoundary, x, y, z)) {
         break;
      }

      Vec3d bval = B.getCell(x, y, z);

      top_flux -= bval[0]*B.dx[yCoord];
      flux[B.dimension[0]->cells * i + x] = top_flux;
   }

   // First, fill the y/z = max - 4 cells
   // Then integrate in -y/z direction
   for(int x = B.dimension[0]->cells - (outerBoundary + 2); x >= outerBoundary; x--) {
      int i = B.dimension[yCoord]->cells - (outerBoundary + 1);
      int y = eqPlane ? i : 0;
      int z = eqPlane ? 0 : i;
      if (isInside(B, innerBoundary, x, y, z)) {
         break;
      }

      Vec3d bval = B.getCell(x, y, z);

      top_flux -= bval[yCoord] * B.dx[0];
      flux[B.dimension[0]->cells * i + x] = top_flux;

      tmp_flux = top_flux;
      for (i--; i >= outerBoundary; i--) {
         y = eqPlane ? i : 0;
         z = eqPlane ? 0 : i;
         if (isInside(B, innerBoundary, x, y, z)) {
            break;
         }

         bval = B.getCell(x, y, z);

         tmp_flux += bval[0] * B.dx[yCoord];
         flux[B.dimension[0]->cells * i + x] = tmp_flux;
      }
   }

   return flux;
}

// Calculate fluxfunction by integrating along -x from the right boundary
std::vector<double> computeFluxLeft(Field& B, int outerBoundary, double innerBoundary) {
   // Create fluxfunction-field to be the same shape as B
   std::vector<double> flux(B.dimension[0]->cells * B.dimension[1]->cells * B.dimension[2]->cells);
   for (int i = 0; i < flux.size(); ++i) {
      flux[i] = NAN;
   }

   bool eqPlane = B.dimension[1]->cells > 1;
   int yCoord = eqPlane ? 1 : 2;

   long double tmp_flux=0.;
   long double right_flux=0.;

   // First calculate flux difference to bottom right corner
   // Now, for each row, integrate in -z-direction.
   for(int i = outerBoundary; i < B.dimension[yCoord]->cells - outerBoundary; i++) {
      int x = B.dimension[0]->cells - (outerBoundary + 1);
      int y = eqPlane ? i : 0;
      int z = eqPlane ? 0 : i;
      if (isInside(B, innerBoundary, x, y, z)) {
         break;
      }

      Vec3d bval = B.getCell(x, y, z);

      right_flux -= bval[0] * B.dx[yCoord];
      flux[B.dimension[0]->cells * i + x] = right_flux;

      tmp_flux = right_flux;
      for(x--; x >= outerBoundary; x--) {
         if (isInside(B, innerBoundary, x, y, z)) {
            break;
         }

         bval = B.getCell(x,y,z);

         tmp_flux -= bval[yCoord] * B.dx[0];
         flux[B.dimension[0]->cells * i + x] = tmp_flux;
      }
   }

   return flux;
}

// Calculate fluxfunction by integrating along -y/z boundary
// Then along the -x boundary
// And finally right in the +x direction
std::vector<double> computeFluxUpRight(Field& B, int outerBoundary, double innerBoundary) {
   // Create fluxfunction-field to be the same shape as B
   std::vector<double> flux(B.dimension[0]->cells * B.dimension[1]->cells * B.dimension[2]->cells);
   for (int i = 0; i < flux.size(); ++i) {
      flux[i] = NAN;
   }

   bool eqPlane = B.dimension[1]->cells > 1;
   int yCoord = eqPlane ? 1 : 2;

   long double tmp_flux=0.;
   long double left_flux=0.;

   // First calculate flux difference from the right edge
   for (int x = B.dimension[0]->cells - (outerBoundary + 1); x >= outerBoundary; x--) {
      int i = outerBoundary;
      int y = eqPlane ? i : 0;
      int z = eqPlane ? 0 : i;
      if (isInside(B, innerBoundary, x, y, z)) {
         break;
      }

      Vec3d bval = B.getCell(x,y,z);

      left_flux -= bval[yCoord] * B.dx[0];
      flux[B.dimension[0]->cells * i + x] = left_flux;
   }

   // Now, for each row, integrate in y/z-direction,
   // Then integrate in +x direction
   for (int i = outerBoundary; i < B.dimension[yCoord]->cells - outerBoundary; i++) {
      int x = outerBoundary;
      int y = eqPlane ? i : 0;
      int z = eqPlane ? 0 : i;
      if (isInside(B, innerBoundary, x, y, z)) {
         break;
      }

      Vec3d bval = B.getCell(x, y, z);

      left_flux -= bval[0] * B.dx[yCoord];
      flux[B.dimension[0]->cells * i + x] = left_flux;

      tmp_flux = left_flux;
      for(x++; x < B.dimension[0]->cells - outerBoundary; x++) {
         if (isInside(B, innerBoundary, x, y, z)) {
            break;
         }

         bval = B.getCell(x,y,z);

         tmp_flux += bval[yCoord] * B.dx[0];
         flux[B.dimension[0]->cells * i + x] = tmp_flux;
      }
   }

   return flux;
}

// Calculate fluxfunction by integrating along +y/z boundary
// Then along the -x boundary
// And finally right in the +x direction
std::vector<double> computeFluxDownRight(Field& B, int outerBoundary, double innerBoundary) {
   // Create fluxfunction-field to be the same shape as B
   std::vector<double> flux(B.dimension[0]->cells * B.dimension[1]->cells * B.dimension[2]->cells);
   for (int i = 0; i < flux.size(); ++i) {
      flux[i] = NAN;
   }

   bool eqPlane = B.dimension[1]->cells > 1;
   int yCoord = eqPlane ? 1 : 2;

   long double tmp_flux=0.;
   long double left_flux=0.;

   // Calculate flux-difference between bottom and top edge
   // of +x boundary (so that values are consistent with computeFluxUp)
   for(int i = outerBoundary; i < B.dimension[yCoord]->cells - outerBoundary; i++) {
      int x = B.dimension[0]->cells - (outerBoundary + 1);
      int y = eqPlane ? i : 0;
      int z = eqPlane ? 0 : i;
      if (isInside(B, innerBoundary, x, y, z)) {
         break;
      }

      Vec3d bval = B.getCell(x, y, z);

      left_flux -= bval[0]*B.dx[yCoord];
      flux[B.dimension[0]->cells * i + x] = left_flux;
   }

   // Then to left edge
   for(int x = B.dimension[0]->cells - (outerBoundary + 2); x >= outerBoundary; x--) {
      int i = B.dimension[yCoord]->cells - (outerBoundary + 1);
      int y = eqPlane ? i : 0;
      int z = eqPlane ? 0 : i;
      if (isInside(B, innerBoundary, x, y, z)) {
         break;
      }

      Vec3d bval = B.getCell(x, y, z);

      left_flux -= bval[yCoord]*B.dx[0];
      flux[B.dimension[0]->cells * i + x] = left_flux;
   }

   // Now, for each row, integrate in y/z-direction,
   // Then integrate in +x direction
   for (int i = B.dimension[yCoord]->cells - (outerBoundary + 2); i >= outerBoundary; i--) {
      int x = outerBoundary;
      int y = eqPlane ? i : 0;
      int z = eqPlane ? 0 : i;
      if (isInside(B, innerBoundary, x, y, z)) {
         break;
      }

      Vec3d bval = B.getCell(x, y, z);

      left_flux += bval[0] * B.dx[yCoord];
      flux[B.dimension[0]->cells * i + x] = left_flux;

      tmp_flux = left_flux;
      for(x++; x < B.dimension[0]->cells - outerBoundary; x++) {
         if (isInside(B, innerBoundary, x, y, z)) {
            break;
         }

         bval = B.getCell(x,y,z);

         tmp_flux += bval[yCoord] * B.dx[0];
         flux[B.dimension[0]->cells * i + x] = tmp_flux;
      }
   }

   return flux;
}

// Get median of vector
double nanMedian(std::vector<double> &v) {
   v.erase(std::remove_if(v.begin(), v.end(), [](const double& value) {return !isfinite(value);}), v.end());
   int n = v.size();
   if (!n) {
      return NAN;
   } else if (n % 2) {
      std::nth_element(v.begin(), v.begin() + n/2, v.end());
      return v[n/2];
   } else {
      std::nth_element(v.begin(), v.begin() + n/2 - 1, v.end()); // Left median
      std::nth_element(v.begin(), v.begin() + n/2, v.end());     // Right median
      return (v[n/2 - 1] + v[n/2])/2;
   }
}

// Get mean of vector
double nanMean(std::vector<double> &v) {
   v.erase(std::remove_if(v.begin(), v.end(), [](const double& value) {return !isfinite(value);}), v.end());
   int n = v.size();
   return n ? std::accumulate(v.begin(), v.end(), 0.0) / n : NAN;
}

int main(int argc, char** argv) {

   // MPI::Init(argc, argv);

   if(argc < 3) {
      cerr << "Syntax: fluxfunction input.vlsv output.bin" << endl;
      cerr << "Output will be two files: output.bin and output.bin.bov." << endl;
      cerr << "Point visit to the BOV file." << endl;
      return 1;
   }
   string inFile(argv[1]);
   string outFile(argv[2]);

   // TODO: Don't uselessly read E, we really only care about B.
   Field E,B,V;
   readfields(inFile.c_str(),E,B,V,false);

   // Make sure we are working with a 2D simulation here.
   if(B.dimension[0]->cells > 1 && B.dimension[1]->cells > 1 && B.dimension[2]->cells > 1) {
         cerr << "This is a 3D simulation output. Flux function calculation only makes sense for 2D data."
            << endl;
         exit(1);
   }

   cerr << "File read, calculating flux function..." << endl;

   double dOuter = 2;
   double rInner = 5*physicalconstants::R_E;  // Default for now
   std::vector<double> fluxUp, fluxDown, fluxLeft, fluxUR, fluxDR;
   fluxUp = computeFluxUp(B, dOuter, rInner);
   fluxDown = computeFluxDown(B, dOuter, rInner);
   fluxLeft = computeFluxLeft(B, dOuter, rInner);
   fluxUR = computeFluxUpRight(B, dOuter, rInner);
   fluxDR = computeFluxDownRight(B, dOuter, rInner);

   for(unsigned int i=0; i<fluxUp.size(); i++) {
      std::vector<double> v {fluxUp[i], fluxDown[i], fluxLeft[i], fluxUR[i], fluxDR[i]};
      fluxUp[i] = nanMedian(v);
      //fluxDown[i] = nanMedian(v);
      //fluxLeft[i] = nanMean(v);
      //fluxUR[i] = abs((fluxDown[i] - fluxLeft[i])/fluxDown[i]);
      //fluxDR[i] = abs((fluxDown[i] - fluxLeft[i])/fluxLeft[i]);
   }

   cerr << "Done. Writing output..." << endl;
   // std::cout << "Mean of medians " << nanMean(fluxDown) << std::endl;
   // std::cout << "Mean of means " << nanMean(fluxLeft) << std::endl;
   // std::cout << "Mean relative difference " << nanMean(fluxUR) << std::endl;
   // std::cout << "Mean relative difference " << nanMean(fluxDR) << std::endl;
   // std::cout << "Maximum difference " << *max_element(fluxUR.begin(), fluxUR.end()) << std::endl;
   // std::cout << "Maximum difference " << *max_element(fluxDR.begin(), fluxDR.end()) << std::endl;

   // Write output as a visit-compatible BOV file
   int fd = open(outFile.c_str(), O_CREAT|O_TRUNC|O_WRONLY, 0644);
   if(!fd || fd == -1) {
      cerr << "Error: cannot open output file " << outFile << ": " << strerror(errno) << endl;
      return 1;
   }
   size_t size=B.dimension[0]->cells*B.dimension[1]->cells*B.dimension[2]->cells*sizeof(double);

   // Write binary blob
   for(ssize_t remain=size; remain > 0; ) {
      remain -= write(fd, ((char*) &(fluxUp[0]))+remain-size, remain);
   }
   close(fd);

   // Write BOV header
   string outBov = outFile + ".bov";
   FILE* f=fopen(outBov.c_str(), "w");
   if(!f) {
      cerr<< "Error: unable to write BOV ascii file " << outBov << ":" << strerror(errno) << endl;
      return 1;
   }
   fprintf(f, "TIME: %lf\n", B.time);
   fprintf(f, "DATA_FILE: %s\n", outFile.c_str());
   fprintf(f, "DATA_SIZE: %i %i %i\n", B.dimension[0]->cells, B.dimension[1]->cells, B.dimension[2]->cells);
   fprintf(f, "DATA_FORMAT: DOUBLE\nVARIABLE: fluxfunction\nDATA_ENDIAN: LITTLE\nCENTERING: zonal\n");
   fprintf(f, "BRICK_ORIGIN: %lf %lf %lf\n", B.dimension[0]->min, B.dimension[1]->min, B.dimension[2]->min);
   fprintf(f, "BRICK_SIZE: %lf %lf %lf\n", B.dimension[0]->max - B.dimension[0]->min, B.dimension[1]->max - B.dimension[1]->min, B.dimension[2]->max - B.dimension[2]->min);
   fprintf(f, "DATA_COMPONENTS: 1\n");

   fclose(f);

   return 0;
}
