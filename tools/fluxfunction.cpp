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
#include "particles/field.h"
#include "particles/readfields.h"

using namespace std;

// Calculate fluxfunction by integrating along -y/z boundary first,
// and then going along y/z-direction.
std::vector<double> computeFluxUp(Field& B, int outerBoundary, int innerBoundaryRE) {
   // Create fluxfunction-field to be the same shape as B
   std::vector<double> flux(B.dimension[0]->cells * B.dimension[1]->cells * B.dimension[2]->cells);
   bool eqPlane = B.dimension[1]->cells > 1;
   int yCoord = eqPlane ? 1 : 2;

   long double tmp_flux=0.;
   long double bottom_flux=0.;

   // First, fill the y/z=3 cells
   // Then integrate in y/z direction
   for(int x = B.dimension[0]->cells - (outerBoundary+1); x > outerBoundary + 1; x--) {
      int y = eqPlane ? outerBoundary : 0;
      int z = eqPlane ? 0 : outerBoundary;
      Vec3d bval = B.getCell(x,y,z);

      bottom_flux -= bval[yCoord] * B.dx[0];
      flux[B.dimension[0]->cells * outerBoundary + x] = tmp_flux;

      tmp_flux = bottom_flux
      for(int i=outerBoundary + 1; i < B.dimension[yCoord]->cells - outerBoundary; i++) {
         y = eqPlane ? i : 0;
         z = eqPlane ? 0 : i;
         bval = B.getCell(x,y,z);

         tmp_flux -= bval[0] * B.dx[yCoord];
         flux[B.dimension[0]->cells * i + x] = tmp_flux;
      }
   }

   return flux;
}

// Calculate fluxfunction by integrating along +y/z boundary first,
// and then going along negative y/z-direction.
std::vector<double> computeFluxDown(Field& B, int outerBoundary, int innerBoundaryRE) {
   // Create fluxfunction-field to be the same shape as B
   std::vector<double> flux(B.dimension[0]->cells * B.dimension[1]->cells * B.dimension[2]->cells);
   bool eqPlane = B.dimension[1]->cells > 1;
   int yCoord = eqPlane ? 1 : 2;

   long double tmp_flux=0.;
   long double top_flux=0.;

   // Calculate flux-difference between bottom and top edge
   // of +x boundary (so that values are consistent with computeFluxUp)
   for(int i = outerBoundary; i < B.dimension[yCoord]->cells - outerBoundary; i++) {
      int y = eqPlane ? i : 0;
      int z = eqPlane ? 0 : i;
      Vec3d bval = B.getCell(B.dimension[0]->cells - (outerBoundary + 1), y, z);

      tmp_flux -= bval[0]*B.dx[yCoord];
   }

   // First, fill the y/z = max - 4 cells
   // Then integrate in -y/z direction
   for(int x = B.dimension[0]->cells - (outerBoundary + 2); x >= outerBoundary; x--) {
      int y = eqPlane ? B.dimension[1]->cells - (outerBoundary + 1) : 0;
      int z = eqPlane ? 0 : B.dimension[2]->cells - (outerBoundary + 1);
      Vec3d bval = B.getCell(x, y, z);

      top_flux -= bval[yCoord] * B.dx[0];
      flux[B.dimension[0]->cells * (B.dimension[yCoord]->cells - (outerBoundary + 1)) + x] = top_flux;

      tmp_flux = top_flux
      for(int i=B.dimension[yCoord]->cells - (outerBoundary + 2); i > 0; i--) {
         y = eqPlane ? i : 0;
         z = eqPlane ? 0 : i;
         bval = B.getCell(x, y, z);

         tmp_flux += bval[0] * B.dx[yCoord];
         flux[B.dimension[0]->cells * i + x] = tmp_flux;
      }
   }

   return flux;
}

// Calculate fluxfunction by integrating along -x from the right boundary
std::vector<double> computeFluxLeft(Field& B, int outerBoundary, int innerBoundaryRE) {
   // Create fluxfunction-field to be the same shape as B
   std::vector<double> flux(B.dimension[0]->cells * B.dimension[1]->cells * B.dimension[2]->cells);
   bool eqPlane = B.dimension[1]->cells > 1;
   int yCoord = eqPlane ? 1 : 2;

   long double tmp_flux=0.;
   long double right_flux=0.;

   // First calculate flux difference to bottom right corner
   // Now, for each row, integrate in -z-direction.
   for(int i=0; i < B.dimension[yCoord]->cells - outerBoundary; i++) {
      int y = eqPlane ? i : 0;
      int z = eqPlane ? 0 : i;
      Vec3d bval = B.getCell(B.dimension[0]->cells - (outerBoundary + 1), y, z);

      right_flux -= bval[0] * B.dx[yCoord];
      flux[B.dimension[0]->cells * i + B.dimension[0]->cells - (outerboundary + 1)] = tmp_flux;

      tmp_flux = right_flux;
      for(int x = B.dimension[0]->cells - (outerBoundary + 2); x >= outerBoundary; x--) {

         bval = B.getCell(x,y,z);

         tmp_flux -= bval[yCoord] * B.dx[0];
         flux[B.dimension[0]->cells * i + x] = tmp_flux;
      }
   }

   return flux;
}

// Get a median of 3 values (branch-free!)
static double median3(double a, double b, double c) {
  return max(min(a,b), min(max(a,b),c));
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


  std::vector<double> fluxUp, fluxDown, fluxLeft;
   fluxUp = computeFluxUp(B, 2, 5);
   fluxDown = computeFluxDown(B, 2, 5);
   fluxLeft = computeFluxLeft(B, 2, 5);

  for(unsigned int i=0; i<fluxUp.size(); i++) {
    // Calc median flux value;
    double a = fluxUp[i];
    double b = fluxDown[i];
    double c = fluxLeft[i];

    fluxUp[i] = median3(a,b,c);
  }

  cerr << "Done. Writing output..." << endl;

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
