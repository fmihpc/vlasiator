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

// Calculate fluxfunction by integrating along -z boundary first,
// and then going along z-direction.
std::vector<double> computeFluxUp(Field& B) {
  // Create fluxfunction-field to be the same shape as B
  std::vector<double> flux(B.cells[0] * B.cells[1] * B.cells[2]);

  long double tmp_flux=0.;

  // First, fill the z=3 cells
  for(int x=B.cells[0]-2; x>0; x--) {
    Vec3d bval = B.getCell(x,0,3);

    tmp_flux -= bval[2] * B.dx[0];
    flux[B.cells[0] * B.cells[1] * 3 + x] = tmp_flux;
  }

  // Now, for each row, integrate in z-direction.
  for(int x=1; x< B.cells[0]-1; x++) {

    tmp_flux = flux[B.cells[0] * B.cells[1] * 3 + x];
    for(int z=4; z< B.cells[2]; z++) {
      Vec3d bval = B.getCell(x,0,z);

      tmp_flux -= bval[0]*B.dx[2];
      flux[B.cells[0] * B.cells[1] * z  +  x] = tmp_flux;
    }
  }

  return flux;
}



// Calculate fluxfunction by integrating along +z boundary first,
// and then going along negative z-direction.
std::vector<double> computeFluxDown(Field& B) {
  // Create fluxfunction-field to be the same shape as B
  std::vector<double> flux(B.cells[0] * B.cells[1] * B.cells[2]);

  long double tmp_flux=0.;

  // Calculate flux-difference between bottom and top edge
  // of +x boundary (so that values are consistent with computeFluxUp)
  for(int z=3; z<B.cells[2]-4; z++) {
    Vec3d bval = B.getCell(B.cells[0]-2,0,z);

    tmp_flux -= bval[0]*B.dx[2];
  }

  // From there, fill the z=max - 4 cells
  for(int x=B.cells[0]-2; x>0; x--) {
    Vec3d bval = B.getCell(x,0,B.cells[2]-4);

    tmp_flux -= bval[2] * B.dx[0];
    flux[B.cells[0] * B.cells[1] * (B.cells[2] - 4) + x] = tmp_flux;
  }

  // Now, for each row, integrate in -z-direction.
  for(int x=1; x< B.cells[0]-1; x++) {

    tmp_flux = flux[B.cells[0] * B.cells[1] * (B.cells[2] - 4) + x];
    for(int z=B.cells[2]-5; z > 0; z--) {
      Vec3d bval = B.getCell(x,0,z);

      tmp_flux += bval[0] * B.dx[2];
      flux[B.cells[0] * B.cells[1] * z  +  x] = tmp_flux;
    }
  }

  return flux;
}



// Calculate fluxfunction by integrating along -x from the right boundary
std::vector<double> computeFluxLeft(Field& B) {
  // Create fluxfunction-field to be the same shape as B
  std::vector<double> flux(B.cells[0] * B.cells[1] * B.cells[2]);

  long double tmp_flux=0.;
  long double bottom_right_flux=0.;

  // First calculate flux difference to bottom right corner
  // Now, for each row, integrate in -z-direction.
  for(int z=0; z < B.cells[2]; z++) {
    Vec3d bval = B.getCell(B.cells[0]-1,0,z);
    bottom_right_flux -= bval[0] * B.dx[2];

    tmp_flux = bottom_right_flux;
    for(int x=B.cells[0]-1; x>0; x--) {

      bval = B.getCell(x,0,z);

      tmp_flux -= bval[2] * B.dx[0];
      flux[B.cells[0] * B.cells[1] * z  +  x] = tmp_flux;
    }
  }

  return flux;
}



// Get a median of 3 values (branch-free!)
static double median3(double a, double b, double c) {
  return max(min(a,b), min(max(a,b),c));
}



int main(int argc, char** argv) {

  MPI::Init(argc, argv);

  if(argc < 3) {
    cerr << "Syntax: fluxfunction input.vlsv output.bin" << endl;
    cerr << "Output will be two files: output.bin and output.bin.bov." << endl;
    cerr << "Point visit to the BOV file." << endl;
    return 1;
  }
  string inFile(argv[1]);
  string outFile(argv[2]);

  // TODO: Don't uselessly read E and V, we really only care about B.
  Field E,B,V;
  readfields(inFile.c_str(),E,B,V);

  cerr << "File read, calculating flux function..." << endl;

  std::vector<double> fluxUp(computeFluxUp(B)), fluxDown(computeFluxDown(B)), fluxLeft(computeFluxLeft(B));

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
  size_t size=B.cells[0]*B.cells[1]*B.cells[2]*sizeof(double);

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
  fprintf(f, "DATA_SIZE: %i %i %i\n", B.cells[0], B.cells[1], B.cells[2]);
  fprintf(f, "DATA_FORMAT: DOUBLE\nVARIABLE: fluxfunction\nDATA_ENDIAN: LITTLE\nCENTERING: zonal\n");
  fprintf(f, "BRICK_ORIGIN: %lf %lf %lf\n", B.min[0], B.min[1], B.min[2]);
  fprintf(f, "BRICK_SIZE: %lf %lf %lf\n", B.max[0] - B.min[0], B.max[1] - B.min[1], B.max[2] - B.min[2]);
  fprintf(f, "DATA_COMPONENTS: 1\n");

  fclose(f);

  return 0;
}
