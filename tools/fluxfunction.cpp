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

int main(int argc, char** argv) {

	MPI::Init(argc, argv);

	if(argc < 3) {
		cerr << "Syntax: fluxfunction input.vlsv output.bin" << endl;
		cerr << "Output will be two files: output.bin and output.bin.bov." << endl;
		cerr << "Point visit to the BOV file." << endl;
		return 1;
	}
	string infile(argv[1]);
	string outfile(argv[2]);

	// TODO: Don't uselessly read E and V, we really only care about V.
	Field E,B,V;
	readfields(infile.c_str(),E,B,V);

	cerr << "File read, calculating flux function..." << endl;

	// Create fluxfunction-field to be the same shape as B
	std::vector<double> flux(B.cells[0] * B.cells[1] * B.cells[2]);

	// Calculate its values.
	// Since B-Values are face-centered, this quantity
	// is volume-centered.
	// TODO: This is assuming we're in the polar plane
	double tmp_flux=0.;

	// First, fill the z=3 cells
	for(int x=B.cells[0]-2; x>0; x--) {
		Vec3d bval = B.getCell(x,0,3);
		
		tmp_flux -= bval[2] * B.dx[0];

		flux[B.cells[0] * B.cells[1] * 3 + x] = tmp_flux;
	}

	// Now, for each row, integrate in x-direction.
	for(int x=1; x< B.cells[0]; x++) {
		for(int z=4; z< B.cells[2]; z++) {
			Vec3d bval = B.getCell(x,0,z);

			// For each cell, check whether integrating in x or in y gives us
			// a smaller absolute flux contribution.
			// (This is to avoid the singularity in earth's centre)
			double up_flux = flux[B.cells[0] * B.cells[1] * (z-1) + x] - bval[0] * B.dx[2];
			//double left_flux = flux[B.cells[0] * B.cells[1] * z + x-1] + bval[2] * B.dx[0];

			//if(fabs(bval[0]) > fabs(bval[2])) {
			//	tmp_flux = up_flux;
			//} else {
			//	tmp_flux = left_flux;
			//}
			
			tmp_flux = up_flux;

			flux[B.cells[0] * B.cells[1] * z  +  x] = tmp_flux;
		}
	}

	cerr << "Done. Writing output..." << endl;

	// Write output as a visit-compatible BOV file
	int fd = open(outfile.c_str(), O_CREAT|O_TRUNC|O_WRONLY, 0644);
	if(!fd || fd == -1) {
		cerr << "Error: cannot open output file " << outfile << ": " << strerror(errno) << endl;
		return 1;
	}
	size_t size=B.cells[0]*B.cells[1]*B.cells[2]*sizeof(double);

	// Write binary blob
	for(ssize_t remain=size; remain > 0; ) {
		remain -= write(fd, ((char*) &(flux[0]))+remain-size, remain);
	}
	close(fd);

	// Write BOV header
	string out_bov = outfile + ".bov";
	FILE* f=fopen(out_bov.c_str(), "w");
	if(!f) {
		cerr<< "Error: unable to write BOV ascii file " << out_bov << ":" << strerror(errno) << endl;
		return 1;
	}
	fprintf(f, "TIME: %lf\n", B.time);
	fprintf(f, "DATA_FILE: %s\n", outfile.c_str());
	fprintf(f, "DATA_SIZE: %i %i %i\n", B.cells[0], B.cells[1], B.cells[2]);
	fprintf(f, "DATA_FORMAT: DOUBLE\nVARIABLE: fluxfunction\nDATA_ENDIAN: LITTLE\nCENTERING: zonal\n");
	fprintf(f, "BRICK_ORIGIN: %lf %lf %lf\n", B.min[0], B.min[1], B.min[2]);
	fprintf(f, "BRICK_SIZE: %lf %lf %lf\n", B.max[0] - B.min[0], B.max[1] - B.min[1], B.max[2] - B.min[2]);
	fprintf(f, "DATA_COMPONENTS: 1\n");

	fclose(f);

	return 0;
}
