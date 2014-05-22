#pragma once
#include <vector>
#include "vectorclass.h"
#include "vector3d.h"

/* A 3D cartesian vector field with suitable interpolation properties for
 * particle pushing */
struct Field {
	/* Coordinate boundaries */
	double min[3];
	double max[3];

	/* Mesh spacing */
	double dx[3];

	/* Mesh cells */
	int cells[3];

	bool periodic[3];

	/* The actual field data */
	std::vector<double> data;

	double* getCellRef(int x, int y, int z) {

		if(cells[2] == 1) {
			return &(data[4*(y*cells[0]+x)]);
		} else {
			return &(data[4*(z*cells[0]*cells[1] + y*cells[0] + x)]);
		}
	}

	Vec3d getCell(int x, int y, int z) {
		double* cell = getCellRef(x,y,z);
		Vec3d retval;
		retval.load_partial(3,cell);
		return retval;
	}

	/* Round-Brace indexing: indexing by physical location, with interpolation */
	Vec3d operator()(double x, double y, double z) {
		Vec3d v(x,y,z);
		Vec3d vmin,vdx;
		vmin.load(min);
		vdx.load(dx);

		v -= vmin;
		v /= vdx;

		int index[3];
		double fract[3];
		round_to_int(v).store(index);
		fraction(v).store(fract);

		/* TODO: Boundary behaviour? */

		/* TODO: 3d? */
		Vec3d interp[4];
		interp[0] = getCell(index[0],index[1],index[2]);
		interp[1] = getCell(index[0]+1,index[1],index[2]);
		interp[2] = getCell(index[0],index[1]+1,index[2]);
		interp[3] = getCell(index[0]+1,index[1]+1,index[2]);

		return fract[0]*(fract[1]*interp[3]+(1.-fract[1])*interp[1])
		+ (1.-fract[0])*(fract[1]*interp[2]+(1.-fract[1])*interp[0]);
	}
	Vec3d operator()(Vec3d pos) {
		return operator()(pos[0],pos[1],pos[2]);
	}

};

