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

	/* The actual field data */
	std::vector<double> data;

	double* getCellRef(int x, int y, int z) {
		/* Bounds checking */
		//if(x > cells[0] || y > cells[1] || z > cells[2] 
		//		|| x <0 || y < 0 || z < 0) {
		//	std::cerr << "Field access out of bounds!" << std::endl;
		//}

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
		/* TODO: Vectorize these */
		x -= min[0];
		y -= min[1];
		z -= min[2];

		x/=dx[0];
		y/=dx[1];
		z/=dx[2];

		/* TODO: Boundary behaviour? */

		/* TODO: interpolation */
		return getCell((int)x,(int)y,(int)z);
	}
	Vec3d operator()(Vec3d pos) {
		return operator()(pos[0],pos[1],pos[2]);
	}

};

