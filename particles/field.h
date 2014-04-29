#pragma once
#include <glm/glm.hpp>
#include <vector>

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
		if(x > cells[0] || y > cells[1] || z > cells[2] 
				|| x <0 || y < 0 || z < 0) {
			std::cerr << "Field access out of bounds!" << std::endl;
		}

		double* cell = &(data[4*(z*cells[0]*cells[1] + y*cells[0] + x)]);
		return cell;
	}

	glm::dvec3 getCell(int x, int y, int z) {
		double* cell = getCellRef(x,y,z);
		glm::dvec3 retval;
		retval.x = cell[0];
		retval.y = cell[1];
		retval.z = cell[2];

		return retval;
	}

	/* Round-Brace indexing: indexing by physical location, with interpolation */
	glm::dvec3 operator()(double x, double y, double z) {
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
	glm::dvec3 operator()(glm::dvec3 pos) {
		return operator()(pos.x,pos.y,pos.z);
	}

};

