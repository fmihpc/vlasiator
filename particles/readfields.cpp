#include <string>
#include <vector>
#include <iostream>
#include <stdint.h>
#include "field.h"
#include "readfields.h"

/* Debugging image output */
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"


/* For debugging purposes - dump a field into a png file
 * We're hardcodedly writing the z=0 plane here. */
void debug_output(Field& F, const char* filename) {
	
	/* Find min and max value */
	glm::dvec3 min, max;

	/* TODO: Darn, this is dirty. */
	min[0] = min[1] = min[2] = 99999999999;
	max[0] = max[1] = max[2] = -99999999999;

	for(int i=0; i<F.cells[0]*F.cells[1]; i++) {
		for(int j=0; j<3; j++) {
			if(F.data[4*i+j] > max[j]) {
				max[j] = F.data[4*i+j];
			}
			if(F.data[4*i+j] < min[j]) {
				min[j] = F.data[4*i+j];
			}
		}
	}

	/* Allocate a rgb-pixel array */
	std::vector<uint8_t> pixels(4*F.cells[0]*F.cells[1]);

	/* And fill it with colors */
	for(int y=0; y<F.cells[1]; y++) {
		for(int x=0; x<F.cells[0]; x++) {

			/* Rescale the field values to lie between 0..255 */
			glm::dvec3 scaled_val = F.getCell(x,y,0);
			scaled_val -= min;
			scaled_val /= (max-min);
			scaled_val *= 255.;

			pixels[4*(y*F.cells[0] + x)] = (uint8_t) scaled_val.x;
			pixels[4*(y*F.cells[0] + x)+1] = (uint8_t) scaled_val.y;
			pixels[4*(y*F.cells[0] + x)+2] = (uint8_t) scaled_val.z;
			pixels[4*(y*F.cells[0] + x)+3] = 255; // Alpha=1
		}
	}

	/* Write it out */
	if(!stbi_write_png(filename, F.cells[0], F.cells[1], 4, pixels.data(), F.cells[0]*4)) {
		std::cerr << "Writing " << filename << " failed: " << strerror(errno) << std::endl;
	}
}

