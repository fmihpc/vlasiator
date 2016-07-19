#include <string>
#include <vector>
#include <iostream>
#include <stdint.h>
#include "field.h"
#include "readfields.h"
#include "vectorclass.h"
#include "vector3d.h"
#include "../definitions.h"

/* Debugging image output */
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

std::string B_field_name;
std::string E_field_name;

/* Read the cellIDs into an array */
std::vector<uint64_t> readCellIds(vlsvinterface::Reader& r) {

   uint64_t arraySize=0;
   uint64_t vectorSize=0;
   uint64_t byteSize=0;
   vlsv::datatype::type dataType;
   std::list<std::pair<std::string,std::string> > attribs;
   attribs.push_back(std::pair<std::string,std::string>("name","CellID"));
   if( r.getArrayInfo("VARIABLE",attribs, arraySize,vectorSize,dataType,byteSize) == false ) {
      std::cerr << "getArrayInfo returned false when trying to read CellID VARIABLE." << std::endl;
      exit(1);
   }

   if(dataType != vlsv::datatype::type::UINT || byteSize != 8 || vectorSize != 1) {
      std::cerr << "Datatype of CellID VARIABLE entries is not uint64_t." << std::endl;
      exit(1);
   }

   /* Allocate memory for the cellIds */
   std::vector<uint64_t> cellIds(arraySize*vectorSize);

   if( r.readArray("VARIABLE",attribs,0,arraySize,(char*) cellIds.data()) == false) {
      std::cerr << "readArray faied when trying to read CellID Variable." << std::endl;
      exit(1);
   }

   return cellIds;
}

/* For debugging purposes - dump a field into a png file
 * We're hardcodedly writing the z=0 plane here. */
void debug_output(Field& F, const char* filename) {
   
   /* Find min and max value */
   Real min[3], max[3];

   /* TODO: ugh, this is an ungly hack */
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

   Vec3d vmin,vmax;

   vmin.load(min);
   vmax.load(max);

   /* Allocate a rgb-pixel array */
   std::vector<uint8_t> pixels(4*F.cells[0]*F.cells[1]);

   /* And fill it with colors */
   for(int y=0; y<F.cells[1]; y++) {
      for(int x=0; x<F.cells[0]; x++) {

         /* Rescale the field values to lie between 0..255 */
         Vec3d scaled_val = F.getCell(x,y,0);
         scaled_val -= vmin;
         scaled_val /= (vmax-vmin);
         scaled_val *= 255.;

         pixels[4*(y*F.cells[0] + x)] = (uint8_t) scaled_val[0];
         pixels[4*(y*F.cells[0] + x)+1] = (uint8_t) scaled_val[1];
         pixels[4*(y*F.cells[0] + x)+2] = (uint8_t) scaled_val[2];
         pixels[4*(y*F.cells[0] + x)+3] = 255; // Alpha=1
      }
   }

   /* Write it out */
   if(!stbi_write_png(filename, F.cells[0], F.cells[1], 4, pixels.data(), F.cells[0]*4)) {
      std::cerr << "Writing " << filename << " failed: " << strerror(errno) << std::endl;
   }
}

