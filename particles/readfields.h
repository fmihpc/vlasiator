#pragma once

#include "vlsv_reader.h"
#include "vlsvreaderinterface.h"
#include "field.h"
#include <algorithm>
#include <vector>
#include <string>
#include <set>

#define DEBUG

extern std::string B_field_name;
extern std::string E_field_name;

/* Read the cellIDs into an array */
std::vector<uint64_t> readCellIds(vlsvinterface::Reader& r);

template <class Reader>
void detect_field_names(Reader& r) {

#ifdef DEBUG
   std::cerr << "Checking for volume-averaged fields... ";
#endif
   std::list<std::string> variableNames;
   std::string gridname("SpatialGrid");

   r.getVariableNames(gridname,variableNames);
   if(find(variableNames.begin(), variableNames.end(), std::string("B_vol"))!=variableNames.end()) {
#ifdef DEBUG
      std::cerr << "yep!" << std::endl;
#endif
      B_field_name = "B_vol";
      E_field_name = "E_vol";
   } else if(find(variableNames.begin(), variableNames.end(), std::string("B"))!=variableNames.end()) {
#ifdef DEBUG
      std::cerr << "Nope!" << std::endl;
#endif
      B_field_name = "B";
      E_field_name = "E";
   } else {
      std::cerr << "No B- or E-fields found! Strange file format?" << std::endl;
      exit(1);
   }
}

/* Read the "raw" field data in file order */
template <class Reader>
std::vector<double> readFieldData(Reader& r, std::string& name, unsigned int numcomponents) {

   uint64_t arraySize=0;
   uint64_t vectorSize=0;
   uint64_t byteSize=0;
   vlsv::datatype::type dataType;
   std::list<std::pair<std::string,std::string> > attribs;
   attribs.push_back(std::pair<std::string,std::string>("name",name));
   if( r.getArrayInfo("VARIABLE",attribs, arraySize,vectorSize,dataType,byteSize) == false ) {
      std::cerr << "getArrayInfo returned false when trying to read VARIABLE \""
         << name << "\"." << std::endl;
      exit(1);
   }

   if(dataType != vlsv::datatype::type::FLOAT || byteSize != 8 || vectorSize != numcomponents) {
      std::cerr << "Datatype of VARIABLE \"" << name << "\" entries is not double." << std::endl;
      exit(1);
   }

   /* Allocate memory for the data */
   std::vector<double> buffer(arraySize*vectorSize);
   
   if( r.readArray("VARIABLE",attribs,0,arraySize,(char*) buffer.data()) == false) {
      std::cerr << "readArray faied when trying to read VARIABLE \"" << name << "\"." << std::endl;
      exit(1);
   }

   return buffer;
}

/* Read the next logical input file. Depending on sign of dt,
 * this may be a numerically larger or smaller file.
 * Return value: true if a new file was read, otherwise false.
 */
template <class Reader>
bool read_next_timestep(const std::string& filename_pattern, double t, int step, Field& E0, Field& E1,
   Field& B0, Field& B1, Field& V, bool doV, int& input_file_counter) {

   char filename_buffer[256];
   bool retval = false;

   while(t < E0.time || t>= E1.time) {
      input_file_counter += step;

      E0=E1;
      B0=B1;
      snprintf(filename_buffer,256,filename_pattern.c_str(),input_file_counter);

      /* Open next file */
      Reader r;
      r.open(filename_buffer);
      double t;
      if(!r.readParameter("time",t)) {
        if(!r.readParameter("t",t)) {
          std::cerr << "Time parameter in file " << filename_buffer << " is neither 't' nor 'time'. Bad file format?" << std::endl;
          exit(1);
        }
      }

      E1.time = t;
      B1.time = t;

      uint64_t cells[3];
      r.readParameter("xcells_ini",cells[0]);
      r.readParameter("ycells_ini",cells[1]);
      r.readParameter("zcells_ini",cells[2]);

      /* Read CellIDs and Field data */
      std::vector<uint64_t> cellIds = readCellIds(r);
      std::string name(B_field_name);
      std::vector<double> Bbuffer = readFieldData(r,name,3u);
      name = E_field_name;
      std::vector<double> Ebuffer = readFieldData(r,name,3u);
			std::vector<double> Vbuffer;
			if(doV) {
				name = "rho_v";
				std::vector<double> rho_v_buffer = readFieldData(r,name,3u);
				name = "rho";
				std::vector<double> rho_buffer = readFieldData(r,name,1u);
				for(int i=0; i<rho_buffer.size(); i++) {
					Vbuffer.push_back(rho_v_buffer[3*i] / rho_buffer[i]);
					Vbuffer.push_back(rho_v_buffer[3*i+1] / rho_buffer[i]);
					Vbuffer.push_back(rho_v_buffer[3*i+2] / rho_buffer[i]);
				}
			}

      /* Assign them, without sanity checking */
      /* TODO: Is this actually a good idea? */
      for(uint i=0; i< cellIds.size(); i++) {
         uint64_t c = cellIds[i];
         int64_t x = c % cells[0];
         int64_t y = (c /cells[0]) % cells[1];
         int64_t z = c /(cells[0]*cells[1]);

         double* Etgt = E1.getCellRef(x,y,z);
         double* Btgt = B1.getCellRef(x,y,z);
         Etgt[0] = Ebuffer[3*i];
         Etgt[1] = Ebuffer[3*i+1];
         Etgt[2] = Ebuffer[3*i+2];
         Btgt[0] = Bbuffer[3*i];
         Btgt[1] = Bbuffer[3*i+1];
         Btgt[2] = Bbuffer[3*i+2];

				 if(doV) {
				   double* Vtgt = V.getCellRef(x,y,z);
					 Vtgt[0] = Vbuffer[3*i];
					 Vtgt[1] = Vbuffer[3*i+1];
					 Vtgt[2] = Vbuffer[3*i+2];
				 }
      }

      r.close();
      retval = true;
   }

   return retval;
}

/* Non-template version, autodetecting the reader type */
static bool read_next_timestep(const std::string& filename_pattern, double t, int step, Field& E0, Field& E1,
   Field& B0, Field& B1, Field& V, bool doV, int& input_file_counter) {

   char filename_buffer[256];
   snprintf(filename_buffer,256,filename_pattern.c_str(),input_file_counter);

   return read_next_timestep<vlsvinterface::Reader>(filename_pattern, t,
       step,E0,E1,B0,B1,V,doV,input_file_counter);
}

/* Read E- and B-Fields as well as velocity field from a vlsv file */
template <class Reader>
void readfields(const char* filename, Field& E, Field& B, Field& V) {
   Reader r;

#ifdef DEBUG
   std::cerr << "Opening " << filename << "...";
#endif
   r.open(filename);
#ifdef DEBUG
   std::cerr <<"ok." << std::endl;
#endif

   /* Check whethere we got volume-centered fields */
   detect_field_names<Reader>(r);

   /* Read the MESH, yielding the CellIDs */
   std::vector<uint64_t> cellIds = readCellIds(r);

   /* Also read the raw field data */
   std::string name= B_field_name;
   std::vector<double> Bbuffer = readFieldData(r,name,3u);
   name = E_field_name;
   std::vector<double> Ebuffer = readFieldData(r,name,3u);

   name = "rho_v";
   std::vector<double> rho_v_buffer = readFieldData(r,name,3u);
   name = "rho";
   std::vector<double> rho_buffer = readFieldData(r,name,1u);

   /* Coordinate Boundaries */
   double min[3], max[3], time;
   uint64_t cells[3];
   r.readParameter("xmin",min[0]);
   r.readParameter("ymin",min[1]);
   r.readParameter("zmin",min[2]);
   r.readParameter("xmax",max[0]);
   r.readParameter("ymax",max[1]);
   r.readParameter("zmax",max[2]);
   r.readParameter("xcells_ini",cells[0]);
   r.readParameter("ycells_ini",cells[1]);
   r.readParameter("zcells_ini",cells[2]);
   if(!r.readParameter("t",time)) {
     r.readParameter("time",time);
   }

   //std::cerr << "Grid is " << cells[0] << " x " << cells[1] << " x " << cells[2] << " Cells, " << std::endl
   //          << " with dx = " << ((max[0]-min[0])/cells[0]) << ", dy = " << ((max[1]-min[1])/cells[1])
   //          << ", dz = " << ((max[2]-min[2])/cells[2]) << "." << std::endl;

   /* Allocate space for the actual field structures */
   E.data.resize(4*cells[0]*cells[1]*cells[2]);
   B.data.resize(4*cells[0]*cells[1]*cells[2]);
   V.data.resize(4*cells[0]*cells[1]*cells[2]);

   /* Sanity-check stored data sizes */
   if(3*cellIds.size() != Bbuffer.size()) {
      std::cerr << "3 * cellIDs.size (" << cellIds.size() << ") != Bbuffer.size (" << Bbuffer.size() << ")!"
         << std::endl;
      exit(1);
   }
   if(3*cellIds.size() != Ebuffer.size()) {
      std::cerr << "3 * cellIDs.size (" << cellIds.size() << ") != Ebuffer.size (" << Ebuffer.size() << ")!"
         << std::endl;
      exit(1);
   }
   if(3*cellIds.size() != rho_v_buffer.size()) {
      std::cerr << "3 * cellIDs.size (" << cellIds.size() << ") != rho_v_buffer.size (" << Ebuffer.size() << ")!"
         << std::endl;
      exit(1);
   }
   if(cellIds.size() != rho_buffer.size()) {
      std::cerr << "cellIDs.size (" << cellIds.size() << ") != rho_buffer.size (" << Ebuffer.size() << ")!"
         << std::endl;
      exit(1);
   }

   /* Set field sizes */
   for(int i=0; i<3;i++) {
      /* Volume-centered values -> shift by half a cell in all directions*/
      E.dx[i] = B.dx[i] = V.dx[i] = (max[i]-min[i])/cells[i];
      double shift = E.dx[i]/2;
      E.min[i] = B.min[i] = V.min[i] = min[i]+shift;
      E.max[i] = B.max[i] = V.max[i] = max[i]+shift;
      E.cells[i] = B.cells[i] = V.cells[i] = cells[i];
   }
   E.time = B.time = V.time = time;

   /* So, now we've got the cellIDs, the mesh size and the field values,
    * we can sort them into place */
   for(uint i=0; i< cellIds.size(); i++) {
      uint64_t c = cellIds[i];
      int64_t x = c % cells[0];
      int64_t y = (c /cells[0]) % cells[1];
      int64_t z = c /(cells[0]*cells[1]);

      double* Etgt = E.getCellRef(x,y,z);
      double* Btgt = B.getCellRef(x,y,z);
      double* Vtgt = V.getCellRef(x,y,z);
      Etgt[0] = Ebuffer[3*i];
      Etgt[1] = Ebuffer[3*i+1];
      Etgt[2] = Ebuffer[3*i+2];
      Btgt[0] = Bbuffer[3*i];
      Btgt[1] = Bbuffer[3*i+1];
      Btgt[2] = Bbuffer[3*i+2];
      Vtgt[0] = rho_v_buffer[3*i] / rho_buffer[i];
      Vtgt[1] = rho_v_buffer[3*i+1] / rho_buffer[i];
      Vtgt[2] = rho_v_buffer[3*i+2] / rho_buffer[i];
   }

   r.close();
}

/* Non-template version, autodetecting the reader type */
static void readfields(const char* filename, Field& E, Field& B, Field& V) {
  readfields<vlsvinterface::Reader>(filename,E,B,V);
}

/* For debugging purposes - dump a field into a png file */
void debug_output(Field& F, const char* filename);
