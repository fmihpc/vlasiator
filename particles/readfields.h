#pragma once

#include "vlsvreader2.h"
#include "vlsv_reader.h"
#include "vlsvreaderinterface.h"
#include "field.h"
#include <vector>

/* Read the cellIDs into an array */
std::vector<uint64_t> readCellIds(oldVlsv::Reader& r);
std::vector<uint64_t> readCellIds(newVlsv::Reader& r);

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

double readDoubleParameter(newVlsv::Reader& r, const char* name);
double readDoubleParameter(oldVlsv::Reader& r, const char* name);

/* Read a single-valued integer parameter */
uint32_t readUintParameter(newVlsv::Reader& r, const char* name);
uint32_t readUintParameter(oldVlsv::Reader& r, const char* name);

/* Read E- and B-Fields as well as velocity field from a vlsv file */
template <class Reader>
void readfields(std::string& filename, Field& E, Field& B, Field& V) {
	Reader r;

	std::cerr << "Opening " << filename << "...";
	r.open(filename);
	std::cerr <<"ok." << std::endl;

	/* Read the MESH, yielding the CellIDs */
	std::vector<uint64_t> cellIds = readCellIds(r);

	/* Also read the raw field data */
	std::string name("B_vol");
	std::vector<double> Bbuffer = readFieldData(r,name,3u);
	name = "E_vol";
	std::vector<double> Ebuffer = readFieldData(r,name,3u);

	name = "rho_v";
	std::vector<double> rho_v_buffer = readFieldData(r,name,3u);
	name = "rho";
	std::vector<double> rho_buffer = readFieldData(r,name,1u);

	/* Coordinate Boundaries */
	double min[3], max[3];
	uint64_t cells[3];
	min[0] = readDoubleParameter(r,"xmin");
	min[1] = readDoubleParameter(r,"ymin");
	min[2] = readDoubleParameter(r,"zmin");
	max[0] = readDoubleParameter(r,"xmax");
	max[1] = readDoubleParameter(r,"ymax");
	max[2] = readDoubleParameter(r,"zmax");
	cells[0] = readUintParameter(r,"xcells_ini");
	cells[1] = readUintParameter(r,"ycells_ini");
	cells[2] = readUintParameter(r,"zcells_ini");

	//std::cerr << "Grid is " << cells[0] << " x " << cells[1] << " x " << cells[2] << " Cells, " << std::endl
	//	<< " with dx = " << ((max[0]-min[0])/cells[0]) << ", dy = " << ((max[1]-min[1])/cells[1])
	//	<< ", dz = " << ((max[2]-min[2])/cells[2]) << "." << std::endl;

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

/* For debugging purposes - dump a field into a png file */
void debug_output(Field& F, const char* filename);
