#include <cstdlib>
#include <map>
#include <limits>
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <assert.h> 
#include <sstream>
#include "vlsreader.h"

using namespace std;
using namespace VlsHeader;
using namespace VlsVariable;

typedef int VtkInt;
typedef float VtkReal;

bool writeBinary = true;
bool writerIsLittleEndian; // If true, vls2vtk is run on LittleEndian computer.

bool almostEqual(const float& a,const float& b,const int& maxDiff) {
   assert(sizeof(float) == sizeof(int));
   if (a == b) return true;
   const int intDiff = abs(*(int*)&a - *(int*)&b);
   if (intDiff <= maxDiff) return true;
   return false;
}

template<typename T> void writeStream(ostream& os,const T& value) {
   const char* ptr = reinterpret_cast<const char*>(&value);
   if (writerIsLittleEndian == true) {
      for (int i=sizeof(value)-1; i>=0; --i) os << ptr[i];
   } else {
      for (int i=0; i<sizeof(value); ++i) os << ptr[i];
   }
}

template<typename REAL> struct NodeCrd {
   static REAL EPS;
   REAL x;
   REAL y;
   REAL z;

   NodeCrd(const REAL& x,const REAL& y,const REAL& z): x(x),y(y),z(z) {
      
   }
   /*
   // Return true if a is less than b:
   static bool operator()(const NodeCrd<REAL>& a,const NodeCrd<REAL>& b) {
      if (a.z < b.z) return true;
      if (a.y < b.y) return true;
      if (a.x < b.x) return true;
      return false;
   }
   */

   bool comp(const NodeCrd<REAL>& n) const {
      REAL EPS1,EPS2,EPS;
      EPS1 = 1.0e-6 * fabs(x);
      EPS2 = 1.0e-6 * fabs(n.x);
      if (x == 0.0) EPS1 = 1.0e-8;
      if (n.x == 0.0) EPS2 = 1.0e-8;
      EPS = max(EPS1,EPS2);
      if (fabs(x - n.x) > EPS) return false;
      //if (almostEqual(x,n.x,5) == false) return false;

      EPS1 = 1.0e-6 * fabs(y);
      EPS2 = 1.0e-6 * fabs(n.y);
      if (y == 0.0) EPS1 = 1.0e-8;
      if (n.y == 0.0) EPS2 = 1.0e-8;
      EPS = max(EPS1,EPS2);
      if (fabs(y - n.y) > EPS) return false;

      EPS1 = 1.0e-6 * fabs(z);
      EPS2 = 1.0e-6 * fabs(n.z);
      if (z == 0.0) EPS1 = 1.0e-8;
      if (n.z == 0.0) EPS2 = 1.0e-8;
      EPS = max(EPS1,EPS2);
      if (fabs(z - n.z) > EPS) return false;
      return true;
   }
};

template<> float NodeCrd<float>::EPS = +1.0e-20;
template<> double NodeCrd<double>::EPS = +1.0e-20;

// Return true if a is to be placed at an earlier position than b, false otherwise.
// Strict weak ordering: objects are equal if both f(a,b) and f(b,a) are false.
struct NodeComp {
   bool operator()(const NodeCrd<float>& a,const NodeCrd<float>& b) const {
      if (a.comp(b) == true) return false;
      float EPS = 0.5e-6 * (fabs(a.z) + fabs(b.z));
      
      if (a.z > b.z + EPS) return false;
      if (a.z < b.z - EPS) return true;
      EPS = 0.5e-6 * (fabs(a.y) + fabs(b.y));
      
      if (a.y > b.y + EPS) return false;
      if (a.y < b.y - EPS) return true;
      EPS = 0.5e-6 * (fabs(a.x) + fabs(b.x));
      
      if (a.x > b.x + EPS) return false;
      if (a.x < b.x - EPS) return true;
      cerr << "ERROR" << endl;
      return false;
   }
};

void printIns(map<NodeCrd<VtkReal>,VtkInt>::iterator& it,const NodeCrd<VtkReal>& crd) {
   if (it->first.comp(crd) == false) cerr << "ERROR" << endl;
   /*
   cout.precision(8);
   cout << showpos << scientific;
   cout << "Insert failed for (";
   cout << crd.x << ',' << crd.y << ',' << crd.z;
   cout << ") because (";
   cout << it->first.x << ',' << it->first.y << ',' << it->first.z;
   cout << ") has the same key" << endl;
    */
}

bool writeVtkHeader(fstream& out) {
   if (out.good() == false) return false;
   out << "# vtk DataFile Version 3.0" << endl;
   if (out.good() == false) return false;
   out << "vlasov sim file" << endl;
   if (out.good() == false) return false;
   if (writeBinary == false) out << "ASCII" << endl;
   else out << "BINARY" << endl;
   if (out.good() == false) return false;
   return true;
}

bool writeVtkPoints(fstream& out,const map<NodeCrd<VtkReal>,VtkInt,NodeComp>& nodes) {
   if (out.good() == false) return false;
   
   if (writeBinary == false) {
      out << "DATASET UNSTRUCTURED_GRID" << endl;
      if (out.good() == false) return false;
      out << "POINTS " << nodes.size() << " float" << endl;
      for (map<NodeCrd<VtkReal>,VtkInt,NodeComp>::const_iterator it=nodes.begin(); it!=nodes.end(); ++it) {
	 out << it->first.x << ' ' << it->first.y << ' ' << it->first.z << endl;
	 if (out.good() == false) return false;
      }
      out << endl;
   } else {
      out << "DATASET UNSTRUCTURED_GRID" << endl;
      out << "POINTS " << nodes.size() << " float" << endl;
      if (out.good() == false) return false;
      
      const char* ptr;
      for (map<NodeCrd<VtkReal>,VtkInt,NodeComp>::const_iterator it=nodes.begin(); it!=nodes.end(); ++it) {
	 //out.write(reinterpret_cast<const char*>(&(it->first.x)),sizeof(VtkReal));
	 //out.write(reinterpret_cast<const char*>(&(it->first.y)),sizeof(VtkReal));
	 //out.write(reinterpret_cast<const char*>(&(it->first.z)),sizeof(VtkReal));
	 writeStream(out,it->first.x);
	 writeStream(out,it->first.y);
	 writeStream(out,it->first.z);
	 if (out.good() == false) return false;
      }
      out << endl;
   }
   return true;
}

bool writeVtkCell(fstream& out,const map<NodeCrd<VtkReal>,VtkInt,NodeComp>& nodes,
		  const VtkReal& x,const VtkReal& y,const VtkReal& z,const VtkReal& dx,const VtkReal& dy,const VtkReal& dz) {
   map<NodeCrd<VtkReal>,VtkInt,NodeComp>::const_iterator it;
   if (it == nodes.end()) return false;
   if (writeBinary == false) {
      out << 8 << ' ';
      it = nodes.find(NodeCrd<VtkReal>(x   ,y   ,z   ));
      if (it != nodes.end()) out << it->second << ' ';
      it = nodes.find(NodeCrd<VtkReal>(x+dx,y   ,z   ));
      if (it != nodes.end()) out << it->second << ' ';
      it = nodes.find(NodeCrd<VtkReal>(x   ,y+dy,z   ));
      if (it != nodes.end()) out << it->second << ' ';
      it = nodes.find(NodeCrd<VtkReal>(x+dx,y+dy,z   ));
      if (it != nodes.end()) out << it->second << ' ';
      it = nodes.find(NodeCrd<VtkReal>(x   ,y   ,z+dz));
      if (it != nodes.end()) out << it->second << ' ';
      it = nodes.find(NodeCrd<VtkReal>(x+dx,y   ,z+dz));
      if (it != nodes.end()) out << it->second << ' ';
      it = nodes.find(NodeCrd<VtkReal>(x   ,y+dy,z+dz));
      if (it != nodes.end()) out << it->second << ' ';
      it = nodes.find(NodeCrd<VtkReal>(x+dx,y+dy,z+dz));
      if (it != nodes.end()) out << it->second << ' ';
      out << endl;
   } else {
      //VtkInt number8 = 8;
      //out.write(reinterpret_cast<char*>(&number8),sizeof(VtkInt));
      writeStream(out,(VtkInt)8);
      
      it = nodes.find(NodeCrd<VtkReal>(x   ,y   ,z   ));
      //if (it != nodes.end()) out.write(reinterpret_cast<char*>(const_cast<VtkInt*>(&(it->second))),sizeof(VtkInt));
      if (it != nodes.end()) writeStream(out,it->second);
      
      it = nodes.find(NodeCrd<VtkReal>(x+dx,y   ,z   ));
      //if (it != nodes.end()) out.write(reinterpret_cast<char*>(const_cast<VtkInt*>(&(it->second))),sizeof(VtkInt));
      if (it != nodes.end()) writeStream(out,it->second);
      
      it = nodes.find(NodeCrd<VtkReal>(x   ,y+dy,z   ));
      //if (it != nodes.end()) out.write(reinterpret_cast<char*>(const_cast<VtkInt*>(&(it->second))),sizeof(VtkInt));
      if (it != nodes.end()) writeStream(out,it->second);
      
      it = nodes.find(NodeCrd<VtkReal>(x+dx,y+dy,z   ));
      //if (it != nodes.end()) out.write(reinterpret_cast<char*>(const_cast<VtkInt*>(&(it->second))),sizeof(VtkInt));
      if (it != nodes.end()) writeStream(out,it->second);
      
      it = nodes.find(NodeCrd<VtkReal>(x   ,y   ,z+dz));
      //if (it != nodes.end()) out.write(reinterpret_cast<char*>(const_cast<VtkInt*>(&(it->second))),sizeof(VtkInt));
      if (it != nodes.end()) writeStream(out,it->second);
      
      it = nodes.find(NodeCrd<VtkReal>(x+dx,y   ,z+dz));
      //if (it != nodes.end()) out.write(reinterpret_cast<char*>(const_cast<VtkInt*>(&(it->second))),sizeof(VtkInt));
      if (it != nodes.end()) writeStream(out,it->second);
      
      it = nodes.find(NodeCrd<VtkReal>(x   ,y+dy,z+dz));
      //if (it != nodes.end()) out.write(reinterpret_cast<char*>(const_cast<VtkInt*>(&(it->second))),sizeof(VtkInt));
      if (it != nodes.end()) writeStream(out,it->second);

      it = nodes.find(NodeCrd<VtkReal>(x+dx,y+dy,z+dz));
      //if (it != nodes.end()) out.write(reinterpret_cast<char*>(const_cast<VtkInt*>(&(it->second))),sizeof(VtkInt));
      if (it != nodes.end()) writeStream(out,it->second);
   }
   return true;
}

bool writeVtkCellTypes(fstream& out,unsigned int N_cells) {
   if (writeBinary == false) {
      out << "CELL_TYPES " << N_cells << endl;
      for (unsigned int i=0; i<N_cells; ++i) out << 11 << endl;
      out << endl;
   } else {
      out << "CELL_TYPES " << N_cells << endl;
      //VtkInt number11 = 11;
      for (unsigned int i=0; i<N_cells; ++i) {
	 //out.write(reinterpret_cast<char*>(&number11),sizeof(VtkInt));
	 writeStream(out,(VtkInt)11);
      }
      out << endl;
   }
   return true;
}

bool writeVtkComponentHeader(fstream& fout,const size_t& varID,VlsReader& vlsReader) {
   if (fout.good() == false) return false;
   switch (vlsReader.getStaticVarType(varID)) {
    case NULLVARIABLE:
      return false;
      break;
    case SCALAR:
      fout << "SCALARS " << vlsReader.getStaticVarName(varID) << " float 1" << endl;
      fout << "LOOKUP_TABLE default" << endl;
      break;
    case VECTOR2:
      fout << "VECTORS " << vlsReader.getStaticVarName(varID) << " float" << endl;
      break;
    case VECTOR3:
      fout << "VECTORS " << vlsReader.getStaticVarName(varID) << " float" << endl;
      break;
    case TENSOR22:
      fout << "TENSORS " << vlsReader.getStaticVarName(varID) << " float" << endl;
      break;
    case TENSOR23:
      fout << "TENSORS " << vlsReader.getStaticVarName(varID) << " float" << endl;
      break;
    case TENSOR32:
      fout << "TENSORS " << vlsReader.getStaticVarName(varID) << " float" << endl;
      break;
    case TENSOR33:
      fout << "TENSORS " << vlsReader.getStaticVarName(varID) << " float" << endl;
      break;
    default:
      return false;
      break;
   }
   return true;
}

bool writeVtkStaticVariable(fstream& fout,const size_t& varID,VlsReader& vlsReader) {
   if (fout.good() == false) return false;
   if (writeBinary == false) {
      switch (vlsReader.getStaticVarType(varID)) {
       case NULLVARIABLE:
	 break;
       case SCALAR:
	 fout << vlsReader.getStaticVar<VtkReal>(varID,0) << endl;
	 break;
       case VECTOR2:
	 fout << vlsReader.getStaticVar<VtkReal>(varID,0) << ' ' << vlsReader.getStaticVar<VtkReal>(varID,1) << ' ' << "0.0" << endl;
	 break;
       case VECTOR3:
	 fout << vlsReader.getStaticVar<VtkReal>(varID,0) << ' ' << vlsReader.getStaticVar<VtkReal>(varID,1) << ' ';
	 fout << vlsReader.getStaticVar<VtkReal>(varID,2) << endl;
	 break;
       case TENSOR22:
	 fout << vlsReader.getStaticVar<VtkReal>(varID,0) << ' ' << vlsReader.getStaticVar<VtkReal>(varID,1) << ' ' << "0.0" << endl;
	 fout << vlsReader.getStaticVar<VtkReal>(varID,2) << ' ' << vlsReader.getStaticVar<VtkReal>(varID,3) << ' ' << "0.0" << endl;
	 fout << "0.0" << ' ' << "0.0" << ' ' << "0.0" << endl;
	 break;
       case TENSOR23:
	 fout << vlsReader.getStaticVar<VtkReal>(varID,0) << ' ' << vlsReader.getStaticVar<VtkReal>(varID,1) << ' ' << "0.0" << endl;
	 fout << vlsReader.getStaticVar<VtkReal>(varID,2) << ' ' << vlsReader.getStaticVar<VtkReal>(varID,3) << ' ' << "0.0" << endl;
	 fout << vlsReader.getStaticVar<VtkReal>(varID,4) << ' ' << vlsReader.getStaticVar<VtkReal>(varID,5) << ' ' << "0.0" << endl;
	 break;
       case TENSOR32:
	 fout << vlsReader.getStaticVar<VtkReal>(varID,0) << ' ' << vlsReader.getStaticVar<VtkReal>(varID,1) << ' ';
	 fout << vlsReader.getStaticVar<VtkReal>(varID,2) << endl;
	 fout << vlsReader.getStaticVar<VtkReal>(varID,3) << ' ' << vlsReader.getStaticVar<VtkReal>(varID,4) << ' ';
	 fout << vlsReader.getStaticVar<VtkReal>(varID,5) << endl;
	 fout << "0.0" << ' ' << "0.0" << ' ' << "0.0" << endl;
	 break;
       case TENSOR33:
	 fout << vlsReader.getStaticVar<VtkReal>(varID,0) << ' ' << vlsReader.getStaticVar<VtkReal>(varID,1) << ' ';
	 fout << vlsReader.getStaticVar<VtkReal>(varID,2) << endl;
	 fout << vlsReader.getStaticVar<VtkReal>(varID,3) << ' ' << vlsReader.getStaticVar<VtkReal>(varID,4) << ' ';
	 fout << vlsReader.getStaticVar<VtkReal>(varID,5) << endl;
	 fout << vlsReader.getStaticVar<VtkReal>(varID,6) << ' ' << vlsReader.getStaticVar<VtkReal>(varID,7) << ' ';
	 fout << vlsReader.getStaticVar<VtkReal>(varID,8) << endl;
	 break;
       default:
	 break;
      }
   } else {
      switch (vlsReader.getStaticVarType(varID)) {
       case NULLVARIABLE:
	 break;
       case SCALAR:
	 writeStream(fout,vlsReader.getStaticVar<VtkReal>(varID,0));
	 break;
       case VECTOR2:
	 writeStream(fout,vlsReader.getStaticVar<VtkReal>(varID,0));
	 writeStream(fout,vlsReader.getStaticVar<VtkReal>(varID,1));
	 writeStream(fout,(VtkReal)0.0);
	 break;
       case VECTOR3:
	 writeStream(fout,vlsReader.getStaticVar<VtkReal>(varID,0));
	 writeStream(fout,vlsReader.getStaticVar<VtkReal>(varID,1));
	 writeStream(fout,vlsReader.getStaticVar<VtkReal>(varID,2));
	 break;
       case TENSOR22:
	 writeStream(fout,vlsReader.getStaticVar<VtkReal>(varID,0));
	 writeStream(fout,vlsReader.getStaticVar<VtkReal>(varID,1));
	 writeStream(fout,(VtkReal)0.0);
	 writeStream(fout,vlsReader.getStaticVar<VtkReal>(varID,2));
	 writeStream(fout,vlsReader.getStaticVar<VtkReal>(varID,3));
	 writeStream(fout,(VtkReal)0.0);
	 writeStream(fout,(VtkReal)0.0);
	 writeStream(fout,(VtkReal)0.0);
	 writeStream(fout,(VtkReal)0.0);
	 break;
       case TENSOR23:
	 writeStream(fout,vlsReader.getStaticVar<VtkReal>(varID,0));
	 writeStream(fout,vlsReader.getStaticVar<VtkReal>(varID,1));
	 writeStream(fout,(VtkReal)0.0);
	 writeStream(fout,vlsReader.getStaticVar<VtkReal>(varID,2));
	 writeStream(fout,vlsReader.getStaticVar<VtkReal>(varID,3));
	 writeStream(fout,(VtkReal)0.0);
	 writeStream(fout,vlsReader.getStaticVar<VtkReal>(varID,4));
	 writeStream(fout,vlsReader.getStaticVar<VtkReal>(varID,5));
	 writeStream(fout,(VtkReal)0.0);
	 break;
       case TENSOR32:
	 writeStream(fout,vlsReader.getStaticVar<VtkReal>(varID,0));
	 writeStream(fout,vlsReader.getStaticVar<VtkReal>(varID,1));
	 writeStream(fout,vlsReader.getStaticVar<VtkReal>(varID,2));
	 writeStream(fout,vlsReader.getStaticVar<VtkReal>(varID,3));
	 writeStream(fout,vlsReader.getStaticVar<VtkReal>(varID,4));
	 writeStream(fout,vlsReader.getStaticVar<VtkReal>(varID,5));
	 writeStream(fout,(VtkReal)0.0);
	 writeStream(fout,(VtkReal)0.0);
	 writeStream(fout,(VtkReal)0.0);
	 break;
       case TENSOR33:
	 writeStream(fout,vlsReader.getStaticVar<VtkReal>(varID,0));
	 writeStream(fout,vlsReader.getStaticVar<VtkReal>(varID,1));
	 writeStream(fout,vlsReader.getStaticVar<VtkReal>(varID,2));
	 writeStream(fout,vlsReader.getStaticVar<VtkReal>(varID,3));
	 writeStream(fout,vlsReader.getStaticVar<VtkReal>(varID,4));
	 writeStream(fout,vlsReader.getStaticVar<VtkReal>(varID,5));
	 writeStream(fout,vlsReader.getStaticVar<VtkReal>(varID,6));
	 writeStream(fout,vlsReader.getStaticVar<VtkReal>(varID,7));
	 writeStream(fout,vlsReader.getStaticVar<VtkReal>(varID,8));
	 break;
      }
   }
   return fout.good();
}

int main(int argn,char* args[]) {
   if (argn < 2) {
      cerr << endl << "USAGE: ./vls2vtk <input file 1> <input file 2> ..." << endl << endl;
      return 1;
   }
   
   // Test if program is run on little endian computer:
     {
	const int number = 1;
	const char* const ptr = reinterpret_cast<const char*>(&number);
	if (ptr[0] == 1) writerIsLittleEndian = true;
	else writerIsLittleEndian = false;
     }
   
   // Init MPI
   int N_processes,myrank;
   if (MPI_Init(&argn,&args) != MPI_SUCCESS) {cerr << "ERROR: MPI failed to init!" << endl; return 1;}
   MPI_Comm_size(MPI_COMM_WORLD,&N_processes);
   MPI_Comm_rank(MPI_COMM_WORLD,&myrank);

   bool success = true;

   for (int arg = 1; arg < argn; arg++) {
	   VlsReader vlsReader;
	   VtkReal x,y,z,dx,dy,dz;
	   map<NodeCrd<VtkReal>,VtkInt,NodeComp> nodes;
	   fstream out;

	   unsigned int N_cells = 0;
	   
	   // Open file for reading:
	   const string infile = args[arg];
	   if (infile.find(".vlsv") == string::npos) {
	      cerr << "ERROR: Input file is not a .vlsv file!" << endl;
	      success == false;      
	   }
	   
	   if (success == true) if (vlsReader.open(infile) == false) {
	      cerr << "ERROR: Could not open file '" << infile << "' for reading!" << endl; 
	      success = false;
	   }
	   
	   // Read file header:
	   if (success == true) if (vlsReader.readHeader() == false) {
	      cerr << "ERROR when reading header!" << endl;
	      vlsReader.close();
	      success = false;
	   }

	   // Read coordinate entries and push all unique node coordinates into map nodes:
	   if (myrank == 0 && success == true) {
	      while (vlsReader.readSpatCellCoordEntry() == true) {
		 // Get coordinate entry values:
		 x = vlsReader.getCrdX<VtkReal>();
		 y = vlsReader.getCrdY<VtkReal>();
		 z = vlsReader.getCrdZ<VtkReal>();
		 dx = vlsReader.getDx<VtkReal>();
		 dy = vlsReader.getDy<VtkReal>();
		 dz = vlsReader.getDz<VtkReal>();

		 // Push nodes into map:
		 pair<map<NodeCrd<VtkReal>,VtkInt>::iterator,bool> pos;
		 pos = nodes.insert( pair<NodeCrd<VtkReal>,VtkInt>(NodeCrd<VtkReal>(x   ,y   ,z   ),numeric_limits<VtkInt>::max()) );
		 pos = nodes.insert( pair<NodeCrd<VtkReal>,VtkInt>(NodeCrd<VtkReal>(x+dx,y   ,z   ),numeric_limits<VtkInt>::max()) );
		 pos = nodes.insert( pair<NodeCrd<VtkReal>,VtkInt>(NodeCrd<VtkReal>(x+dx,y+dy,z   ),numeric_limits<VtkInt>::max()) );
		 pos = nodes.insert( pair<NodeCrd<VtkReal>,VtkInt>(NodeCrd<VtkReal>(x   ,y+dy,z   ),numeric_limits<VtkInt>::max()) );
		 pos = nodes.insert( pair<NodeCrd<VtkReal>,VtkInt>(NodeCrd<VtkReal>(x   ,y   ,z+dz),numeric_limits<VtkInt>::max()) );
		 pos = nodes.insert( pair<NodeCrd<VtkReal>,VtkInt>(NodeCrd<VtkReal>(x+dx,y   ,z+dz),numeric_limits<VtkInt>::max()) );
		 pos = nodes.insert( pair<NodeCrd<VtkReal>,VtkInt>(NodeCrd<VtkReal>(x+dx,y+dy,z+dz),numeric_limits<VtkInt>::max()) );
		 pos = nodes.insert( pair<NodeCrd<VtkReal>,VtkInt>(NodeCrd<VtkReal>(x   ,y+dy,z+dz),numeric_limits<VtkInt>::max()) );
		 ++N_cells;
	      }
	   }
	   vlsReader.close();
	   //cerr << "nodes.size() = " << nodes.size() << endl;
	   
	   // Give each node a unique ID:
	   if (myrank == 0 && success == true) {
	      unsigned int counter = 0;
	      for (map<NodeCrd<VtkReal>,VtkInt,NodeComp>::iterator it=nodes.begin(); it!=nodes.end(); ++it) {
		 it->second = counter;
		 ++counter;
	      }
	   }
	   
	   // Open Vtk file for writing:
	   string fileout;
	   if (myrank == 0 && success == true) {
	      fileout = infile;
	      size_t pos = fileout.rfind(".vlsv");
	      if (pos != string::npos) fileout.replace(pos,5,".vtk");
	      //cerr << "out file is '" << fileout << "'" << endl;
	      
	      if (writeBinary == false) out.open(fileout.c_str(),fstream::out);
	      else out.open(fileout.c_str(),fstream::out);
	      if (out.good() == false) {
		 cerr << "ERROR: Could not open '" << fileout << "' for writing!" << endl;
		 success = false;
	      }
	   }
	   
	   // Write Vtk header:
	   if (myrank == 0 && success == true) {
	      if (writeVtkHeader(out) == false) {
		 cerr << "ERROR: Failed to write Vtk header!" << endl;
		 success = false;
	      }
	   }

	   // Write node coordinates into Vtk file:
	   if (myrank == 0 && success == true) {
	      if (writeVtkPoints(out,nodes) == false) {
		 cerr << "ERROR: Failed to write Vtk file!" << endl;
		 success = false;
	      }
	   }

	   // Open one output file per variable:
	   size_t N_staticVars = vlsReader.getNumberOfStaticVars();
	   vector<fstream*> varFiles(N_staticVars);
	   for (size_t i=0; i<N_staticVars; ++i) {
	      stringstream ss;
	      ss << "var";
	      ss.width(5);
	      ss.fill('0');
	      ss << i << ".vtk";
	      string fname;
	      ss >> fname;

	      varFiles[i] = new fstream;
	      varFiles[i]->open(fname.c_str(), fstream::out);
	      if (varFiles[i]->good() == false) {
		 cerr << "ERROR: Failed to create outfile for variable '" << vlsReader.getStaticVarName(i) << "'" << endl;
	      }
	      
	      //*(varFiles[i]) << "SCALARS " << vlsReader.getStaticVarName(i) << " float 1" << endl;
	      //*(varFiles[i]) << "LOOKUP_TABLE default" << endl;
	      writeVtkComponentHeader(*(varFiles[i]),i,vlsReader);
	   }

	   // Write cells into Vtk file:
	   if (success == true) if (vlsReader.open(infile) == false) success = false;
	   if (success == true) if (vlsReader.readHeader() == false) success = false;
	   if (myrank == 0 && success == true) {
	      out << "CELLS " << N_cells << ' ' << N_cells*9 << endl;
	      
	      while (vlsReader.readSpatCellCoordEntry() == true) {
		 x = vlsReader.getCrdX<VtkReal>();
		 y = vlsReader.getCrdY<VtkReal>();
		 z = vlsReader.getCrdZ<VtkReal>();
		 dx = vlsReader.getDx<VtkReal>();
		 dy = vlsReader.getDy<VtkReal>();
		 dz = vlsReader.getDz<VtkReal>();

		 writeVtkCell(out,nodes,x,y,z,dx,dy,dz);
		 
		 for (size_t i=0; i<N_staticVars; ++i) {
		    if (writeVtkStaticVariable(*(varFiles[i]),i,vlsReader) == false) {
		       cerr << "vls2vtk: An ERROR has occurred while writing component '" << vlsReader.getStaticVarName(i) << "'" << endl;
		       break;
		    }
		 }
	      }
	      
	      out << endl;
	      for (size_t i=0; i<N_staticVars; ++i) *(varFiles[i]) << endl;
	   }
	   vlsReader.close();

	   // Write cell types into Vtk file:
	   if (success == true && myrank == 0) {
	      if (writeVtkCellTypes(out,N_cells) == false) {
		 cerr << "ERROR: Failed to write cell types!" << endl;
		 success = false;
	      }
	   }
	   
	   // Append each component file to the output vtk file:
	   if (myrank == 0) {
	      out << "CELL_DATA " << N_cells << endl;
	      out.close();
	      
	      for (size_t i=0; i<N_staticVars; ++i) {
		 stringstream ss;
		 ss << "var";
		 ss.width(5);
		 ss.fill('0');
		 ss << i << ".vtk";
		 
		 string command = "cat " + ss.str();
		 command = command + " >> " + fileout;
		 //cout << "Command is '" << command << "'" << endl;
		 if (system(command.c_str()) != 0) {
		    cerr << "An ERROR occurred while appending variable file '" << ss.str() << "'" << endl;
		 }
	      }
	   }
	   
	   // Close files:
	   out.close();
	   for (size_t i=0; i<varFiles.size(); ++i) {
	      varFiles[i]->close();
	      delete varFiles[i];
	      varFiles[i] = NULL;
	   }
	   
	   // Delete temporary variable files:
	   if (myrank == 0) {
	      for (size_t i=0; i<N_staticVars; ++i) {
		 stringstream ss;
		 ss << "var";
		 ss.width(5);
		 ss.fill('0');
		 ss << i << ".vtk";
		 string command = "rm " + ss.str();
		 //cout << "Command is '" << command << "'" << endl;
		 if (system(command.c_str()) != 0) {
		    cerr << "An ERROR occurred while removing temporary variable file '" << ss.str() << "'" << endl;
		 }
	      }
	   }
   }
   
   // Exit program:
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();   
   return 0;
}




