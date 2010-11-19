#include <cmath>
#include <limits>
#include "parameters.h"

using namespace std;

typedef Parameters P;

// Define static members:
Real P::xmin = NAN;
Real P::xmax = NAN;
Real P::ymin = NAN;
Real P::ymax = NAN;
Real P::zmin = NAN;
Real P::zmax = NAN;
Real P::dx_ini = NAN;
Real P::dy_ini = NAN;
Real P::dz_ini = NAN;

Real P::vxmin = NAN;
Real P::vxmax = NAN;
Real P::vymin = NAN;
Real P::vymax = NAN;
Real P::vzmin = NAN;
Real P::vzmax = NAN;

uint P::xcells_ini = numeric_limits<uint>::max();
uint P::ycells_ini = numeric_limits<uint>::max();
uint P::zcells_ini = numeric_limits<uint>::max();
uint P::vxblocks_ini = numeric_limits<uint>::max();
uint P::vyblocks_ini = numeric_limits<uint>::max();
uint P::vzblocks_ini = numeric_limits<uint>::max();

Real P::dt = NAN;
uint P::tstep = 0;
uint P::tsteps = 0;
uint P::saveInterval = numeric_limits<uint>::max();
uint P::diagnInterval = numeric_limits<uint>::max();

uint P::transmit = 0;

Parameters::Parameters() {
   xmin = -1.2;
   xmax = +1.2;
   ymin = -1.2;
   ymax = +1.2;
   zmin = +0.0;
   zmax = +0.6;

   vxmin = -1.0;
   vxmax = +1.0;
   vymin = -1.0;
   vymax = +1.0;
   vzmin = -1.0;
   vzmax = +1.0;
   
   xcells_ini = 20;
   ycells_ini = 10;
   zcells_ini = 10;
   vxblocks_ini = 10;
   vyblocks_ini = 10;
   vzblocks_ini = 10;
   
   dx_ini = (xmax-xmin)/xcells_ini;
   dy_ini = (ymax-ymin)/ycells_ini;
   dz_ini = (zmax-zmin)/zcells_ini;

   dt = 0.025;
   tsteps = 50;
   diagnInterval = 5;
}







