#include <cmath>
#include <limits>
#include "parameters.h"

using namespace std;

typedef Parameters P;

real P::xmin = NAN;
real P::xmax = NAN;
real P::ymin = NAN;
real P::ymax = NAN;
real P::zmin = NAN;
real P::zmax = NAN;

real P::vxmin = NAN;
real P::vxmax = NAN;
real P::vymin = NAN;
real P::vymax = NAN;
real P::vzmin = NAN;
real P::vzmax = NAN;

uint P::xcells_ini = numeric_limits<uint>::max();
uint P::ycells_ini = numeric_limits<uint>::max();
uint P::zcells_ini = numeric_limits<uint>::max();
uint P::vxblocks_ini = numeric_limits<uint>::max();
uint P::vyblocks_ini = numeric_limits<uint>::max();
uint P::vzblocks_ini = numeric_limits<uint>::max();

uint P::tsteps = 0;
uint P::saveInterval = numeric_limits<uint>::max();
uint P::diagnInterval = numeric_limits<uint>::max();

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
   
   xcells_ini = 34;
   ycells_ini = 7;
   zcells_ini = 1;
   vxblocks_ini = 5;
   vyblocks_ini = 5;
   vzblocks_ini = 5;
   
   tsteps = 21;
   diagnInterval = 10;
}







