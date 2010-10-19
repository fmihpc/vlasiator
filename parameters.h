#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "definitions.h"

cuint MAX_SPA_CELLS = 240;
cuint MAX_VEL_BLOCKS = 30000;

struct Parameters {
   static real xmin;
   static real xmax;
   static real ymin;
   static real ymax;
   static real zmin;
   static real zmax;
   
   static real vxmin;
   static real vxmax;
   static real vymin;
   static real vymax;
   static real vzmin;
   static real vzmax;
   
   static uint xcells_ini;
   static uint ycells_ini;
   static uint zcells_ini;
   static uint vxblocks_ini;
   static uint vyblocks_ini;
   static uint vzblocks_ini;

   static uint tsteps;
   static uint saveInterval;
   static uint diagnInterval;
   
   Parameters();
};

#endif
