/*
This file is part of Vlasiator.

Copyright 2015 Finnish Meteorological Institute

*/

#include "fs_common.h"
#include "fs_cache.h"
#include "ldz_hall.hpp"

#ifndef NDEBUG
   #define DEBUG_FSOLVER
#endif

using namespace std;

// X
template<typename REAL> inline
REAL JXBX_000_100(
                  const REAL* const pC,
                  creal BGBY,
                  creal BGBZ,
                  creal dx,
                  creal dy,
                  creal dz
                 ) {
   using namespace Rec;
   return -(pC[a_zz]*BGBZ)/dz+(pC[a_z]*BGBZ)/dz-(pC[a_yz]*BGBZ)/(2*dz)-(pC[c_xzz]*BGBZ)/(6*dx)+(pC[c_xz]*BGBZ)/(2*dx)-(pC[c_xyz]*BGBZ)/(4*dx)+(pC[c_xy]*BGBZ)/(2*dx)-(pC[c_x]*BGBZ)/dx-(pC[a_yz]*BGBY)/(2*dy)-(pC[a_yy]*BGBY)/dy+(pC[a_y]*BGBY)/dy+(pC[b_xz]*BGBY)/(2*dx)-(pC[b_xyz]*BGBY)/(4*dx)-(pC[b_xyy]*BGBY)/(6*dx)+(pC[b_xy]*BGBY)/(2*dx)-(pC[b_x]*BGBY)/dx-(pC[a_zz]*pC[c_zz])/(6*dz)+(pC[a_z]*pC[c_zz])/(6*dz)-
     (pC[a_yz]*pC[c_zz])/(12*dz)+(pC[a_zz]*pC[c_z])/(2*dz)-(pC[a_z]*pC[c_z])/(2*dz)+(pC[a_yz]*pC[c_z])/(4*dz)-(pC[a_zz]*pC[c_yz])/(4*dz)+(pC[a_z]*pC[c_yz])/(4*dz)-(pC[a_yz]*pC[c_yz])/(8*dz)+(pC[a_zz]*pC[c_y])/(2*dz)-(pC[a_z]*pC[c_y])/(2*dz)+(pC[a_yz]*pC[c_y])/(4*dz)+(pC[a_xzz]*pC[c_xz])/(24*dz)-(pC[a_xz]*pC[c_xz])/(24*dz)+(pC[a_xyz]*pC[c_xz])/(48*dz)-(pC[a_xzz]*pC[c_x])/(12*dz)+(pC[a_xz]*pC[c_x])/(12*dz)-(pC[a_xyz]*pC[c_x])/(24*dz)-(pC[a_zz]*pC[c_0])/dz+(pC[a_z]*pC[c_0])/dz-(pC[a_yz]*pC[c_0])/(2*dz)+(pC[a_yz]*pC[b_z])/(4*dy)+
     (pC[a_yy]*pC[b_z])/(2*dy)-(pC[a_y]*pC[b_z])/(2*dy)-(pC[a_yz]*pC[b_yz])/(8*dy)-(pC[a_yy]*pC[b_yz])/(4*dy)+(pC[a_y]*pC[b_yz])/(4*dy)-(pC[a_yz]*pC[b_yy])/(12*dy)-(pC[a_yy]*pC[b_yy])/(6*dy)+(pC[a_y]*pC[b_yy])/(6*dy)+(pC[a_yz]*pC[b_y])/(4*dy)+(pC[a_yy]*pC[b_y])/(2*dy)-(pC[a_y]*pC[b_y])/(2*dy)+(pC[a_xyz]*pC[b_xy])/(48*dy)+(pC[a_xyy]*pC[b_xy])/(24*dy)-(pC[a_xy]*pC[b_xy])/(24*dy)-(pC[a_xyz]*pC[b_x])/(24*dy)-(pC[a_xyy]*pC[b_x])/(12*dy)+(pC[a_xy]*pC[b_x])/(12*dy)-(pC[a_yz]*pC[b_0])/(2*dy)-(pC[a_yy]*pC[b_0])/dy+(pC[a_y]*pC[b_0])/dy-
     (pC[c_xzz]*pC[c_zz])/(36*dx)+(pC[c_xz]*pC[c_zz])/(12*dx)-(pC[c_xyz]*pC[c_zz])/(24*dx)+(pC[c_xy]*pC[c_zz])/(12*dx)-(pC[c_x]*pC[c_zz])/(6*dx)+(pC[c_xzz]*pC[c_z])/(12*dx)-(pC[c_xz]*pC[c_z])/(4*dx)+(pC[c_xyz]*pC[c_z])/(8*dx)-(pC[c_xy]*pC[c_z])/(4*dx)+(pC[c_x]*pC[c_z])/(2*dx)-(pC[c_xzz]*pC[c_yz])/(24*dx)+(pC[c_xz]*pC[c_yz])/(8*dx)-(pC[c_xyz]*pC[c_yz])/(16*dx)+(pC[c_xy]*pC[c_yz])/(8*dx)-(pC[c_x]*pC[c_yz])/(4*dx)+(pC[c_xzz]*pC[c_y])/(12*dx)-(pC[c_xz]*pC[c_y])/(4*dx)+(pC[c_xyz]*pC[c_y])/(8*dx)-(pC[c_xy]*pC[c_y])/(4*dx)+(pC[c_x]*pC[c_y])/(2*dx)-
     (pC[c_0]*pC[c_xzz])/(6*dx)-(pC[c_xxz]*pC[c_xz])/(24*dx)+(pC[c_xx]*pC[c_xz])/(12*dx)+(pC[c_0]*pC[c_xz])/(2*dx)-(pC[c_0]*pC[c_xyz])/(4*dx)+(pC[c_0]*pC[c_xy])/(2*dx)+(pC[c_x]*pC[c_xxz])/(12*dx)-(pC[c_x]*pC[c_xx])/(6*dx)-(pC[c_0]*pC[c_x])/dx-(pC[b_xz]*pC[b_z])/(4*dx)+(pC[b_xyz]*pC[b_z])/(8*dx)+(pC[b_xyy]*pC[b_z])/(12*dx)-(pC[b_xy]*pC[b_z])/(4*dx)+(pC[b_x]*pC[b_z])/(2*dx)+(pC[b_xz]*pC[b_yz])/(8*dx)-(pC[b_xyz]*pC[b_yz])/(16*dx)-(pC[b_xyy]*pC[b_yz])/(24*dx)+(pC[b_xy]*pC[b_yz])/(8*dx)-(pC[b_x]*pC[b_yz])/(4*dx)+(pC[b_xz]*pC[b_yy])/(12*dx)-
     (pC[b_xyz]*pC[b_yy])/(24*dx)-(pC[b_xyy]*pC[b_yy])/(36*dx)+(pC[b_xy]*pC[b_yy])/(12*dx)-(pC[b_x]*pC[b_yy])/(6*dx)-(pC[b_xz]*pC[b_y])/(4*dx)+(pC[b_xyz]*pC[b_y])/(8*dx)+(pC[b_xyy]*pC[b_y])/(12*dx)-(pC[b_xy]*pC[b_y])/(4*dx)+(pC[b_x]*pC[b_y])/(2*dx)+(pC[b_0]*pC[b_xz])/(2*dx)-(pC[b_0]*pC[b_xyz])/(4*dx)-(pC[b_0]*pC[b_xyy])/(6*dx)-(pC[b_xxy]*pC[b_xy])/(24*dx)+(pC[b_xx]*pC[b_xy])/(12*dx)+(pC[b_0]*pC[b_xy])/(2*dx)+(pC[b_x]*pC[b_xxy])/(12*dx)-(pC[b_x]*pC[b_xx])/(6*dx)-(pC[b_0]*pC[b_x])/dx;
}

template<typename REAL> inline
REAL JXBX_010_110(
                  const REAL* const pC,
                  creal BGBY,
                  creal BGBZ,
                  creal dx,
                  creal dy,
                  creal dz
                 ) {
   using namespace Rec;
   return -(pC[a_zz]*BGBZ)/dz+(pC[a_z]*BGBZ)/dz+(pC[a_yz]*BGBZ)/(2*dz)-(pC[c_xzz]*BGBZ)/(6*dx)+(pC[c_xz]*BGBZ)/(2*dx)+(pC[c_xyz]*BGBZ)/(4*dx)-(pC[c_xy]*BGBZ)/(2*dx)-(pC[c_x]*BGBZ)/dx-(pC[a_yz]*BGBY)/(2*dy)+(pC[a_yy]*BGBY)/dy+(pC[a_y]*BGBY)/dy+(pC[b_xz]*BGBY)/(2*dx)+(pC[b_xyz]*BGBY)/(4*dx)-(pC[b_xyy]*BGBY)/(6*dx)-(pC[b_xy]*BGBY)/(2*dx)-(pC[b_x]*BGBY)/dx-(pC[a_zz]*pC[c_zz])/(6*dz)+(pC[a_z]*pC[c_zz])/(6*dz)+
   (pC[a_yz]*pC[c_zz])/(12*dz)+(pC[a_zz]*pC[c_z])/(2*dz)-(pC[a_z]*pC[c_z])/(2*dz)-(pC[a_yz]*pC[c_z])/(4*dz)+(pC[a_zz]*pC[c_yz])/(4*dz)-(pC[a_z]*pC[c_yz])/(4*dz)-(pC[a_yz]*pC[c_yz])/(8*dz)-(pC[a_zz]*pC[c_y])/(2*dz)+(pC[a_z]*pC[c_y])/(2*dz)+(pC[a_yz]*pC[c_y])/(4*dz)+(pC[a_xzz]*pC[c_xz])/(24*dz)-(pC[a_xz]*pC[c_xz])/(24*dz)-(pC[a_xyz]*pC[c_xz])/(48*dz)-(pC[a_xzz]*pC[c_x])/(12*dz)+(pC[a_xz]*pC[c_x])/(12*dz)+(pC[a_xyz]*pC[c_x])/(24*dz)-(pC[a_zz]*pC[c_0])/dz+(pC[a_z]*pC[c_0])/dz+(pC[a_yz]*pC[c_0])/(2*dz)+(pC[a_yz]*pC[b_z])/(4*dy)-
   (pC[a_yy]*pC[b_z])/(2*dy)-(pC[a_y]*pC[b_z])/(2*dy)+(pC[a_yz]*pC[b_yz])/(8*dy)-(pC[a_yy]*pC[b_yz])/(4*dy)-(pC[a_y]*pC[b_yz])/(4*dy)-(pC[a_yz]*pC[b_yy])/(12*dy)+(pC[a_yy]*pC[b_yy])/(6*dy)+(pC[a_y]*pC[b_yy])/(6*dy)-(pC[a_yz]*pC[b_y])/(4*dy)+(pC[a_yy]*pC[b_y])/(2*dy)+(pC[a_y]*pC[b_y])/(2*dy)-(pC[a_xyz]*pC[b_xy])/(48*dy)+(pC[a_xyy]*pC[b_xy])/(24*dy)+(pC[a_xy]*pC[b_xy])/(24*dy)-(pC[a_xyz]*pC[b_x])/(24*dy)+(pC[a_xyy]*pC[b_x])/(12*dy)+(pC[a_xy]*pC[b_x])/(12*dy)-(pC[a_yz]*pC[b_0])/(2*dy)+(pC[a_yy]*pC[b_0])/dy+(pC[a_y]*pC[b_0])/dy-
   (pC[c_xzz]*pC[c_zz])/(36*dx)+(pC[c_xz]*pC[c_zz])/(12*dx)+(pC[c_xyz]*pC[c_zz])/(24*dx)-(pC[c_xy]*pC[c_zz])/(12*dx)-(pC[c_x]*pC[c_zz])/(6*dx)+(pC[c_xzz]*pC[c_z])/(12*dx)-(pC[c_xz]*pC[c_z])/(4*dx)-(pC[c_xyz]*pC[c_z])/(8*dx)+(pC[c_xy]*pC[c_z])/(4*dx)+(pC[c_x]*pC[c_z])/(2*dx)+(pC[c_xzz]*pC[c_yz])/(24*dx)-(pC[c_xz]*pC[c_yz])/(8*dx)-(pC[c_xyz]*pC[c_yz])/(16*dx)+(pC[c_xy]*pC[c_yz])/(8*dx)+(pC[c_x]*pC[c_yz])/(4*dx)-(pC[c_xzz]*pC[c_y])/(12*dx)+(pC[c_xz]*pC[c_y])/(4*dx)+(pC[c_xyz]*pC[c_y])/(8*dx)-(pC[c_xy]*pC[c_y])/(4*dx)-(pC[c_x]*pC[c_y])/(2*dx)-
   (pC[c_0]*pC[c_xzz])/(6*dx)-(pC[c_xxz]*pC[c_xz])/(24*dx)+(pC[c_xx]*pC[c_xz])/(12*dx)+(pC[c_0]*pC[c_xz])/(2*dx)+(pC[c_0]*pC[c_xyz])/(4*dx)-(pC[c_0]*pC[c_xy])/(2*dx)+(pC[c_x]*pC[c_xxz])/(12*dx)-(pC[c_x]*pC[c_xx])/(6*dx)-(pC[c_0]*pC[c_x])/dx-(pC[b_xz]*pC[b_z])/(4*dx)-(pC[b_xyz]*pC[b_z])/(8*dx)+(pC[b_xyy]*pC[b_z])/(12*dx)+(pC[b_xy]*pC[b_z])/(4*dx)+(pC[b_x]*pC[b_z])/(2*dx)-(pC[b_xz]*pC[b_yz])/(8*dx)-(pC[b_xyz]*pC[b_yz])/(16*dx)+(pC[b_xyy]*pC[b_yz])/(24*dx)+(pC[b_xy]*pC[b_yz])/(8*dx)+(pC[b_x]*pC[b_yz])/(4*dx)+(pC[b_xz]*pC[b_yy])/(12*dx)+
   (pC[b_xyz]*pC[b_yy])/(24*dx)-(pC[b_xyy]*pC[b_yy])/(36*dx)-(pC[b_xy]*pC[b_yy])/(12*dx)-(pC[b_x]*pC[b_yy])/(6*dx)+(pC[b_xz]*pC[b_y])/(4*dx)+(pC[b_xyz]*pC[b_y])/(8*dx)-(pC[b_xyy]*pC[b_y])/(12*dx)-(pC[b_xy]*pC[b_y])/(4*dx)-(pC[b_x]*pC[b_y])/(2*dx)+(pC[b_0]*pC[b_xz])/(2*dx)+(pC[b_0]*pC[b_xyz])/(4*dx)-(pC[b_0]*pC[b_xyy])/(6*dx)-(pC[b_xxy]*pC[b_xy])/(24*dx)-(pC[b_xx]*pC[b_xy])/(12*dx)-(pC[b_0]*pC[b_xy])/(2*dx)-(pC[b_x]*pC[b_xxy])/(12*dx)-(pC[b_x]*pC[b_xx])/(6*dx)-(pC[b_0]*pC[b_x])/dx;
}

template<typename REAL> inline
REAL JXBX_001_101(
                  const REAL* const pC,
                  creal BGBY,
                  creal BGBZ,
                  creal dx,
                  creal dy,
                  creal dz
) {
   using namespace Rec;
   return(pC[a_zz]*BGBZ)/dz+(pC[a_z]*BGBZ)/dz-(pC[a_yz]*BGBZ)/(2*dz)-(pC[c_xzz]*BGBZ)/(6*dx)-(pC[c_xz]*BGBZ)/(2*dx)+(pC[c_xyz]*BGBZ)/(4*dx)+(pC[c_xy]*BGBZ)/(2*dx)-(pC[c_x]*BGBZ)/dx+(pC[a_yz]*BGBY)/(2*dy)-(pC[a_yy]*BGBY)/dy+(pC[a_y]*BGBY)/dy-(pC[b_xz]*BGBY)/(2*dx)+(pC[b_xyz]*BGBY)/(4*dx)-(pC[b_xyy]*BGBY)/(6*dx)+(pC[b_xy]*BGBY)/(2*dx)-(pC[b_x]*BGBY)/dx+(pC[a_zz]*pC[c_zz])/(6*dz)+(pC[a_z]*pC[c_zz])/(6*dz)-
   (pC[a_yz]*pC[c_zz])/(12*dz)+(pC[a_zz]*pC[c_z])/(2*dz)+(pC[a_z]*pC[c_z])/(2*dz)-(pC[a_yz]*pC[c_z])/(4*dz)-(pC[a_zz]*pC[c_yz])/(4*dz)-(pC[a_z]*pC[c_yz])/(4*dz)+(pC[a_yz]*pC[c_yz])/(8*dz)-(pC[a_zz]*pC[c_y])/(2*dz)-(pC[a_z]*pC[c_y])/(2*dz)+(pC[a_yz]*pC[c_y])/(4*dz)+(pC[a_xzz]*pC[c_xz])/(24*dz)+(pC[a_xz]*pC[c_xz])/(24*dz)-(pC[a_xyz]*pC[c_xz])/(48*dz)+(pC[a_xzz]*pC[c_x])/(12*dz)+(pC[a_xz]*pC[c_x])/(12*dz)-(pC[a_xyz]*pC[c_x])/(24*dz)+(pC[a_zz]*pC[c_0])/dz+(pC[a_z]*pC[c_0])/dz-(pC[a_yz]*pC[c_0])/(2*dz)+(pC[a_yz]*pC[b_z])/(4*dy)-
   (pC[a_yy]*pC[b_z])/(2*dy)+(pC[a_y]*pC[b_z])/(2*dy)-(pC[a_yz]*pC[b_yz])/(8*dy)+(pC[a_yy]*pC[b_yz])/(4*dy)-(pC[a_y]*pC[b_yz])/(4*dy)+(pC[a_yz]*pC[b_yy])/(12*dy)-(pC[a_yy]*pC[b_yy])/(6*dy)+(pC[a_y]*pC[b_yy])/(6*dy)-(pC[a_yz]*pC[b_y])/(4*dy)+(pC[a_yy]*pC[b_y])/(2*dy)-(pC[a_y]*pC[b_y])/(2*dy)-(pC[a_xyz]*pC[b_xy])/(48*dy)+(pC[a_xyy]*pC[b_xy])/(24*dy)-(pC[a_xy]*pC[b_xy])/(24*dy)+(pC[a_xyz]*pC[b_x])/(24*dy)-(pC[a_xyy]*pC[b_x])/(12*dy)+(pC[a_xy]*pC[b_x])/(12*dy)+(pC[a_yz]*pC[b_0])/(2*dy)-(pC[a_yy]*pC[b_0])/dy+(pC[a_y]*pC[b_0])/dy-
   (pC[c_xzz]*pC[c_zz])/(36*dx)-(pC[c_xz]*pC[c_zz])/(12*dx)+(pC[c_xyz]*pC[c_zz])/(24*dx)+(pC[c_xy]*pC[c_zz])/(12*dx)-(pC[c_x]*pC[c_zz])/(6*dx)-(pC[c_xzz]*pC[c_z])/(12*dx)-(pC[c_xz]*pC[c_z])/(4*dx)+(pC[c_xyz]*pC[c_z])/(8*dx)+(pC[c_xy]*pC[c_z])/(4*dx)-(pC[c_x]*pC[c_z])/(2*dx)+(pC[c_xzz]*pC[c_yz])/(24*dx)+(pC[c_xz]*pC[c_yz])/(8*dx)-(pC[c_xyz]*pC[c_yz])/(16*dx)-(pC[c_xy]*pC[c_yz])/(8*dx)+(pC[c_x]*pC[c_yz])/(4*dx)+(pC[c_xzz]*pC[c_y])/(12*dx)+(pC[c_xz]*pC[c_y])/(4*dx)-(pC[c_xyz]*pC[c_y])/(8*dx)-(pC[c_xy]*pC[c_y])/(4*dx)+(pC[c_x]*pC[c_y])/(2*dx)-
   (pC[c_0]*pC[c_xzz])/(6*dx)-(pC[c_xxz]*pC[c_xz])/(24*dx)-(pC[c_xx]*pC[c_xz])/(12*dx)-(pC[c_0]*pC[c_xz])/(2*dx)+(pC[c_0]*pC[c_xyz])/(4*dx)+(pC[c_0]*pC[c_xy])/(2*dx)-(pC[c_x]*pC[c_xxz])/(12*dx)-(pC[c_x]*pC[c_xx])/(6*dx)-(pC[c_0]*pC[c_x])/dx-(pC[b_xz]*pC[b_z])/(4*dx)+(pC[b_xyz]*pC[b_z])/(8*dx)-(pC[b_xyy]*pC[b_z])/(12*dx)+(pC[b_xy]*pC[b_z])/(4*dx)-(pC[b_x]*pC[b_z])/(2*dx)+(pC[b_xz]*pC[b_yz])/(8*dx)-(pC[b_xyz]*pC[b_yz])/(16*dx)+(pC[b_xyy]*pC[b_yz])/(24*dx)-(pC[b_xy]*pC[b_yz])/(8*dx)+(pC[b_x]*pC[b_yz])/(4*dx)-(pC[b_xz]*pC[b_yy])/(12*dx)+
   (pC[b_xyz]*pC[b_yy])/(24*dx)-(pC[b_xyy]*pC[b_yy])/(36*dx)+(pC[b_xy]*pC[b_yy])/(12*dx)-(pC[b_x]*pC[b_yy])/(6*dx)+(pC[b_xz]*pC[b_y])/(4*dx)-(pC[b_xyz]*pC[b_y])/(8*dx)+(pC[b_xyy]*pC[b_y])/(12*dx)-(pC[b_xy]*pC[b_y])/(4*dx)+(pC[b_x]*pC[b_y])/(2*dx)-(pC[b_0]*pC[b_xz])/(2*dx)+(pC[b_0]*pC[b_xyz])/(4*dx)-(pC[b_0]*pC[b_xyy])/(6*dx)-(pC[b_xxy]*pC[b_xy])/(24*dx)+(pC[b_xx]*pC[b_xy])/(12*dx)+(pC[b_0]*pC[b_xy])/(2*dx)+(pC[b_x]*pC[b_xxy])/(12*dx)-(pC[b_x]*pC[b_xx])/(6*dx)-(pC[b_0]*pC[b_x])/dx ;
}

template<typename REAL> inline
REAL JXBX_011_111(
                  const REAL* const pC,
                  creal BGBY,
                  creal BGBZ,
                  creal dx,
                  creal dy,
                  creal dz
                 ) {
   using namespace Rec;
   return (pC[a_zz]*BGBZ)/dz+(pC[a_z]*BGBZ)/dz+(pC[a_yz]*BGBZ)/(2*dz)-(pC[c_xzz]*BGBZ)/(6*dx)-(pC[c_xz]*BGBZ)/(2*dx)-(pC[c_xyz]*BGBZ)/(4*dx)-(pC[c_xy]*BGBZ)/(2*dx)-(pC[c_x]*BGBZ)/dx+(pC[a_yz]*BGBY)/(2*dy)+(pC[a_yy]*BGBY)/dy+(pC[a_y]*BGBY)/dy-(pC[b_xz]*BGBY)/(2*dx)-(pC[b_xyz]*BGBY)/(4*dx)-(pC[b_xyy]*BGBY)/(6*dx)-(pC[b_xy]*BGBY)/(2*dx)-(pC[b_x]*BGBY)/dx+(pC[a_zz]*pC[c_zz])/(6*dz)+(pC[a_z]*pC[c_zz])/(6*dz)+
   (pC[a_yz]*pC[c_zz])/(12*dz)+(pC[a_zz]*pC[c_z])/(2*dz)+(pC[a_z]*pC[c_z])/(2*dz)+(pC[a_yz]*pC[c_z])/(4*dz)+(pC[a_zz]*pC[c_yz])/(4*dz)+(pC[a_z]*pC[c_yz])/(4*dz)+(pC[a_yz]*pC[c_yz])/(8*dz)+(pC[a_zz]*pC[c_y])/(2*dz)+(pC[a_z]*pC[c_y])/(2*dz)+(pC[a_yz]*pC[c_y])/(4*dz)+(pC[a_xzz]*pC[c_xz])/(24*dz)+(pC[a_xz]*pC[c_xz])/(24*dz)+(pC[a_xyz]*pC[c_xz])/(48*dz)+(pC[a_xzz]*pC[c_x])/(12*dz)+(pC[a_xz]*pC[c_x])/(12*dz)+(pC[a_xyz]*pC[c_x])/(24*dz)+(pC[a_zz]*pC[c_0])/dz+(pC[a_z]*pC[c_0])/dz+(pC[a_yz]*pC[c_0])/(2*dz)+(pC[a_yz]*pC[b_z])/(4*dy)+
   (pC[a_yy]*pC[b_z])/(2*dy)+(pC[a_y]*pC[b_z])/(2*dy)+(pC[a_yz]*pC[b_yz])/(8*dy)+(pC[a_yy]*pC[b_yz])/(4*dy)+(pC[a_y]*pC[b_yz])/(4*dy)+(pC[a_yz]*pC[b_yy])/(12*dy)+(pC[a_yy]*pC[b_yy])/(6*dy)+(pC[a_y]*pC[b_yy])/(6*dy)+(pC[a_yz]*pC[b_y])/(4*dy)+(pC[a_yy]*pC[b_y])/(2*dy)+(pC[a_y]*pC[b_y])/(2*dy)+(pC[a_xyz]*pC[b_xy])/(48*dy)+(pC[a_xyy]*pC[b_xy])/(24*dy)+(pC[a_xy]*pC[b_xy])/(24*dy)+(pC[a_xyz]*pC[b_x])/(24*dy)+(pC[a_xyy]*pC[b_x])/(12*dy)+(pC[a_xy]*pC[b_x])/(12*dy)+(pC[a_yz]*pC[b_0])/(2*dy)+(pC[a_yy]*pC[b_0])/dy+(pC[a_y]*pC[b_0])/dy-
   (pC[c_xzz]*pC[c_zz])/(36*dx)-(pC[c_xz]*pC[c_zz])/(12*dx)-(pC[c_xyz]*pC[c_zz])/(24*dx)-(pC[c_xy]*pC[c_zz])/(12*dx)-(pC[c_x]*pC[c_zz])/(6*dx)-(pC[c_xzz]*pC[c_z])/(12*dx)-(pC[c_xz]*pC[c_z])/(4*dx)-(pC[c_xyz]*pC[c_z])/(8*dx)-(pC[c_xy]*pC[c_z])/(4*dx)-(pC[c_x]*pC[c_z])/(2*dx)-(pC[c_xzz]*pC[c_yz])/(24*dx)-(pC[c_xz]*pC[c_yz])/(8*dx)-(pC[c_xyz]*pC[c_yz])/(16*dx)-(pC[c_xy]*pC[c_yz])/(8*dx)-(pC[c_x]*pC[c_yz])/(4*dx)-(pC[c_xzz]*pC[c_y])/(12*dx)-(pC[c_xz]*pC[c_y])/(4*dx)-(pC[c_xyz]*pC[c_y])/(8*dx)-(pC[c_xy]*pC[c_y])/(4*dx)-(pC[c_x]*pC[c_y])/(2*dx)-
   (pC[c_0]*pC[c_xzz])/(6*dx)-(pC[c_xxz]*pC[c_xz])/(24*dx)-(pC[c_xx]*pC[c_xz])/(12*dx)-(pC[c_0]*pC[c_xz])/(2*dx)-(pC[c_0]*pC[c_xyz])/(4*dx)-(pC[c_0]*pC[c_xy])/(2*dx)-(pC[c_x]*pC[c_xxz])/(12*dx)-(pC[c_x]*pC[c_xx])/(6*dx)-(pC[c_0]*pC[c_x])/dx-(pC[b_xz]*pC[b_z])/(4*dx)-(pC[b_xyz]*pC[b_z])/(8*dx)-(pC[b_xyy]*pC[b_z])/(12*dx)-(pC[b_xy]*pC[b_z])/(4*dx)-(pC[b_x]*pC[b_z])/(2*dx)-(pC[b_xz]*pC[b_yz])/(8*dx)-(pC[b_xyz]*pC[b_yz])/(16*dx)-(pC[b_xyy]*pC[b_yz])/(24*dx)-(pC[b_xy]*pC[b_yz])/(8*dx)-(pC[b_x]*pC[b_yz])/(4*dx)-(pC[b_xz]*pC[b_yy])/(12*dx)-
   (pC[b_xyz]*pC[b_yy])/(24*dx)-(pC[b_xyy]*pC[b_yy])/(36*dx)-(pC[b_xy]*pC[b_yy])/(12*dx)-(pC[b_x]*pC[b_yy])/(6*dx)-(pC[b_xz]*pC[b_y])/(4*dx)-(pC[b_xyz]*pC[b_y])/(8*dx)-(pC[b_xyy]*pC[b_y])/(12*dx)-(pC[b_xy]*pC[b_y])/(4*dx)-(pC[b_x]*pC[b_y])/(2*dx)-(pC[b_0]*pC[b_xz])/(2*dx)-(pC[b_0]*pC[b_xyz])/(4*dx)-(pC[b_0]*pC[b_xyy])/(6*dx)-(pC[b_xxy]*pC[b_xy])/(24*dx)-(pC[b_xx]*pC[b_xy])/(12*dx)-(pC[b_0]*pC[b_xy])/(2*dx)-(pC[b_x]*pC[b_xxy])/(12*dx)-(pC[b_x]*pC[b_xx])/(6*dx)-(pC[b_0]*pC[b_x])/dx;
}

// Y
template<typename REAL> inline
REAL JXBY_000_010(
                  const REAL* const pC,
                  creal BGBX,
                  creal BGBZ,
                  creal dx,
                  creal dy,
                  creal dz
                 ) {
   using namespace Rec;
   return -(pC[b_zz]*BGBZ)/dz+(pC[b_z]*BGBZ)/dz-(pC[b_xz]*BGBZ)/(2*dz)-(pC[c_yzz]*BGBZ)/(6*dy)+(pC[c_yz]*BGBZ)/(2*dy)-(pC[c_y]*BGBZ)/dy-(pC[c_xyz]*BGBZ)/(4*dy)+(pC[c_xy]*BGBZ)/(2*dy)+(pC[a_yz]*BGBX)/(2*dy)-(pC[a_y]*BGBX)/dy-(pC[a_xyz]*BGBX)/(4*dy)+(pC[a_xy]*BGBX)/(2*dy)-(pC[a_xxy]*BGBX)/(6*dy)-(pC[b_xz]*BGBX)/(2*dx)-(pC[b_xx]*BGBX)/dx+(pC[b_x]*BGBX)/dx-(pC[b_zz]*pC[c_zz])/(6*dz)+(pC[b_z]*pC[c_zz])/(6*dz)-
     (pC[b_xz]*pC[c_zz])/(12*dz)+(pC[b_zz]*pC[c_z])/(2*dz)-(pC[b_z]*pC[c_z])/(2*dz)+(pC[b_xz]*pC[c_z])/(4*dz)+(pC[b_yzz]*pC[c_yz])/(24*dz)-(pC[b_yz]*pC[c_yz])/(24*dz)+(pC[b_xyz]*pC[c_yz])/(48*dz)-(pC[b_yzz]*pC[c_y])/(12*dz)+(pC[b_yz]*pC[c_y])/(12*dz)-(pC[b_xyz]*pC[c_y])/(24*dz)-(pC[b_zz]*pC[c_xz])/(4*dz)+(pC[b_z]*pC[c_xz])/(4*dz)-(pC[b_xz]*pC[c_xz])/(8*dz)+(pC[b_zz]*pC[c_x])/(2*dz)-(pC[b_z]*pC[c_x])/(2*dz)+(pC[b_xz]*pC[c_x])/(4*dz)-(pC[b_zz]*pC[c_0])/dz+(pC[b_z]*pC[c_0])/dz-(pC[b_xz]*pC[c_0])/(2*dz)-(pC[c_yzz]*pC[c_zz])/(36*dy)+
     (pC[c_yz]*pC[c_zz])/(12*dy)-(pC[c_y]*pC[c_zz])/(6*dy)-(pC[c_xyz]*pC[c_zz])/(24*dy)+(pC[c_xy]*pC[c_zz])/(12*dy)+(pC[c_yzz]*pC[c_z])/(12*dy)-(pC[c_yz]*pC[c_z])/(4*dy)+(pC[c_y]*pC[c_z])/(2*dy)+(pC[c_xyz]*pC[c_z])/(8*dy)-(pC[c_xy]*pC[c_z])/(4*dy)-(pC[c_xz]*pC[c_yzz])/(24*dy)+(pC[c_x]*pC[c_yzz])/(12*dy)-(pC[c_0]*pC[c_yzz])/(6*dy)-(pC[c_yyz]*pC[c_yz])/(24*dy)+(pC[c_yy]*pC[c_yz])/(12*dy)+(pC[c_xz]*pC[c_yz])/(8*dy)-(pC[c_x]*pC[c_yz])/(4*dy)+(pC[c_0]*pC[c_yz])/(2*dy)+(pC[c_y]*pC[c_yyz])/(12*dy)-(pC[c_y]*pC[c_yy])/(6*dy)-(pC[c_xz]*pC[c_y])/(4*dy)+
     (pC[c_x]*pC[c_y])/(2*dy)-(pC[c_0]*pC[c_y])/dy-(pC[c_xyz]*pC[c_xz])/(16*dy)+(pC[c_xy]*pC[c_xz])/(8*dy)+(pC[c_x]*pC[c_xyz])/(8*dy)-(pC[c_0]*pC[c_xyz])/(4*dy)-(pC[c_x]*pC[c_xy])/(4*dy)+(pC[c_0]*pC[c_xy])/(2*dy)-(pC[a_yz]*pC[a_z])/(4*dy)+(pC[a_y]*pC[a_z])/(2*dy)+(pC[a_xyz]*pC[a_z])/(8*dy)-(pC[a_xy]*pC[a_z])/(4*dy)+(pC[a_xxy]*pC[a_z])/(12*dy)+(pC[a_xz]*pC[a_yz])/(8*dy)+(pC[a_xx]*pC[a_yz])/(12*dy)-(pC[a_x]*pC[a_yz])/(4*dy)+(pC[a_0]*pC[a_yz])/(2*dy)-(pC[a_y]*pC[a_yy])/(6*dy)+(pC[a_xy]*pC[a_yy])/(12*dy)-(pC[a_xz]*pC[a_y])/(4*dy)+
     (pC[a_xyy]*pC[a_y])/(12*dy)-(pC[a_xx]*pC[a_y])/(6*dy)+(pC[a_x]*pC[a_y])/(2*dy)-(pC[a_0]*pC[a_y])/dy-(pC[a_xyz]*pC[a_xz])/(16*dy)+(pC[a_xy]*pC[a_xz])/(8*dy)-(pC[a_xxy]*pC[a_xz])/(24*dy)-(pC[a_xx]*pC[a_xyz])/(24*dy)+(pC[a_x]*pC[a_xyz])/(8*dy)-(pC[a_0]*pC[a_xyz])/(4*dy)-(pC[a_xy]*pC[a_xyy])/(24*dy)+(pC[a_xx]*pC[a_xy])/(12*dy)-(pC[a_x]*pC[a_xy])/(4*dy)+(pC[a_0]*pC[a_xy])/(2*dy)-(pC[a_xx]*pC[a_xxy])/(36*dy)+(pC[a_x]*pC[a_xxy])/(12*dy)-(pC[a_0]*pC[a_xxy])/(6*dy)+(pC[a_z]*pC[b_xz])/(4*dx)-(pC[a_xz]*pC[b_xz])/(8*dx)-
     (pC[a_xx]*pC[b_xz])/(12*dx)+(pC[a_x]*pC[b_xz])/(4*dx)-(pC[a_0]*pC[b_xz])/(2*dx)-(pC[a_y]*pC[b_xyz])/(24*dx)+(pC[a_xy]*pC[b_xyz])/(48*dx)+(pC[a_y]*pC[b_xy])/(12*dx)-(pC[a_xy]*pC[b_xy])/(24*dx)-(pC[a_y]*pC[b_xxy])/(12*dx)+(pC[a_xy]*pC[b_xxy])/(24*dx)+(pC[a_z]*pC[b_xx])/(2*dx)-(pC[a_xz]*pC[b_xx])/(4*dx)-(pC[a_xx]*pC[b_xx])/(6*dx)+(pC[a_x]*pC[b_xx])/(2*dx)-(pC[a_0]*pC[b_xx])/dx-(pC[a_z]*pC[b_x])/(2*dx)+(pC[a_xz]*pC[b_x])/(4*dx)+(pC[a_xx]*pC[b_x])/(6*dx)-(pC[a_x]*pC[b_x])/(2*dx)+(pC[a_0]*pC[b_x])/dx;
}

template<typename REAL> inline
REAL JXBY_100_110(
                  const REAL* const pC,
                  creal BGBX,
                  creal BGBZ,
                  creal dx,
                  creal dy,
                  creal dz
                 ) {
   using namespace Rec;
   return -(pC[b_zz]*BGBZ)/dz+(pC[b_z]*BGBZ)/dz+(pC[b_xz]*BGBZ)/(2*dz)-(pC[c_yzz]*BGBZ)/(6*dy)+(pC[c_yz]*BGBZ)/(2*dy)-(pC[c_y]*BGBZ)/dy+(pC[c_xyz]*BGBZ)/(4*dy)-(pC[c_xy]*BGBZ)/(2*dy)+(pC[a_yz]*BGBX)/(2*dy)-(pC[a_y]*BGBX)/dy+(pC[a_xyz]*BGBX)/(4*dy)-(pC[a_xy]*BGBX)/(2*dy)-(pC[a_xxy]*BGBX)/(6*dy)-(pC[b_xz]*BGBX)/(2*dx)+(pC[b_xx]*BGBX)/dx+(pC[b_x]*BGBX)/dx-(pC[b_zz]*pC[c_zz])/(6*dz)+(pC[b_z]*pC[c_zz])/(6*dz)+
     (pC[b_xz]*pC[c_zz])/(12*dz)+(pC[b_zz]*pC[c_z])/(2*dz)-(pC[b_z]*pC[c_z])/(2*dz)-(pC[b_xz]*pC[c_z])/(4*dz)+(pC[b_yzz]*pC[c_yz])/(24*dz)-(pC[b_yz]*pC[c_yz])/(24*dz)-(pC[b_xyz]*pC[c_yz])/(48*dz)-(pC[b_yzz]*pC[c_y])/(12*dz)+(pC[b_yz]*pC[c_y])/(12*dz)+(pC[b_xyz]*pC[c_y])/(24*dz)+(pC[b_zz]*pC[c_xz])/(4*dz)-(pC[b_z]*pC[c_xz])/(4*dz)-(pC[b_xz]*pC[c_xz])/(8*dz)-(pC[b_zz]*pC[c_x])/(2*dz)+(pC[b_z]*pC[c_x])/(2*dz)+(pC[b_xz]*pC[c_x])/(4*dz)-(pC[b_zz]*pC[c_0])/dz+(pC[b_z]*pC[c_0])/dz+(pC[b_xz]*pC[c_0])/(2*dz)-(pC[c_yzz]*pC[c_zz])/(36*dy)+
     (pC[c_yz]*pC[c_zz])/(12*dy)-(pC[c_y]*pC[c_zz])/(6*dy)+(pC[c_xyz]*pC[c_zz])/(24*dy)-(pC[c_xy]*pC[c_zz])/(12*dy)+(pC[c_yzz]*pC[c_z])/(12*dy)-(pC[c_yz]*pC[c_z])/(4*dy)+(pC[c_y]*pC[c_z])/(2*dy)-(pC[c_xyz]*pC[c_z])/(8*dy)+(pC[c_xy]*pC[c_z])/(4*dy)+(pC[c_xz]*pC[c_yzz])/(24*dy)-(pC[c_x]*pC[c_yzz])/(12*dy)-(pC[c_0]*pC[c_yzz])/(6*dy)-(pC[c_yyz]*pC[c_yz])/(24*dy)+(pC[c_yy]*pC[c_yz])/(12*dy)-(pC[c_xz]*pC[c_yz])/(8*dy)+(pC[c_x]*pC[c_yz])/(4*dy)+(pC[c_0]*pC[c_yz])/(2*dy)+(pC[c_y]*pC[c_yyz])/(12*dy)-(pC[c_y]*pC[c_yy])/(6*dy)+(pC[c_xz]*pC[c_y])/(4*dy)-
     (pC[c_x]*pC[c_y])/(2*dy)-(pC[c_0]*pC[c_y])/dy-(pC[c_xyz]*pC[c_xz])/(16*dy)+(pC[c_xy]*pC[c_xz])/(8*dy)+(pC[c_x]*pC[c_xyz])/(8*dy)+(pC[c_0]*pC[c_xyz])/(4*dy)-(pC[c_x]*pC[c_xy])/(4*dy)-(pC[c_0]*pC[c_xy])/(2*dy)-(pC[a_yz]*pC[a_z])/(4*dy)+(pC[a_y]*pC[a_z])/(2*dy)-(pC[a_xyz]*pC[a_z])/(8*dy)+(pC[a_xy]*pC[a_z])/(4*dy)+(pC[a_xxy]*pC[a_z])/(12*dy)-(pC[a_xz]*pC[a_yz])/(8*dy)+(pC[a_xx]*pC[a_yz])/(12*dy)+(pC[a_x]*pC[a_yz])/(4*dy)+(pC[a_0]*pC[a_yz])/(2*dy)-(pC[a_y]*pC[a_yy])/(6*dy)-(pC[a_xy]*pC[a_yy])/(12*dy)+(pC[a_xz]*pC[a_y])/(4*dy)-
     (pC[a_xyy]*pC[a_y])/(12*dy)-(pC[a_xx]*pC[a_y])/(6*dy)-(pC[a_x]*pC[a_y])/(2*dy)-(pC[a_0]*pC[a_y])/dy-(pC[a_xyz]*pC[a_xz])/(16*dy)+(pC[a_xy]*pC[a_xz])/(8*dy)+(pC[a_xxy]*pC[a_xz])/(24*dy)+(pC[a_xx]*pC[a_xyz])/(24*dy)+(pC[a_x]*pC[a_xyz])/(8*dy)+(pC[a_0]*pC[a_xyz])/(4*dy)-(pC[a_xy]*pC[a_xyy])/(24*dy)-(pC[a_xx]*pC[a_xy])/(12*dy)-(pC[a_x]*pC[a_xy])/(4*dy)-(pC[a_0]*pC[a_xy])/(2*dy)-(pC[a_xx]*pC[a_xxy])/(36*dy)-(pC[a_x]*pC[a_xxy])/(12*dy)-(pC[a_0]*pC[a_xxy])/(6*dy)+(pC[a_z]*pC[b_xz])/(4*dx)+(pC[a_xz]*pC[b_xz])/(8*dx)-
     (pC[a_xx]*pC[b_xz])/(12*dx)-(pC[a_x]*pC[b_xz])/(4*dx)-(pC[a_0]*pC[b_xz])/(2*dx)-(pC[a_y]*pC[b_xyz])/(24*dx)-(pC[a_xy]*pC[b_xyz])/(48*dx)+(pC[a_y]*pC[b_xy])/(12*dx)+(pC[a_xy]*pC[b_xy])/(24*dx)+(pC[a_y]*pC[b_xxy])/(12*dx)+(pC[a_xy]*pC[b_xxy])/(24*dx)-(pC[a_z]*pC[b_xx])/(2*dx)-(pC[a_xz]*pC[b_xx])/(4*dx)+(pC[a_xx]*pC[b_xx])/(6*dx)+(pC[a_x]*pC[b_xx])/(2*dx)+(pC[a_0]*pC[b_xx])/dx-(pC[a_z]*pC[b_x])/(2*dx)-(pC[a_xz]*pC[b_x])/(4*dx)+(pC[a_xx]*pC[b_x])/(6*dx)+(pC[a_x]*pC[b_x])/(2*dx)+(pC[a_0]*pC[b_x])/dx;
}

template<typename REAL> inline
REAL JXBY_001_011(
                  const REAL* const pC,
                  creal BGBX,
                  creal BGBZ,
                  creal dx,
                  creal dy,
                  creal dz
                 ) {
   using namespace Rec;
   return (pC[b_zz]*BGBZ)/dz+(pC[b_z]*BGBZ)/dz-(pC[b_xz]*BGBZ)/(2*dz)-(pC[c_yzz]*BGBZ)/(6*dy)-(pC[c_yz]*BGBZ)/(2*dy)-(pC[c_y]*BGBZ)/dy+(pC[c_xyz]*BGBZ)/(4*dy)+(pC[c_xy]*BGBZ)/(2*dy)-(pC[a_yz]*BGBX)/(2*dy)-(pC[a_y]*BGBX)/dy+(pC[a_xyz]*BGBX)/(4*dy)+(pC[a_xy]*BGBX)/(2*dy)-(pC[a_xxy]*BGBX)/(6*dy)+(pC[b_xz]*BGBX)/(2*dx)-(pC[b_xx]*BGBX)/dx+(pC[b_x]*BGBX)/dx+(pC[b_zz]*pC[c_zz])/(6*dz)+(pC[b_z]*pC[c_zz])/(6*dz)-
     (pC[b_xz]*pC[c_zz])/(12*dz)+(pC[b_zz]*pC[c_z])/(2*dz)+(pC[b_z]*pC[c_z])/(2*dz)-(pC[b_xz]*pC[c_z])/(4*dz)+(pC[b_yzz]*pC[c_yz])/(24*dz)+(pC[b_yz]*pC[c_yz])/(24*dz)-(pC[b_xyz]*pC[c_yz])/(48*dz)+(pC[b_yzz]*pC[c_y])/(12*dz)+(pC[b_yz]*pC[c_y])/(12*dz)-(pC[b_xyz]*pC[c_y])/(24*dz)-(pC[b_zz]*pC[c_xz])/(4*dz)-(pC[b_z]*pC[c_xz])/(4*dz)+(pC[b_xz]*pC[c_xz])/(8*dz)-(pC[b_zz]*pC[c_x])/(2*dz)-(pC[b_z]*pC[c_x])/(2*dz)+(pC[b_xz]*pC[c_x])/(4*dz)+(pC[b_zz]*pC[c_0])/dz+(pC[b_z]*pC[c_0])/dz-(pC[b_xz]*pC[c_0])/(2*dz)-(pC[c_yzz]*pC[c_zz])/(36*dy)-
     (pC[c_yz]*pC[c_zz])/(12*dy)-(pC[c_y]*pC[c_zz])/(6*dy)+(pC[c_xyz]*pC[c_zz])/(24*dy)+(pC[c_xy]*pC[c_zz])/(12*dy)-(pC[c_yzz]*pC[c_z])/(12*dy)-(pC[c_yz]*pC[c_z])/(4*dy)-(pC[c_y]*pC[c_z])/(2*dy)+(pC[c_xyz]*pC[c_z])/(8*dy)+(pC[c_xy]*pC[c_z])/(4*dy)+(pC[c_xz]*pC[c_yzz])/(24*dy)+(pC[c_x]*pC[c_yzz])/(12*dy)-(pC[c_0]*pC[c_yzz])/(6*dy)-(pC[c_yyz]*pC[c_yz])/(24*dy)-(pC[c_yy]*pC[c_yz])/(12*dy)+(pC[c_xz]*pC[c_yz])/(8*dy)+(pC[c_x]*pC[c_yz])/(4*dy)-(pC[c_0]*pC[c_yz])/(2*dy)-(pC[c_y]*pC[c_yyz])/(12*dy)-(pC[c_y]*pC[c_yy])/(6*dy)+(pC[c_xz]*pC[c_y])/(4*dy)+
     (pC[c_x]*pC[c_y])/(2*dy)-(pC[c_0]*pC[c_y])/dy-(pC[c_xyz]*pC[c_xz])/(16*dy)-(pC[c_xy]*pC[c_xz])/(8*dy)-(pC[c_x]*pC[c_xyz])/(8*dy)+(pC[c_0]*pC[c_xyz])/(4*dy)-(pC[c_x]*pC[c_xy])/(4*dy)+(pC[c_0]*pC[c_xy])/(2*dy)-(pC[a_yz]*pC[a_z])/(4*dy)-(pC[a_y]*pC[a_z])/(2*dy)+(pC[a_xyz]*pC[a_z])/(8*dy)+(pC[a_xy]*pC[a_z])/(4*dy)-(pC[a_xxy]*pC[a_z])/(12*dy)+(pC[a_xz]*pC[a_yz])/(8*dy)-(pC[a_xx]*pC[a_yz])/(12*dy)+(pC[a_x]*pC[a_yz])/(4*dy)-(pC[a_0]*pC[a_yz])/(2*dy)-(pC[a_y]*pC[a_yy])/(6*dy)+(pC[a_xy]*pC[a_yy])/(12*dy)+(pC[a_xz]*pC[a_y])/(4*dy)+
     (pC[a_xyy]*pC[a_y])/(12*dy)-(pC[a_xx]*pC[a_y])/(6*dy)+(pC[a_x]*pC[a_y])/(2*dy)-(pC[a_0]*pC[a_y])/dy-(pC[a_xyz]*pC[a_xz])/(16*dy)-(pC[a_xy]*pC[a_xz])/(8*dy)+(pC[a_xxy]*pC[a_xz])/(24*dy)+(pC[a_xx]*pC[a_xyz])/(24*dy)-(pC[a_x]*pC[a_xyz])/(8*dy)+(pC[a_0]*pC[a_xyz])/(4*dy)-(pC[a_xy]*pC[a_xyy])/(24*dy)+(pC[a_xx]*pC[a_xy])/(12*dy)-(pC[a_x]*pC[a_xy])/(4*dy)+(pC[a_0]*pC[a_xy])/(2*dy)-(pC[a_xx]*pC[a_xxy])/(36*dy)+(pC[a_x]*pC[a_xxy])/(12*dy)-(pC[a_0]*pC[a_xxy])/(6*dy)+(pC[a_z]*pC[b_xz])/(4*dx)-(pC[a_xz]*pC[b_xz])/(8*dx)+
     (pC[a_xx]*pC[b_xz])/(12*dx)-(pC[a_x]*pC[b_xz])/(4*dx)+(pC[a_0]*pC[b_xz])/(2*dx)+(pC[a_y]*pC[b_xyz])/(24*dx)-(pC[a_xy]*pC[b_xyz])/(48*dx)+(pC[a_y]*pC[b_xy])/(12*dx)-(pC[a_xy]*pC[b_xy])/(24*dx)-(pC[a_y]*pC[b_xxy])/(12*dx)+(pC[a_xy]*pC[b_xxy])/(24*dx)-(pC[a_z]*pC[b_xx])/(2*dx)+(pC[a_xz]*pC[b_xx])/(4*dx)-(pC[a_xx]*pC[b_xx])/(6*dx)+(pC[a_x]*pC[b_xx])/(2*dx)-(pC[a_0]*pC[b_xx])/dx+(pC[a_z]*pC[b_x])/(2*dx)-(pC[a_xz]*pC[b_x])/(4*dx)+(pC[a_xx]*pC[b_x])/(6*dx)-(pC[a_x]*pC[b_x])/(2*dx)+(pC[a_0]*pC[b_x])/dx;
}

template<typename REAL> inline
REAL JXBY_101_111(
                  const REAL* const pC,
                  creal BGBX,
                  creal BGBZ,
                  creal dx,
                  creal dy,
                  creal dz
                 ) {
   using namespace Rec;
   return (pC[b_zz]*BGBZ)/dz+(pC[b_z]*BGBZ)/dz+(pC[b_xz]*BGBZ)/(2*dz)-(pC[c_yzz]*BGBZ)/(6*dy)-(pC[c_yz]*BGBZ)/(2*dy)-(pC[c_y]*BGBZ)/dy-(pC[c_xyz]*BGBZ)/(4*dy)-(pC[c_xy]*BGBZ)/(2*dy)-(pC[a_yz]*BGBX)/(2*dy)-(pC[a_y]*BGBX)/dy-(pC[a_xyz]*BGBX)/(4*dy)-(pC[a_xy]*BGBX)/(2*dy)-(pC[a_xxy]*BGBX)/(6*dy)+(pC[b_xz]*BGBX)/(2*dx)+(pC[b_xx]*BGBX)/dx+(pC[b_x]*BGBX)/dx+(pC[b_zz]*pC[c_zz])/(6*dz)+(pC[b_z]*pC[c_zz])/(6*dz)+
     (pC[b_xz]*pC[c_zz])/(12*dz)+(pC[b_zz]*pC[c_z])/(2*dz)+(pC[b_z]*pC[c_z])/(2*dz)+(pC[b_xz]*pC[c_z])/(4*dz)+(pC[b_yzz]*pC[c_yz])/(24*dz)+(pC[b_yz]*pC[c_yz])/(24*dz)+(pC[b_xyz]*pC[c_yz])/(48*dz)+(pC[b_yzz]*pC[c_y])/(12*dz)+(pC[b_yz]*pC[c_y])/(12*dz)+(pC[b_xyz]*pC[c_y])/(24*dz)+(pC[b_zz]*pC[c_xz])/(4*dz)+(pC[b_z]*pC[c_xz])/(4*dz)+(pC[b_xz]*pC[c_xz])/(8*dz)+(pC[b_zz]*pC[c_x])/(2*dz)+(pC[b_z]*pC[c_x])/(2*dz)+(pC[b_xz]*pC[c_x])/(4*dz)+(pC[b_zz]*pC[c_0])/dz+(pC[b_z]*pC[c_0])/dz+(pC[b_xz]*pC[c_0])/(2*dz)-(pC[c_yzz]*pC[c_zz])/(36*dy)-
     (pC[c_yz]*pC[c_zz])/(12*dy)-(pC[c_y]*pC[c_zz])/(6*dy)-(pC[c_xyz]*pC[c_zz])/(24*dy)-(pC[c_xy]*pC[c_zz])/(12*dy)-(pC[c_yzz]*pC[c_z])/(12*dy)-(pC[c_yz]*pC[c_z])/(4*dy)-(pC[c_y]*pC[c_z])/(2*dy)-(pC[c_xyz]*pC[c_z])/(8*dy)-(pC[c_xy]*pC[c_z])/(4*dy)-(pC[c_xz]*pC[c_yzz])/(24*dy)-(pC[c_x]*pC[c_yzz])/(12*dy)-(pC[c_0]*pC[c_yzz])/(6*dy)-(pC[c_yyz]*pC[c_yz])/(24*dy)-(pC[c_yy]*pC[c_yz])/(12*dy)-(pC[c_xz]*pC[c_yz])/(8*dy)-(pC[c_x]*pC[c_yz])/(4*dy)-(pC[c_0]*pC[c_yz])/(2*dy)-(pC[c_y]*pC[c_yyz])/(12*dy)-(pC[c_y]*pC[c_yy])/(6*dy)-(pC[c_xz]*pC[c_y])/(4*dy)-
     (pC[c_x]*pC[c_y])/(2*dy)-(pC[c_0]*pC[c_y])/dy-(pC[c_xyz]*pC[c_xz])/(16*dy)-(pC[c_xy]*pC[c_xz])/(8*dy)-(pC[c_x]*pC[c_xyz])/(8*dy)-(pC[c_0]*pC[c_xyz])/(4*dy)-(pC[c_x]*pC[c_xy])/(4*dy)-(pC[c_0]*pC[c_xy])/(2*dy)-(pC[a_yz]*pC[a_z])/(4*dy)-(pC[a_y]*pC[a_z])/(2*dy)-(pC[a_xyz]*pC[a_z])/(8*dy)-(pC[a_xy]*pC[a_z])/(4*dy)-(pC[a_xxy]*pC[a_z])/(12*dy)-(pC[a_xz]*pC[a_yz])/(8*dy)-(pC[a_xx]*pC[a_yz])/(12*dy)-(pC[a_x]*pC[a_yz])/(4*dy)-(pC[a_0]*pC[a_yz])/(2*dy)-(pC[a_y]*pC[a_yy])/(6*dy)-(pC[a_xy]*pC[a_yy])/(12*dy)-(pC[a_xz]*pC[a_y])/(4*dy)-(pC[a_xyy]*pC[a_y])/(12*dy)-
     (pC[a_xx]*pC[a_y])/(6*dy)-(pC[a_x]*pC[a_y])/(2*dy)-(pC[a_0]*pC[a_y])/dy-(pC[a_xyz]*pC[a_xz])/(16*dy)-(pC[a_xy]*pC[a_xz])/(8*dy)-(pC[a_xxy]*pC[a_xz])/(24*dy)-(pC[a_xx]*pC[a_xyz])/(24*dy)-(pC[a_x]*pC[a_xyz])/(8*dy)-(pC[a_0]*pC[a_xyz])/(4*dy)-(pC[a_xy]*pC[a_xyy])/(24*dy)-(pC[a_xx]*pC[a_xy])/(12*dy)-(pC[a_x]*pC[a_xy])/(4*dy)-(pC[a_0]*pC[a_xy])/(2*dy)-(pC[a_xx]*pC[a_xxy])/(36*dy)-(pC[a_x]*pC[a_xxy])/(12*dy)-(pC[a_0]*pC[a_xxy])/(6*dy)+(pC[a_z]*pC[b_xz])/(4*dx)+(pC[a_xz]*pC[b_xz])/(8*dx)+(pC[a_xx]*pC[b_xz])/(12*dx)+(pC[a_x]*pC[b_xz])/(4*dx)+
     (pC[a_0]*pC[b_xz])/(2*dx)+(pC[a_y]*pC[b_xyz])/(24*dx)+(pC[a_xy]*pC[b_xyz])/(48*dx)+(pC[a_y]*pC[b_xy])/(12*dx)+(pC[a_xy]*pC[b_xy])/(24*dx)+(pC[a_y]*pC[b_xxy])/(12*dx)+(pC[a_xy]*pC[b_xxy])/(24*dx)+(pC[a_z]*pC[b_xx])/(2*dx)+(pC[a_xz]*pC[b_xx])/(4*dx)+(pC[a_xx]*pC[b_xx])/(6*dx)+(pC[a_x]*pC[b_xx])/(2*dx)+(pC[a_0]*pC[b_xx])/dx+(pC[a_z]*pC[b_x])/(2*dx)+(pC[a_xz]*pC[b_x])/(4*dx)+(pC[a_xx]*pC[b_x])/(6*dx)+(pC[a_x]*pC[b_x])/(2*dx)+(pC[a_0]*pC[b_x])/dx;
}

// Z
template<typename REAL> inline
REAL JXBZ_000_001(
                  const REAL* const pC,
                  creal BGBX,
                  creal BGBY,
                  creal dx,
                  creal dy,
                  creal dz
                 ) {
   using namespace Rec;
   return -(pC[b_z]*BGBY)/dz+(pC[b_yz]*BGBY)/(2*dz)-(pC[b_yyz]*BGBY)/(6*dz)+(pC[b_xz]*BGBY)/(2*dz)-(pC[b_xyz]*BGBY)/(4*dz)-(pC[c_yy]*BGBY)/dy+(pC[c_y]*BGBY)/dy-(pC[c_xy]*BGBY)/(2*dy)-(pC[a_z]*BGBX)/dz+(pC[a_yz]*BGBX)/(2*dz)+(pC[a_xz]*BGBX)/(2*dz)-(pC[a_xyz]*BGBX)/(4*dz)-(pC[a_xxz]*BGBX)/(6*dz)-(pC[c_xy]*BGBX)/(2*dx)-(pC[c_xx]*BGBX)/dx+(pC[c_x]*BGBX)/dx-(pC[b_z]*pC[b_zz])/(6*dz)+(pC[b_yz]*pC[b_zz])/(12*dz)+
     (pC[b_yzz]*pC[b_z])/(12*dz)-(pC[b_yy]*pC[b_z])/(6*dz)+(pC[b_y]*pC[b_z])/(2*dz)-(pC[b_xy]*pC[b_z])/(4*dz)+(pC[b_x]*pC[b_z])/(2*dz)-(pC[b_0]*pC[b_z])/dz-(pC[b_yz]*pC[b_yzz])/(24*dz)+(pC[b_yy]*pC[b_yz])/(12*dz)-(pC[b_y]*pC[b_yz])/(4*dz)+(pC[b_xy]*pC[b_yz])/(8*dz)-(pC[b_x]*pC[b_yz])/(4*dz)+(pC[b_0]*pC[b_yz])/(2*dz)-(pC[b_yy]*pC[b_yyz])/(36*dz)+(pC[b_y]*pC[b_yyz])/(12*dz)-(pC[b_xy]*pC[b_yyz])/(24*dz)+(pC[b_x]*pC[b_yyz])/(12*dz)-(pC[b_0]*pC[b_yyz])/(6*dz)+(pC[b_xz]*pC[b_yy])/(12*dz)-(pC[b_xyz]*pC[b_yy])/(24*dz)-(pC[b_xz]*pC[b_y])/(4*dz)+
     (pC[b_xyz]*pC[b_y])/(8*dz)+(pC[b_xy]*pC[b_xz])/(8*dz)-(pC[b_x]*pC[b_xz])/(4*dz)+(pC[b_0]*pC[b_xz])/(2*dz)-(pC[b_xy]*pC[b_xyz])/(16*dz)+(pC[b_x]*pC[b_xyz])/(8*dz)-(pC[b_0]*pC[b_xyz])/(4*dz)-(pC[a_z]*pC[a_zz])/(6*dz)+(pC[a_xz]*pC[a_zz])/(12*dz)+(pC[a_y]*pC[a_z])/(2*dz)+(pC[a_xzz]*pC[a_z])/(12*dz)-(pC[a_xy]*pC[a_z])/(4*dz)-(pC[a_xx]*pC[a_z])/(6*dz)+(pC[a_x]*pC[a_z])/(2*dz)-(pC[a_0]*pC[a_z])/dz-(pC[a_y]*pC[a_yz])/(4*dz)+(pC[a_xy]*pC[a_yz])/(8*dz)+(pC[a_xx]*pC[a_yz])/(12*dz)-(pC[a_x]*pC[a_yz])/(4*dz)+(pC[a_0]*pC[a_yz])/(2*dz)-
     (pC[a_xz]*pC[a_y])/(4*dz)+(pC[a_xyz]*pC[a_y])/(8*dz)+(pC[a_xxz]*pC[a_y])/(12*dz)-(pC[a_xz]*pC[a_xzz])/(24*dz)+(pC[a_xy]*pC[a_xz])/(8*dz)+(pC[a_xx]*pC[a_xz])/(12*dz)-(pC[a_x]*pC[a_xz])/(4*dz)+(pC[a_0]*pC[a_xz])/(2*dz)-(pC[a_xy]*pC[a_xyz])/(16*dz)-(pC[a_xx]*pC[a_xyz])/(24*dz)+(pC[a_x]*pC[a_xyz])/(8*dz)-(pC[a_0]*pC[a_xyz])/(4*dz)-(pC[a_xxz]*pC[a_xy])/(24*dz)-(pC[a_xx]*pC[a_xxz])/(36*dz)+(pC[a_x]*pC[a_xxz])/(12*dz)-(pC[a_0]*pC[a_xxz])/(6*dz)+(pC[b_z]*pC[c_yz])/(12*dy)-(pC[b_yz]*pC[c_yz])/(24*dy)-(pC[b_z]*pC[c_yyz])/(12*dy)+
     (pC[b_yz]*pC[c_yyz])/(24*dy)-(pC[b_yy]*pC[c_yy])/(6*dy)+(pC[b_y]*pC[c_yy])/(2*dy)-(pC[b_xy]*pC[c_yy])/(4*dy)+(pC[b_x]*pC[c_yy])/(2*dy)-(pC[b_0]*pC[c_yy])/dy+(pC[b_yy]*pC[c_y])/(6*dy)-(pC[b_y]*pC[c_y])/(2*dy)+(pC[b_xy]*pC[c_y])/(4*dy)-(pC[b_x]*pC[c_y])/(2*dy)+(pC[b_0]*pC[c_y])/dy-(pC[b_z]*pC[c_xyz])/(24*dy)+(pC[b_yz]*pC[c_xyz])/(48*dy)-(pC[b_yy]*pC[c_xy])/(12*dy)+(pC[b_y]*pC[c_xy])/(4*dy)-(pC[b_xy]*pC[c_xy])/(8*dy)+(pC[b_x]*pC[c_xy])/(4*dy)-(pC[b_0]*pC[c_xy])/(2*dy)+(pC[a_z]*pC[c_xz])/(12*dx)-(pC[a_xz]*pC[c_xz])/(24*dx)-
     (pC[a_z]*pC[c_xyz])/(24*dx)+(pC[a_xz]*pC[c_xyz])/(48*dx)+(pC[a_y]*pC[c_xy])/(4*dx)-(pC[a_xy]*pC[c_xy])/(8*dx)-(pC[a_xx]*pC[c_xy])/(12*dx)+(pC[a_x]*pC[c_xy])/(4*dx)-(pC[a_0]*pC[c_xy])/(2*dx)-(pC[a_z]*pC[c_xxz])/(12*dx)+(pC[a_xz]*pC[c_xxz])/(24*dx)+(pC[a_y]*pC[c_xx])/(2*dx)-(pC[a_xy]*pC[c_xx])/(4*dx)-(pC[a_xx]*pC[c_xx])/(6*dx)+(pC[a_x]*pC[c_xx])/(2*dx)-(pC[a_0]*pC[c_xx])/dx-(pC[a_y]*pC[c_x])/(2*dx)+(pC[a_xy]*pC[c_x])/(4*dx)+(pC[a_xx]*pC[c_x])/(6*dx)-(pC[a_x]*pC[c_x])/(2*dx)+(pC[a_0]*pC[c_x])/dx;
}

template<typename REAL> inline
REAL JXBZ_100_101(
                  const REAL* const pC,
                  creal BGBX,
                  creal BGBY,
                  creal dx,
                  creal dy,
                  creal dz
) {
   using namespace Rec;
   return -(pC[b_z]*BGBY)/dz+(pC[b_yz]*BGBY)/(2*dz)-(pC[b_yyz]*BGBY)/(6*dz)-(pC[b_xz]*BGBY)/(2*dz)+(pC[b_xyz]*BGBY)/(4*dz)-(pC[c_yy]*BGBY)/dy+(pC[c_y]*BGBY)/dy+(pC[c_xy]*BGBY)/(2*dy)-(pC[a_z]*BGBX)/dz+(pC[a_yz]*BGBX)/(2*dz)-(pC[a_xz]*BGBX)/(2*dz)+(pC[a_xyz]*BGBX)/(4*dz)-(pC[a_xxz]*BGBX)/(6*dz)-(pC[c_xy]*BGBX)/(2*dx)+(pC[c_xx]*BGBX)/dx+(pC[c_x]*BGBX)/dx-(pC[b_z]*pC[b_zz])/(6*dz)+(pC[b_yz]*pC[b_zz])/(12*dz)+
     (pC[b_yzz]*pC[b_z])/(12*dz)-(pC[b_yy]*pC[b_z])/(6*dz)+(pC[b_y]*pC[b_z])/(2*dz)+(pC[b_xy]*pC[b_z])/(4*dz)-(pC[b_x]*pC[b_z])/(2*dz)-(pC[b_0]*pC[b_z])/dz-(pC[b_yz]*pC[b_yzz])/(24*dz)+(pC[b_yy]*pC[b_yz])/(12*dz)-(pC[b_y]*pC[b_yz])/(4*dz)-(pC[b_xy]*pC[b_yz])/(8*dz)+(pC[b_x]*pC[b_yz])/(4*dz)+(pC[b_0]*pC[b_yz])/(2*dz)-(pC[b_yy]*pC[b_yyz])/(36*dz)+(pC[b_y]*pC[b_yyz])/(12*dz)+(pC[b_xy]*pC[b_yyz])/(24*dz)-(pC[b_x]*pC[b_yyz])/(12*dz)-(pC[b_0]*pC[b_yyz])/(6*dz)-(pC[b_xz]*pC[b_yy])/(12*dz)+(pC[b_xyz]*pC[b_yy])/(24*dz)+(pC[b_xz]*pC[b_y])/(4*dz)-
     (pC[b_xyz]*pC[b_y])/(8*dz)+(pC[b_xy]*pC[b_xz])/(8*dz)-(pC[b_x]*pC[b_xz])/(4*dz)-(pC[b_0]*pC[b_xz])/(2*dz)-(pC[b_xy]*pC[b_xyz])/(16*dz)+(pC[b_x]*pC[b_xyz])/(8*dz)+(pC[b_0]*pC[b_xyz])/(4*dz)-(pC[a_z]*pC[a_zz])/(6*dz)-(pC[a_xz]*pC[a_zz])/(12*dz)+(pC[a_y]*pC[a_z])/(2*dz)-(pC[a_xzz]*pC[a_z])/(12*dz)+(pC[a_xy]*pC[a_z])/(4*dz)-(pC[a_xx]*pC[a_z])/(6*dz)-(pC[a_x]*pC[a_z])/(2*dz)-(pC[a_0]*pC[a_z])/dz-(pC[a_y]*pC[a_yz])/(4*dz)-(pC[a_xy]*pC[a_yz])/(8*dz)+(pC[a_xx]*pC[a_yz])/(12*dz)+(pC[a_x]*pC[a_yz])/(4*dz)+(pC[a_0]*pC[a_yz])/(2*dz)+(pC[a_xz]*pC[a_y])/(4*dz)
       -(pC[a_xyz]*pC[a_y])/(8*dz)+(pC[a_xxz]*pC[a_y])/(12*dz)-(pC[a_xz]*pC[a_xzz])/(24*dz)+(pC[a_xy]*pC[a_xz])/(8*dz)-(pC[a_xx]*pC[a_xz])/(12*dz)-(pC[a_x]*pC[a_xz])/(4*dz)-(pC[a_0]*pC[a_xz])/(2*dz)-(pC[a_xy]*pC[a_xyz])/(16*dz)+(pC[a_xx]*pC[a_xyz])/(24*dz)+(pC[a_x]*pC[a_xyz])/(8*dz)+(pC[a_0]*pC[a_xyz])/(4*dz)+(pC[a_xxz]*pC[a_xy])/(24*dz)-(pC[a_xx]*pC[a_xxz])/(36*dz)-(pC[a_x]*pC[a_xxz])/(12*dz)-(pC[a_0]*pC[a_xxz])/(6*dz)+(pC[b_z]*pC[c_yz])/(12*dy)-(pC[b_yz]*pC[c_yz])/(24*dy)-(pC[b_z]*pC[c_yyz])/(12*dy)+(pC[b_yz]*pC[c_yyz])/(24*dy)-
     (pC[b_yy]*pC[c_yy])/(6*dy)+(pC[b_y]*pC[c_yy])/(2*dy)+(pC[b_xy]*pC[c_yy])/(4*dy)-(pC[b_x]*pC[c_yy])/(2*dy)-(pC[b_0]*pC[c_yy])/dy+(pC[b_yy]*pC[c_y])/(6*dy)-(pC[b_y]*pC[c_y])/(2*dy)-(pC[b_xy]*pC[c_y])/(4*dy)+(pC[b_x]*pC[c_y])/(2*dy)+(pC[b_0]*pC[c_y])/dy+(pC[b_z]*pC[c_xyz])/(24*dy)-(pC[b_yz]*pC[c_xyz])/(48*dy)+(pC[b_yy]*pC[c_xy])/(12*dy)-(pC[b_y]*pC[c_xy])/(4*dy)-(pC[b_xy]*pC[c_xy])/(8*dy)+(pC[b_x]*pC[c_xy])/(4*dy)+(pC[b_0]*pC[c_xy])/(2*dy)+(pC[a_z]*pC[c_xz])/(12*dx)+(pC[a_xz]*pC[c_xz])/(24*dx)-(pC[a_z]*pC[c_xyz])/(24*dx)-
     (pC[a_xz]*pC[c_xyz])/(48*dx)+(pC[a_y]*pC[c_xy])/(4*dx)+(pC[a_xy]*pC[c_xy])/(8*dx)-(pC[a_xx]*pC[c_xy])/(12*dx)-(pC[a_x]*pC[c_xy])/(4*dx)-(pC[a_0]*pC[c_xy])/(2*dx)+(pC[a_z]*pC[c_xxz])/(12*dx)+(pC[a_xz]*pC[c_xxz])/(24*dx)-(pC[a_y]*pC[c_xx])/(2*dx)-(pC[a_xy]*pC[c_xx])/(4*dx)+(pC[a_xx]*pC[c_xx])/(6*dx)+(pC[a_x]*pC[c_xx])/(2*dx)+(pC[a_0]*pC[c_xx])/dx-(pC[a_y]*pC[c_x])/(2*dx)-(pC[a_xy]*pC[c_x])/(4*dx)+(pC[a_xx]*pC[c_x])/(6*dx)+(pC[a_x]*pC[c_x])/(2*dx)+(pC[a_0]*pC[c_x])/dx;
}

template<typename REAL> inline
REAL JXBZ_010_011(
                  const REAL* const pC,
                  creal BGBX,
                  creal BGBY,
                  creal dx,
                  creal dy,
                  creal dz
                 ) {
   using namespace Rec;
   return -(pC[b_z]*BGBY)/dz-(pC[b_yz]*BGBY)/(2*dz)-(pC[b_yyz]*BGBY)/(6*dz)+(pC[b_xz]*BGBY)/(2*dz)+(pC[b_xyz]*BGBY)/(4*dz)+(pC[c_yy]*BGBY)/dy+(pC[c_y]*BGBY)/dy-(pC[c_xy]*BGBY)/(2*dy)-(pC[a_z]*BGBX)/dz-(pC[a_yz]*BGBX)/(2*dz)+(pC[a_xz]*BGBX)/(2*dz)+(pC[a_xyz]*BGBX)/(4*dz)-(pC[a_xxz]*BGBX)/(6*dz)+(pC[c_xy]*BGBX)/(2*dx)-(pC[c_xx]*BGBX)/dx+(pC[c_x]*BGBX)/dx-(pC[b_z]*pC[b_zz])/(6*dz)-(pC[b_yz]*pC[b_zz])/(12*dz)-
     (pC[b_yzz]*pC[b_z])/(12*dz)-(pC[b_yy]*pC[b_z])/(6*dz)-(pC[b_y]*pC[b_z])/(2*dz)+(pC[b_xy]*pC[b_z])/(4*dz)+(pC[b_x]*pC[b_z])/(2*dz)-(pC[b_0]*pC[b_z])/dz-(pC[b_yz]*pC[b_yzz])/(24*dz)-(pC[b_yy]*pC[b_yz])/(12*dz)-(pC[b_y]*pC[b_yz])/(4*dz)+(pC[b_xy]*pC[b_yz])/(8*dz)+(pC[b_x]*pC[b_yz])/(4*dz)-(pC[b_0]*pC[b_yz])/(2*dz)-(pC[b_yy]*pC[b_yyz])/(36*dz)-(pC[b_y]*pC[b_yyz])/(12*dz)+(pC[b_xy]*pC[b_yyz])/(24*dz)+(pC[b_x]*pC[b_yyz])/(12*dz)-(pC[b_0]*pC[b_yyz])/(6*dz)+(pC[b_xz]*pC[b_yy])/(12*dz)+(pC[b_xyz]*pC[b_yy])/(24*dz)+(pC[b_xz]*pC[b_y])/(4*dz)+
     (pC[b_xyz]*pC[b_y])/(8*dz)-(pC[b_xy]*pC[b_xz])/(8*dz)-(pC[b_x]*pC[b_xz])/(4*dz)+(pC[b_0]*pC[b_xz])/(2*dz)-(pC[b_xy]*pC[b_xyz])/(16*dz)-(pC[b_x]*pC[b_xyz])/(8*dz)+(pC[b_0]*pC[b_xyz])/(4*dz)-(pC[a_z]*pC[a_zz])/(6*dz)+(pC[a_xz]*pC[a_zz])/(12*dz)-(pC[a_y]*pC[a_z])/(2*dz)+(pC[a_xzz]*pC[a_z])/(12*dz)+(pC[a_xy]*pC[a_z])/(4*dz)-(pC[a_xx]*pC[a_z])/(6*dz)+(pC[a_x]*pC[a_z])/(2*dz)-(pC[a_0]*pC[a_z])/dz-(pC[a_y]*pC[a_yz])/(4*dz)+(pC[a_xy]*pC[a_yz])/(8*dz)-(pC[a_xx]*pC[a_yz])/(12*dz)+(pC[a_x]*pC[a_yz])/(4*dz)-(pC[a_0]*pC[a_yz])/(2*dz)+(pC[a_xz]*pC[a_y])/(4*dz)
       +(pC[a_xyz]*pC[a_y])/(8*dz)-(pC[a_xxz]*pC[a_y])/(12*dz)-(pC[a_xz]*pC[a_xzz])/(24*dz)-(pC[a_xy]*pC[a_xz])/(8*dz)+(pC[a_xx]*pC[a_xz])/(12*dz)-(pC[a_x]*pC[a_xz])/(4*dz)+(pC[a_0]*pC[a_xz])/(2*dz)-(pC[a_xy]*pC[a_xyz])/(16*dz)+(pC[a_xx]*pC[a_xyz])/(24*dz)-(pC[a_x]*pC[a_xyz])/(8*dz)+(pC[a_0]*pC[a_xyz])/(4*dz)+(pC[a_xxz]*pC[a_xy])/(24*dz)-(pC[a_xx]*pC[a_xxz])/(36*dz)+(pC[a_x]*pC[a_xxz])/(12*dz)-(pC[a_0]*pC[a_xxz])/(6*dz)+(pC[b_z]*pC[c_yz])/(12*dy)+(pC[b_yz]*pC[c_yz])/(24*dy)+(pC[b_z]*pC[c_yyz])/(12*dy)+(pC[b_yz]*pC[c_yyz])/(24*dy)
         +(pC[b_yy]*pC[c_yy])/(6*dy)+(pC[b_y]*pC[c_yy])/(2*dy)-(pC[b_xy]*pC[c_yy])/(4*dy)-(pC[b_x]*pC[c_yy])/(2*dy)+(pC[b_0]*pC[c_yy])/dy+(pC[b_yy]*pC[c_y])/(6*dy)+(pC[b_y]*pC[c_y])/(2*dy)-(pC[b_xy]*pC[c_y])/(4*dy)-(pC[b_x]*pC[c_y])/(2*dy)+(pC[b_0]*pC[c_y])/dy-(pC[b_z]*pC[c_xyz])/(24*dy)-(pC[b_yz]*pC[c_xyz])/(48*dy)-(pC[b_yy]*pC[c_xy])/(12*dy)-(pC[b_y]*pC[c_xy])/(4*dy)+(pC[b_xy]*pC[c_xy])/(8*dy)+(pC[b_x]*pC[c_xy])/(4*dy)-(pC[b_0]*pC[c_xy])/(2*dy)+(pC[a_z]*pC[c_xz])/(12*dx)-(pC[a_xz]*pC[c_xz])/(24*dx)+(pC[a_z]*pC[c_xyz])/(24*dx)-
     (pC[a_xz]*pC[c_xyz])/(48*dx)+(pC[a_y]*pC[c_xy])/(4*dx)-(pC[a_xy]*pC[c_xy])/(8*dx)+(pC[a_xx]*pC[c_xy])/(12*dx)-(pC[a_x]*pC[c_xy])/(4*dx)+(pC[a_0]*pC[c_xy])/(2*dx)-(pC[a_z]*pC[c_xxz])/(12*dx)+(pC[a_xz]*pC[c_xxz])/(24*dx)-(pC[a_y]*pC[c_xx])/(2*dx)+(pC[a_xy]*pC[c_xx])/(4*dx)-(pC[a_xx]*pC[c_xx])/(6*dx)+(pC[a_x]*pC[c_xx])/(2*dx)-(pC[a_0]*pC[c_xx])/dx+(pC[a_y]*pC[c_x])/(2*dx)-(pC[a_xy]*pC[c_x])/(4*dx)+(pC[a_xx]*pC[c_x])/(6*dx)-(pC[a_x]*pC[c_x])/(2*dx)+(pC[a_0]*pC[c_x])/dx;
}

template<typename REAL> inline
REAL JXBZ_110_111(
                  const REAL* const pC,
                  creal BGBX,
                  creal BGBY,
                  creal dx,
                  creal dy,
                  creal dz
) {
   using namespace Rec;
   return -(pC[b_z]*BGBY)/dz-(pC[b_yz]*BGBY)/(2*dz)-(pC[b_yyz]*BGBY)/(6*dz)-(pC[b_xz]*BGBY)/(2*dz)-(pC[b_xyz]*BGBY)/(4*dz)+(pC[c_yy]*BGBY)/dy+(pC[c_y]*BGBY)/dy+(pC[c_xy]*BGBY)/(2*dy)-(pC[a_z]*BGBX)/dz-(pC[a_yz]*BGBX)/(2*dz)-(pC[a_xz]*BGBX)/(2*dz)-(pC[a_xyz]*BGBX)/(4*dz)-(pC[a_xxz]*BGBX)/(6*dz)+(pC[c_xy]*BGBX)/(2*dx)+(pC[c_xx]*BGBX)/dx+(pC[c_x]*BGBX)/dx-(pC[b_z]*pC[b_zz])/(6*dz)-(pC[b_yz]*pC[b_zz])/(12*dz)-
     (pC[b_yzz]*pC[b_z])/(12*dz)-(pC[b_yy]*pC[b_z])/(6*dz)-(pC[b_y]*pC[b_z])/(2*dz)-(pC[b_xy]*pC[b_z])/(4*dz)-(pC[b_x]*pC[b_z])/(2*dz)-(pC[b_0]*pC[b_z])/dz-(pC[b_yz]*pC[b_yzz])/(24*dz)-(pC[b_yy]*pC[b_yz])/(12*dz)-(pC[b_y]*pC[b_yz])/(4*dz)-(pC[b_xy]*pC[b_yz])/(8*dz)-(pC[b_x]*pC[b_yz])/(4*dz)-(pC[b_0]*pC[b_yz])/(2*dz)-(pC[b_yy]*pC[b_yyz])/(36*dz)-(pC[b_y]*pC[b_yyz])/(12*dz)-(pC[b_xy]*pC[b_yyz])/(24*dz)-(pC[b_x]*pC[b_yyz])/(12*dz)-(pC[b_0]*pC[b_yyz])/(6*dz)-(pC[b_xz]*pC[b_yy])/(12*dz)-(pC[b_xyz]*pC[b_yy])/(24*dz)-(pC[b_xz]*pC[b_y])/(4*dz)-
     (pC[b_xyz]*pC[b_y])/(8*dz)-(pC[b_xy]*pC[b_xz])/(8*dz)-(pC[b_x]*pC[b_xz])/(4*dz)-(pC[b_0]*pC[b_xz])/(2*dz)-(pC[b_xy]*pC[b_xyz])/(16*dz)-(pC[b_x]*pC[b_xyz])/(8*dz)-(pC[b_0]*pC[b_xyz])/(4*dz)-(pC[a_z]*pC[a_zz])/(6*dz)-(pC[a_xz]*pC[a_zz])/(12*dz)-(pC[a_y]*pC[a_z])/(2*dz)-(pC[a_xzz]*pC[a_z])/(12*dz)-(pC[a_xy]*pC[a_z])/(4*dz)-(pC[a_xx]*pC[a_z])/(6*dz)-(pC[a_x]*pC[a_z])/(2*dz)-(pC[a_0]*pC[a_z])/dz-(pC[a_y]*pC[a_yz])/(4*dz)-(pC[a_xy]*pC[a_yz])/(8*dz)-(pC[a_xx]*pC[a_yz])/(12*dz)-(pC[a_x]*pC[a_yz])/(4*dz)-(pC[a_0]*pC[a_yz])/(2*dz)-(pC[a_xz]*pC[a_y])/(4*dz)-
     (pC[a_xyz]*pC[a_y])/(8*dz)-(pC[a_xxz]*pC[a_y])/(12*dz)-(pC[a_xz]*pC[a_xzz])/(24*dz)-(pC[a_xy]*pC[a_xz])/(8*dz)-(pC[a_xx]*pC[a_xz])/(12*dz)-(pC[a_x]*pC[a_xz])/(4*dz)-(pC[a_0]*pC[a_xz])/(2*dz)-(pC[a_xy]*pC[a_xyz])/(16*dz)-(pC[a_xx]*pC[a_xyz])/(24*dz)-(pC[a_x]*pC[a_xyz])/(8*dz)-(pC[a_0]*pC[a_xyz])/(4*dz)-(pC[a_xxz]*pC[a_xy])/(24*dz)-(pC[a_xx]*pC[a_xxz])/(36*dz)-(pC[a_x]*pC[a_xxz])/(12*dz)-(pC[a_0]*pC[a_xxz])/(6*dz)+(pC[b_z]*pC[c_yz])/(12*dy)+(pC[b_yz]*pC[c_yz])/(24*dy)+(pC[b_z]*pC[c_yyz])/(12*dy)+(pC[b_yz]*pC[c_yyz])/(24*dy)+
     (pC[b_yy]*pC[c_yy])/(6*dy)+(pC[b_y]*pC[c_yy])/(2*dy)+(pC[b_xy]*pC[c_yy])/(4*dy)+(pC[b_x]*pC[c_yy])/(2*dy)+(pC[b_0]*pC[c_yy])/dy+(pC[b_yy]*pC[c_y])/(6*dy)+(pC[b_y]*pC[c_y])/(2*dy)+(pC[b_xy]*pC[c_y])/(4*dy)+(pC[b_x]*pC[c_y])/(2*dy)+(pC[b_0]*pC[c_y])/dy+(pC[b_z]*pC[c_xyz])/(24*dy)+(pC[b_yz]*pC[c_xyz])/(48*dy)+(pC[b_yy]*pC[c_xy])/(12*dy)+(pC[b_y]*pC[c_xy])/(4*dy)+(pC[b_xy]*pC[c_xy])/(8*dy)+(pC[b_x]*pC[c_xy])/(4*dy)+(pC[b_0]*pC[c_xy])/(2*dy)+(pC[a_z]*pC[c_xz])/(12*dx)+(pC[a_xz]*pC[c_xz])/(24*dx)+(pC[a_z]*pC[c_xyz])/(24*dx)+
     (pC[a_xz]*pC[c_xyz])/(48*dx)+(pC[a_y]*pC[c_xy])/(4*dx)+(pC[a_xy]*pC[c_xy])/(8*dx)+(pC[a_xx]*pC[c_xy])/(12*dx)+(pC[a_x]*pC[c_xy])/(4*dx)+(pC[a_0]*pC[c_xy])/(2*dx)+(pC[a_z]*pC[c_xxz])/(12*dx)+(pC[a_xz]*pC[c_xxz])/(24*dx)+(pC[a_y]*pC[c_xx])/(2*dx)+(pC[a_xy]*pC[c_xx])/(4*dx)+(pC[a_xx]*pC[c_xx])/(6*dx)+(pC[a_x]*pC[c_xx])/(2*dx)+(pC[a_0]*pC[c_xx])/dx+(pC[a_y]*pC[c_x])/(2*dx)+(pC[a_xy]*pC[c_x])/(4*dx)+(pC[a_xx]*pC[c_x])/(6*dx)+(pC[a_x]*pC[c_x])/(2*dx)+(pC[a_0]*pC[c_x])/dx;
}

void calculateEdgeHallTermXComponents(Real* cp,Real* derivs,const Real* const perturbedCoefficients,cint& RKCase) {
   #warning Particles (charge) assumed to be protons here

   Real hallRho;
   Real By, Bz, EXHall = 0.0;
   
   switch (Parameters::ohmHallTerm) {
    case 0:
      cerr << __FILE__ << __LINE__ << "You shouldn't be in a Hall term function if Parameters::ohmHallTerm == 0." << endl;
      break;
      
    case 1:
      if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
         By = cp[CellParams::PERBY]+cp[CellParams::BGBY];
         Bz = cp[CellParams::PERBZ]+cp[CellParams::BGBZ];
         hallRho =  (cp[CellParams::RHO] <= Parameters::hallMinimumRho ) ? Parameters::hallMinimumRho : cp[CellParams::RHO] ;
         EXHall = (Bz*((derivs[fieldsolver::dBGBxdz]+derivs[fieldsolver::dPERBxdz])/cp[CellParams::DZ] -
                       (derivs[fieldsolver::dBGBzdx]+derivs[fieldsolver::dPERBzdx])/cp[CellParams::DX]) -
                   By*((derivs[fieldsolver::dBGBydx]+derivs[fieldsolver::dPERBydx])/cp[CellParams::DX]-
                       ((derivs[fieldsolver::dBGBxdy]+derivs[fieldsolver::dPERBxdy])/cp[CellParams::DY])))
                / (physicalconstants::MU_0*hallRho*physicalconstants::CHARGE);
      }
      if (RKCase == RK_ORDER2_STEP1) {
         By = cp[CellParams::PERBY_DT2]+cp[CellParams::BGBY];
         Bz = cp[CellParams::PERBZ_DT2]+cp[CellParams::BGBZ];
         hallRho =  (cp[CellParams::RHO_DT2] <= Parameters::hallMinimumRho ) ? Parameters::hallMinimumRho : cp[CellParams::RHO_DT2] ;
         EXHall = (Bz*((derivs[fieldsolver::dBGBxdz]+derivs[fieldsolver::dPERBxdz])/cp[CellParams::DZ] -
                       (derivs[fieldsolver::dBGBzdx]+derivs[fieldsolver::dPERBzdx])/cp[CellParams::DX]) -
                   By*((derivs[fieldsolver::dBGBydx]+derivs[fieldsolver::dPERBydx])/cp[CellParams::DX]-
                       ((derivs[fieldsolver::dBGBxdy]+derivs[fieldsolver::dPERBxdy])/cp[CellParams::DY])))
                / (physicalconstants::MU_0*hallRho*physicalconstants::CHARGE);
      }

      cp[CellParams::EXHALL_000_100] =
      cp[CellParams::EXHALL_010_110] =
      cp[CellParams::EXHALL_001_101] =
      cp[CellParams::EXHALL_011_111] = EXHall;
      break;
    case 2:
      if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
         hallRho =  (cp[CellParams::RHO] <= Parameters::hallMinimumRho ) ? Parameters::hallMinimumRho : cp[CellParams::RHO] ;
         cp[CellParams::EXHALL_000_100] = JXBX_000_100(perturbedCoefficients, cp[CellParams::BGBY], cp[CellParams::BGBZ], 
                                                       cp[CellParams::DX], cp[CellParams::DY], cp[CellParams::DZ]) 
                                        / (physicalconstants::MU_0*hallRho*physicalconstants::CHARGE);
         cp[CellParams::EXHALL_010_110] = JXBX_010_110(perturbedCoefficients, cp[CellParams::BGBY], cp[CellParams::BGBZ], 
                                                       cp[CellParams::DX], cp[CellParams::DY], cp[CellParams::DZ]) 
                                        / (physicalconstants::MU_0*hallRho*physicalconstants::CHARGE);
         cp[CellParams::EXHALL_001_101] = JXBX_001_101(perturbedCoefficients, cp[CellParams::BGBY], cp[CellParams::BGBZ], 
                                                       cp[CellParams::DX], cp[CellParams::DY], cp[CellParams::DZ]) 
                                        / (physicalconstants::MU_0*hallRho*physicalconstants::CHARGE);
         cp[CellParams::EXHALL_011_111] = JXBX_011_111(perturbedCoefficients, cp[CellParams::BGBY], cp[CellParams::BGBZ], 
                                                       cp[CellParams::DX], cp[CellParams::DY], cp[CellParams::DZ]) 
                                        / (physicalconstants::MU_0*hallRho*physicalconstants::CHARGE);
      }
      if (RKCase == RK_ORDER2_STEP1) {
         hallRho =  (cp[CellParams::RHO_DT2] <= Parameters::hallMinimumRho ) ? Parameters::hallMinimumRho : cp[CellParams::RHO_DT2] ;
         cp[CellParams::EXHALL_000_100] = JXBX_000_100(perturbedCoefficients, cp[CellParams::BGBY], cp[CellParams::BGBZ], 
                                                       cp[CellParams::DX], cp[CellParams::DY], cp[CellParams::DZ]) 
                                        / (physicalconstants::MU_0*hallRho*physicalconstants::CHARGE);
         cp[CellParams::EXHALL_010_110] = JXBX_010_110(perturbedCoefficients, cp[CellParams::BGBY], cp[CellParams::BGBZ], 
                                                       cp[CellParams::DX], cp[CellParams::DY], cp[CellParams::DZ]) 
                                        / (physicalconstants::MU_0*hallRho*physicalconstants::CHARGE);
         cp[CellParams::EXHALL_001_101] = JXBX_001_101(perturbedCoefficients, cp[CellParams::BGBY], cp[CellParams::BGBZ], 
                                                       cp[CellParams::DX], cp[CellParams::DY], cp[CellParams::DZ]) 
                                        / (physicalconstants::MU_0*hallRho*physicalconstants::CHARGE);
         cp[CellParams::EXHALL_011_111] = JXBX_011_111(perturbedCoefficients, cp[CellParams::BGBY], cp[CellParams::BGBZ], 
                                                       cp[CellParams::DX], cp[CellParams::DY], cp[CellParams::DZ]) 
                                        / (physicalconstants::MU_0*hallRho*physicalconstants::CHARGE);
      }
      break;

    default:
      cerr << __FILE__ << ":" << __LINE__ << "You are welcome to code higher-order Hall term correction terms." << endl;
      break;
   }
}

void calculateEdgeHallTermYComponents(Real* cp,Real* derivs,const Real* const perturbedCoefficients,cint& RKCase) {
   #warning Particles (charge) assumed to be protons here

   Real hallRho;
   Real Bx, Bz, EYHall = 0.0;
   
   switch (Parameters::ohmHallTerm) {
    case 0:
      cerr << __FILE__ << __LINE__ << "You shouldn't be in a Hall term function if Parameters::ohmHallTerm == 0." << endl;
      break;
      
    case 1:
      if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
         Bx = cp[CellParams::PERBX]+cp[CellParams::BGBX];
         Bz = cp[CellParams::PERBZ]+cp[CellParams::BGBZ];
         hallRho =  (cp[CellParams::RHO] <= Parameters::hallMinimumRho ) ? Parameters::hallMinimumRho : cp[CellParams::RHO] ;
         EYHall = (Bx*((derivs[fieldsolver::dBGBydx]+derivs[fieldsolver::dPERBydx])/cp[CellParams::DX] -
                       (derivs[fieldsolver::dBGBxdy]+derivs[fieldsolver::dPERBxdy])/cp[CellParams::DY]) -
                   Bz*((derivs[fieldsolver::dBGBzdy]+derivs[fieldsolver::dPERBzdy])/cp[CellParams::DY] -
                       ((derivs[fieldsolver::dBGBydz]+derivs[fieldsolver::dPERBydz])/cp[CellParams::DZ] )))
                / (physicalconstants::MU_0*hallRho*physicalconstants::CHARGE);
      }
      if (RKCase == RK_ORDER2_STEP1) {
         Bx = cp[CellParams::PERBX_DT2]+cp[CellParams::BGBX];
         Bz = cp[CellParams::PERBZ_DT2]+cp[CellParams::BGBZ];
         hallRho =  (cp[CellParams::RHO_DT2] <= Parameters::hallMinimumRho ) ? Parameters::hallMinimumRho : cp[CellParams::RHO_DT2] ;
         EYHall = (Bx*((derivs[fieldsolver::dBGBydx]+derivs[fieldsolver::dPERBydx])/cp[CellParams::DX] -
                       (derivs[fieldsolver::dBGBxdy]+derivs[fieldsolver::dPERBxdy])/cp[CellParams::DY]) -
                   Bz*((derivs[fieldsolver::dBGBzdy]+derivs[fieldsolver::dPERBzdy])/cp[CellParams::DY] -
                       ((derivs[fieldsolver::dBGBydz]+derivs[fieldsolver::dPERBydz])/cp[CellParams::DZ] )))
                / (physicalconstants::MU_0*hallRho*physicalconstants::CHARGE);
      }

      cp[CellParams::EYHALL_000_010] =
      cp[CellParams::EYHALL_100_110] =
      cp[CellParams::EYHALL_101_111] =
      cp[CellParams::EYHALL_001_011] = EYHall;
      break;
      
    case 2:
      if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
         hallRho =  (cp[CellParams::RHO] <= Parameters::hallMinimumRho ) ? Parameters::hallMinimumRho : cp[CellParams::RHO] ;
         cp[CellParams::EYHALL_000_010] = JXBY_000_010(perturbedCoefficients, cp[CellParams::BGBX], cp[CellParams::BGBZ], 
						       cp[CellParams::DX], cp[CellParams::DY], cp[CellParams::DZ]) 
	                                / (physicalconstants::MU_0*hallRho*physicalconstants::CHARGE);
         cp[CellParams::EYHALL_100_110] = JXBY_100_110(perturbedCoefficients, cp[CellParams::BGBX], cp[CellParams::BGBZ], 
						       cp[CellParams::DX], cp[CellParams::DY], cp[CellParams::DZ]) 
	                                / (physicalconstants::MU_0*hallRho*physicalconstants::CHARGE);
         cp[CellParams::EYHALL_001_011] = JXBY_001_011(perturbedCoefficients, cp[CellParams::BGBX], cp[CellParams::BGBZ], 
						       cp[CellParams::DX], cp[CellParams::DY], cp[CellParams::DZ]) 
	   / (physicalconstants::MU_0*hallRho*physicalconstants::CHARGE);
         cp[CellParams::EYHALL_101_111] = JXBY_101_111(perturbedCoefficients, cp[CellParams::BGBX], cp[CellParams::BGBZ], 
						       cp[CellParams::DX], cp[CellParams::DY], cp[CellParams::DZ]) 
	   / (physicalconstants::MU_0*hallRho*physicalconstants::CHARGE);
      }
      if (RKCase == RK_ORDER2_STEP1) {
         hallRho =  (cp[CellParams::RHO_DT2] <= Parameters::hallMinimumRho ) ? Parameters::hallMinimumRho : cp[CellParams::RHO_DT2] ;
         cp[CellParams::EYHALL_000_010] = JXBY_000_010(perturbedCoefficients, cp[CellParams::BGBX], cp[CellParams::BGBZ], 
						       cp[CellParams::DX], cp[CellParams::DY], cp[CellParams::DZ]) 
	   / (physicalconstants::MU_0*hallRho*physicalconstants::CHARGE);
         cp[CellParams::EYHALL_100_110] = JXBY_100_110(perturbedCoefficients, cp[CellParams::BGBX], cp[CellParams::BGBZ], 
						       cp[CellParams::DX], cp[CellParams::DY], cp[CellParams::DZ]) 
	   / (physicalconstants::MU_0*hallRho*physicalconstants::CHARGE);
         cp[CellParams::EYHALL_001_011] = JXBY_001_011(perturbedCoefficients, cp[CellParams::BGBX], cp[CellParams::BGBZ], 
						       cp[CellParams::DX], cp[CellParams::DY], cp[CellParams::DZ]) 
	   / (physicalconstants::MU_0*hallRho*physicalconstants::CHARGE);
         cp[CellParams::EYHALL_101_111] = JXBY_101_111(perturbedCoefficients, cp[CellParams::BGBX], cp[CellParams::BGBZ], 
						       cp[CellParams::DX], cp[CellParams::DY], cp[CellParams::DZ]) 
	   / (physicalconstants::MU_0*hallRho*physicalconstants::CHARGE);
      }
      break;

    default:
      cerr << __FILE__ << ":" << __LINE__ << "You are welcome to code higher-order Hall term correction terms." << endl;
      break;
   }
}

void calculateEdgeHallTermZComponents(Real* cp,Real* derivs,const Real* const perturbedCoefficients,cint& RKCase) {
  #warning Particles (charge) assumed to be protons here

   Real hallRho;
   Real Bx, By, EZHall=0.0;

   switch (Parameters::ohmHallTerm) {
   case 0:
     cerr << __FILE__ << __LINE__ << "You shouldn't be in a Hall term function if Parameters::ohmHallTerm == 0." << endl;
     break;

   case 1:
     if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
        Bx = cp[CellParams::PERBX]+cp[CellParams::BGBX];
        By = cp[CellParams::PERBY]+cp[CellParams::BGBY];
        hallRho =  (cp[CellParams::RHO] <= Parameters::hallMinimumRho ) ? Parameters::hallMinimumRho : cp[CellParams::RHO] ;
        EZHall = (By*((derivs[fieldsolver::dBGBzdy]+derivs[fieldsolver::dPERBzdy])/cp[CellParams::DY] -
                      (derivs[fieldsolver::dBGBydz]+derivs[fieldsolver::dPERBydz])/cp[CellParams::DZ]) -
                  Bx*((derivs[fieldsolver::dBGBxdz]+derivs[fieldsolver::dPERBxdz])/cp[CellParams::DZ] -
                      ((derivs[fieldsolver::dBGBzdx]+derivs[fieldsolver::dPERBzdx])/cp[CellParams::DX])))
          / (physicalconstants::MU_0*hallRho*physicalconstants::CHARGE);
     }
     if (RKCase == RK_ORDER2_STEP1) {
        Bx = cp[CellParams::PERBX_DT2]+cp[CellParams::BGBX];
        By = cp[CellParams::PERBY_DT2]+cp[CellParams::BGBY];
        hallRho =  (cp[CellParams::RHO_DT2] <= Parameters::hallMinimumRho ) ? Parameters::hallMinimumRho : cp[CellParams::RHO_DT2] ;
        EZHall = (By*((derivs[fieldsolver::dBGBzdy]+derivs[fieldsolver::dPERBzdy])/cp[CellParams::DY] -
                      (derivs[fieldsolver::dBGBydz]+derivs[fieldsolver::dPERBydz])/cp[CellParams::DZ]) -
                  Bx*((derivs[fieldsolver::dBGBxdz]+derivs[fieldsolver::dPERBxdz])/cp[CellParams::DZ] -
                      ((derivs[fieldsolver::dBGBzdx]+derivs[fieldsolver::dPERBzdx])/cp[CellParams::DX])))
          / (physicalconstants::MU_0*hallRho*physicalconstants::CHARGE);
     }

      cp[CellParams::EZHALL_000_001] =
      cp[CellParams::EZHALL_100_101] =
      cp[CellParams::EZHALL_110_111] =
      cp[CellParams::EZHALL_010_011] = EZHall;
     break;

   case 2:
      if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
         hallRho =  (cp[CellParams::RHO] <= Parameters::hallMinimumRho ) ? Parameters::hallMinimumRho : cp[CellParams::RHO] ;
         cp[CellParams::EZHALL_000_001] = JXBZ_000_001(perturbedCoefficients, cp[CellParams::BGBX], cp[CellParams::BGBY], 
                                                       cp[CellParams::DX], cp[CellParams::DY], cp[CellParams::DZ]) 
           / (physicalconstants::MU_0*hallRho*physicalconstants::CHARGE);
         cp[CellParams::EZHALL_100_101] = JXBZ_100_101(perturbedCoefficients, cp[CellParams::BGBX], cp[CellParams::BGBY], 
                                                       cp[CellParams::DX], cp[CellParams::DY], cp[CellParams::DZ]) 
           / (physicalconstants::MU_0*hallRho*physicalconstants::CHARGE);
         cp[CellParams::EZHALL_010_011] = JXBZ_010_011(perturbedCoefficients, cp[CellParams::BGBX], cp[CellParams::BGBY], 
                                                       cp[CellParams::DX], cp[CellParams::DY], cp[CellParams::DZ]) 
           / (physicalconstants::MU_0*hallRho*physicalconstants::CHARGE);
         cp[CellParams::EZHALL_110_111] = JXBZ_110_111(perturbedCoefficients, cp[CellParams::BGBX], cp[CellParams::BGBY], 
                                                       cp[CellParams::DX], cp[CellParams::DY], cp[CellParams::DZ]) 
           / (physicalconstants::MU_0*hallRho*physicalconstants::CHARGE);
      }
      if (RKCase == RK_ORDER2_STEP1) {
         hallRho =  (cp[CellParams::RHO_DT2] <= Parameters::hallMinimumRho ) ? Parameters::hallMinimumRho : cp[CellParams::RHO_DT2] ;
         cp[CellParams::EZHALL_000_001] = JXBZ_000_001(perturbedCoefficients, cp[CellParams::BGBX], cp[CellParams::BGBY], 
                                                       cp[CellParams::DX], cp[CellParams::DY], cp[CellParams::DZ]) 
           / (physicalconstants::MU_0*hallRho*physicalconstants::CHARGE);
         cp[CellParams::EZHALL_100_101] = JXBZ_100_101(perturbedCoefficients, cp[CellParams::BGBX], cp[CellParams::BGBY], 
                                                       cp[CellParams::DX], cp[CellParams::DY], cp[CellParams::DZ]) 
           / (physicalconstants::MU_0*hallRho*physicalconstants::CHARGE);
         cp[CellParams::EZHALL_010_011] = JXBZ_010_011(perturbedCoefficients, cp[CellParams::BGBX], cp[CellParams::BGBY], 
                                                       cp[CellParams::DX], cp[CellParams::DY], cp[CellParams::DZ]) 
           / (physicalconstants::MU_0*hallRho*physicalconstants::CHARGE);
         cp[CellParams::EZHALL_110_111] = JXBZ_110_111(perturbedCoefficients, cp[CellParams::BGBX], cp[CellParams::BGBY], 
                                                       cp[CellParams::DX], cp[CellParams::DY], cp[CellParams::DZ]) 
           / (physicalconstants::MU_0*hallRho*physicalconstants::CHARGE);
      }
      break;
      
    default:
      cerr << __FILE__ << ":" << __LINE__ << "You are welcome to code higher-order Hall term correction terms." << endl;
      break;
   }
}

/** Calculate the Hall term on all given cells.
 * @param sysBoundaries System boundary condition functions.
 * @param cache Cache for local cells.
 * @param cells Local IDs of calculated cells, one of the vectors in fs_cache::CacheContainer.
 * @param RKCase
 */
void calculateHallTerm(SysBoundary& sysBoundaries,std::vector<fs_cache::CellCache>& cache,
                       const std::vector<uint16_t>& cells,cint& RKCase) {

   #pragma omp parallel for
   for (size_t c=0; c<cells.size(); ++c) { // DO_NOT_COMPUTE cells already removed
      const uint16_t localID = cells[c];

      #ifdef DEBUG_FSOLVER
      if (localID >= cache.size()) {
         cerr << "local index out of bounds in " << __FILE__ << ":" << __LINE__ << endl;
         exit(1);
      }
      if (cache[localID].cells[fs_cache::calculateNbrID(1,1,1)] == NULL) {
         cerr << "NULL pointer in " << __FILE__ << ":" << __LINE__ << endl;
         exit(1);
      }
      if (cache[localID].cells[fs_cache::calculateNbrID(1,1,1)]->parameters == NULL) {
         cerr << "NULL cell parameters in " << __FILE__ << ":" << __LINE__ << endl;
         exit(1);
      }
      if (cache[localID].cells[fs_cache::calculateNbrID(1,1,1)]->derivatives == NULL) {
         cerr << "NULL derivatives in " << __FILE__ << ":" << __LINE__ << endl;
         exit(1);
      }
      #endif

      cuint fieldSolverSysBoundaryFlag = cache[localID].existingCellsFlags;
      cuint cellSysBoundaryFlag        = cache[localID].sysBoundaryFlag;
      cuint cellSysBoundaryLayer       = cache[localID].cells[fs_cache::calculateNbrID(1,1,1)]->sysBoundaryLayer;

      Real perturbedCoefficients[Rec::N_REC_COEFFICIENTS];

      reconstructionCoefficients(cache[localID],
                                 perturbedCoefficients,
                                 3, // Reconstruction order of the fields after Balsara 2009, 2 used for general B, 3 used here for 2nd-order Hall term
                                 RKCase
                                );

      if ((fieldSolverSysBoundaryFlag & CALCULATE_EX) == CALCULATE_EX) {
         if ((cellSysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY) &&
             (cellSysBoundaryLayer != 1)) {
            sysBoundaries.getSysBoundary(cellSysBoundaryFlag)->fieldSolverBoundaryCondHallElectricField(cache[localID],RKCase,0);
         } else {
            Real* cp     = cache[localID].cells[fs_cache::calculateNbrID(1,1,1)]->parameters;
            Real* derivs = cache[localID].cells[fs_cache::calculateNbrID(1,1,1)]->derivatives;
            calculateEdgeHallTermXComponents(cp,derivs,perturbedCoefficients,RKCase);
         }
      }
      if ((fieldSolverSysBoundaryFlag & CALCULATE_EY) == CALCULATE_EY) {
         if ((cellSysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY) &&
             (cellSysBoundaryLayer != 1)) {
            sysBoundaries.getSysBoundary(cellSysBoundaryFlag)->fieldSolverBoundaryCondHallElectricField(cache[localID],RKCase,1);
         } else {
            Real* cp     = cache[localID].cells[fs_cache::calculateNbrID(1,1,1)]->parameters;
            Real* derivs = cache[localID].cells[fs_cache::calculateNbrID(1,1,1)]->derivatives;
            calculateEdgeHallTermYComponents(cp,derivs,perturbedCoefficients,RKCase);
         }
      }
      if ((fieldSolverSysBoundaryFlag & CALCULATE_EZ) == CALCULATE_EZ) {
         if ((cellSysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY) &&
             (cellSysBoundaryLayer != 1)) {
            sysBoundaries.getSysBoundary(cellSysBoundaryFlag)->fieldSolverBoundaryCondHallElectricField(cache[localID],RKCase,2);
         } else {
            Real* cp     = cache[localID].cells[fs_cache::calculateNbrID(1,1,1)]->parameters;
            Real* derivs = cache[localID].cells[fs_cache::calculateNbrID(1,1,1)]->derivatives;
            calculateEdgeHallTermZComponents(cp,derivs,perturbedCoefficients,RKCase);
         }
      }
   }
}

void calculateHallTermSimple(
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   SysBoundary& sysBoundaries,
   const vector<CellID>& localCells,
   cint& RKCase,
   const bool communicateDerivatives
) {

   namespace fs = fieldsolver;
   int timer;
   size_t N_cells;
   phiprof::start("Calculate Hall term");
   if(communicateDerivatives) {
      SpatialCell::set_mpi_transfer_type(Transfer::CELL_DERIVATIVES);

      fs_cache::CacheContainer& cacheContainer = fs_cache::getCache();

      timer=phiprof::initializeTimer("Start communication of derivatives","MPI");
      phiprof::start(timer);
      mpiGrid.start_remote_neighbor_copy_updates(FIELD_SOLVER_NEIGHBORHOOD_ID);
      phiprof::stop(timer);

      // Calculate Hall term on inner cells
      timer=phiprof::initializeTimer("Compute inner cells");
      phiprof::start(timer);
      calculateHallTerm(sysBoundaries,cacheContainer.localCellsCache,cacheContainer.cellsWithLocalNeighbours,RKCase);
      phiprof::stop(timer,cacheContainer.cellsWithLocalNeighbours.size(),"Spatial Cells");

      timer=phiprof::initializeTimer("Wait for receives","MPI","Wait");
      phiprof::start(timer);
      mpiGrid.wait_remote_neighbor_copy_update_receives(FIELD_SOLVER_NEIGHBORHOOD_ID);
      phiprof::stop(timer);
      
      // Calculate Hall term on boundary cells:
      timer=phiprof::initializeTimer("Compute boundary cells");
      phiprof::start(timer);
      calculateHallTerm(sysBoundaries,cacheContainer.localCellsCache,cacheContainer.cellsWithRemoteNeighbours,RKCase);
      phiprof::stop(timer,cacheContainer.cellsWithRemoteNeighbours.size(),"Spatial Cells");

      timer=phiprof::initializeTimer("Wait for sends","MPI","Wait");
      phiprof::start(timer);
      mpiGrid.wait_remote_neighbor_copy_update_sends();
      phiprof::stop(timer);

      N_cells = cacheContainer.cellsWithRemoteNeighbours.size()
      + cacheContainer.cellsWithLocalNeighbours.size();
   } else {
      fs_cache::CacheContainer& cacheContainer = fs_cache::getCache();
      calculateHallTerm(sysBoundaries,cacheContainer.localCellsCache,cacheContainer.local_NOT_DO_NOT_COMPUTE,RKCase);
      N_cells = cacheContainer.local_NOT_DO_NOT_COMPUTE.size();
   }

   phiprof::stop("Calculate Hall term",N_cells,"Spatial Cells");
}
