/*
This file is part of Vlasiator.

Copyright 2010, 2011, 2012, 2013, 2014 Finnish Meteorological Institute

*/

#include "fs_limiters.h"

/*
Real limiter(creal& left,creal& cent,creal& rght) {
   //const Real limited = minmod(left,cent,rght);
   //const Real limited = MClimiter(left,cent,rght);
   const Real limited = vanLeer(left,cent,rght);
   
   #ifdef DEBUG_SOLVERS
   if (limited != limited
      || limited * 0 != 0) {
      std::cerr << __FILE__ << ":" << __LINE__
                  << " Limiter returned an invalid value " << limited
                  << " with left, center, right: " << left << ", " << cent << ", " << rght
                  << std::endl;
      abort();
   }
   #endif
   
   return limited;
}
*/

