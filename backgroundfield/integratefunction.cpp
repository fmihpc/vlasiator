/*
Background magnetic field class of Vlasiator.

Copyright 1997, 1998, 1999, 2000, 2001, 2010, 2011, 2012 Finnish Meteorological Institute
*/

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>

#include "../common.h"
#include "integratefunction.hpp"
#include "functions.hpp"
#include "quadr.hpp"


double LineAverage(const T3DFunction& f1, coordinate line, double accuracy,
		 const double r1[3], double L)
{
  //subroutines not threadsafe 
#pragma omp critical
  {
    const double norm = 1/L;
    const double acc = accuracy*L;
    const double a = r1[line];
    const double b = r1[line] + L;
    
    switch (line) {
    case X: 
      {
	T3D_fix23 f(f1,r1[1],r1[2]); 
	return Romberg(f,a,b,acc)*norm;
      }
      break;
    case Y: 
      {
	T3D_fix13 f(f1,r1[0],r1[2]); 
	return Romberg(f,a,b,acc)*norm;
      }
      break;
    case Z: 
      {
	T3D_fix12 f(f1,r1[0],r1[1]); 
	return Romberg(f,a,b,acc)*norm;
      }
      break;
    default:
      cerr << "*** lineAverage  is bad\n";
      exit(1);
      break;
    }
  }
}


double surfaceAverage(const T3DFunction& f1, 
		    coordinate face, double accuracy,
		    const double r1[3],
		    double L1,
		    double L2) {
#pragma omp critical
  {
    const double acc = accuracy*L1*L2;
    const double norm = 1/(L1*L2);
    double value;
    switch (face) {
    case X:
      {
	T3D_fix1 f(f1,r1[0]);
	value = Romberg(f, r1[1],r1[1]+L1, r1[2],r1[2]+L2, acc)*norm;
      }
      break;
    case Y:
      {
	T3D_fix2 f(f1,r1[1]);
	value = Romberg(f, r1[0],r1[0]+L1, r1[2],r1[2]+L2, acc)*norm; 
      }
      break;
    case Z:
      {
	T3D_fix3 f(f1,r1[2]);
	value = Romberg(f, r1[0],r1[0]+L1, r1[1],r1[1]+L2, acc)*norm;
      }
      break;
    default:
      cerr << "*** SurfaceAverage  is bad\n";
      exit(1);
    }
    return value;
  }
}


double volumeAverage(const T3DFunction& f1, double accuracy,
		   const double r1[3], const double r2[3])
{
#pragma omp critical
  {
    const double acc = accuracy*(r2[0]-r1[0])*(r2[1]-r1[1])*(r2[2]-r1[2]);
    const double norm = 1.0/((r2[0]-r1[0])*(r2[1]-r1[1])*(r2[2]-r1[2]));
    return  Romberg(f1, r1[0],r2[0], r1[1],r2[1], r1[2],r2[2], acc)*norm;
  }
}

