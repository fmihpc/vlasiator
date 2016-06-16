
#ifndef NOISE_H
#define NOISE_H

#pragma once
// #include "simplexnoise1234.h"
#include "vectorclass.h"
#include "vector3d.h"
#include <vector>
#include <math.h>

Vec3d divergence_free_noise(Vec3d p, Vec3d scales, int num_octaves);

//static unsigned char perm[512];
//static double grad2lut[8][2];
//static double grad3lut[16][3];
//static double grad4lut[32][4];
//static unsigned char simplex[64][4];

//#pragma omp threadprivate(perm,grad2lut,grad3lut,grad4lut,simplex)

double sdnoise1( double x, double *dnoise_dx);
double sdnoise2( double x, double y, double *dnoise_dx, double *dnoise_dy );
double sdnoise3( double x, double y, double z,
                double *dnoise_dx, double *dnoise_dy, double *dnoise_dz );
double sdnoise4( double x, double y, double z, double w,
                double *dnoise_dx, double *dnoise_dy,
                double *dnoise_dz, double *dnoise_dw);

#endif
