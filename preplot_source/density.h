#ifndef __DENSITY_H__
#define __DENSITY_H__

#include "include.h"
#include "helper.h"
#include "io.h"
#include <omp.h>
#include <algorithm>

extern int id;
extern float slen;
extern float exrad;

int doDensityNormalization (double ** &data, int num_lines, int num_fields);

int doGaussianSmoothDensity3D (double ** & halo_data, int num_lines, int num_fields, int x, int y);

int doGaussianSmoothDensity (double ** & halo_data, int num_lines, int num_fields, int x, int y);

int doBuildDensityCatalog (double ** & halo_data, int num_lines, int num_fields, int x, int y);

int doCalcVoxelDensityDistributions (int x, int y, int num_bins, bool NORM);

int doCombineHlistFiles (double ** & halo_data, int num_lines, int num_fields);

int doDensityPercentileEvolution (double ** &data, int num_lines, int num_fields, int x);

#endif
