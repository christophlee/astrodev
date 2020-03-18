#ifndef __FILEOP_H__
#define __FILEOP_H__

#include "include.h"
#include "helper.h"
#include "io.h"
#include "stats.h"
#include <algorithm>

#define THREEQ_BINNING 1
#define PERCENT_BINNING 2

extern double doBin_data_min;
extern double doBin_data_max;
extern bool binnedHist2D;
extern bool binnedMedian;
extern bool binnedHist1D;

int doBin(double ** &data, int num_lines, int num_fields, int x, int num_bins, threshold_object &thresholds, int threshold_type, bool binary_write);

int doBinnedAnalysis (double ** &data, int num_lines, int num_fields, int binx, int num_bins, threshold_object &thresholds, int threshold_type, int x, int y, int num_bins_x, int num_bins_y, double bin_width_x, double bin_width_y, bool * LOG, bool NORM, bool SMOOTH, bool binary_write, bool autobin, int auto_tnum, bool use_percentiles_x, bool use_percentiles_y, bool RANK, bool CDF, int use_rw, bool MASSFUNC, int POINTS_PER_CELL, bool SCATTER);

int doAddFieldToCatalog (double ** &data, int num_lines, int num_fields, bool binary_write);

int doStackCatalogs ();

int doReweightMassDistribution (double ** &data, int num_lines, int num_fields, int x);

int doReorderFields (double ** &data, int num_lines, int num_fields);
int doReorderFields2 (double ** &data, int num_lines, int num_fields);

int doCatalogMerge (double ** &data, int num_lines, int num_fields, std::string mergeTarget);
int doCatalogMerge (double ** &data, int num_lines, int num_fields, double ** &data2, int num_lines2, int num_fields2);

#endif
