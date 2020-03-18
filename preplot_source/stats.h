#ifndef __STATS_H__
#define __STATS_H__

#include "include.h"
#include "helper.h"
#include "io.h"
#include <omp.h>
#include <limits>
#include <random>
#include <algorithm>
#include <functional>

//extern std::string output;
extern double doBin_data_min;
extern double doBin_data_max;

struct spline{
    double a;
    double b;
    double c;
    double d;
    double x;
};

double pdf_integrate (double lbound, double rbound, int x);

int doAvg(double ** &data, int num_lines, int x, int y, int num_bins);

int doMedian(double ** &data, int num_lines, int x, int y, int num_bins, bool autobin, int auto_tnum, bool use_percentiles_x = false, bool use_percentiles_y = false);

int doHist1D(double ** &data, int num_lines, int x, int num_bins, bool NORM, bool LOGZ, bool CDF, bool autobin, int auto_tnum, bool MASSFUNC, double bin_width);

int doHist2D(double ** &data, int num_lines, int x, int y, int num_bins_x, int num_bins_y, double bin_width_x, double bin_width_y, bool LOGZ, bool NORM, bool SMOOTH, bool RANK, int POINTS_PER_CELL, bool SCATTER);

int doRankCorrelation (double ** &data, int num_lines, int num_fields, int x, int y);

int doGenPlotStats (double ** &data, int num_lines, int num_fields, int x, int y);

std::vector<spline> cubicSpline (std::vector<double> &x, std::vector<double> &y);


/*****************
 * vavg
 * ---------------------------------------------------------------------------------------
 * Finds and returns the average value of vector passed in.
 * ---------------------------------------------------------------------------------------
 * v        sorted vector of data (flexible type)
 * returns  average value
 *****************/
template <typename T>
T vavg (std::vector<T> &v) {

    int len = v.size();
    T avg;
    for (int i = 0; i < len; i++) avg += v[i];

    return avg/len;
}

/*****************
 * vmedian
 * ---------------------------------------------------------------------------------------
 * Finds and returns the median value of vector passed in.
 * ---------------------------------------------------------------------------------------
 * v        sorted vector of data (flexible type)
 * returns  median value
 *****************/
template <typename T>
T vmedian (std::vector<T> &v) {

    int len = v.size();

    if (len % 2) return v[(len-1)/2];
    else return 0.5*(v[(len-2)/2] + v[(len)/2]);
}

/*****************
 * vCIlo
 * ---------------------------------------------------------------------------------------
 * Finds and returns the low bound of 95% CI on median of vector passed in.
 * ---------------------------------------------------------------------------------------
 * v        sorted vector of data (flexible type)
 * returns  low bound of 95% CI
 *****************/
template <typename T>
T vCIlo (std::vector<T> &v) {

    int len = v.size();
    double lo = len*.5-1.96*sqrt(len*.25)+1;
    return v[(int)floor(lo)];
}

/*****************
 * vCIhi
 * ---------------------------------------------------------------------------------------
 * Finds and returns the high bound of 95% CI on median of vector passed in.
 * ---------------------------------------------------------------------------------------
 * v        sorted vector of data (flexible type)
 * returns  high bound of 95% CI
 *****************/
template <typename T>
T vCIhi (std::vector<T> &v) {

    int len = v.size();
    double hi = len*.5+1.96*sqrt(len*.25)-1.;
    return v[(int)floor(hi)];
}

/*****************
 * vdisp
 * ---------------------------------------------------------------------------------------
 * Finds and returns the high bound of 95% CI on median of vector passed in.
 * ---------------------------------------------------------------------------------------
 * v        sorted vector of data (flexible type)
 * p        percentile value to return (range 0 - 1) (note: passing in p = 1 will not work)
 * returns  high bound of 95% CI
 *****************/
template <typename T>
T vdisp (std::vector<T> &v, double p) {

    int len = v.size();
    if (p >= 0.5) return v[(int)ceil(p*len)];
    else return v[(int)floor(p*len)];
}

/*****************
 * bootstrap
 * ---------------------------------------------------------------------------------------
 * Uses resampling with replacement (bootstrapping) to determine confidence interval on
 * the median of vector v.  Confidence interval specified by ci (e.g. 0.95 = 95% CI)
 * ---------------------------------------------------------------------------------------
 * v        sorted vector of data (flexible type)
 * ci       confidence interval to return (0-1)
 * ci_low   lower value of confidence interval (will be updated when determiend)
 * ci_high  upper value of confidence interval (will be updated when determined)
 * returns  success/failure
 *****************/
template <typename T>
int bootstrap (std::vector<T> &v, double ci, double &ci_low, double &ci_high) {

    std::random_device rnd_device;
    std::mt19937 mersenne_engine(rnd_device());

    int n = 1000;    // number of samples to take
    int m = v.size();     // number of draws per sample
    std::uniform_int_distribution<int> dist(0,v.size()-1);
    auto gen = std::bind(dist,mersenne_engine);
    std::vector<int> rand_ind(m*n);
    std::vector<T> sample(m);
    std::vector<T> median;

    std::cout << "Called bootstrap method. Vector has " << v.size() << " elements." << std::endl;
    //std::cout << "Random indices to sample are: " << std::endl;
    std::cout << "Generating medians..." << std::endl;

    std::generate(rand_ind.begin(), rand_ind.end(), gen);

    for (int i = 0; i < n; i++) {

        //if (i<20) {
        //    for (int j = 0; j < m; j++) {
        //        std::cout << rand_ind[i*m+j] << " ";
        //    }
        //}

        for (int j = 0; j < m; j++) {
            sample[j] = v[rand_ind[i*m+j]];
        }
        std::sort(sample.begin(),sample.end());
        median.push_back(vmedian(sample));
        //if (i<20) {
        //    std::cout << median.back() << " ";
        //    std::cout << std::endl;
        //}
    }

    std::sort(median.begin(),median.end());
    
    std::cout << "Median of original sample was: " << vmedian(v) << ", median of medians is: " << vmedian(median) << std::endl;

    ci_low = median[(int)floor(n*0.025)];
    ci_high = median[(int)ceil(n*0.975)];

    std::cout << "95% Confidence interval is (" << ci_low << ", " << ci_high << ")" << std::endl;

    return -1;
}

#endif
