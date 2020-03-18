#ifndef __PROFILE_H__
#define __PROFILE_H__

#include "include.h"
#include "helper.h"
#include "io.h"
#include "stats.h"
#include <unordered_map>
#include <algorithm>
#include <queue>
#include <ctime>

#define FORCE_RES 1.0

struct new_profile_fields {
    double rs1, rs2, rs3;
    double c02, c03;
    double gamma3;
    double GOF1, GOF2, GOF3;
    double chi2_1, chi2_2, chi2_3;
    double np_rvir;
    double nb_rvir;
};

struct halo_profile {
    long int id;
    long int parent_id;
    double mass;
    double radius;
    double vmax;
    double x, y, z;
    double vx, vy, vz;
    int np;
    int n_rbins;
    int ppbin;

    // these must be listed last and in this order, since above fields will
    // be written/read as a block
    new_profile_fields newf;

    std::vector<double> bins;
    std::vector<double> weights;
};

struct bound {
    double lower;
    double upper;
};

int doStackProfiles (double ** &data, int num_lines, int num_fields);

int doProfileAnalysis (std::string name, bool binary_read);

#endif

