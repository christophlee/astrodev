#ifndef __INCLUDE_H__
#define __INCLUDE_H__

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <vector>

struct BP_CAT_DEF {
    int a = 0;
    int id = 1;
    int mvir = 10;
    int rvir = 11;
    int rs = 12;
    int almm = 15;
    int vmax = 16;
    int x = 17;
    int y = 18;
    int z = 19;
    int spinp = 26;
    int orig_id = 30;
    int rsk = 34;
    int xoff = 40;
    int voff = 41;
    int spinb = 42;
    int bta = 43;
    int cta = 44;
    int ax = 45;
    int ay = 46;
    int az = 47;
    int bta500 = 48;
    int cta500 = 49;
    int ax500 = 50;
    int ay500 = 51;
    int az500 = 52;
    int tu = 53;
    int mpeak = 57;
    int amhalf = 60;
    int mar_inst = 61;
    int rhocic = 72;
    int rhohalf = 73;
    int rho1 = 74;
    int rho2 = 75;
    int rho4 = 76;
    int rho8 = 77;
    int rho16 = 78;
    int cnfw = 79;
    int bs_ratio = 80;
    int pvir = 82;
    int p500 = 83;
    int smar_inst = 84;
    int smar_dyn = 85;
    int cklyp = 88;
    int tf = 90;
    int tf_dyn = 92;

// tree analyzed fields.
// these are appended to catalog only after tree-analyzed and merged
    //int mpeak = 93;
    int a_peak = 94;
    int mdot_dyn = 95;
    int npt = 96;
    int high_tf = 97;
    int maxtf_peak = 98;
    int a_maxtf_peak = 99;
    int min_bsr_peak = 100;
    int a_min_bsr_peak = 101;
    int min_bsr_lmm = 102;
    int a_min_bsr_lmm = 103;
    int delta_a_peak = 104;
    int delta_a_lmm = 105;
};

struct BP_TREE_DEF {
    int a = 0;
    int id = 1;
    int desc_id = 3;
    int num_prog = 4;
    int upid = 6;
    int mvir = 10;
    int rvir = 11;
    int rs = 12;
    int almm = 15;
    int vmax = 16;
    int x = 17;
    int y = 18;
    int z = 19;
    int spinp = 26;
    int orig_id = 30;
    int tf = 35;
    int rsk = 37;
    int xoff = 43;
    int voff = 44;
    int spinb = 45;
    int bta = 46;
    int cta = 47;
    int ax = 48;
    int ay = 49;
    int az = 50;
    int bta500 = 51;
    int cta500 = 52;
    int ax500 = 53;
    int ay500 = 54;
    int az500 = 55;
    int tu = 56;
    int mpeak = 59;
    int a_peak = 60;
    int mdot_dyn = 61;
    int npt = 62;
    int high_tf = 63;
    int maxtf_peak = 64;
    int a_maxtf_peak = 65;
    int min_bsr_peak = 66;
    int a_min_bsr_peak = 67;
    int min_bsr_lmm = 68;
    int a_min_bsr_lmm = 69;
    int spin_peak = 70;
    int tu_peak = 71;
    int mass_lost_as_sh = 72;  // mass lost while a subhalo since mpeak

    // the below are not being used for now, but are included so tree.cpp will compile
    int almms = 60;
    int mlumpy = 61;
    int mlumpy_frac = 62;
    int msmooth = 63;
    int msmooth_frac = 64;
    int merger_count = 65;
    int merger_mass_avg = 66;
    //int high_tf = 67;
    //int mdot_dyn = 68;
    int mdot_dyn_lumpy = 69;
    int mdot_dyn_smooth = 70;
};

struct threshold_object {
    double * minValues = NULL;
    double * maxValues = NULL;
    bool range = false;
    bool percent = false;
};

extern std::string output;
extern BP_TREE_DEF BP_TREE;
extern BP_CAT_DEF BP_CAT;

#endif
