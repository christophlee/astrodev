#ifndef __TREE_H__
#define __TREE_H__

#include "include.h"
#include "helper.h"
#include "io.h"
#include "fileop.h"
#include <unordered_map>
#include <algorithm>
#include <functional>
#include <cstdio>

//extern std::string output;

int doTreeWalk (double ** &data, int num_lines, int num_fields, bool mmp, std::string outfile_name = "");

int doProgenitorHistory (double ** &data, int num_lines, int num_fields);

int doSelectZ (double ** &data, int num_lines, int num_fields, int x, std::string outfile_name = "");

int doConvertMergerTracks (double ** &data, int num_lines, int num_fields);

int doTreeAnalysis (double ** &data, int num_lines, int num_fields, std::string outfile_name = "");

int doIndvHaloTracks (double ** &data, int num_lines, int num_fields, int x);

int doFullTreeAnalysis (double ** &data, int num_lines, int num_fields, bool mmp); 

#endif
