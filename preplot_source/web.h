#ifndef __WEB_H__
#define __WEB_H__

#include "include.h"
#include "helper.h"
#include "io.h"
#include "stats.h"
#include <unordered_map>
#include <algorithm>

int doCombineSpineCatalog (double ** &data, int num_lines, int num_fields, bool binary_write);

int doSpineCatalogAnalysis (double ** &data, int num_lines, int num_fields, int x = 0);

int doSpineTest (double ** &data, int num_lines, int num_fields);

#endif
