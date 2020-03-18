#ifndef __HELPER_H__
#define __HELPER_H__

#include "include.h"
#include <limits>

//extern std::string output;
extern std::string fileName;

void free (double ** &data, int &num_fields);
void free (double ** & data, long long & num_fields);
void free (double *** &data, int dim1, int dim2);
void free (float *** & data, int dim1, int dim2);
void free (long double *** &data, int dim1, int dim2);

std::string str (int value);

void toggle_stream (std::ostream & stream);
void off_stream (std::ostream & stream);
void on_stream (std::ostream & stream);

bool isBinaryExt (std::string fileName);
bool isBinary2Ext (std::string fileName);

void parseOutputFileName (std::string inptFile);

void handleSameInputOutputName();

void setOutputFileName(int x, int y);
void setOutputFileName(int x);
void setOutputFileName(bool b);
void setOutputFileName();

int handleRangeExpression (std::string expr, double &minVal, double &maxVal, std::string &errmsg);

#endif
