#ifndef __IO_H__
#define __IO_H__

#include "include.h"
#include "helper.h"
#include <omp.h>
#include <sstream>
#include <limits>
#include <iomanip>
#include <algorithm>

#define MAX_WORD 1048576    // max number of bytes in character buffer
                            // for reading files
#define MAX_LINE 1048576    // same as above but for whole lines.
                            // 1024 * 1024 bytes = 1MB

int readFile(std::string name, int &num_lines, int &num_fields, double ** &data, char delim = ' ');

int readBinaryFileParallel (std::string fileName, int & num_lines, int & num_fields, double ** & data, bool peek = false);

int readBinaryFile (std::string fileName, int &num_lines, int &num_fields, double ** & data, bool peek = false);

int readBinary2File (std::string fileName, int &num_lines, int &num_fields, double ** & data, bool peek = false);

int writeOutputFileBinary (std::string fileName, double ** &data, int num_lines, int num_fields, int type=0);

int writeOutputFileBinary (std::string fileName, double ** &data, int num_lines, int num_fields, int num_new_fields, double ** &newdata, int * &insert_index);

int writeOutputFile (double ** &data, int num_lines, int num_fields, int type = 0);

int writeOutputFile (double ** &data, int num_lines, int num_fields, int type, int num_new_fields, double ** & newFieldVals, int * &newFieldInds);


//void doBinaryRead3D (std::string fileName, long double *** & data, int dim, int header = 0) {

// NOTE: header is number of bytes in header
template <typename T>
int doBinaryRead3D (std::string fileName, T *** & data, int dim, int header = 0) {//, std::string type) {

    std::cout << "Reading binary file: \"" << fileName << "\"" << std::endl;

    std::ifstream inpt(fileName,std::ios::in);

    int dot_freq = int(ceil(dim/10.));

    // skip header (if there is one)
    inpt.seekg(header);

    //data = new void ** [dim];
    data = (T ***) malloc (dim * sizeof (T **));
    for (int i = 0; i < dim; i++) {
        if (i % dot_freq == 0) {
            std::cout << ".";
            std::cout.flush();
        }
        data[i] = (T **)malloc(dim * sizeof(T *));
        //data[i] = new void * [dim];
        for (int j = 0; j < dim; j++) {

            //int size = sizeof(T);

            //if (type == "long double") {
            //    size = sizeof(long double);
            //}
            //else if (type == "int") {
            //    size = sizeof(int);
            //}
            //else if (type == "char") {
            //    size = sizeof(char);
            //}
            //else if (type == "double") {
            //    size = sizeof(double);
            //}

            data[i][j] = (T *)malloc(dim * sizeof(T));
            memset(data[i][j],0,dim*sizeof(T));
            //for (int k = 0; k < dim; k++) data[i][j][k] = 0;
    
            inpt.read((char *)data[i][j],sizeof(T)*dim);
            //if (type == "double")
            //    inpt.read((char *)data[i][j],sizeof(double)*dim);
            //else if (type == "long double")
            //    inpt.read((char *)data[i][j],sizeof(long double)*dim);
            //else if (type == "float")
            //    inpt.read((char *)data[i][j],sizeof(float)*dim);
            //else if (type == "int")
            //    inpt.read((char *)data[i][j],sizeof(int)*dim);
            //else if (type == "char")
            //    inpt.read((char *)data[i][j],sizeof(char)*dim);
            //else
            //    std::cerr << "doBinaryRead3D error: unrecognized data type: " << type << std::endl;

        }   
    }

    inpt.close();

    std::cout << " complete." << std::endl;

    return 0;
}

#endif
