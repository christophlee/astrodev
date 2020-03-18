#include "helper.h"

void free (double ** & data, long long & num_fields) {

	// delete and free memory
	for (long long j = 0; j < num_fields; j++) {
		if (data[j]) delete[] data[j];
	}
	delete[] data;
}

void free (double ** & data, int & num_fields) {

	// delete and free memory
	for (int j = 0; j < num_fields; j++) {
		if (data[j]) delete[] data[j];
	}
	delete[] data;
}

void free (double *** & data, int dim1, int dim2) {

    for (int i = 0; i < dim1; i++) {
        for (int j = 0; j < dim2; j++) {
            delete [] data[i][j];
        }
        delete [] data[i];
    }
    delete [] data;
}

void free (float *** & data, int dim1, int dim2) {

    for (int i = 0; i < dim1; i++) {
        for (int j = 0; j < dim2; j++) {
            delete [] data[i][j];
        }
        delete [] data[i];
    }
    delete [] data;
}

void free (long double *** & data, int dim1, int dim2) {

    for (int i = 0; i < dim1; i++) {
        for (int j = 0; j < dim2; j++) {
            delete [] data[i][j];
        }
        delete [] data[i];
    }
    delete [] data;
}

std::string str (int value) {return std::to_string((long long)value);}

void toggle_stream (std::ostream & stream) {
    if (stream.fail() && !stream.bad()) stream.clear();
    else stream.setstate(std::ios_base::failbit);
}

void off_stream (std::ostream & stream) {
    stream.setstate(std::ios_base::failbit);
}

void on_stream (std::ostream & stream) {
    stream.clear();
}

bool isBinaryExt (std::string fileName) {

    int ext_pos = fileName.find_last_of('.');
    if (ext_pos != std::string::npos)
        if (fileName.substr(ext_pos+1) == "bin") return true;

    return false;
}

bool isBinary2Ext (std::string fileName) {

    int ext_pos = fileName.find_last_of('.');
    if (ext_pos != std::string::npos)
        if (fileName.substr(ext_pos+1) == "bin2") return true;

    return false;
}

void parseOutputFileName (std::string inptFile) {

    if (output.find('*') != std::string::npos) {
        // get name of input file without extension
        int ext_pos = inptFile.find_last_of('.');
        if (ext_pos != std::string::npos)
            output = inptFile.substr(0,ext_pos)+output.substr(output.find('*')+1);
    }
}

void handleSameInputOutputName() {
    char cin_char;
    if (output != fileName) return;
    std::cout << "Warning: output file would overwrite input file, proceed (y/n)?:";
    std::cin >> cin_char;
    if (cin_char != 'y') {
        std::cout << "Understood, appending \".out\" to output file name." << std::endl;
        output = fileName + ".out";
    }
    else {
        std::cout << "Understood, proceeding to overwrite input file." << std::endl;
    }
}

void setOutputFileName (int x, int y) {

	if (!output.length()) {
		output = fileName+".out_"+std::to_string((long long)x)+"_"+std::to_string((long long)y);
	}
    else handleSameInputOutputName();
}

void setOutputFileName (int x) {

	if (!output.length()) {
		output = fileName+".out_"+std::to_string((long long)x);
	}
    else handleSameInputOutputName();
}

// overloaded function which does not call handleSameInputOutputName.
// used in doBin.
void setOutputFileName (bool b) {

    if (!output.length()) {
        output = fileName+".out";
    }
}

void setOutputFileName () {

    if (!output.length()) {
        output = fileName+".out";
    }
    else handleSameInputOutputName();
}

int handleRangeExpression (std::string expr, double &minVal, double &maxVal, std::string &errmsg) {

    std::string arg1, arg2;

    //std::cout << "In handle range expr with expr: " << expr << std::endl;

    // parse expression

    // find colon
    int pos = expr.find(':');
    arg1 = expr.substr(0,expr.find(':'));
    if (pos == expr.length()) arg2 = "";
    else arg2 = expr.substr(expr.find(':')+1);

    //std::cout << arg1 << " , " << arg2 << std::endl;

    if (!arg1.length()) minVal = -std::numeric_limits<double>::infinity();
    else minVal = std::strtod(arg1.c_str(),NULL);

    if (!arg2.length()) maxVal = std::numeric_limits<double>::infinity();
    else maxVal = std::strtod(arg2.c_str(),NULL);

    return 0;
}

