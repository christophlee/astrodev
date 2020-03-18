#include "io.h"

int readFile(std::string name, int &num_lines, int &num_fields, double ** &data, char delim) {
	std::ifstream inpt;
	std::string tmp_str;
	//char tmp_word_c[std::numeric_limits<int>::max()];
	//char tmp_line_c[std::numeric_limits<int>::max()];
	char tmp_word_c[MAX_WORD];
	char tmp_line_c[MAX_LINE];
    std::string tmp_line(tmp_line_c);
    std::string tmp_word(tmp_word_c);
	double tmp_d;

    std::cout << std::numeric_limits<int>::max() << std::endl;

    std::cout << "Reading file: \"" << name << "\"" << std::endl;

    inpt.open(name.c_str(),std::ifstream::in);

    int startline = 1;

	// lets check if this input file has a header
	// and skip it if it does
	int position = inpt.tellg();
	while (getline(inpt,tmp_str)) {
        startline++;
		if (tmp_str[0] != '#') {
            startline--;
			inpt.seekg(position);
			break;
		}
		position = inpt.tellg();
	}

	// lets find out how many columns of data we are dealing with
	bool at_delim = false;
	getline(inpt,tmp_str);
	for (int i = 0; i < tmp_str.length(); i++) {
		if ((tmp_str.at(i) == ' ') || (tmp_str.at(i) == '\t') || (tmp_str.at(i) == '\n') ||
            (tmp_str.at(i) == '\v') || (tmp_str.at(i) == '\f') || (tmp_str.at(i) == '\r') || (tmp_str.at(i) == ',')) {
			if (!at_delim) {
				at_delim = true;
			}
		}
		else {
			if (i == 0) num_fields++;
			else if (at_delim) {
				at_delim = false;
				num_fields++;
			}
		}
	}
	inpt.seekg(position);

    //std::cout << "Startline: " << startline << std::endl;
	std::cout << name + " num_fields: " << num_fields << std::endl;
		
	inpt.seekg(position);

	// lets scan through file and see how many lines of data we have
	inpt.seekg(0,inpt.end);
	long int file_eof = inpt.tellg();
    //std::cout << "file_eof: " << file_eof << std::endl;
	inpt.seekg(position);
	while (true) {
		if (inpt.tellg() == file_eof) break;
		getline(inpt,tmp_str);
		num_lines++;
        //if (num_lines % 1000000 == 0) std::cout << num_lines/1000000 << std::endl; //<< ": "+tmp_str << std::endl;
	}
	inpt.seekg(position);

	std::cout << name + " num_lines: " << num_lines << std::endl;

	// now allocate space for that data
    try {
	    data = new double*[num_fields];

    	for (int i = 0; i < num_fields; i++) {
            // note: using num_lines+1 seems to fix some seg faults etc on certain files
		    data[i] = new double[num_lines+1];
        }
    }
    catch (std::bad_alloc & ba) {
        std::cerr << "ERROR: bad_alloc exception caught: " << ba.what() << std::endl;
        return -1;
	}

	// finally, read file into memory
	int i = 0;
    int skipped_lines = 0;
    int dot_freq = int(ceil(num_lines/10));

    //std::istringstream line_stream;

    try {
	    while(inpt.good()) {
	    	//for (int j = 0; j < num_fields; j++) {
	    	//	inpt >> tmp_word;
            //    if (std::string(tmp_word).find_first_not_of("0123456789.e+-") != std::string::npos) {
            //        std::cerr << "Found bad value in data file: \n Line: " << i << "\tField: " << j << "\nEntry: " << tmp_word << std::endl;
            //        return -1;
            //    }
	    	//	tmp_d = std::strtod(tmp_word,NULL);
	    	//	data[j][i] = tmp_d;
	    	//}

            std::getline(inpt, tmp_line);
            std::istringstream line_stream(tmp_line);

            if (tmp_line.find_first_not_of("0123456789.eE+-, infa") != std::string::npos) {
                std::cerr << "Found bad value in data file: \nLine: " << i+startline+skipped_lines << "\nEntry:\n" << tmp_line << std::endl;
                num_lines--;
                skipped_lines++;
                continue;
                //return -1;
            }

            int j = 0;
            bool cont = false;
            //for (int j = 0; j < num_fields; j++) {
            while (getline(line_stream,tmp_word,delim)) {//line_stream >> tmp_word) {
                //line_stream >> tmp_word;
                if (!tmp_word.size()) continue;
                if (j >= num_fields) {
                    std::cerr << "Found too many fields at line: " << i+startline+skipped_lines << std::endl;
                    std::cerr << "Data:\n" << tmp_line << std::endl;
                    std::cerr << tmp_d << std::endl;
                    std::cerr << j << std::endl;
                    num_lines--;
                    skipped_lines++;
                    cont = true;
                    //break;
                    return -1;
                }

	    		tmp_d = std::strtod(tmp_word.c_str(),NULL);
	    		data[j++][i] = tmp_d;
            }

            if (cont) continue;
         
            if (j < num_fields) {
                if (inpt.good()) {
                    std::cerr << "Found too few fields at line: " << i+startline+skipped_lines << ".  Found " << j << ", expected " << num_fields << std::endl;
                    std::cerr << "Data:\n" << tmp_line << std::endl;
                    num_lines--;
                    skipped_lines++;
                    continue;
                    //return -1;
                }
            }

            if (dot_freq && i%dot_freq == 0) {
                std::cout << ".";
                std::cout.flush();
            }
	    	i++;
	    }
    }
    catch (std::exception & e) {
        std::cerr << "Exception caught: " << e.what() << std::endl;
        return -1;
    }

	inpt.close();

    std::cout << " complete." << std::endl;
    //std::cout << "num_lines: " << num_lines << std::endl;
    std::cout.flush();

    return 0;
}

//int readBinaryFile2 (std::string fileName, int &num_lines, int &num_fields, double ** & data) {
//
//    std::cout << "Reading binary file: \"" << fileName << "\"" << std::endl;
//
//    std::string tmp_str;
//
//    std::ifstream inpt(fileName,std::ios::in);
//
//    FILE * inpt = fopen(fileName.c_str(),'rb')
//
//	// lets check if this input file has a header
//	// and skip it if it does
//	int position = inpt.tellg();
//	while (getline(inpt,tmp_str)) {
//		if (tmp_str[0] != '#') {
//			inpt.seekg(position);
//			break;
//		}
//        else {
//            if (tmp_str.find("NUM_LINES") != -1) {
//                int start = tmp_str.find("NUM_LINES");
//                num_lines = std::strtod((tmp_str.substr(tmp_str.find(' ',start))).c_str(),NULL);
//                std::cout << "num_lines: " << num_lines << std::endl;
//            }
//            else if (tmp_str.find("NUM_FIELDS") != -1) {
//                int start = tmp_str.find("NUM_FIELDS");
//                num_fields = std::strtod((tmp_str.substr(tmp_str.find(' ',start))).c_str(),NULL);
//                std::cout << "num_fields: " << num_fields << std::endl;
//            }
//        }
//		position = inpt.tellg();
//	}
//
//	data = new double*[num_fields];
//
//    int dot_freq = int(ceil(num_fields/10.));
//
//    for (int i = 0; i < num_fields; i++) {
//        if (i % dot_freq == 0) {
//            std::cout << ".";
//            std::cout.flush();
//        }
//		data[i] = new double[num_lines];        
//        inpt.read((char *)data[i],sizeof(double)*num_lines);   
//    }
//
//    inpt.close();
//
//    std::cout << " complete." << std::endl;
//
//    /*for (int i = 0; i < num_fields; i++) {
//        std::cout << data[i][0] << std::endl;
//    }*/
//
//    return 0;
//}


int readBinaryFileParallel (std::string fileName, int & num_lines, int & num_fields, double ** & data, bool peek) {

    std::cout << "Reading binary file (in parallel): \"" << fileName << "\"" << std::endl;

    std::string tmp_str;

    std::ifstream inpt(fileName,std::ios::in);

	// lets check if this input file has a header
	// and skip it if it does
	int position = inpt.tellg();
	while (getline(inpt,tmp_str)) {
		if (tmp_str[0] != '#') {
			inpt.seekg(position);
			break;
		}
        else {
            if (tmp_str.find("NUM_LINES") != -1) {
                int start = tmp_str.find("NUM_LINES");
                num_lines = std::strtod((tmp_str.substr(tmp_str.find(' ',start))).c_str(),NULL);
                std::cout << "num_lines: " << num_lines << std::endl;
            }
            else if (tmp_str.find("NUM_FIELDS") != -1) {
                int start = tmp_str.find("NUM_FIELDS");
                num_fields = std::strtod((tmp_str.substr(tmp_str.find(' ',start))).c_str(),NULL);
                std::cout << "num_fields: " << num_fields << std::endl;
            }
        }
		position = inpt.tellg();
	}

    if (peek) return 0;

    inpt.close();

    data = new double*[num_fields];

    for (int i = 0; i < num_fields; i++) data[i] = new double[num_lines];

    int dot_freq = int(ceil(num_fields/10.));

    int counter = 0;

#pragma omp parallel shared(num_lines,num_fields,fileName,data,dot_freq,counter,position)
    {

    std::ifstream inpt_parallel(fileName,std::ios::in);

    #pragma omp for schedule(dynamic,4) nowait
    for (int i = 0; i < num_fields; i++) {
        
        // determine file offset based on iteration number
        inpt_parallel.seekg(position + sizeof(double)*num_lines*i);
        inpt_parallel.read((char *)data[i],sizeof(double)*num_lines);

        //if (omp_get_thread_num() == 0)
        //    std::cout << position << " " << inpt_parallel.tellg() << std::endl;

        //#pragma omp critical
        //if (++counter % dot_freq == 0) {
        //    std::cout << ".";
        //    std::cout.flush();
        //}
    }

    inpt_parallel.close();

    }

    std::cout << " complete." << std::endl;

    return 0;
}

int readBinaryFile (std::string fileName, int &num_lines, int &num_fields, double ** & data, bool peek) {

    std::cout << "Reading binary file: \"" << fileName << "\"" << std::endl;

    std::string tmp_str;

    std::ifstream inpt(fileName,std::ios::in);

	// lets check if this input file has a header
	// and skip it if it does
	int position = inpt.tellg();
	while (getline(inpt,tmp_str)) {
		if (tmp_str[0] != '#') {
			inpt.seekg(position);
			break;
		}
        else {
            if (tmp_str.find("NUM_LINES") != -1) {
                int start = tmp_str.find("NUM_LINES");
                num_lines = std::strtod((tmp_str.substr(tmp_str.find(' ',start))).c_str(),NULL);
                std::cout << "num_lines: " << num_lines << std::endl;
            }
            else if (tmp_str.find("NUM_FIELDS") != -1) {
                int start = tmp_str.find("NUM_FIELDS");
                num_fields = std::strtod((tmp_str.substr(tmp_str.find(' ',start))).c_str(),NULL);
                std::cout << "num_fields: " << num_fields << std::endl;
            }
        }
		position = inpt.tellg();
	}

    if (peek) return 0;

	data = new double*[num_fields];

    int dot_freq = int(ceil(num_fields/10.));

    for (int i = 0; i < num_fields; i++) {
        if (i % dot_freq == 0) {
            std::cout << ".";
            std::cout.flush();
        }
		data[i] = new double[num_lines];        
        inpt.read((char *)data[i],sizeof(double)*num_lines);
        //std::cout << inpt.tellg() << std::endl;
    }

    inpt.close();

    std::cout << " complete." << std::endl;

    /*for (int i = 0; i < num_fields; i++) {
        std::cout << data[i][0] << std::endl;
    }*/

    return 0;
}

int readBinary2File (std::string fileName, int &num_lines, int &num_fields, double ** & data, bool peek) {

    std::cout << "Reading binary2 file: \"" << fileName << "\"" << std::endl;

    std::string tmp_str;

    std::ifstream inpt(fileName,std::ios::in);

    long int header[2];

    // read header
    inpt.read((char *)header, 2*sizeof(long int));

    num_lines = header[0];
    num_fields = header[1];

    std::cout << "num_lines: " << num_lines << std::endl;
    std::cout << "num_fields: " << num_fields << std::endl;

    if (peek) return 0;

	data = new double*[num_lines];

    int dot_freq = int(ceil(num_lines/10.));

    for (int i = 0; i < num_lines; i++) {
        if (dot_freq && i % dot_freq == 0) {
            std::cout << ".";
            std::cout.flush();
        }
		data[i] = new double[num_fields];        
        inpt.read((char *)data[i],sizeof(double)*num_fields);
    }

    inpt.close();

    std::cout << " complete." << std::endl;

    //for (int i = 0; i < num_fields; i++) {
    //    std::cout << data[0][i] << std::endl;
    //}

    return 0;
}

//void doBinaryRead3D (std::string fileName, long double *** & data, int dim, int header = 0) {
//    void *** data_v;
//    doBinaryRead3D (fileName,data_v,dim,header,"long double");
//    data = (long double ***)data_v;
//}
//void doBinaryRead3D (std::string fileName, double *** & data, int dim, int header = 0) {
//    void *** data_v;
//    doBinaryRead3D (fileName,data_v,dim,header,"double");
//    data = (double ***)data_v;
//}
//void doBinaryRead3D (std::string fileName, char *** & data, int dim, int header = 0) {
//    void *** data_v;
//    doBinaryRead3D (fileName,data_v,dim,header,"char");
//    data = (char ***)data_v;
//}
//void doBinaryRead3D (std::string fileName, int *** & data, int dim, int header = 0) {
//    void *** data_v;
//    doBinaryRead3D (fileName,data_v,dim,header,"int");
//    data = (int ***)data_v;
//}

// type = 0 means data is in regular format.  type = 1 means data was read and store in binary2 format, so
// we need to reoder before writing.
int writeOutputFileBinary (std::string fileName, double ** &data, int num_lines, int num_fields, int type) {

    std::ofstream outfile;
    std::string file;
    
    if (!fileName.length()) {
        setOutputFileName();
    
        outfile.open(output.c_str(),std::ofstream::out);
        file = output;
    }

    else {
        outfile.open(fileName.c_str(),std::ofstream::out);
        file = fileName;
    }

    // check if we need to reorder data
    if (type == 1) {
        double ** data2;
        data2 = (double **) malloc (num_fields * sizeof(double *));
        for (int i = 0; i < num_fields; i++) {
            data2[i] = (double *) malloc (num_lines * sizeof(double));
        }

        // copy data over
        for (int i = 0; i < num_lines; i++) {
            for (int j = 0; j < num_fields; j++) {
                data2[j][i] = data[i][j];
            }
        }
        
        // now free data
        free (data, num_lines);

        // and redirect pointer to new data
        data = data2;
    }   

    std::cout << "Writing output file \"" << file << "\"" << std::endl;
    std::cout << "num_lines: " << num_lines << std::endl;
    std::cout << "num_fields: " << num_fields << std::endl;

    outfile << "# NUM_LINES " << num_lines << std::endl;
    outfile << "# NUM_FIELDS " << num_fields << std::endl;

    int dot_freq = int(ceil(num_fields/10.));

    for (int i = 0; i < num_fields; i++) {
        if (dot_freq && i % dot_freq == 0) {
            std::cout << ".";
            std::cout.flush();
        }
        outfile.write((const char *)data[i],num_lines*sizeof(double));
    }

    std::cout << " complete." << std::endl;

    outfile.close();
}


/********************************
* writeOutputFileBinary:  Writes output file in binary format, while inserting new fields at specified columns in the data.
*
* fileName:      File name of output file to be written.  If empty string, use default.
* data:          2D data to write to output file.  Consists of num_lines (rows) * num_fields (cols) entries.
* num_lines:     Number of rows in data.
* num_fields:    Number of cols in data.
* num_new_fields:Number of new fields to insert into data while writing.
* newdata:       2D array containing num_lines (rows) * num_new_fields (cols) to be inserted.
* insert_index:  1D array containing num_new_fields elements, where each element corresponds to the insertion index for the corresponding field in newdata.  The insertion index refers to the position in the unmodified data (before inserting any fields), and is assumed to be in order of insertion.
* *********************************/
int writeOutputFileBinary (std::string fileName, double ** &data, int num_lines, int num_fields, int num_new_fields, double ** &newdata, int * &insert_index) {

    std::ofstream outfile;
    std::string file;
    
    if (!fileName.length()) {
        setOutputFileName();
    
        outfile.open(output.c_str(),std::ofstream::out);
        file = output;
    }

    else {
        outfile.open(fileName.c_str(),std::ofstream::out);
        file = fileName;
    }

    std::cout << "Writing output file \"" << file << "\"" << std::endl;
    std::cout << "num_lines: " << num_lines << std::endl;
    std::cout << "num_fields: " << num_fields+num_new_fields << std::endl;

    outfile << "# NUM_LINES " << num_lines << std::endl;
    outfile << "# NUM_FIELDS " << num_fields+num_new_fields << std::endl;

    int dot_freq = int(ceil((num_fields+num_new_fields)/10.));

    int next_insert_field = 0;

    for (int i = 0; i < num_fields; i++) {
        if ((i+next_insert_field) % dot_freq == 0) {
            std::cout << ".";
            std::cout.flush();
        }

        // check if we have a field to insert at this index
        if (next_insert_field < num_new_fields) {
            if (i == insert_index[next_insert_field]) {
                outfile.write((const char *)newdata[next_insert_field++],num_lines*sizeof(double));
                // loop back with same i value to see if we have another new field
                // to insert at this index
                i--;
                continue;
            }
        }

        // if we did not insert any new fields, then write out next field from data
        outfile.write((const char *)data[i],num_lines*sizeof(double));
    }
    
    // write out any remaining new fields to be inserted at the end of the file
    while (next_insert_field < num_new_fields) {
        if (insert_index[next_insert_field] == num_fields) {
            outfile.write((const char *)newdata[next_insert_field++],num_lines*sizeof(double));
        }
    }

    std::cout << " complete." << std::endl;

    outfile.close();
}

int writeOutputFile (double ** &data, int num_lines, int num_fields, int type) {

    std::ofstream of;
    
    setOutputFileName();
    
    of.open(output.c_str(),std::ofstream::out);

    of << std::setprecision(10);

    std::cout << "Writing output file" << std::endl;

    int dot_freq = int(ceil(num_lines/10.));

    for (int i = 0; i < num_lines; i++) {

        if (dot_freq && i % dot_freq == 0) {
            std::cout << ".";
            std::cout.flush();
        }

        for (int j = 0; j < num_fields; j++) {

            switch (type) {

                case 0:

                    if (!j) of << data[0][i];

                    else {
                        // are we outputing this hlist style? (full halo id)
                        if (j == 1) of << " " << (unsigned int)data[j][i];
                        //else of << " " << data[j][i];

                        // I added the below specifically for reformatting the density augmented catalogs for public release
                        else {
                            if (j == 34) {
                                of << " " << data[89][i];
                                of << " " << data[90][i];
                                of << " " << data[91][i];
                                of << " " << data[j][i];
                            }
                            else if (j == 72) {
                                of << " " << data[92][i];
                                of << " " << data[j][i];
                            }
                            else if (j < 79) of << " " << data[j][i];
                        }
                        
                        // correction for T/|U| field
                        //else if (j == 53)
                        //    of << " " << 2*data[j][i];
                        // correction for nfw goodness of fit proxy as
                        // (r_s,nfw - r_s,klypin) / r_s,nfw
                        //else if (j == 82)
                        //    of << " " << (data[12][i] - data[34][i])/data[12][i];
                    }

                    break;

                case 1:

                    if (!j) of << data[i][0];

                    else {
                        // are we outputing this hlist style? (full halo id)
                        if (j == 1) of << " " << (unsigned int)data[i][j];
                        else of << " " << data[i][j];
                    }

                    break;
            }
        }

    of << std::endl;
    }

    of.close();

    std::cout << " complete." << std::endl;

    return 0;
}

int writeOutputFile (double ** &data, int num_lines, int num_fields, int type, int num_new_fields, double ** &newFieldVals, int * &newFieldInds) {

    std::ofstream of;
    
    setOutputFileName();
    
    of.open(output.c_str(),std::ofstream::out);

    of << std::setprecision(10);

    std::cout << "Writing output file" << std::endl;

    int dot_freq = int(ceil(num_lines/10.));
    int insert_index = 0;
    bool increment_insert_index = true;

    for (int i = 0; i < num_lines; i++) {

        if (dot_freq && i % dot_freq == 0) {
            std::cout << ".";
            std::cout.flush();
        }

        for (int j = 0; j < num_fields; j++) {

            // check if we got passed some data to insert
            if (newFieldVals) {
                for (int k = 0; k < num_new_fields; k++) {
                    if (newFieldInds[k] == j) {
                        increment_insert_index = false;
                        of << " " << newFieldVals[k][i];
                    }
                }
                if (increment_insert_index) insert_index++;
                else increment_insert_index = true;
            }

            if (!j) {
                of << data[0][i];
            }

            else {
                // are we outputing this hlist style? (full halo id)
                if (type == 0) {
                    if (j == 1)
                        of << " " << (unsigned int)data[j][i];
                    else
                        of << " " << data[j][i];
                }
                // correction for T/|U| field
                //else if (j == 53)
                //    of << " " << 2*data[j][i];
                // correction for nfw goodness of fit proxy as
                // (r_s,nfw - r_s,klypin) / r_s,nfw
                //else if (j == 82)
                //    of << " " << (data[12][i] - data[34][i])/data[12][i];
                else
                    of << " " << data[j][i];
            }
        }

    of << std::endl;
    }

    of.close();

    std::cout << " complete." << std::endl;

    return 0;
}
