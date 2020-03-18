#include "fileop.h"

int doBin(double ** &data, int num_lines, int num_fields, int x, int num_bins, threshold_object &thresholds, int threshold_type, bool binary_write) {

    // first we need to take all of the data and sort it
    std::vector<double> data_x;

    // populate vector with data
    for (int i = 0; i < num_lines; i++) {
        if (std::isfinite(data[x][i])) {
            data_x.push_back(data[x][i]);
        }
    }

    std::sort(data_x.begin(), data_x.end());

    std::cout << "Data min: " << data_x[0] << std::endl; 
    std::cout << "Data max: " << data_x[data_x.size()-1] << std::endl; 

    doBin_data_min = data_x[0];
    doBin_data_max = data_x[data_x.size()-1];

    setOutputFileName(false);

    // now lets create our output files
    std::ofstream * outputs;
    outputs = new std::ofstream[num_bins];

    // determine output file names and open
    std::string outfileNames[num_bins];

    if ((num_bins == 1) && (output == fileName)) {
        char cin_char;
        std::cout << "Warning: output file would overwrite input file, proceed (y/n)?:";
        std::cin >> cin_char;
        if (cin_char != 'y') return -1;
    }

    for (int i = 0; i < num_bins; i++) {
        int ind = output.find_last_of('.');
        if (num_bins > 1) outfileNames[i] = output.substr(0,ind)+"_"+std::to_string((long long)i)+output.substr(ind);
        else outfileNames[i] = output;
    }

    if (!binary_write) {

        for (int i = 0; i < num_bins; i++) {
            outputs[i].open(outfileNames[i],std::ofstream::out);
        }
    }

    std::cout << "Threshold_type: " << threshold_type << std::endl;

    //std::cout << "num_lines: " << num_lines << ", num_fields: " << num_fields << ", num_bins: " << num_bins << ", x: " << x << std::endl;

    std::cout << "Thresholds: ";

    // find binning thresholds
    int num_per_bin = (data_x.size() - (data_x.size() % num_bins))/num_bins;
    if (!thresholds.minValues) {
        thresholds.minValues = new double[num_bins];
        thresholds.maxValues = new double[num_bins];

        for (int i = 0; i < num_bins; i++) {
            thresholds.minValues[i] = data_x[i*num_per_bin];
            if (i != num_bins-1) thresholds.maxValues[i] = data_x[(i+1)*num_per_bin];
            else thresholds.maxValues[i] = data_x[data_x.size()-1];
        }
    }

    else {
        if (threshold_type == THREEQ_BINNING) {
            thresholds.minValues[0] = data_x[0];
            thresholds.minValues[1] = data_x[floor(data_x.size()*.25)];
            thresholds.minValues[2] = data_x[floor(data_x.size()*.75)];
            thresholds.maxValues[0] = thresholds.minValues[1];
            thresholds.maxValues[1] = thresholds.minValues[2];
            thresholds.maxValues[2] = data_x[data_x.size()-1];
        }

        else if (thresholds.range) {
            // go through and fix any extreme values (the inf values are used as placeholders to
            // indicate that they need to be updated to the proper min/max values

            //  UPDATE: the +- inf values work fine if those are the intended boundaries, will just leave them be

            // we'll handle ranged percentiles slightly differently than non-ranged here
            if (threshold_type == PERCENT_BINNING) {

                // this is a hard coded option for using the whole box density values rather than
                // the provided catalog data to determine actual bin values.
                // this will overwrite the min/max thresholds already set above.
                //bool doWholeBoxPercentileConversion = true;
                bool doWholeBoxPercentileConversion = false;
    
                if (doWholeBoxPercentileConversion) {
    
                    for (int i = 0; i < num_bins; i++) {
                        if (std::isinf(thresholds.minValues[i])) thresholds.minValues[i] = data_x[0];//pdf_integrate(0.,0.,x);
                        else thresholds.minValues[i] = pdf_integrate(thresholds.minValues[i],0.,x);
                        if (std::isinf(thresholds.maxValues[i])) thresholds.maxValues[i] = data_x[data_x.size()-1];//pdf_integrate(100.,0.,x);
                        else thresholds.maxValues[i] = pdf_integrate(thresholds.maxValues[i],0.,x);
                    }
                }

                // this is the normal route, when not doing the special binning above
                else {
                    for (int i = 0; i < num_bins; i++) {
                        if (std::isinf(thresholds.minValues[i])) thresholds.minValues[i] = data_x[0];
                        else thresholds.minValues[i] = data_x[floor(data_x.size()*(thresholds.minValues[i]/100.))];
                        if (std::isinf(thresholds.maxValues[i])) thresholds.maxValues[i] = data_x[data_x.size()-1];
                        else thresholds.maxValues[i] = data_x[floor(data_x.size()*(thresholds.maxValues[i]/100.))];
                    }
                }
            }

            else {
                for (int i = 0; i < num_bins; i++) {
                    //if (std::isinf(thresholds.minValues[i])) thresholds.minValues[i] = data_x[0];
                    //if (std::isinf(thresholds.maxValues[i])) thresholds.maxValues[i] = data_x[data_x.size()-1];
                    continue;
                }
            }


        }

        // handle non-ranged percentiles
        else if (threshold_type == PERCENT_BINNING) {
            for (int i = 0; i < num_bins; i++) {
                if (!i) {
                    thresholds.minValues[0] = data_x[0];
                }
                else {
                    thresholds.minValues[i] = data_x[floor(data_x.size()*(thresholds.minValues[i]/100.))];
                    thresholds.maxValues[i-1] = thresholds.minValues[i];
                }

                if (i == num_bins-1) thresholds.maxValues[i] = data_x[data_x.size()-1];
            }
        }



        // thresholds must be standard, value specified thresholds
        else {
            for (int i = 0; i < num_bins; i++) {
                if (!i) {
                    thresholds.minValues[i] = data_x[0];
                }
                else {
                    thresholds.maxValues[i-1] = thresholds.minValues[i];
                }
                if (i == num_bins-1) thresholds.maxValues[i] = data_x[data_x.size()-1];
            }
        }
    }

    for (int i = 0; i < num_bins; i++) std::cout << std::setprecision(4) << std::fixed << "\t[" << thresholds.minValues[i] << ":" << thresholds.maxValues[i] << "]";
    std::cout << std::endl;

    // finally, determine which bin each galaxy falls into, and output to appropriate file
    /*for (int i = 0; i < num_lines; i++) {
        int j;
        for (j = 0; j < num_bins; j++) {
                if ((j < (num_bins-1)) && (thresholds[j+1] > data[x][i])) break;
                else if (j == num_bins-1) break;
                else continue;
        }

        for (int k = 0; k < num_fields; k++) {
            if (!k) outputs[j] << data[k][i];
            else outputs[j] << " " << data[k][i];
        }
        outputs[j] << std::endl;
    }*/

    int dot_freq = int(ceil(num_lines/10.));

    if (!binary_write) {
        std::cout << "Writing output files" << std::endl;
        for (int i = 0; i < num_lines; i++) {

            if (i % dot_freq == 0) {
                std::cout << ".";
                std::cout.flush();
            }

            int j;
            for (j = 0; j < num_bins; j++) {
                if ((thresholds.minValues[j] < data[x][i]) && (thresholds.maxValues[j] >= data[x][i])) {
                        for (int k = 0; k < num_fields; k++) {
                            if (!k) outputs[j] << data[k][i];
                            else outputs[j] << " " << data[k][i];
                        }
                        outputs[j] << std::endl;
                }
            }
        }
        std::cout << std::endl;
    }

    else {
        double ** binnedData = new double*[num_fields];
        for (int i = 0; i < num_fields; i++) {
            binnedData[i] = new double[num_lines];
        }
        for (int j = 0; j < num_bins; j++) {
            int counter = 0;
            for (int i = 0; i < num_lines; i++) {
                if ((thresholds.minValues[j] < data[x][i]) && (thresholds.maxValues[j] >= data[x][i])) {
                    for (int k = 0; k < num_fields; k++) {
                        binnedData[k][counter] = data[k][i];
                    }
                    counter++;
                }
            }
            writeOutputFileBinary(outfileNames[j],binnedData,counter,num_fields);
        }

        for (int i = 0; i < num_fields; i++)
            delete [] binnedData[i];
        delete [] binnedData;
    }

    // cleanup
    if (!binary_write) {
        for (int i = 0; i < num_bins; i++) {
            outputs[i].close();
        }

        delete [] outputs;
    }

    delete [] thresholds.minValues;
    delete [] thresholds.maxValues;

    return 0;
}

int doBinnedAnalysis (double ** &data, int num_lines, int num_fields, int binx, int num_bins, threshold_object &thresholds, int threshold_type, int x, int y, int num_bins_x, int num_bins_y, double bin_width_x, double bin_width_y, bool * LOG, bool NORM, bool SMOOTH, bool binary_write, bool autobin, int auto_tnum, bool use_percentiles_x, bool use_percentiles_y, bool RANK, bool CDF, int use_rw, bool MASSFUNC, int POINTS_PER_CELL, bool SCATTER) {

    // first lets just bin this puppy. if we are using reweighted files, the binning has been precomputed, so
    // don't do it here also.
    if (use_rw == -1) {
        doBin (data, num_lines, num_fields, binx, num_bins, thresholds, threshold_type, binary_write);
    }

    // now those files should be done, just pass them through the hist2D routine
    std::string orig_output = output;

    std::cout << "doing binned analysis" << std::endl;

    for (int i = 0; i < num_bins; i++) {
        int ind = orig_output.find_last_of('.');
    
        std::string file_name;
        std::string file_name_rw;

        // normal naming convention. just insert bin number before file suffix
        if (use_rw == -1) {
            if (num_bins == 1) {
                file_name = orig_output;
            }
            else {
                file_name = orig_output.substr(0,ind)+"_"+std::to_string((long long)i)+orig_output.substr(ind);
            }
        }

        // special naming convention for reweighted mass bins. insert "_rw$i" after hmbins,
        // where i is rw number, in addition to bin number before suffix.
        // NOTE: since we are not binning, there is no need to first read non-binned data.  Thus, we should
        // pass in -file=None or equivalent, and -output=bp_z0_centrals_hmbins_cnfw.bin or equivalent.
        else {
            int ind2 = orig_output.find("hmbins")+6;
            // this will be used to pass in to fileReaders, so should just be equiv. bp_z0_centrals_hmbins_0_rw1.bin
            file_name_rw = orig_output.substr(0,ind2)+"_"+str(i)+"_rw"+str(use_rw)+orig_output.substr(ind);
            // this will be used for output files
            file_name = orig_output.substr(0,ind2)+"_rw"+str(use_rw)+orig_output.substr(ind2);
            ind = file_name.find_last_of('.');
            file_name = file_name.substr(0,ind)+"_"+str(i)+file_name.substr(ind);
        }

        ind = file_name.find_last_of('.');

        // change global output for each file so hist2D outputs correctly
        if (binnedHist2D) output = file_name.substr(0,ind)+"_hist2D"+file_name.substr(ind);
        else if (binnedHist1D && !CDF) output = file_name.substr(0,ind)+"_hist1D"+file_name.substr(ind);
        else if (binnedHist1D && CDF) output = file_name.substr(0,ind)+"_cdf"+file_name.substr(ind);
        else if (binnedMedian) output = file_name.substr(0,ind)+"_med"+file_name.substr(ind);

        double ** data_2;
        int num_lines_2=0, num_fields_2=0;

        // read in each file and pass appropraite fields to hist2D
        // note: if these files were written in binary, we must read back in binary
        if (use_rw != -1) file_name = file_name_rw;
        if (binary_write) readBinaryFile (file_name,num_lines_2,num_fields_2,data_2);
        else readFile (file_name,num_lines_2,num_fields_2,data_2);

        // if using reweighted files, need to handle log of input at this stage, since didn't have a chance
        // to do it on pre-binned catalog.
        if (use_rw != -1) {
            if (LOG[0]) for (int j = 0; j < num_lines_2; j++) data_2[x][j] = log10(data_2[x][j]);
            if (LOG[1]) for (int j = 0; j < num_lines_2; j++) data_2[y][j] = log10(data_2[y][j]);
        }

        if (binnedHist2D) doHist2D (data_2, num_lines_2, x, y, num_bins_x, num_bins_y, bin_width_x, bin_width_y, LOG[2], NORM, SMOOTH,RANK, POINTS_PER_CELL, SCATTER);
        else if (binnedHist1D) doHist1D (data_2, num_lines_2, x, num_bins_x, NORM, LOG[2], CDF, autobin, auto_tnum, MASSFUNC, bin_width_x);
        else if (binnedMedian) doMedian (data_2, num_lines_2, x, y, num_bins_x, autobin, auto_tnum,use_percentiles_x,use_percentiles_y);

        free (data_2, num_fields_2);
    }

    return 0;
}

int doAddFieldToCatalog (double ** &data, int num_lines, int num_fields, bool binary_write) {

    //std::cout << "ADD FIELD TO CATALOG" << std::endl;
    //std::cout.flush();
    
    //int numNewFields = 10;
    int numNewFields = 4;

    double ** newFields  = new double * [numNewFields];
    for (int i = 0; i < numNewFields; i++)
        newFields[i] = new double[num_lines];

    int * newFieldInds = new int[numNewFields];

    double mcoll = 0;
    double a = data[0][0];

    if      (a == 1.00231) mcoll = 12.703;
    else if (a == 0.67325) mcoll = 11.9711;
    else if (a == 0.50112) mcoll = 11.2092;
    else if (a == 0.33406) mcoll = 9.82278;
    else if (a == 0.24800) mcoll = 8.67282;
    else if (a == 0.20244) mcoll = 7.71916;


    auto time = [](double a)->double {
        return 11.5570*asinh(1.50244*pow(a,3./2.));
    };

    for (int i = 0; i < num_lines; i++) {
        double ca500 = data[49][i];
        double cavir = data[44][i];
        double value = 0;
        //double sfr_ave = data[24][i];
        //double mstar = data[9][i];
        double accrateinst = data[61][i];
        double accratedyn = data[63][i];
        double accrate2dyn = data[64][i];
        double rs_nfw = data[12][i];
        double rs_klypin = data[34][i];
        double rvir = data[11][i];
        double mpeak = data[57][i];
        double mvir = data[10][i];
        double am[3] = {data[23][i],data[24][i],data[25][i]};
        double shape[3] = {data[45][i],data[46][i],data[47][i]};
        double s = data[44][i];
        double q = data[43][i];
        double s500 = data[49][i];
        double q500 = data[48][i];
        double cost = am[0]*shape[0]+am[1]*shape[1]+am[2]*shape[2];
        cost /= sqrt(am[0]*am[0]+am[1]*am[1]+am[2]*am[2]);
        cost /= sqrt(shape[0]*shape[0]+shape[1]*shape[1]+shape[2]*shape[2]);
        cost = fabs(cost);
        double Prvir = 1-sqrt(s*s+q*q)/sqrt(2);//(1-q*q)/(1-s*s);
        double Pr500 = 1-sqrt(s500*s500+q500*q500)/sqrt(2);//(1-q500*q500)/(1-s500*s500);
        //value = log10(sfr_ave/pow(10.,mstar));
        double mar_inst = accrateinst/mvir;
        double mar_dyn = accratedyn/mvir;
    
        newFields[0][i] = data[BP_CAT.a_min_bsr_peak][i] - data[BP_CAT.a_peak][i];
        newFields[1][i] = data[BP_CAT.a_min_bsr_lmm][i] - data[BP_CAT.almm][i];
        newFields[2][i] = time(data[BP_CAT.a_min_bsr_peak][i]) - time(data[BP_CAT.a_peak][i]);
        newFields[3][i] = time(data[BP_CAT.a_min_bsr_lmm][i]) - time(data[BP_CAT.almm][i]);

        //// C_NFW
        //newFields[0][i] = rvir/rs_nfw;
        //// BS ratio
        //newFields[1][i] = mvir/mpeak;
        //// cos(theta)
        //newFields[2][i] = cost;
        //// P_rvir
        //newFields[3][i] = Prvir;
        //// P_r500
        //newFields[4][i] = Pr500;
        //// MAR_inst
        //newFields[5][i] = mar_inst;
        //// MAR_dyn
        //newFields[6][i] = mar_dyn;
        //// MAR_2dyn
        //newFields[7][i] = accrate2dyn/mvir;
        //// log M/M*
        //newFields[8][i] = log10(mvir)/mcoll;
        //// C_klypin
        //newFields[9][i] = rvir/rs_klypin;

    }

    // insertion indeces for new fields
    for (int j = 0; j < numNewFields; j++) newFieldInds[j] = num_fields;

    if (binary_write)
        writeOutputFileBinary ("", data, num_lines, num_fields, numNewFields, newFields, newFieldInds);

    else
        writeOutputFile (data, num_lines, num_fields, 0/*type: 0==hlist*/, numNewFields, newFields, newFieldInds);

}

int doStackCatalogs () {

    int num_stacks = 5;     // number of bins
    int num_per_stack = 5;  // number of redshifts

    std::string z[] = {"0","0.5","1","2","3"};

    for (int i = 0; i < num_stacks; i++) {

        std::cout << "Stacking bin " << i << "... " << std::endl;

        int total_num_lines = 0;
        int total_num_fields = 0;

        for (int j = 0; j < num_per_stack; j++) {

            double ** hlist_data;
            int num_lines_2=0, num_fields_2=0;

            std::string hlist_file = "/pfs/chtlee/workspace/bp_z"+z[j]+"_centrals_mmstar_"+std::to_string((long long)i) + ".bin";

            readBinaryFile(hlist_file,num_lines_2,num_fields_2,hlist_data,true);

            total_num_lines += num_lines_2;
            total_num_fields = num_fields_2;
        }

        double stacked[total_num_fields][total_num_lines];

        long counter = 0;

        for (int j = 0; j < num_per_stack; j++) {

            double ** hlist_data;
            int num_lines_2=0, num_fields_2=0;

            std::string hlist_file = "/pfs/chtlee/workspace/bp_z"+z[j]+"_centrals_mmstar_"+std::to_string((long long)i) + ".bin";

            readBinaryFile(hlist_file,num_lines_2,num_fields_2,hlist_data);

            // lets scale some of our data.
            // scale mar by (1+z)-(5/2.)

            double z_d = std::strtod(z[j].c_str(),NULL);

            for (int k = 0; k < num_lines_2; k++) 
                hlist_data[85][k] /= pow(1.+z_d,5./2.);

            for (int k = 0; k < num_lines_2; k++) 
                hlist_data[86][k] /= pow(1.+z_d,5./2.);

            // scale c by (1+z)
            for (int k = 0; k < num_lines_2; k++) 
                hlist_data[79][k] *= 1.+z_d;

            for (int k = 0; k < num_fields_2; k++) {
                for (int l = 0; l < num_lines_2; l++) {
                    stacked[k][l+counter] = hlist_data[k][l];
                }
            }

            counter += num_lines_2;

            free(hlist_data, num_fields_2);
        }

        std::ofstream outfile;
        int ind = output.find_last_of('.');
        std::string new_output = output.substr(0,ind)+"_"+std::to_string((long long)i)+output.substr(ind);
        outfile.open(new_output.c_str(),std::ofstream::out);

        std::cout << "Writing output file \"" << new_output << "\"" << std::endl;
        std::cout << "num_lines: " << total_num_lines << std::endl;
        std::cout << "num_fields: " << total_num_fields << std::endl;

        outfile << "# NUM_LINES " << total_num_lines << std::endl;
        outfile << "# NUM_FIELDS " << total_num_fields << std::endl;

        int dot_freq = int(ceil(total_num_fields/10.));

        for (int i = 0; i < total_num_fields; i++) {
            if (i % dot_freq == 0) {
                std::cout << ".";
                std::cout.flush();
            }
            outfile.write((const char *)stacked[i],total_num_lines*sizeof(double));
        }

        std::cout << std::endl;
    }
}

//====================================================
// doReweightMassDistribution:
// Takes in a mass binned catalog and smoothes the distribution relative
// to the specified column (x).
// data:
// num_lines:
//          number of lines in the catalog.
// num_fields:
//          number of fields for each line (halo).
//====================================================
int doReweightMassDistribution (double ** &data, int num_lines, int num_fields, int x) {

    // object to hold halo fields, as well as allowing sorting by x and y fields
    struct halo {
        std::vector<double> fields;
        double x;
        double y;
        int rankx;
        int ranky;
    };

    // sorting comparator for rank ordering by field y
    struct sort_comp_y {
        bool operator()(const halo &left, const halo &right) {
            return (left.y < right.y);
        }
    };

    // sorting comparator for rank ordering by field x
    struct sort_comp_x {
        bool operator()(const halo &left, const halo &right) {
            return (left.x < right.x);
        }
    };

    // container for binning all the data
    std::vector<std::vector<std::vector<halo>>> bins;

    // container for sorting all the raw data
    std::vector<halo> vdata;

    for (int i = 0; i < num_lines; i++) {
        halo h;
        h.x = data[x][i];
        h.y = data[10][i];  // y is mvir for now
        for (int j = 0; j < num_fields; j++) {
            h.fields.push_back(data[j][i]); // push all halos fields to that halo
        }
        vdata.push_back(h);
    }

    std::sort(vdata.begin(),vdata.end(),sort_comp_y());

    for (int i = 0; i < vdata.size(); i++) {
        vdata[i].ranky = i;
    }

    std::sort(vdata.begin(),vdata.end(),sort_comp_x());

    for (int i = 0; i < vdata.size(); i++) {
        vdata[i].rankx = i;
    }

    // determine number of bins dynamically
    int num_bins_x = (int)floor(3*pow(2,log10(vdata.size())-2));
    double bin_size_x = 100./(double)num_bins_x;
    std::cout << "Num x bins determined: " << num_bins_x << std::endl;

    int num_bins_y = 2*(int)floor(3*pow(2,log10(vdata.size()/(double)num_bins_x)-2));
    double bin_size_y = 100./(double)num_bins_y;
    std::cout << "Num y bins determined: " << num_bins_y << std::endl;

    // create bin vectors
    for (int i = 0; i < num_bins_x; i++) {

        std::vector<std::vector<halo>> xbin;
        bins.push_back(xbin);

        for (int j = 0; j < num_bins_y; j++) {

            std::vector<halo> ybin;
            bins[i].push_back(ybin);
        }
    }

    // calculate percentile rank of current position in vdata and
    // use to compute index for binning
    for (int i = 0; i < vdata.size(); i++) {

        int index_x = (int) floor ((vdata[i].rankx/(double)vdata.size())*100. / bin_size_x);
        if (index_x < 0) index_x = 0;
        else if (index_x >= num_bins_x) index_x = num_bins_x - 1;

        int index_y = (int) floor ((vdata[i].ranky/(double)vdata.size())*100. / bin_size_y);
        if (index_y < 0) index_y = 0;
        else if (index_y >= num_bins_y) index_y = num_bins_y - 1;

        bins[index_x][index_y].push_back(vdata[i]);
    }

    // now we need to go through each ybin and randomly remove entries until we reach global
    // minimum entry count for each corresponding ybin

    // let's find them minimum counts for each ybin accross all corresponding ybins
    std::vector<double> min_count;
    
    for (int i = 0; i < num_bins_y; i++) {

        int min = bins[0][i].size();

        for (int j = 0; j < num_bins_x; j++) {

            if (bins[j][i].size() < min) {
                min = bins[j][i].size();
            }
        }

        min_count.push_back(min);
    }

    srand(time(NULL));

    // now lets equilibrate all ybins to have same count as corresponding min_count
    for (int i = 0; i < num_bins_x; i++) {
        for (int j = 0; j < num_bins_y; j++) {

            while (bins[i][j].size() > min_count[j]) {

                int remove_index = rand() % bins[i][j].size();
                bins[i][j].erase(bins[i][j].begin()+remove_index);
            }
        }
    }

    // now we need to output our new subset catalog
    int num_halos = 0;
    for (int i = 0; i < num_bins_x; i++) {
        for (int j = 0; j < num_bins_y; j++) {
            num_halos += bins[i][j].size();
        }
    }

    double ** subset = new double * [num_fields];
    for (int i = 0; i < num_fields; i++) subset[i] = new double [num_halos];

    int halo_count = 0;

    for (int i = 0; i < num_bins_x; i++) {
        for (int j = 0; j < num_bins_y; j++) {
            for (int k = 0; k < bins[i][j].size(); k++) {
                for (int l = 0; l < num_fields; l++) {
                    subset[l][halo_count] = bins[i][j][k].fields[l];
                }
                halo_count++;
            }
        }
    }

    writeOutputFileBinary ("", subset, num_halos, num_fields);

    free(subset,num_fields);

    return 0;

    // the below code implements a median forcing approach, and was replaced by
    // the above distribution forcing approach.

    /*

    // ADJUSTABLE PARAMETERS=============
    int miss_limit = 100;
    double min_retainment_fraction = 0.9;
    //===================================

    // container for binnning all the data
    std::vector<std::vector<halo>> bins;

    // keeps track of the worst retainment fraction out of all bins
    // after reweighting procedure has completed a full iteration.
    // initialize to 0 so that we can enter while loop
    double worst_retainment = 0;


    // iterate the full reweighting procedure until minimal changes are made
    while (worst_retainment < min_retainment_fraction) {

        std::cout << "Beginning of Iteration" << std::endl;

        // need to reset for next iteration
        worst_retainment = 1;

        // container for sorting all the raw data
        std::vector<halo> vdata;

        // will be true if this is not the first iteration.  if one iteration
        // has already happened, we should copy halo data from our previous
        // binning container rather than in raw catalog data 
        if (bins.size()) {
            for (int i = 0; i < bins.size(); i++) {
                for (int j = 0; j < bins[i].size(); j++) {
                    vdata.push_back(bins[i][j]);
                }
            }
        }

        else {
            for (int i = 0; i < num_lines; i++) {
                halo h;
                h.x = data[x][i];
                h.y = data[10][i];  // y is mvir for now
                for (int j = 0; j < num_fields; j++) {
                    h.fields.push_back(data[j][i]); // push all halos fields to that halo
                }
                vdata.push_back(h);
            }
        }

        bins.clear();


        std::sort(vdata.begin(),vdata.end(),sort_comp_x());

        // determine number of bins dynamically
        int num_bins = (int)floor(3*pow(2,log10(num_lines)-2));
        double bin_size = 100./(double)num_bins;
        std::cout << "Num bins determined: " << num_bins << std::endl;

        // create bin vectors
        for (int i = 0; i < num_bins; i++) {
            std::vector<halo> bin;
            bins.push_back(bin);
        }

        // calculate percentile rank of current position in vdata and
        // use to compute index for binning
        for (int i = 0; i < vdata.size(); i++) {
            int index = (int) floor (((double)i/(double)vdata.size())*100. / bin_size);
            if (index < 0) index = 0;
            else if (index >= num_bins) index = num_bins - 1;
            bins[index].push_back(vdata[i]);
        }

        // now go through and compute average median mvir
        double median = 0;
        for (int i = 0; i < num_bins; i++) {
            std::sort(bins[i].begin(),bins[i].end(),sort_comp_y());
            int s = bins[i].size();
            if (s % 2) median += bins[i][(s-1)/2].y;
            else median += 0.5*(bins[i][(s-2)/2].y + bins[i][s/2].y);
        }

        median /= (double)num_bins;

        std::cout << "Median Mvir on first pass is: " << log10(median) << std::endl;

        srand(time(NULL));

        // now go through again and iterate on each bin, forcing it towards the avg median
        for (int i = 0; i < num_bins; i++) {

            double median_i = 0;
            int s = bins[i].size();
            if (s % 2) median_i = bins[i][(s-1)/2].y;
            else median_i = 0.5*(bins[i][(s-2)/2].y + bins[i][s/2].y);

            bool forcing_down = median_i > median;
            int miss_count = 0;
            int init_size = s;

            while (miss_count < miss_limit) {

                int test_index = rand() % bins[i].size();

                if ((forcing_down  && (bins[i][test_index].y > median_i)) ||
                    (!forcing_down && (bins[i][test_index].y < median_i)))
                {

                    //if (bins[i].size()-1 < (int)floor(min_retainment_fraction*init_size)) {
                    //    std::cout << "bin " << i << " dropped below minimum retainment fraction. quitting..." << std::endl;
                    //    return -1;
                    //}

                    miss_count = 0;

                    halo h = bins[i][test_index];

                    bins[i].erase(bins[i].begin()+test_index);

                    // check if the bin is empty.  this can happen if the min retainment threshold
                    // is set to zero.
                    if (bins[i].size() == 0) {
                        continue;
                    }

                    s = bins[i].size();
                    if (s % 2) median_i = bins[i][(s-1)/2].y;
                    else median_i = 0.5*(bins[i][(s-2)/2].y + bins[i][s/2].y);
                    

                    // we've over shot
                    if ((   forcing_down && (median_i < median) ) ||
                        (  !forcing_down && (median_i > median) ))
                    {
                        // let's randomly insert halo that over-shot median so that we have
                        // equal over shoots and under shoots
                        if (rand() % 2) bins[i].push_back(h);

                        break;
                    }
                }
                else miss_count++;
            }

            std::cout << "median_" << i << ": " << log10(median_i) << ", retainment fraction: ";
            std::cout << ((double)bins[i].size())/(double)init_size << std::endl;
            if ((double)bins[i].size()/(double)init_size < worst_retainment) {
                worst_retainment = (double)bins[i].size()/(double)init_size;
            }
        }
    }

    // now we need to output our new subset catalog
    int num_halos = 0;
    for (int i = 0; i < bins.size(); i++) {
        num_halos += bins[i].size();
    }

    double ** subset = new double * [num_fields];
    for (int i = 0; i < num_fields; i++) subset[i] = new double [num_halos];

    int halo_count = 0;

    for (int i = 0; i < bins.size(); i++) {
        for (int j = 0; j < bins[i].size(); j++) {
            for (int k = 0; k < num_fields; k++) {
                subset[k][halo_count] = bins[i][j].fields[k];
            }
            halo_count++;
        }
    }


    writeOutputFileBinary ("", subset, num_halos, num_fields);

    free(subset,num_fields);

    return 0; */
}

int doReorderFields (double ** &data, int num_lines, int num_fields) {

    std::ofstream outfile;
    
    setOutputFileName();

    outfile.open(output.c_str(),std::ofstream::out);

    std::cout << "Writing output file \"" << output << "\"" << std::endl;
    std::cout << "num_lines: " << num_lines << std::endl;
    std::cout << "num_fields: " << num_fields << std::endl;

    outfile << "# NUM_LINES " << num_lines << std::endl;
    outfile << "# NUM_FIELDS " << num_fields << std::endl;

    int dot_freq = int(ceil(num_fields/10.));

    for (int i = 0; i < num_fields; i++) {
        if (i % dot_freq == 0) {
            std::cout << ".";
            std::cout.flush();
        }
        // take these columns and move to end.  These were fields Peter added to catalog after I had started using them,
        // so keeping them inserted in middle of fields would require me to redo all my column reference, which would be a pain
        // in the but.
        if (i == 34 ||
            i == 35 ||
            i == 36 ||
            i == 75   ) {
            continue;
        }
        outfile.write((const char *)data[i],num_lines*sizeof(double));
    }

    outfile.write((const char *)data[34],num_lines*sizeof(double));
    outfile.write((const char *)data[35],num_lines*sizeof(double));
    outfile.write((const char *)data[36],num_lines*sizeof(double));
    outfile.write((const char *)data[75],num_lines*sizeof(double));

    std::cout << " complete." << std::endl;

    outfile.close();

    return 0;
}

int doReorderFields2 (double ** &data, int num_lines, int num_fields) {

    std::ofstream outfile;
    
    setOutputFileName();

    outfile.open(output.c_str(),std::ofstream::out);

    std::cout << "Writing output file \"" << output << "\"" << std::endl;
    std::cout << "num_lines: " << num_lines << std::endl;
    std::cout << "num_fields: " << num_fields << std::endl;

    outfile << "# NUM_LINES " << num_lines << std::endl;
    outfile << "# NUM_FIELDS " << num_fields << std::endl;

    int dot_freq = int(ceil(num_fields/10.));

    for (int i = 0; i < num_fields; i++) {
        if (i % dot_freq == 0) {
            std::cout << ".";
            std::cout.flush();
        }
        // take these columns and move to end.  These were fields Peter added to catalog after I had started using them,
        // so keeping them inserted in middle of fields would require me to redo all my column reference, which would be a pain
        // in the but.
        if (i == 79 ||
            i == 80 ||
            i == 81 ||
            i == 82   ) {
            continue;
        }
        outfile.write((const char *)data[i],num_lines*sizeof(double));
    }

    outfile.write((const char *)data[79],num_lines*sizeof(double));
    outfile.write((const char *)data[80],num_lines*sizeof(double));
    outfile.write((const char *)data[81],num_lines*sizeof(double));
    outfile.write((const char *)data[82],num_lines*sizeof(double));

    std::cout << " complete." << std::endl;

    outfile.close();

    return 0;
}

// this is the standard form called from command line
int doCatalogMerge (double ** &data, int num_lines, int num_fields, std::string mergeTarget) {

    double ** data2;
    int num_lines2=0, num_fields2=0;

    // we've already read in one file.  now we need to read in the merge target.
    if (isBinaryExt (mergeTarget)) {
        if (readBinaryFile(mergeTarget,num_lines2,num_fields2,data2)) {
            return -1;
        }
    }
    else if (isBinary2Ext (mergeTarget)) {
        std::cout << "Warning: mergeTarget is binary2 file.  This is not yet supported." << std::endl;
        if (readBinary2File(mergeTarget,num_lines2,num_fields,data2)) {
            return -1;
        }
    }
    else {
        if (readFile(mergeTarget,num_lines2,num_fields,data2)) {
            return -1;
        }
    }

    doCatalogMerge (data, num_lines, num_fields, data2, num_lines2, num_fields2);

    return 0;
}

// merge data to target catalog, and save in specified output file. Right now we're assuming data is treeAnalysis output
// (z=0 from merger tree data) and data2 is rockstar augmented halo catalog (mergeTarget)
int doCatalogMerge (double ** &data, int num_lines, int num_fields, double ** &data2, int num_lines2, int num_fields2) {

    if (num_lines != num_lines2) {
        std::cout << "Error: number of lines in mergeTarget not equal to input file." << std::endl;
        return -1;
    }

    struct index {
        int i = -1;
        long int id = -1;
        int merge_i = -1;
    };

    // preserve ordering of merge target
    std::vector<index> dataVec;
    std::vector<index> data2Vec;

    for (int i = 0; i < num_lines; i++) {
        index in;
        in.i = i;
        in.id = data[BP_CAT.id][i];
        dataVec.push_back(in);
    }

    for (int i = 0; i < num_lines2; i++) {
        index in;
        in.i = i;
        in.id = data2[BP_CAT.id][i];
        data2Vec.push_back(in);
    }

    std::sort(dataVec.begin(),dataVec.end(),[](index left, index right){return left.id < right.id;});
    std::sort(data2Vec.begin(),data2Vec.end(),[](index left, index right){return left.id < right.id;});

    for (int i = 0; i < num_lines; i++) {
        if (dataVec[i].id != data2Vec[i].id) {
            std::cout << "ID mismatch at index " << i << std::endl;
        }
        else {
            data2Vec[i].merge_i = dataVec[i].i;
        }
    }

    std::sort(data2Vec.begin(),data2Vec.end(),[](index left, index right){return left.i < right.i;});

    std::ofstream of;

    setOutputFileName();

    of.open(output.c_str(),std::ofstream::out);

    int num_fields_tot = num_fields2 + 14;

    std::cout << "Writing output file \"" << output << "\"" << std::endl;
    std::cout << "num_lines: " << num_lines2 << std::endl;
    std::cout << "num_fields: " << num_fields_tot << std::endl;

    of << "# NUM_LINES " << num_lines2 << std::endl;
    of << "# NUM_FIELDS " << num_fields_tot << std::endl;

    int dot_freq = int(ceil(num_fields_tot/10.));

    // first output all fields in merge target
    for (int i = 0; i < num_fields2; i++) {
        if (i % dot_freq == 0) {
            std::cout << ".";
            std::cout.flush();
        }
        of.write((const char *)data2[i],num_lines2*sizeof(double));
    }

    // now collect new fields one at a time
    for (int i = 0; i < num_lines2; i++) {
        of.write((const char *)&data[BP_TREE.mpeak][data2Vec[i].merge_i],sizeof(double));
        if (i != data2Vec[i].i) std::cout << "index mismatch" << std::endl;
    }

    for (int i = 0; i < num_lines2; i++) of.write((const char *)&data[BP_TREE.a_peak][data2Vec[i].merge_i],sizeof(double));
    for (int i = 0; i < num_lines2; i++) of.write((const char *)&data[BP_TREE.mdot_dyn][data2Vec[i].merge_i],sizeof(double));
    for (int i = 0; i < num_lines2; i++) of.write((const char *)&data[BP_TREE.npt][data2Vec[i].merge_i],sizeof(double));
    for (int i = 0; i < num_lines2; i++) of.write((const char *)&data[BP_TREE.high_tf][data2Vec[i].merge_i],sizeof(double));
    for (int i = 0; i < num_lines2; i++) of.write((const char *)&data[BP_TREE.maxtf_peak][data2Vec[i].merge_i],sizeof(double));
    for (int i = 0; i < num_lines2; i++) of.write((const char *)&data[BP_TREE.a_maxtf_peak][data2Vec[i].merge_i],sizeof(double));
    for (int i = 0; i < num_lines2; i++) of.write((const char *)&data[BP_TREE.min_bsr_peak][data2Vec[i].merge_i],sizeof(double));
    for (int i = 0; i < num_lines2; i++) of.write((const char *)&data[BP_TREE.a_min_bsr_peak][data2Vec[i].merge_i],sizeof(double));
    for (int i = 0; i < num_lines2; i++) of.write((const char *)&data[BP_TREE.min_bsr_lmm][data2Vec[i].merge_i],sizeof(double));
    for (int i = 0; i < num_lines2; i++) of.write((const char *)&data[BP_TREE.a_min_bsr_lmm][data2Vec[i].merge_i],sizeof(double));
    for (int i = 0; i < num_lines2; i++) of.write((const char *)&data[BP_TREE.spin_peak][data2Vec[i].merge_i],sizeof(double));
    for (int i = 0; i < num_lines2; i++) of.write((const char *)&data[BP_TREE.tu_peak][data2Vec[i].merge_i],sizeof(double));
    for (int i = 0; i < num_lines2; i++) of.write((const char *)&data[BP_TREE.mass_lost_as_sh][data2Vec[i].merge_i],sizeof(double));

    std::cout << " complete." << std::endl;

    of.close();

    return 0;
}
