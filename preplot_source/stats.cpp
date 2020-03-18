#include "stats.h"

// these are used to set the min/max values of the binning parameter
// so that it can be used while computing veffective of each bin.
double doBin_data_min = 0;
double doBin_data_max = 0;

/*****************
 * getIndex
 * ---------------------------------------------------------------------------------------
 * Finds and returns the middle position (percentile) of target value in the sorted
 * vector vdata.  As of now only called from median function when using y data percentiles.
 * ---------------------------------------------------------------------------------------
 * vdata        vector containing sorted y value data
 * target       value for which to look up position in vdata
 * returns      percentile value of target
 *****************/
double getIndex (std::vector<double> &vdata, double target) {
    
    if (!vdata.size()) return -1;

    //std::cout << "getIndex called with target (" << target << ") and vdata range: " << vdata[0] << " - " << vdata[vdata.size()-1] << std::endl;

    int lbound = 0, rbound = vdata.size();
    int pivot = int(floor(vdata.size()/2.));

    if (target == vdata[0]) return 0.;
    if (target == vdata[vdata.size()-1]) return 100.;

    while (true) {

        if (target == vdata[pivot]) return (100.*pivot)/double(vdata.size());
        else if (target < vdata[pivot]) rbound = pivot;
        else lbound = pivot;

        if ((rbound-lbound) == 1) {
            //std::cerr << "getIndex: target (" << target << ") not found. Closest (lbound,rbound): (" << vdata[lbound] << ", " << vdata[rbound] << ")" << std::endl;
            //return -2;
            return (100*(lbound+0.5))/double(vdata.size());
        }

        pivot = int(floor((rbound+lbound)/2.));
    }
}

double pdf_integrate (double lbound, double rbound, std::string pdfFileName) {

    // let's first read the file in
    double ** data;
    int num_fields=0, num_lines=0;

    // read in pdf.  first column should be x values, 2nd column normalized pdf.
    // 3rd column is total bin width of given bin -- use this to compute area
    toggle_stream(std::cout);
    if (readFile(pdfFileName,num_lines,num_fields,data)) {
        return -1;
    }
    toggle_stream(std::cout);

    double result = 0;
    double total = 0;
    bool percentile = false;

    if (rbound < lbound) {
        percentile = true;
        std::cout << "Determining value at percentile: " << lbound << std::endl;
    }

    // now we need to step from lbound to rbound and add up area
    for (int i = 0; i < num_lines; i++) {
        double area = data[2][i]*data[1][i];
        total += area;
        if (data[0][i] > lbound && data[0][i] < rbound) {
            result += area;
        }
    }

    double pvalue = 0;

    if (percentile) {
        double total2 = 0;
        for (int i = 0; i < num_lines; i++) {
            if (lbound == 100.) {
                pvalue = data[0][num_lines-1];
                break;
            }
            double area = data[2][i]*data[1][i];
            total2 += area;
            if (100.*total2/total >= lbound) {
                pvalue = data[0][i];
                break;
            }
        }

        pvalue = pow(10.,pvalue);

        std::cout << "total integral: " << total << std::endl;
        std::cout << "density at " << lbound << "%: " << pvalue << std::endl;
    }

    else {
        std::cout << "total integral: " << total << std::endl;
        std::cout << "fraction with bounds: (" << lbound << ", " << rbound << "): " << result/total << std::endl;
    }

    free (data, num_fields);

    if (percentile) return pvalue;

    return result/total;
}

// overloaded function calls pdf_integrate(double,double,string) once appropriate
// pdf file name has been determined using the binning variable x
double pdf_integrate (double lbound, double rbound, int x) {

    std::string pdfFileName = "../bolshoi_plank/DougsCICs/BP_";
    
    std::vector<std::string> z = {"_z0_","_z0.5_","_z1_","_z2_"};
    std::vector<std::string> snap = {"0214","0170","0136","0103"};
    
    // let's check beginning of intput file name to see if we can get the redshift
    std::string::size_type n;
    int i = 0;
    for (i = 0; i < z.size(); i++) {
        n = fileName.find(z[i]);
        if (n == std::string::npos) continue;
        else break;
    }
    if (n == std::string::npos) {
        std::cerr << "Error: compute_veff could not determine redshift of input." << std::endl;
        return -1;
    }
    
    pdfFileName += snap[i]+"_densities_1024_";

    switch (x) {
        case 75:    pdfFileName += "2"; break;
        case 76:    pdfFileName += "4"; break;
        case 77:    pdfFileName += "8"; break;
        case 78:    pdfFileName += "16"; break;
    }

    pdfFileName += "_pdf.dat";

    if (rbound > 100) std::cout << "Reading from " << pdfFileName << std::endl;

    return pdf_integrate (lbound,rbound,pdfFileName);
}

/*****************
 * compute_veff
 * ---------------------------------------------------------------------------------------
 * Computes the effective volume of a given density range by computing the fraction of 
 * voxels in the smoothed density cubes that occupy that density range.  Need to deduce
 * the redshift and density range from input file name and output file name.
 * ---------------------------------------------------------------------------------------
 * veffective       the effective volume of the current density range in use
 * returns          success/fail.  veffective updated through reference.
 *****************/
int compute_veff (double & veffective) {

    // we need to know some stuff in order to do this computation.
    // first thing to figure out is what file should we look at to find
    // the appropriate pdf.

    std::string pdfFileName = "../bolshoi_plank/DougsCICs/BP_";

    std::vector<std::string> z = {"_z0_","_z0.5_","_z1_","_z2_"};
    std::vector<std::string> snap = {"0214","0170","0136","0103"};

    // let's check beginning of intput file name to see if we can get the redshift
    std::string::size_type n;
    int i = 0;
    for (i = 0; i < z.size(); i++) {
        n = fileName.find(z[i]);
        if (n == std::string::npos) continue;
        else break;
    }
    if (n == std::string::npos) {
        std::cerr << "Error: compute_veff could not determine redshift of input." << std::endl;
        return -1;
    }
    
    pdfFileName += snap[i]+"_densities_1024_";

    // now we need to determine the density range to use.
    // let's define the possible ranges here, then deduce which bin we should select based
    // on the output file number of this current execution.
    //std::vector<double> bin_edge {-0.5,-0.175,-0.05,0.05,0.175,0.5};
    //std::vector<double> bin_edge {30,45,55,65,80};
    std::vector<double> bin_edge {10,30,50,70,90};
    bool percentiles = true;

    n = output.find("_hist1D");
    if (n == std::string::npos) {
        std::cerr << "Error: compute_veff could not determine output file number." << std::endl;
        return -1;
    }

    // Note: we are assuming output file numbers 0-9 only, two digit output file numbers not handled.
    int ofnum = output[n-1]-'0';//std::strtol(fileName[n-1].c_str(),NULL,10);

    std::cout << "Detected redshift " << z[i] << ", ofnum " << ofnum << std::endl;

    double bin_edge_l, bin_edge_r;

    // the doBin_data_min/max come from the doBin function, and are set during
    // the inital binning phase, when all of the prebinned data is available.
    // this allows us to set bounds on what the minimum/maximum halo density sampling
    // from the density pdf should be.
    if (ofnum == 0) {
        //bin_edge_l = -std::numeric_limits<double>::infinity();
        if (percentiles) bin_edge_l = 0.;
        else bin_edge_l = doBin_data_min;
        bin_edge_r = bin_edge[ofnum];
    }
    else if (ofnum == bin_edge.size()) {
        bin_edge_l = bin_edge[ofnum-1];
        //bin_edge_r = std::numeric_limits<double>::infinity();
        if (percentiles) bin_edge_r = 100.;
        else bin_edge_r = doBin_data_max;
    }
    // this is the last output file, which is for all densities
    else if (ofnum == bin_edge.size()+1) {
        //bin_edge_l = -std::numeric_limits<double>::infinity();
        //bin_edge_r = std::numeric_limits<double>::infinity();
        if (percentiles) {
            bin_edge_l = 0.;
            bin_edge_r = 100.;
        }
        else {
            bin_edge_l = doBin_data_min;
            bin_edge_r = doBin_data_max;
        }
    }
    else {
        bin_edge_l = bin_edge[ofnum-1];
        bin_edge_r = bin_edge[ofnum];
    }

    // now we have our density range and we know what redshift we are at.
    // just need to compute veff for the scale we are interested in.

    // let's start with scale 4mpc/h. we could always infer the scale from the output
    // file name in the future as well, so long as we include it in the output file name.
    pdfFileName += "4_pdf.dat";

    // now pass to helper function to do the integratation.  this will return a fraction
    // (integral divided by total integral, though total should be 1 if normalized anyways).
    if (percentiles) veffective = (bin_edge_r-bin_edge_l)/100.;
    else veffective = pdf_integrate (bin_edge_l, bin_edge_r, pdfFileName);

    veffective *= (250*250*250.);

    std::cout << "Veff: " << veffective << ", Vtotal: " << (250*250*250.) << std::endl;

    return 0;
}

int doAvg(double ** &data, int num_lines, int x, int y, int num_bins) {

	double minX, maxX;
	double bin_size;
	int index;
	
	for (int i = 0; i < num_lines; i++) {
		if (i == 0) {
			minX = data[x][i];
			maxX = data[x][i];
		}
		if (data[x][i] < minX) minX = data[x][i];
		if (data[x][i] > maxX) maxX = data[x][i];
	}

	std::cout << "minX: " << minX << std::endl << "maxX: " << maxX << std::endl;

	bin_size = (maxX-minX)/(double)num_bins;

	double bins[num_bins][3];

	for (int i = 0; i < num_bins; i++) {
		bins[i][0] = minX + i * bin_size;
		bins[i][1] = 0.0;
		bins[i][2] = 0.0;
	}

	for (int i = 0; i < num_lines; i++) {
		index = (int) floor ((data[x][i] - minX) / bin_size);
		if (index < 0) index = 0;
		else if (index >= num_bins) index = num_bins - 1;
		if (std::isfinite(data[y][i])) {
			bins[index][1] += 1.0;
			bins[index][2] += data[y][i];
		}
	}

	for (int i = 0; i < num_bins; i++) {
		if (bins[i][2] == 0.0) {
			if (!i) {
				bins[i][1] = bins[i-1][1];
			}
		}
		else {
			bins[i][1] = bins[i][2] / bins[i][1];
		}
	}

	// now write to output file
	std::ofstream of_hist;

	setOutputFileName(x,y);

	of_hist.open(output.c_str(),std::ofstream::out);

	for (int i = 0; i < num_bins; i++) {
		of_hist << bins[i][0] << " " << bins[i][1] << std::endl;
	}
	
	of_hist.close();

	return 0;
}

int doMedian(double ** &data, int num_lines, int x, int y, int num_bins, bool autobin, int auto_tnum, bool use_percentiles_x, bool use_percentiles_y) {

	double minX, maxX;
	double bin_size;
	int index;

    // NOTE: this is a hack to allow us to compute some custom quantities that I'm to lazy to write a proper routine to implement.
    // we will just overwrite the provided y field (could be anything, so chose somewhat carefully at input) with our custom quantity for now.
    //std::cout << "Using hack to compute rs_jiang values..." << std::endl;
    //for (int i = 0; i < num_lines; i++) {
    //    data[y][i] = data[BP_CAT.rs][i]*pow((data[BP_CAT.rvir][i]/data[BP_CAT.rs][i])/7.,0.4);
    //}

    typedef struct {
        std::vector<double> elements;
        double width = 0;
        double center = 0;
        double median = 0;
        double disp80 = 0;
        double disp20 = 0;
        double mad = 0;
        double ci_95_lo = 0;
        double ci_95_hi = 0;
    } bin;
    
    typedef struct {
        double x = 0;
        double y = 1;
        double bsr = 0;
    } element;

    // check if we need to set bin size
    if (!num_bins) {
        num_bins = (int)floor(3*pow(2,log10(num_lines)-2));//(int)ceil(log(num_lines));
        std::cout << "Num bins determined: " << num_bins << std::endl;
    }

	for (int i = 0; i < num_lines; i++) {
		if (i == 0) {
			minX = data[x][i];
			maxX = data[x][i];
		}
		if (data[x][i] < minX) minX = data[x][i];
		if (data[x][i] > maxX) maxX = data[x][i];
	}

	std::cout << "minX: " << minX << std::endl << "maxX: " << maxX << std::endl;

    // if using percentiles for the x-axis, each bin will by definition have
    // the same number of halos.  ajust min/max for initializing bin centers etc.
    if (use_percentiles_x) {
        bin_size = 100./(double)num_bins;
        minX = 0.;
        maxX = 100.;
    }

	else bin_size = (maxX-minX)/(double)num_bins;

    std::vector<bin> vbins (num_bins);

    for (int i = 0; i < num_bins; i++) {
        vbins[i].width = bin_size/2.;
        vbins[i].center = minX + i * bin_size + bin_size/2.;
    }

    //std::vector<std::pair<double,double>> vdata_rank_x;
    std::vector<element> vdata_rank_x;
    std::vector<double> vdata_rank_y;

    if (use_percentiles_x || use_percentiles_y) {

        //struct sort_comp {
        //    bool operator()(const std::pair<double,double> &left, const std::pair<double,double> &right) {
        //        return left.first < right.first;
        //    }
        //};
        
        struct sort_comp {
            bool operator()(const element &left, const element &right) {
                return left.x < right.x;
            }
        };

        for (int i = 0; i < num_lines; i++) {

            //vdata_rank_x.push_back(std::make_pair(data[x][i],data[y][i]));
            element e;
            e.x = data[x][i];
            e.y = data[y][i];
            e.bsr = data[80][i];
            vdata_rank_x.push_back(e);
            vdata_rank_y.push_back(data[y][i]);
        }

        std::sort(vdata_rank_x.begin(),vdata_rank_x.end(),sort_comp());
        std::sort(vdata_rank_y.begin(),vdata_rank_y.end());
    }

    if (use_percentiles_x) {

        //for (int i = 0; i < num_lines; i++) {

        int i = 0;

        for (std::vector<element>::iterator it = vdata_rank_x.begin(); it != vdata_rank_x.end(); it++) {

            // exclude stripped halos for this run
            if (vdata_rank_x[i].bsr > 0.98) {
                i++;
                continue;
            }

            // calculate percentile rank of current position in vdata and
            // use to compute index for binning
            index = (int) floor (((double)i/(double)num_lines)*100. / bin_size);
            if (index < 0) index = 0;
            else if (index >= num_bins) index = num_bins - 1;
            if (std::isfinite(vdata_rank_x[i].y)) {
                vbins[index].elements.push_back(vdata_rank_x[i].y);
            }

            i++;
        }
    }

    else {

	    for (int i = 0; i < num_lines; i++) {
	    	index = (int) floor ((data[x][i] - minX) / bin_size);
	    	if (index < 0) index = 0;
	    	else if (index >= num_bins) index = num_bins - 1;
	    	if (std::isfinite(data[y][i])) {
                vbins[index].elements.push_back(data[y][i]);
	    	}
	    }
    }

    if (autobin) {

        if (auto_tnum == -1) {
            auto_tnum = num_bins;
        }

        std::cout << "Auto rebinning..." << std::endl;
        std::cout.flush();

        // rebin auto
        bool merged = false;
        std::vector<bin>::iterator it = vbins.begin();
        while (true) {

            // reached end of bins vector
            if (it == vbins.end()) {
                if (merged) {
                    it = vbins.begin();
                    merged = false;
                }
                else break;
            }

            //std::cout << "center: " << it->center << "\t" << it->elements.size() << std::endl;
            //std::cout.flush();

            // check if num elements is less than threshold
            if (it->elements.size() < auto_tnum) {
                // merge unless right neighbor is smaller -- if so, skip
                if ((it+1) != vbins.end()) {
                    if ((it+1)->elements.size() < it->elements.size()) {
                        it++;
                        continue;
                    }
                }

                merged = true;
                int j = 0; // merge direction
                // merge to smallest neighbor
                if (it == vbins.begin()) j = 1;
                else if ((it+1) == vbins.end()) j = -1;
                else if ((it-1)->elements.size() < (it+1)->elements.size()) j = -1;
                else j = 1;
                //std::string arrow;
                //if (j==1) arrow = " -> ";
                //else arrow = " <- ";
                //std::cout << "merging centers (" << j << ")\t" << it->center << arrow << (it+j)->center;
                //std::cout << "\tcounts: (" << it->elements.size() << "," << (it+j)->elements.size() << ")";

                (it+j)->elements.insert((it+j)->elements.end(),it->elements.begin(),it->elements.end());
                (it+j)->center = it->center + j * (it+j)->width;
                (it+j)->width += it->width;
                // remove from bin list
                vbins.erase(it);
                //if (it == vbins.end()) std::cout << "\t after erasing, the end." << std::endl;
                //else std::cout << "\t after erasing, count is " << it->elements.size() << std::endl;
                //std::cout.flush();
                continue;
            }
            it++;
        }

        num_bins = vbins.size();

        std::cout << "Number of bins determined from auto-rebinning: " << num_bins << std::endl;
    }

    for (std::vector<bin>::iterator it = vbins.begin(); it != vbins.end(); it++) {
        if (it->elements.size()) {
            std::sort(it->elements.begin(),it->elements.end());
            //it->median = it->elements[(int)floor(it->elements.size()/2.)];
            int s = it->elements.size();
            if (s % 2)  it->median = it->elements[(s-1)/2];
            else        it->median = 0.5*(it->elements[(s-2)/2]+it->elements[s/2]);
            it->disp20 = it->elements[(int)floor(s*0.2)];
            it->disp80 = it->elements[(int)ceil(s*0.8)];
        }
    }


    // compute median absolute deviation
    for (int i = 0; i < vbins.size(); i++) {
        std::vector<double> dev;
        if (vbins[i].elements.size()) {
            for (std::vector<double>::iterator it = vbins[i].elements.begin(); it < vbins[i].elements.end(); it++) {
                dev.push_back(fabs(*it-vbins[i].median));
            }
            std::sort(dev.begin(), dev.end());
            vbins[i].mad = dev[(int)floor(dev.size()*.5)];
        }
    }

    // compute confidence interval for median
    for (std::vector<bin>::iterator it = vbins.begin(); it != vbins.end(); it++) {
        // for low n, compute quartiles of binomial distribution
        //if (it->elements.size() < 100) {
        //    int k = 1+3; //pass
        //}
        //// use large n approximation
        //else {
        if (it->elements.size()) {
            int n = it->elements.size();
            double j = n*.5-1.96*sqrt(n*.25)+1;
            double k = n*.5+1.96*sqrt(n*.25)-1.;
            it->ci_95_lo = it->elements[(int)floor(j)];
            it->ci_95_hi = it->elements[(int)ceil(k)];
            double bs_hi, bs_lo;
            // not going to use bootstrapping to determine confidence interval for now, since the approximation above
            // produces nearly identical results
            //bootstrap(it->elements,0.95,bs_hi,bs_lo);
            //std::cout << "95% Confidence Interval approximation is (" << it->ci_95_lo << ", " << it->ci_95_hi << ")" << std::endl;
        }
        //}
    }

    // determine percentiles of y data relative all ydata in plot, rather than separately
    // in each bin
    if (use_percentiles_y) {
        for (std::vector<bin>::iterator it = vbins.begin(); it != vbins.end(); it++) {
            if (it->elements.size()) {
                it->median = getIndex(vdata_rank_y,it->median);
                it->disp20 = getIndex(vdata_rank_y,it->disp20);
                it->disp80 = getIndex(vdata_rank_y,it->disp80);
                // can't use percentiles for mad, unless we add median +- mad
                //it->mad = getIndex(vdata_rank_y,it->mad);
                it->ci_95_lo = getIndex(vdata_rank_y,it->ci_95_lo);
                it->ci_95_hi = getIndex(vdata_rank_y,it->ci_95_hi);
            }
        }
    }

	// now write to output file
	std::ofstream of_hist;

	setOutputFileName(x,y);

	of_hist.open(output.c_str(),std::ofstream::out);

	for (int i = 0; i < num_bins; i++) {
        of_hist << vbins[i].center << " " << vbins[i].median << " " << vbins[i].disp20 << " " << vbins[i].disp80 << " " << vbins[i].mad << " " << vbins[i].ci_95_lo << " " << vbins[i].ci_95_hi << " " << vbins[i].center-vbins[i].width << " " << vbins[i].center+vbins[i].width << " " << vbins[i].elements.size() << std::endl;
	}
	
	of_hist.close();

	return 0;

}

int doHist1D(double ** &data, int num_lines, int x, int num_bins, bool NORM, bool LOGZ, bool CDF, bool autobin, int auto_tnum, bool MASSFUNC, double bin_width) {

	// now lets begin with the actual histogram
	double minX, maxX;
	int index;
	
    bool doingBacksplashFraction = false;
    //bool doingBacksplashFraction = true;

    bool doing_zspace = false;

    bool notSet = 1;

    // object to hold bin properties
    struct bin {
        double pdf = 0;         // pdf value for this bin
        double cdf = 0;         // cdf value for this bin
        int count = 0;          // actual number of halos in this bin
        double width = 0;       // width (radius) of bin
        double center = 0;      // center of bin
        bin ** self = NULL;     // pointer to self pointer anchor
        bin ** left = NULL;     // pointer to bin to the left
        bin ** right = NULL;    // pointer to bin to the right
        bin ** lowest_n = NULL; // neighboring bin with lowest count
        int lowest_ncount = 0;  // count of neighboring with lowest count
    };

    std::cout << "Using hack to compute rs_jiang values..." << std::endl;
    for (int i = 0; i < num_lines; i++) {
        data[x][i] = log10(data[BP_CAT.rs][i]*pow((data[BP_CAT.rvir][i]/data[BP_CAT.rs][i])/7.,0.4));
    }

	for (int i = 0; i < num_lines; i++) {
        if (std::isfinite(data[x][i])) {
		    if (notSet) {
		    	minX = data[x][i];
		    	maxX = data[x][i];
                notSet = 0;
		    }
		    if (data[x][i] < minX) minX = data[x][i];
		    if (data[x][i] > maxX) maxX = data[x][i];
        }
	}

	std::cout << "minX: " << minX << std::endl << "maxX: " << maxX << std::endl;

    // check if we need to set bin width
    if (!num_bins) {
        // adjusting this to seem about right for pdfs.  makes sense that this should
        // be more bins for a given sample size that the median analysis
        num_bins = (int)floor(10*pow(3.,log10(num_lines)-2));//(int)ceil(log(num_lines));
    }

    if (bin_width == 0.0) bin_width = (maxX-minX)/(double)num_bins;
    else num_bins = ceil((maxX-minX)/bin_width);


    std::cout << "Num bins determined: " << num_bins << std::endl;

    // set initial bin size
	bin_width = (maxX-minX)/(double)num_bins;

    std::vector<bin> vbins (num_bins);      // vector containing all current bins
    std::vector<bin *> anchor (num_bins);   // vector containing pointers to current bins

    // NOTE: the idea with the anchor vector is to have a static array of pointers to the actual
    // bins.  Since the bins are being sorted repeatedly, their addresses will change, but we
    // still want to be able to point to their neighbors, etc.  To achieve this, every bin must
    // update its new address in the anchor vector after being sorted.  All subsequent calls
    // to neighboring bins, etc, are then routed through their respective anchor points before
    // being directed to the current location of the actual bins.

    // almm
    if (x == 15 || x == 71 || x == 72) {
        doing_zspace = true;
    }

    std::vector<double> vdata;

    if (doing_zspace) {
        for (int i = 0; i < num_lines; i++) {
            vdata.push_back(data[x][i]);
        }

        std::sort(vdata.begin(),vdata.end());

        int j = 0;
        for (int i = 0; i < vdata.size(); i++) {
            if (i > 0) {
                if (vdata[i] != vdata[i-1]) {
                    j++;
                    if (j >= vbins.size()) {
                        bin b;
                        vbins.push_back (b);
                        bin * bp;
                        anchor.push_back (bp);
                    }
                    // handle zero-bin width
                    if (j == 1) {
                        vbins[j-1].width = (vdata[i]-vbins[j-1].center)/2.;
                    }
                    else {
                        vbins[j-1].width = (vdata[i]-vbins[j-2].center)/4.;
                    }
                }
            }
            vbins[j].center = vdata[i];
            vbins[j].count++;

            // set bin width for final bin
            if (i+1 == vdata.size()) {
                vbins[j].width = (vdata[i] - vbins[j-1].width)/2.;

                // erase any remaining unused bins from bin vector
                if (j < vbins.size()-1) {
                    vbins.erase(vbins.begin()+j+1,vbins.end());
                }
            }
        }

        std::cout << "Total number of bins determined from z-space: " << vbins.size() << std::endl;
    }

    // establish inital bin sizes, locations, pointers
    for (int i = 0; i < vbins.size(); i++) {

        if (!doing_zspace) {
            vbins[i].width = bin_width/2.;
            vbins[i].center = minX + i * bin_width + bin_width/2.;
        }
        
        vbins[i].self = &(anchor[i]);
        anchor[i] = &(vbins[i]);
    }

    for (int i = 0; i < vbins.size(); i++) {
        if (i > 0) {
            vbins[i].left = vbins[i-1].self;
            //std::cout << "left bin: " << vbins[i].left->center << std::endl;
            //std::cout << "left bin: " << (*(vbins[i].self))->left->center << std::endl;
        }
        if (i < vbins.size()-1) {
            vbins[i].right = vbins[i+1].self;
        }
    }

	double bins[num_bins][2];
    double bins_bs[num_bins][5];

	for (int i = 0; i < num_bins; i++) {
		bins[i][0] = minX + i * bin_width;
		bins[i][1] = 0.0;
	}

    //000000000000000000000 for bs fraction only
    if (doingBacksplashFraction) {
        for (int i = 0; i < num_bins; i++) {
            bins_bs[i][0] = minX + i * bin_width;
            for (int j = 1; j < 5; j++) bins_bs[i][j] = 0.0;
        }
    }

    int outofbounds1 = 0, outofbounds2 = 0;

	for (int i = 0; i < num_lines; i++) {
        if (doing_zspace) break;
		index = (int) floor ((data[x][i] - minX) / bin_width);
		if (index < 0) {
            outofbounds1++;
            index = 0;
        }
		else if (index >= num_bins) {
            outofbounds2++;
            index = num_bins - 1;
        }
		bins[index][1] += 1.0;
        vbins[index].count++;

        //00000000 for bs fraction only
        if (doingBacksplashFraction) {
            if (data[80][i] > 0.8)
                bins_bs[index][1] += 1.0;
            else {
                bins_bs[index][2] += 1.0;
                if (data[80][i] <= 0.6)
                    bins_bs[index][3] += 1.0;
                if (data[80][i] <= 0.4)
                    bins_bs[index][4] += 1.0;
            }
        }
	}

    std::cout << "outofbounds < 0: " << outofbounds1 << ", outofbounds > : " << outofbounds2 << std::endl;

    // turn pdf into cumulative distribution function, normalized to 1 at the maximum.
    // Deprecated use of bins. Current implementation found after autobinning section.
    if (CDF) {
        for (int i = 1; i < num_bins; i++) bins[i][1] += bins[i-1][1];
        for (int i = 0; i < num_bins; i++) bins[i][1] /= (double)num_lines;
    }

    // now lets take care of autobinning if requested. this will merge low count bins together
    // to provide a pdf that adjusts in resolution with available statistics in a particular region
    if (autobin) {

        if (auto_tnum == -1) {
            auto_tnum = num_bins;
        }

        std::cout << "Autobinning..." << std::endl;
        std::cout << "Using bin threshold number: " << auto_tnum << std::endl;

        struct sort_comp {
            bool operator()(const bin &left, const bin &right) {
                return (left.count == right.count) ? left.lowest_ncount < right.lowest_ncount : left.count < right.count;
            }
        };

        double min_bin_width = vbins[0].width;

        // determine lowest neighbor counts
        for (int i = 0; i < vbins.size(); i++) {
            if (vbins[i].left == NULL) {
                vbins[i].lowest_n = vbins[i].right;
            }
            else if (vbins[i].right == NULL) {
                vbins[i].lowest_n = vbins[i].left;
            }
            else {
                vbins[i].lowest_n = (*(vbins[i].left))->count < (*(vbins[i].right))->count ? vbins[i].left : vbins[i].right;
            }
            if (vbins[i].width < min_bin_width) {
                min_bin_width = vbins[i].width;
            }
            vbins[i].lowest_ncount = (*(vbins[i].lowest_n))->count;
            vbins[i].pdf = vbins[i].count;
        }

        if (min_bin_width < bin_width) {
            for (int i = 0; i < vbins.size(); i++) {
                vbins[i].pdf /= (vbins[i].width/min_bin_width);
            }
        }

        //std::cout << "min_bin_width: " << min_bin_width << std::endl;

        //for (int i = 0; i < vbins.size(); i++) {
        //    std::cout << "Bin " << i << ": " << vbins[i].count << ", " << (*(vbins[i].lowest_n))->count << ", " << vbins[i].pdf << std::endl;
        //}

        int j = 0;

        // now we iterate through the bins vector and merge bins until done
        while (true) {

            //std::cout << "Iteration " << j++ << std::endl;

            //sort the bins vector
            std::sort(vbins.begin(),vbins.end(),sort_comp());

            // now lets update all the anchor pointers
            for (int i = 0; i < vbins.size(); i++) {
                *(vbins[i].self) = &(vbins[i]);
            }

            std::vector<bin>::iterator it = vbins.begin();

            // don't merge bins with zero count into bins with nonzero count (other way around is ok though).
            // this is to discourage outliers from heavily affecting other parts of histogram.
            // also enforce a maximum bin size (right now set to 10 times the initially determined bin_width).
            // so this would be a maximum of 10 original bins merged together.

            // while ((it->width > bin_width*5 && it->count == 0) ||
            //       ((*(it->lowest_n))->width > bin_width*5 && (*(it->lowest_n))->count == 0)) {
            while ( (it->count == 0) ||
                  ( (it->width > bin_width*5) || ((*(it->lowest_n))->width > bin_width*5) ) ) {

                //std::cout << "Can't merge large bin at " << it->center << " since count is zero." << std::endl;

                // go to next priority bin
                it++;

                // we have no mergeable bins remaining
                if (it == vbins.end()) break;
            }

            // break out of both while loops
            if (it == vbins.end()) break;

            // merge until lowest bin count is above auto_tnum.
            if (it->count > auto_tnum) {
                break;
            }

            //for (int i = 0; i < vbins.size(); i++) {
            //    std::cout << "Bin " << i << ": " << vbins[i].count << ", " << (*(vbins[i].lowest_n))->count << ", " << vbins[i].pdf << std::endl;
            //}

            bin * target = *(it->lowest_n);

            //std::cout << "count: " << it->count << "->" << target->count << std::endl;
            //std::cout << "center: " << it->center << "->" << target->center << std::endl;
            //std::cout << "width: " << it->width << "->" << target->width << std::endl;
            //std::cout << "pdf: " << it->pdf << "->" << target->pdf << std::endl;

            // merge lowest bin with lowest count neighbor.
            target->pdf = target->pdf * target->width;
            target->pdf += it->pdf * it->width;
            target->pdf /= target->width + it->width;

            target->count += it->count;

            if (target->center > it->center) {
                target->center -= it->width;
            }
            else {
                target->center += it->width;
            }

            target->width += it->width;

            //std::cout << "new count: " << target->count << std::endl;
            //std::cout << "new center: " << target->center << std::endl;
            //std::cout << "new width: " << target->width << std::endl;
            //std::cout << "new pdf: " << target->pdf << std::endl;

            // with only 2 bins left, there is no point in doing the below steps, since
            // only 1 bin will remain after this iteration anyways
            if (vbins.size() == 2) {
                vbins.erase(it);
                break;
            }

            // now need to nullify old low count bin and update
            // neighboring bin pointers and lowest_n/counts
            if (it->left) {

                // reassign target to left neighbor (to reuse the below code)
                target = (*(it->left));

                // skip bin that got merged
                target->right = it->right;

                // now need to update lowest_n for left neighbor
                if (target->left == NULL) {
                    target->lowest_n = target->right;
                }
                else if (target->right == NULL) {
                    target->lowest_n = target->left;
                }
                else {
                    target->lowest_n = (*(target->left))->count < (*(target->right))->count ? target->left : target->right;
                }
                target->lowest_ncount = (*(target->lowest_n))->count;

            }
            // same procedure as above for right neighbor
            if (it->right) {

                target = (*(it->right));

                target->left = it->left;

                // now need to update lowest_n for right neighbor
                if (target->left == NULL) {
                    target->lowest_n = target->right;
                }
                else if (target->right == NULL) {
                    target->lowest_n = target->left;
                }
                else {
                    target->lowest_n = (*(target->left))->count < (*(target->right))->count ? target->left : target->right;
                }
                target->lowest_ncount = (*(target->lowest_n))->count;
            }

            // nullify anchor pointer
            *(it->self) = NULL;

            // remove merged bin from bin vector
            vbins.erase(it);  
        }

        struct sort_comp_final {
            bool operator()(const bin &left, const bin &right) {
                return left.center < right.center;
            }
        };

        // now need to resort the final vector in order

        std::sort(vbins.begin(),vbins.end(),sort_comp_final());

        std::cout << "Number of bins remaining after autobinning: " << vbins.size() << std::endl;
    }

    // turn pdf into cumulative distribution function, normalized to 1 at the maximum.
    // I've moved this to after the autobinning code since cdfs would need to be updated
    // at this stage anyways (if bins are changed, removed, etc), so it makes more sense
    // to just compute the cdf values here.
    if (CDF) {

        // using new vector bins. store cdf in separate field.
        for (int i = 0; i < vbins.size(); i++) {
            if (!i) vbins[i].cdf = vbins[i].count;
            else vbins[i].cdf = vbins[i].count + vbins[i-1].cdf;
        }
        for (int i = 0; i < vbins.size(); i++) vbins[i].cdf /= (double)num_lines;

        // now, convert
    }

    // lets normalize the histogram if specified
    if (NORM) {
        double area=0.0;
        std::cout << "normalizing" << std::endl;
        // first lets find the area of the histogram
        for (int i = 0; i < vbins.size(); i++) {
            //area += bins[i][1];
            area += vbins[i].pdf * 2 * vbins[i].width;
        }
        //area *= bin_width;
        std::cout << area << std::endl;
        // now normalize

        double volume = pow((250./.7),3);
        double veffective = 0;
        double area2=0.0;

        if (MASSFUNC) {
            std::cout << "Normalizing to SMF Limits" << std::endl;

            // compute effective volume
            if (compute_veff (veffective)) return -1;
        }

        if (area > 0.0) {
            for (int i = 0; i < vbins.size(); i++) {
                //bins[i][1] /= area;
                vbins[i].pdf /= area;
                if (MASSFUNC) {
                    vbins[i].pdf *= (num_lines/veffective);
                    //bins[i][1] *= (num_lines/volume);
                }
                //area2 += bins[i][1];
                area2 += vbins[i].pdf * 2 * vbins[i].width;
            }
            //area2 *= bin_width;
            std::cout << area2 << std::endl;
        }
    }

	// now write 1d histogram to output file
	std::ofstream of_hist;

	setOutputFileName(x);

	of_hist.open(output.c_str(),std::ofstream::out);

    // total/cumm area used for computing density cdf - whole box cdf comparison
    //double total_area=1;
    //double total_area = pdf_integrate (-std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),x);

	for (int i = 0; i < vbins.size(); i++) {
        if (doingBacksplashFraction) {
            of_hist << bins_bs[i][0] << " " << bins_bs[i][1]/(bins_bs[i][1]+bins_bs[i][2]) << " " << bins_bs[i][2]/(bins_bs[i][1]+bins_bs[i][2]) << " " << bins_bs[i][3]/(bins_bs[i][1]+bins_bs[i][2]) << " " << bins_bs[i][4]/(bins_bs[i][1]+bins_bs[i][2]) << std::endl;
        }
        else {
            if (LOGZ) {
		        of_hist << bins[i][0] << " " << log10(bins[i][1]) << std::endl;
            }
            else {
		        //of_hist << bins[i][0] << " " << bins[i][1] << std::endl;
                //double cumm_area=1;
                //double cumm_area = pdf_integrate (-std::numeric_limits<double>::infinity(),vbins[i].center,x);

		        of_hist << vbins[i].center << " " << vbins[i].pdf << " " << 2 * vbins[i].width << " " << vbins[i].cdf << " " << vbins[i].count << std::endl;
		        //of_hist << vbins[i].center << " " << vbins[i].pdf << " " << 2 * vbins[i].width << " " << vbins[i].cdf << " " << cumm_area/total_area << std::endl;
            }
        }
	}
	
	of_hist.close();

	return 0;
}

int doHist2D(double ** &data, int num_lines, int x, int y, int num_bins_x, int num_bins_y, double bin_width_x, double bin_width_y, bool LOGZ, bool NORM, bool SMOOTH, bool RANK, int POINTS_PER_CELL, bool SCATTER) {

	// now lets begin with the actual histogram
	double minX=-1, minY=-1, maxX=-1, maxY=-1;
	int indexX, indexY;

    typedef struct {
        double x;
        double y;
        double prank_x = 0;
        double prank_y = 0;
        double z = 0;
    } element;

    std::vector<element> vdata;

    // choose any selection criteria here
    for (int i = 0; i < num_lines; i++) {
       element e;
       e.x = data[x][i];
       e.y = data[y][i];
       vdata.push_back(e);
    }

    std::vector<element *> scatter;

    if (RANK) {

        //std::cout << "vector built, sorting x..." << std::endl;

        struct sort_comp_x {
            bool operator()(const element &left, const element &right) {
                return left.x < right.x;
            }
        };

        struct sort_comp_y {
            bool operator()(const element &left, const element &right) {
                return left.y < right.y;
            }
        };

        std::sort(vdata.begin(),vdata.end(),sort_comp_x());

        int i = 0;

        for (std::vector<element>::iterator it = vdata.begin(); it != vdata.end(); it++) {

            it->prank_x = (i++)/(double)num_lines*100.;
        }

        //std::cout << "sorting y..." << std::endl;

        std::sort(vdata.begin(),vdata.end(),sort_comp_y());

        i = 0;

        for (std::vector<element>::iterator it = vdata.begin(); it != vdata.end(); it++) {

            it->prank_y = (i++)/(double)num_lines*100.;
        }

	    minX = 0.;
	    maxX = 100.;
	    minY = 0.;
	    maxY = 100.;
    }

    else {

        bool notSet = true;

        for (int i = 0; i < num_lines; i++) {
            if (std::isfinite(vdata[i].x) && std::isfinite(vdata[i].y)) {
	    	    if (notSet) {
	    	    	minX = vdata[i].x;
	    	    	maxX = vdata[i].x;
	    	    	minY = vdata[i].y;
	    	    	maxY = vdata[i].y;
                    notSet = false;
	    	    }

	    	    if (vdata[i].x < minX) minX = vdata[i].x;
	    	    if (vdata[i].x > maxX) maxX = vdata[i].x;
	    	    if (vdata[i].y < minY) minY = vdata[i].y;
	    	    if (vdata[i].y > maxY) maxY = vdata[i].y;
            }
	    }
    }

	std::cout << "minX: " << minX << std::endl << "maxX: " << maxX << std::endl;
	std::cout << "minY: " << minY << std::endl << "maxY: " << maxY << std::endl;

    // I'm adding in the option of specifying the numerical bin widths.
    // If these are specified on the command line, they will take precedent to
    // any specification number of bins.  If only x_bin_width is specified at
    // the command line, num_bins_x will be determined here, bin_width_y will be fixed to
    // bin_width_x, and num_bins_y also determined here.
	if (bin_width_x == 0.0) bin_width_x = (maxX-minX)/(double)num_bins_x;
    else num_bins_x = ceil((maxX-minX)/bin_width_x);

    if (bin_width_y == 0.0) {
        if (!num_bins_y) {
            bin_width_y = bin_width_x;
            num_bins_y = ceil((maxY-minY)/bin_width_y);
        }

        else bin_width_y = (maxY-minY)/(double)num_bins_y;
    }

    else num_bins_y = ceil((maxY-minY)/bin_width_y);

    typedef struct {
        double x;
        double y;
        double z;
        double z_smooth;
        std::vector<element *> points;
    } cell;

    std::vector<cell> bins;

	for (int i = 0; i < num_bins_x; i++) {
		for (int j = 0; j < num_bins_y; j++) {

            cell c;
            c.x = minX + i * bin_width_x;
            c.y = minY + j * bin_width_y;
            c.z = 0.;
            bins.push_back(c);
		}
	}

	for (int i = 0; i < num_lines; i++) {
		indexX = (int) floor (( (RANK?vdata[i].prank_x:vdata[i].x) - minX) / bin_width_x);
		if (indexX < 0) indexX = 0;
		else if (indexX >= num_bins_x) indexX = num_bins_x - 1;
		indexY = (int) floor (( (RANK?vdata[i].prank_y:vdata[i].y) - minY) / bin_width_y);
		if (indexY < 0) indexY = 0;
		else if (indexY >= num_bins_y) indexY = num_bins_y - 1;

        bins[indexX * num_bins_y + indexY].z++;
        bins[indexX * num_bins_y + indexY].points.push_back(&vdata[i]);
	}

    double single_hit_value = 1.0;

    // lets normalize the histogram if specified
    if (NORM) {
        double area=0.0;
        std::cout << "normalizing" << std::endl;
        
        // first lets find the area of the histogram
        for (int i = 0; i < num_bins_x; i++) {
            for (int j = 0; j <num_bins_y; j++) {
                area += bins[i*num_bins_y + j].z;
            }
        }
        area *= bin_width_x * bin_width_y;
        std::cout << area << std::endl;

        single_hit_value = 1.0 / area;
        
        if (LOGZ) {
            std::cout << "cb range min: " << log10(single_hit_value) << std::endl;
        }

        // now normalize
        double area2=0.0;
        double max=0.0;
        if (area > 0.0) {
            for (int i = 0; i < num_bins_x; i++) {
                for (int j = 0; j < num_bins_y; j++) {
                    bins[i*num_bins_y + j].z /= area;
                    area2 += bins[i*num_bins_y + j].z;
                }
            }
            area2 *= bin_width_x * bin_width_y;
            std::cout << area2 << std::endl;
        }
    }

	// now write 2d histogram to output file
	std::ofstream of_hist;

	setOutputFileName(x,y);

	of_hist.open(output.c_str(),std::ofstream::out);

	for (int i = 0; i < num_bins_x * num_bins_y; i++) {
		if (LOGZ) {
			if (bins[i].z == 0) bins[i].z = -10.0; //-1.0;
			else bins[i].z = log10(bins[i].z);
		}
    }

    

    // do smoothing if requested
    if (SMOOTH) {
        double kernel[3][3] = {{1,2,1},{2,4,2},{1,2,1}};
        for (int i = 1; i < num_bins_x - 1; i++) {
            for (int j = 1; j < num_bins_y - 1; j++) {
                double newVal = 0.0;
                for (int k = -1; k < 2; k++) {
                    for (int l = -1; l < 2; l++) {
                        newVal += bins[(i+k)*num_bins_y+(j+l)].z*kernel[k+1][l+1];
                    }
                }
                newVal /= 16.0;
                int index = i * num_bins_y + j;
                bins[index].z_smooth = newVal;
            }
        }
        // update z to be consistent with other flags/options
        for (int i = 0; i < num_bins_x * num_bins_y; i++) bins[i].z = bins[i].z_smooth;
    }

    // do scatter plot halo selection if chosen
    if (SCATTER) {

        for (int i = 0; i < num_bins_x * num_bins_y; i++) {
            
            int count = 0;

            while ((count < bins[i].points.size()) && (count < POINTS_PER_CELL)) {

                scatter.push_back(bins[i].points[count++]);
                scatter.back()->z = bins[i].z;
            }
        }

        for (int i = 0; i < scatter.size(); i++) {

            of_hist << scatter[i]->x << " " << scatter[i]->y << " " << scatter[i]->z << std::endl;
        }
    }

    // if not a scatter plot then do hist2D as usual
    else {
        for (int i = 0; i < num_bins_x * num_bins_y; i++) {
            of_hist << bins[i].x << " " << bins[i].y << " " << bins[i].z << std::endl;
	    }
    }

	of_hist.close();


	return 0;
}

int doRankCorrelation (double ** &data, int num_lines, int num_fields, int x, int y) {

    typedef struct {
        double x;
        double y;
        int rank_x = 0;
        int rank_y = 0;
    } element;

    std::vector<element> vdata;

    // choose any selection criteria here
    for (int i = 0; i < num_lines; i++) {
       element e;
       e.x = data[x][i];
       e.y = data[y][i];
       vdata.push_back(e);
    }

    std::cout << "vector built, sorting x..." << std::endl;

    struct sort_comp_x {
        bool operator()(const element &left, const element &right) {
            return left.x < right.x;
        }
    };

    struct sort_comp_y {
        bool operator()(const element &left, const element &right) {
            return left.y < right.y;
        }
    };

    std::sort(vdata.begin(),vdata.end(),sort_comp_x());

    int i = 0;

    for (std::vector<element>::iterator it = vdata.begin(); it != vdata.end(); it++) {

        it->rank_x = i++;
    }

    std::cout << "sorting y..." << std::endl;

    std::sort(vdata.begin(),vdata.end(),sort_comp_y());

    i = 0;

    for (std::vector<element>::iterator it = vdata.begin(); it != vdata.end(); it++) {

        it->rank_y = i++;
    }

    std::cout << "computing r..." << std::endl;

    // now compute spearman rank correlation coefficient
    double r = 0;
    for (std::vector<element>::iterator it = vdata.begin(); it != vdata.end(); it++) {
        r += pow((it->rank_x - it->rank_y),2.);
    }

    r = 1 - ( 6*r / ( vdata.size() * (pow(vdata.size(),2.) - 1.) ) );

    std::cout << "Spearman rank correlation coefficient is: " << r << std::endl;

    return 0;
}

//====================================================
// doGenPlotStats:
// Computes several statistics to be used in plotting data involving field x and y,
// such as median x, median y, and best fit lines.
//====================================================
int doGenPlotStats (double ** &data, int num_lines, int num_fields, int x, int y) {

    std::vector<double> vx;
    std::vector<double> vy;

    for (int i = 0; i < num_lines; i++) {
        vx.push_back(data[x][i]);
        vy.push_back(data[y][i]);
    }

    std::sort(vx.begin(),vx.end());
    std::sort(vy.begin(),vy.end());

    double x_median, y_median;

    // if odd in length, take middle element
    if (vx.size() % 2) {
        x_median = vx[(vx.size()-1)/2];
    }
    // else average two middle elements
    else {
        x_median = vx[vx.size()/2]/2. + vx[(vx.size()/2)-1]/2.;
    }

    // if odd in length, take middle element
    if (vy.size() % 2) {
        y_median = vy[(vy.size()-1)/2];
    }
    // else average two middle elements
    else {
        y_median = vy[vy.size()/2]/2. + vy[(vy.size()/2)-1]/2.;
    }

    std::cout << "x_median: " << x_median << std::endl;
    std::cout << "y_median: " << y_median << std::endl;

    // tablulate median x and y lines
    /*int num_points = 100;
    double x_median_line[num_points], y_median_line[num_points];

    double interval = (vy[vy.size()-1] - vy[0])/(double)num_points;
    double value = vy[0];

    for (int i = 0; i < num_points; i++) {
    */

    return 0;
}

std::vector<spline> cubicSpline (std::vector<double> &x, std::vector<double> &y) {

    // Natural cubic spline has form
    // S_i(x) = a_i + b_i * (x - x_i) + c_i * (x - x_i)^2 + d_i * (x - x_i)^3
    // where we can replace (x - x_i) with another parameter h_i
    // S_i(x) = a_i + b_i * h_i + c_i * h_i^2 + d_i * h_i^3

    // Boundary conditions on interior nodes are (from n = 1 to n):
    // S_i(0) = S_i-1(h_i-1)            (continuous 0th order)
    // S_i'(0) = S_i-1'(h_i-1)          (continuous 1st order)
    // S_i''(0) = S_i-1''(h_i-1)        (continuous 2nd order)

    // this provides 3*(n-1) equations

    // S_i(0) = y_i

    // provides another (n+1) equations

    // For end nodes, we also require
    // S_0''(0) = S_n-1''(h_n-1) = 0

    // altogether we have 3(n-1) + (n+1) + 2 = 4n equations and 4n unknowns

    // trivially, we have that:
    // S_i(0) = a_i = y_i

    //std::vector<double> x(5);
    //std::vector<double> y(5);

    ////choose some random coordinates
    //x[0] = 1.;
    //x[1] = 3.;
    //x[2] = 3.4;
    //x[3] = 4.4;
    //x[4] = 6.1;

    //y[0] = -2;
    //y[1] = 1;
    //y[2] = 1.1;
    //y[3] = 3.4;
    //y[4] = 2.8;

    // establish n, where number of input coords is n+1
    int n = x.size()-1;

    // initialize vectors for computing splines
    std::vector<double> a(n+1);
    std::vector<double> b(n);
    std::vector<double> d(n);
    std::vector<double> h(n);
    std::vector<double> alpha(n);
    std::vector<double> c(n+1);
    std::vector<double> l(n+1);
    std::vector<double> mu(n+1);
    std::vector<double> z(n+1);

    // solve trivially for a 
    for (int i = 0; i < n+1; i++) {
        a[i] = y[i];
    }

    // determine h
    for (int i = 0; i < n; i++) {
        h[i] = x[i+1] - x[i];
    }

    for (int i = 1; i < n; i++) {
        alpha[i] = 3*(a[i+1]-a[i])/h[i] - 3*(a[i]-a[i-1])/h[i-1];
    }

    l[0] = 1;
    mu[0] = 0;
    z[0] = 0;

    for (int i = 1; i < n; i++) {
        l[i] = 2*(x[i+1]-x[i-1])-h[i-1]*mu[i-1];
        mu[i] = h[i]/l[i];
        z[i] = (alpha[i]-h[i-1]*z[i-1])/l[i];
    }

    l[n] = 1;
    z[n] = 0;
    c[n] = 0;

    for (int j = n-1; j >= 0; j--) {
        c[j] = z[j] - mu[j]*c[j+1];
        b[j] = (a[j+1]-a[j])/h[j] - h[j]*(c[j+1]+2*c[j])/3.;
        d[j] = (c[j+1]-c[j])/(3.*h[j]);
    }

    std::vector<spline> result(n);

    for (int i = 0; i < n; i++) {
        spline s;
        s.a = a[i];
        s.b = b[i];
        s.c = c[i];
        s.d = d[i];
        s.x = x[i];
        result[i] = s;
    }

    // now output a trajectory
    //std::ofstream of;
    //setOutputFileName();
    //of.open(output.c_str(),std::ofstream::out);

    //int line = 0;
    //int steps = 100;
    //double stepsize = 0;
    //double p = 0;
    //double q = 0;

    //for (int i = 0; i < x.size()-1; i++) {
    //    for (int j = 0; j < steps; j++) {
    //        double stepsize = (x[i+1]-x[i])/steps;
    //        p = stepsize*j;
    //        q = a[i]+b[i]*p+c[i]*p*p+d[i]*p*p*p;
    //        of << x[i]+p << " " << q;
    //        if (line < x.size()) {
    //            of << " " << x[line] << " " << y[line];
    //            line++;
    //        }
    //        of << std::endl;
    //    }
    //}

    //// make sure to connect to last node
    //of << x.back() << " " << y.back();

    //of.close();

    return result;
}
