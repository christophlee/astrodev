#include "include.h"
#include "helper.h"
#include "io.h"
#include "stats.h"
#include "density.h"
#include "tree.h"
#include "fileop.h"
#include "profiles.h"
#include "web.h"
#include "viz.h"
#include "sam_am.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <string.h>
#include <cmath>
#include <algorithm>
#include <vector>
#include <limits>
#include <iterator>

bool hist1D = false;
bool hist2D = false;
bool avg = false;
bool median = false;
bool binnedHist2D = false;
bool binnedMedian = false;
bool binnedHist1D = false;
bool bin = false;
bool toBin = false;
bool toAscii = false;
bool treeWalk = false;
bool catalogMerge = false;
bool other = false;
bool LOG[3] = {false,false,false};
std::string fileName;
std::string output;

int id = 0;

BP_CAT_DEF BP_CAT;

BP_TREE_DEF BP_TREE;

int doOther (double ** &data, int num_lines, int num_fields, int num_bins, int x, int y, bool NORM, bool binary_write, std::string otherparam, bool mmp);

int setType(std::string arg) {

	if (arg == "hist1D") {
		hist1D = true;
	}
	else if (arg == "hist2D") {
		hist2D = true;
	}
	else if (arg == "avg") {
		avg = true;
	}
	else if (arg == "median") {
		median = true;
	}
    else if (arg == "bin") {
        bin = true;
    }
    else if (arg == "binnedHist2D") {
        binnedHist2D = true;
    }
    else if (arg == "binnedHist1D") {
        binnedHist1D = true;
    }
    else if (arg == "binnedMedian") {
        binnedMedian = true;
    }
    else if (arg == "toBin") {
        toBin = true;
    }
    else if (arg == "toAscii") {
        toAscii = true;
    }
    else if (arg == "treeWalk") {
        treeWalk = true;
    }
    else if (arg == "merge") {
        catalogMerge = true;
    }
    else if (arg == "other") {
        other = true;
    }
	else {
		std::cerr << "Error: type \"" << arg << "\" not recognized." << std::endl;
		return -1;
	}

	return 0;
}

bool isTypeSet() {
    return (hist1D || hist2D || avg || median || bin || binnedHist2D || binnedHist1D || binnedMedian || toBin || toAscii || treeWalk || catalogMerge || other);
}

int main (int argc, char ** argv) {

	int x=0, y=1, num_bins=1, num_bins_x=0, num_bins_y=0, binx=0, auto_tnum=-1, use_rw=-1, POINTS_PER_CELL = 3;
    bool NORM = false, SMOOTH = false, logbinx = false, binary_read = false, binary_write = false;
    bool binary2_read = false, autobin = true, use_percentiles_x = false, use_percentiles_y = false;
    bool RANK = false, CDF= false, mmp=false, MASSFUNC = false, PROFILES = false, SCATTER = false;
    double bin_width_x = 0.0, bin_width_y = 0.0, mstar = -1.;
    //double * thresholds = NULL;
    threshold_object thresholds;
    int threshold_type = 0;
    std::string otherparam = "", mergeTarget = "";

	// lets parse our command line arguments
	for (int i = 1; i < argc; i++) {
		
		// all arguments should be composed like:
		// -file=filename -x=2 -y=3 etc
		if (argv[i][0] != '-') {
			std::cerr << "Error: command line arguments malformed, quitting." << std::endl;
			return -1;
		}
		else {
			std::string arg = std::string(argv[i]);
			std::string flag = arg.substr(1,arg.find('=')-1);
			arg = arg.substr(arg.find('=')+1);
			
			if (flag == "file") {
                // set input file name
				fileName = arg;
			}
            else if (flag == "br") {
                // indicate we are reading a binary file
                binary_read = true;
            }
            else if (flag == "bw") {
                // indicate we are writing any binned catalogs in binary format
                binary_write = true;
            }
			else if (flag == "type") {
                // set type of analysis
				if (setType(arg)) return -1;
			}
			else if (flag == "x") {
                // set x column of data for analysis
				x = std::strtol(arg.c_str(),NULL,10);
			}
			else if (flag == "y") {
                // set y column of data for analysis
				y = std::strtol(arg.c_str(),NULL,10);
			}
            else if (flag == "binx") {
                // set x column of data for file binning
                binx = std::strtol(arg.c_str(),NULL,10);
            }
            else if (flag == "nbins") {
                // can specify number of bins or just the thresholds
                // to specify thresholds, use:
                // -nbins=[10.2,10.6,11.2]
                // etc, which will give 4 bins with said thresholds
                // OR can use
                // -nbins=%[10,25,75,90]
                // which would give 5 bins of said percentiles,
                // bin 0 would be 0-10% rank ordered halos in binning
                // variable, bin 1 would be 10-25%, etc
                if (arg[0] == '%' && (arg[1] == '[' || (arg[1] == 'r' && arg[2] == '['))) {
                    threshold_type = PERCENT_BINNING;
                    arg = arg.substr(1);
                }

                if (arg[0] == 'r' && arg[1] == '[') {
                    arg = arg.substr(1);
                    thresholds.range = true;
                }

                if (arg[0] == '[') {

                    // set the thresholds
                    int ind = 1, count = 1;

                    //std::cout << arg << std::endl;

                    // do this twice, on the first pass, determine how many bins there are
                    // and on the second pass, allocate space and store them
                    for (int i = 0; i < 2; i++) {
    
                        // handle second pass differently
                        if (i) {
                            num_bins = count;
                            //std::cout << "num bins: " << num_bins << std::endl;
                            ind = 1;
                            count = 1;

                            // range specified bins will have one less than standard, since
                            // the bin range is completely specified in each expression, rather
                            // than including the low and high extreme bins (extra count)
                            if (thresholds.range) num_bins--;

                            //std::cout << "num_bins: " << num_bins << std::endl;

                            thresholds.minValues = new double[num_bins];
                            thresholds.maxValues = new double[num_bins];
                        }

                        while (true) {
                            int nextComma = arg.find(',',ind);
                            if (nextComma == -1) nextComma = arg.length()-1;
                            //std::cout << "nextComma: " << nextComma << ": ";
                            if (i) {
                                if (thresholds.range) {
                                    std::string errmsg = "";
                                    // note that "count" skips the zeroth index for non range bins, so that it can be filled in later, so
                                    // here we need to reduce it because we want to specify the zeroth index here
                                    if (handleRangeExpression(arg.substr(ind,nextComma-ind),thresholds.minValues[count-1],thresholds.maxValues[count-1],errmsg) == -1) {
                                        std::cerr << errmsg << std::endl;
                                        return -1;
                                    }

                                    //std::cout << thresholds.minValues[count-1] << " " << thresholds.maxValues[count-1] << std::endl;
                                }

                                else {
                                    thresholds.minValues[count] = std::strtod((arg.substr(ind,nextComma-ind)).c_str(),NULL);
                                }
                            }
                            //std::cout << std::strtod((arg.substr(ind,nextComma-ind)).c_str(),NULL) << " ";
                            ind = nextComma+1;
                            count++;
                            if (ind == arg.length()) break;
                        }
                    }
                }
                else if (arg == "3Q") {
                    thresholds.minValues = new double[3];
                    thresholds.maxValues = new double[3];
                    num_bins = 3;
                    // just indicate for later (by assigning -1 to each threshold)
                    // that we want the binning function to use the 3 quartile
                    // method
                    //for (int i = 0; i < 3; i++) thresholds[i] = -1;
                    threshold_type = THREEQ_BINNING;
                }
                else {
                    // set num bins for file binning
                    num_bins = std::strtol(arg.c_str(),NULL,10);
                }
            }
			else if (flag == "nbinsx") {
                // set num bins on x axis for analysis
				num_bins_x = std::strtol(arg.c_str(),NULL,10);
			}
			else if (flag == "nbinsy") {
                // set num bins on y axis for analysis
				num_bins_y = std::strtol(arg.c_str(),NULL,10);
			}
            else if (flag == "binwidthx") {
                // set numerical width of x bins
                bin_width_x = std::strtod(arg.c_str(),NULL);
            }
            else if (flag == "binwidthy") {
                // set numerical width of y bins
                bin_width_y = std::strtod(arg.c_str(),NULL);
            }
            else if (flag == "autobin") {
                // turn on auto-rebinning
                if (arg == "false") autobin = false;
            }
            else if (flag == "autotnum") {
                // set threshold for minimum bin count
                auto_tnum = std::strtol(arg.c_str(),NULL,10);
            }
            else if (flag == "logx") {
                // log x data before using for analysis
                LOG[0] = true;
            }
            else if (flag == "logy") {
                // log y data before using for analysis
                LOG[1] = true;
            }
			else if (flag == "logz") {
                // log z data (contour level / colorbar) before outputting
				LOG[2] = true;
			}
            else if (flag == "logbinx") {
                logbinx = true;
            }
            else if (flag == "norm") {
                // normalize result (hist1D or hist2D)
                NORM = true;
            }
            else if (flag == "smooth") {
                // smooth result with gaussian kernel
                SMOOTH = true;
            }
			else if (flag == "output") {
                // set output file name
				output = arg;
			}
            else if (flag == "exrad") {
                // exclusion radius for local density calculation
                exrad = std::strtod(arg.c_str(),NULL);
            }
            else if (flag == "smoothing_length" || flag == "slen") {
                // smoothing radius for local density calculation
                slen = std::strtod(arg.c_str(),NULL);
            }
            else if (flag == "id") {
                id = std::strtol(arg.c_str(),NULL,10);
            }
            else if (flag == "otherparam") {
                otherparam = arg;
            }
            else if (flag == "mstar") {
                mstar = std::strtod(arg.c_str(),NULL);
            }
            else if (flag == "use_percentiles_x") {
                use_percentiles_x = true;
            }
            else if (flag == "use_percentiles_y") {
                use_percentiles_y = true;
            }
            else if (flag == "rank") {
                // use percentile rank for x-y hist2d plot
                RANK = true;
            }
            else if (flag == "cdf") {
                // produce a cdf instead of a pdf for 1d distributions
                CDF = true;
            }
            else if (flag == "mmp") {
                mmp = true;
            }
            else if (flag == "use_rw") {
                // run analysis over reweighted mass bins (must be precomputed)
                use_rw = std::strtol(arg.c_str(),NULL,10);
            }
            else if (flag == "massfunc") {
                MASSFUNC = true;
            }
            else if (flag == "mergeTo") {
                mergeTarget = arg;
            }
            else if (flag == "profiles") {
                PROFILES = true;
            }
            else if (flag == "scatter") {
                SCATTER = true;
            }
            else if (flag == "points_per_cell") {
                POINTS_PER_CELL = std::strtol(arg.c_str(),NULL,10);
            }
			else {
				std::cerr << "Error: flag not recognized: " << flag << std::endl;
			}
		}
	}

    if (!isTypeSet()) {
        std::cerr << "Error: type of analysis must be specified." << std::endl;
    }

	if (!fileName.length()) {
		std::cerr << "Error: no file name specified." << std::endl;
		return -1;
	}

    if (output.length()) parseOutputFileName(fileName);

    if (!binary_read) binary_read = isBinaryExt (fileName);
    if (!binary_write) binary_write = isBinaryExt (output);

    binary2_read = isBinary2Ext (fileName);

    //std::cout << "output file name: " << output << std::endl;

	//std::cout << x << " " << y << std::endl;

    double ** data = NULL;
    int num_fields=0, num_lines=0;

    if (PROFILES) {

        doProfileAnalysis(fileName, binary_read);
        //if (binary_read) {
        //    if (readProfileBinaryFile(fileName)) {
        //        return -1;
        //    }
        //}
        //else {
        //    if (readProfileFile(fileName)) {
        //        return -1;
        //    }
        //}

        //doFitProfiles(halos,"nfw");
    }
    
    else {
        if (binary_read) {
            if (readBinaryFile(fileName,num_lines,num_fields,data)) {
                return -1;
            }
        }
        else if (binary2_read) {
            if (readBinary2File(fileName,num_lines,num_fields,data)) {
                return -1;
            }
        }
        else {
            if (readFile(fileName,num_lines,num_fields,data)) {
                return -1;
            }
        }
    }

    // log data for specified columns
    if (LOG[0]) for (int i = 0; i < num_lines; i++) data[x][i] = log10(data[x][i]);
    if (LOG[1]) for (int i = 0; i < num_lines; i++) data[y][i] = log10(data[y][i]);
    if (logbinx) for (int i = 0; i < num_lines; i++) data[binx][i] = log10(data[binx][i]);

	// now we choose which type of analysis to do
	if 	    (hist1D) 	    doHist1D(data,num_lines,x,num_bins_x,NORM,LOG[2],CDF,autobin,auto_tnum,MASSFUNC,bin_width_x);
	else if	(hist2D)	    doHist2D(data,num_lines,x,y,num_bins_x,num_bins_y,bin_width_x,bin_width_y,LOG[2],NORM,SMOOTH,RANK,POINTS_PER_CELL,SCATTER);
	else if	(avg)		    doAvg(data,num_lines,x,y,num_bins_x);
	else if	(median)	    doMedian(data,num_lines,x,y,num_bins_x,autobin,auto_tnum,use_percentiles_x,use_percentiles_y);
    else if (bin)           doBin(data,num_lines,num_fields,binx,num_bins,thresholds,threshold_type,binary_write);
    else if (binnedHist2D || binnedMedian || binnedHist1D)
                            doBinnedAnalysis(data,num_lines,num_fields,binx,num_bins,thresholds,threshold_type,x,y,num_bins_x,num_bins_y,bin_width_x,bin_width_y,LOG,NORM,SMOOTH,binary_write,autobin,auto_tnum,use_percentiles_x,use_percentiles_y,RANK,CDF,use_rw,MASSFUNC,POINTS_PER_CELL,SCATTER);
    else if (toBin) {        if (binary2_read) writeOutputFileBinary ("", data, num_lines, num_fields, 1);
                            else              writeOutputFileBinary ("", data, num_lines, num_fields); }
    else if (toAscii) {      if (binary2_read) writeOutputFile (data, num_lines, num_fields, 1);
                            else              writeOutputFile (data, num_lines, num_fields); }
    else if (treeWalk)      doTreeWalk (data, num_lines, num_fields, mmp);
    else if (catalogMerge)  doCatalogMerge (data, num_lines, num_fields, mergeTarget);
    else if (other)         doOther(data,num_lines,num_fields,num_bins_x,x,y,NORM,binary_write,otherparam,mmp);

    if (data) {
        if (binary2_read) free (data, num_lines);
        else free (data, num_fields);
    }

    return 0;
}

int doOther (double ** &data, int num_lines, int num_fields, int num_bins, int x, int y, bool NORM, bool binary_write, std::string otherparam, bool mmp) {

    std::string desc = std::string("sam centrals");
    std::string nullstr = std::string("");

    if (otherparam == "") {

        //doUnmatchedHalos (data,num_lines);

        //doCreateLimitedCatalog (data, num_lines, num_fields, x, 12.1, -1, nullstr);

        //doCalcQuiescentFractions (data, num_lines, num_fields, x, y);

        //doCalcStellarMassFractions (data, num_lines, num_fields, x, y);

        //doGaussianSmoothDensity (data, num_lines, num_fields, x, y);
        //doGaussianSmoothDensity3D (data, num_lines, num_fields, x, y);

        //doBuildDensityCatalog (data, num_lines, num_fields, x, y);

        //doCalcVoxelDensityDistributions (x, y, num_bins, NORM);

        //doCombineHlistFiles (data, num_lines, num_fields);

        // doAddFieldToCatalog (data, num_lines, num_fields, binary_write);

        //doCombineSpineCatalog (data, num_lines, num_fields, binary_write);

        //doDensityNormalization (data, num_lines, num_fields);

        // do nothing
        std::cout << "Doing nothing... " << std::endl;
    }

    else if (otherparam == "gaussianSmoothDensity") {
        doGaussianSmoothDensity (data, num_lines, num_fields, x, y);
    }

    else if (otherparam == "doCalcVoxelDensityDistributions") {
        doCalcVoxelDensityDistributions (x, y, num_bins, NORM);
    }

    else if (otherparam == "combineHlistFiles") {
        doCombineHlistFiles (data, num_lines, num_fields);
    }
    else if (otherparam == "buildDensityCatalog") {
        doBuildDensityCatalog (data, num_lines, num_fields, x, y);
    }
    else if (otherparam == "addFieldToCatalog") {
        doAddFieldToCatalog (data, num_lines, num_fields, binary_write);
    }
    else if (otherparam == "stackCatalogs") {
        doStackCatalogs();
    }
    else if (otherparam == "rankCorrelation") {
        doRankCorrelation (data, num_lines, num_fields, x, y);
    }
    else if (otherparam == "genPlotStats") {
        doGenPlotStats (data, num_lines, num_fields, x, y);
    }
    else if (otherparam == "doProgenitorHistory") {
        doProgenitorHistory (data, num_lines, num_fields);
    }
    else if (otherparam == "doDensityPercentileEvolution") {
        doDensityPercentileEvolution (data, num_lines, num_fields, x);
    }
    else if (otherparam == "doReweightMassDistribution") {
        doReweightMassDistribution (data, num_lines, num_fields, x);
    }
    else if (otherparam == "doStackProfiles") {
        doStackProfiles (data, num_lines, num_fields);
    }
    else if (otherparam == "doSelectZ") {
        doSelectZ (data, num_lines, num_fields, x);
    }
    else if (otherparam == "doConvertMergerTracks") {
        doConvertMergerTracks (data, num_lines, num_fields);
    }
    else if (otherparam == "doTreeAnalysis") {
        doTreeAnalysis (data, num_lines, num_fields);
    }
    else if (otherparam == "doReorderFields") {
        doReorderFields (data, num_lines, num_fields);
    }
    else if (otherparam == "doReorderFields2") {
        doReorderFields2 (data, num_lines, num_fields);
    }
    else if (otherparam == "doSpineTest") {
        doSpineTest (data, num_lines, num_fields);
    }
    else if (otherparam == "doSpineCatalogAnalysis") {
        doSpineCatalogAnalysis (data, num_lines, num_fields, x);
    }
    else if (otherparam == "doConvertToVizFormat") {
        doConvertToVizFormat (data, num_lines, num_fields);
    }
    else if (otherparam == "doIndvHaloTracks") {
        doIndvHaloTracks (data, num_lines, num_fields, x);
    }
    else if (otherparam == "doFullTreeAnalysis") {
        doFullTreeAnalysis (data, num_lines, num_fields, mmp);
    }
    //else if (otherparam == "cubicSpline") {
    //    cubicSpline ();
    //}
    else {
        std::cerr << "other param: \"" << otherparam << "\" not recognized, exiting..." << std::endl;
        return -1;
    }

    return 0;
}
