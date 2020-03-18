#include "density.h"

float exrad = 0, slen = 0;


int doDensityNormalization (double ** &data, int num_lines, int num_fields) {

    for (int i = 0; i < num_lines; i++) {
        for (int j = 72; j < 78; j++) {
            data[j][i] /= 8.0;
        }
    }

    writeOutputFileBinary ("", data, num_lines, num_fields);

    return 0;
}

int doHist1DVoxels(double *** &data, int dim, int num_bins, bool NORM) {

	// now lets begin with the actual histogram
	double minX, maxX;
	double bin_size;
	int index;

    bool use_log_log = true;

    //for (int i = 0; i < dim; i++) {
    //    for (int j = 0; j < dim; j++) {
    //        for (int k = 0; k < dim; k++) {
    //            data[i][j][k] /= 8.0;
    //        }
    //    }
    //}

    if (use_log_log) {
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                for (int k = 0; k < dim; k++) {
                    data[i][j][k] = log10(data[i][j][k]);
                }
            }
        }
    }

	for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            for (int k = 0; k < dim; k++) {
		        if ((i+j+k) == 0) {
			        minX = data[i][j][k];
			        maxX = data[i][j][k];
		        }
		        else {
                    if (data[i][j][k] < minX && std::isfinite(data[i][j][k])) minX = data[i][j][k];
		            if (data[i][j][k] > maxX && std::isfinite(data[i][j][k])) maxX = data[i][j][k];
                }
            }
        }
	}

	std::cout << "minX: " << minX << std::endl << "maxX: " << maxX << std::endl;

	bin_size = (maxX-minX)/(double)num_bins;

    std::cout << "bin_size: " << bin_size << std::endl;

	double bins[num_bins][2];

	for (int i = 0; i < num_bins; i++) {
		bins[i][0] = minX + i * bin_size;
		bins[i][1] = 0.0;
	}

	for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            for (int k = 0; k < dim; k++) {
		        index = (int) floor ((data[i][j][k] - minX) / bin_size);
		        if (index < 0) index = 0;
		        else if (index >= num_bins) index = num_bins - 1;
		        bins[index][1] += 1.0;
            }
        }
	}
    

    // lets normalize the histogram if specified
    if (NORM) {
        double area=0.0;
        std::cout << "normalizing" << std::endl;
        // first lets find the area of the histogram
        for (int i = 0; i < num_bins; i++) {
            area += bins[i][1];
        }
        area *= bin_size;
        std::cout << area << std::endl;
        // now normalize

        double area2=0.0;
        if (area > 0.0) {
            for (int i = 0; i < num_bins; i++) {
                bins[i][1] /= area;
                area2 += bins[i][1];
            }
            area2 *= bin_size;
            std::cout << area2 << std::endl;
        }
    }

	// now write 1d histogram to output file
	std::ofstream of_hist;

	setOutputFileName();

	of_hist.open(output.c_str(),std::ofstream::out);

    // adjust bin positions to reflect centers of bins, also output bin width (for pdf_integrate).
	for (int i = 0; i < num_bins; i++) {
		of_hist << bins[i][0]+bin_size/2. << " " << bins[i][1] << " " << bin_size << std::endl;
	}
	
	of_hist.close();

	return 0;
}

int genGaussianKernel (double * & kernel, int dim, double FWHM, double threshold = -1.0) {

    int length = ceil((float)dim/2.0);

    double sigma = FWHM/(2.0*sqrt(2*log(2.0)));

    kernel = new double [length];

    for (int i = 0; i < length; i++) kernel[i] = 0.0;

    for (int i = 0; i < length; i++) {
        double value = exp(-i*i/(2.0*sigma*sigma));
        value /= sqrt(2*M_PI*sigma*sigma);
        if (threshold < 0) threshold = value/1.e5;
        if (value >= threshold) kernel[i] = value;
        else {
            kernel[i] = 0.0;
            break;
        }

        std::cout << kernel[i] << " ";
    }

    std::cout << std::endl;
    
    double sum = 0.0;

    // normalize
    for (int i = (-length+1); i < length; i++) {
        if (i < 0) sum += kernel[-i];
        else sum += kernel [i];
    }

    for (int i = 0; i < length; i++) kernel[i] /= sum;

    std::cout << "Kernel normalization factor: " << sum << std::endl;

    //std::ofstream kernel_function ("kernel_function_gaussian.out",std::ios::out);

    //for (int i = -length+1; i < length; i++) {
    //    for (int j = -length+1; j < length; j++) {
    //        kernel_function << i << " " << j << " " << kernel[(int)fabs(i)]*kernel[(int)fabs(j)] << std::endl;
    //    }
    //}

    //kernel_function.close();

    return 0;
}

int genGaussianKernelExcludeCenter (double * & kernel, int dim, double FWHM, double threshold, double exclusion_radius) {

    int length = ceil((float)dim/2.0);

    double sigma = FWHM/(2.0*sqrt(2*log(2.0)));

    kernel = new double [length];

    for (int i = 0; i < length; i++) {
        if (i > exclusion_radius) {
            double value = exp(-i*i/(2.0*sigma*sigma));
            value /= sqrt(2*M_PI*sigma*sigma);
            if (value >= threshold) kernel[i] = value;
            else kernel[i] = 0.0;
        }
        else kernel[i] = 0.0;

        std::cout << kernel[i] << " ";
    }

    std::cout << std::endl;

    double sum = 0.0;

    // normalize
    for (int i = (-length+1); i < length; i++) {
        if (i < 0) sum += kernel[-i];
        else sum += kernel [i];
    }

    std::cout << "Kernel normalization factor: " << sum << std::endl;

    for (int i = 0; i < length; i++) {
        kernel[i] /= sum;
    }

    std::ofstream kernel_function ("kernel_function_ex.out",std::ios::out);

    for (int i = -length+1; i < length; i++) {
        for (int j = -length+1; j < length; j++) {
            kernel_function << i << " " << j << " " << kernel[(int)fabs(i)]*kernel[(int)fabs(j)] << std::endl;
        }
    }

    kernel_function.close();

    return 0;
}

int genGaussianKernelExcludeCenter3D (double *** & kernel, int dim, double FWHM, double threshold, double exclusion_radius) {

    int length = ceil((float)dim/2.0);
    double radius = 0.0, value = 0.0;

    double sigma = FWHM/(2.0*sqrt(2*log(2.0)));

    kernel = new double ** [length];
    for (int i = 0; i < length; i++) {
        kernel[i] = new double * [length];
        for (int j = 0; j < length; j++) {
            kernel[i][j] = new double [length];
            for (int k = 0; k < length; k++) {
                radius = sqrt(i*i+j*j+k*k);
                if (radius <= exclusion_radius) kernel[i][j][k] = 0.0;
                else {
                    value = exp(-radius*radius/(2*sigma*sigma));
                    value /= pow(2*M_PI,3./2.)*pow(sigma,3.0);
                    if (value < threshold) kernel[i][j][k] = 0.0;
                    else kernel[i][j][k] = value;
                }
            }
        }
    }

    double sum = 0.0;

    // normalize
    for (int i = -length+1; i < length; i++) {
        for (int j = -length+1; j < length; j++) {
            for (int k = -length+1; k < length; k++) {
                sum += kernel[(int)fabs(i)][(int)fabs(j)][(int)fabs(k)];
            }
        }
    }

    std::cout << "Kernel Normalization Factor: " << sum << std::endl;

    for (int i = 0; i < length; i++) {
        for (int j = 0; j < length; j++) {
            for (int k = 0; k < length; k++) {
                kernel[i][j][k] /= sum;
            }
        }
    }
    
    std::ofstream kernel_function ("kernel_function_ex3d.out",std::ios::out);

    for (int i = -length+1; i < length; i++) {
        for (int j = -length+1; j < length; j++) {
            kernel_function << i << " " << j << " " << kernel[(int)fabs(i)][(int)fabs(j)][0] << std::endl;
        }
    }

    kernel_function.close();

    return 0;
}

int doGaussianSmoothDensity3D (double ** & halo_data, int num_lines, int num_fields, int x, int y) {

    int rank=3, dim = 1024;

    double box_size = 250.0;    

    long double *** data;

    doBinaryRead3D ("bolshoi_plank/DougsCICs/bolshoi_plank_densities.dat",data,dim);

//    long double *** data2 = new long double ** [dim];
//    for (int i = 0; i < dim; i++) {
//        data2[i] = new long double * [dim];
//        for (int j = 0; j < dim; j++) {
//            data2[i][j] = new long double [dim];
//            for (int k = 0; k < dim; k++) data2[i][j][k] = 0.0;
//        }   
//    }

    double densities[num_lines];

    double *** kernel;

    // specify kernel properties: smoothing length, exclusion radius, cutoff
    genGaussianKernelExcludeCenter3D (kernel, 1024, slen, 1e-8, exrad);

    int kernel_length = 0;

    bool out_of_center = false;

    for (int i = 0; i < (int)ceil((float)dim/2.0); i++) {
        if (!kernel[i][0][0]) {
            if (out_of_center) {
                kernel_length = i;
                break;
            }
        }
        else out_of_center = true;
    }

    std::cout << "kernel_length: " << kernel_length << std::endl;

    int num_threads = 1;

#pragma omp parallel
{
    if (omp_get_thread_num() == 0) {
        std::cout << "Number of omp threads/procs: " << omp_get_num_threads() << "/" << omp_get_num_procs() << std::endl;
        num_threads = omp_get_num_threads();
    }
}

    long double sum = 0.0;

    int length_i = 0, length_j = 0, length_k = 0;
    int l2 = 0, m2 = 0, n2 = 0;


    int chunk_size = (int) ceil (num_lines / num_threads);

    std::cout << "Chunk size: " << chunk_size << std::endl;

    std::ofstream progress("progress",std::ios::out);

    // now compute the chunck allocation
    // lets just make 16 blocks, diving i, j into 4x4 squares, and for all of k

#pragma omp parallel for schedule(dynamic,chunk_size) private(sum,length_i,length_j,length_k,l2,m2,n2) shared(data)//,data2)

    for (int i = 0; i < num_lines; i++) {

        if (!omp_get_thread_num()) {
            if (i % chunk_size == 0) progress << "i range: " << i << " - " << i+chunk_size << std::endl;
            if (i % 1000 == 0) {
                progress << i/1000 << std::endl;
            }
            progress.flush();
        }

        // get halo position
        int position[3];
        //std::cout << i << ": " << "(" << halo_data[17][i] << " , " << halo_data[18][i] << " , " << halo_data[19][i] << ") : (";
        for (int j = 0; j < 3; j++) {
            position[j] = (int)floor(halo_data[17+j][i]/(box_size/(float)dim));
            if (position[j] >= dim || position[j] < 0) {
                std::cout << "Out of bounds("<<i<<"): " << position[j] << std::endl;
                if (position[j] >= dim) position[j] = dim-1;
                else if (position[j] < 0) position[j] = 0;
            }
            //std::cout << position[j] << " , ";
        }

        sum = 0.0;

        for (int l = position[0] - kernel_length; l < position[0] + kernel_length; l++) {
            for (int m = position[1] - kernel_length; m < position[1] + kernel_length; m++) {
                for (int n = position[2] - kernel_length; n < position[2] + kernel_length; n++) {
                    length_i = (int)fabs(l-position[0]);
                    length_j = (int)fabs(m-position[1]);
                    length_k = (int)fabs(n-position[2]);

                    l2 = (l < 0) ? (dim + l) : ( (l >= dim) ? (l - dim) : l);
                    m2 = (m < 0) ? (dim + m) : ( (m >= dim) ? (m - dim) : m);
                    n2 = (n < 0) ? (dim + n) : ( (n >= dim) ? (n - dim) : n);

                    sum += data[l2][m2][n2]*kernel[length_i][length_j][length_k];
                }
            }
        } 

        densities[i] = sum;

    }

    // update catalog (assuming input file is old catalog)
    //std::ofstream outfile("hlist_1.00230.list.testnew",std::ios::out);
    std::ofstream outfile("hlist_1.00230.list.testnew_"+std::to_string((long long)id),std::ios::out);

    outfile << "# " << id << " " << slen << " " << exrad << std::endl;

    bool pass = false;

    for (int i = 0; i < num_lines; i++) {

        //for (int j = 0; j < num_fields; j++) {
        //    if (!j) {
        //        outfile << halo_data[0][i];
        //    }
        //    else {
        //        if (j == 1)
        //            outfile << " " << (unsigned int)halo_data[j][i];
        //        //else if (j == 77 && !pass) {
        //        //    outfile << " " << densities[i];
        //        //    j--;
        //        //    pass = true;
        //        //}
        //        else
        //            outfile << " " << halo_data[j][i];
        //    }
        //}
        outfile << (unsigned int)halo_data[1][i];
        outfile << " " << densities[i];
        outfile << std::endl;
    }

    outfile.close();
    progress.close();

//    for (int b = 0; b < 16; b++) {
//
//        // determine which chunk you got
//        int i_start = b/4;
//        int j_start = b%4;
//        int length_i = 0, length_j = 0, length_k = 0, radius=0;
//        int l2 = 0, m2 = 0, n2 = 0;
//
//        if (!omp_get_thread_num()) std::cout << "thread " << omp_get_thread_num() << ": " << b << " (" << i_start << "," << j_start << ")" << std::endl;
//
//        for (int i = i_start*dim/4; i < (i_start+1)*dim/4; i++) {
//            for (int j = j_start*dim/4; j < (j_start+1)*dim/4; j++) {
//                for (int k = 0; k < dim; k++) {
//    
//                    if (!omp_get_thread_num() && !k) {
//                        progress << "0: " << i << " " << j << std::endl;
//                        progress.flush();
//                        std::cout << i << " " << j << std::endl;
//                    }

//                    sum = 0.0;

                    // now i, j, k is our origin point -> convolve with all points around it
//                    for (int l = 0; l < dim; l++) {
//                        for (int m = 0; m < dim; m++) {
//                            for (int n = 0; n < dim; n++) {
//                                length_i = (int)fabs(l-i);
//                                length_j = (int)fabs(m-j);
//                                length_k = (int)fabs(n-k);
//
//                                // handle periodic boundary
//                                if (length_i > dim/2) length_i = dim - length_i;
//                                if (length_j > dim/2) length_j = dim - length_j;
//                                if (length_k > dim/2) length_k = dim - length_k;
//
//                                radius = sqrt(length_i*length_i+length_j*length_j+length_k*length_k);
//
//                                if (radius <= kernel_length) sum += data[l][m][n]*kernel[length_i][length_j][length_k];
//                            }
//                        }
//                    }
//
//                    for (int l = i - kernel_length; l < i + kernel_length; l++) {
//                        for (int m = j - kernel_length; m < j + kernel_length; m++) {
//                            for (int n = k - kernel_length; n < k + kernel_length; n++) {
//                                length_i = (int)fabs(l-i);
//                                length_j = (int)fabs(m-j);
//                                length_k = (int)fabs(n-k);
//
//                                l2 = (l < 0) ? (dim + l) : ( (l >= dim) ? (l - dim) : l);
//                                m2 = (m < 0) ? (dim + m) : ( (m >= dim) ? (m - dim) : m);
//                                n2 = (n < 0) ? (dim + n) : ( (n >= dim) ? (n - dim) : n);
//
//                                sum += data[l2][m2][n2]*kernel[length_i][length_j][length_k];
//                            }
//                        }
//                    }           
//
//                    data2[i][j][k] = sum;
//                }
//            }
//        }
//        
//      //if (!omp_get_thread_num()) progress << std::endl;
//    }
//
//    progress << "writing output file" << std::endl;
//
//    std::ofstream outfile("BolshoiP_density_4_ex1_mpc2.dat",std::ios::out);
//    for (int i = 0; i < dim; i++) {
//        for (int j = 0; j < dim; j++) {
//            outfile.write((const char *)data2[i][j],dim*sizeof(long double));
//        }
//    }
//
//    outfile.close();
//
//    progress.close();

    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            delete [] data[i][j];
         //   delete [] data2[i][j];
        }
        delete [] data[i];
       // delete [] data2[i];
    }
    delete [] data;
   // delete [] data2;

    for (int i = 0; i < dim/2; i++) {
        for (int j = 0; j < dim/2; j++) {
            delete [] kernel[i][j];
        }
        delete [] kernel[i];
    }
    delete [] kernel;

    return 0;
}

int doGaussianSmoothDensity (double ** & halo_data, int num_lines, int num_fields, int x, int y) {
 
    int rank = 3, dim = 1024;

    double *** data;

    std::string snap = std::to_string((long long)id);

    if (id >= 100) snap = "0"+snap;
    if (id < 100) snap = "00"+snap;
    if (id < 10)  snap = "0"+snap;

    std::string file_name = "/pfs/chtlee/bolshoi_plank/DougsCICs/BP_"+snap+"_densities_1024";

    doBinaryRead3D (file_name+".bin",data,dim);
    
    for (int i = 0; i < 20; i++) {
        std::cout << data[0][0][i] << " ";
    }
    std::cout << std::endl;

    double *** data2 = new double ** [dim];
    for (int i = 0; i < dim; i++) {
        data2[i] = new double * [dim];
        for (int j = 0; j < dim; j++) {
            data2[i][j] = new double [dim];
            for (int k = 0; k < dim; k++) data2[i][j][k] = 0.0;
        }   
    } 

    double * kernel;

    double HWHM = 0.5; //units of Mpc/h (smoothing radius)

    // use commandline arg for smoothing length if specified
    if (slen) HWHM = slen;

    // 1 Mpc h-1 FWHM
    genGaussianKernel (kernel, 1024, 2*HWHM*1024./250.);

    // this will exlude a box of length 1 mpc/h ^3 around each point.
    // i think the remaining area will follow a gaussian contribution
    // with a FWHM of 4 mpc/h and a cutof threshold of 1e-8
    //genGaussianKernelExcludeCenter (kernel, 1024, 16.384, 1e-8, 4.096);

    int kernel_length = 0;

    // find kernel relevance length
    for (int i = 0; i < (int)ceil((float)dim/2.0); i++) {
        if (kernel[i] == 0.0) {
            kernel_length = i;
            break;
        }
    }

    std::cout << "Starting Convolution" << std::endl;
    std::cout << "kernel length: " << kernel_length << std::endl;

#pragma omp parallel
{
    if (omp_get_thread_num() == 0) {
        std::cout << "Number of omp threads/procs: " << omp_get_num_threads() << "/" << omp_get_num_procs() << std::endl;
    }
}

    double sum = 0.0;
    int length;

    std::ofstream progress("progress_"+std::to_string((long long)id)+"_"+std::to_string((long long)slen),std::ios::out);

#pragma omp parallel for schedule(dynamic) reduction (+:sum)
    for (int i = 0; i <dim; i++) {
        for (int j = 0; j < dim; j++) {
            for (int k = 0; k < dim; k++) {
                sum += data[i][j][k];
            }
        }
    }

    double avg = sum/double(dim*dim*dim);

    std::cout << "The average density of the universe is: " << sum << "/" << dim*dim*dim << " = " << avg << std::endl;

    // do the first dimension of convolution
#pragma omp parallel for schedule(dynamic) private(sum,length) shared(data,data2)
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            if (!omp_get_thread_num()) {
                if (!j) {
                    progress << "1: " << i << " ";
                    progress.flush();
                }
            }
            for (int k = 0; k < dim; k++) {
                sum = 0.0;
                for (int l = 0; l < dim; l++) {

                    length = abs(k-l);

                    // handle periodic boundary
                    if (length > dim/2) length = dim - length;

                    if (length <= kernel_length) sum += data[i][j][l]*kernel[length];
                }
                data2[k][i][j] = sum;
            }
        }
        if (!omp_get_thread_num()) progress << std::endl;
    }

    // do the second dimension of convolution
#pragma omp parallel for schedule(dynamic) private(sum,length) shared(data,data2)
    for (int k = 0; k < dim; k++) {
        for (int i = 0; i < dim; i++) {
            if (!omp_get_thread_num()) {
                if (!i) {
                    progress << "2: " << k << " ";
                    progress.flush();
                }
            }
            for (int j = 0; j < dim; j++) {
                sum = 0.0;
                for (int l = 0; l < dim; l++) {

                    length = abs(j-l);

                    // handle periodic boundary
                    if (length > dim/2) length = dim - length;

                    if (length <= kernel_length) sum += data2[k][i][l]*kernel[length];
                }
                data[j][k][i] = sum;
            }
        }
        if (!omp_get_thread_num()) progress << std::endl;
    }

    // do the third dimension of convolution
#pragma omp parallel for schedule(dynamic) private(sum,length) shared(data,data2)
    for (int j = 0; j < dim; j++) {
        for (int k = 0; k < dim; k++) {
            if (!omp_get_thread_num()) {
                if (!k) {
                    progress << "3: " << j << " ";
                    progress.flush();
                }
            }
            for (int i = 0; i < dim; i++) {
                sum = 0.0;
                for (int l = 0; l < dim; l++) {

                    length = abs(i-l);

                    // handle periodic boundary
                    if (length > dim/2) length = dim - length;

                    if (length <= kernel_length) sum += data[j][k][l]*kernel[length];
                }
                // divide out average density of universe
                data2[i][j][k] = sum/avg;
            }
        }
        if (!omp_get_thread_num()) progress << std::endl;
    }

    progress << "writing output file" << std::endl;

    std::ofstream outfile(file_name+"_"+std::to_string((long long)HWHM)+".bin",std::ios::out);
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            outfile.write((const char *)data2[i][j],dim*sizeof(double));
        }
    }

    outfile.close();

    progress.close();

    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            delete [] data[i][j];
            delete [] data2[i][j];
        }
        delete [] data[i];
        delete [] data2[i];
    }
    delete [] data;
    delete [] data2;

    delete [] kernel;

    return 0;
}

int doBuildDensityCatalog (double ** & halo_data, int num_lines, int num_fields, int x, int y) {

    double *** data_bin;

    double ** data;
    int num_lines_2=0, num_fields_2=0;

    int dim = 1024;
    double box_size = 250.0;

    std::string snap = std::to_string((long long)id);

    if (id >= 100) snap = "0"+snap;
    if (id < 100) snap = "00"+snap;
    if (id < 10)  snap = "0"+snap;

    std::string file_name = "/zang/chtlee/pfs_backup/bolshoi_plank/DougsCICs/BP_"+snap+"_densities_1024";

    if (slen) file_name += "_" + std::to_string((long long) slen); 

    doBinaryRead3D (file_name+".bin",data_bin,dim);

    //readFile("hlist_1.00230.list.updated.old",num_lines_2,num_fields_2,data);

    // lets output a slice
    std::ofstream slice(file_name+"_slice.dat",std::ios::out);

    for (int i = 0; i < 1024; i++) {
        for (int j = 0; j < 1024; j++) {
            slice << i << " " << j << " " << data_bin[i][j][0] << std::endl;
        }
    }

    slice.close();
    
    // update catalog (assuming input file is old catalog)
    //std::ofstream outfile("hlist_1.00230.list.testnew",std::ios::out);
    std::ofstream outfile(file_name+"_hlist.bin");

    std::cout << "Building c_nfw - density relation" << std::endl;

    unsigned int central_num = 0;
    double density[num_lines];
    long chunk_size = 100000;

#pragma omp parallel for schedule(dynamic,chunk_size) shared(halo_data,data_bin,num_lines,box_size,dim)
    // lets start running through the halo data
    for (int i = 0; i < num_lines; i++) {

        if (!omp_get_thread_num()) {
            if (i % 100000 == 0) std::cout << i/100000 << std::endl;
        }

        // lets start by just looking at centrals
        int pid = halo_data[5][i];
        if (true/*pid == -1*/) {
            double c_nfw = halo_data[11][i]/halo_data[12][i];

            // find local density
            int position[3];
            //std::cout << i << ": " << "(" << halo_data[17][i] << " , " << halo_data[18][i] << " , " << halo_data[19][i] << ") : (";
            for (int j = 0; j < 3; j++) {
                position[j] = (int)floor(halo_data[17+j][i]/(box_size/(float)dim));
                if (position[j] >= dim || position[j] < 0) {
                    std::cout << "Out of bounds("<<i<<"): " << position[j] << std::endl;
                    if (position[j] >= dim) position[j] = dim-1;
                    else if (position[j] < 0) position[j] = 0;
                }
                //std::cout << position[j] << " , ";
            }
            //std::cout << ")" << std::endl;

            //double density = data_bin[position[0]][position[1]][position[2]];
            density[i] = data_bin[position[0]][position[1]][position[2]];
        }
    }

    std::cout << "Writing output file" << std::endl;

    outfile << "# NUM_LINES " << num_lines << std::endl;
    //outfile << "# NUM_FIELDS " << num_fields+1 << std::endl;
    outfile << "# NUM_FIELDS " << 1 << std::endl;


//    for (int i = 0; i < num_fields; i++) {
//        std::cout << i << std::endl;
//        outfile.write((const char *)halo_data[i],num_lines*sizeof(double));
//    }

    outfile.write((const char *)density,num_lines*sizeof(double));

    outfile.close();
/*
    for (int i = 0; i < num_lines; i++) {

        //unsigned int halo_id = halo_data[1][i];
        //double mvir = halo_data[10][i];
        //outfile << halo_id << " " << mvir << " " << c_nfw << " " << density_half << " " << density_1 << " " << density_2 << " " << density_4 << std::endl;
        if (i % 100000 == 0) std::cout << i/100000 << std::endl;

        outfile << std::setprecision(10);
        for (int j = 0; j < num_fields; j++) {
            if (!j) {
                outfile << halo_data[0][i];
            }
            else {
                if (j == 1)
                    outfile << " " << (unsigned int)halo_data[j][i];
                else
                    outfile << " " << halo_data[j][i];
            }
        }
        //outfile << " " << c_nfw;
        //outfile << " " << density << std::endl;
        //for (int j = 63; j < 70; j++) {
        //    outfile << " " << data[j][i]/8.0;
        //}
        outfile << " " << density[i];
        outfile << std::endl;
        //outfile << halo_id << " " << mvir << " " << c_nfw << " " << data[3][central_num] << " " << data[4][central_num] << " " << density_2 << " " << density_4 << std::endl;
        //central_num++;
    }

    outfile.close();
*/
    //free (data, num_fields_2);
    free (data_bin, dim, dim);

    return 0;
}

int doCalcVoxelDensityDistributions (int x, int y, int num_bins, bool NORM) {

    double *** data_bin;
    int dim = 1024;

    std::string d_str = std::string("");
    std::string z_str = std::string("");

    switch (x) {
        case 0:
            d_str += "0";
            break;
        case 1:
            d_str += "1";
            break;
        case 2:
            d_str += "2";
            break;
        case 3:
            d_str += "4";
            break;
        case 4:
            d_str += "8";
            break;
        case 5:
            d_str += "16";
            break;
        // 32 no longer used
        case 6:
            d_str += "32";
            break;
    }

    switch(y) {
        case 0:
            z_str += "214";
            break;
        case 1:
            z_str += "170";
            break;
        case 2:
            z_str += "136";
            break;
        case 3:
            z_str += "103";
            break;
    }

    std::string inputfile = "../bolshoi_plank/DougsCICs/BP_0" + z_str + "_densities_1024_" + d_str + ".bin";

    if (output.length()) parseOutputFileName(inputfile); 

    doBinaryRead3D (inputfile,data_bin,dim);

//    double sum = 0;
//
//#pragma omp parallel for schedule(dynamic) reduction (+:sum)
//    for (int i = 0; i <dim; i++) {
//        for (int j = 0; j < dim; j++) {
//            for (int k = 0; k < dim; k++) {
//                sum += data_bin[i][j][k];
//            }
//        }
//    }
//
//    double avg = sum/double(dim*dim*dim);
//
//    std::cout << "The average density of the universe is: " << sum << "/" << dim*dim*dim << " = " << avg << std::endl;

    doHist1DVoxels (data_bin, dim, num_bins, NORM);

    // for now, we're going to cannibalize this function to compute the effective volume of
    // regions of different density
//    std::vector<double> bin_edge {-0.175,-0.05,0.05,0.175};
//    std::vector<int> bin_counts (bin_edge.size()+1, 0);
//
//    #pragma omp parallel shared(bin_edge,bin_counts)
//    {
//
//        std::vector<int> bin_counts_priv (bin_counts.begin(),bin_counts.end());
//
//        #pragma omp for schedule(dynamic)
//        for (int i = 0; i < dim; i++) {
//            if (!omp_get_thread_num()) std::cout << i << " ";
//            for (int j = 0; j < dim; j++) {
//                for (int k = 0; k < dim; k++) {
//                    for (int l = 0; l < bin_edge.size(); l++) {
//
//                        // check from left most bin to right most bin
//                        if (log10(data_bin[i][j][k]) < bin_edge[l]) {
//                            bin_counts_priv[l]++;
//                            break;
//                        }
//                        // if we made it past all edges, must be highest bin (no right edge provided)
//                        if (l+1 == bin_edge.size()) bin_counts_priv[l+1]++;
//                    }
//                }
//            }
//        }
//
//        #pragma omp critical
//        {
//            for (int i = 0; i < bin_counts.size(); i++) {
//                bin_counts[i] += bin_counts_priv[i];
//            }
//        }
//
//    } // end parallel
//
//    std::cout << std::endl;
//
//    std::cout << "Volume fraction in regions with density:" << std::endl;
//    for (int i = 0; i < bin_edge.size(); i++) {
//        if (!i) std::cout << "< ";
//        else std::cout << bin_edge[i-1] << " - ";
//        std::cout << bin_edge[i] << ": " << ((double)bin_counts[i])/(1024*1024*1024.) << std::endl;
//    }
//    std::cout << "> " << bin_edge.back() << ": " << ((double)bin_counts.back())/(1024*1024*1024.) << std::endl;

    return 0;
}

int doCombineHlistFiles (double ** & halo_data, int num_lines, int num_fields) {

    int num_files = 7;

    double new_data[num_files][num_lines];

    std::string snap = std::to_string((long long)id);

    if (id >= 100) snap = "0"+snap;
    if (id < 100) snap = "00"+snap;
    if (id < 10)  snap = "0"+snap;

    std::string hlist_file = "/zang/chtlee/pfs_backup/bolshoi_plank/DougsCICs/BP_"+snap+"_densities_1024";

    for (int i = 0; i < num_files; i++) {

        double ** hlist_data;
        int num_lines_2=0, num_fields_2=0;
        //std::string hlist_file = std::string("hlist_1.00230.list.testnew_"+std::to_string((long long)i));

        std::string suffix("");
        if (i) suffix = "_" + std::to_string((long long )pow(2.,(i-2.)));
        suffix += "_hlist.bin";

        readBinaryFile(hlist_file+suffix,num_lines_2,num_fields_2,hlist_data);

        if (num_lines_2 != num_lines)
            std::cerr << "num_lines_2 != num_lines in combinedHlistFiles at i = " << i << std::endl;

        for (int j = 0; j < num_lines_2; j++) {
            // check if halo id is the same
            //if (halo_data[1][j] != hlist_data[0][j]) {
            //    std::cerr << "halo id mismatch at j = " << j << std::endl;
            //}

            new_data[i][j] = hlist_data[0][j];
        }

        free(hlist_data, num_fields_2);
    }

    std::cout << "new_data test:" << std::endl;
    std::cout << new_data[0][0] << " " << new_data[1][100] << " " << new_data[5][12031] << std::endl;

    std::ofstream outfile;
    outfile.open(output.c_str(),std::ofstream::out);

    std::cout << "Writing output file \"" << output << "\"" << std::endl;
    std::cout << "num_lines: " << num_lines << std::endl;
    std::cout << "num_fields: " << num_fields+num_files << std::endl;

    outfile << "# NUM_LINES " << num_lines << std::endl;
    outfile << "# NUM_FIELDS " << num_fields+num_files << std::endl;

    int dot_freq = int(ceil((num_fields+num_files)/10.));

    for (int i = 0; i < num_fields+num_files; i++) {
        if (i % dot_freq == 0) {
            std::cout << ".";
            std::cout.flush();
        }
        if (i < num_fields)
            outfile.write((const char *)halo_data[i],num_lines*sizeof(double));
        else
            outfile.write((const char *)new_data[i-num_fields],num_lines*sizeof(double));
    }

    std::cout << std::endl;

    //setOutputFileName();

    //std::ofstream outfile("hlist_1.00230.list.testnew",std::ios::out);

    //for (int i = 0; i < num_lines; i++) {

    //        for (int j = 0; j < num_fields; j++) {
    //            if (!j) {
    //                outfile << halo_data[0][i];
    //            }
    //            else {
    //                if (j == 1)
    //                    outfile << " " << (unsigned int)halo_data[j][i];
    //                else
    //                    outfile << " " << halo_data[j][i];
    //            }
    //        }

    //        for (int j = 0; j < num_files; j++) {
    //            outfile << " " << new_data[j][i];
    //        }

    //        outfile << std::endl;
    //}

    //outfile.close();
}

int doDensityPercentileEvolution (double ** &data, int num_lines, int num_fields, int x) {

    // data is MMPB merger tree
    // x is 0-4 corresponding to redshifts 0,0.5,1,2,3,4

    std::string catalog;

    // determine which redshift we want to output based on x parameter
    switch (x) {
        case 0: catalog = "bp_z0_centrals_cut10.bin";     break;
        case 1: catalog = "bp_z0.5_centrals_cut10.bin";   break;
        case 2: catalog = "bp_z1_centrals_cut10.bin";     break;
        case 3: catalog = "bp_z2_centrals_cut10.bin";     break;
        case 4: catalog = "bp_z3_centrals_cut10.bin";     break;
        case 5: catalog = "bp_z4_centrals_cut10.bin";     break;
        default:
                std::cout << "WTF is x??" << std::endl;
    }

    // read in desired redshift catalog
    double ** cat_data;
    int num_fields_cat=0, num_lines_cat=0;

    readBinaryFile(catalog,num_lines_cat,num_fields_cat,cat_data);

    // let's now collect the halo ids from the merger trees for the chosen redshift
    double a = cat_data[0][0];  // scale factor of chosen catalog. use to identify
                                // correct halos to extract from merger trees

    std::cout << "scale factor is: " << a << std::endl;

    std::vector<long int> ids;  // container to hold ids of progenitor halos from chosen redshift.
    std::vector<double> mass;   // container to hold masses of progenitor halos

    for (int i = 0; i < num_lines; i++) {
        if (data[i][0] == 0.67325) {              // remember these are bin2 format, row major order.
            ids.push_back(data[i][1]);
            mass.push_back(data[i][10]);
        }
    }

    std::cout << "Found " << ids.size() << " halos in merger tree file at this redshift." << std::endl;

    // now that we know something about the halos we're interested in, let's
    // make a mass cut on the catalog data to simplify things.

    // range of mass percentiles to include from merger tree halos, symmetric about
    // the median (i.e. 0.9 would be the 5th - 95th percetiles used).
    double mass_cut_percentile_range = 0.9;

    std::sort(mass.begin(),mass.end());

    std::cout << "Median mass is " << log10(mass[(int)mass.size()/2.]) << std::endl;
    std::cout << "Full mass range is " << log10(mass[0]) << " - " << log10(mass[mass.size()-1]) << std::endl;

    double low_mass_cut = mass[(int)floor((1.-mass_cut_percentile_range)/2.*mass.size())];
    double high_mass_cut = mass[(int)ceil((mass_cut_percentile_range/2.+0.5)*mass.size())];

    std::vector<std::vector<double>> subset;

    for (int i = 0; i < num_lines_cat; i++) {
        
        if ((long int)cat_data[1][i] == ids[0]) std::cout << "found 0th id" << std::endl;

        // find all halos within this mass range
        if (cat_data[10][i] > low_mass_cut && cat_data[10][i] < high_mass_cut) {

            // once we've found one, copy all its properties and insert into new container
            std::vector<double> h;

            for (int j = 0; j < num_fields_cat; j++) h.push_back(cat_data[j][i]);

            subset.push_back(h);
        }
    }

    std::cout << "Found " << subset.size() << " halos in mass range (" << log10(low_mass_cut) << ", " << log10(high_mass_cut) << ")" << std::endl;

    struct sort_comp {
        bool operator()(const std::vector<double> &left, const std::vector<double> &right) {
            return (left[1] < right[1]);
        }
    };

    std::cout << "Sorting...";  std::cout.flush();

    std::sort(subset.begin(),subset.end(),sort_comp());

    std::cout << " done." << std::endl;

    return 0;
}
