#include "profiles.h"

int readProfileFile(std::string name, std::vector<halo_profile> & halos);
int readProfileBinaryFile(std::string name, std::vector<halo_profile> & halos);
int writeProfileBinaryFile (std::vector<halo_profile> & halos);
int doFitProfiles (std::vector<halo_profile> & halos);

int doProfileAnalysis (std::string name, bool binary_read) {

    std::vector<halo_profile> halos;

    if (binary_read) {
        if (readProfileBinaryFile(fileName, halos)) {
            return -1;
        }
    }
    else {
        if (readProfileFile(fileName, halos)) {
            return -1;
        }

        writeProfileBinaryFile(halos);
    }
    

    // subset these halos
    //std::vector<halo_profile *> hpp;
    //
    //for (int i = 0; i < halos.size(); i++) {

    //    hpp.push_back(&halos[i]);
    //}

    //std::sort(hpp.begin(),hpp.end(),[](halo_profile * left, halo_profile * right){return left->mass > right->mass;});

    //std::vector<halo_profile> output_halos;
    //int num_output_halos = 20;

    //for (int i = 0; i < num_output_halos; i++) {
    //    output_halos.push_back(*(hpp[(int)floor((halos.size()-1)/(double)num_output_halos)*i]));
    //}

    //writeProfileBinaryFile(output_halos);

    doFitProfiles(halos);
};

int readProfileFile(std::string name, std::vector<halo_profile> & halos) {
	std::ifstream inpt;
	std::string tmp_str;
	char tmp_word_c[MAX_WORD];
	char tmp_line_c[MAX_LINE];
    std::string tmp_line(tmp_line_c);
    std::string tmp_word(tmp_word_c);
	double tmp_d;

    std::cout << "Reading file: \"" << name << "\"" << std::endl;

    inpt.open(name.c_str(),std::ifstream::in);

    int startline = 1;

	// lets check if this input file has a header
	// and skip it if it does
	int position = inpt.tellg();
	while (getline(inpt,tmp_str)) {
        startline++;
		if (tmp_str[0] != '#' && tmp_str[0] != 'M') {
            startline--;
			inpt.seekg(position);
			break;
		}
		position = inpt.tellg();
	}

	// lets scan through file and see how many lines of data we have
	inpt.seekg(0,inpt.end);
	long int file_eof = inpt.tellg();
	inpt.seekg(position);

	// finally, read file into memory
	int i = 0;
    int skipped_lines = 0;
    int dot_freq = int(ceil((file_eof-position)/10.));
    //std::cout << file_eof << " " << file_eof-position << " " << dot_freq << std::endl;

    int progress_bar = -1;

    //std::vector<halo_profile> halos;

    while (inpt.good()) {

        double current_progress = floor(((long int)inpt.tellg()-position)/(double)dot_freq);
        if (current_progress > progress_bar) {
            std::cout << (int)++progress_bar << ".";
            std::cout.flush();
        }

        halo_profile hp;

        std::getline(inpt, tmp_line);
        std::istringstream line_stream(tmp_line);
        
        line_stream >> hp.id >> hp.parent_id >> hp.mass >> hp.radius >> hp.vmax >> hp.x >> hp.y >> hp.z >> hp.vx >> hp.vy >> hp.vz >> hp.np >> hp.n_rbins >> hp.ppbin;

        //if (i < 20) std::cout << "Reading profile of halo " << hp.id << " with mass " << log10(hp.mass) << " and with " << hp.n_rbins << " radial bins." << std::endl;

        int j = 0;
        while (j++ < hp.n_rbins) {
            double bin;
            line_stream >> bin;
            hp.bins.push_back(bin);
        }

        halos.push_back(hp);
        
        i++;
    }

    std::cout << " complete." << std::endl;

    std::cout << "Read in " << halos.size() << " halos." << std::endl;

    inpt.close();

    return 0;
}

int readProfileBinaryFile(std::string name, std::vector<halo_profile> & halos) {

    std::cout << "Reading binary file: \"" << fileName << "\"" << std::endl;

    std::ifstream inpt(fileName,std::ios::in);

    long int header[1] = {0};

    // first read all halo profile fields

    inpt.read((char *)header,sizeof(long int));

    std::cout << "Reading " << header[0] << " halo profiles." << std::endl;

    int static_profile_size = sizeof(halo_profile)-2*sizeof(std::vector<double>)-sizeof(new_profile_fields);

    for (int i = 0; i < header[0]; i++) {

        halo_profile hp;

        inpt.read((char *)&hp, static_profile_size);

        halos.push_back(hp);
    }

    double buffer[header[0]];

    long int position = inpt.tellg();
	inpt.seekg(0,inpt.end);
	long int file_eof = inpt.tellg();
	inpt.seekg(position);
    int dot_freq = int(ceil((file_eof-position)/10.));

    int progress_bar = -1;

    while (true) {
        
        double current_progress = floor(((long int)inpt.tellg()-position)/(double)dot_freq);
        if (current_progress > progress_bar) {
            std::cout << (int)++progress_bar << ".";
            std::cout.flush();
        }

        inpt.read((char *)header,sizeof(long int));

        if (header[0] == 0) break;

        inpt.read((char *)buffer,header[0]*sizeof(double));

        for (int i = 0; i < header[0]; i++) {

            halos[i].bins.push_back(buffer[i]);
        }
    }

    std::cout << " complete." << std::endl;

    return 0;
}

int writeProfileBinaryFile (std::vector<halo_profile> & halos) {

    std::ofstream of;
    
    setOutputFileName();
    
    of.open(output.c_str(),std::ofstream::out);

    std::cout << "Writing output file " << output << std::endl;

    // first, sort profiles by number of radial bins

    std::vector<halo_profile *> hpp;
    
    for (int i = 0; i < halos.size(); i++) {

        hpp.push_back(&halos[i]);
    }

    std::sort(hpp.begin(),hpp.end(),[](halo_profile * left, halo_profile * right){return left->n_rbins > right->n_rbins;});

    long int header[1] = {(long int)halos.size()};

    // write how many halos we will be reading
    of.write((const char *)header, sizeof(long int));

    std::cout << "Writing halo profile fields... ";
    std::cout.flush();

    int static_profile_size = sizeof(halo_profile)-2*sizeof(std::vector<double>)-sizeof(new_profile_fields);

    // now go through all halos and output data starting with fields, then radial bins
    for (int i = 0; i < hpp.size(); i++) {

        // write all fields except for vector of bins for each halo profile
        of.write((const char *)hpp[i],static_profile_size);
    }

    std::cout << "complete" << std::endl;

    int current_bin = 0;
    
    double buffer[halos.size()];

    //long int header_pos = of.tellp();

    bool done = false;

    // fill in header junk data -- will be rewritten with correct value once determined
    //of.write((const char *)header,sizeof(long int));

    std::cout << "Writing radial bins ";
    std::cout.flush();

    int dot_freq = (int)ceil(hpp[0]->bins.size()/10.);

    while (true) {

        if (current_bin % dot_freq == 0) {
            std::cout << current_bin/dot_freq << ".";
            std::cout.flush();
        }

        // reset to zero and count how many halos have appropriate number of bins
        header[0] = 0;

        int i;

        for (i = 0; i < hpp.size(); i++) {

            // remember, halos already sorted by bin number, so when we find first one
            // with fewer bins than the current output bin number, we know we're done searching
            if (hpp[i]->bins.size()-1 >= current_bin) {

                buffer[i] = hpp[i]->bins[current_bin];
                header[0]++;
            }

            else {
                if (i == 0) done = true;
                break;
            }
        }

        if (done) {
            // write 0 to let us know we are at the end of the file:
            header[0] = 0;
            of.write((const char *)header, sizeof(long int));
            break;
        }

        of.write((const char *)header, sizeof(long int));
        of.write((const char *)buffer, header[0]*sizeof(double));
        current_bin++;
    }

    std::cout << " complete." << std::endl;

    of.close();

    return 0;
}

//====================================================
// doStackProfiles:
// 
// data:
// num_lines:
//          number of lines in the catalog.
// num_fields:
//          number of fields for each line (halo).
//====================================================
int doStackProfiles (double ** &data, int num_lines, int num_fields) {

    // data is the subset halo catalog for us to look up and stack profiles

    // profile fields
    int RVIR_P = 3, NBINS = 12, PPB = 13, R1 = 15;          

    // halo catalog fields
    int A = 0, MVIR = 10, RVIR = 11, RS = 12, ORIG_ID = 30, ID = 1;    

    // what is the largest halo we have?
    std::vector<double> rvir;
    std::vector<double> rs;
    std::vector<double> cnfw;
    std::vector<double> rmin;
    std::vector<double> rmax;

    // let's first see what redshift we are working with
    std::string z("");

    if      (data[0][0] == 1.00230)     z += '0';
    else if (data[0][0] == 0.50112)     z += '1';

    //std::string profiles_path("/zang/chtlee/pfs_backup/bolshoi_plank/profiles/profiles_z"+z+".mod.bin");
    std::string profiles_path("/zang/chtlee/pfs_backup/bolshoi_plank/profiles/extended_profiles.mod.bin");

    int num_lines_p=0, num_fields_p=0;
    double ** data_p;

    readBinaryFile (profiles_path,num_lines_p,num_fields_p,data_p);

    // create a map so we can look up halos by orig halo id
    std::unordered_map <long int, long int> halo_map;
    std::unordered_map<long int, long int>::const_iterator fetch;

    // for extended profiles, first field is nbody id
    for (int i = 0; i < num_lines_p; i++) {
        halo_map[data_p[0][i]] =  i;
    }

    // now lets get a sense for max r
    for (int i = 0; i < num_lines; i++) {

        int j = 0;

        // find halo id in locations data
        fetch = halo_map.find(data[ID][i]);

        if (fetch == halo_map.end()) {
            //std::cerr << "Halo " << data[ID][j] << " not found at j " << j << std::endl;
            continue;
        }

        else j = fetch->second;

        rmax.push_back(data_p[R1+(int)data_p[NBINS][j]-1][j]);
        rmin.push_back(data_p[R1][j]);
    }

    double a = data[A][0];

    for (int i = 0; i < num_lines; i++) {
        rvir.push_back(a*data[RVIR][i]);
        rs.push_back(a*data[RS][i]);
        cnfw.push_back(data[RVIR][i]/data[RS][i]);
    }

    std::sort(rvir.begin(),rvir.end());
    std::sort(rs.begin(),rs.end());
    std::sort(cnfw.begin(),cnfw.end());
    std::sort(rmin.begin(),rmin.end());
    std::sort(rmax.begin(),rmax.end());

    double rvirmax = rvir.back();
    double rmaxmax = rmax.back();

    std::cout << "median rvir: " << vmedian(rvir) <<  " " << rvir[(int)(0.2*rvir.size())] << " " << rvir[(int)(0.8*rvir.size())] << std::endl;
    std::cout << "median rs: " << vmedian(rs) << " " << rs[(int)(0.2*rs.size())] << " " << rs[(int)(0.8*rs.size())] << std::endl;
    std::cout << "median cnfw: " << vmedian(cnfw) << std::endl;
    std::cout << "median rmin: " << vmedian(rmin) <<  " " << rmin[(int)(0.2*rmin.size())] << " " << rmin[(int)(0.8*rmin.size())] << std::endl;
    std::cout << "median rmax: " << vmedian(rmax) <<  " " << rmax[(int)(0.2*rmax.size())] << " " << rmax[(int)(0.8*rmax.size())] << std::endl;

    // how many bins to we want to use for profiles?
    int nbins = 2000;
    //double bin_min = log10(rmaxmax/nbins);
    double bin_min = log10(rmin.front());
    double bin_max = log10(rmaxmax);
    double bin_size = (bin_max-bin_min)/nbins;
    //double bin_size = (log10(rmaxmax)-bin_min)/nbins;

    // recenter bin_min
    bin_min += bin_size/2.;

    std::cout << "rvirmax: " << rvirmax << std::endl;
    std::cout << "bin_min: " << bin_min << ", bin_max: " << bin_max << std::endl;
    std::cout << "rmaxmax: " << rmaxmax << std::endl;
    std::cout << "bin_size: " << bin_size << std::endl;

    std::vector<std::vector<double>> profiles;

    for (int i = 0; i < nbins; i++) {
        std::vector<double> bin;
        profiles.push_back(bin);
    }

    // now lets go through and look up profiles
    // keep track of halo id index in locations data
    int i = 0;

    //bool output_single_halo = true;
    bool output_single_halo = false;

    std::ofstream of;
    setOutputFileName();
    of.open(output.c_str(),std::ofstream::out);

    int dot_freq = (int) ceil(num_lines/10.);

    int miss_count = 0;

    int it_num = -1;

    // loop over z = 0 halos in input catalog
    for (int j = 0; j < num_lines; j++) {

        if (output_single_halo) {

            it_num++;

            if (it_num >= 1 && miss_count < it_num) break;
        }

        if (j % dot_freq == 0) {
            std::cout << ".";
            std::cout.flush();
        }

        // find halo id in locations data
        fetch = halo_map.find(data[ID][j]);

        if (fetch == halo_map.end()) {
            //std::cerr << "Halo " << data[ID][j] << " not found at j " << j << std::endl;
            miss_count++;
            continue;
        }

        else i = fetch->second;

        // build vector of first radial bin position for each halo
        //rmin.push_back(data_p[R1][i]);

        //for (int k = 0; k < num_fields_p; k++) {
        //    std::cout << data_p[k][i] << " ";
        //}

        //std::cout << std::endl;

        std::ofstream of2;
        std::ofstream of3;
        if (output_single_halo) {
            of2.open((output+"2").c_str(),std::ofstream::out);
            of3.open((output+"3").c_str(),std::ofstream::out);

            // lets print some info about the halo we found
            std::cout << "Halo Mass: " << log10(data[MVIR][j]) << std::endl;
            std::cout << "Rvir: " << log10(data[RVIR][i]) << ", RS: " << log10(data[RS][i]) << std::endl;
            std::cout << "MAX R: " << data_p[R1+(int)data_p[NBINS][i]-1][i] << std::endl;
        }

        // natural cubic spline interpolation implementation

        // determine x,y node vectors
        std::vector<double> x;
        std::vector<double> y;
        double center = 0;
        double rho = 0;

        int k = 0;

        x.push_back(0);
        y.push_back(0);

        for (k = R1; k < num_fields_p; k++) {
            if (data_p[k][i] == 0.0) break;

            // fill up log R - M node vectors 
            x.push_back(a*data_p[k][i]);
            y.push_back(y.back()+data_p[PPB][i]);
            if (output_single_halo) of2 << x.back() << " " << y.back() << std::endl;
        }

        // still need to do final node
        //x.push_back(a*data_p[RVIR_P][i]);
        //y.push_back(y.back()+data_p[PPB][i]);
        //if (output_single_halo) of2 << x.back() << " " << y.back() << std::endl;

        if (output_single_halo) of2.close();

        // now get the spline coefficients
        std::vector<spline> s = cubicSpline(x,y);

        int scount = 0;
        double prev_m = 0;
        bool output3 = false;

        //std::cout << "max r: " << data_p[R1+(int)data_p[NBINS][i]-1][i] << std::endl;

        for (k = 0; k < nbins; k++) {

            // bin position
            double pos = bin_min+k*bin_size;

            int nbins_p = data_p[NBINS][i];

            // comparing (log scale) bin position to profile r
            if (pos >= log10(a*data_p[R1+nbins_p-1][i])) break;

            // need to convert bin position back to linear scaling for spline comparison
            pos = pow(10.,pos);

            if (scount < s.size()-1) {
                while (pos > s[scount+1].x) {
                    scount++;
                    output3 = true;

                    // break if we reach end of spline
                    if (scount >= s.size()-1) break;
                }
            }

            // distance from bin position to spine node
            double h = pos-s[scount].x;

            //rho = s[scount].a + s[scount].b * h + s[scount].c * h*h + s[scount].d * h*h*h;
            double m = s[scount].a + s[scount].b * h + s[scount].c * h*h + s[scount].d * h*h*h;

            // need to probe below bin_min to get proper prev_m value
            if (k == 0) {
                double h0 = pow(10.,bin_min-bin_size);
                prev_m = s[scount].a + s[scount].b * h0 + s[scount].c * h0*h0 + s[scount].d * h0*h0*h0;
            }

            rho = (m-prev_m)/((4.*M_PI/3.)*(pow(pos,3.)-pow(pow(10.,bin_min+(k-1)*bin_size),3.)));

            // keeps track of previous mass value, for computing rho
            prev_m = m;

            if (output3 && output_single_halo) {
                output3 = false;
                of3 << s[scount].x << " " << rho << std::endl;
            }

            profiles[k].push_back(rho);
        }

        //for (k = 15; k < num_fields_p; k++) {
        //    if (data_p[k][i] == 0.0) break;

        //    // fill up log R - log Rho node vectors 
        //    center = (a*data_p[k][i]+a*data_p[k-1][i])/2.;
        //    x.push_back(log10(center));
        //    // assume r^-1 profile for first bin
        //    if (k == 15) {
        //        rho = data_p[13][i]/(2.*M_PI*pow(a*data_p[k][i],2.));
        //    }
        //    else {
        //        rho = data_p[13][i]/((4.*M_PI/3.)*(pow(a*data_p[k][i],3.)-pow(a*data_p[k-1][i],3.)));
        //    }
        //    y.push_back(log10(rho));
        //    if (output_single_halo) of2 << log10(center) << " " << log10(rho) << std::endl;
        //}

        //// still need to do final node
        //center = (a*data_p[3][i]+a*data_p[k-1][i])/2.;
        //x.push_back(log10(center));
        //rho = data_p[13][i]/((4.*M_PI/3.)*(pow(a*data_p[3][i],3.)-pow(a*data_p[k-1][i],3.)));
        //y.push_back(log10(rho));     
        //if (output_single_halo) of2 << log10(center) << " " << log10(rho) << std::endl;

        //if (output_single_halo) of2.close();

        //// now get the spline coefficients
        //std::vector<spline> s = cubicSpline(x,y);

        //int scount = 0;

        //// now fill in profile data
        //for (k = 0; k < nbins; k++) {

        //    double pos = bin_min+k*bin_size;

        //    if (pos >= log10(a*data_p[3][i])) break;

        //    // initially, we want to follow 0th spline from before node
        //    // all the way to next node
        //    if (scount == 0) {
        //        if (pos > s[1].x) {
        //            scount++;
        //        }
        //    }

        //    else if (scount < s.size()-1) {
        //        if (pos > s[scount+1].x) {
        //            scount++;
        //        }
        //    }

        //    double h = pos-s[scount].x;

        //    // if before first node, only use linear extrapolation from
        //    // 0th spline
        //    if (pos < s[0].x) {
        //        rho = s[0].a + s[0].b*h;
        //    }

        //    else {
        //        rho = s[scount].a + s[scount].b * h + s[scount].c * h*h + s[scount].d * h*h*h;
        //    }

        //    profiles[k].push_back(rho);
        //}

        // original implementation below

        //// initializing k here so final value can be saved to use after loop
        //int k = 0;

        //int rindex = 0;
        //int lindex = 0;
        //double norm = 0;

        //for (k = 15; k < num_fields_p; k++) {

        //    if (data_p[k][i] == 0.0) break;

        //    rindex = (int) floor((log10(a*data_p[k][i])-bin_min)/bin_size);

        //    if (data_p[k-1][i] == 0.0) lindex = 0;
        //    else lindex = (int) floor((log10(a*data_p[k-1][i])-bin_min)/bin_size);

        //    // num particles divided by shell volume
        //    norm = data_p[13][i]/((4.*M_PI/3.)*(pow(a*data_p[k][i],3.)-pow(a*data_p[k-1][i],3.)));

        //    //std::cout << norm << " ";
        //    //std::cout << rindex << " " << lindex << " " << data_p[13][i] << " " << a << " ";
        //    //std::cout.flush();

        //    if (lindex < 0 || rindex > nbins) {
        //        std::cout << "PROBLEM: L/R index out of bounds: " << lindex << ", " << rindex << std::endl;
        //    }

        //    for (int l = lindex; l < rindex; l++) {
        //        profiles[l].push_back(norm);
        //    }
        //}

        //// still need to do last bin
        //rindex = (int) floor((log10(a*data_p[3][i])-bin_min)/bin_size);
        //lindex = (int) floor((log10(a*data_p[k-1][i])-bin_min)/bin_size);

        //if (rindex > nbins-1) rindex = nbins-1;

        //// if we get messed up indices here, lets just skip to next halo
        //if (lindex < 0 || rindex > nbins) {
        //    std::cout << "PROBLEM: L/R index out of bounds: " << lindex << ", " << rindex << std::endl;
        //    std::cout << data_p[k-1][i] << std::endl;
        //    continue;
        //}

        //norm = data_p[13][i]/((4.*M_PI/3.)*(pow(a*data_p[3][i],3.)-pow(a*data_p[k-1][i],3.)));

        //for (int l = lindex; l < rindex; l++) {
        //    profiles[l].push_back(norm);
        //}
    }

    std::cout << "Number of halos not found: " << miss_count << std::endl;

    //std::sort(rmin.begin(),rmin.end());
    //std::cout << "median rmin: " << vmedian(rmin) <<  " " << rmin[(int)(0.2*rmin.size())] << " " << rmin[(int)(0.8*rmin.size())] << std::endl;

    for (int j = 0; j < nbins; j++) {
        if (profiles[j].size()) {

            std::sort(profiles[j].begin(),profiles[j].end());

            of << pow(10,bin_min+j*bin_size) << " " << vmedian(profiles[j]) << " ";
            of << vCIlo(profiles[j]) << " " << vCIhi(profiles[j]) << " ";
            of << vdisp(profiles[j],0.2) << " " << vdisp(profiles[j],0.8) << " ";
            of << profiles[j].size() << std::endl;
        }
    }

    of.close();

}

// mass enclosed by nfw profile 0->r, not including constant multiplier
double nfw_menc (double r, double rs) {

    return (log((rs+r)/rs) - (r/(rs+r)));
}

// full mass enclosed by nfw profile, in units of particle mass
double nfw_menc (double r, double rs, double c0) {

    return c0*pow(rs,3.)*nfw_menc(r,rs);
}

// Expecting parameters vector to have rs as first element and c0 as second element;
double chi2_nfw1 (halo_profile & hp, std::vector<double> & p) {

    if (p.size() != 1) {
        std::cout << "chi2_nfw passed " << p.size() << " parameters instead of 1." << std::endl;
    }

    double rs = p[0];
    double chi2 = 0;
    double last_menc = 0;
    double bin_mass = nfw_menc(hp.bins[hp.newf.nb_rvir-1.],rs)/(hp.newf.nb_rvir-1.);

    // note start from i = 0 because first bin will be at r = 0
    for (int i = 1; i < hp.bins.size(); i++) {

        // only want to fit up to the virial radius
        if (hp.bins[i] > hp.radius) break;
        
        double menc = nfw_menc(hp.bins[i],rs);
        double dx = hp.weights[i]*(menc-last_menc-bin_mass)/bin_mass;
        chi2 += dx*dx;
        last_menc = menc;
    }

    return chi2;
}

// Expecting parameters vector to have rs as first element and c0 as second element;
double chi2_nfw2 (halo_profile & hp, std::vector<double> & p) {

    if (p.size() != 2) {
        std::cout << "chi2_nfw passed " << p.size() << " parameters instead of 2." << std::endl;
    }

    double rs = p[0];
    double c0 = p[1];
    double chi2 = 0;
    double last_menc = 0;
    double bin_mass;// = nfw_menc(hp.bins[hp.newf.nb_rvir-1.],rs,c0)/(hp.newf.nb_rvir-1.);

    // note start from i = 0 because first bin will be at r = 0
    for (int i = 1; i < hp.bins.size(); i++) {

        // only want to fit up to the virial radius
        if (hp.bins[i] > hp.radius) break;
        
        double menc = nfw_menc(hp.bins[i],rs,c0);
        double dx = hp.weights[i]*(menc-last_menc-hp.ppbin)/hp.ppbin;
        chi2 += dx*dx;
        last_menc = menc;
    }

    return chi2;
}

// Expecting parameters vector to have rs as first element, c0 as second element, and gamma as third;
double chi2_nfw3 (halo_profile & hp, std::vector<double> & p) {

    if (p.size() != 3) {
        std::cout << "chi2_nfw passed " << p.size() << " parameters instead of 2." << std::endl;
    }

    double rs = p[0];
    double c0 = p[1];
    double gamma = p[2];
    double chi2 = 0;
    double last_menc = 0;
    double bin_mass;// = nfw_menc(hp.bins[hp.newf.nb_rvir-1.],rs,c0)/(hp.newf.nb_rvir-1.);

    // note start from i = 0 because first bin will be at r = 0
    for (int i = 1; i < hp.bins.size(); i++) {

        double r = hp.bins[i], menc;

        // only want to fit up to the virial radius
        if (r > hp.radius) break;

        // we have different forms for n = 1 and n = 2 exact solutions since general function is indeterminate there.
        // these are included for completeness but in practice are not really necessary (since could just evaluate
        // general function at gamma = 2.000001 etc
        if (gamma == 1.0) {
            menc = c0 * rs*rs * (r - rs * log (1.+r/rs));
        }
        else if (gamma == 2.0) {
            menc = c0 * rs*rs*rs * (log (1.+r/rs) - r/(r+rs));
        }
        else {
            menc = c0 * rs*rs * (rs - pow(1.+r/rs,1.-gamma) * (rs + (gamma-1.)*r));
            menc /= (gamma-1.)*(gamma-2.);
        }

        double dx = hp.weights[i]*(menc-last_menc-hp.ppbin)/hp.ppbin;
        chi2 += dx*dx;
        last_menc = menc;
    }

    return chi2;
}

//double chi2_nfw2 (halo_profile & hp, std::vector<double> & p) {
//
//    if (p.size() != 2) {
//        std::cout << "chi2_nfw passed " << p.size() << " parameters instead of 2." << std::endl;
//    }
//
//    double rs = p[0];
//    double c0 = p[1];
//    double chi2 = 0;
//    double menc = 0;
//
//    // note start from i = 0 because first bin will be at r = 0
//    for (int i = 1; i < hp.bins.size(); i++) {
//
//        // only want to fit up to the virial radius
//        if (hp.bins[i] > hp.radius) break;
//    
//        menc += hp.ppbin;
//        chi2 += pow (hp.weights[i]*(menc-nfw_menc(hp.bins[i],rs,c0)),2.);
//    }
//
//    return chi2;
//}

double chi2_nfw (halo_profile & hp, double rs, double c0) {

    double chi2 = 0;
    double menc = 0;

    // note start from i = 0 because first bin will be at r = 0
    for (int i = 1; i < hp.bins.size(); i++) {

        // only want to fit up to the virial radius
        if (hp.bins[i] > hp.radius) break;
    
        menc += hp.ppbin;
        chi2 += pow (hp.weights[i]*(menc-nfw_menc(hp.bins[i],rs,c0)),2.);
    }

    return chi2;
}

// a += b
int vadd (std::vector<double> & a, std::vector<double> & b) {

    for (int i = 0; i < a.size(); i++) a[i] += b[i];
    return 0;
}

// a[] *= b[]
int vmult (std::vector<double> & a, std::vector<double> & b) {
    for (int i = 0; i < a.size(); i++) a[i] *= b[i];
    return 0;
}

// a[] *= b
int vmult (std::vector<double> & a, double b) {
    for (int i = 0; i < a.size(); i++) a[i] *= b;
    return 0;
}

// a[] += b[] * c
int vaddmult (std::vector<double> & a, std::vector<double> & b, double c) {
    for (int i = 0; i < a.size(); i++) a[i] += b[i] * c;
    return 0;
}

// assumes m is nxn square matrix
int invert (std::vector<std::vector<double>> & m) {

    if (m.size() == 0) {
        double size;
        std::cout << "Scanning matrix from input. Enter matrix size (nxn):";
        std::cin >> size;
        std::cout << "Now enter each element of the matrix:";
        std::vector<double> row (size);
        for (int i = 0; i < row.size(); i++) {
            for (int j = 0; j < row.size(); j++) {
                std::cin >> row[j];
            }
            m.push_back(row);
        }
        std::cout << "done." << std::endl;
    }

    //std::cout << "Input matrix m is: " << std::endl;
    //for (int i = 0; i < m.size(); i++) {
    //    for (int j = 0; j < m.size(); j++) {

    //        std::cout << m[i][j] << "\t\t";
    //    }
    //    std::cout << std::endl;
    //}

    // use Gauss Jordan method to find inverse of n x n matrix
    std::vector<std::vector<double>> inv (m);

    // initialize to identity matrix
    for (int i = 0; i < m.size(); i++) {
        for (int j = 0; j < m.size(); j++) inv[i][j] = i==j ? 1. : 0.;
    }

    // perform row reductions.
    // n-1 different reductions to apply top down
    for (int i = 0; i < m.size()-1; i++) {

        for (int j = i+1; j < m.size(); j++) {

            //vmult (m[i][i]m[j]
            //m[j] += -m[i]*(m[j][i]/m[i][i])
            if (m[i][i] == 0) {
                std::cerr << "Matrix not invertible." << std::endl;
                return -1;
            }
            else if (m[j][i] == 0) continue;
            double factor = -m[j][i]/m[i][i];
            vaddmult (inv[j], inv[i], factor);
            vaddmult (m[j], m[i], factor);
        }
    }

    // now n-1 more reductions to apply bottom up
    for (int i = m.size()-1; i >= 0; i--) {

        for (int j = i-1; j >= 0; j--) {

            if (m[i][i] == 0) {
                std::cerr << "Matrix not invertible." << std::endl;
                return -1;
            }
            else if (m[j][i] == 0) continue;

            double factor = -m[j][i]/m[i][i];
            vaddmult (inv[j], inv[i], factor);
            vaddmult (m[j], m[i], factor);
        }
    }

    // finally, normalize diagonal elements to unity
    for (int i = 0; i < m.size(); i++) {
        if (m[i][i] == 0. || m[i][i] == 1.) continue;
        vmult (inv[i],1./m[i][i]);
        vmult (m[i],1./m[i][i]);
    }

    //std::cout << "Original matrix m is: " << std::endl;
    //for (int i = 0; i < m.size(); i++) {
    //    for (int j = 0; j < m.size(); j++) {

    //        std::cout << m[i][j] << "\t\t";
    //    }
    //    std::cout << std::endl;
    //}

    //std::cout << "Inverse is: " << std::endl;
    //for (int i = 0; i < m.size(); i++) {
    //    for (int j = 0; j < m.size(); j++) {

    //        std::cout << inv[i][j] << "\t\t";
    //    }
    //    std::cout << std::endl;
    //}

    // now change m to match inv
    for (int i = 0; i < m.size(); i++) {
        for (int j = 0; j < m.size(); j++) {
            m[i][j] = inv[i][j];
        }
    }

    return 0;
}

// ADADELTA implementation
int doFitADA (halo_profile * &h, std::vector<double> & p, double (*chi2f) (halo_profile &, std::vector<double> &), std::vector<bound> & bounds) {

    bool debug = !(std::cout.fail());
    bool write_of = true;

    int nump = p.size();

    if (nump == 0) {
        if (debug) std::cout << "Need to provide at least one parameter initialization" << std::endl;
        return -1;
    }

    double chi2 = chi2f (*h, p);
    double last_chi2 = 1e30;
    double global_min_chi2 = last_chi2;

    // revert to this at end of fitting if different from final value
    std::vector<double> global_min_p (nump);

    std::vector<double> dp_step (nump), dp (nump), move (nump), new_p(nump);

    std::vector<double> J (nump);

    std::vector<double> G (nump);
    std::vector<double> E (nump);
    std::vector<double> update (nump);

    double lambda = 0.9, epsilon = 1.e-8;

    for (int i = 0; i < nump; i++) {
        dp_step[i] = 10000.;
        dp[i] = 1e-8; //p[i]/dp_step[i];
        move[i] = 0.;
        new_p[i] = 0.;
        G[i] = 0.;
        E[i] = 0.;
    }

    std::ofstream of;

    if (write_of) of.open("error_trajectory_ada.dat",std::ofstream::out);

    int it = 0;
    int total_its = -1;
    bool converged = false;
    std::queue<double> change;
    int change_queue_size = 50;
    double sum_change = 0;
    
    if (debug) std::cout << "\tBeginning fitting procedure..." << std::endl;
    
    while (true) {
    
        total_its++;
    
        for (int i = 0; i < nump; i++) {
            //dp[i] = p[i]/dp_step[i];
            update[i] = p[i];
        }

        if (it && (change.size() == change_queue_size) && (fabs(sum_change/change.size()-chi2) < 1.e-6*chi2)) {

            //std::cout << "Converged with parameters rs: " << rs << ", c0: " << c0 << std::endl;
            if (debug) {
                std::cout << "\tConverged after " << it << " iterations with p[0]: " << p[0];
                for (int i = 1; i < nump; i++) std::cout << ", p[" << i << "]: " << p[i];
                std::cout << std::endl;
            }
            converged = true;
            break;
        }

        //if (it > 10000) {
        //    break;
        //}
   
        if (debug) {
            if (change.size()) std::cout << "\tAvg change: " << sum_change/change.size() << std::endl;
            std::cout << "\titer " << it << ", chi2: " << chi2 << ", %E:" << fabs(last_chi2-chi2)/chi2;
            for (int i = 0; i < nump; i++) std::cout << ", p[" << i << "]: " << p[i] << ", dp[" << i << "]: " << dp[i];
            std::cout << std::endl;
        }
    
        // compute Jacobian (gradient in this case, since we have a scalar field)
        for (int i = 0; i < nump; i++) {

            if (bounds.size()) {
                if      (p[i]-dp[i] < bounds[i].lower) dp[i] = 0.9*(p[i]-bounds[i].lower);
                else if (p[i]+dp[i] > bounds[i].upper) dp[i] = 0.9*(bounds[i].upper-p[i]);
            }

            p[i] += dp[i];
            J[i] =  chi2f(*h, p);
            p[i] -= 2*dp[i];
            J[i] -= chi2f(*h, p);
            p[i] += dp[i];
            J[i] /= 2*dp[i];
        }

        double min_chi2 = chi2;

        for (int i = 0; i < nump; i++) {

            // for each parameter, update G. Previously stored value was G_t-1 from last iteration.
            //double logJ = log10(J[i]);

            G[i] = lambda * G[i] + (1-lambda) * (J[i]*J[i]);

            //double move_scale = fabs(log10(J[i]));

            //move[i] = - pow((E[i]+epsilon)/(G[i]+epsilon),0.5) * J[i];
            move[i] = - J[i];


            if (bounds.size()) {
                if      (p[i]+move[i] < bounds[i].lower) move[i] = 0.9*(bounds[i].lower-p[i]);
                else if (p[i]+move[i] > bounds[i].upper) move[i] = 0.9*(bounds[i].upper-p[i]);
            }

            //if      (move[i] > 3.*p[i])     move[i] = 3.*p[i];
            //else if (move[i] < -0.75*p[i])  move[i] = -0.75*p[i];

            new_p[i] = p[i] + move[i];

            if (debug) std::cout << "\t\tp["<<i<<"]: move: " << move[i] << ", J: " << J[i] << ", G: " << G[i] << ", E: " << E[i] << std::endl;

        }

        double new_chi2 = chi2f(*h,new_p);

        if (debug) std::cout << "\t\tnew_chi2: " << new_chi2 << std::endl;

        // if we are adopting a worse solution, then check global mins
        if (new_chi2 > chi2) {
            if (chi2 < global_min_chi2) {
                for (int i = 0; i < nump; i++) global_min_p[i] = p[i];
                global_min_chi2 = chi2;
            }
        }

        for (int i = 0; i < nump; i++) {
            //dp[i] = new_p[i]/dp_step[i];
            p[i] = new_p[i];
        }

        if (debug) {
            std::cout << "\t\tDecided to keep values: ";
            for (int i = 0; i < nump; i++) {
                if (i) std::cout << ", "; 
                std::cout << "p" << i << ": " << p[i];
            }
            std::cout << ", with new chi2: " << new_chi2 << ", an improvement of " << (chi2-new_chi2)/new_chi2*100. << "%" << std::endl;
        }

        if (write_of) {
            of << p[0];
            for (int i = 1; i < nump; i++) of << " " << p[i];
            of << " " << total_its << std::endl;
        }
    
        // want to only update last_chi2 if we found new values, otherwise we would be starting new iteration
        // with no change in chi2 -- a problem.
        last_chi2 = chi2;
        chi2 = new_chi2;

        for (int i = 0; i < nump; i++) {

            update[i] = p[i] - update[i];

            // this is now setting E_t (for this iteration), so next iteration when used it will be E_t-1
            E[i] = lambda * E[i] + (1.-lambda) * (update[i] * update[i]);
        }

        change.push(last_chi2);
        sum_change += change.back();

        if (change.size() > change_queue_size) {
            sum_change -= change.front();
            change.pop();
        }

        it++;
    }

    if (write_of) of.close();

    if (global_min_chi2 < chi2) {
        if (debug) {
            std::cout << "\tAdopting previous minimum chi2 solution, with p[0]: " << global_min_p[0];
            for (int i = 1; i < nump; i++) std::cout << ", p[" << i << "]: " << global_min_p[i];
            std::cout << std::endl;
        }
        for (int i = 0; i < nump; i++) p[i] = global_min_p[i];
    }

    if (debug) std::cout << "\tIterated in total " << total_its << " times (including overshooting)." << std::endl;

    return total_its;
}

// Levenburg - Marquardt implementation
int doFitLM (halo_profile * &h, std::vector<double> & p, double (*chi2f) (halo_profile &, std::vector<double> &), std::vector<bound> & bounds) {

    bool debug = !(std::cout.fail());
    bool write_of = true;

    int nump = p.size();

    if (nump == 0) {
        if (debug) std::cout << "Need to provide at least one parameter initialization" << std::endl;
        return -1;
    }

    double chi2 = chi2f (*h, p);
    double last_chi2 = 1e30;
    double global_min_chi2 = last_chi2;

    // revert to this at end of fitting if different from final value
    std::vector<double> global_min_p (nump);

    std::vector<double> dp_step (nump), dp (nump), move (nump), new_p(nump);

    std::vector<double> J (nump);
    std::vector<std::vector<double>> H;

    for (int i = 0; i < nump; i++) {
        std::vector<double> hi (nump);
        H.push_back(hi);
    }

    double lambda = 0.001;

    for (int i = 0; i < nump; i++) {
        dp_step[i] = 10000.;
        dp[i] = p[i]/dp_step[i];
        move[i] = 0.;
        new_p[i] = 0.;
    }

    std::ofstream of;

    if (write_of) of.open("error_trajectory_lm.dat",std::ofstream::out);

    int it = 0;
    int total_its = -1;
    bool converged = false;
    std::queue<double> change;
    double sum_change = 0;
    int change_queue_size = 10;
    
    if (debug) std::cout << "\tBeginning fitting procedure..." << std::endl;
    
    while (true) {
    
        total_its++;
    
        for (int i = 0; i < nump; i++) {
            dp[i] = p[i]/dp_step[i];
        }

        if (it && (change.size() == change_queue_size) && (fabs(sum_change/change.size()-chi2) < 0.00005*chi2)) {

            //std::cout << "Converged with parameters rs: " << rs << ", c0: " << c0 << std::endl;
            if (debug) {
                std::cout << "\tConverged after " << it << " iterations with p[0]: " << p[0];
                for (int i = 1; i < nump; i++) std::cout << ", p[" << i << "]: " << p[i];
                std::cout << std::endl;
            }
            converged = true;
            break;
        }

        //if (it > 10000) {
        //    break;
        //}
   
        if (debug) {
            if (change.size()) std::cout << "\tAvg change: " << sum_change/change.size() << std::endl;
            std::cout << "\titer " << it << ", chi2: " << chi2 << ", %E:" << fabs(last_chi2-chi2)/chi2;
            for (int i = 0; i < nump; i++) std::cout << ", p[" << i << "]: " << p[i] << ", dp[" << i << "]: " << dp[i];
            std::cout << std::endl;
        }
    
        // compute Jacobian (gradient in this case, since we have a scalar field)
        // and Hessian
        for (int i = 0; i < nump; i++) {

            if (bounds.size()) {
                if      (p[i]-dp[i] < bounds[i].lower) dp[i] = 0.9*(p[i]-bounds[i].lower);
                else if (p[i]+dp[i] > bounds[i].upper) dp[i] = 0.9*(bounds[i].upper-p[i]);
            }


            //if (debug) std::cout << "\t\tComputation of J[" << i << "] = ";
            p[i] += dp[i];
            J[i] =  chi2f(*h, p);
            //if (debug) std::cout << J[i] << " - ";
            p[i] -= 2*dp[i];
            //if (debug) std::cout << chi2f(*h,p) << " / ";
            J[i] -= chi2f(*h, p);
            p[i] += dp[i];
            //if (debug) std::cout << 2*dp[i] << std::endl;
            J[i] /= 2*dp[i];

            for (int j = 0; j < nump; j++) {

                //p[i] += dp[i];
                //p[j] += dp[j];
                //H[i][j] = chi2f(*h,p);  // +f(i+di,j+dj)
                //p[i] -= 2*dp[i];
                //H[i][j] -= chi2f(*h,p); // -f(i-di,j+dj)
                //p[j] -= 2*dp[j];
                //H[i][j] += chi2f(*h,p); // +f(i-di,j-dj)
                //p[i] += 2*dp[i];
                //H[i][j] -= chi2f(*h,p); // -f(i+di,j-dj)
                //p[i] -= dp[i];
                //p[j] += dp[j];
                //H[i][j] /= 4*dp[i]*dp[j];

                H[i][j] = J[i]*J[j];

                if (i == j) H[i][j] *= (1.+lambda);
            }
        }

        if (invert(H) < 0) {
            return -1;
        }

        double min_chi2 = chi2;

        for (int i = 0; i < nump; i++) {

            double mult;

            for (int j = 0; j < nump; j++) mult += H[i][j]*J[j];

            // Newton-Raphson update
            move[i] = -mult;

            if (bounds.size()) {
                if      (p[i]+move[i] < bounds[i].lower) move[i] = 0.9*(bounds[i].lower-p[i]);
                else if (p[i]+move[i] > bounds[i].upper) move[i] = 0.9*(bounds[i].upper-p[i]);
            }

            //if      (move[i] > 3.*p[i])     move[i] = 3.*p[i];
            //else if (move[i] < -0.75*p[i])  move[i] = -0.75*p[i];

            new_p[i] = p[i] + move[i];

            if (debug) {
                std::cout << "\t\tp["<<i<<"]: move: " << move[i] << ", J: " << J[i];
                for (int j = 0; j < nump; j++) std::cout << ", H["<<i<<"]["<<j<<"]: " << H[i][j];
                std::cout << ", dp: " << dp[i] << std::endl;
            }
        }

        double new_chi2 = chi2f(*h,new_p);

        if (debug) std::cout << "\t\tnew_chi2: " << new_chi2 << std::endl;

        // reject values
        if (new_chi2 > chi2) {
            lambda *= 10.;

            if (debug) {
                std::cout << "\t\tDecided to reject values. New lambda: " << lambda << std::endl;;
            }

            continue;
        }

        else {
            lambda /= 10.;
        }

        min_chi2 = new_chi2;

        for (int i = 0; i < nump; i++) {
            dp[i] = new_p[i]/dp_step[i];
            p[i] = new_p[i];
        }
        
        if (debug) {
            std::cout << "\t\tDecided to keep values: ";
            for (int i = 0; i < nump; i++) {
                if (i) std::cout << ", "; 
                std::cout << "p" << i << ": " << p[i];
            }
            std::cout << ", with new lambda: " << lambda << ", and chi2: " << min_chi2 << ", an improvement of " << (chi2-min_chi2)/min_chi2*100. << "%" << std::endl;
        }

        if (write_of) {
            of << p[0];
            for (int i = 1; i < nump; i++) of << " " << p[i];
            of << " " << total_its << std::endl;
        }
    
        // want to only update last_chi2 if we found new values, otherwise we would be starting new iteration
        // with no change in chi2 -- a problem.
        last_chi2 = chi2;
        chi2 = min_chi2;

        change.push(last_chi2);
        sum_change += change.back();

        if (change.size() > change_queue_size) {
            sum_change -= change.front();
            change.pop();
        }

        it++;
    }

    if (write_of) of.close();

    if (global_min_chi2 < chi2) {
        if (debug) {
            std::cout << "\tAdopting previous minimum chi2 solution, with p[0]: " << global_min_p[0];
            for (int i = 1; i < nump; i++) std::cout << ", p[" << i << "]: " << global_min_p[i];
            std::cout << std::endl;
        }
        for (int i = 0; i < nump; i++) p[i] = global_min_p[i];
    }

    if (debug) std::cout << "\tIterated in total " << total_its << " times (including overshooting)." << std::endl;

    return total_its;
}

// Vanilla Batch Gradient Descent
int doFitBGD (halo_profile * &h, std::vector<double> & p, double (*chi2f) (halo_profile &, std::vector<double> &), std::vector<bound> & bounds) {

    bool debug = !(std::cout.fail());
    bool write_of = false;//true;

    int nump = p.size();

    if (nump == 0) {
        if (debug) std::cout << "Need to provide at least one parameter initialization" << std::endl;
        return -1;
    }

    double chi2 = chi2f (*h, p);
    double last_chi2 = 1e30;
    double global_min_chi2 = last_chi2;

    // revert to this at end of fitting if different from final value
    std::vector<double> global_min_p (nump);

    std::vector<double> dp_step (nump), dp (nump), lr (nump), move (nump), new_p(nump), new_chi2_p(nump);

    std::vector<double> J (nump);

    double lr_up = 2., lr_down = 5.;

    for (int i = 0; i < nump; i++) {
        dp_step[i] = 10000.;
        dp[i] = p[i]/dp_step[i];
        lr[i] = 0.01;
        move[i] = 0.;
        new_p[i] = 0.;
    }

    std::ofstream of;

    if (write_of) of.open("error_trajectory_nr.dat",std::ofstream::out);

    int it = 0;
    int total_its = -1;
    bool converged = false;
    std::queue<double> change;
    double sum_change = 0;
    int change_queue_size = 5;
    
    if (debug) std::cout << "\tBeginning fitting procedure..." << std::endl;
    
    while (true) {
    
        total_its++;
    
        for (int i = 0; i < nump; i++) {
            dp[i] = p[i]/dp_step[i];
        }

        if (lr[0] < 1e-20 || (it && (change.size() == change_queue_size) && (fabs(sum_change/change.size()-chi2) < 1e-8*chi2))) {

            //std::cout << "Converged with parameters rs: " << rs << ", c0: " << c0 << std::endl;
            if (debug) {
                std::cout << "\tConverged after " << it << " iterations with p[0]: " << p[0];
                for (int i = 1; i < nump; i++) std::cout << ", p[" << i << "]: " << p[i];
                std::cout << std::endl;
            }
            converged = true;
            break;
        }

        if (it > 50000) {
            break;
        }
   
        if (debug) {
            if (change.size()) std::cout << "\tAvg change: " << sum_change/change.size() << std::endl;
            std::cout << "\titer " << it << ", chi2: " << chi2 << ", %E:" << fabs(last_chi2-chi2)/chi2;
            for (int i = 0; i < nump; i++) std::cout << ", p[" << i << "]: " << p[i] << ", dp[" << i << "]: " << dp[i];
            std::cout << std::endl;
        }
    
        // compute Jacobian (gradient in this case, since we have a scalar field)
        for (int i = 0; i < nump; i++) {

            if (bounds.size()) {
                if      (p[i]-dp[i] < bounds[i].lower) dp[i] = 0.9*(p[i]-bounds[i].lower);
                else if (p[i]+dp[i] > bounds[i].upper) dp[i] = 0.9*(bounds[i].upper-p[i]);
            }

            p[i] += dp[i];
            J[i] =  chi2f(*h, p);
            p[i] -= 2*dp[i];
            J[i] -= chi2f(*h, p);
            p[i] += dp[i];
            J[i] /= 2*dp[i];
        }

        double min_chi2 = chi2;

        for (int i = 0; i < nump; i++) {

            // Basic gradient descent update
            move[i] = -lr[i] * J[i];

            if (bounds.size()) {
                if      (p[i]+move[i] < bounds[i].lower) move[i] = 0.9*(bounds[i].lower-p[i]);
                else if (p[i]+move[i] > bounds[i].upper) move[i] = 0.9*(bounds[i].upper-p[i]);
            }

            //if      (move[i] > 3.*p[i])     move[i] = 3.*p[i];
            //else if (move[i] < -0.75*p[i])  move[i] = -0.75*p[i];

            new_p[i] = p[i] + move[i];

            p[i] += move[i];
            new_chi2_p[i] = chi2f(*h, p);
            p[i] -= move[i];

            if (debug) std::cout << "\t\tp["<<i<<"]: learning rate: " << lr[i] << ", move: " << -lr[i]*J[i] << ", move(bounded): " << move[i] << ", new_chi2: " << new_chi2_p[i] << ", min_chi2? " << ((new_chi2_p[i] < min_chi2) ? "Yes":"No") << std::endl;

            if (new_chi2_p[i] < min_chi2) min_chi2 = new_chi2_p[i];
        }

        double new_chi2 = chi2f(*h,new_p);

        if (debug) std::cout << "\t\tnew_chi2: " << new_chi2 << std::endl;

        if (new_chi2 < min_chi2) min_chi2 = new_chi2;

        for (int i = 0; i < nump; i++) {
            if (new_chi2_p[i] < chi2) lr[i] *= lr_up;
            else lr[i] /= lr_down;
        }

        // let's just keep both values and see what happens
        if (new_chi2 == min_chi2) {
            for (int i = 0; i < nump; i++) {
                dp[i] = new_p[i]/dp_step[i];
                p[i] = new_p[i];
            }

        }

        else {

            bool cont = false;

            for (int i = 0; i < nump; i++) {
                if (new_chi2_p[i] == min_chi2) {
                    dp[i] = new_p[i]/dp_step[i];
                    p[i] = new_p[i];
                    break;
                }

                // if none of our new parameter guesses got a new min_chi2, then reject and try again
                if (i == nump -1) {
                    if (debug) std::cout << "\t\tDecided to reject new values." << std::endl;
                    cont = true;
                }
            }

            if (cont) continue;
        }

        if (debug) {
            std::cout << "\t\tDecided to keep values: ";
            for (int i = 0; i < nump; i++) {
                if (i) std::cout << ", "; 
                std::cout << "p" << i << ": " << p[i];
            }
            std::cout << ", with new chi2: " << min_chi2 << ", an improvement of " << (chi2-min_chi2)/min_chi2*100. << "%" << std::endl;
        }

        if (write_of) {
            of << p[0];
            for (int i = 1; i < nump; i++) of << " " << p[i];
            of << " " << total_its << std::endl;
        }
    
        // want to only update last_chi2 if we found new values, otherwise we would be starting new iteration
        // with no change in chi2 -- a problem.
        last_chi2 = chi2;
        chi2 = min_chi2;

        change.push(last_chi2);
        sum_change += change.back();

        if (change.size() > change_queue_size) {
            sum_change -= change.front();
            change.pop();
        }

        it++;
    }

    if (write_of) of.close();

    if (global_min_chi2 < chi2) {
        if (debug) {
            std::cout << "\tAdopting previous minimum chi2 solution, with p[0]: " << global_min_p[0];
            for (int i = 1; i < nump; i++) std::cout << ", p[" << i << "]: " << global_min_p[i];
            std::cout << std::endl;
        }
        for (int i = 0; i < nump; i++) p[i] = global_min_p[i];
    }

    if (debug) std::cout << "\tIterated in total " << total_its << " times (including overshooting)." << std::endl;

    return total_its;
}

int doFitGD (halo_profile * &h, std::vector<double> & p, double (*chi2f) (halo_profile &, std::vector<double> &), std::vector<bound> & bounds, int method) {

    bool debug = !(std::cout.fail());
    bool write_of = true;

    bool BGD = false, NR = false, ADA = false;

    if (method > 0) {
        switch (method) {
            case 1: BGD = true; if (debug) std::cout << "Using Batch Gradient Descent optimization." << std::endl; break;
            case 2: NR = true; if (debug) std::cout << "Using Newton-Raphson optimization." << std::endl; break;
            case 3: ADA = true; if (debug) std::cout << "Using Adadelta Gradient Descent optimization." << std::endl; break;
        }
    }

    // default
    else {
        ADA = true;
    }

    int nump = p.size();

    if (nump == 0) {
        if (debug) std::cout << "Need to provide at least one parameter initialization" << std::endl;
        return -1;
    }

    double chi2 = chi2f (*h, p);
    double last_chi2 = 1e30;
    double global_min_chi2 = last_chi2;

    // revert to this at end of fitting if different from final value
    std::vector<double> global_min_p (nump);

    std::vector<double> dp_step (nump), dp (nump), lr (nump), move (nump), new_p(nump), new_chi2_p(nump);

    std::vector<double> J (nump);
    std::vector<std::vector<double>> H;

    for (int i = 0; i < nump; i++) {
        std::vector<double> hi (nump);
        H.push_back(hi);
    }

    std::vector<double> G (nump);
    std::vector<double> E (nump);
    std::vector<double> update (nump);

    double lambda = 0.5, epsilon = 1.e-8;

    for (int i = 0; i < nump; i++) {
        dp_step[i] = 10000.;
        dp[i] = p[i]/dp_step[i];
        lr[i] = 0.01;
        move[i] = 0.;
        new_p[i] = 0.;
        G[i] = 0.;
        E[i] = 0.;
    }

    std::ofstream of;

    if (write_of) of.open("error_trajectory_nr.dat",std::ofstream::out);

    int it = 0;
    int total_its = -1;
    bool converged = false;
    std::queue<double> change;
    double sum_change = 0;
    int change_queue_size = 10;
    
    if (debug) std::cout << "\tBeginning fitting procedure..." << std::endl;
    
    while (true) {
    
        total_its++;
    
        for (int i = 0; i < nump; i++) {
            dp[i] = p[i]/dp_step[i];
            update[i] = p[i];
        }

        if (it && (change.size() == change_queue_size) && (fabs(sum_change/change.size()-chi2) < 1e-8*chi2)) {

            //std::cout << "Converged with parameters rs: " << rs << ", c0: " << c0 << std::endl;
            if (debug) {
                std::cout << "\tConverged after " << it << " iterations with p[0]: " << p[0];
                for (int i = 1; i < nump; i++) std::cout << ", p[" << i << "]: " << p[i];
                std::cout << std::endl;
            }
            converged = true;
            break;
        }

        if (it > 10000) {
            break;
        }
   
        if (debug) {
            if (change.size()) std::cout << "\tAvg change: " << sum_change/change.size() << std::endl;
            std::cout << "\titer " << it << ", chi2: " << chi2 << ", %E:" << fabs(last_chi2-chi2)/chi2;
            for (int i = 0; i < nump; i++) std::cout << ", p[" << i << "]: " << p[i] << ", dp[" << i << "]: " << dp[i];
            std::cout << std::endl;
        }
    
        // compute Jacobian (gradient in this case, since we have a scalar field)
        // and Hessian
        for (int i = 0; i < nump; i++) {

            if (bounds.size()) {
                if      (p[i]-dp[i] < bounds[i].lower) dp[i] = 0.9*(p[i]-bounds[i].lower);
                else if (p[i]+dp[i] > bounds[i].upper) dp[i] = 0.9*(bounds[i].upper-p[i]);
            }

            p[i] += dp[i];
            J[i] =  chi2f(*h, p);
            p[i] -= 2*dp[i];
            J[i] -= chi2f(*h, p);
            p[i] += dp[i];
            J[i] /= 2*dp[i];

            if (NR) {
                //std::vector<double> hi (nump);
                //H.push_back(hi);

                for (int j = 0; j < nump; j++) {

                    p[i] += dp[i];
                    p[j] += dp[j];
                    H[i][j] = chi2f(*h,p);  // +f(i+di,j+dj)
                    p[i] -= 2*dp[i];
                    H[i][j] -= chi2f(*h,p); // -f(i-di,j+dj)
                    p[j] -= 2*dp[j];
                    H[i][j] += chi2f(*h,p); // +f(i-di,j-dj)
                    p[i] += 2*dp[i];
                    H[i][j] -= chi2f(*h,p); // -f(i+di,j-dj)
                    p[i] -= dp[i];
                    p[j] += dp[j];
                    H[i][j] /= 4*dp[i]*dp[j];
                }
            }
        }

        if (NR) {
            if (invert(H) < 0) {
                return -1;
            }
        }

        double min_chi2 = chi2;

        for (int i = 0; i < nump; i++) {

            if (NR) {
                double mult;

                for (int j = 0; j < nump; j++) mult += H[i][j]*J[j];

                // Newton-Raphson update
                move[i] = -mult;//;-lr[i] * mult;
            }

            if (BGD) {
                // Basic gradient descent update
                move[i] = -lr[i] * J[i];
            }

            if (ADA) {
                // Adadelta update
                G[i] = lambda * G[i] + (1-lambda) * (J[i]*J[i]);
                move[i] = - sqrt(E[i]+epsilon)/sqrt(G[i]+epsilon) * J[i];
            }

            if (bounds.size()) {
                if      (p[i]+move[i] < bounds[i].lower) move[i] = 0.9*(bounds[i].lower-p[i]);
                else if (p[i]+move[i] > bounds[i].upper) move[i] = 0.9*(bounds[i].upper-p[i]);
            }

            //if      (move[i] > 3.*p[i])     move[i] = 3.*p[i];
            //else if (move[i] < -0.75*p[i])  move[i] = -0.75*p[i];

            new_p[i] = p[i] + move[i];

            if (BGD) {
                p[i] += move[i];
                new_chi2_p[i] = chi2f(*h, p);
                p[i] -= move[i];
            }

            if (debug) std::cout << "\t\tp["<<i<<"]: learning rate: " << lr[i] << ", move: " << move[i] << ", new_chi2: " << new_chi2_p[i] << ", min_chi2? " << ((new_chi2_p[i] < min_chi2) ? "Yes":"No") << std::endl;

            if (BGD) {
                if (new_chi2_p[i] < min_chi2) min_chi2 = new_chi2_p[i];
            }
        }

        double new_chi2 = chi2f(*h,new_p);

        if (debug) std::cout << "\t\tnew_chi2: " << new_chi2 << std::endl;


        if (BGD) {

            if (new_chi2 < min_chi2) min_chi2 = new_chi2;

            for (int i = 0; i < nump; i++) {
                if (new_chi2_p[i] < chi2) lr[i] *= 2;
                else lr[i] /= 2;
            }

            // let's just keep both values and see what happens
            if (new_chi2 == min_chi2) {
                for (int i = 0; i < nump; i++) {
                    dp[i] = new_p[i]/dp_step[i];
                    p[i] = new_p[i];
                }

            }

            else {


                
                bool cont = false;

                for (int i = 0; i < nump; i++) {
                    if (new_chi2_p[i] == min_chi2) {
                        dp[i] = new_p[i]/dp_step[i];
                        p[i] = new_p[i];
                        break;
                    }

                    // if none of our new parameter guesses got a new min_chi2, then reject and try again
                    if (i == nump -1) {
                        if (debug) std::cout << "\t\tDecided to reject new values." << std::endl;
                        cont = true;
                    }
                }

                if (cont) continue;
            }
        }

        if (NR) {

            min_chi2 = new_chi2;

            for (int i = 0; i < nump; i++) {
                //if (new_chi2 < chi2) lr[i] *= 2;
                //else lr[i] /= 2.;
                dp[i] = new_p[i]/dp_step[i];
                p[i] = new_p[i];
            }
        }

        if (ADA) {

            min_chi2 = new_chi2;

            // if we are adopting a worse solution, then check global mins
            if (min_chi2 > chi2) {
                if (chi2 < global_min_chi2) {
                    for (int i = 0; i < nump; i++) global_min_p[i] = p[i];
                    global_min_chi2 = chi2;
                }
            }

            for (int i = 0; i < nump; i++) {
                dp[i] = new_p[i]/dp_step[i];
                p[i] = new_p[i];
            }
        }

        if (debug) {
            std::cout << "\t\tDecided to keep values: ";
            for (int i = 0; i < nump; i++) {
                if (i) std::cout << ", "; 
                std::cout << "p" << i << ": " << p[i];
            }
            std::cout << ", with new chi2: " << min_chi2 << ", an improvement of " << (chi2-min_chi2)/min_chi2*100. << "%" << std::endl;
        }



        if (write_of) {
            of << p[0];
            for (int i = 1; i < nump; i++) of << " " << p[i];
            of << " " << total_its << std::endl;
        }
    
        // want to only update last_chi2 if we found new values, otherwise we would be starting new iteration
        // with no change in chi2 -- a problem.
        last_chi2 = chi2;
        chi2 = min_chi2;

        if (ADA) {
            // ADADELTA
            for (int i = 0; i < nump; i++) {

                update[i] = p[i] - update[i];

                E[i] = lambda * E[i] + (1.-lambda) * (update[i] * update[i]);
            }
        }

        change.push(last_chi2);
        sum_change += change.back();

        if (change.size() > change_queue_size) {
            sum_change -= change.front();
            change.pop();
        }

        it++;
    }

    if (write_of) of.close();

    if (global_min_chi2 < chi2) {
        if (debug) {
            std::cout << "\tAdopting previous minimum chi2 solution, with p[0]: " << global_min_p[0];
            for (int i = 1; i < nump; i++) std::cout << ", p[" << i << "]: " << global_min_p[i];
            std::cout << std::endl;
        }
        for (int i = 0; i < nump; i++) p[i] = global_min_p[i];
    }

    if (debug) std::cout << "\tIterated in total " << total_its << " times (including overshooting)." << std::endl;

    return total_its;
}

int doFitGD (halo_profile * &h, double & rs, double & c0) {

    double chi2 = chi2_nfw (*h, rs, c0);
    double last_chi2 = 1e30;
    double drs_step = 1000.;
    double dc0_step = 1000.;
    double drs = rs/drs_step;
    double dc0 = c0/dc0_step;
    int it = 0;
    int total_its = -1;
    bool converged = false;
    double lr_rs = 0.01, lr_c0 = 0.01;
    std::queue<double> change;
    double sum_change = 0;
    
    std::cout << "\tBeginning fitting procedure..." << std::endl;
    
    while (true) {
    
        total_its++;
    
        drs = rs/drs_step;
        dc0 = c0/dc0_step;

        //if (fabs(last_chi2-chi2) < 0.00005*chi2) {
        if (it && fabs(sum_change/change.size()-chi2) < 0.00005*chi2) {
    
            //std::cout << "Converged with parameters rs: " << rs << ", c0: " << c0 << std::endl;
            std::cout << "\tConverged after " << it << " iterations with rs: " << rs << " and c0: " << c0 << ", cnfw: " << h->radius/rs << std::endl;
            converged = true;
            break;
        }
    
        //if (it >= 100) { 
        //    std::cout << "Did not converge after " << it << " iterations, breaking..." << std::endl;
        //    break;
        //}
   
        std::cout << "\tAvg change: " << sum_change/change.size() << std::endl;
        std::cout << "\titer " << it << ", chi2: " << chi2 << ", %E:" << fabs(last_chi2-chi2)/chi2 << ", rs:" << rs << ", drs:" << drs << ", c0:" << c0 << ", dc0:" << dc0 << std::endl;
    
        //last_chi2 = chi2;
   
        double chi2_left_rs;
        double chi2_right_rs;
        double dx_rs;
        double new_rs = rs;
        //double old_
        
        double chi2_left_c0;
        double chi2_right_c0;
        double dx_c0;
        double new_c0 = c0;
    
        double move;
        double new_chi2;
    
        double new_chi2_rs = last_chi2;
        double new_chi2_c0 = last_chi2;
        double min_chi2 = chi2;
    
        // first, lets compute finite differences for each parameter
        if (c0-dc0 < 0) dc0 = 0.1*c0;
    
        chi2_left_c0 = chi2_nfw(*h, rs, c0-dc0);
        chi2_right_c0 = chi2_nfw(*h, rs, c0+dc0);
        dx_c0 = (chi2_right_c0-chi2_left_c0)/(2.*dc0);
    
        // choose new test values for rs
        if (rs+drs > h->radius) drs = 0.1*(h->radius-rs);
        else if (rs-drs < h->bins[0]) drs = 0.1*(rs-h->bins[0]);
    
        chi2_left_rs = chi2_nfw(*h, rs-drs, c0);
        chi2_right_rs = chi2_nfw(*h, rs+drs, c0);
        dx_rs = (chi2_right_rs-chi2_left_rs)/(2.*drs);
    
        move = 0;
    
        std::cout << "\t\tlr_rs: " << lr_rs << ", lr_c0: " << lr_c0 << std::endl;
    
        move = -lr_c0*dx_c0;
        if (move > 3.*c0) move = 3.*c0;
        else if (move < -0.75*c0) move = -0.75*c0;
        new_c0 = c0+move;
        if (new_c0 < 0) new_c0 = 0.1*c0;
        new_chi2_c0 = chi2_nfw(*h,rs,new_c0);
        //if (chi2_right_c0 < new_chi2_c0) {
        //    new_chi2_c0 = chi2_right_c0;
        //    new_c0 = c0+dc0;
        //}
        //if (chi2_left_c0 < new_chi2_c0) {
        //    new_chi2_c0 = chi2_left_c0;
        //    new_c0 = c0-dc0;
        //}
    
        std::cout << "\t\tc0_move: " << move << ", dc0: " << dc0 << std::endl;
    
        std::cout << "\t\tchi2_left: " << chi2_left_c0 << ", chi2_right: " << chi2_right_c0 << ", new_chi2_c0:" << new_chi2_c0 << std::endl;
    
        if (new_chi2_c0 != chi2) {
            std::cout << "\t\tChose new c0: " << new_c0 << " dc0: " << dc0 << " new_chi2 " << new_chi2_c0 << " " << fabs(chi2-new_chi2_c0)/new_chi2_c0 << std::endl;
        }
        else {
            std::cout << "\t\tc0 remains the same." << std::endl;
        }
    
        if (new_chi2_c0 < min_chi2) min_chi2 = new_chi2_c0; 
    
        move = -lr_rs*dx_rs;///dx2_rs;//*mag_dx/mag_dx2;
        if (move > 3.*rs) move = 3.*rs;
        else if (move < -0.75*rs) move = -0.75*rs;
        new_rs = rs+move;
        if (new_rs > h->radius) new_rs = rs + 0.8*(h->radius-rs);
        else if (new_rs < h->bins[0]) new_rs = rs - 0.8*(rs-h->bins[0]);
        new_chi2_rs = chi2_nfw(*h,new_rs,c0);
        //if (chi2_right_rs < new_chi2_rs) {
        //    new_chi2_rs = chi2_right_rs;
        //    new_rs = rs+drs;
        //}
        //if (chi2_left_rs < new_chi2_rs) {
        //    new_chi2_rs = chi2_left_rs;
        //    new_rs = rs-drs;
        //}
    
        std::cout << "\t\trs_move: " << move << ", drs: " << drs << std::endl;
    
        std::cout << "\t\tchi2_left: " << chi2_left_rs << ", chi2_right: " << chi2_right_rs << ", new_chi2_rs: " << new_chi2_rs << std::endl;
    
        if (new_chi2_rs != chi2) {
            std::cout << "\t\tChose new rs: " << new_rs << " drs: " << drs << " new_chi2 " << new_chi2_rs << " " << fabs(chi2-new_chi2_rs)/new_chi2_rs << std::endl;
        }
        else {
            std::cout << "\t\tRs remains the same." << std::endl;
        }
    
        if (new_chi2_rs < min_chi2) min_chi2 = new_chi2_rs;
    
    
        // we now have new_rs and new_c0, so check together
        new_chi2 = chi2_nfw(*h,new_rs,new_c0);
    
        std::cout << "\t\tnew_chi2: " << new_chi2 << std::endl;
    
        if (new_chi2 < min_chi2) min_chi2 = new_chi2;
    
        std::cout << "\t\tmin_chi2: " << min_chi2 << std::endl;

        if (new_chi2_rs > chi2) lr_rs /= 10.;
        else lr_rs *= 10.;

        if (new_chi2_c0 > chi2) lr_c0 /= 10.;
        else lr_c0 *= 10.;

        // let's just keep both values and see what happens
        if (new_chi2 == min_chi2) {
            drs = new_rs/drs_step;
            dc0 = new_c0/dc0_step;
            rs = new_rs;
            c0 = new_c0;
        }

        else if (new_chi2_rs == min_chi2) {
            drs = new_rs/drs_step;
            rs = new_rs;
        }

        else if (new_chi2_c0 == min_chi2) {
            dc0 = new_c0/dc0_step;
            c0 = new_c0;
        }

        else {
            std::cout << "\t\tDecided to reject new values." << std::endl;
            continue;
        }

        // no improved parameter set found -- decrease step size and try again
        //if (chi2 == min_chi2) {
        //    drs_step *= 2.;
        //    dc0_step *= 2.;
    
        //    std::cout << "\t\tDecided to reject new values." << std::endl;
        //    continue;
        //}
    
        //else if (new_chi2_rs == min_chi2) {
        //    //drs = fabs (rs-new_rs)/1.;
        //    //drs = fabs (0.01*new_rs);
        //    drs = new_rs/drs_step;
        //    rs = new_rs;
        //    std::cout << "\t\tKeeping only rs" << std::endl;
        //}
        //else if (new_chi2_c0 == min_chi2) {
        //    //dc0 = fabs (c0-new_c0)/1.;
        //    //dc0 = fabs (0.01*new_c0);
        //    dc0 = new_c0/dc0_step;
        //    c0 = new_c0;
        //    std::cout << "\t\tKeeping only c0" << std::endl;
        //}
        //else if (new_chi2 == min_chi2) {
        //    //drs = fabs (rs-new_rs)/1.;
        //    //dc0 = fabs (c0-new_c0)/1.;
        //    //drs = fabs (new_rs*0.01);
        //    //dc0 = fabs (new_c0*0.01);
        //    drs = new_rs/drs_step;
        //    dc0 = new_c0/dc0_step;
        //    rs = new_rs;
        //    c0 = new_c0;
        //    std::cout << "\t\tKeeping both new values" << std::endl;
        //}
    
        //std::cout << "\t\tDecided to keep values: c0: " << c0 << ", rs: " << rs << ", with new chi2: " << new_chi2 << ", an improvement of " << (chi2-new_chi2)/new_chi2*100. << "%" << std::endl;
        std::cout << "\t\tDecided to keep values: c0: " << c0 << ", rs: " << rs << ", with new chi2: " << min_chi2 << ", an improvement of " << (chi2-min_chi2)/min_chi2*100. << "%" << std::endl;
    
        // want to only update last_chi2 if we found new values, otherwise we would be starting new iteration
        // with no change in chi2 -- a problem.
        last_chi2 = chi2;
        chi2 = min_chi2;

        change.push(last_chi2);
        sum_change += change.back();

        if (change.size() > 5) {
            sum_change -= change.front();
            change.pop();
        }

        it++;
    }

    std::cout << "Iterated in total " << total_its << " times (including overshooting)." << std::endl;

    return 1;
}

int doFitProfiles (std::vector<halo_profile> & halos) {

    std::cout << "Fitting halo profiles... " << std::endl;

    //while (true) {
    //    std::vector<std::vector<double>> m;
    //    invert (m);
    //    std::cout << "Invert another matrix?";
    //    std::string cin_str = "y";
    //    std::cin >> cin_str;
    //    if (cin_str == "n") break;
    //}

    // sort halos by bin number
    std::vector<halo_profile *> hpp;
    
    for (int i = 0; i < halos.size(); i++) {

        hpp.push_back(&halos[i]);
    }

    std::sort(hpp.begin(),hpp.end(),[](halo_profile * left, halo_profile * right){return left->n_rbins > right->n_rbins;});

    // first, we may want to go through and find how many halos we even can fit
    std::vector<int> bin_frequency (hpp[0]->n_rbins+1,0);
    std::vector<double> bin_mass_max (hpp[0]->n_rbins+1,0);

    //for (int i = 0; i < hpp.size(); i++) {

    //    bin_frequency[hpp[i]->n_rbins]++;
    //}

    //for (int i = 0; i < bin_frequency.size(); i++) {

    //    std::cout << "N_rbins: " << i << ", count: " << bin_frequency[i] << std::endl;
    //}

    // now check how many bins halos have within the virial radius
    for (int i = 0; i < hpp.size(); i++) {

        for (int j = 0; j < hpp[i]->bins.size(); j++) {

            if (hpp[i]->bins[j] >= hpp[i]->radius) {

                bin_frequency[j]++;
                if (hpp[i]->mass > bin_mass_max[j]) bin_mass_max[j] = hpp[i]->mass;
                break;
            }
        }
    }

    //for (int i = 0; i < bin_frequency.size(); i++) {

    //    std::cout << "N_rbins: " << i << ", count: " << bin_frequency[i] << ", max mass: " << log10(bin_mass_max[i]) << std::endl;
    //}

    //std::sort(hpp.begin(),hpp.end(),[](halo_profile * left, halo_profile * right){return left->n_rbins < right->n_rbins;});
    std::sort(hpp.begin(),hpp.end(),[](halo_profile * left, halo_profile * right){return left->mass > right->mass;});

    std::vector<halo_profile *> fitted_profiles;

    int skipped = 0;
    double df = 0;
    double cnfw = 0;
    int cnfw_count = 0;
    int of_num = 0;
    int break_after = -1;//1000;//10000;
    clock_t total_t = 0;

    bool debug = false;//true; 

    if (!debug) off_stream(std::cout);

    int dot_freq;
    if (break_after > -1) dot_freq = (int)ceil(break_after/1000.);
    else dot_freq = (int)ceil(hpp.size()/1000.);
    std::cout << hpp.size() << ", " << dot_freq << std::endl;
    int counter = 0;

    for (int i = 0; i < hpp.size(); i++) {

        halo_profile * h = hpp[i];

        //if (h->mass < 10e13) continue;

        // check if halo has enough particles
        int np_rvir = 0;
        int rvir_bin = -1;
        int n_downweights = 0;

        h->weights = std::vector<double> (h->n_rbins,1.0);

        for (int j = 0; j < h->n_rbins; j++) {
            
            // check if bins are within 3*force resolution of halo center, and downweight if so
            if (h->bins[j] < 3*FORCE_RES) {
                h->weights[j] = 0.1;
                n_downweights++;
            }

            // add up total number of particles out to rvir
            if (j > 0) {
                if (h->bins[j] <= h->radius) np_rvir+=h->ppbin;
                else {
                    if (rvir_bin == -1) {
                        rvir_bin = j;
                        h->newf.np_rvir = np_rvir;
                        h->newf.nb_rvir = rvir_bin;
                    }
                }
            }

            // break if bins no longer downweighted and rvir has been reached
            if (h->weights[j] == 1 && h->bins[j] > h->radius) break;
        }

        // let's not fit halos with less than a certain number of particles
        // or those with 2 or less bins within rvir (rvir_bin is next bin position outside of rvir)
        if (np_rvir < 100 || rvir_bin < 4) {

            skipped++;
            continue;
        }

        if (!debug) {
            if (break_after > -1) {
                if ((i-skipped) % dot_freq == 0) {
                    std::cout << counter++ << ".";
                    std::cout.flush();
                }
            }

            else if (i % dot_freq == 0) {
                std::cout << counter++ << ".";
                std::cout.flush();
            }
        }

        // we will fit this profile, so save pointer to it in a new vector (so we can easily access all fit profiles later)
        fitted_profiles.push_back(h);

        // now we need to provide an initial guess for rs
        double rs=h->radius/12.;
        int rs_bin = -1;
        double c0;
        double gamma = 2.;

        //for (int j = 1; j < h->n_rbins; j++) {

        //    if (j*h->ppbin > np_rvir/20.) {
        //        rs = h->bins[j];
        //        rs_bin = j;

        //        // don't want rs = 0 -- need to include more particles here
        //        if (rs == 0) continue;
        //        break;
        //    }
        //}

        // and an initial guess for constant c_0 = 16*pi*rho_0
        c0 = np_rvir / (rs*rs*rs * nfw_menc (h->radius, rs));
        
        // actual c0 at r_s
        double c0_check = 0.;//12 * h->ppbin / (pow(.5*(h->bins[rs_bin]+h->bins[rs_bin+1]),3.)-pow(.5*(h->bins[rs_bin-1]+h->bins[rs_bin]),3.));

        if (debug) std::cout << "n_rbins: " << h->n_rbins << ", downweighted bins: " << n_downweights << ", rvir(" << rvir_bin << "): " << h->radius << ", rs_init(" << rs_bin << "): " << rs << ", c0 init: " << c0 << ", and computed c0: " << c0_check << ", difference factor: " << c0_check/c0 << std::endl;

        //df += c0_check/c0;

        // now, we'll compute the fits
        // we have out first guesses for our parameters (e.g. rs, c0, gamma).
        // basic procedure will be to compute chi^2, then chose new values and compute chi^2 again.  repeat until convergence reached.

        std::vector<double> p1(1), p2(2), p3(3);
        std::vector<bound> bounds1(1), bounds2(2), bounds3(3);
        p1[0] = rs; p2[0] = rs; p3[0] = rs;
        p2[1] = c0; p3[1] = c0;
        p3[2] = gamma;

        bounds1[0].lower = 0; bounds1[0].upper = h->radius;
        bounds2[0].lower = 0; bounds2[0].upper = h->radius;
        bounds3[0].lower = 0; bounds3[0].upper = h->radius;

        bounds2[1].lower = 0; bounds2[1].upper = 1000.;
        bounds3[1].lower = 0; bounds3[1].upper = 1000.;

        bounds3[2].lower = 0; bounds3[2].upper = 100.;

        if (break_after == -10) {

            std::ofstream of;

	        of.open("error_hist2D.dat",std::ofstream::out);

            double rs_min = 100, rs_max = 400;
            double c0_min = 0.01, c0_max = 2;
            double steps = 1000.;
            double rs_step = (rs_max-rs_min)/steps;
            double c0_step = (c0_max-c0_min)/steps;

            for (double j = rs_min; j < rs_max; j+=rs_step) {
                for (double k = c0_min; k < c0_max; k+=c0_step) {
                    of << j << " " << k << " " << chi2_nfw(*h,j,k) << std::endl;
                }
            }
            
            of.close();
        }

        int converged;
        clock_t begin_t, end_t;
        double final_chi2;

        ///////////////////////////////////////////////////////////////

        if (debug) std::cout << "\tStarting clock... " << std::endl;

        off_stream(std::cout);

        begin_t = clock();
        
        //converged = doFitGD (h, rs, c0);

        converged = doFitBGD (h, p1, chi2_nfw1, bounds1);

        end_t = clock();
        
        on_stream(std::cout);
        
        total_t += end_t - begin_t;

        if (debug) std::cout << "\tStopped clock. Elapsed CPU time: " << (double(end_t-begin_t)) / CLOCKS_PER_SEC << std::endl;

        rs = p1[0];

        final_chi2 = chi2_nfw1 (*h, p1);
        
        if (debug) std::cout << "\trs_fit: " << rs << ", chi2: " << final_chi2 << ", cnfw: " << h->radius/rs << ", total iterations: " << converged << std::endl;

        h->newf.rs1 = rs;
        h->newf.chi2_1 = final_chi2;
        h->newf.GOF1 = final_chi2 / (h->newf.nb_rvir - 1.);

        ///////////////////////////////////////////////////////////////

        // update guesses based on previous fitting results
        p2[0] = rs;
        p2[1] = np_rvir / (rs*rs*rs * nfw_menc (h->radius, rs));

        if (debug) std::cout << "\tStarting clock... " << std::endl;

        off_stream(std::cout);

        begin_t = clock();
        
        //converged = doFitGD (h, rs, c0);

        converged = doFitBGD (h, p2, chi2_nfw2, bounds2);

        end_t = clock();
        
        on_stream(std::cout);
        
        total_t += end_t - begin_t;

        if (debug) std::cout << "\tStopped clock. Elapsed CPU time: " << (double(end_t-begin_t)) / CLOCKS_PER_SEC << std::endl;

        rs = p2[0];
        c0 = p2[1];

        final_chi2 = chi2_nfw2 (*h, p2);
        
        if (debug) std::cout << "\trs_fit: " << rs << ", c0_fit: " << c0 << ", chi2: " << final_chi2 << ", cnfw: " << h->radius/rs << ", total iterations: " << converged << std::endl;

        h->newf.rs2 = rs;
        h->newf.c02 = c0;
        h->newf.chi2_2 = final_chi2;
        h->newf.GOF2 = final_chi2 / (h->newf.nb_rvir - 2.);

        ///////////////////////////////////////////////////////////////

        // update guesses based on previous fitting results
        p3[0] = rs;
        p3[1] = c0;
        p3[2] = 2.0;

        if (debug) std::cout << "\tStarting clock... " << std::endl;

        off_stream(std::cout);

        begin_t = clock();
        
        //converged = doFitGD (h, rs, c0);

        converged = doFitBGD (h, p3, chi2_nfw3, bounds3);

        end_t = clock();
        
        on_stream(std::cout);
        
        total_t += end_t - begin_t;

        if (debug) std::cout << "\tStopped clock. Elapsed CPU time: " << (double(end_t-begin_t)) / CLOCKS_PER_SEC << std::endl;

        rs = p3[0];
        c0 = p3[1];
        gamma = p3[2];

        final_chi2 = chi2_nfw3 (*h, p3);
        
        if (debug) std::cout << "\trs_fit: " << rs << ", c0_fit: " << c0 << ", gamma_fit: " << gamma << ", chi2: " << final_chi2 << ", cnfw: " << h->radius/rs << ", total iterations: " << converged << std::endl;

        h->newf.rs3 = rs;
        h->newf.c03 = c0;
        h->newf.gamma3 = gamma;
        h->newf.chi2_3 = final_chi2;
        h->newf.GOF3 = final_chi2 / (h->newf.nb_rvir - 3.);

        ///////////////////////////////////////////////////////////////

        //if (converged && of_num < 20) {

        //    std::ofstream of;
        //    setOutputFileName();
        //    std::string of_name = output.substr(0,output.find_last_of('.'))+"_"+str(of_num)+output.substr(output.find_last_of('.'));
        //    of.open(of_name.c_str(),std::ofstream::out);

        //    // now lets output some of these profiles to check
        //    for (int j = 0; j < h->newf.nb_rvir; j++) {

        //        double rho_calc, rho_fit;
        //        if (j == 0) { rho_calc = 0; rho_fit = 0;}
        //        else {
        //            rho_calc = 3 * h->ppbin / (4.*M_PI) / (pow(.5*(h->bins[j]+h->bins[j+1]),3.)-pow(.5*(h->bins[j-1]+h->bins[j]),3.));
        //            rho_fit = c0 / (4.*M_PI) / (h->bins[j]/rs) / pow(1.+h->bins[j]/rs,2.);
        //        }
        //        of << h->bins[j] << " " << h->ppbin*j << " " << nfw_menc(h->bins[j],rs,c0) << " " << rho_calc << " " << rho_fit;
        //        if (j == 0) of << " " << h->radius << " " << rs;
        //        of << std::endl;
        //    }
        //    
        //    of.close();

        //    of_num++;
        //}

        if (rs>0) {
            cnfw += h->radius/rs;
            cnfw_count++;
        }

        if ((break_after > -1) && (i+1-skipped >= break_after)) {
            break;
        }
    }

    if (!debug) std::cout << std::endl;

    // now write output file with new data
    std::ofstream of;
    setOutputFileName();
    of.open(output.c_str(),std::ofstream::out);

    for (int i = 0; i < fitted_profiles.size(); i++) {

        halo_profile * h = fitted_profiles[i];

        of << h->id << " " << h->mass << " " << h->radius << " " << h->newf.rs1 << " " << h->newf.rs2 << " " << " " << h->newf.rs3 << " " << h->newf.c02 << " " << h->newf.c03 << " " << h->newf.gamma3 << " " << h->newf.chi2_1 << " " << h->newf.GOF1 << " " << h->newf.chi2_2 << " " << h->newf.GOF2 << " " << h->newf.chi2_3 << " " << h->newf.GOF3 << " " << h->radius / h->newf.rs2 << " " << h->newf.nb_rvir << std::endl;
    }

    of.close();

    if (!debug) on_stream(std::cout);

    std::cout << "Skipped " << skipped << " halos out of " << (break_after == -1 ? hpp.size() : skipped+break_after) << " checked so far." << std::endl;
    std::cout << "Total fitting CPU time: " << (double(total_t)) / CLOCKS_PER_SEC << std::endl;
    //std::cout << "average difference factor: " << df/200. << std::endl;
    std::cout << "average concentration: " << cnfw/cnfw_count << std::endl;

    return 0;
}


