#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <map>
#include <cmath>
#include <random>
#include <chrono>

///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
////////////////////////////////////////////////////
//////////////////////////////////////////////
/////////////////////////////////////////
/////////////////////////////////////
// DEFINE CONSTANTS //////////////
////////////////////////////////
///////////////////////////////
#define GALPROP 0
#define SAM 1
#define MOCK 2

#define galprop_num_fields 39
#define halos_num_fields 25
#define mock_num_fields 13

// halos.dat //
#define halos_halo_id_index 0
#define halos_halo_id_nbody_index 1
#define halos_m_vir_index 3
#define halos_c_NFW_index 7
#define halos_spin_index 9
#define halos_mstar_index 11
#define halos_sfr_index 14

// mock //
#define mock_halo_id_nbody_index 0
#define mock_x_index 1
#define mock_y_index 2
#define mock_z_index 3
#define mock_mvir_index 7
#define mock_m_star_index 9
#define mock_gr_color_index 10
#define mock_ssfr_index 10
#define mock_mhost_index 11

// galprop.dat //
#define galprop_halo_id_index 0
#define galprop_gal_id_index 1
#define galprop_mhalo_index 6
#define galprop_mstar_index 9
#define galprop_mbh_index 15
#define galprop_sfr_ave_index 24
#define galprop_x_index 32
#define galprop_y_index 33
#define galprop_z_index 34

typedef struct {
	double stellar_mass = 0;
	double sfr = 0;
	double sfr_ave = 0;
    double spin = 0;
    double c_NFW = 0;
	double m_vir = 0;
    double mbh = 0;
	double gr_color = 0;
	double x = 0;
	double y = 0;
	double z = 0;
} model;

typedef struct {
	model am;
	model sam;
} halo_object;

double unit_random () {

    return (double) std::rand() / (double) RAND_MAX;
}

double random_normal (double mean, double stdDev) {

    double w = 1.0, x1, x2, y1, y2;

    while (w >= 1.0) {
        x1 = 2*unit_random()-1.0;
        x2 = 2*unit_random()-1.0;
        w = x1*x1 + x2*x2;
    }

    w = sqrt ( (-2.0 * log(w)) / w);
    y1 = x1*w;
    y2 = x2*w;

    // take the first number and adjust std dev and mean as necessary;
    y1 *= stdDev;
    y1 += mean;

    return y1;
}

// binning function to take in a 1D stream of data and return a 2D binned array for outputing
double ** bin_data (double * data, unsigned int size, int num_bins) {

	double min, max, bin_size;

	// set max/min
	for (int i = 0; i < size; i++) {
		if (!i) {
			min = data[i];
			max = data[i];
		}
		else {
			if (data[i] < min && isfinite(data[i])) min = data[i];
			if (data[i] > max && isfinite(data[i])) max = data[i];
		}
	}

	// for now, we'll allow bin_size to be a double (non integer)
	bin_size = (max-min)/(double)num_bins;
	
	std::cout << min << " " << max << " " << bin_size << std::endl;

	// allocate space for the binned array
	double ** binned_data = new double * [2];
	for (int i = 0; i < 2; i++) binned_data[i] = new double [size];

	// initialize
	// here, the value of the bin means where the bin starts, so the hits for that bin fall in the
	// range of value + bin_size
	double sum = min;
	for (int i = 0; i < num_bins; i++) {
		binned_data[0][i] = sum;
		binned_data[1][i] = 0;
		sum += bin_size;
	}

	// populate bins, catching bad values
	for (int i = 0; i < size; i++) {
		int index = floor((data[i]-min)/bin_size);
		if (index > num_bins) index = num_bins-1;
		if (index < 0) index = 0;
		binned_data[1][index]++;
	}

	return binned_data;
}

// binning function to take in a 1D stream of data and write binned output to file, returning success/fail
int write_histogram_data (std::string file_name, double * data, unsigned int size, int num_bins) {

	std::ofstream of;

	of.open (file_name.c_str(), std::ofstream::out);
	
	// function returns new array with binned data
	double ** binned_data = bin_data (data, size, num_bins);

	// written as "x y" data columns --> x is the bins, y the abundances
	for (int i = 0; i < num_bins; i++)
		of << binned_data[0][i] << " " << binned_data[1][i] << std::endl;
	
	of.close();

	// free data
	for (int i = 0; i < 2; i++) delete [] binned_data[i];
	delete [] binned_data;

	return 1;
}

int main(int argc, char ** argv) {
	
	std::ifstream * datafiles;
	std::string tmp_str;
	int numFiles = 3;
	std::string runs_dir = "/pfs/chtlee/runs/sam_bp_8_5_15";

	datafiles = new std::ifstream[numFiles];

	if (argc > 1) runs_dir = "/pfs/chtlee/runs/sam_bp_8_5_15/" + std::string(argv[1]);

	datafiles[GALPROP].open((runs_dir + "/galprop.dat").c_str(),std::ifstream::in);
	datafiles[SAM].open((runs_dir + "/halos.dat").c_str(),std::ifstream::in);
	//datafiles[MOCK].open("Mr19_age_distribution_matching_mock.dat",std::ifstream::in);
	//datafiles[MOCK].open("sm_gr_fiducial_mock.dat",std::ifstream::in);
	datafiles[MOCK].open("/pfs/chtlee/catalogs/age_matching/sm9.8_age_matching_SFR_mock.dat",std::ifstream::in);
		
	std::cout << "Now reading input file headers " << std::endl;

	for (int file = 0; file < numFiles; file++) {
		int position=0;
		//std::cout << "File " << file << std::endl;	
		while(getline(datafiles[file],tmp_str)) {
			//std::cout << tmp_str << std::endl;
			if (tmp_str[0] != '#') {
				datafiles[file].seekg(position);
				break;
			}
			position = datafiles[file].tellg();
		}
	}

	std::map<unsigned int, halo_object*>* halos = new std::map<unsigned int, halo_object*>();
	std::map<unsigned int, unsigned int>* halo_id_map = new std::map<unsigned int, unsigned int>();
	double *** data = new double**[numFiles];
	int num_fields[numFiles];
	int num_halos[numFiles];

	num_fields[GALPROP] = galprop_num_fields;	// galprop.dat
	num_fields[SAM] = halos_num_fields;		// halos.dat
	num_fields[MOCK] = mock_num_fields;		// mock catalog

	// scan through files to check how many halos are in each
	for (int i = 0; i < numFiles; i++) {
		num_halos[i] = 0;
		int position = datafiles[i].tellg();
		datafiles[i].seekg (0, datafiles[i].end);
		int file_eof = datafiles[i].tellg();
		datafiles[i].seekg (position);
		while (true) {
			if (datafiles[i].tellg() == file_eof) break;
			getline(datafiles[i],tmp_str);
			num_halos[i]++;
		}
		datafiles[i].seekg(position);
	}

	for (int j = 0; j < numFiles; j++) {
		data[j] = new double*[num_fields[j]];
		for (int i = 0; i < num_fields[j]; i++)
			data[j][i] = new double[num_halos[j]];
	}

	std::cout << "Now reading and storing data into arrays" << std::endl;

	for (int k = 0; k < numFiles; k++) {
		//std::cout << "File " << k << ": " << datafiles[k].tellg() << std::endl;
		int i = 0;
		while (datafiles[k].good()) {
			if (i < 20) std::cout << "Halo " << i << ": ";
			for (int j = 0; j < num_fields[k]; j++) {
				double tmp_d;
				char tmp_word[256];
				datafiles[k] >> tmp_word;
				tmp_d = std::strtod(tmp_word,NULL);
				data[k][j][i] = tmp_d;
				if (i < 20) std::cout << tmp_d << " ";
			}
			if (i < 20) std::cout << std::endl;
			i++;
		}
	}

	// at this point all the data is stored in the array "data"

	// step through data and match up halo properties in halo map

	// sam properties
	for (int i = 0; i < num_halos[SAM]; i++) {
			halo_object * halo;

			halo = (*halos)[data[SAM][halos_halo_id_nbody_index][i]];
	
			// create new object to store in map if it doesn't already exist
			if (!halo) halo = new halo_object;

			// update values
			halo->sam.m_vir = data[SAM][halos_m_vir_index][i];
			halo->sam.sfr = data[SAM][halos_sfr_index][i];
            halo->sam.c_NFW = data[SAM][halos_c_NFW_index][i];
            halo->sam.spin = data[SAM][halos_spin_index][i];

			(*halos)[data[SAM][halos_halo_id_nbody_index][i]] = halo;

			// update halo_id_map to map halo_id_nbody to halo_id
			(*halo_id_map)[data[SAM][halos_halo_id_index][i]] = data[SAM][halos_halo_id_nbody_index][i];
	}

	// mock properties
/*	for (int i = 0; i < num_halos[MOCK]; i++) {

			halo_object * halo;

			halo = (*halos)[data[MOCK][mock_halo_id_nbody_index][i]];

			if (!halo) halo = new halo_object;

			halo->am.m_vir = log10(data[MOCK][mock_mvir_index][i]);
			halo->am.x = data[MOCK][mock_x_index][i];
			halo->am.y = data[MOCK][mock_y_index][i];
			halo->am.z = data[MOCK][mock_z_index][i];
			halo->am.stellar_mass = data[MOCK][mock_m_star_index][i];
			//halo->am.gr_color = data[MOCK][mock_gr_color_index][i];
			//halo->am.sfr = pow(10.,data[MOCK][mock_ssfr_index][i])*pow(10.,halo->am.stellar_mass);
			halo->am.sfr = data[MOCK][mock_ssfr_index][i];

			(*halos)[data[MOCK][mock_halo_id_nbody_index][i]] = halo;
	}*/

    int satellites = 0;
    int centrals = 0;

    int output_freq = ceil(num_halos[GALPROP]/10.);

	// galprop properties
	for (int i = 0; i < num_halos[GALPROP]; i++) {
			halo_object * halo;
			unsigned int halo_id_nbody;

            if (i % output_freq == 0) std::cout << i/output_freq << std::endl;

			// check if galaxy is a central, other wise let's skip it
			if ((long)data[GALPROP][galprop_gal_id_index][i] != 1) continue;

			// retrieve nbody id
			halo_id_nbody = (*halo_id_map)[data[GALPROP][galprop_halo_id_index][i]];

			if (!halo_id_nbody) {
				std::cerr << "ERROR: halo_id returned null halo_id_nbody from map: " << data[GALPROP][galprop_halo_id_index][i] << std::endl;
				continue;
			}

			halo = (*halos)[halo_id_nbody];

			if (!halo) {
				std::cerr << "ERROR: halo_id_nbody returned null halo from map when it should not" << std::endl;
				continue;
			}

			if (halo->sam.x != 0) {
				std::cerr << "ERROR: it appears properties of central galaxy in halo have already been set.  Proceeding to override..." << std::endl;
                std::cerr << "halo_id_nbody: " << (unsigned int) halo_id_nbody << std::endl;
                std::cerr << "sam.cnfw: " << halo->sam.c_NFW << std::endl;
                std::cerr << "sam.x: " << halo->sam.x << std::endl;
                std::cerr << "sam.stellar_mass " << halo->sam.stellar_mass << std::endl;
                std::cerr << "new sam.x: " << data[GALPROP][galprop_x_index][i] << std::endl;
                std::cerr << "new sam.stellar_mass: " << data[GALPROP][galprop_mstar_index][i] << std::endl;
                std::cerr << "i: " << i << std::endl;
			}

			halo->sam.x = data[GALPROP][galprop_x_index][i];
			halo->sam.y = data[GALPROP][galprop_y_index][i];
			halo->sam.z = data[GALPROP][galprop_z_index][i];
			halo->sam.sfr_ave = data[GALPROP][galprop_sfr_ave_index][i];
			halo->sam.stellar_mass = data[GALPROP][galprop_mstar_index][i];
            halo->sam.mbh = data[GALPROP][galprop_mbh_index][i];

		//	// using temporaily for sSFR
		//	if (halo->sam.sfr_ave == 0) {
        //        halo->sam.sfr = -15.0;
        //        halo->sam.sfr_ave = 0.00000001;
        //    }
		//	else halo->sam.sfr = log10 (halo->sam.sfr_ave / pow(10.,halo->sam.stellar_mass));
	}

	std::map<unsigned int, halo_object*>::iterator it = halos->begin();

	// create output file to write sfr ///////////////////////
/*	std::ofstream of_sfr;

	of_sfr.open (("comparison_sfr"+std::string(argv[1]).substr(7)+".dat").c_str(),std::ofstream::out);

	it = halos->begin();

	while (it != halos->end()) {
		if (it->second->sam.m_vir && it->second->am.m_vir)
			of_sfr << it->second->sam.sfr << " " << it->second->am.sfr << " " << it->second->sam.sfr_ave << " " << pow(10.,it->second->am.sfr)*pow(10.,it->second->am.stellar_mass) << std::endl;
		it++;
	}

	of_sfr.close();*/
	//////////////////////////////////////////////////////////

	// create output file to write SAM - AM comparisons //////
	std::ofstream of_comp;

	of_comp.open (("comparison"+std::string(argv[1]).substr(7)+".dat").c_str(),std::ofstream::out);

	it = halos->begin();

	int count = 0;
    double sSFR, SFR;

	while (it != halos->end()) {
		if ((it->second->sam.stellar_mass != 0) && true) {//(it->second->am.m_vir != 0)) {
		// 	nbody_halo_id		sam mstar	am mstar	sam mvir	am mvir	
			of_comp << it->first << " " << it->second->sam.stellar_mass << " " /*<< it->second->am.stellar_mass << " "*/ << it->second->sam.m_vir/* << " " << it->second->am.m_vir*/;

            SFR = it->second->sam.sfr_ave;

            // calculate sam sSFR and "corrected" SFR
            if (it->second->sam.sfr_ave == 0.0) {
                sSFR = pow (10., random_normal (-12.0,0.2));
                SFR = sSFR * pow (10., it->second->sam.stellar_mass);
            }

            else {
                sSFR = SFR / pow (10., it->second->sam.stellar_mass);
            }

            // calculate sSFR correlation agreement parameter
            /*int agreement_param = 0;
            double sSFR_AM = it->second->am.sfr;
            if (log10(sSFR) < -11.0) {
                if (sSFR_AM < -11.0)    agreement_param = 1;
                else                    agreement_param = 4;
            }
            else {
                if (sSFR_AM < -11.0)    agreement_param = 2;
                else                    agreement_param = 3;
            }*/

		//	sam sfr	    sam "corrected" SFR     sam sSFR	am sSFR     c_NFW   a_form  mbh     spin    agreement_param
			of_comp  << " " << it->second->sam.sfr_ave << " " << SFR << " " << log10(sSFR) << " " << /*it->second->am.sfr << " " <<*/ it->second->sam.c_NFW << " " << 4.1/it->second->sam.c_NFW << " " << it->second->sam.mbh << " " << it->second->sam.spin /*<< " " << agreement_param */<< std::endl;
		}
		it++;
	}

	of_comp.close();

    /*
	//////////////////////////////////////////////////////////

	// create output file for halo mass vs stellar mass //////
	// fraction for all centrals vs satalites
	std::ofstream of_mh_ms_cen;
	std::ofstream of_mh_ms_sat;

	of_mh_ms_cen.open (("mh_ms_cen"+std::string(argv[1]).substr(7)+".dat").c_str(),std::ofstream::out);
	of_mh_ms_sat.open (("mh_ms_sat"+std::string(argv[1]).substr(7)+".dat").c_str(),std::ofstream::out);

	for (int i = 0; i < num_halos[GALPROP]; i++) {
		int halo_id_nbody = (*halo_id_map)[data[GALPROP][galprop_halo_id_index][i]];

		// get full halo mass
		double m_halo = (*halos)[halo_id_nbody]->sam.m_vir;

		// centrals
		//if (data[GALPROP][galprop_gal_id_index][i] == 1) {
		//	of_mh_ms_cen << data[GALPROP][galprop_mhalo_index][i] << " " << data[GALPROP][galprop_mstar_index][i]/data[GALPROP][galprop_mhalo_index][i] << std::endl;
		//}
		//// satelites
		//else {
		//	of_mh_ms_sat << data[GALPROP][galprop_mhalo_index][i] << " " << data[GALPROP][galprop_mstar_index][i]/data[GALPROP][galprop_mhalo_index][i] << std::endl;
		//}

		// centrals
		if (data[GALPROP][galprop_gal_id_index][i] == 1) {
			of_mh_ms_cen << m_halo << " " << data[GALPROP][galprop_mstar_index][i]/m_halo << std::endl;
		}
		// satelites
		else {
			of_mh_ms_sat << m_halo << " " << data[GALPROP][galprop_mstar_index][i]/m_halo << std::endl;
		}
	}

	of_mh_ms_cen.close();
	of_mh_ms_sat.close();
	
	//////////////////////////////////////////////////////////

	// write output file SAM galaxies / halos ////////////////

	std::ofstream of_sam;

	of_sam.open(("sam_prop"+std::string(argv[1]).substr(7)+".dat").c_str(),std::ofstream::out);

	for (int i = 0; i < num_halos[SAM]; i++) {

		// lets skip the halos with mstar == -99 for now
		if (data[SAM][halos_mstar_index][i] < 0) {
			//std::cout << std::string(argv[1]).substr(7) << ": mstar = " << data[SAM][halos_mstar_index][i] << std::endl;
 			continue;
		}

		// writing:
		// 	m_vir	SFR	sSFR	a_form	   m_star
		of_sam << data[SAM][halos_m_vir_index][i] << " " << data[SAM][halos_sfr_index][i] << " " << log10(data[SAM][halos_sfr_index][i]/pow(10.,data[SAM][halos_mstar_index][i])) << " " << 4.1/data[SAM][halos_c_NFW_index][i] << " " << data[SAM][halos_mstar_index][i] << std::endl;
	}
	
	of_sam.close();

	//////////////////////////////////////////////////////////

	// write output file for galaxies ////////////////

	std::ofstream of_gal;

	of_gal.open(("gal_prop"+std::string(argv[1]).substr(7)+".dat").c_str(),std::ofstream::out);

	for (int i = 0; i < num_halos[GALPROP]; i++) {
        // do a few calculations on sfr to get numbers more usable
        double sfr, sSFR;
        sfr = data[GALPROP][galprop_sfr_ave_index][i];
        if (sfr == 0) {
            //sSFR = pow (10.0, -11.99);
            sSFR = pow (10., random_normal (-12.0,0.2));
            sfr = sSFR * pow (10., data[GALPROP][galprop_mstar_index][i]);
        }

        else {
            sSFR = sfr / pow (10., data[GALPROP][galprop_mstar_index][i]);
        }

		// mhalo    mstar 	sfr_ave     sSFR      gal_id
		of_gal << data[GALPROP][galprop_mhalo_index][i] << " " << data[GALPROP][galprop_mstar_index][i] << " " << log10(sfr) << " " << log10(sSFR) << " " << data[GALPROP][galprop_gal_id_index][i] << std::endl;
	}

	of_gal.close();
*/

	//////////////////////////////////////////////////////////

	// write output file for mock ////////////////////////////
//
//	std::ofstream of_mock;
//
//	of_mock.open(("mock_prop"+std::string(argv[1]).substr(7)+".dat").c_str(),std::ofstream::out);
//
//	for (int i = 0; i < num_halos[MOCK]; i++) {
//		// m_vir	mstar		SFR		sSFR		central?1:0
//		of_mock << log10(data[MOCK][mock_mvir_index][i]) << " " << data[MOCK][mock_m_star_index][i] << " ";
//		of_mock << pow(10.,data[MOCK][mock_ssfr_index][i])*pow(10.,data[MOCK][mock_m_star_index][i]);
//		of_mock << " " << data[MOCK][mock_ssfr_index][i] << " ";
//		// print "1" if there is a matching sam halo, else "0"
//		halo_object * halo = (*halos)[data[MOCK][mock_halo_id_nbody_index][i]];
//		if (halo->sam.m_vir != 0.0) of_mock << "1";
//		else	of_mock << "0";
//		of_mock << std::endl;
//	}
//
//	of_mock.close();
//
	//////////////////////////////////////////////////////////	

	// write a sSFR histogram output file ////////////////////
	//double * sSFR = new double[num_halos[GALPROP]];

	//for (int i = 0; i < num_halos[GALPROP]; i++)
	//	sSFR[i] = log10(data[GALPROP][galprop_sfr_ave_index][i] / pow(10.,data[GALPROP][galprop_mstar_index][i]));

	//write_histogram_data (std::string("sSFR_hist.dat"),sSFR,num_halos[GALPROP],50);

	//delete [] sSFR;
	//////////////////////////////////////////////////////////

	// clean up and free memory

	(*halo_id_map).clear();
	delete halo_id_map;

	(*halos).clear();
	delete halos;

	for (int i = 0; i < numFiles; i++) datafiles[i].close();

	for (int j = 0; j < numFiles; j++) {
		for (int i = 0; i < num_fields[j]; i++)
			delete [] data[j][i];
		delete [] data[j];
	}
	delete [] data;
	delete [] datafiles;

	return 1;
}

