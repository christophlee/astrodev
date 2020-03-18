#include "sam_am.h"

int doUnmatchedHalos (double ** & data, int num_lines) {

    // lets write some code to output information for each of the unmatched "centrals" from the mock
    // now for indeces:
    int halo_id = 0, x = 1, y = 2, z = 3, mvir = 7, mhost = 11;

    typedef struct {
        double mvir;
        double x;
        double y;
        double z;
    } halo;

    std::vector<halo> unmatched;
    std::vector<unsigned int> matched;

    double ** matched_data;
    int num_lines_2=0, num_fields_2=0;
    std::string matched_file = std::string("matched_halos.dat");

    readFile(matched_file,num_lines_2,num_fields_2,matched_data);

    for (int i = 0; i < num_lines_2; i++) {
        matched.push_back(matched_data[0][i]);
    }

    free(matched_data, num_fields_2);

    int count = 0;
    for (int i = 0; i < num_lines; i++) {
        if (data[mvir][i] == data[mhost][i]) {
            std::vector<unsigned int>::iterator it = find (matched.begin(),matched.end(),data[halo_id][i]);
            if (it == matched.end()) {
                count++;
                halo newhalo;
                newhalo.mvir = data[mvir][i];
                newhalo.x = data[x][i];
                newhalo.y = data[y][i];
                newhalo.z = data[z][i];
                unmatched.push_back(newhalo);
            }
        }
    }

    std::cout << "Number of unmatched \"centrals\" is " << count << std::endl;

    std::ofstream of_mock;

    setOutputFileName();

    of_mock.open(output.c_str(),std::ofstream::out);

    for (int i = 0; i < unmatched.size(); i++) {
        of_mock << log10(unmatched[i].mvir) << " " << unmatched[i].x << " " << unmatched[i].y << " " << unmatched[i].z << std::endl;
    }

    of_mock.close();

    return 0;
}

int doCreateLimitedCatalog (double ** & data, int num_lines, int num_fields, int x, double minMass, double maxMass, std::string descriptor) {

    std::ofstream of_lim;

    setOutputFileName(x);

    of_lim.open(output.c_str(),std::ofstream::out);

    for (int i = 0; i < num_lines; i++) {

        // here is where we decide how the catalog is being limited
        //
        // check mass limits
        if (minMass >= 0) {
            if (data[x][i] < minMass) continue;
        }
        if (maxMass >= 0) {
            if (data[x][i] > maxMass) continue;
        }
        if (descriptor.length()) {
            if (descriptor == "sam centrals") {

                // if this is sam output, we need to check gal id
                if (data[4][i] != 1) continue;
            }
            else if (descriptor == "am centrals") {

                // for age matching, compare mvir to mhost
                if (data[7][i] != data[11][i]) continue;
            }
            else if (descriptor == "sam sats") {

                if (data[4][i] == 1) continue;
            }
            else if (descriptor == "am sats") {

                if (data[7][i] == data[11][i]) continue;
            }
            else {

                std::cerr << "Error: unrecognized descriptor in doCreateLimitedCatalog." << std::endl;
            }
        }
     
        // output to file
        for (int j = 0; j < num_fields; j++) {
            if (j) of_lim << " " << data[j][i];
            else of_lim << data[j][i];
        }
        of_lim << std::endl;
    }

    of_lim.close();

    return 0;
}

int doCalcQuiescentFractions (double ** & data, int num_lines, int num_fields, int x, int y) {

    double sSFR_SAM, sSFR_AM;
    int both_quiescent = 0, both_sf = 0, sam_quiescent = 0, sam_sf = 0;

    std::ofstream * of_agree = new std::ofstream;
    std::ofstream * of_disagree = new std::ofstream;
    std::ofstream * temp;

    // modified function to instead output updated catalog with
    // an agreement parameter for each halo rather than separate
    // catalogs                                               //
    std::ofstream * of_all = new std::ofstream;               //
    (*of_all).open(output.c_str(),std::ofstream::out);        //
    int agreement_param = 0;                                  //
    ////////////////////////////////////////////////////////////

    std::string output_a = std::string("matched_ssfr_agree.dat");
    std::string output_d = std::string("matched_ssfr_disagree.dat");

    std::cout << "starting calculations" << std::endl;
    (*of_agree).open(output_a.c_str(),std::ofstream::out);
    (*of_disagree).open(output_d.c_str(),std::ofstream::out);

    std::cout << "starting calculations" << std::endl;

    // lets use the SAM data as our starting point
    // determine whether a galaxy is star forming or quiescent using the criterior
    // that a quiescent galaxy has sSFR < -11
    for (int i = 0; i < num_lines; i++) {
        sSFR_SAM = data[x][i];
        sSFR_AM = data[y][i];
        if (sSFR_SAM < -11.0) {
            if (sSFR_AM < -11.0) {
                both_quiescent++;
                agreement_param = 1;
                //temp = of_agree;
            }
            else {
                sam_quiescent++;
                agreement_param = 4;
                //temp = of_disagree;
            }
        }
        else {
            if (sSFR_AM < -11.0) {
                agreement_param = 2;
                sam_sf++;
                //temp = of_disagree;
            }
            else {
                agreement_param = 3;
                both_sf++;
                //temp = of_agree;
            }
        }

        temp = of_all;

        for (int j = 0; j < num_fields; j++) {
            if (!j) (*temp) << data[j][i];
            else (*temp) << " " << data[j][i];
        }
        (*temp) << " " << agreement_param << std::endl;
    }

    (*of_all).close();

    (*of_agree).close();
    (*of_disagree).close();

    delete of_all;
    delete of_agree;
    delete of_disagree;

    std::cout << "Quiescent Fractions are as follows:" << std::endl;
    std::cout << "Both Star Forming: " << (double) both_sf / (double) num_lines << std::endl;
    std::cout << "Both Quiescent: " << (double) both_quiescent / (double) num_lines << std::endl;
    std::cout << "SAM Star Forming, Age Matching Quiescent: " << (double) sam_sf / (double) num_lines << std::endl;
    std::cout << "SAM Quiescent, Age Matching Star Forming: " << (double) sam_quiescent / (double) num_lines << std::endl;

    return 0;
}

int doCalcStellarMassFractions (double ** & data, int num_lines, int num_fields, int x, int y) {

    std::ofstream of;

    setOutputFileName();

    of.open(output.c_str(),std::ofstream::out);

    double totalx = 0.0, totaly = 0.0;

    for (int i = 0; i < num_lines; i++) {
        totalx += data[x][i];
        totaly += data[y][i];
    }

    of << totalx/(double)num_lines << " " << totaly/(double)num_lines << std::endl;

    /*int skipped = 0;

    for (int i = 0; i < num_lines; i++) {
        int pid = data[5][i];
        // only do centrals
        if (pid == -1) {
            double c_nfw = data[63][i];
            unsigned int halo_id = data[1][i];
            double mvir = log10(data[10][i]);
            double mpeak = log10(data[55][i]);
            double lambda = data[26][i];
            if (pow(10.,mvir)/pow(10.,mpeak) > 0.6) {
                skipped++;
                continue;
            }
            of << (unsigned int)halo_id << " " << mvir << " " << lambda << " " << c_nfw;
            for (int j = 0; j < 17; j++) {
                of << " " << data[64+j][i];
            }
            of << std::endl;
        }
        //of << data[x][i] << " " << pow(10.,data[y][i])/pow(10.,data[x][i]) << std::endl;
        //of << data[x][i] << " " << 4.1/data[y][i] << std::endl;

    }

    std::cout << "Skipped halos due to MVIR / MPEAK < 0.8: " << skipped << std::endl;
    */
    of.close();

    return 0;

}

