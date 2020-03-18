#include "tree.h"

int num_lines_loc = 0, num_fields_loc = 0;
double ** loc_data = NULL;
std::unordered_map<long int, long int> halo_id;

//====================================================
// doFullTreeAnalysis:
// Takes halo catalog as input, produces tree-analyzed halo
// catalog as output, with new fields appended to original catalog
// data:    halo catalog to analyze
// num_lines:
//          number of lines in the halo catalog
// num_fields:
//          number of fields each halo has in the catalog
// mmp:
//          true indicates to only analyze the most massive progenitor branch.
//          if false, all halos in merger tree will be considered for analysis
//====================================================
int doFullTreeAnalysis (double ** &data, int num_lines, int num_fields, bool mmp) {

    std::cout << "Running full tree analysis pipeline..." << std::endl;

    // need to determine chunk size of halos to assign to each processor initially
    unsigned int CHUNK_SIZE = (int)floor(num_lines/100.);

    // keep CHUNK_SIZE greater than a minimum threshold and less than a max threshold
    if (CHUNK_SIZE < 100) CHUNK_SIZE = 100;
    if (CHUNK_SIZE > 10000) CHUNK_SIZE = 10000;

    int max_i = (int)ceil(num_lines/(double)CHUNK_SIZE);
    
    double dot_freq = max_i/100.;

    int chunks_completed = 0, counter = 0;

    std::cout << "Dividing catalog into " << max_i << " chunks of roughly " << CHUNK_SIZE << " halos each." << std::endl;
    std::cout << "Number of cores available: " << omp_get_num_procs() << "/" << omp_get_max_threads() << std::endl;

    // do this so we can modify output file name for sending to treewalk etc
    setOutputFileName();

    std::string loc_path = "/zang/chtlee/pfs_backup/workspace/locations.mod.bin";

    // load locations into memory only once
    readBinaryFile(loc_path,num_lines_loc,num_fields_loc,loc_data);
   
    off_stream (std::cout);

    for (int i = 0; i < num_lines_loc; i++) halo_id[loc_data[0][i]] = i;

#pragma omp parallel for schedule(dynamic) shared(data,num_lines,num_fields,mmp,CHUNK_SIZE,max_i,loc_data,halo_id,chunks_completed,counter)
    for (int i = 0; i < max_i; i++) {

        // now need to create a new data structure to hold temporary chunk for each
        // core to work on and pass to treeWalk, etc.
       
        double ** data_chunk;
        int num_lines_chunk = CHUNK_SIZE;

        data_chunk = (double **) malloc (num_fields * sizeof(double *));

        if (CHUNK_SIZE*(i+1) > num_lines) num_lines_chunk = num_lines-i*CHUNK_SIZE; 

        for (int j = 0; j < num_fields; j++) {

            data_chunk[j] = (double *) malloc (num_lines_chunk * sizeof (double));
        }

        for (int j = 0; j < num_fields; j++) {
            for (int k = 0; k < num_lines_chunk; k++) {
                data_chunk[j][k] = data[j][i*CHUNK_SIZE+k];
            }
        }

        // create new output file unique to this iteration for each step of pipeline
        int ind = output.find_last_of('.');
        
        std::string outfile_name = output.substr(0,ind)+"_"+str(i)+".bin2";

        //if (!omp_get_thread_num()) std::cerr << "Calling treeWalk..." << std::endl;

        // find merger trees for current chunk and write to unique bin2 output file
        doTreeWalk (data_chunk, num_lines_chunk, num_fields, mmp, outfile_name);

        double ** tree_data = NULL; int num_fields_i = 0, num_lines_i = 0;

        //if (!omp_get_thread_num()) std::cerr << "Reading treeWalk output file..." << std::endl;

        // now we need to read in the unique bin2 treeWalk output file
        readBinary2File (outfile_name, num_lines_i, num_fields_i, tree_data);

        // now we can pass this merger tree on to be tree analyzed and output again
        // as unique output file (but we can overwrite treewalk output file)
        doTreeAnalysis (tree_data, num_lines_i, num_fields_i, outfile_name);

        // free tree_data so we can use it again
        free (tree_data, num_lines_i);

        // now read tree analysis to pass to selectZ
        readBinary2File (outfile_name, num_lines_i, num_fields_i, tree_data);

        outfile_name = output.substr(0,ind)+"_"+str(i)+".bin";

        // now do selectZ on tree-analyzed bin2 catalog. overwrite file to become
        // z=0 bin catalog with all tree properties
        doSelectZ (tree_data, num_lines_i, num_fields_i, 0, outfile_name);

        free (tree_data, num_lines_i);
        free (data_chunk, num_fields);

        // cleanup intermediate output file
        std::remove ((output.substr(0,ind)+"_"+str(i)+".bin2").c_str());

        // check if stdout is turned off and if so do thread by thread progress
        if (std::cout.fail()) {

            #pragma omp critical
            {
            chunks_completed++;

            while (chunks_completed > counter*dot_freq) {

                std::cerr << counter << ".";
                std::cerr.flush();
                counter++;
            }
            }
        }
    }

    if (std::cout.fail()) std::cerr << " complete." << std::endl;

    // free up locations data now since we're done with it
    free (loc_data,num_fields_loc);

    on_stream (std::cout);

    double ** data2;
    int num_lines2 = num_lines, num_fields2 = 0;

    // now we need to loop through and merge all new z=0 catalog fragments
    for (int i = 0; i < max_i; i++) {

        double ** data_chunk; int num_lines_chunk = 0, num_fields_chunk = 0;

        int ind = output.find_last_of('.');

        std::string infile_name = output.substr(0,ind)+"_"+str(i)+".bin";

        readBinaryFile (infile_name, num_lines_chunk, num_fields_chunk, data_chunk);

        // cleanup file once its been read in and we no longer need it
        std::remove (infile_name.c_str());

        // still need to set num fields for data2 and allocate memory
        if (i == 0) {

            num_fields2 = num_fields_chunk;

            data2 = (double **) malloc (num_fields2 * sizeof (double *));

            for (int j = 0; j < num_fields2; j++) {

                data2[j] = (double *) malloc (num_lines2 * sizeof(double));
            }
        }

        for (int j = 0; j < num_fields_chunk; j++) {
            for (int k = 0; k < num_lines_chunk; k++) {

                data2[j][i*CHUNK_SIZE+k] = data_chunk[j][k];
            }
        }
    }

    doCatalogMerge (data2, num_lines2, num_fields2, data, num_lines, num_fields);

    free (data2, num_fields2);

    return 0;
}

//====================================================
// doTreeWalk:
// Takes a halo catalog and looks up merger trees for each halo
// in the catalog.  The merger trees listed one after another in one
// output file.  The output file can be written in ASCII or binary format.
// data:    halo catalog data of halos to look up merger trees for
// num_lines:
//          number of lines in the halo catalog
// num_fields:
//          number of fields each halo has in the catalog
// mmp:
//          true indicates to only output the most massive progenitor branch.
//          if false, all halos in merger tree will be output
//====================================================
int doTreeWalk (double ** &data, int num_lines, int num_fields, bool mmp, std::string outfile_name) {

    // we're going to need to read in the locations file
    std::string loc_path = "/zang/chtlee/pfs_backup/workspace/locations.mod.bin";

    // true unless locations were read externally, in which case we don't want
    // to free before returning from treewalk
    bool free_loc = true;

    //int num_lines_loc = 0, num_fields_loc = 0;
    //double ** loc_data = NULL;

    // check if locations have already been loaded into memory (could be case if called
    // from fullTreeAnalysis)
    if (loc_data == NULL)   readBinaryFile(loc_path,num_lines_loc,num_fields_loc,loc_data);
    else free_loc = false;

    if (!halo_id.size())    for (int i = 0; i < num_lines_loc; i++) halo_id[loc_data[0][i]] = i;

    // loc_data[0][i] = halo_id
    // loc_data[1][i] = file_id
    // loc_data[2][i] = offset
    // loc_data[3][i] = compressed filename (e.g. tree_0_1_3.dat -> 013)    

    // let's make use of a map to look up our halos in the loc data
    //std::unordered_map<long int, long int> halo_id;
    //for (int i = 0; i < num_lines_loc; i++) halo_id[loc_data[0][i]] = i;

    std::ifstream mtree;    // merger tree file to read tree from
    std::string base = std::string("/zang/pbehrooz/Bolshoi_Planck/trees/");

    // how many fields per halo in merger trees?
    mtree.open((base+"tree_0_0_0.dat").c_str(), std::ifstream::in);

    // skip header
    std::string tmp_str;
    getline(mtree,tmp_str);
    while (tmp_str[0] == '#') getline(mtree,tmp_str);
    
    // we need to get 2 more lines to get to a halo line
    getline(mtree,tmp_str);
    getline(mtree,tmp_str);

    // now we have first halo data in tmp_str. find # fields.
    std::istringstream line(tmp_str);
    int num_fields_mtree = 0;
    while(line >> tmp_str) num_fields_mtree++;
    mtree.close();

    std::cout << "Num fields in merger trees is " << num_fields_mtree << std::endl;
    std::cout.flush();

    // now let's allocate some space for these trees.
    // let's start with 200 lines per merger tree. if this space is used up, we
    // can write it to disk and start refilling the array
    int CHUNKSIZE;
    if (mmp) {
        CHUNKSIZE = 200*num_lines;
        // however, keep it limited to 1GB if the above would go over
        if (CHUNKSIZE*num_fields_mtree*sizeof(double) > 1024*1024*1024) {
            CHUNKSIZE = (1024*1024*1024)/(sizeof(double)*num_fields_mtree);
        }
    }

    else CHUNKSIZE = (1024*1024*512)/(sizeof(double)*num_fields_mtree);

    //if (!omp_get_thread_num()) {
    //    std::cerr << "CHUNKSIZE = " << CHUNKSIZE << std::endl;
    //    std::cerr << "num_fields_mtree = " << num_fields_mtree << std::endl;
    //    std::cerr << "num_lines = " << num_lines << std::endl;
    //}

    long long mtree_size = CHUNKSIZE;

    double ** mtree_data = (double **) malloc (mtree_size * sizeof(double *));
    
    for (int i = 0; i < mtree_size; i++) mtree_data[i] = (double *) malloc (num_fields_mtree * sizeof(double));

    //if (!omp_get_thread_num()) std::cerr << "Just reached problem region..." << std::endl;

    // let's also set up our output file
    std::ofstream outfile;

    //setOutputFileName();
    if (!outfile_name.size()) outfile_name = output;
    
    outfile.open(outfile_name.c_str(),std::ofstream::out);

    // next, we'll need to start reading in trees and storing them in memory

    int dot_freq = int(ceil(num_lines/100.));

    std::cout << "Number of halo trees to lookup: " << num_lines << std::endl;
    std::cout << "CHUNKSIZE: " << CHUNKSIZE << std::endl;

    // iterator used to store lookups to halo_id map
    std::unordered_map<long int, long int>::const_iterator fetch;

    // keep track halo id index in locations data
    int i = 0;

    // keep track of lines in output file data
    long long of_count = 0;

    // keep track of total number of lines stored
    long long global_of_count = 0;

    // loop over z = 0 halos in input catalog
    for (int j = 0; j < num_lines; j++) {

        if (j % dot_freq == 0) {
            std::cout << ceil(j/dot_freq) << ".";
            std::cout.flush();
        }

        // find halo id in locations data
        fetch = halo_id.find(data[1][j]);

        if (fetch == halo_id.end()) {
            std::cerr << "Halo " << data[1][j] << " not found." << std::endl;
            continue;
        }

        else i = fetch->second;

        // extract file name digits from file_id
        int k = (int)loc_data[1][i];
        std::stringstream s;
        s << base << "tree_" << k/25 << "_" << (k%25)/5 << "_" << (k%25)%5 << ".dat";

        // open corresponding merger tree file and find start of specified tree
        mtree.open(s.str(), std::ifstream::in);
        mtree.seekg(loc_data[2][i]);

        // lets read in the whole tree and output the
        // most massive progenitor branch
        long int curr_mmp = -1;     // id of current most massive progenitor
        long int desc = data[1][j]; // id of descendant halo
        double max_mvir = 0;        // maximum mass of halos in redshift interval
        double curr_z = -1;         // global redshift interval
        int z = 0, id = 1, desc_id = 3, mvir = 10;
        double curr_mmp_data[num_fields_mtree]; // holder for mmp data until we permanently store it
        double tmp_data[num_fields_mtree];      // holder for data in mtree file as it is read in

        // keep track of line number in merger tree
        k = -1;

        // scan through merger tree until we reach end of MMPB or next tree starts
        while (true) {

            // we may need to reallocate mtree_data if we encounter a tree that is too large
            if (of_count >= mtree_size) {

                std::cerr << "treeWalk chunk trying to read in more lines than there is space for in mtree_data!" << std::endl;
                std::cerr << "of_count: " << of_count << ", mtree_size: " << mtree_size << std::endl;
            }

            // increment merger tree line counter
            k++;

            // read in one line of data
            getline(mtree,tmp_str);

            // check if we've reached the end of the merger tree
            if (tmp_str[0] == '#') break;

            // parse each field
            std::istringstream line_stream(tmp_str);
            tmp_str = "";
            int l = 0;

            while (line_stream >> tmp_str) {
                tmp_data[l++] = std::strtod(tmp_str.c_str(),NULL);
            }
            
            // only output most massive progenitor branch?
            if (mmp) {

                // always store output from final timestep
                if (!k) {

                    for (l = 0; l < num_fields_mtree; l++) mtree_data[of_count][l] = tmp_data[l];
                    of_count++;

                    continue;
                }

                // second line in tree. set global redshift.
                else if (k == 1) curr_z = tmp_data[z];

                // when we reach a new timestep, output mmp from previous timestep
                if (curr_z != tmp_data[z]) {

                    // if mmp never gets set after first redshift
                    if (curr_mmp == -1) {
                        std::cout << "No MMP found for z = " << curr_z << std::endl;
                        break;
                    }

                    // we haven't found a new mmp, so this must be the end of the MMPB.
                    if (desc == curr_mmp) {
                        break;
                    }
                    
                    // update and output mmp from prev timestep.
                    // descendant is now mmp from previous z.
                    desc = curr_mmp;
                    max_mvir = 0;

                    // store mmp data
                    for (l = 0; l < num_fields_mtree; l++) mtree_data[of_count][l] = curr_mmp_data[l];

                    of_count++;

                    curr_z = tmp_data[z];
                    
                }

                // determine if this halo is mmp.
                // make sure it has correct descendant.
                if (tmp_data[desc_id] == desc) {

                    // only choose most massive
                    if (tmp_data[mvir] > max_mvir) {

                        // nice.. found most massive progenitor so far in this z
                        curr_mmp = (long int)tmp_data[id];
                        max_mvir = tmp_data[mvir];

                        for (l = 0; l < num_fields_mtree; l++) curr_mmp_data[l] = tmp_data[l];
                    }
                }
            }

            // just output all halos in merger tree
            else {
                for (l = 0; l < num_fields_mtree; l++) mtree_data[of_count][l] = tmp_data[l];
                of_count++;
            }
        }

        mtree.close();

        // once we fill up our data structure, let's write it to disk.
        // we should keep a buffer of ~500 lines so we don't go over.
        // if the buffer was too small, we would continue to the next iteration
        // and overfill our data structure.

        int bufferamount;
        if (mmp) bufferamount = 500;
        else bufferamount = 800000;

        // also do this if we are on last iteration (done collecting data)
        // or if whole merger trees are being output rather than just mmp branches
        if (of_count > (CHUNKSIZE - bufferamount) || (j == num_lines - 1)) {

            // finish writing progress dots 'complete' for j loop
            if (j == num_lines - 1) std::cout << " complete.";

            std::cout << std::endl << "Writing chunk to output file \"" << outfile_name << "\"" << std::endl;
            std::cout << "num_lines: " << of_count << std::endl;
            std::cout << "num_fields: " << num_fields_mtree << std::endl;

            // update total line count
            global_of_count += of_count;

            long int header[2] = {global_of_count, num_fields_mtree};

            // save current file position. (-1 for unset).
            long int pos = -1;
            
            // if file hasn't been written to yet, don't save position.  if we did,
            // we would overwrite our data.
            if (outfile.tellp() != outfile.beg) pos = outfile.tellp();

            // write the file header with updated info
            outfile.seekp(outfile.beg);
            outfile.write((const char *)header,2*sizeof(long int));

            // place back at end
            if (pos != -1) outfile.seekp(pos);

            int dot_freq_of = int(ceil(of_count/10.));

            // write data
            for (k = 0; k < of_count; k++) {
                if (k % dot_freq_of == 0) {
                    std::cout << ".";
                    std::cout.flush();
                }
                outfile.write((const char *)mtree_data[k],num_fields_mtree*sizeof(double));
            }

            std::cout <<" complete." << std::endl;

            of_count = 0;
        }
    }

    outfile.close();

    free (mtree_data, mtree_size);

    if (free_loc) free(loc_data,num_fields_loc);

    return 0;
}

//====================================================
// doProgenitorHistory:
// Produces median halo progenitor histories for the provided
// group of halo merger trees.
// data:    the merger trees of all desired halos, with each
//          complete merger tree listed before the next tree starts.  This
//          should be provided as an input file in standard binary format,
//          consisting of all the merger trees listed one after the other.
// num_lines:
//          number of lines in the entire stacked merger tree file
// num_fields:
//          number of fields each halo has at each snapshot

//====================================================
int doProgenitorHistory (double ** &data, int num_lines, int num_fields) {

    // let's first set up our data structures
    struct halo {
        halo * parent = NULL;   // pointer to parent halo at higher z
        halo * child = NULL;    // pointer to child halo at lower z
        bool alive = true;      // indicates whether tree has been dropped or not
        bool rez = true;        // indicates whether halo should be resurrected for next property
        double z = 0;           // current redshift
        double mvir = 0;        // virial mass (will be normalized by final mass)
        double cnfw = 0;        // nfw concentration
        double lambdap = 0;     // lambda prime: bullock spin parameter
        double dlambdap = 0;    // halo torque (dlambdap/dt)
        double highspin = 1;    // keeps track of if halo has had high spin throughout history
        int pop = -1;           // keeps track of which initial population halo belonged to (for spin mobility, for ex)
        double tf = 0;          // tidal force
        double rs = 0;          // scale radius (NFW) in comoving coordinates
        double rsp = 0;         // scale radius (NFW) in physical coordinates
        double rvir = 0;        // virial radius in comoving coordinates
        double rvirp = 0;       // virial radius in physical coordinates
        double smar = 0;        // specific mass accretion rate (normalized by final mass)
        double vmax = 0;        // maximum circular velocity
        double xoff = 0;        // offset between average position and peak density
        double tu = 0;          // virial ratio T/|U|
        double e = 0;           // elongation
        double e500 = 0;        // elonogation measured at r500
        double eratio = 0;      // ratio of virial to r500 elongation
        double lambda = 0;      // peebles spin parameter
        double cklypin = 0;     // concentration using klypin's scale radius
        double almm = 0;        // scale of last major merger
        double almmZ0 = 0;      // scale of last major merger of final halo
        double almms = 0;       // scale of last stellar major merger
        double almmZ0s = 0;     // scale of last stellar major merger of final halo
        double temp = 0;        // used as a placeholder to swap other properties in/out of
        double mlumpy = 0;      // total mass of dark matter accreted in lumps
        double msmooth = 0;     // total mass of dark matter accreted smoothly
        double merger_mass_avg = 0; // average mass of mergers at each timestep
        double merger_ratio_avg = 0;// average merger ratio at each timestep
        double merger_count = 0;   // merger count at each timestep
        double mdot_dyn = 0;    // dynamically time averaged mass accretion rate (requires treeAnalysis)
        double smdot_dyn = 0;   // specific mass accretion rate (t_dyn) (requires treeAnalysis)
    };

    struct mobility_rank {
        double p20[3] = {0,0,0};
        double p40[3] = {0,0,0};
        double p60[3] = {0,0,0};
        double p80[3] = {0,0,0};
    };

    struct halo_stats {
        halo median;
        halo disp20;
        halo disp80;
        halo ci_95_lo;
        halo ci_95_hi;
        mobility_rank spin_rank;
        mobility_rank tf_rank;
    };

    int num_trees = 1;          // number of merger trees in input file

    double a0 = data[0][0];     // initial scale factor in file. recurrences of this
                                // scale factor will tell us how mnay merger trees are present.

    int SMA_N = 20;             // number of timesteps (progenitors) to include in the SMA

    // let's determine how many merger trees were provided in the input file
    for (int i = 1; i < num_lines; i++) {
        if (data[i][0] == a0) num_trees++;
    }

    std::cout << "Processing " << num_trees << " merger trees..." << std::endl;

    // ADJUSTABLE PARAMETERS =========================================================

    bool output_z_list = false;     // output a list of redshifts we need time data for

    int tree_limit = num_trees;     // process only this many merger trees from file.
                                    // could be set to any number <= num_trees

    double mvir_lower_limit = 10.0; // lower threshold on mass before tree gets dropped

    int halo_count_lower_limit = 5; // once alive halo population drops below this number
                                    // we will stop recording statistics.

    std::string tin_path = "t_list.dat";    // file path for time data corresponding to
                                            // each redshift in the merger trees.
    std::string zout_path = "z_list.dat";   // file path to write output file containing
                                            // all redshifts in the merger trees we need
                                            // time data for.

    //================================================================================

    std::vector<std::vector<halo>> z;

    int j = 0;                  // z vector index.
    halo * curr_child = NULL;   // pointer to current child.
    bool alive = true;          // keeps track of whether tree has be declared 
                                // dead or not for ancestors.

    int tree_count = 0;

    std::cout << "Populating data structure with merger trees...";
    std::cout.flush();

    // it will be necessary to reindex the data to row major order in the future to improve
    // performance.
    // now let's populate our data structure with all the halos and their progenitors
    for (int i = 0; i < num_lines; i++) {

        // check if we've reached the start of a new merger tree. (a > 1.0)
        // don't do this if we are on the first tree.
        if ((data[i][0] > 1.0) && i) {

            // if previous tree is dead, continue with new one
            if (!alive) {
                curr_child = NULL;
                j = 0;
                alive = true;
                if (++tree_count >= tree_limit) break;
            }

            //if not, we need to add dead parent in there
            else {
                //std::cout << "Encountered alive tree when next tree started. Inserting dead parent." << std::endl;
                alive = false;
                i--;
            }
        }

        halo h;
        h.z = 1./data[i][BP_TREE.a] - 1.;
        h.mvir = data[i][BP_TREE.mvir];
        h.rvir = data[i][BP_TREE.rvir];
        h.rvirp = data[i][BP_TREE.rvir]*data[i][BP_TREE.a];
        h.rs = data[i][BP_TREE.rs];
        h.rsp = data[i][BP_TREE.rs]*data[i][BP_TREE.a];
        h.cnfw = h.rvir / h.rs;
        h.cklypin = h.rvir / data[i][BP_TREE.rsk];
        h.vmax = data[i][BP_TREE.vmax];
        h.tf = data[i][BP_TREE.tf];
        h.xoff = data[i][BP_TREE.xoff]/h.rvir;
        h.tu = data[i][BP_TREE.tu];
        h.lambdap = data[i][BP_TREE.spinb];
        h.lambda = data[i][BP_TREE.spinp];
        h.e = 1-sqrt(data[i][BP_TREE.bta]*data[i][BP_TREE.bta]+data[i][BP_TREE.cta]*data[i][BP_TREE.cta])/sqrt(2.);
        h.e500 = 1-sqrt(data[i][BP_TREE.bta500]*data[i][BP_TREE.bta500]+data[i][BP_TREE.cta500]*data[i][BP_TREE.cta500])/sqrt(2.);
        h.eratio = h.e/h.e500;
        h.almm = data[i][BP_TREE.almm];

        // the below quantities rely upon having run treeAnalysis on the merger trees first. check if this is the case before using them.
        if (num_fields-1 >= BP_TREE.merger_mass_avg) {
            h.almms = data[i][BP_TREE.almms];
            h.mlumpy = data[i][BP_TREE.mlumpy]/data[i][BP_TREE.mpeak]; // normalize by mpeak
            h.msmooth = data[i][BP_TREE.msmooth]/data[i][BP_TREE.mpeak];
            h.merger_count = data[i][BP_TREE.merger_count];
            h.merger_mass_avg = data[i][BP_TREE.merger_mass_avg];
            h.merger_ratio_avg = h.merger_mass_avg/h.mvir;
            h.mdot_dyn = data[i][BP_TREE.mdot_dyn];
            h.smdot_dyn = h.mdot_dyn/h.mvir;
        }

        // once a tree has died, make sure all ancestors also die
        if (!alive) {
            h.alive = false;
            h.rez = false;
        }

        // drop trees with masses below a threshold mass
        if (log10(h.mvir) < mvir_lower_limit) {
            h.alive = false;
            h.rez = false;
            alive = false;
        }

        // update pointer from this "parent" to child (previous halo in tree)
        if (curr_child) {
            h.child = curr_child;
        }

        // check if this redshift exists already. if not, initialize a vector.
        // reserve is necessary so that we don't get undesireable reallocation
        // of the vectors as they grow, invalidating all of our pointers.
        if (j >= z.size()) {
            std::vector<halo> new_z;
            z.push_back(new_z);
            z.back().reserve(num_trees);
        }

        z[j].push_back(h);

        // update pointer from child (previous halo in tree) to this "parent"
        if (curr_child) {
            curr_child->parent = &(z[j].back());
        }

        // use this to update pointers to/from this "child" to parent (next halo in tree)
        curr_child = &(z[j].back());

        // increment redshift counter
        j++;

        // check if we are on the last line in the file
        if (i+1 == num_lines) {

            // if we have a live tree, we need to insert dead parent
            // before terminating
            if (alive) {
                alive = false;
                i--;
            }
        }

    }

    int z_counter = 0;  // counts the number of redshifts we have alive halos at

    // let's output a list of redshifts
    std::ofstream zout;
    if (output_z_list) zout.open(zout_path,std::ofstream::out);

    for (j = 0; j < z.size(); j++) {
        for (int i = 0; i < z[j].size(); i++) {
            if (z[j][i].alive) {
                if (output_z_list) zout << z[j][i].z << std::endl;
                z_counter++;
                break;
            }
        }
    }

    if (output_z_list) {

        zout.close();

        std::cout << "Z list output to file: \"" << zout_path << "\"" << std::endl;

        // return since we need to create t_list file before continuing
        return 0;
    }

    std::cout << "Expecting T list file to be found at: \"" << tin_path << "\"" << std::endl;
    std::cout.flush();

    double t[z_counter];

    // the following was used to create t_list file:
    // $ rm -f t_list.dat; cat z_list.dat | awk '{print "python cosmocalc.py "$1" 67.8 0.307 0.693"}' | /bin/bash | awk '{print $1}' >> t_list.dat

    // now, let's read in and store our time data (to correspond to each redshift
    // in the merger trees).
    // note, we are assuming the file containing time data can be found at the path below.
    std::ifstream tin (tin_path,std::ifstream::in);

    for (int i = 0; i < z_counter; i++) {
        tin >> t[i];
    }

    tin.close();

    std::cout << "Updating time relevant halo data...";
    std::cout.flush();

    // now we'll go through and update all halo properties that need time information
    for (j = 0; j < z.size()-1; j++) {
        for (int i = 0; i < z[j].size(); i++) {
            halo * h1 = &(z[j][i]);
            halo * h2 = h1->parent;

            // only compute this for halos that are still alive.
            if (h1->alive) {

                // if the two z's are the same it is because of an
                // artificial (dead) ancestor -- which won't have proper data. might as
                // well kill the child to avoid any issues.
                if (h1->z == h2->z) {
                    h1->alive = false;
                }

                else {
                    // compute relevant quantities here
                    h1->smar = (h1->mvir - h2->mvir)/(t[j]-t[j+1]);
                    h1->dlambdap = (h1->lambdap - h2->lambdap)/(t[j]-t[j+1]);

                    // set final halos last major merger scale for all descendants
                    // so that we can easily refer to it in the future to check if
                    // that progenitor is past the last major merger scale or not
                    if (j == 0) {
                        // for halo mass
                        h1->almmZ0 = h1->almm;
                        double almmZ0 = h1->almmZ0;

                        // for stellar mass
                        h1->almmZ0s = h1->almms;
                        double almmZ0s = h1->almmZ0s;

                        halo * p = h2;
                        while (p->alive) {
                            p->almmZ0 = almmZ0;
                            p->almmZ0s = almmZ0s;
                            p = p->parent;
                        }
                    }
                }
            }
        }
    }

    std::cout << " done." << std::endl;

    std::cout << "Normalizing relevant quantities...";
    std::cout.flush();

    // now let's go through and normalize necessary quantities
    for (int i = 0; i < z[0].size(); i++) {
        double mvir0 = z[0][i].mvir;

        // normalize final snapshot
        //z[0][i].mvir /= mvir0;
        z[0][i].smar /= mvir0;

        // normalize all alive ancestors
        halo * p = z[0][i].parent;
        while(p->alive) {
            //p->mvir /= mvir0;
            p->smar /= mvir0;
            p = p->parent;
        }
    }

    std::cout << " done." << std::endl;

    std::cout << "Handling high spin determination...";
    std::cout.flush();

    // parameters for high spin function
    double lambda_min = 0.035;
    double lambda_cut = 0.035;
    double zt = 2.0;

    // function returning lambdap cut at different redshifts.
    // if a halo has lambdap higher than the value returned by this function,
    // it would be considered a high spin halo.
    auto highspin_f = [&lambda_min,&lambda_cut,&zt](double z)->double {
        if (z > zt) return lambda_cut;
        else return (lambda_min+(z/zt)*(lambda_cut-lambda_min));
    };
    //----------------------------------

    //bool output_track = true;
    bool output_track = false;

    std::ofstream output_track_1;
    std::ofstream output_track_2;

    // Handling high spin determination
    for (int i = 0; i < z[0].size(); i++) {

        // shouldn't happen
        if (!z[0][i].alive) continue;

        if (output_track) {
            output_track_1.open(output+".track1_"+str(i),std::ofstream::out);
            output_track_2.open(output+".track2_"+str(i),std::ofstream::out);
        }

        // this is the halo at the back of the SMA group (highest z, most recently added)
        halo * h1 = &(z[0][i]);

        // this is the lowest z halo at the front of the SMA group that
        //will be popped to "move" the average
        halo * SMA_pop = h1;

        // this is the halo at the center of the SMA group, used to determine
        // appropriate scale for SMA
        halo * SMA_center = h1;

        // simple moving average
        double SMA = 0;

        // counter to let us know for how many iterations so far we've been at the
        // end of the merge tree
        int end_of_tree = 0;

        // iteration counter
        int j = 0;

        while (true) {

            // if less than SMA_N halos since starting, scale appropriately
            if (j < SMA_N) {

                // this could happen if SMA_N is large enough or the tree is short enough, but
                // its just a weird edge case we should have to deal with in general.
                if (end_of_tree) {
                    //std::cout << "This is a problem.. we've reached end of tree before having SMA_N = ";
                    //std::cout << SMA_N << " halos in the moving average." << std::endl;
                    break;
                }

                SMA *= j/(j+1.);
                SMA += h1->lambdap / (j+1);
            }

            // if we've reached end of tree, scale SMA accordingly as halos pop off front
            else if (end_of_tree) {
                SMA -= SMA_pop->lambdap;
                SMA *= (SMA_N+1.-end_of_tree)/(SMA_N-end_of_tree);
                SMA_pop = SMA_pop->parent;
            }

            // this is the normal case, where SMA is moving along, one halo pushed on back, one popped
            // off front
            else {
                SMA += (h1->lambdap - SMA_pop->lambdap) / SMA_N;
                SMA_pop = SMA_pop->parent;
            }

            // update center halo if necessary
            if (j > (int)floor(SMA_N/2.)) {

                // we only should check the SMA condition when we are actually going to update the
                // SMA_center.  otherwise, we may not be checking the actual SMA values that we
                // want (for example at the beginning, when the SMA group is small and we are still
                // adding halos to the back.
                if (SMA < highspin_f(SMA_center->z)) {

                    SMA_center->highspin = 0;

                    // falsify ancestors as well
                    halo * p = SMA_center->parent;
                    while(p->alive) {
                        p->highspin = 0;
                        p = p->parent;
                    }
                    break;
                }

                // output to individual halo tracks if desired
                if (output_track) {
                    output_track_1 << SMA_center->z << " " << SMA_center->lambdap << std::endl;
                    output_track_2 << SMA_center->z << " " << SMA << std::endl;
                }

                // update center halo
                SMA_center = SMA_center->parent;
            }

            // update h1 unless we've reached end of tree.  in that case,
            // h1 remains the same, but SMA_pop continues until all halos
            // are popped off
            if (h1->parent->alive) h1 = h1->parent;

            else {

                // end_of_tree increments each time we continue to compute SMA up to the last halo
                // in the tree.  So the value of end_of_tree is by how much the SMA group has been
                // reduced.
                end_of_tree++;

                // this would indicate SMA_center has reached the last alive halo in the tree, so
                // we can stop here.
                if (SMA_center == h1) {
                    break;
                }
            }

            j++;
        }

        // print remaining halos in output track if desired (if any are left). if we dropped below high spin
        // threshold, then there should be halos left.  if we made it to end of tree, there might be one halo
        // left (the last one in the tree)
        halo * p = SMA_center;
        int remaining = 0;
        while(p->alive) {

            remaining++;

            if (output_track) {
                output_track_1 << p->z << " " << p->lambdap << std::endl;
            }

            p = p->parent;
        }

        // only want to output a few tracks
        if (i>20) {
            output_track = false;
        }

        //else std::cout << "remaining: " << remaining << std::endl;

        if (output_track) {
            output_track_1.close();
            output_track_2.close();
        }
    }

    std::cout << " done." << std::endl;    

    // define comparator for sorting halos by alive/dead as well as temp halo property.
    // the ordering should reflect alive/small, alive/large, dead/small, dead/large
    // (where small/large refers to whichever property is stored in halo.temp).
    struct sort_comp {
        bool operator()(const halo &left, const halo &right) {
            return (left.alive == right.alive) ? left.temp < right.temp : left.alive > right.alive ;
        }
    };

    int highspincount = 0; // number of high spin halos at a given redshift
    int highspincountZ0 = 0; // number of high spin halos at z = 0
    int majormergercount = 0; // number of halos in high spin group that have had a major merger since z = 0

    // keeps track of what halo property we are following 
    int field_counter = 0;

    // number of halo progenitor fields we're going to track
    int num_fields_prog = 28;

    // vector containing halo statistics for each redshift
    std::vector<halo_stats> stats;

    std::cout << "Processing halo property number ";
    std::cout.flush();

    while (true) {

        // now place temp property into appropriate halo field.
        // do this only if we've already iterated through one property.
        // note: using field_counter-1 b/c it is next interation, but we
        // still need to clean up from previous iteration.
        if (field_counter) {
            for (j = 0; j < stats.size(); j++) {
                switch (field_counter-1) {
                    case 0: stats[j].median.mvir                = stats[j].median.temp;
                            stats[j].disp20.mvir                = stats[j].disp20.temp;
                            stats[j].disp80.mvir                = stats[j].disp80.temp;
                            stats[j].ci_95_lo.mvir              = stats[j].ci_95_lo.temp;
                            stats[j].ci_95_hi.mvir              = stats[j].ci_95_hi.temp; break;
                    case 1: stats[j].median.cnfw                = stats[j].median.temp;
                            stats[j].disp20.cnfw                = stats[j].disp20.temp;
                            stats[j].disp80.cnfw                = stats[j].disp80.temp;
                            stats[j].ci_95_lo.cnfw              = stats[j].ci_95_lo.temp;
                            stats[j].ci_95_hi.cnfw              = stats[j].ci_95_hi.temp; break;
                    case 2: stats[j].median.rvir                = stats[j].median.temp;
                            stats[j].disp20.rvir                = stats[j].disp20.temp;
                            stats[j].disp80.rvir                = stats[j].disp80.temp;
                            stats[j].ci_95_lo.rvir              = stats[j].ci_95_lo.temp;
                            stats[j].ci_95_hi.rvir              = stats[j].ci_95_hi.temp; break;
                    case 3: stats[j].median.rs                  = stats[j].median.temp;
                            stats[j].disp20.rs                  = stats[j].disp20.temp;
                            stats[j].disp80.rs                  = stats[j].disp80.temp;
                            stats[j].ci_95_lo.rs                = stats[j].ci_95_lo.temp;
                            stats[j].ci_95_hi.rs                = stats[j].ci_95_hi.temp; break;
                    case 4: stats[j].median.smar                = stats[j].median.temp;
                            stats[j].disp20.smar                = stats[j].disp20.temp;
                            stats[j].disp80.smar                = stats[j].disp80.temp;
                            stats[j].ci_95_lo.smar              = stats[j].ci_95_lo.temp;
                            stats[j].ci_95_hi.smar              = stats[j].ci_95_hi.temp; break;
                    case 5: stats[j].median.lambdap             = stats[j].median.temp;
                            stats[j].disp20.lambdap             = stats[j].disp20.temp;
                            stats[j].disp80.lambdap             = stats[j].disp80.temp;
                            stats[j].ci_95_lo.lambdap           = stats[j].ci_95_lo.temp;
                            stats[j].ci_95_hi.lambdap           = stats[j].ci_95_hi.temp; break;
                    case 6: stats[j].median.dlambdap            = stats[j].median.temp;
                            stats[j].disp20.dlambdap            = stats[j].disp20.temp;
                            stats[j].disp80.dlambdap            = stats[j].disp80.temp;
                            stats[j].ci_95_lo.dlambdap          = stats[j].ci_95_lo.temp;
                            stats[j].ci_95_hi.dlambdap          = stats[j].ci_95_hi.temp; break;
                    case 7: stats[j].median.tf                  = stats[j].median.temp;
                            stats[j].disp20.tf                  = stats[j].disp20.temp;
                            stats[j].disp80.tf                  = stats[j].disp80.temp;
                            stats[j].ci_95_lo.tf                = stats[j].ci_95_lo.temp;
                            stats[j].ci_95_hi.tf                = stats[j].ci_95_hi.temp; break;
                    case 8: stats[j].median.vmax                = stats[j].median.temp;
                            stats[j].disp20.vmax                = stats[j].disp20.temp;
                            stats[j].disp80.vmax                = stats[j].disp80.temp;
                            stats[j].ci_95_lo.vmax              = stats[j].ci_95_lo.temp;
                            stats[j].ci_95_hi.vmax              = stats[j].ci_95_hi.temp; break;
                    case 9: stats[j].median.xoff                = stats[j].median.temp;
                            stats[j].disp20.xoff                = stats[j].disp20.temp;
                            stats[j].disp80.xoff                = stats[j].disp80.temp;
                            stats[j].ci_95_lo.xoff              = stats[j].ci_95_lo.temp;
                            stats[j].ci_95_hi.xoff              = stats[j].ci_95_hi.temp; break;
                    case 10:stats[j].median.e                   = stats[j].median.temp;
                            stats[j].disp20.e                   = stats[j].disp20.temp;
                            stats[j].disp80.e                   = stats[j].disp80.temp;
                            stats[j].ci_95_lo.e                 = stats[j].ci_95_lo.temp;
                            stats[j].ci_95_hi.e                 = stats[j].ci_95_hi.temp; break;
                    case 11:stats[j].median.rsp                 = stats[j].median.temp;
                            stats[j].disp20.rsp                 = stats[j].disp20.temp;
                            stats[j].disp80.rsp                 = stats[j].disp80.temp;
                            stats[j].ci_95_lo.rsp               = stats[j].ci_95_lo.temp;
                            stats[j].ci_95_hi.rsp               = stats[j].ci_95_hi.temp; break;
                    case 12:stats[j].median.rvirp               = stats[j].median.temp;
                            stats[j].disp20.rvirp               = stats[j].disp20.temp;
                            stats[j].disp80.rvirp               = stats[j].disp80.temp;
                            stats[j].ci_95_lo.rvirp             = stats[j].ci_95_lo.temp;
                            stats[j].ci_95_hi.rvirp             = stats[j].ci_95_hi.temp; break;
                    case 13:stats[j].median.e500                = stats[j].median.temp;
                            stats[j].disp20.e500                = stats[j].disp20.temp;
                            stats[j].disp80.e500                = stats[j].disp80.temp;
                            stats[j].ci_95_lo.e500              = stats[j].ci_95_lo.temp;
                            stats[j].ci_95_hi.e500              = stats[j].ci_95_hi.temp; break;
                    case 14:stats[j].median.lambda              = stats[j].median.temp;
                            stats[j].disp20.lambda              = stats[j].disp20.temp;
                            stats[j].disp80.lambda              = stats[j].disp80.temp;
                            stats[j].ci_95_lo.lambda            = stats[j].ci_95_lo.temp;
                            stats[j].ci_95_hi.lambda            = stats[j].ci_95_hi.temp; break;
                    case 15:stats[j].median.cklypin             = stats[j].median.temp;
                            stats[j].disp20.cklypin             = stats[j].disp20.temp;
                            stats[j].disp80.cklypin             = stats[j].disp80.temp;
                            stats[j].ci_95_lo.cklypin           = stats[j].ci_95_lo.temp;
                            stats[j].ci_95_hi.cklypin           = stats[j].ci_95_hi.temp; break;
                    case 16:stats[j].median.highspin            = stats[j].median.temp;
                            stats[j].disp20.highspin            = stats[j].disp20.temp;
                            stats[j].disp80.highspin            = stats[j].disp80.temp;
                            stats[j].ci_95_lo.highspin          = stats[j].ci_95_lo.temp;
                            stats[j].ci_95_hi.highspin          = stats[j].ci_95_hi.temp; break;
                    case 17:    // spin mobility rank
                    case 18:    // tf mobility rank
                                                                                          break;
                    case 19:stats[j].median.mlumpy              = stats[j].median.temp;
                            stats[j].disp20.mlumpy              = stats[j].disp20.temp;
                            stats[j].disp80.mlumpy              = stats[j].disp80.temp;
                            stats[j].ci_95_lo.mlumpy            = stats[j].ci_95_lo.temp;
                            stats[j].ci_95_hi.mlumpy            = stats[j].ci_95_hi.temp; break;
                    case 20:stats[j].median.msmooth             = stats[j].median.temp;
                            stats[j].disp20.msmooth             = stats[j].disp20.temp;
                            stats[j].disp80.msmooth             = stats[j].disp80.temp;
                            stats[j].ci_95_lo.msmooth           = stats[j].ci_95_lo.temp;
                            stats[j].ci_95_hi.msmooth           = stats[j].ci_95_hi.temp; break;
                    case 21:stats[j].median.merger_count        = stats[j].median.temp;
                            stats[j].disp20.merger_count        = stats[j].disp20.temp;
                            stats[j].disp80.merger_count        = stats[j].disp80.temp;
                            stats[j].ci_95_lo.merger_count      = stats[j].ci_95_lo.temp;
                            stats[j].ci_95_hi.merger_count      = stats[j].ci_95_hi.temp; break;
                    case 22:stats[j].median.merger_mass_avg     = stats[j].median.temp;
                            stats[j].disp20.merger_mass_avg     = stats[j].disp20.temp;
                            stats[j].disp80.merger_mass_avg     = stats[j].disp80.temp;
                            stats[j].ci_95_lo.merger_mass_avg   = stats[j].ci_95_lo.temp;
                            stats[j].ci_95_hi.merger_mass_avg   = stats[j].ci_95_hi.temp; break;
                    case 23:stats[j].median.merger_ratio_avg    = stats[j].median.temp;
                            stats[j].disp20.merger_ratio_avg    = stats[j].disp20.temp;
                            stats[j].disp80.merger_ratio_avg    = stats[j].disp80.temp;
                            stats[j].ci_95_lo.merger_ratio_avg  = stats[j].ci_95_lo.temp;
                            stats[j].ci_95_hi.merger_ratio_avg  = stats[j].ci_95_hi.temp; break;
                    case 24:stats[j].median.eratio              = stats[j].median.temp;
                            stats[j].disp20.eratio              = stats[j].disp20.temp;
                            stats[j].disp80.eratio              = stats[j].disp80.temp;
                            stats[j].ci_95_lo.eratio            = stats[j].ci_95_lo.temp;
                            stats[j].ci_95_hi.eratio            = stats[j].ci_95_hi.temp; break;
                    case 25:stats[j].median.tu                  = stats[j].median.temp;
                            stats[j].disp20.tu                  = stats[j].disp20.temp;
                            stats[j].disp80.tu                  = stats[j].disp80.temp;
                            stats[j].ci_95_lo.tu                = stats[j].ci_95_lo.temp;
                            stats[j].ci_95_hi.tu                = stats[j].ci_95_hi.temp; break;
                    case 26:stats[j].median.mdot_dyn            = stats[j].median.temp;
                            stats[j].disp20.mdot_dyn            = stats[j].disp20.temp;
                            stats[j].disp80.mdot_dyn            = stats[j].disp80.temp;
                            stats[j].ci_95_lo.mdot_dyn          = stats[j].ci_95_lo.temp;
                            stats[j].ci_95_hi.mdot_dyn          = stats[j].ci_95_hi.temp; break;
                    case 27:stats[j].median.smdot_dyn           = stats[j].median.temp;
                            stats[j].disp20.smdot_dyn           = stats[j].disp20.temp;
                            stats[j].disp80.smdot_dyn           = stats[j].disp80.temp;
                            stats[j].ci_95_lo.smdot_dyn         = stats[j].ci_95_lo.temp;
                            stats[j].ci_95_hi.smdot_dyn         = stats[j].ci_95_hi.temp; break;

                }

                // re-initialize so it doesn't contaminate future analysis
                stats[j].median.temp    = 0;
                stats[j].disp20.temp    = 0;
                stats[j].disp80.temp    = 0;
                stats[j].ci_95_lo.temp  = 0;
                stats[j].ci_95_hi.temp  = 0;
            }
        }

        // let's finish this loop if we've analyzed all desired properties
        if (field_counter >= num_fields_prog) break;

        // let's choose what halo property we want to track.
        for (j = 0; j < z.size(); j++) {
            for (int i = 0; i < z[j].size(); i++) {

                // we are copying the desired property to the temp field
                // so that we can use the same sorting comparator for all
                // fields.
                switch (field_counter) {
                    case 0: z[j][i].temp = z[j][i].mvir;            break;
                    case 1: z[j][i].temp = z[j][i].cnfw;            break;
                    case 2: z[j][i].temp = z[j][i].rvir;            break;
                    case 3: z[j][i].temp = z[j][i].rs;              break;
                    case 4: z[j][i].temp = z[j][i].smar;            break;
                    case 5: z[j][i].temp = z[j][i].lambdap;         break;
                    case 6: z[j][i].temp = z[j][i].dlambdap;        break;
                    case 7: z[j][i].temp = z[j][i].tf;              break;
                    case 8: z[j][i].temp = z[j][i].vmax;            break;
                    case 9: z[j][i].temp = z[j][i].xoff;            break;
                    case 10:z[j][i].temp = z[j][i].e;               break;
                    case 11:z[j][i].temp = z[j][i].rsp;             break;
                    case 12:z[j][i].temp = z[j][i].rvirp;           break;
                    case 13:z[j][i].temp = z[j][i].e500;            break;
                    case 14:z[j][i].temp = z[j][i].lambda;          break;
                    case 15:z[j][i].temp = z[j][i].cklypin;         break;
                    case 16:z[j][i].temp = z[j][i].lambdap;         break; // highspin
                    case 17:z[j][i].temp = z[j][i].lambdap;         break; // spin mobility
                    case 18:z[j][i].temp = z[j][i].tf;              break; // tf mobility
                    case 19:z[j][i].temp = z[j][i].mlumpy;          break;
                    case 20:z[j][i].temp = z[j][i].msmooth;         break;
                    case 21:z[j][i].temp = z[j][i].merger_count;    break;
                    case 22:z[j][i].temp = z[j][i].merger_mass_avg; break;
                    case 23:z[j][i].temp = z[j][i].merger_ratio_avg;break;
                    case 24:z[j][i].temp = z[j][i].eratio;          break;
                    case 25:z[j][i].temp = z[j][i].tu;              break;
                    case 26:z[j][i].temp = z[j][i].mdot_dyn;        break;
                    case 27:z[j][i].temp = z[j][i].smdot_dyn;       break;
                }
            }
        }

        std::cout << field_counter << ".";
        std::cout.flush();

        // increment field counter for next iteration
        field_counter++;

        // now that our data structure is populated with all the merger trees, let's
        // iterate through each redshift and sort / find desired median properties
        for (j = 0; j < z.size(); j++) {

            // sort redshift vector by alive/dead and by halo property. this will group all
            // alive trees together and all dead trees together. each alive/dead group
            // will be further sorted by the specified halo property
            std::sort(z[j].begin(),z[j].end(),sort_comp());

            int sep = 0;    // alive/dead separation index

            // find where the alive/dead group divide is in the vector
            for (int i = z[j].size()-1; i >= 0; i--) {
                if (z[j][i].alive) {
                    sep = i+1;
                    break;
                }
            }

            // determine high spin divide percentile
            if (field_counter == 17) {
                highspincount = 0;
                majormergercount = 0;
                for (int i = 0; i < sep; i++) {
                    if (z[j][i].highspin == 1) {
                        highspincount++;

                        // this assumes we are using final_ms 1:10 type comparison, where
                        // a -1 would indicate that halo never had a stellar mm with significant
                        // mass of final galaxy
                        if (z[j][i].almmZ0s != -1) majormergercount++;
                        
                        //double zlmm = (1./z[j][i].almmZ0-1.);
                        //double zlmm = (1./z[j][i].almmZ0s-1.);
                        //if (z[j][i].z > zlmm) {
                        //    majormergercount++;
                        //}
                    }
                }
                // use this to normalize other redshifts to be relative to z0 fraction
                if (j == 0) highspincountZ0 = highspincount;
            }

            // set rank mobility populations
            if (field_counter == 18 || field_counter == 19) {

                // let's determine which populations these halos belong to
                if (j == 0) {
                    for (int i = 0; i < sep; i++) {

                        z[j][i].pop = -1; // reset since could have changed in previous iteration

                        // pop 0 is < 20th percentile at z = 0
                        if (i < (int)(0.2*sep)) {
                            z[j][i].pop = 0;
                        }
                        // pop 1 is 40-60th percentile at z = 0
                        else if (i > (int)(0.4*sep) && i < (int)(0.6*sep)) {
                            z[j][i].pop = 1;
                        }
                        // pop 2 is > 80th percentile at z = 0
                        else if (i > (int)(0.8*sep)) {
                            z[j][i].pop = 2;
                        }

                        // propagate z = 0 pops to all progenitors
                        halo * p = z[j][i].parent;
                        while (p->alive) {
                            p->pop = z[j][i].pop;
                            p = p->parent;
                        }
                    }
                }
            }


            //std::cout << "Alive/Dead Separation at index: " << sep << "/" << z[j].size() << ", z = " << z[j][0].z << std::endl;

            // now, let's go through and find trees with dead parents
            // and drop corresponding trees with opposite ranks
            for (int i = 0; i < sep; i++) { 

                if (!z[j][i].parent->alive) {

                    //std::cout << "Found dead parent. Opp ranked tree at index " << sep-1-i << std::endl;

                    // find tree with opposing rank
                    halo * drop = &(z[j][sep-1-i]);
                    
                    // mark all ancestors as dead.
                    // if ancestors are already dead, we can stop.
                    while (drop->parent) {
                        if (drop->parent->alive) {
                            drop->parent->alive = false;
                            drop = drop->parent;
                        }
                        else break;
                    }
                }
            }

            // now we need to update parent/child relationships
            for (int i = 0; i < z[j].size(); i++) {
                halo * curr = &(z[j][i]);
                if (z[j][i].parent) z[j][i].parent->child = curr;
                if (z[j][i].child) z[j][i].child->parent = curr;
            }

            // finally, let's store the halo statistics (only if alive halos remain at
            // this redshift). let's limit this to the case where we have at least as many
            // halos as specified by our lower limit parameter.
            if (sep >= halo_count_lower_limit) {
                if (j >= stats.size()) {
                    halo_stats stat;
                    stats.push_back(stat);
                }

                // handle high spin fraction separately from normal properties
                if (field_counter == 17) { // high spin
                    //stats[j].median.temp = (sep-highspincount)/(double)sep;

                    // fraction of halos w/ high spin at z = 0 that still have high spin at this redshift
                    stats[j].median.temp = highspincount/(double)highspincountZ0; 
                    stats[j].disp20.temp = highspincount/(double)sep;
                    stats[j].disp80.temp = majormergercount/(double)highspincount;
                }

                else if (field_counter == 18 || field_counter == 19) {

                    int pop_count[3] = {0,0,0};

                    // now let's determine where the quintile cuts are at this redshift.
                    // first, how many of each pop are still alive?
                    for (int i = 0; i < sep; i++) {
                        switch (z[j][i].pop) {
                            case -1:                 break;
                            case  0: pop_count[0]++; break;
                            case  1: pop_count[1]++; break;
                            case  2: pop_count[2]++; break;
                        }
                    }

                    int pop_index[3] = {0,0,0};

                    // now go back through and record cuts
                    for (int i = 0; i < sep; i++) {
                        switch (z[j][i].pop) {
                            case -1:                 break;
                            case  0: pop_index[0]++; break;
                            case  1: pop_index[1]++; break;
                            case  2: pop_index[2]++; break;
                        }
                        
                        double percentile = 0;          // percentile w/i remaining population
                        double percentileplusone = 0;   // percentile of next ranking halo w/i pop
                        double total_rank = 0;          // percentile of halo wrt all alive halos

                        for (int k = 0; k < 3; k++) {
                            percentile = ((double)pop_index[k])/pop_count[k];
                            percentileplusone = ((double)(pop_index[k]+1))/pop_count[k];
                            total_rank = ((double)(i+1))/sep;

                            if (field_counter == 18) {

                                if (percentileplusone > 0.2 && percentile <= 0.2) {
                                    stats[j].spin_rank.p20[k] = z[j][i].temp; //total_rank;
                                }
                                else if (percentileplusone > 0.4 && percentile <= 0.4) {
                                    stats[j].spin_rank.p40[k] = z[j][i].temp; //total_rank;
                                }
                                else if (percentileplusone > 0.6 && percentile <= 0.6) {
                                    stats[j].spin_rank.p60[k] = z[j][i].temp; //total_rank;
                                }
                                else if (percentileplusone > 0.8 && percentile <= 0.8) {
                                    stats[j].spin_rank.p80[k] = z[j][i].temp; //total_rank;
                                }
                            }

                            else if (field_counter == 19) {

                                if (percentileplusone > 0.2 && percentile <= 0.2) {
                                    stats[j].tf_rank.p20[k] = z[j][i].temp; //total_rank;
                                }
                                else if (percentileplusone > 0.4 && percentile <= 0.4) {
                                    stats[j].tf_rank.p40[k] = z[j][i].temp; //total_rank;
                                }
                                else if (percentileplusone > 0.6 && percentile <= 0.6) {
                                    stats[j].tf_rank.p60[k] = z[j][i].temp; //total_rank;
                                }
                                else if (percentileplusone > 0.8 && percentile <= 0.8) {
                                    stats[j].tf_rank.p80[k] = z[j][i].temp; //total_rank;
                                }
                            }

                        }
                    }

                }

                // merger count, merger mass, merger ratio
                else if (field_counter == 22 || field_counter == 23 || field_counter == 24) {

                    // NOTE: while the below code uses median, disp20/80, ci_95_lo/hi, what is actually
                    // being computed in those variables is average, avg+-1sigma, avg+-std_mean

                    // take averages for these parameters since there are so many zero values
                    double sum = 0;
                    for (int i = 0; i < sep; i++) {
                        sum += z[j][i].temp;
                    }
                    stats[j].median.temp = sum/sep;

                    double stdv = 0;
                    for (int i = 0; i < sep; i++) {
                        stdv += pow(z[j][i].temp-stats[j].median.temp,2.);
                    }

                    stdv /= sep;
                    stdv = sqrt(stdv);

                    stats[j].disp20.temp = stats[j].median.temp-stdv;
                    stats[j].disp80.temp = stats[j].median.temp+stdv;

                    double stdv_mean = stats[j].median.temp/sqrt(sep);

                    stats[j].ci_95_lo.temp = stats[j].median.temp-stdv_mean;
                    stats[j].ci_95_hi.temp = stats[j].median.temp+stdv_mean;
                }

                // handle normal halo properties
                else {

                    //if (field_counter == 20 || field_counter == 21) {
                    //    std::cout << "Field counter: " << field_counter << ", 1+z = " << 1.+z[j][0].z << std::endl;
                    //    std::cout << "\tsep: " << sep << ", median sep: ";
                    //    if (sep % 2) std::cout << (sep-1)/2;
                    //    else std::cout << 0.5*((sep-2)/2 + sep/2);
                    //    std::cout << "\tdisp20: " << (int)floor(0.2*sep) << "\t disp80: " << (int)ceil(0.8*sep);
                    //    std::cout << "\n\tci_95_lo: " << (int)floor(sep*.5-1.96*sqrt(sep*.25)+1);
                    //    std::cout << "\tci_95_hi: " << (int)ceil(sep*.5+1.96*sqrt(sep*.25)-1.) << std::endl;
                    //}

                    // store medians
                    if (sep % 2) stats[j].median.temp = z[j][(sep-1)/2].temp;
                    else stats[j].median.temp = 0.5*(z[j][(sep-2)/2].temp + z[j][(sep)/2].temp);

                    // store z for medians only (no need to redundantly store it anywhere else)
                    stats[j].median.z = z[j][0].z;

                    // store disp20/80
                    stats[j].disp20.temp = z[j][(int)floor(0.2*sep)].temp;
                    stats[j].disp80.temp = z[j][(int)ceil(0.8*sep)].temp;

                    // store ci_95_lo/ci_95_hi
                    double lo = sep*.5-1.96*sqrt(sep*.25)+1;
                    double hi = sep*.5+1.96*sqrt(sep*.25)-1.;
                    stats[j].ci_95_lo.temp = z[j][(int)floor(lo)].temp;
                    stats[j].ci_95_hi.temp = z[j][(int)ceil(hi)].temp;
                }
            }
        }

        // let's restore the data structure so that we can run another halo property on it.
        // we can do this by just ressurecting halos marked for rezzing (those that are dead
        // due to rank dropping rather than low masses or tree death)
        for (j = 0; j < z.size(); j++) {
            for (int i = 0; i < z[j].size(); i++) {
                if (z[j][i].rez) z[j][i].alive = true;
            }
        }
    }

    std::cout << " done." << std::endl;

    //std::cout << "Medians:" << std::endl;
    //for (int j = 0; j < medians.size(); j++) {
    //    std::cout << "z = " << medians[j].z << ", M_vir = " << medians[j].mvir;
    //    std::cout << ", C_nfw = " << medians[j].cnfw << ", R_s = " << medians[j].rs;
    //    std::cout << ", R_vir = " << medians[j].rvir << ", MAR = " << medians[j].smar << std::endl;
    //}

    std::cout << "Writing output file...";

    // write the median data to output file
    std::ofstream of;
    setOutputFileName();
    of.open(output.c_str(),std::ofstream::out);

    for (j = 0; j < stats.size(); j++) {

        // column 1
        of << stats[j].median.z;

        // median  disp20  disp80  ci_95_lo  ci_95_hi
        // 2,3,4,5,6
        of << " " << stats[j].median.mvir       << " " << stats[j].disp20.mvir;
        of << " " << stats[j].disp80.mvir       << " " << stats[j].ci_95_lo.mvir;
        of << " " << stats[j].ci_95_hi.mvir;

        // 7,8,9,10,11
        of << " " << stats[j].median.smar       << " " << stats[j].disp20.smar;
        of << " " << stats[j].disp80.smar       << " " << stats[j].ci_95_lo.smar;
        of << " " << stats[j].ci_95_hi.smar;

        // 12,13,14,15,16
        of << " " << stats[j].median.cnfw       << " " << stats[j].disp20.cnfw;
        of << " " << stats[j].disp80.cnfw       << " " << stats[j].ci_95_lo.cnfw;
        of << " " << stats[j].ci_95_hi.cnfw;

        // 17,18,19,20,21
        of << " " << stats[j].median.rvir       << " " << stats[j].disp20.rvir;
        of << " " << stats[j].disp80.rvir       << " " << stats[j].ci_95_lo.rvir;
        of << " " << stats[j].ci_95_hi.rvir;

        // 22,23,24,25,26
        of << " " << stats[j].median.rs         << " " << stats[j].disp20.rs;
        of << " " << stats[j].disp80.rs         << " " << stats[j].ci_95_lo.rs;
        of << " " << stats[j].ci_95_hi.rs;

        // 27,28,29,30,31
        of << " " << stats[j].median.lambdap    << " " << stats[j].disp20.lambdap;
        of << " " << stats[j].disp80.lambdap    << " " << stats[j].ci_95_lo.lambdap;
        of << " " << stats[j].ci_95_hi.lambdap;
    
        // 32,33,34,35,36
        of << " " << stats[j].median.dlambdap   << " " << stats[j].disp20.dlambdap;
        of << " " << stats[j].disp80.dlambdap   << " " << stats[j].ci_95_lo.dlambdap;
        of << " " << stats[j].ci_95_hi.dlambdap;

        // 37,38,39,40,41
        of << " " << stats[j].median.tf         << " " << stats[j].disp20.tf;
        of << " " << stats[j].disp80.tf         << " " << stats[j].ci_95_lo.tf;
        of << " " << stats[j].ci_95_hi.tf;

        // 42,43,44,45,46
        of << " " << stats[j].median.vmax       << " " << stats[j].disp20.vmax;
        of << " " << stats[j].disp80.vmax       << " " << stats[j].ci_95_lo.vmax;
        of << " " << stats[j].ci_95_hi.vmax;

        // 47,48,49,50,51
        of << " " << stats[j].median.xoff       << " " << stats[j].disp20.xoff;
        of << " " << stats[j].disp80.xoff       << " " << stats[j].ci_95_lo.xoff;
        of << " " << stats[j].ci_95_hi.xoff;

        // 52,53,54,55,56
        of << " " << stats[j].median.e          << " " << stats[j].disp20.e;
        of << " " << stats[j].disp80.e          << " " << stats[j].ci_95_lo.e;
        of << " " << stats[j].ci_95_hi.e;

        // 57,58,59,60,61
        of << " " << stats[j].median.rsp        << " " << stats[j].disp20.rsp;
        of << " " << stats[j].disp80.rsp        << " " << stats[j].ci_95_lo.rsp;
        of << " " << stats[j].ci_95_hi.rsp;

        // 62,63,64,65,66
        of << " " << stats[j].median.rvirp      << " " << stats[j].disp20.rvirp;
        of << " " << stats[j].disp80.rvirp      << " " << stats[j].ci_95_lo.rvirp;
        of << " " << stats[j].ci_95_hi.rvirp;

        // 67,68,69,70,71
        of << " " << stats[j].median.e500       << " " << stats[j].disp20.e500;
        of << " " << stats[j].disp80.e500       << " " << stats[j].ci_95_lo.e500;
        of << " " << stats[j].ci_95_hi.e500;

        // 72,73,74,75,76
        of << " " << stats[j].median.lambda     << " " << stats[j].disp20.lambda;
        of << " " << stats[j].disp80.lambda     << " " << stats[j].ci_95_lo.lambda;
        of << " " << stats[j].ci_95_hi.lambda;

        // 77,78,79,80,81
        of << " " << stats[j].median.cklypin    << " " << stats[j].disp20.cklypin;
        of << " " << stats[j].disp80.cklypin    << " " << stats[j].ci_95_lo.cklypin;
        of << " " << stats[j].ci_95_hi.cklypin;

        // 82, 83, 84
        of << " " << stats[j].median.highspin   << " " << stats[j].disp20.highspin;
        of << " " << stats[j].disp80.highspin;

        // 85, 86, 87, 88
        of << " " << stats[j].spin_rank.p20[0] << " " << stats[j].spin_rank.p40[0];
        of << " " << stats[j].spin_rank.p60[0] << " " << stats[j].spin_rank.p80[0];

        // 89, 90, 91, 92
        of << " " << stats[j].spin_rank.p20[1] << " " << stats[j].spin_rank.p40[1];
        of << " " << stats[j].spin_rank.p60[1] << " " << stats[j].spin_rank.p80[1];

        // 93, 94, 95, 96
        of << " " << stats[j].spin_rank.p20[2] << " " << stats[j].spin_rank.p40[2];
        of << " " << stats[j].spin_rank.p60[2] << " " << stats[j].spin_rank.p80[2];

        // 97, 98, 99, 100 
        of << " " << stats[j].tf_rank.p20[0] << " " << stats[j].tf_rank.p40[0];
        of << " " << stats[j].tf_rank.p60[0] << " " << stats[j].tf_rank.p80[0];

        // 101, 102, 103, 104
        of << " " << stats[j].tf_rank.p20[1] << " " << stats[j].tf_rank.p40[1];
        of << " " << stats[j].tf_rank.p60[1] << " " << stats[j].tf_rank.p80[1];

        // 105, 106, 107, 108
        of << " " << stats[j].tf_rank.p20[2] << " " << stats[j].tf_rank.p40[2];
        of << " " << stats[j].tf_rank.p60[2] << " " << stats[j].tf_rank.p80[2];

        // 109,110,111,112,113
        of << " " << stats[j].median.mlumpy    << " " << stats[j].disp20.mlumpy;
        of << " " << stats[j].disp80.mlumpy    << " " << stats[j].ci_95_lo.mlumpy;
        of << " " << stats[j].ci_95_hi.mlumpy;

        // 114,115,116,117,118
        of << " " << stats[j].median.msmooth    << " " << stats[j].disp20.msmooth;
        of << " " << stats[j].disp80.msmooth    << " " << stats[j].ci_95_lo.msmooth;
        of << " " << stats[j].ci_95_hi.msmooth;

        // 119,120,121,122,123
        of << " " << stats[j].median.merger_count    << " " << stats[j].disp20.merger_count;
        of << " " << stats[j].disp80.merger_count    << " " << stats[j].ci_95_lo.merger_count;
        of << " " << stats[j].ci_95_hi.merger_count;

        // 124,125,126,127,128
        of << " " << stats[j].median.merger_mass_avg    << " " << stats[j].disp20.merger_mass_avg;
        of << " " << stats[j].disp80.merger_mass_avg    << " " << stats[j].ci_95_lo.merger_mass_avg;
        of << " " << stats[j].ci_95_hi.merger_mass_avg;

        // 129,130,131,132,133
        of << " " << stats[j].median.merger_ratio_avg    << " " << stats[j].disp20.merger_ratio_avg;
        of << " " << stats[j].disp80.merger_ratio_avg    << " " << stats[j].ci_95_lo.merger_ratio_avg;
        of << " " << stats[j].ci_95_hi.merger_ratio_avg;

        // 134,135,136,137,138
        of << " " << stats[j].median.eratio    << " " << stats[j].disp20.eratio;
        of << " " << stats[j].disp80.eratio    << " " << stats[j].ci_95_lo.eratio;
        of << " " << stats[j].ci_95_hi.eratio;

        // 139,140,141,142,143
        of << " " << stats[j].median.tu    << " " << stats[j].disp20.tu;
        of << " " << stats[j].disp80.tu    << " " << stats[j].ci_95_lo.tu;
        of << " " << stats[j].ci_95_hi.tu;

        // 144,145,146,147,148
        of << " " << stats[j].median.mdot_dyn    << " " << stats[j].disp20.mdot_dyn;
        of << " " << stats[j].disp80.mdot_dyn    << " " << stats[j].ci_95_lo.mdot_dyn;
        of << " " << stats[j].ci_95_hi.mdot_dyn;

        // 149,150,151,152,153
        of << " " << stats[j].median.smdot_dyn    << " " << stats[j].disp20.smdot_dyn;
        of << " " << stats[j].disp80.smdot_dyn    << " " << stats[j].ci_95_lo.smdot_dyn;
        of << " " << stats[j].ci_95_hi.smdot_dyn;

        of << std::endl;
    }

    of.close();

    std::cout << " done." << std::endl;

    // let's print out what we've got so far
    //for (j = z.size()-20; j < z.size(); j++) {
    //for (j = 0; j < 3; j++) {
    //    for (int i = 0; i < z[j].size(); i++) {
    //        if (!i) {
    //            std::cout << "z = " << z[j][i].z << ":" << std::endl;
    //        }
    //        std::cout << "M_vir = " << z[j][i].mvir;
    //        if (z[j][i].alive) std::cout << ", alive";
    //        else std::cout << ", dead";
    //        std::cout << std::endl;
    //    }
    //}

    //char quit = 'y';

    //while (quit == 'n') {

    //    int nhalo;
    //    std::cout << "choose a halo merger tree (0-" << z[0].size()-1 << "):";
    //    std::cin >> nhalo;

    //    // let's print out a merger tree to see if the pointers are working
    //    halo * h = &(z[0][nhalo]);
    //    std::cout << "Printing merger tree for random halo:" << std::endl;
    //    while (true) {
    //        std::cout << "M_vir = " << h->mvir  << ", z = " << h->z << ", C_nfw = " << h->cnfw;
    //        if (h->alive) std::cout << ", alive";
    //        else std::cout << ", dead";
    //        std::cout << std::endl;

    //        if (h->parent) h = h->parent;
    //        else break;
    //    }

    //    std::cout << "quit? (y/n)";
    //    std::cin >> quit;
    //}

    return 0;
}

// assuming data in bin2 format
int doSelectZ (double ** &data, int num_lines, int num_fields, int x, std::string outfile_name) {

    double a = 0;
    int count = 0;

    switch (x) {
        case 0:     a = 1.00231; break;
        case 1:     a = 0.50112; break;
    }

    for (int i = 0; i < num_lines; i++) {

        if (data[i][0] == a) count++;
    }

    std::cout << "Found " << count << " halos at a = " << a << std::endl;

    double ** subset = new double *[num_fields];
    for (int i = 0; i < num_fields; i++) subset[i] = new double[count];

    count = 0;

    for (int i = 0; i < num_lines; i++) {

        if (data[i][0] == a) {

            for (int j = 0; j < num_fields; j++) {
                
                subset[j][count] = data[i][j];
            }
            count++;
        }
    }

    writeOutputFileBinary (outfile_name, subset, count, num_fields);

}

//====================================================
// doConvertMergerTrack:
// 
// data:        selected merger tree catalog in bin2 format (i.e. data[num_lines][num_fields])
// num_lines:
//              number of lines in the catalog.
// num_fields:
//              number of fields for each line (halo).
//====================================================
int doConvertMergerTracks (double ** &data, int num_lines, int num_fields) {

    struct track {
        long int id = -1;
        long int descid = -1;
        double mvir = -1;
        double mar = 0;
        bool mmp = false;
        int track_num = 0;
    };

    // sorting comparator for rank ordering by mvir in track vector
    struct sort_comp_vtrack {
        bool operator()(const track &left, const track &right) {
            return (left.mvir > right.mvir);
        }
    };

    std::unordered_map <long int, track> track_map;
    std::unordered_map<long int, track>::const_iterator fetch;

    std::unordered_map <int, double> mpeak_map;
    std::unordered_map<int, double>::const_iterator fetch_mpeak;

    // determine what final scale of merger trees is
    double final_a = data[0][0];

    // we want to reorganize the merger trees by mass (i.e. list most massive tree first, then less massive).
    // vector to store pair of halo root line number, halo root mass.
    // this will allow us to look up tree root in data after sorting by mass.
    std::vector<std::pair<double,long int>> root_mvir;

    int list = 1;
    double min_a = final_a;

    // for each root halo, save mass and line number
    for (int i = 0; i < num_lines; i++) {
        if (data[i][0] == final_a) {
            root_mvir.push_back(std::make_pair(data[i][10],i));
        }
        // list reflects the total number of timesteps/redshifts in the provided merger trees
        if (data[i][0] < min_a) {
            min_a = data[i][0];
            list++;
        }
    }

    // now sort by mvir
    std::sort(root_mvir.begin(),root_mvir.end());

    std::ofstream of;
    setOutputFileName();
    of.open(output.c_str(),std::ofstream::out);

    for (int i = root_mvir.size()-1; i >= 0; i--) {
        std::cout << (long int)data[root_mvir[i].second][1] << ": " << log10(root_mvir[i].first) << " , " << i << std::endl;
    }

    int tree_count = 1;

    // create datastructure to store xrot, zrot, list_count, tree_count
    struct extra_params {
        double xrot = 0;
        double zrot = 0;
        long int list_count = 0;
        long int tree_count = 0;
        double mpeak = 0;
    };

    extra_params extras[num_lines];

    // now lets output first tree in halo track format
    for (int i = root_mvir.size()-1; i >= 0; i--) {

        std::vector<track> vtrack;

        int track_num = 1;
        int list_count = list;
        int j = root_mvir[i].second;
        double curr_a = final_a;

        // loop over all halos in this merger tree
        while (true) {

            // check if we've reached a new scale or end of file
            if (j == num_lines || data[j][0] != curr_a) {
                if (j != num_lines) {
                    curr_a = data[j][0];    // set new current scale
                    list_count--;           // decrement to next scale
                }

                // now we need to go through and set tracks

                // sort by mvir
                std::sort(vtrack.begin(),vtrack.end(),sort_comp_vtrack());

                // for root halo, just insert right away
                if (track_num == 1) {
                    vtrack[0].track_num = track_num++;
                    track_map[vtrack[0].id] = vtrack[0];

                    mpeak_map[track_num] = vtrack[0].mvir; 
                }

                // for other halos, we need to check descendent to determine track
                else {

                    // loop through all halos at this scale and look up descendents in map
                    for (int k = 0; k < vtrack.size(); k++) {
                                
                        fetch = track_map.find(vtrack[k].descid);

                        if (fetch == track_map.end()) {
                            std::cout << "Error: halo descendent id not found in track map!" << std::endl;
                            return 0;
                        }

                        else {
                            track t = fetch->second;
                            
                            // if mmp has already been found, need to start new track
                            if (t.mmp) {

                                vtrack[k].track_num = track_num++;
                                track_map[vtrack[k].id] = vtrack[k];
                            }

                            // this halo is mmp of descendent
                            else {

                                // update descendent in map 
                                t.mmp = true;
                                t.mar = t.mvir-vtrack[k].mvir; // this is not currently implemented correctly.
                                                               // need to do proper Mdot
                                track_map[t.id] = t;

                                // insert this halo with same track as desc
                                vtrack[k].track_num = t.track_num;
                                track_map[vtrack[k].id] = vtrack[k];
                            }

                            // update mpeak for all halos at this scale
                            fetch_mpeak = mpeak_map.find(vtrack[k].track_num);

                            if (fetch_mpeak == mpeak_map.end()) {
                                mpeak_map[vtrack[k].track_num] = vtrack[k].mvir;
                            }

                            else {
                                if (vtrack[k].mvir > fetch_mpeak->second) {
                                    mpeak_map[vtrack[k].track_num] = vtrack[k].mvir;
                                }
                            }
                        }
                    }
                }

                // clear track vector for next scale
                vtrack.clear();         
            }

            // include condition to stop when next tree is reached
            // or when end of file is reached
            if (j == num_lines) break;
            if (j > root_mvir[i].second && data[j][0] == final_a) break;


            track t;
            t.id = (long int)data[j][1];
            t.descid = (long int)data[j][3];
            t.mvir = data[j][10];
            t.mmp = false;
            vtrack.push_back(t);

            // now we need to print additional info like rotation components, track, list
            double ax = data[j][48]; double ay = data[j][49]; double az = data[j][50];
            double xrot = 0., zrot = 0., ar;
            ar = sqrt(ax*ax + ay*ay + az*az);
            if (ar != 0) {
                xrot = acos (az/ar) * 180./ M_PI;
                if (xrot > 90) xrot = -(xrot-90);
                if (xrot < -90) xrot = -(xrot+90);

                zrot = atan (ay/ax) * 180./ M_PI;
                if (zrot > 90) zrot = -(zrot-90);
                if (zrot < -90) zrot = -(zrot+90);
            }

            extras[j].xrot = xrot;
            extras[j].zrot = zrot;
            extras[j].list_count = list_count;
            extras[j].tree_count = tree_count;

            j++;
        }

        // save end position, then reiterate through all halos in tree
        int j_end = j;

        // for each halo in tree, update mpeak value
        for (j = root_mvir[i].second; j < j_end; j++) {

            // find appropriate track num for halo
            fetch = track_map.find(data[j][1]);
            track t;

            if (fetch == track_map.end()) {
                std::cout << "Error: halo id not found in track map!" << std::endl;
                return 0;
            }

            else {
                t = fetch->second;
            }

            // look up corresponding mpeak for that track
            fetch_mpeak = mpeak_map.find(t.track_num);

            if (fetch_mpeak == mpeak_map.end()) {
                std::cout << "Error: track num not found in mpeak map!" << std::endl;
                return 0;
            }

            // save in extra params for that halo
            else {
                extras[j].mpeak = fetch_mpeak->second;
            }
        }

        // reset j (unnecessary?)
        j = j_end;

        // clear map for next tree, since track nums will be reused
        mpeak_map.clear();

        tree_count++;
    }

    // loop over trees until all scales have been output
    for (int l = list; l > 0; l--) {

        // loop over all trees and output scale l
        for (int i = root_mvir.size()-1; i >= 0; i--) {

            int j = root_mvir[i].second;

            // if this tree has finished, skip
            if (j == -1) continue;

            double curr_a = data[j][0];

            // loop over all halos at scale l in this tree
            while (true) {

                // include condition to stop when end of data is reached
                if (j == num_lines) {
                    root_mvir[i].second = -1;
                    break;
                }

                // include condition to stop when next scale is reached            
                if (j > root_mvir[i].second && data[j][0] == final_a) {
                    root_mvir[i].second = -1;
                    break;
                }

                // stop when next scale is reached
                if (curr_a != data[j][0]) {
                    root_mvir[i].second = j;
                    break;
                }

                // loop over fields in each halo
                for (int k = 0; k < num_fields; k++) {
                    switch (k) {
                        case 0:     // scale
                        case 10:    // mvir
                        case 11:    // rvir
                        case 12:    // rs
                        case 15:    // almm
                        case 16:    // vmax
                        case 17:    // x
                        case 18:    // y
                        case 19:    // z
                        case 20:    // vx
                        case 21:    // vy
                        case 22:    // vz
                        case 26:    // lambda
                        case 35:    // tf
                        case 37:    // rs_klypin
                        case 43:    // xoff
                        case 44:    // voff
                        case 45:    // lambadp
                        case 46:    // b/a
                        case 47:    // c/a
                        case 51:    // b/a500
                        case 52:    // c/a500
                        case 56:    // T/|U|
                            of << data[j][k] << " ";
                            break;
                        case 1:     // id
                        case 3:     // descid
                            of << (long int)data[j][k] << " ";
                            break;
                        default:
                            break;
                    }
                }

                fetch = track_map.find(data[j][1]);

                track t;

                if (fetch == track_map.end()) {
                    std::cout << "Error: halo id not found in track map!" << std::endl;
                    return 0;
                }

                else {
                    t = fetch->second;
                }

                of << extras[j].xrot << " " << extras[j].zrot << " " << extras[j].list_count << " ";
                of << extras[j].tree_count << " " << t.track_num << " " << t.mar << std::endl;

                j++;
            }
        }
    }
            
    return 0;
}

//====================================================
// doTreeAnalysis:
//              This is a tree walking implementation designed to be able to
//              easily answer questions concerning merger trees.  It expects an input file
//              containing each merger tree listed consecutively.
// data:        selected merger tree catalog in bin2 format (i.e. data[num_lines][num_fields])
// num_lines:
//              number of lines in the catalog.
// num_fields:
//              number of fields for each line (halo).
//====================================================
int doTreeAnalysis (double ** &data, int num_lines, int num_fields, std::string outfile_name) {

    //struct time_avg_prop {

    //    // for each quantity saved, we need vector of property
    //    double mdot_dyn = -1;
    //    double mdot_dyn_lumpy = -1;
    //    double mdot_dyn_smooth = -1;
    //    double lambda_dyn = -1;
    //    double lambda_dyn_disp = -1;
    //    double tdyn = -1;
    //    void * 
    //};

    //num_lines = 1000000;

    // we need a data structure for each node on the tree
    struct node {
        // ids
        long int id = -1;
        long int desc_id = -1;
        long int upid = -1;

        bool root = false;  // is this a root node?
        long int i = -1;    // index to location in raw data

        // pointer to descendent node
        node * child;

        // pointer to host node (upid) (same scale factor)
        node * host;

        // vector containing pointers to parent nodes
        std::vector<node *> parents;

        // vector containing pointers to subhalos (those that list current
        // halo as upid)
        std::vector<node *> subhalos;

        int num_prog = -1;

        // halo properties from merger trees
        bool   mmp = false;
        double scale = -1;
        double mvir = -1;
        double rs = -1;
        double rvir = -1;
        double almm = -1;
        double vmax = -1;
        double lambda = -1;
        double lambdap = -1;
        double tf = -1;
        double rs_klypin = -1;
        double xoff = -1;
        double b_a = -1;
        double c_a = -1;
        double b_a500 = -1;
        double c_a500 = -1;
        double tu = -1;

        // custom halo properties (computed here)
        double z = -1;
        double mpeak = -1;
        double a_mpeak = -1;
        double cnfw = -1;
        double e = -1;
        double e500 = -1;
        double xrot = -1;
        double zrot = -1;
        long int track = -1;
        long int list_count = -1;
        long int tree_count = -1;
        double almms = -1;
        double almm2 = -1;      // my determination
        double alminm = -1;
        double msmooth = -1;
        double msmooth_frac = -1;
        double mlumpy = -1;
        double mlumpy_frac = -1;
        double merger_count = 0;
        double merger_mass_avg = 0;
        double high_tf = 0;
        bool ongoing_mm = false;
        bool ongoing_minm = false;
        double ongoing_mm_scale = -1;
        double ongoing_minm_scale = -1;
        double maxtf_peak = -1;
        double a_maxtf_peak = -1;
        double min_bsr_peak = -1;
        double a_min_bsr_peak = -1;
        double min_bsr_lmm = -1;
        double a_min_bsr_lmm = -1;

        double mdot_dyn = -1;
        double mdot_dyn_lumpy = -1;
        double mdot_dyn_smooth = -1;
        double lambda_dyn = -1;
        double lambda_dyn_disp = -1;
        //double tdyn = -1;

        double spin_peak;   // peak lambda_B since mpeak
        double tu_peak;     // peak T/|U| since mpeak

        double npt = 0; // number of pass-throughs since mpeak
        double mass_lost_as_sh = 0; // total mass loss since mpeak that
                                    // occured while halo was a subhalo
        double mass_lost_mpeak = 0; // total mass loss since mpeak
    };

    // we need a top level vector to point to each tree root.
    // i'm using pointers here so that when we sort the vector, any pointers
    // in use will still be valid.
    std::vector<node *> forest;

    // vector to store all nodes in all trees
    std::vector<node> vnode;
    vnode.reserve(num_lines);

    std::unordered_map<long int, node *> id_map;
    std::unordered_map<long int, node *>::const_iterator id_fetch;

    double final_a = data[0][0];

    int list = 1;
    double min_a = final_a;
    int tree_count = 1;

    // now, we need iterate through the merger tree input and transform into a forest
    for (int i = 0; i < num_lines; i++) {

        node n;

        // now push this node onto the node vector and make further changes from there
        vnode.push_back(n);
        node * np = &vnode.back();

        np->id = data[i][BP_TREE.id];
        np->desc_id = data[i][BP_TREE.desc_id];
        np->upid = data[i][BP_TREE.upid];

        np->num_prog = data[i][BP_TREE.num_prog];

        id_map[np->id] = np;

        np->i = i;
        
        np->scale = data[i][BP_TREE.a];
        np->mvir = data[i][BP_TREE.mvir];
        np->rvir = data[i][BP_TREE.rvir];
        np->rs = data[i][BP_TREE.rs];
        np->almm = data[i][BP_TREE.almm];
        np->vmax = data[i][BP_TREE.vmax];
        np->lambda = data[i][BP_TREE.spinp];
        np->tf = data[i][BP_TREE.tf];
        np->rs_klypin = data[i][BP_TREE.rsk];
        np->xoff = data[i][BP_TREE.xoff];
        np->lambdap = data[i][BP_TREE.spinb];
        np->b_a = data[i][BP_TREE.bta];
        np->c_a = data[i][BP_TREE.cta];
        np->b_a500 = data[i][BP_TREE.bta500];
        np->c_a500 = data[i][BP_TREE.cta500];
        np->tu = data[i][BP_TREE.tu];

        np->z = 1./np->scale - 1.;
        np->cnfw = np->rvir/np->rs;
        np->e = 1.-sqrt(np->b_a*np->b_a+np->c_a*np->c_a)/sqrt(2.);
        np->e500 = 1.-sqrt(np->b_a500*np->b_a500 + np->c_a500*np->c_a500)/sqrt(2.);

        // now we need to compute rotation components
        double ax = data[i][BP_TREE.ax]; double ay = data[i][BP_TREE.ay]; double az = data[i][BP_TREE.az];
        double xrot = 0., zrot = 0., ar;
        ar = sqrt(ax*ax + ay*ay + az*az);
        if (ar != 0) {
            xrot = acos (az/ar) * 180./ M_PI;
            if (xrot > 90) xrot = -(xrot-90);
            if (xrot < -90) xrot = -(xrot+90);

            zrot = atan (ay/ax) * 180./ M_PI;
            if (zrot > 90) zrot = -(zrot-90);
            if (zrot < -90) zrot = -(zrot+90);
        }

        np->xrot = xrot;
        np->zrot = zrot;

        // list reflects the total number of timesteps/redshifts in the provided merger trees
        if (np->scale < min_a) {
            min_a = np->scale;
            list++;
        }

        // this is a root node
        if (np->scale == final_a) {
            //if (forest.size() > 100) break;
            np->root = true;
            np->mmp = true;
            forest.push_back(np);
            
            //std::cout << "Root " << forest.size() << ": " << log10(forest.back()->mvir) << ", " << forest.back()->scale << std::endl;
        }

        // this is not a root node, so we need to find child (and update child's parents vector)
        else {
            id_fetch = id_map.find(np->desc_id);

            if (id_fetch == id_map.end()) {
                std::cout << "Error: descendent halo " << np->desc_id << " not found!" << std::endl;
            }

            else {
                np->child = id_fetch->second;
                //std::cout << np->child << " " << &vnode[vnode.size()-2] << " " << np->id << " " << np->desc_id << std::endl;
                np->child->parents.push_back(np);

                //std::cout << "Inserting node at scale " << np->scale << ", desc scale " << np->child->scale << ", " << log10(np->mvir) << "/" << log10(np->child->mvir) << std::endl;
                //if (np->child->parents.back()->mvir != np->mvir) {
                //    std::cout << "problem" << std::endl;
                //}
            }
        }
    }

    int err_count = 0;
    for (int i = 0; i < vnode.size(); i++) {

        node * np = &vnode[i];

        // check if halo has host halo and update subhalo pointers if so
        if (np->upid != -1) {

            id_fetch = id_map.find(np->upid);

            if (id_fetch == id_map.end()) {
                err_count++;
                //std::cout << "Error: most massive host halo " << np->upid << " not found!" << std::endl;
            }

            else {
                np->host = id_fetch->second;
                np->host->subhalos.push_back(np);
            }
        }
    }

    //int count  = 0;
    //double avg_ns = 0;
    //for (int i = 0; i < vnode.size(); i++) {

    //    if (vnode[i].subhalos.size() > 0) {
    //        count++;
    //        avg_ns+=vnode[i].subhalos.size();
    //        if (count < 5) {
    //            std::cout << "Num subhalos: " << vnode[i].subhalos.size() << " and size of node " << sizeof(vnode[i]) << std::endl;
    //        std::cout << "before: " << &vnode[i+1]-&vnode[i] << std::endl;
    //        for (int j = 0; j < 100; j++) {
    //            vnode[i].subhalos.push_back(&vnode[i+1]);
    //        }
    //        std::cout << "after: " << &vnode[i+1]-&vnode[i] << " and size of node " << sizeof(vnode[i]) << std::endl;
    //        }
    //    }
    //}

    //std::cout << "Node size: " << sizeof(node) << std::endl;
    //std::cout << "Number of halos with subhalos: " << count << std::endl;
    //std::cout << "Avg number of subhalos among halos with subhalos: " << avg_ns/count << std::endl;
    //std::cout << "Number of halos with missing hosts: " << err_count << std::endl;
    //std::cout << "Total number of halos processed: " << vnode.size() << std::endl;

    //std::cout << "# Subhalos " << forest[0]->subhalos.size() << std::endl;
    //std::cout << "# Subhalos (a<1) " << forest[0]->parents[0]->subhalos.size() << std::endl;

    std::cout << "list: " << list << std::endl;

    std::function<double(node*,int,int)> mpeak;
    mpeak = [&mpeak](node * n, int tree_count, int list)->double {
        double max_mpeak = 0;
        double a_max_mpeak = -1;

        double maxtf_peak = -1;
        double a_maxtf_peak = -1;
        double min_bsr_peak = -1;
        double a_min_bsr_peak = -1;

        double min_bsr_lmm = -1;
        double a_min_bsr_lmm = -1;

        double spin_peak = -1;
        double tu_peak = -1;

        double mass_lost_as_sh = -1;
        double mass_lost_mpeak = -1;

        n->list_count = list;           // set list count (timestep num, 1 being earliest timestep in whole file)
        n->tree_count = tree_count;
        if (n->parents.size() == 0) {
            n->mpeak = n->mvir;
            n->a_mpeak = n->scale;
            n->maxtf_peak = n->tf;
            n->a_maxtf_peak = n->scale;
            n->min_bsr_peak = 1.;
            n->a_min_bsr_peak = n->scale;
            n->min_bsr_lmm = 1.;
            n->a_min_bsr_lmm = n->scale;
            n->spin_peak = n->lambdap;
            n->tu_peak = n->tu;
            n->mass_lost_as_sh = 0;
            n->mass_lost_mpeak = 0;

            return n->mpeak;
        }

        // sort by halo mass, mmp first
        std::sort(n->parents.begin(),n->parents.end(),[](node * left, node * right){return left->mvir > right->mvir;});

        // set mmp for all tree branches (I don't think this is correct, should only be for actual MMPB)
        n->parents.front()->mmp = true;

        // find greatest peak mass among progenitors
        for (int i = 0; i < n->parents.size(); i++) {

            n->parents[i]->mpeak = mpeak(n->parents[i],tree_count,list-1);
            //std::cout << "mpeak at " << n->parents[i]->scale << ", " << log10(n->parents[i]->mvir) << ", " << log10(n->parents[i]->mpeak) << std::endl;;
            if (n->parents[i]->mpeak > max_mpeak) {
                max_mpeak = n->parents[i]->mpeak;
                a_max_mpeak = n->parents[i]->a_mpeak;

                maxtf_peak = n->parents[i]->maxtf_peak;
                a_maxtf_peak = n->parents[i]->a_maxtf_peak;
                min_bsr_peak = n->parents[i]->min_bsr_peak;
                a_min_bsr_peak = n->parents[i]->a_min_bsr_peak;
                min_bsr_lmm = n->parents[i]->min_bsr_lmm;
                a_min_bsr_lmm = n->parents[i]->a_min_bsr_lmm;

                spin_peak = n->parents[i]->spin_peak;
                tu_peak = n->parents[i]->tu_peak;

                mass_lost_as_sh = n->parents[i]->mass_lost_as_sh;
                mass_lost_mpeak = n->parents[i]->mass_lost_mpeak;
            }
        }

        // check if current timestep is mpeak
        if (n->mvir > max_mpeak) {
            max_mpeak = n->mvir;
            a_max_mpeak = n->scale;

            // if we're at mpeak, then reset 'since mpeak' quantities
            maxtf_peak = n->tf;
            a_maxtf_peak = n->scale;
            min_bsr_peak = 1.;
            a_min_bsr_peak = n->scale;

            mass_lost_as_sh = 0;
            mass_lost_mpeak = 0;
        }

        // if not at mpeak, then check if current timestep values should be used rather than
        // max progenitor values for 'since mpeak' quantities
        else {
            if (n->tf > maxtf_peak) {
                maxtf_peak = n->tf;
                a_maxtf_peak = n->scale;
            }

            if (n->mvir/max_mpeak < min_bsr_peak) {
                min_bsr_peak = n->mvir/max_mpeak;
                a_min_bsr_peak = n->scale;
            }

            double dm = n->parents[0]->mvir - n->mvir;
            if (dm > 0) {
                if (n->upid != -1 && n->parents[0]->upid != -1) {
                    mass_lost_as_sh += dm;
                }
                mass_lost_mpeak += dm;
            }
        }

        // now do similar for the 'since almm' properties
        if (n->scale == n->almm) {
            min_bsr_lmm = 1.;
            a_min_bsr_lmm = n->scale;

            spin_peak = n->lambdap;
            tu_peak = n->tu;
        }

        else {
            
            if (n->mvir/max_mpeak < min_bsr_lmm) {
                min_bsr_lmm = n->mvir/max_mpeak;
                a_min_bsr_lmm = n->scale;
            }

            
            if (n->lambdap > spin_peak) {
                spin_peak = n->lambdap;
            }

            if (n->tu > tu_peak) {
                tu_peak = n->tu;
            }
        }

        // need to set this here since not returning like max_mpeak
        n->a_mpeak = a_max_mpeak;

        n->maxtf_peak = maxtf_peak;
        n->a_maxtf_peak = a_maxtf_peak;
        n->min_bsr_peak = min_bsr_peak;
        n->a_min_bsr_peak = a_min_bsr_peak;

        n->min_bsr_lmm = min_bsr_lmm;
        n->a_min_bsr_lmm = a_min_bsr_lmm;

        n->spin_peak = spin_peak;
        n->tu_peak = tu_peak;

        n->mass_lost_as_sh = mass_lost_as_sh;
        n->mass_lost_mpeak = mass_lost_mpeak;

        return max_mpeak;
    };

    // sort by halo mass
    std::sort(forest.begin(),forest.end(),[](node * left, node * right){return left->mvir > right->mvir;});

    std::cout << "Computing Mpeak..." << std::endl;

    // compute mpeak for each halo
    for (int i = 0; i < forest.size(); i++) {
        forest[i]->mpeak = mpeak(forest[i],i+1,list);

        //std::cout << log10(forest[i]->mvir) << " " << log10(forest[i]->mpeak) << std::endl;
    }

    std::function<double(node*)> tfpeak;
    tfpeak = [&tfpeak](node * n)->double {

        if (n->parents.size() == 0) {
            //n->high_tf = 0.;
            return 0.;
        }

        for (int i = 0; i < n->parents.size(); i++) {

            n->parents[i]->high_tf = tfpeak(n->parents[i]);
            //std::cout << "mpeak at " << n->parents[i]->scale << ", " << log10(n->parents[i]->mvir) << ", " << log10(n->parents[i]->mpeak) << std::endl;
        }
        // return value will set high_tf for this halo.
        // always return true if tf > 1
        if (n->tf > 1.) return 1.;
        // if this halo is at mpeak, then parents don't matter. return false
        else if (n->mvir == n->mpeak) return 0.;
        // not mpeak, so parent's value (mmp) is what we want
        else return n->parents[0]->high_tf;
    };

    std::cout << "Computing TF peaks..." << std::endl;

    // compute mpeak for each halo
    for (int i = 0; i < forest.size(); i++) {
        forest[i]->high_tf = tfpeak(forest[i]);

        //std::cout << log10(forest[i]->mvir) << " " << log10(forest[i]->mpeak) << std::endl;
    }

    std::cout << "Computing tracks..." << std::endl;

    // now set track num for each
    for (int i = 0; i < forest.size(); i++) {

        int track = 1;

        // set track for root, increment so next track set will be 2
        forest[i]->track = track++;

        // vector for storing parent nodes at each scale
        std::vector<node *> gather;

        // populate initial gather with root nodes parents
        for (int j = 0; j < forest[i]->parents.size(); j++) gather.push_back(forest[i]->parents[j]);

        while (true) {

            // since we will be pushing to gather, need to keep end of loop value constant
            int scale_size = gather.size();

            // loop over each parent, set track num, and push all grandparents onto gather
            for (int j = 0; j < scale_size; j++) {
                if (gather[j]->mmp) gather[j]->track = gather[j]->child->track;
                else gather[j]->track = track++;

                if (gather[j]->parents.size()) {
                    for (int k = 0; k < gather[j]->parents.size(); k++) gather.push_back(gather[j]->parents[k]);
                }
            }

            // no more parents added, so we've reached the end of the tree
            if (gather.size() == scale_size) {
                break;
            }

            // reduce gather vector to just newly added nodes
            gather.erase(gather.begin(),gather.begin()+scale_size);
        }
    }

    // now set stellar major merger scale for each.
    // to determine the stellar major merger, we need to assume some stellar to halo mass ratio.
    // this does not evolve too drastically at low redshift, so we can assume a constant ratio for now.
    // using this ratio M*(Mh), check each merger (most massive vs 2nd most massive) halo progenitor
    // and indicate those which stellar mass ratio greater than 1:3.

    // helper function for mstellar(). computes f(x) in peters fitting formula.
    auto mstellar_f = [](double x)->double {
        double alpha = -1.84754, delta = 3.85767, gamma = 0.44557;
        double r1 = -log10(pow(10.,alpha*x)+1.);
        double r2 = delta*pow(log10(1.+exp(x)),gamma);
        r2 /= 1.+exp(pow(10.,-x));
        return r1+r2;
    };

    // lambda function for converting halo mass to stellar mass. expects input as either log10(mvir)
    // or Mvir (will detect whether its logged or now) with mvir in units of Msol/h.
    // returns M* in units Msol/h (not logged)
    auto mstellar = [&mstellar_f](double mhalo)->double {
        
        // check if mhalo is logged already or not. we want it not logged for this function.
        if (mhalo<100.) mhalo = pow(10.,mhalo);

        // using Peter's fitting forumla (Eq. 3 from Behroozi et al 2013).
        // expects input Mvir in units of Msol/h, and output is log10(M*)
        // with M* in units of Msol/h^2
        double le = -1.959715, lm1 = 11.483698, h = 0.678;
        double ret = le+lm1+mstellar_f(log10(mhalo)-lm1)-mstellar_f(0.);

        // peter's fitting formula returns M* in units of Msun/h^2, but we want Msun/h
        return pow(10.,ret)/h;
    };

    std::cout << "Computing almms... " << std::endl;

    std::function<double(node*,double)> almms;
    almms = [&almms,&mstellar](node * n, double final_ms)->double {

        // base case #1
        // we'll treat end of tree as major merger, since peter seems to do that for almm
        if (n->parents.size() == 0) {
            if (final_ms == -1) {
                // we'll treat end of tree as major merger, since peter seems to do that for almm
                return n->scale;
            }
            // but if we're doing comparison to final_ms, we want to be able to select halos that
            // never had stellar mm, so let's signify that with -1
            else {
                return -1;
            }
        }

        // sort by halo mass, mmp first.
        // would want to do this generally, though we've already sorted these by mass in the
        // mpeak function above, no need to do so again.
        //std::sort(n->parents.begin(),n->parents.end(),[](node * left, node * right){return left->mvir > right->mvir;});

        // check if a stellar major merger happened to this halo's parents
        if (n->parents.size() > 1) {

            double ms_ratio = 0;

            // -1 means to do regular major merger determination by comparing two most massive progenitors
            if (final_ms == -1) {
                // since we are using mpeak to determine the stellar mass, but parents are rank ordered
                // by mvir, it is possible to have a merger ratio less than one (2nd mmp had higher mpeak
                // than mmp).
                ms_ratio = mstellar(n->parents[0]->mpeak) / mstellar(n->parents[1]->mpeak);

                // using mass ratio of 1:3 to define major merger. only comparing most massive
                // to second most massive.

                // this is a major merger
                if (ms_ratio < 3.) {
                    n->almms = n->scale;
                }
            }

            // positive value of final_ms means to use final_ms as the stellar mass comparison.  this means
            // always comparing progenitor stellar masses to the same final halo stellar mass.
            else {
                // we're now determining merger ratio relative to (mpeak) stellar mass of final halo
                ms_ratio = final_ms / mstellar(n->parents[1]->mpeak);

                // this is a major merger
                if (ms_ratio < 10.) {
                    n->almms = n->scale;
                }
            }
        }

        // even if we set almms here, still need to recurse down through parents to 
        // set almms for them.
        for (int i = 0; i < n->parents.size(); i++) {

            n->parents[i]->almms = almms(n->parents[i],final_ms);

            //std::cout << "almms at " << n->parents[i]->scale << ", " << log10(n->parents[i]->mvir) << ", almm " << n->parents[i]->almm << ", " << n->parents[i]->almms << std::endl;
        }

        // we already set n->almms, so this is redundant, but it will be set again
        // using this return value.
        if (n->almms != -1) return n->almms;

        // otherwise, we are only interested in the most massive progenitor branch
        // for return value.
        else return n->parents[0]->almms;
    };

    // sort by halo mass
    std::sort(forest.begin(),forest.end(),[](node * left, node * right){return left->mvir > right->mvir;});

    // compute almms for each halo
    for (int i = 0; i < forest.size(); i++) {
        //forest[i]->almms = almms(forest[i],-1);
        forest[i]->almms = almms(forest[i],mstellar(forest[i]->mpeak));

        //std::cout << log10(forest[i]->mvir) << " " << log10(forest[i]->mpeak) << std::endl;
    }

    std::cout << "Computing almm/alminm... " << std::endl;

    double me_gt_pb = 0, me_lt_pb = 0, me_eq_pb = 0;
    double total_count = 0;

    std::function<void(node*)> almm;
    almm = [&almm,&me_gt_pb,&me_lt_pb,&me_eq_pb,&total_count](node * n)->void {

        bool flag = false;

        // base case #1
        // we'll treat end of tree as major merger, since peter seems to do that for almm
        if (n->parents.size() == 0) {

                if (n->almm != n->scale) {
                    //std::cout << "Problem, expected EOL MM at a = " << n->scale << ", but catalog says a = " << n->almm << std::endl;
                    //std::cout.flush();
                    flag = true;
                }
                n->almm2 = n->scale;
                n->alminm = n->scale;

                // set number of pass-throughs
                n->npt = 0;

                if (flag) {
                    //std::cout << "Ended up setting almm2 to " << n->almm2 << std::endl;
                    //std::cout.flush();

                    if (n->almm2>n->almm) me_gt_pb += 1.;
                    if (n->almm2<n->almm) me_lt_pb += 1.;
                    if (n->almm2==n->almm) me_eq_pb += 1.;
                    total_count += 1.;
                }

                return;
        }

        // still need to recurse down through parents to set almm for them.
        for (int i = 0; i < n->parents.size(); i++) {

            almm(n->parents[i]);

        }

        // for subhalos, set almm to current a
        //if ((n->upid != -1) && (n->parents[0]->upid == -1)) {

        //    n->almm2 = n->scale;
        //    n->alminm = n->scale;

        //    if (n->almm2>n->almm) me_gt_pb += 1.;
        //    if (n->almm2<n->almm) me_lt_pb += 1.;
        //    if (n->almm2==n->almm) me_eq_pb += 1.;
        //    total_count += 1.;
        //}

        // check if a major merger is underway
        if (n->subhalos.size() > 0) {

            // sort by halo mass, mmp first.
            // would want to do this generally, though we've already sorted these by mass in the
            // mpeak function above, no need to do so again.
            std::sort(n->subhalos.begin(),n->subhalos.end(),[](node * left, node * right){return left->mvir > right->mvir;});
            std::sort(n->parents.begin(),n->parents.end(),[](node * left, node * right){return left->mvir > right->mvir;});

            double mratio = 0;

            //mratio = n->parents[1]->mvir / n->parents[0]->mvir;
            //int i = 1; bool upid_flag = false;
            //while (n->parents[i]->upid != n->parents[0]->id) {
            //    i++;
            //    if (i == n->parents.size()-1) {
            //        upid_flag = true;
            //        break;
            //    }
            //}

            //mratio = n->parents[1]->mvir / (n->mvir - n->parents[1]->mvir);
            mratio = n->subhalos[0]->mvir / (n->mvir - n->subhalos[0]->mvir);

            // using mass ratio of 0.3:1 to define major merger, 0.1:1 for minor merger. only comparing most massive
            // to second most massive.

            // this is at least a minor merger
            if (mratio > 0.2) {

                //if (!n->parents[0]->ongoing_minm) {
                //    //n->alminm = n->scale;
                //    n->ongoing_minm_scale = n->scale;
                //}
                //else n->ongoing_minm_scale = n->parents[0]->ongoing_minm_scale;

                //n->ongoing_minm = true;

                n->alminm = n->scale;

                // this is a major merger
                if (mratio > 0.3) {

                    //if (!n->parents[0]->ongoing_mm) {

                    //    //if (n->almm != n->scale) {
                    //    //    //std::cout << "Problem, expected 0.3:1 MM at a = " << n->scale << ", but catalog says a = " << n->almm << std::endl;
                    //    //    //std::cout << "MMP mass is " << n->parents[0]->mvir << ", and 2nd MMP mass is " << n->parents[1]->mvir << std::endl;
                    //    //    //std::cout.flush();
                    //    //    flag = true;
                    //    //}
                    //    //n->almm2 = n->scale;
                    //    n->ongoing_mm_scale = n->scale;
                    //}
                    //else n->ongoing_mm_scale = n->parents[0]->ongoing_mm_scale;

                    //n->ongoing_mm = true;

                    n->almm2 = n->scale;
                }

                //else if (n->parents[0]->ongoing_mm) {

                //    node * np = n->parents[0];

                //    while (np->ongoing_mm) {
                //        np->almm2 = np->ongoing_mm_scale;
                //        if (!np->parents.size()) break;
                //        np = np->parents[0];
                //    }
                //}

                //else n->ongoing_mm = false;
            }

            //else if (n->parents[0]->ongoing_minm) {

            //    //n->parents[0]->alminm = n->parents[0]->ongoing_minm_scale; 

            //    node * np = n->parents[0];

            //    while (np->ongoing_minm) {
            //        np->alminm = np->ongoing_minm_scale;
            //        if (!np->parents.size()) break;
            //        np = np->parents[0];
            //    }
            //}

            //else n->ongoing_minm = false;
        }

        //else if (n->parents[0]->ongoing_minm) n->parents[0]->alminm = n->parents[0]->ongoing_minm_scale;
        //else if (n->parents[0]->ongoing_mm) n->parents[0]->almm2 = n->parents[0]->ongoing_mm_scale;

        // if not set here, we want to set to the most massive progenitor branch merger scales
        if (n->almm2 == -1) n->almm2 = n->parents[0]->almm2;
        if (n->alminm == -1) n->alminm = n->parents[0]->alminm;

        // inherit number of passthroughs from parent halo
        n->npt = n->parents[0]->npt;

        // want to only report number of pass throughs since mpeak. so if we're at mpeak, set to zero
        if (n->mvir == n->mpeak) n->npt = 0;

        // look for transition from subhalo to central to increment number of pass-throughs
        if ((n->upid == -1) && (n->parents[0]->upid != -1)) n->npt += 1;

        if (flag) {
            //std::cout << "Ended up setting almm2 to " << n->almm2 << std::endl;
            //std::cout.flush();
        }

        if (n->almm2>n->almm) me_gt_pb += 1.;
        if (n->almm2<n->almm) me_lt_pb += 1.;
        if (n->almm2==n->almm) me_eq_pb += 1.;
        total_count += 1.;

        return;
    };


    // sort by halo mass
    std::sort(forest.begin(),forest.end(),[](node * left, node * right){return left->mvir > right->mvir;});

    // compute almm/alminm for each halo
    for (int i = 0; i < forest.size(); i++) {
        almm(forest[i]);
    }

    std::cout << "Found " << me_gt_pb << " halos where a_lmm2 > a_lmm out of " << total_count << " total halos (" << me_gt_pb/total_count*100. << "%)" << std::endl;
    std::cout << "Found " << me_lt_pb << " halos where a_lmm2 < a_lmm out of " << total_count << " total halos (" << me_lt_pb/total_count*100. << "%)" << std::endl;
    std::cout << "Found " << me_eq_pb << " halos where a_lmm2 == a_lmm out of " << total_count << " total halos (" << me_eq_pb/total_count*100. << "%)" << std::endl;

    std::cout << "Computing smooth/lumpy accretion... " << std::endl;

    struct acc_type {
        double mlumpy = 0;
        double msmooth = 0;
        double merger_count = 0;
        double merger_mass_avg = 0;
    };
    
    std::function<acc_type(node*)> accretion_type;
    accretion_type = [&accretion_type](node * n)->acc_type {

        // base case #1
        // we'll treat end of tree as pure smooth accretion, no lumpy, since we have no knowledge
        // of previous mergers
        if (n->parents.size() == 0) {
            acc_type ret;
            ret.msmooth = n->mvir;
            ret.mlumpy = 0;
            ret.merger_count = 0;
            ret.merger_mass_avg = 0;
            return ret;
        }

        double mlumpy_inst = 0;
        double msmooth_inst = 0;
        double merger_count = 0;

        // check if a major merger happened to this halo's parents
        // iterate through progenitors and accumulate instantaneous lumpy vs smooth accretion
        if (n->parents.size() > 1) {

            // this defines the merger ratio at which accreted material would be
            // treated as lumpy vs smooth
            //double lumpy_merger_ratio = 10.; 
            double lumpy_merger_ratio = 1e10; 

            for (int i = 1; i < n->parents.size(); i++) {

                if ((n->parents[0]->mvir / n->parents[i]->mvir) < lumpy_merger_ratio) {
                    mlumpy_inst += n->parents[i]->mvir;
                    merger_count += 1;
                }
            }
        }

        msmooth_inst = n->mvir - n->parents[0]->mvir - mlumpy_inst;

        // at the moment, we're only interested in positive accumulation, so ignore negative accretion
        if (msmooth_inst < 0) msmooth_inst = 0;

        for (int i = 0; i < n->parents.size(); i++) {

            acc_type progenitor = accretion_type(n->parents[i]);

            n->parents[i]->mlumpy = progenitor.mlumpy;
            n->parents[i]->msmooth = progenitor.msmooth;
            n->parents[i]->merger_count = progenitor.merger_count;
            n->parents[i]->merger_mass_avg = progenitor.merger_mass_avg;
            double acc_tot = n->parents[i]->mlumpy + n->parents[i]->msmooth;
            n->parents[i]->mlumpy_frac = n->parents[i]->mlumpy / acc_tot;
            n->parents[i]->msmooth_frac = n->parents[i]->msmooth / acc_tot;

            //std::cout << "almms at " << n->parents[i]->scale << ", " << log10(n->parents[i]->mvir) << ", almm " << n->parents[i]->almm << ", " << n->parents[i]->almms << std::endl;
        }

        acc_type ret;
        ret.msmooth = n->parents[0]->msmooth + msmooth_inst;
        ret.mlumpy = n->parents[0]->mlumpy + mlumpy_inst;
        ret.merger_count = n->parents[0]->merger_count + merger_count;
        if (ret.merger_count) ret.merger_mass_avg = ret.mlumpy/ret.merger_count;

        return ret;
    };

    // compute accretion type for each halo
    for (int i = 0; i < forest.size(); i++) {
        acc_type root = accretion_type(forest[i]);
        forest[i]->mlumpy = root.mlumpy;
        forest[i]->msmooth = root.msmooth;
        forest[i]->merger_count = root.merger_count;
        forest[i]->merger_mass_avg = root.merger_mass_avg;
        double acc_tot = forest[i]->mlumpy + forest[i]->msmooth;
        forest[i]->mlumpy_frac = forest[i]->mlumpy / acc_tot;
        forest[i]->msmooth_frac = forest[i]->msmooth / acc_tot;

        //std::cout << log10(forest[i]->mvir) << " " << log10(forest[i]->mpeak) << std::endl;
    }

    std::cout << "Computing time-averaged quantities..." << std::endl;

    std::string tin_path = "t_list.dat";    // file path for time data corresponding to
                                            // each redshift in the merger trees.

    std::cout << "Expecting T list file to be found at: \"" << tin_path << "\"" << std::endl;
    std::cout.flush();

    double t[list];

    // the following was used to create t_list file:
    // $ rm -f t_list.dat; cat z_list.dat | awk '{print "python cosmocalc.py "$1" 67.8 0.307 0.693"}' | /bin/bash | awk '{print $1}' >> t_list.dat

    // now, let's read in and store our time data (to correspond to each redshift
    // in the merger trees).
    // note, we are assuming the file containing time data can be found at the path below.
    std::ifstream tin (tin_path,std::ifstream::in);

    for (int i = 0; i < list; i++) {
        tin >> t[i];
    }

    tin.close();

    // computes dynamical time in units of billions of years
    auto t_dyn = [](double a)->double {

        double y = 0.307*pow(1./a,3.);
        double Omega = y/(y+0.693);
        double H = 0.06934*sqrt(y+0.693);
        double x = Omega-1.;
        double delta_vir = (18.*M_PI*M_PI + 82.*x - 39.*x*x)/Omega;
        return pow(delta_vir*3.*H*H/(8.*M_PI),-1./2.);
    };

    auto update_tavg = [&t,&list,&t_dyn](node * n1, node * n2)->int {

        // not exactly sure how to handle this case at the moment
        if (n1 == n2) return 0;

        double tdyn = t_dyn (n1->scale);
        double t1 = t[list-n1->list_count];
        double t2 = t[list-n2->list_count];
        double t2c = t[list-n2->child->list_count];

        //std::cout << "Updating time-averaged quantities (" << n1->track << "," << list-n1->list_count << "; " << n2->child->track << "," << n2->track << "," << list-n2->list_count << ")" << std::endl;
        //std::cout << "tdyn: " << tdyn << std::endl;

        // most basic quantities are those for which we only need to know values at endpoints (separated by t_dyn).
        // handle those first

        // check if we've passed dynamical time threshold.
        // but only do this once (don't repeat again for this n1 for future n2's).
        // also, only want direct descendants of halos here
        if ( ((t1-t2) >= tdyn) && ((t1-t2c) < tdyn) && (n1->track == n2->track) ) {

            // if so, we need to interpolate between n2 and n2's child to get approximate value
            double interp_frac = (tdyn - (t1 - t2c)) / (t2c - t2);

            //std::cout << "reached tdyn epoch, t1:" << t1 << ", t2: " << t2 << ", t2c: " << t2c << ", interp_frac: " << interp_frac << std::endl;

            n1->mdot_dyn = (n1->mvir - (n2->child->mvir * (1. - interp_frac) + n2->mvir * interp_frac))/(tdyn*1.e9);
            n1->mdot_dyn_lumpy = (n1->mlumpy - (n2->child->mlumpy * (1. - interp_frac) + n2->mlumpy * interp_frac))/(tdyn*1.e9);
            n1->mdot_dyn_smooth = (n1->msmooth - (n2->child->msmooth * (1. - interp_frac) + n2->msmooth * interp_frac))/(tdyn*1.e9);

            //std::cout << "mvir1: " << n1->mvir << ", mvir2: " << n2->mvir << ", mvir2c: " << n2->child->mvir << ", mdot_dyn: " << n1->mdot_dyn << std::endl;
        }

        // however, if t2 is last ancestor in tree, update halos that have not yet reached tdyn scale with t2
        else if ( (n2->parents.size() == 0) && ((t1-t2) < tdyn) && (n1->track == n2->track) ) {

            //std::cout << "reached end of track, t1:" << t1 << ", t2: " << t2 << std::endl;

            n1->mdot_dyn = (n1->mvir - n2->mvir)/((t1-t2)*1.e9);
            n1->mdot_dyn_lumpy = (n1->mlumpy - n2->mlumpy)/((t1-t2)*1.e9);
            n1->mdot_dyn_smooth = (n1->msmooth - n2->msmooth)/((t1-t2)*1.e9);

            //std::cout << "mvir1: " << n1->mvir << ", mvir2: " << n2->mvir << ", mdot_dyn: " << n1->mdot_dyn << std::endl;
        }

        // let calling function know that this is not correct scale for tdyn
        else {
            return -1;
        }

        // next handle integrated quantities (basically, those that can fluctuate more dramatically)
        n1->lambda_dyn = -1;

        n1->lambda_dyn_disp = -1;

        return 0;
    };

    int dot_freq = int(ceil(forest.size()/100.));
    int counter = 0;

    for (int i = 0; i < forest.size(); i++) {

        if (i % dot_freq == 0) {
            std::cout << ++counter << ".";
            std::cout.flush();
        }

        std::vector<std::vector<node *>> scale;
        
        while(true) {

            std::vector<node *> gather;

            // base case for root node
            if (scale.size() == 0) {
                gather.push_back(forest[i]);
                update_tavg(gather.back(),gather.back());
                scale.push_back(gather);
                continue;
            }

            // gather all parents of last scale
            else {
                for (int j = 0; j < scale.back().size(); j++) {
                    for (int k = 0; k < scale.back()[j]->parents.size(); k++) {
                        gather.push_back(scale.back()[j]->parents[k]);
                    }
                }
            }

            //std::cout << "Found " << gather.size() << " halos at next scale." << std::endl;

            // no more halos left
            if (gather.size() == 0) break;

            // sort gathered halos by parent's track
            std::sort(gather.begin(),gather.end(),[](node * left, node * right){return left->child->track < right->child->track;});
    
            // only update halo if child of next_gen halo is of same track. should be only one match per scale.
            for (int j = 0; j < scale.size(); j++) {

                int k = 0; // scale index
                int l = 0; // gather index

                // adam is vector for containing all nodes from gather that are at end of line and need
                // to update halos of same track which will not ever reach tdyn scale
                std::vector<node *> adam;

                // pgather is vector pointer to switch between gather and adam in cases
                // where adam must be used instead of gather
                std::vector<node *> * pgather;

                pgather = &gather;

                bool first_it = true;

                // step through each scale in parallel with gathered halos
                while (l < pgather->size()) {

                    if (first_it) {

                        first_it = false;

                        // check if this scale will be used
                        double tdyn = t_dyn (scale[j][k]->scale);
                        double t1 = t[list-scale[j][k]->list_count];
                        double t2 = t[list-gather[l]->list_count];
                        double t2c = t[list-gather[l]->child->list_count];

                        if ( ! (((t1-t2) >= tdyn) && ((t1-t2c) < tdyn)) ) {
                            
                            if ((t1-t2) < tdyn) {
                            
                                // find any gathered halos that are final halos (adam) of their track b/c we want to run these.
                                for (int m = 0; m < gather.size(); m++) {
                                    if (gather[m]->parents.size() == 0) {
                                        adam.push_back(gather[m]);
                                    }
                                }

                                // if we didn't find any, then yes progress to next scale
                                if (adam.size() == 0) break;

                                // switch to using adam rather than gather
                                else {
                                    //std::cout << "Switching to adam with size " << adam.size() << std::endl;
                                    pgather = &adam;
                                }
                            }

                            else break;
                        }
                    }

                    //if (adam.size()) std::cout << "k: " << k << ", l: " << l << ", j: " << j << std::endl;

                    if ((*pgather)[l]->child->track > scale[j][k]->track) {
                        if (k+1 < scale[j].size()) k++;
                        else {
                            //std::cout << "incremented k to end without finding match for gather." << std::endl;
                            break;
                        }
                    }

                    // they must be equal, because we should not skip gather halos -- all parents tracks must be present and in order
                    else {

                        //std::cout << "updating tavg" << std::endl;
                        int result = update_tavg(scale[j][k],(*pgather)[l]);

                        l++;
                    }
                }

                // update time avged quantities for end of line halos also
                if (adam.size()) {
                    for (int m = 0; m < adam.size(); m++) {

                        adam[m]->mdot_dyn = adam[m]->child->mdot_dyn;
                        adam[m]->mdot_dyn_lumpy = adam[m]->child->mdot_dyn_lumpy;
                        adam[m]->mdot_dyn_smooth = adam[m]->child->mdot_dyn_smooth;
                    }
                }

            }

            // now sort gathered halos by their own tracks
            std::sort(gather.begin(),gather.end(),[](node * left, node * right){return left->track < right->track;});

            // and add all halos gathered from current scale onto scale vector
            scale.push_back(gather);
        }
    
        // REMOVE
        //break;
    }

    node * np = forest.front();

    std::vector<node *> backp;

    // interactive tree walking.
    // to disable this, set cin_char = 'y' initially
    int cin_int = -1;
    //char cin_char = 0;
    char cin_char = 'y';
    std::string cin_str = "";
    bool quit = false;

    auto isInt = [](std::string & str, int i = 0)->bool {
        if (!str.size()) return false;
        if (!isdigit(str[i]) && (str[i] != '-')) return false;

        char * p;
        strtol(str.substr(i).c_str(),&p,10);
        return (*p == 0);
    };

    auto printHaloProperties = [&mstellar](node * np) {
        std::cout << log10(np->mvir) << ", mpeak " << log10(np->mpeak) << " at scale " << np->a_mpeak << ", track " << np->track << ", almm " << np->almm << ", almm2 " << np->almm2 << ", alminm " << np->alminm << ", almms " << np->almms << ", mstar " << log10(mstellar(np->mpeak)) << ", mlumpy " << log10(np->mlumpy) << ", msmooth " << log10(np->msmooth) << ", merger count " << np->merger_count << ", Avg merger mass " << log10(np->merger_mass_avg) << ", mdot_dyn " << np->mdot_dyn/np->mvir << ", TF " << np->tf << ", TF > 1 since Mpeak? " << bool(np->high_tf) << ", max TF since mpeak " << np->maxtf_peak << " at scale " << np->a_maxtf_peak << ", npt " << np->npt << ", min BSR since mpeak of " << np->min_bsr_peak << " at scale " << np->a_min_bsr_peak << ", min BSR since lmm of " << np->min_bsr_lmm << " at scale " << np->a_min_bsr_lmm << std::endl;
    };

    while (cin_char != 'y' && cin_char != 'Y') {

        while (cin_int < 0) {
            std::cout << "Select tree in forest (0-" << forest.size()-1 << "): ";
            std::cin >> cin_str;
            if (isInt(cin_str)) cin_int = std::stoi(cin_str,nullptr,10);
            if (cin_int < 0 || cin_int > forest.size()-1) {
                std::cout << "Invalid selection." << std::endl;
                cin_int = -1;
            }
        }

        np = forest[cin_int];
        std::cout << "Tree at scale " << np->scale << " with root mass ";
        printHaloProperties (np);
        //<< log10(np->mvir) << ", mpeak " << log10(np->mpeak) << ", track " << np->track << ", almm " << np->almm << ", almms " << np->almms << ", mstar " << mstellar(np->mpeak) << ", mlumpy " << log10(np->mlumpy) << ", msmooth " << log10(np->msmooth) << std::endl;

        bool back = false;
        bool getSubhalo = false;
        bool getParent = false;

        while (true) {
            cin_int = -1;
            while (cin_int < 0 && back == false) {

                std::cout << "Select progenitor (p0-" << np->parents.size()-1 << ")";
                if (np->subhalos.size()) std::cout << " or subhalo (s0-" << (int)np->subhalos.size()-1 << ")";
                std::cout << " or -1 to go back or 'q' to quit: ";
                std::cin >> cin_str;

                // set to other arbitrary negative number (-1 signifies to go back)
                cin_int = -2;

                // now parse cin_str
                if (cin_str.size() > 0) {
                    if (cin_str[0] == 's' || cin_str[0] == 'p') {
                        if (cin_str[0] == 's') getSubhalo = true;
                        else if (cin_str[0] == 'p') getParent = true;
                        if (isInt(cin_str,1)) cin_int = std::stoi(cin_str.substr(1),nullptr,10);
                    }
                    else if (cin_str == "q") {
                        quit = true;
                        break;
                    }
                    else if (cin_str == "-1") cin_int = -1;
                }
                
                if (cin_int < 0 || (getParent && cin_int > np->parents.size()-1) || (getSubhalo && cin_int > np->subhalos.size()-1)) {
                    if (cin_int == -1 && np->child != NULL) {
                        back = true;
                    }
                    else {
                        std::cout << "Invalid selection." << std::endl;
                        cin_int = -1;
                    }
                }
            }
    
            if (quit) {
                quit = false;
                break;
            }

            if (back) {
                //np = np->child;
                if (backp.size()) {
                    np = backp.back();
                    backp.pop_back();
                }
                back = false;
            }
            else {
                backp.push_back(np);
                if (getParent) {
                    np = np->parents[cin_int];
                    getParent = false;
                }
                else if (getSubhalo) {
                    np = np->subhalos[cin_int];
                    getSubhalo = false;
                }
            }
            if (cin_int == -1) std::cout << "Descendent ";
            else std::cout << "Progenitor " << cin_int;
            std::cout << " at scale " << np->scale << ", with mass ";
            printHaloProperties(np);
            // << log10(np->mvir) << ", mpeak " << log10(np->mpeak) << ", track " << np->track << ", almm " << np->almm << ", almms " << np->almms << ", mstar " << mstellar(np->mpeak) << ", mlumpy " << log10(np->mlumpy) << ", msmooth " << log10(np->msmooth) << std::endl;
            if (np->parents.size() == 0) {
                std::cout << "Reached end of merger tree." << std::endl;
                break;
            }
        }

        std::cout << "Quit? (y/n) ";
        std::cin >> cin_char;
        cin_int = -1;
    }

    //std::sort(forest.begin(),forest.end(),[](node * left, node * right){return left->mvir/left->mpeak < right->mvir/right->mpeak;});

    //np = forest.front();

    //for (int i = 0; i < forest.size(); i++) {
    //    std::cout << "mvir: " << log10(forest[i]->mvir) << " " << "mpeak: " << log10(forest[i]->mpeak) << " " << forest[i]->mvir/forest[i]->mpeak << std::endl;
    //}

    //// print random tree branch
    //while (true) {
    //    std::cout << "Scale: " << np->scale << ", mvir: " << log10(np->mvir) << ", " << log10(np->mpeak) << std::endl;
    //    if (np->parents.size() == 0) break;
    //    np = np->parents.front();
    //    //std::cout << np << std::endl;
    //}

    
    // ok, now its time to write these trees into output files.
    // we can either do this in the same format as we read them in, with additional fields tagged onto the end
    // of the already provided fields,
    // or we can use the format I've already established in converMergerTracks, which is what the vizlab
    // code is set up to use.
    std::ofstream outfile;
    
    setOutputFileName();

    if (!outfile_name.size()) outfile_name = output;

    outfile.open(outfile_name.c_str(),std::ofstream::out);

    std::cout << "Writing output file " << outfile_name << std::endl;

    // same number of lines, number of fields increased by number of new fields
    long int header[2] = {num_lines, num_fields+14};

    outfile.write((const char *)header,2*sizeof(long int));

    dot_freq = int(ceil(num_lines/10.));

    //------------------------------------------------------------------------------
    // write output file in binary 2 format with only mmps, starting from most massive
    // to least massive final halos.
    int num_lines_mmp = 0;

    dot_freq = int(ceil(forest.size()/10.));

    for (int i = 0; i < forest.size(); i++) {

        if (i % dot_freq == 0) {
            std::cout << ".";
            std::cout.flush();
        }

        // np = node *, already declared so just reuse.
        np = forest[i];

        // REMOVE
        //if (np->almm == np->almm2) continue;

        while (true) {

            // NOTE: when adding new fields to output, be sure to update header[1] above
            num_lines_mmp++;
            outfile.write((const char *)data[np->i],num_fields*sizeof(double));
            outfile.write((const char *)(&(np->mpeak)),sizeof(double));
            outfile.write((const char *)(&(np->a_mpeak)),sizeof(double));
            //outfile.write((const char *)(&(np->almms)),sizeof(double));
            //outfile.write((const char *)(&(np->mlumpy)),sizeof(double));
            //outfile.write((const char *)(&(np->mlumpy_frac)),sizeof(double));
            //outfile.write((const char *)(&(np->msmooth)),sizeof(double));
            //outfile.write((const char *)(&(np->msmooth_frac)),sizeof(double));
            //outfile.write((const char *)(&(np->merger_count)),sizeof(double));
            //outfile.write((const char *)(&(np->merger_mass_avg)),sizeof(double));
            outfile.write((const char *)(&(np->mdot_dyn)),sizeof(double));
            //outfile.write((const char *)(&(np->mdot_dyn_lumpy)),sizeof(double));
            //outfile.write((const char *)(&(np->mdot_dyn_smooth)),sizeof(double));
            //outfile.write((const char *)(&(np->almm2)),sizeof(double));
            //outfile.write((const char *)(&(np->alminm)),sizeof(double));
            outfile.write((const char *)(&(np->npt)),sizeof(double));
            outfile.write((const char *)(&(np->high_tf)),sizeof(double));
            outfile.write((const char *)(&(np->maxtf_peak)),sizeof(double));
            outfile.write((const char *)(&(np->a_maxtf_peak)),sizeof(double));
            outfile.write((const char *)(&(np->min_bsr_peak)),sizeof(double));
            outfile.write((const char *)(&(np->a_min_bsr_peak)),sizeof(double));
            outfile.write((const char *)(&(np->min_bsr_lmm)),sizeof(double));
            outfile.write((const char *)(&(np->a_min_bsr_lmm)),sizeof(double));

            // change mass_lost_mpeak to fraction of mass lost as a subhalo since mpeak
            if (np->mass_lost_mpeak > 0) np->mass_lost_mpeak = np->mass_lost_as_sh / np->mass_lost_mpeak;

            np->spin_peak = np->lambdap / np->spin_peak;
            np->tu_peak = np->tu / np->tu_peak;

            outfile.write((const char *)(&(np->spin_peak)),sizeof(double));
            outfile.write((const char *)(&(np->tu_peak)),sizeof(double));
            outfile.write((const char *)(&(np->mass_lost_mpeak)),sizeof(double));

            // only write z = 0 halos for now
            //break;

            if (np->parents.size()) np = np->parents[0];
            else break;
        }
    }

    // now update header with correct number of lines
    header[0] = num_lines_mmp;
    outfile.seekp(outfile.beg);
    outfile.write((const char *)header,2*sizeof(long int));

    //-------------------------------------------------------------------------------
    // write output file in binary 2 format with whole trees, in same order as input
    // merger tree file.
    //for (int i = 0; i < vnode.size(); i++) {

    //    if (i % dot_freq == 0) {
    //        std::cout << ".";
    //        std::cout.flush();
    //    }

    //    outfile.write((const char *)data[vnode[i].i],num_fields*sizeof(double));
    //    outfile.write((const char *)(&(vnode[i].mpeak)),sizeof(double));
    //    outfile.write((const char *)(&(vnode[i].almms)),sizeof(double));
    //}
    //-------------------------------------------------------------------------------


    std::cout << " complete." << std::endl;

    outfile.close();

    return 0;
}

//int doConvertMergerTracks (double ** &data, int num_lines, int num_fields) {
//
//    struct track {
//        long int id = -1;
//        long int descid = -1;
//        double mvir = -1;
//        double mar = 0;
//        bool mmp = false;
//        int track_num = 0;
//    };
//
//    // sorting comparator for rank ordering by mvir in track vector
//    struct sort_comp_vtrack {
//        bool operator()(const track &left, const track &right) {
//            return (left.mvir > right.mvir);
//        }
//    };
//
//    std::unordered_map <long int, track> track_map;
//    std::unordered_map<long int, track>::const_iterator fetch;
//
//    std::unordered_map <int, double> mpeak_map;
//    std::unordered_map<int, double>::const_iterator fetch_mpeak;
//
//    // determine what final scale of merger trees is
//    double final_a = data[0][0];
//
//    // we want to reorganize the merger trees by mass (i.e. list most massive tree first, then less massive).
//    // vector to store pair of halo root line number, halo root mass.
//    // this will allow us to look up tree root in data after sorting by mass.
//    std::vector<std::pair<double,long int>> root_mvir;
//
//    int list = 1;
//    double min_a = final_a;
//
//    // for each root halo, save mass and line number
//    for (int i = 0; i < num_lines; i++) {
//        if (data[i][0] == final_a) {
//            root_mvir.push_back(std::make_pair(data[i][10],i));
//        }
//        // list reflects the total number of timesteps/redshifts in the provided merger trees
//        if (data[i][0] < min_a) {
//            min_a = data[i][0];
//            list++;
//        }
//    }
//
//    // now sort by mvir
//    std::sort(root_mvir.begin(),root_mvir.end());
//
//    std::ofstream of;
//    setOutputFileName();
//    of.open(output.c_str(),std::ofstream::out);
//
//    for (int i = root_mvir.size()-1; i >= 0; i--) {
//        std::cout << (long int)data[root_mvir[i].second][1] << ": " << log10(root_mvir[i].first) << " , " << i << std::endl;
//    }
//
//    int tree_count = 1;
//
//    // create datastructure to store xrot, zrot, list_count, tree_count
//    struct extra_params {
//        double xrot = 0;
//        double zrot = 0;
//        long int list_count = 0;
//        long int tree_count = 0;
//        double mpeak = 0;
//    };
//
//    extra_params extras[num_lines];
//
//    // now lets output first tree in halo track format
//    for (int i = root_mvir.size()-1; i >= 0; i--) {
//
//        std::vector<track> vtrack;
//
//        int track_num = 1;
//        int list_count = list;
//        int j = root_mvir[i].second;
//        double curr_a = final_a;
//
//        // loop over all halos in this merger tree
//        while (true) {
//
//            // check if we've reached a new scale or end of file
//            if (j == num_lines || data[j][0] != curr_a) {
//                if (j != num_lines) {
//                    curr_a = data[j][0];    // set new current scale
//                    list_count--;           // decrement to next scale
//                }
//
//                // now we need to go through and set tracks
//
//                // sort by mvir
//                std::sort(vtrack.begin(),vtrack.end(),sort_comp_vtrack());
//
//                // for root halo, just insert right away
//                if (track_num == 1) {
//                    vtrack[0].track_num = track_num++;
//                    track_map[vtrack[0].id] = vtrack[0];
//
//                    mpeak_map[track_num] = vtrack[0].mvir; 
//                }
//
//                // for other halos, we need to check descendent to determine track
//                else {
//
//                    // loop through all halos at this scale and look up descendents in map
//                    for (int k = 0; k < vtrack.size(); k++) {
//                                
//                        fetch = track_map.find(vtrack[k].descid);
//
//                        if (fetch == track_map.end()) {
//                            std::cout << "Error: halo descendent id not found in track map!" << std::endl;
//                            return 0;
//                        }
//
//                        else {
//                            track t = fetch->second;
//                            
//                            // if mmp has already been found, need to start new track
//                            if (t.mmp) {
//
//                                vtrack[k].track_num = track_num++;
//                                track_map[vtrack[k].id] = vtrack[k];
//                            }
//
//                            // this halo is mmp of descendent
//                            else {
//
//                                // update descendent in map 
//                                t.mmp = true;
//                                t.mar = t.mvir-vtrack[k].mvir; // this is not currently implemented correctly.
//                                                               // need to do proper Mdot
//                                track_map[t.id] = t;
//
//                                // insert this halo with same track as desc
//                                vtrack[k].track_num = t.track_num;
//                                track_map[vtrack[k].id] = vtrack[k];
//                            }
//
//                            // update mpeak for all halos at this scale
//                            fetch_mpeak = mpeak_map.find(vtrack[k].track_num);
//
//                            if (fetch_mpeak == mpeak_map.end()) {
//                                mpeak_map[vtrack[k].track_num] = vtrack[k].mvir;
//                            }
//
//                            else {
//                                if (vtrack[k].mvir > fetch_mpeak->second) {
//                                    mpeak_map[vtrack[k].track_num] = vtrack[k].mvir;
//                                }
//                            }
//                        }
//                    }
//                }
//
//                // clear track vector for next scale
//                vtrack.clear();         
//            }
//
//            // include condition to stop when next tree is reached
//            // or when end of file is reached
//            if (j == num_lines) break;
//            if (j > root_mvir[i].second && data[j][0] == final_a) break;
//
//
//            track t;
//            t.id = (long int)data[j][1];
//            t.descid = (long int)data[j][3];
//            t.mvir = data[j][10];
//            t.mmp = false;
//            vtrack.push_back(t);
//
//            // now we need to print additional info like rotation components, track, list
//            double ax = data[j][48]; double ay = data[j][49]; double az = data[j][50];
//            double xrot = 0., zrot = 0., ar;
//            ar = sqrt(ax*ax + ay*ay + az*az);
//            if (ar != 0) {
//                xrot = acos (az/ar) * 180./ M_PI;
//                if (xrot > 90) xrot = -(xrot-90);
//                if (xrot < -90) xrot = -(xrot+90);
//
//                zrot = atan (ay/ax) * 180./ M_PI;
//                if (zrot > 90) zrot = -(zrot-90);
//                if (zrot < -90) zrot = -(zrot+90);
//            }
//
//            extras[j].xrot = xrot;
//            extras[j].zrot = zrot;
//            extras[j].list_count = list_count;
//            extras[j].tree_count = tree_count;
//
//            j++;
//        }
//
//        // save end position, then reiterate through all halos in tree
//        int j_end = j;
//
//        // for each halo in tree, update mpeak value
//        for (j = root_mvir[i].second, j < j_end; j++) {
//
//            // find appropriate track num for halo
//            fetch = track_map.find(data[j][1]);
//            track t;
//
//            if (fetch == track_map.end()) {
//                std::cout << "Error: halo id not found in track map!" << std::endl;
//                return 0;
//            }
//
//            else {
//                t = fetch->second;
//            }
//
//            // look up corresponding mpeak for that track
//            fetch_mpeak = mpeak_map.find(t.track_num);
//
//            if (fetch_mpeak == mpeak_map.end()) {
//                std::cout << "Error: track num not found in mpeak map!" << std::endl;
//                return 0;
//            }
//
//            // save in extra params for that halo
//            else {
//                extras[j].mpeak = fetch_mpeak->second;
//            }
//        }
//
//        // reset j (unnecessary?)
//        j = j_end;
//
//        // clear map for next tree, since track nums will be reused
//        mpeak_map.clear();
//
//        tree_count++;
//    }
//
//    // loop over trees until all scales have been output
//    for (int l = list; l > 0; l--) {
//
//        // loop over all trees and output scale l
//        for (int i = root_mvir.size()-1; i >= 0; i--) {
//
//            int j = root_mvir[i].second;
//
//            // if this tree has finished, skip
//            if (j == -1) continue;
//
//            double curr_a = data[j][0];
//
//            // loop over all halos at scale l in this tree
//            while (true) {
//
//                // include condition to stop when end of data is reached
//                if (j == num_lines) {
//                    root_mvir[i].second = -1;
//                    break;
//                }
//
//                // include condition to stop when next scale is reached            
//                if (j > root_mvir[i].second && data[j][0] == final_a) {
//                    root_mvir[i].second = -1;
//                    break;
//                }
//
//                // stop when next scale is reached
//                if (curr_a != data[j][0]) {
//                    root_mvir[i].second = j;
//                    break;
//                }
//
//                // loop over fields in each halo
//                for (int k = 0; k < num_fields; k++) {
//                    switch (k) {
//                        case 0:     // scale
//                        case 10:    // mvir
//                        case 11:    // rvir
//                        case 12:    // rs
//                        case 15:    // almm
//                        case 16:    // vmax
//                        case 17:    // x
//                        case 18:    // y
//                        case 19:    // z
//                        case 20:    // vx
//                        case 21:    // vy
//                        case 22:    // vz
//                        case 26:    // lambda
//                        case 35:    // tf
//                        case 37:    // rs_klypin
//                        case 43:    // xoff
//                        case 44:    // voff
//                        case 45:    // lambadp
//                        case 46:    // b/a
//                        case 47:    // c/a
//                        case 51:    // b/a500
//                        case 52:    // c/a500
//                        case 56:    // T/|U|
//                            of << data[j][k] << " ";
//                            break;
//                        case 1:     // id
//                        case 3:     // descid
//                            of << (long int)data[j][k] << " ";
//                            break;
//                        default:
//                            break;
//                    }
//                }
//
//                fetch = track_map.find(data[j][1]);
//
//                track t;
//
//                if (fetch == track_map.end()) {
//                    std::cout << "Error: halo id not found in track map!" << std::endl;
//                    return 0;
//                }
//
//                else {
//                    t = fetch->second;
//                }
//
//                of << extras[j].xrot << " " << extras[j].zrot << " " << extras[j].list_count << " ";
//                of << extras[j].tree_count << " " << t.track_num << " " << t.mar << std::endl;
//
//                j++;
//            }
//        }
//    }
//            
//    return 0;
//}

// assume data is read in as bin2 mmp tree
int doIndvHaloTracks (double ** &data, int num_lines, int num_fields, int x) {

    double final_a = data[0][0];

    std::ofstream of;

    setOutputFileName();

    int of_count = 0;

    std::cout << "Writing individual halo track files... " << std::endl;
    std::cout << output.substr(0,output.find_last_of('.'))+"_"+str(0)+"-"+(x>0?str(x-1):"*")+output.substr(output.find_last_of('.')) << std::endl;

    // want to write ascii data.
    // write x files, unless x < num_trees
    for (int i = 0; i < num_lines; i++) {

        if (data[i][0] == final_a) {

            std::cout << of_count << ".";
            std::cout.flush();

            if (of_count > 0) {
                of.close();
            }

            // default is x = 0, which means print all tracks in file
            if (x > 0 && of_count > x) break;

            int ind = output.find_last_of('.');
            std::string new_output = output.substr(0,ind)+"_"+str(of_count++)+output.substr(ind);
            of.open(new_output.c_str(),std::ofstream::out);
        }

        for (int j = 0; j < num_fields-1; j++) {

            if (j == 1) of << (unsigned int)data[i][j] << " ";
            else of << data[i][j] << " ";
        }

        of << data[i][num_fields-1] << std::endl;

        if (i == num_lines-1) {
            of.close();
        }
    }

    std::cout << "done." << std::endl;

    return 0;
}
