#include "web.h"

int doCombineSpineCatalog (double ** &data, int num_lines, int num_fields, bool binary_write) {

    // lets read in a spine catalog, say the 4mpc smoothed walls
    std::cout << "Starting spine catalog reading..." << std::endl;
    std::cout.flush();
    int *** wall4;
    int dim = 1024;
    std::string src = "/home/maragon/SPINE_0.1/SPINE/";
    std::string wall4_fn = "0273.GAU10.Y-0-1-2-3.iwall";
    doBinaryRead3D (src+wall4_fn, wall4, dim);

    std::unordered_map <int, std::vector <int *> > wall4_map;
    int id = 0;
    int dot_freq = int(ceil(dim/10.));

    for (int i = 0; i < dim; i++) {
        if (i % dot_freq == 0) {
            std::cout << ".";
            std::cout.flush();
        }
        for (int j = 0; j < dim; j++) {
            for (int k = 0; k < dim; k++) {
                id = wall4[i][j][k];
                if (!wall4_map.count(id)) wall4_map.insert(std::make_pair(id, std::vector<int *>()));
                int * pos = new int[3];
                pos[0] = i; pos[1] = j; pos[2] = k;
                wall4_map.at(id).push_back(pos);
            }
        }
    }

    std::cout << "Number of walls: " << wall4_map.size() << std::endl;

    std::ofstream of;
    setOutputFileName();
    of.open(output.c_str(),std::ofstream::out);

    for (auto it = wall4_map.begin(); it != wall4_map.end(); it++) {
        of << it->first << "\t" << it->second.size();
        // avg position
        int pos[3] = {0,0,0};
        for (auto vit = it->second.begin(); vit != it->second.end(); vit++) {
            for (int i = 0; i < 3; i++) pos[i] += (*vit)[i];
        }
        for (int i = 0; i < 3; i++) pos[i]/=it->second.size();

        of << " " << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
        
    }

    of.close();

    return -1;
}

int doSpineCatalogAnalysis (double ** &data, int num_lines, int num_fields, int x) {

    std::cout << "Reading spine catalog..." << std::endl;
    std::cout.flush();
    int *** spine_data;
    int dim = 1024;
    std::string src = "/zang/chtlee/pfs_backup/miguelac/";
    std::string file_name;
    
    switch (x) {
        case 0: file_name = "Particles.0273.0000.GAU-2.0.Y-0-1-2.ifila"; break;
        case 1: file_name = "Particles.0273.0000.GAU-2.0.Y-0-1-2.iwall"; break;
        case 2: file_name = "Particles.0273.0000.GAU-2.0.Y-0-1.ifila"; break;
        case 3: file_name = "Particles.0273.0000.GAU-2.0.Y-0-1.iwall"; break;
        case 4: file_name = "Particles.0273.0000.GAU-2.0.Y-0.ifila"; break;
        case 5: file_name = "Particles.0273.0000.GAU-2.0.Y-0.iwall"; break;
    }

    doBinaryRead3D (src+file_name, spine_data, dim, 256);

//-----------------------------------------------------------
//// minimal code to read spine catalog for elliot.
//
//// file name of input file
//std::string fileName = "spine.fila.0.0-128.bin";
//
//// dimension of input file (e.g. 128^3 or 1024^3)
//dim = 128;
//std::cout << "Reading binary file: \"" << fileName << "\"" << std::endl;
//
//// open input file stream (may need <iostream> <fstream>)
//std::ifstream inpt(fileName,std::ios::in);
//
//// create variable to store data. allocate dynamically so
//// that its easier to change dimensions if desired.
//int *** data2 = (int ***) malloc (dim * sizeof (int **));
//
//// keeps track of progress reading file
//int dot_freq2 = int(ceil(dim/10.));
//
//// allocate space for data2 and read in one dimension at a time
//for (int i = 0; i < dim; i++) {
//    if (i % dot_freq2 == 0) {
//        std::cout << ".";
//        std::cout.flush();
//    }
//
//    data2[i] = (int **)malloc(dim * sizeof(int *));
//
//    for (int j = 0; j < dim; j++) {
//
//        data2[i][j] = (int *)malloc(dim * sizeof(int));
//        memset(data2[i][j],0,dim * sizeof(int));
//
//        inpt.read((char *)data2[i][j],dim * sizeof(int));
//        if (!inpt.good()) {
//            std::cerr << "Stream bad after reading i = " << i << ", j = " << j << std::endl;
//            for (int k = 0; k < dim; k++) std::cout << data2[i][j][k] << "  ";
//            std::cout << std::endl;
//            return 0;
//        }
//    }   
//}
//
//inpt.close();
//
//std::cout << " complete." << std::endl;
//
//int count = 0;
//
//for (int i = 0; i < dim; i++) {
//    for (int j = 0; j < dim; j++) {
//        for (int k = 0; k < dim; k++) {
//            if (data2[i][j][k] != -1) count++;
//        }
//    }
//}
//
//std::cout << "Found " << count << " filament voxels (" << count/pow(128.,3.)*100 << "%)" << std::endl;
//
//return 0;
//------------------------------------------------------

    struct loc {
        double x;
        double y;
        double z;
    };

    struct web_substructure {

        int sector_id = -1;

        std::vector<loc> voxels;

        loc min;
        loc max;
    };

    struct web_element {

        // structure id
        int id = -1;

        // to distinguish between walls, filaments, voids, etc
        int type = 0;

        // box sector(s) -- points to web substructures within each sector
        std::vector<web_substructure> sectors;

        // member particles/voxels
        std::vector<loc> voxels;

        // bounds
        loc min;
        loc max;

        // position
        loc med_center;
        loc avg_center;

        double volume;
        double radius;  // used for walls
        double length;  // used for filaments
    };

    std::unordered_map <int, web_element> web_map;
    std::unordered_map<int, web_element>::const_iterator fetch;

    double cell_size = 250./dim;
    int sectors_per_side = 64;
    int sectors_per_slice = sectors_per_side*sectors_per_side;

    std::vector<web_element> web_v;
    loc l;

    // create sector vector. each element is a vector containing pointers to web_elements in that sector.
    // currently using 125 sectors (five per side).
    std::vector<std::vector<web_element *>> all_sectors;
    for (int i = 0; i < (int)pow(sectors_per_side,3.); i++) {
        std::vector<web_element *> init;
        all_sectors.push_back(init);
    }

    // vector of web pointers -- better for doing custom sorting, etc. without messing up
    // pointers to individual web structrues
    std::vector<web_element *> pweb_v;

    std::cout << "Building web data structures and collecting voxels..." << std::endl;

    int dot_freq = int(ceil(dim/10.));

    for (int i = 0; i < dim; i++) {
        if (i % dot_freq == 0) {
            std::cout << ".";
            std::cout.flush();
        }
        for (int j = 0; j < dim; j++) {
            for (int k = 0; k < dim; k++) {
                if (spine_data[i][j][k] > -1) {

                    int id = spine_data[i][j][k];
                    while (id >= web_v.size()) {
                        web_element new_element;
                        new_element.id = web_v.size();
                        web_v.push_back(new_element);
                    }

                    // reversing the ordering here, so that x = k, z = i, due
                    // to the way the cube was stored
                    l.x = k; l.y = j; l.z = i;

                    web_v[id].voxels.push_back(l);
                }
            }
        }
    }

    std::cout << " complete." << std::endl;

//    std::sort(web_v.begin(),web_v.end(),[](web_element left, web_element right){return left.id < right.id;});
    auto vmin = [](std::vector<loc> &v)->loc {
        loc ret;
        ret.x = v[0].x; ret.y = v[0].y; ret.z = v[0].z;
        for (int i = 1; i < v.size(); i++) {
            if (v[i].x < ret.x) ret.x = v[i].x;
            if (v[i].y < ret.y) ret.y = v[i].y;
            if (v[i].z < ret.z) ret.z = v[i].z;
        }
        return ret;
    };

    auto vmax = [](std::vector<loc> &v)->loc {
        loc ret;
        ret.x = v[0].x; ret.y = v[0].y; ret.z = v[0].z;
        for (int i = 1; i < v.size(); i++) {
            if (v[i].x > ret.x) ret.x = v[i].x;
            if (v[i].y > ret.y) ret.y = v[i].y;
            if (v[i].z > ret.z) ret.z = v[i].z;
        }
        return ret;
    };

    auto vavg = [](std::vector<loc> &v)->loc {
        loc ret;
        for (int i = 0; i < v.size(); i++) {
            ret.x += v[i].x;
            ret.y += v[i].y;
            ret.z += v[i].z;
        }
        ret.x = ceil(((double)ret.x)/v.size()-0.5);
        ret.y = ceil(((double)ret.y)/v.size()-0.5);
        ret.z = ceil(((double)ret.z)/v.size()-0.5);
        return ret;
    };

    auto vmed = [](std::vector<loc> &v)->loc {
        loc ret;
        int len = v.size();
        std::sort(v.begin(),v.end(),[](loc left,loc right){return left.x < right.x;});
        if (len % 2) ret.x = v[(len-1)/2].x;
        else ret.x = 0.5*(v[(len-2)/2].x + v[(len)/2].x);
        std::sort(v.begin(),v.end(),[](loc left,loc right){return left.y < right.y;});
        if (len % 2) ret.y = v[(len-1)/2].y;
        else ret.y = 0.5*(v[(len-2)/2].y + v[(len)/2].y);
        std::sort(v.begin(),v.end(),[](loc left,loc right){return left.z < right.z;});
        if (len % 2) ret.z = v[(len-1)/2].z;
        else ret.z = 0.5*(v[(len-2)/2].z + v[(len)/2].z);
        return ret;
    };

    struct sector_voxel {
        int sector_id = -1;
        loc voxel;
    };

    std::cout << "Computing structure properties... " << std::endl;

    dot_freq = int(ceil(web_v.size()/100.));
    int counter = 0;

    for (int i = 0; i < web_v.size(); i++) {
        if (i % dot_freq == 0) {
            std::cout << ++counter << ".";
            std::cout.flush();
        }
        if (web_v[i].voxels.size() == 0) continue;

        web_v[i].min = vmin(web_v[i].voxels);
        web_v[i].max = vmax(web_v[i].voxels);

        web_v[i].avg_center = vavg(web_v[i].voxels);

        // compute sector membership before correcting for periodic boundary

        // sector indeces in each dimension
        int nx, ny, nz;

        // vector of sector for each voxel in structure
        std::vector<sector_voxel> vsectors(web_v[i].voxels.size());

        // determine sector of each voxel in structure
        for (int j = 0; j < web_v[i].voxels.size(); j++) {

            nx = floor(web_v[i].voxels[j].x/(double)(dim/sectors_per_side));
            ny = floor(web_v[i].voxels[j].y/(double)(dim/sectors_per_side));
            nz = floor(web_v[i].voxels[j].z/(double)(dim/sectors_per_side));
            vsectors[j].sector_id = nx+sectors_per_side*ny+sectors_per_slice*nz;

            // adjust voxel positions to centers of voxels and in units of mpc/h
            vsectors[j].voxel.x = 250./dim*(web_v[i].voxels[j].x+0.5);
            vsectors[j].voxel.y = 250./dim*(web_v[i].voxels[j].y+0.5);
            vsectors[j].voxel.z = 250./dim*(web_v[i].voxels[j].z+0.5);
        }

        // sort vsectors so that we can build sector-voxel substructure
        std::sort(vsectors.begin(),vsectors.end(),[](sector_voxel left, sector_voxel right){return left.sector_id < right.sector_id;});

        for (int j = 0; j < vsectors.size(); j++) {
            if ((j == 0) || (vsectors[j].sector_id != vsectors[j-1].sector_id)) {

                web_substructure ws;
                ws.sector_id = vsectors[j].sector_id;
                ws.voxels.push_back(vsectors[j].voxel);

                web_v[i].sectors.push_back(ws);
            }
            else {
                web_v[i].sectors.back().voxels.push_back(vsectors[j].voxel);
            }
        }

        // now check for cases that cross periodic boundary of simulation.
        // if that happens, adjust internal positions to be +-512 by subtracting 1024 from values >= 512

        int x_adjust = 0, y_adjust = 0, z_adjust = 0;
        if (web_v[i].max.x - web_v[i].min.x == dim-1) {
            if (web_v[i].avg_center.x >= dim/2) x_adjust = dim;
            else x_adjust = -dim;
        }
        if (web_v[i].max.y - web_v[i].min.y == dim-1) {
            if (web_v[i].avg_center.y >= dim/2) y_adjust = dim;
            else y_adjust = -dim;
        }
        if (web_v[i].max.z - web_v[i].min.z == dim-1) {
            if (web_v[i].avg_center.z >= dim/2) z_adjust = dim;
            else z_adjust = -dim;
        }

        if (x_adjust || y_adjust || z_adjust) {

            for (int j = 0; j < web_v[i].voxels.size(); j++) {
                if ((web_v[i].voxels[j].x >= dim/2 && x_adjust < 0) ||
                    (web_v[i].voxels[j].x < dim/2 && x_adjust > 0)) web_v[i].voxels[j].x += x_adjust;
                if ((web_v[i].voxels[j].y >= dim/2 && y_adjust < 0) ||
                    (web_v[i].voxels[j].y < dim/2 && y_adjust > 0)) web_v[i].voxels[j].y += y_adjust;
                if ((web_v[i].voxels[j].z >= dim/2 && z_adjust < 0) ||
                    (web_v[i].voxels[j].z < dim/2 && z_adjust > 0)) web_v[i].voxels[j].z += z_adjust;
            }

            // now recompute avg position
            web_v[i].avg_center = vavg(web_v[i].voxels);

            web_v[i].min = vmin(web_v[i].voxels);
            web_v[i].max = vmax(web_v[i].voxels);
        }

        web_v[i].med_center = vmed(web_v[i].voxels);

        // volume of min/max bounds
        web_v[i].volume = pow(cell_size,3.)*(web_v[i].max.x-web_v[i].min.x)*(web_v[i].max.y-web_v[i].min.y)*(web_v[i].max.z-web_v[i].min.z);

        // effective radius determined assuming circular distribution of area
        web_v[i].radius = cell_size*sqrt(web_v[i].voxels.size())/M_PI;

        // also compute min/max for each sector substructure
        for (int j = 0; j < web_v[i].sectors.size(); j++) {
            web_v[i].sectors[j].min = vmin(web_v[i].sectors[j].voxels);
            web_v[i].sectors[j].max = vmax(web_v[i].sectors[j].voxels);
        }
    }

    std::cout << " complete." << std::endl;

    for (int i = 0; i < web_v.size(); i++) {
        pweb_v.push_back(&(web_v[i]));
        for (int j = 0; j < web_v[i].sectors.size(); j++) {
            all_sectors[web_v[i].sectors[j].sector_id].push_back(pweb_v.back());
        }
    }

    std::sort(pweb_v.begin(),pweb_v.end(),[](web_element * left, web_element * right){return left->volume < right->volume;});

    std::cout << "Min volume: " << pweb_v.front()->volume << ", max volume: " << pweb_v.back()->volume << std::endl;

    std::sort(pweb_v.begin(),pweb_v.end(),[](web_element * left, web_element * right){return left->radius < right->radius;});

    std::cout << "Min radius: " << pweb_v.front()->radius << ", max radius: " << pweb_v.back()->radius << std::endl;

    int cin_int = -1;
    //char cin_char = 0;
    char cin_char = 'y';
    std::string cin_str = "";
    bool quit = false;

    auto isInt = [](std::string & str)->bool {
        if (!str.size()) return false;
        if (!isdigit(str[0]) && (str[0] != '-')) return false;

        char * p;
        strtol(str.c_str(),&p,10);
        return (*p == 0);
    };

    while (cin_char != 'y' && cin_char != 'Y') {

        while (cin_int < 0) {
            std::cout << "Select: structures (a), sectors (b): ";
            std::cin >> cin_str;
            cin_int = 0;
            if (cin_str != "a" && cin_str != "b") {
                std::cout << "Invalid selection." << std::endl;
                cin_int = -1;
            }
        }

        cin_int = -1;

        if (cin_str == "a") {

            while (cin_int < 0) {
                std::cout << "Select structure id (0-" << web_v.size()-1 << "): ";
                std::cin >> cin_str;
                if (isInt(cin_str)) cin_int = std::stoi(cin_str,nullptr,10);
                if (cin_int < 0 || cin_int > web_v.size()-1) {
                    std::cout << "Invalid selection." << std::endl;
                    cin_int = -1;
                }
            }

            std::cout << "Structure " << cin_int << " has " << web_v[cin_int].voxels.size() << " member voxels." << std::endl;
            std::cout << "bounds: (" << web_v[cin_int].min.x << "," << web_v[cin_int].min.y << "," << web_v[cin_int].min.z << ") - (" << web_v[cin_int].max.x << "," << web_v[cin_int].max.y << "," << web_v[cin_int].max.z << ")" << std::endl;
            std::cout << "avg_center: (" << web_v[cin_int].avg_center.x << "," << web_v[cin_int].avg_center.y << "," << web_v[cin_int].avg_center.z << ")" << std::endl;
            std::cout << "med_center: (" << web_v[cin_int].med_center.x << "," << web_v[cin_int].med_center.y << "," << web_v[cin_int].med_center.z << ")" << std::endl;
            std::cout << "volume (mpc/h)^3: " << web_v[cin_int].volume << ", radius (mpc/h): " << web_v[cin_int].radius << std::endl;
            std::cout << "Sector membership:";
            for (int i = 0; i < web_v[cin_int].sectors.size(); i++) {
                std::cout << " " << web_v[cin_int].sectors[i].sector_id;
                std::cout << " ("<<web_v[cin_int].sectors[i].voxels.size()/(double)web_v[cin_int].voxels.size()*100.<<"%) ";
            }
            std::cout << std::endl;
        }

        else if (cin_str == "b") {

            while (cin_int < 0) {
                std::cout << "Select sector id (0-" << all_sectors.size()-1 << "): ";
                std::cin >> cin_str;
                if (isInt(cin_str)) cin_int = std::stoi(cin_str,nullptr,10);
                if (cin_int < 0 || cin_int > all_sectors.size()-1) {
                    std::cout << "Invalid selection." << std::endl;
                    cin_int = -1;
                }
            }

            std::cout << "Sector " << cin_int << " has " << all_sectors[cin_int].size() << " member structures." << std::endl;
            std::cout << "Structure IDs: ";
            for (int i = 0; i < all_sectors[cin_int].size(); i++) std::cout << " " << all_sectors[cin_int][i]->id;
            std::cout << std::endl;
        }

        std::cout << "Quit? (y/n) ";
        std::cin >> cin_char;
        cin_int = -1;
    }

    int nx, ny, nz;
    double hx, hy, hz, hx_adjust, hy_adjust, hz_adjust;
    double sector_size = 250./(double)sectors_per_side, dist = 0.;//, dilation;
    int halo_sector, curr_sector, structures_checked, all_structures;
    long int voxels_checked = 0;

    std::vector<double> min_dist(num_lines);
    std::vector<int> id_dist(num_lines);
    std::vector<loc> voxel_dist(num_lines);

    std::cout << "Now computing distance from halos to nearest structure... " << std::endl;

    // build vector of sector offsets to iterate over
    std::vector<loc> sector_offset(729);

    sector_offset[0].x = 0;
    sector_offset[0].y = 0;
    sector_offset[0].z = 0;
    int so_count = 1;

    for (int sox = -1; sox < 2; sox++) {
        for (int soy = -1; soy < 2; soy++) {
            for (int soz = -1; soz < 2; soz++) {
                if (sox != 0 || soy != 0 || soz != 0) {
                    sector_offset[so_count].x = sox;
                    sector_offset[so_count].y = soy;
                    sector_offset[so_count].z = soz;
                    so_count++;
                }
            }
        }
    }

    for (int sox = -2; sox < 3; sox++) {
        for (int soy = -2; soy < 3; soy++) {
            for (int soz = -2; soz < 3; soz++) {
                if (sox*sox==4 || soy*soy==4 || soz*soz==4) {
                    sector_offset[so_count].x = sox;
                    sector_offset[so_count].y = soy;
                    sector_offset[so_count].z = soz;
                    so_count++;
                }
            }
        }
    }

    for (int sox = -3; sox < 4; sox++) {
        for (int soy = -3; soy < 4; soy++) {
            for (int soz = -3; soz < 4; soz++) {
                if (sox*sox==9 || soy*soy==9 || soz*soz==9) {
                    sector_offset[so_count].x = sox;
                    sector_offset[so_count].y = soy;
                    sector_offset[so_count].z = soz;
                    so_count++;
                }
            }
        }
    }

    for (int sox = -4; sox < 5; sox++) {
        for (int soy = -4; soy < 5; soy++) {
            for (int soz = -4; soz < 5; soz++) {
                if (sox*sox==16 || soy*soy==16 || soz*soz==16) {
                    sector_offset[so_count].x = sox;
                    sector_offset[so_count].y = soy;
                    sector_offset[so_count].z = soz;
                    so_count++;
                }
            }
        }
    }

    //std::cout << "final sector_offset count: "<<so_count-1<<std::endl;

    // number of halos with no nearby structure found
    int not_found = 0;

    int level_four_count = 0, level_three_count = 0, level_two_count = 0;

    dot_freq = int(ceil(num_lines/100.));
    counter = 0;

    // now, begin the process of computing nearest web structure to halos in input catalog
    for (int i = 0; i < num_lines; i++) {

        if (i % dot_freq == 0) {
            std::cout << ++counter << ".";
            std::cout.flush();
        }

        // number of structures we've had to check for this halo
        structures_checked = 0;
        all_structures = 0;
        voxels_checked = 0;

        // initialize distance and id vectors
        min_dist[i] = -1;
        id_dist[i] = -1;

        //dilation = 4*250./(double)sectors_per_side;

        hx = data[BP_CAT.x][i];
        hy = data[BP_CAT.y][i];
        hz = data[BP_CAT.z][i];

        // compute halo sector
        nx = floor(hx/sector_size);
        ny = floor(hy/sector_size);
        nz = floor(hz/sector_size);

        //hx = (int)floor(hx*dim/250.);
        //hy = (int)floor(hy*dim/250.);
        //hz = (int)floor(hz*dim/250.);

        halo_sector = nx+sectors_per_side*ny+sectors_per_slice*nz;

        //std::cout << "Halo " << (long int)data[BP_CAT.id][i] << " at ("<<hx<<","<<hy<<","<<hz<<") in sector " << halo_sector << std::endl;

        // now iterate over x,y,z += -1,0,1
        //for (int sx = -1; sx < 2; sx++) {
        //    for (int sy = -1; sy < 2; sy++) {
        //        for (int sz = -1; sz < 2; sz++) {

        int sx, sy, sz;
        bool skip_level_two = false, skip_level_three = false, skip_level_four = false;
        for (int offset = 0; offset < sector_offset.size(); offset++) {

            // if we find a structure closer than the sector size in 0th or 1st level, then no need to
            // check sectors on the second level (e.g. 0,0,2), since they must be farther
            if (min_dist[i] > -1) {
                if (offset < 27 && min_dist[i] < 250./sector_size) skip_level_two = true;
                else if (offset < 125 && min_dist[i] < 2*250./sector_size) skip_level_three = true;
                else if (offset < 343 && min_dist[i] < 3*250./sector_size) skip_level_four = true;
            }

            if (skip_level_two && offset >= 27) break;
            else if (skip_level_three && offset >= 125) break;
            else if (skip_level_four && offset >= 343) break;


            if (offset == 27) level_two_count++;
            if (offset == 125) level_three_count++;
            if (offset == 343) level_four_count++;

            // adjust halo position to be consistent with current sector, if necessary.
            // remember, substructure voxels are not corrected already.
            hx_adjust = hx;
            hy_adjust = hy;
            hz_adjust = hz;

            sx = sector_offset[offset].x;
            sy = sector_offset[offset].y;
            sz = sector_offset[offset].z;

            curr_sector = halo_sector;

            // handle edge cases appropriately
            if (nx + sx < 0) {
                curr_sector += (sectors_per_side+sx);
                hx_adjust = hx + 250.;
            }
            else if (nx + sx > (sectors_per_side-1)) {
                curr_sector += -(sectors_per_side-sx);
                hx_adjust = hx - 250.;
            }
            else curr_sector += sx;

            if (ny + sy < 0) {
                curr_sector += (sectors_per_side+sy)*sectors_per_side;
                hy_adjust = hy + 250.;
            }
            else if (ny + sy > (sectors_per_side-1)) {
                curr_sector += -(sectors_per_side-sy)*sectors_per_side;
                hy_adjust = hy - 250.;
            }
            else curr_sector += sy*sectors_per_side;

            if (nz + sz < 0) {
                curr_sector += (sectors_per_side+sz)*sectors_per_slice;
                hz_adjust = hz + 250.;
            }
            else if (nz + sz > (sectors_per_side-1)) {
                curr_sector += -(sectors_per_side-sz)*sectors_per_slice;
                hz_adjust = hz - 250.;
            }
            else curr_sector += sz*sectors_per_slice;

            //if (curr_sector > sectors_per_side*sectors_per_slice) {
            //    std::cout << "curr_sector: " << curr_sector<<std::endl;
            //    std::cout << sx << " "<<sy<<" "<<sz<<std::endl;
            //    std::cout << nx<<" "<<ny<<" "<<nz<<std::endl;
            //}

            //std::cout << "\tchecking nearby sector ("<<sx<<","<<sy<<","<<sz<<") "<<curr_sector;
            //std::cout << " with " << all_sectors[curr_sector].size() << " member structures."<<std::endl;

            all_structures += all_sectors[curr_sector].size();

            // now eliminate structures with incompatible bounds
            for (int j = 0; j < all_sectors[curr_sector].size(); j++) {

                web_substructure * ws;

                // find componenet of structure that exists in this sector
                for (int k = 0; k < all_sectors[curr_sector][j]->sectors.size(); k++) {
            
                    if (all_sectors[curr_sector][j]->sectors[k].sector_id == curr_sector) {

                        ws = &(all_sectors[curr_sector][j]->sectors[k]);
                        break;
                    }
                }

                // first correct bounds across periodic boundary to match up with structure, if necessary
                //if (hx-ws->med_center.x > 512) hx -= 1024;
                //else if (ws->med_center.x-hx > 512) hx += 1024;

                //if (hy-ws->med_center.y > 512) hy -= 1024;
                //else if (ws->med_center.y-hy > 512) hy += 1024;

                //if (hz-ws->med_center.z > 512) hz -= 1024;
                //else if (ws->med_center.z-hz > 512) hz += 1024;

                // now adjust bounds checking to be limited by closest structure found so far
                //if (min_dist[i] > -1) dilation = min_dist[i];

                // now just check if bounds are ok, otherwise skip.
                // check that bounds (dilated by 16 voxels) contain halo location
                if (min_dist[i] > -1) {
                    if (ws->min.x-hx_adjust > min_dist[i] ||
                        hx_adjust-ws->max.x > min_dist[i]) continue;

                    if (ws->min.y-hy_adjust > min_dist[i] ||
                        hy_adjust-ws->max.y > min_dist[i]) continue;

                    if (ws->min.z-hz_adjust > min_dist[i] ||
                        hz_adjust-ws->max.z > min_dist[i]) continue;
                }

                //std::cout << "\t\tstructure "<<all_sectors[curr_sector][j]->id<<" has acceptable bounds"<<std::endl;

                structures_checked++;
                voxels_checked += ws->voxels.size();

                // so now we've checked that this structure has bounds potentially compatible with the halo.
                // now need to iterate through voxels and compute distances.
                for (int k = 0; k < ws->voxels.size(); k++) {

                    // compute distance
                    dist = pow(ws->voxels[k].x-hx_adjust,2.);
                    dist += pow(ws->voxels[k].y-hy_adjust,2.);
                    dist += pow(ws->voxels[k].z-hz_adjust,2.);
                    dist = sqrt(dist);

                    if (min_dist[i] == -1) {
                        min_dist[i] = dist;
                        id_dist[i] = all_sectors[curr_sector][j]->id;
                        voxel_dist[i] = ws->voxels[k];
                        //std::cout << "\t\tInitialized min distance ("<<dist<<") to voxel ("<<ws->voxels[k].x<<","<<ws->voxels[k].y<<","<<ws->voxels[k].z<<") in structure "<<all_sectors[curr_sector][j]->id<<std::endl;
                    }

                    else if (min_dist[i] > -1 && dist <= min_dist[i]) {
                        //std::cout << "\t\tcomparable or smaller distance ("<<dist<<") to voxel ("<<ws->voxels[k].x<<","<<ws->voxels[k].y<<","<<ws->voxels[k].z<<") in structure "<<all_sectors[curr_sector][j]->id<<std::endl;
                        if (dist == min_dist[i]) {
                            // might want to update this to just compare exact min distance instead of structure size
                            if (ws->voxels.size() > web_v[id_dist[i]].voxels.size()) {
                                id_dist[i] = all_sectors[curr_sector][j]->id;
                                voxel_dist[i] = ws->voxels[k];
                            }
                        }
                        else {
                            min_dist[i] = dist;
                            id_dist[i] = all_sectors[curr_sector][j]->id;
                            voxel_dist[i].x = ws->voxels[k].x;
                            voxel_dist[i].y = ws->voxels[k].y;
                            voxel_dist[i].z = ws->voxels[k].z;
                        }
                    }
                }
            }
        }

        //if (i%1000==0) {
        //    std::cout << "("<<i<<") checked " << structures_checked << " out of "<<all_structures<<" complete structures and found min distance "<<min_dist[i]<<" to structure "<<id_dist[i]<<std::endl;
        //    std::cout << "\tcomputed distances to "<<voxels_checked<<" voxels ("<<voxels_checked/(pow(dim/(double)sectors_per_side,3)*125)*100<<"%)"<<std::endl;
        //}

        //if (min_dist[i] == -1) {
        //    not_found++;
        //    std::cout << "Halo " << (long int)data[BP_CAT.id][i] << " not near any structures.  Voxels checked " << voxels_checked << std::endl;
        //}
    }

    std::cout << " complete." << std::endl;

    std::cout << "Could not find nearby structures for " << not_found << " halos." << std::endl;
    std::cout << "Had to search on level four for " << level_four_count << " halos." << std::endl;
    std::cout << "Had to search on level three for " << level_three_count << " halos." << std::endl;
    std::cout << "Had to search on level two for " << level_two_count << " halos." << std::endl;

    // now lets output some of these results
    std::ofstream outfile;

    setOutputFileName();

    outfile.open(output.c_str(),std::ofstream::out);

    std::cout << "Writing output file \"" << output << "\" ";

    dot_freq = int(ceil(num_lines/10.));

    for (int i = 0; i < num_lines; i++) {
        if (i % dot_freq == 0) {
            std::cout << ".";
            std::cout.flush();
        }
        outfile << (long int)data[BP_CAT.id][i] << " " << data[BP_CAT.x][i] << " " << data[BP_CAT.y][i] << " " << data[BP_CAT.z][i] << " ";
        outfile << min_dist[i] << " ";
        outfile << voxel_dist[i].x*1024./250.-0.5 << " " << voxel_dist[i].y*1024./250.-0.5 << " " << voxel_dist[i].z*1024./250.-0.5 << " ";
        outfile << id_dist[i];
        outfile << std::endl;
    }

    std::cout << " complete." << std::endl;

    outfile.close();

    //for (int i = 0; i < 100; i++) {





    return -1;
}

/*
    double ** spine_data;
    int num_lines_2=0, num_fields_2=0;
    std::string spine_file = std::string("spine_output.bin");
    
    readBinaryFile(spine_file,num_lines_2,num_fields_2,spine_data);

    if (num_lines != num_lines_2) {
        std::cerr << "Number of halos does not match." << std::endl;
        return -1;
    }

    std::ofstream outfile;
    std::string file;
    
    setOutputFileName();
    
    outfile.open(output.c_str(),std::ofstream::out);
    file = output;

    std::cout << "Writing output file \"" << file << "\"" << std::endl;
    std::cout << "num_lines: " << num_lines << std::endl;
    std::cout << "num_fields: " << num_fields+num_fields_2 << std::endl;

    outfile << "# NUM_LINES " << num_lines << std::endl;
    outfile << "# NUM_FIELDS " << num_fields+num_fields_2 << std::endl;

    int dot_freq = int(ceil((num_fields+num_fields_2)/10.));

    for (int i = 0; i < num_fields; i++) {
        if (i % dot_freq == 0) {
            std::cout << ".";
            std::cout.flush();
        }
        outfile.write((const char *)data[i],num_lines*sizeof(double));
    }

    for (int i = 0; i < num_fields_2; i++) {
        if (i % dot_freq == 0) {
            std::cout << ".";
            std::cout.flush();
        }
        outfile.write((const char *)spine_data[i],num_lines*sizeof(double));
    }

    std::cout << " complete." << std::endl;

    outfile.close();
}
*/

int doSpineTest (double ** &data, int num_lines, int num_fields) {

    // we have z = 0 catalog read into data.  now read in a given spine field

    std::string spine_path = "/pfs/chtlee/shared/spine_in_fila2.bin";

    int num_lines_2=0, num_fields_2=0;
    double ** spine;

    readBinaryFile(spine_path,num_lines_2,num_fields_2,spine);

    struct halo {
        int x = 0;
        int y = 0;
        int z = 0;
    };

    std::vector<int> vspine;
    std::vector<halo> vspinehalos;

    for (int i = 0; i < num_lines_2; i++) {
        vspine.push_back((int)(spine[0][i]));
        int j = vspine.back();
        int x = data[BP_CAT.x][j];
        int y = data[BP_CAT.y][j];
        int z = data[BP_CAT.z][j];
        halo h;
        h.x = x;
        h.y = y;
        h.z = z;
        vspinehalos.push_back(h);
    }

    std::ofstream outfile;

    setOutputFileName();

    outfile.open(output.c_str(),std::ofstream::out);

    for (int i = 0; i < vspinehalos.size(); i++) {
        if (vspinehalos[i].z < 10) {
            outfile << vspinehalos[i].x << " " << vspinehalos[i].y << std::endl;
        }
    }

    std::cout << vmedian(vspine) << " ";
    std::cout << vCIlo(vspine) << " " << vCIhi(vspine) << " ";
    std::cout << vdisp(vspine,0.2) << " " << vdisp(vspine,0.8) << " ";
    std::cout << vspine.size() << std::endl;    

    //std::cout << vmedian(vspinehalos) << " ";
    //std::cout << vCIlo(vspinehalos) << " " << vCIhi(vspinehalos) << " ";
    //std::cout << vdisp(vspinehalos,0.2) << " " << vdisp(vspinehalos,0.8) << " ";
    //std::cout << vspinehalos.size() << std::endl;    

    return 0;
}
