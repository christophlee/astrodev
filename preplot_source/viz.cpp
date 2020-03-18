#include "viz.h"

int doConvertToVizFormat (double ** &data, int num_lines, int num_fields) {

    std::ofstream of;

    setOutputFileName();

    of.open(output.c_str(),std::ofstream::out);

    for (int i = 0; i < num_lines; i++) {

        // now we need to print additional info like rotation component
        double ax = data[48][i]; double ay = data[49][i]; double az = data[50][i];
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
                     of << data[k][i] << " ";
                     break;
                 case 1:     // id
                 case 3:     // descid
                     of << (long int)data[k][i] << " ";
                     break;
                 default:
                     break;
             }
         }

         of << xrot << " " << zrot << std::endl;
    }

    of.close();

    return 0;

}
