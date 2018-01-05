#ifndef INPUTS_H 
#define INPUTS_H

#include <string>

using namespace std;

struct fingerprintProperties {

    string type;

    string calculate_derivatives;
    int nderivatives;
    int direction;

    string strategy;
    string weight_type;

    double cutoff;
    double inner_cutoff;


    bool is_box_size_provided;
    double *box_size;


    int netas;
    double *etas;

    int netas2;
    double *etas2;

    int ngammas;
    double *gammas;

    int nzetas;
    double *zetas;

    int natomtypes;
    string *atomtypes;

    int nmax;  //Used for Zernike fingerprint calculator
    int jmax;  //Used for Bispectrum fingerprint calculator

    string output_file;
    string mode;
};

fingerprintProperties read_prop_file(string);



#endif
