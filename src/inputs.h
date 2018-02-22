#ifndef INPUTS_H 
#define INPUTS_H

#include <string>

using namespace std;


//! Structure to hold fingerprint parameters as provided by the user in the input file
struct fingerprintProperties {

    //! The type of fingerprint to compute (e.g. AGNI, Zernike, Behler-Parinello, etc.)
    string type;

    //! Values: yes or no 
    string calculate_derivatives;

    //! The number of derivatives to calculate 
    /*! Derivatives are calculated with respect to the central atom first, then wrt neighbor atoms 
        ordered by increasing distance.
    */
    int nderivatives;

    //! How many directions (x,y,z) to calculate derivatives for
    int ndirections;

    //! Directions to calculate derivatives for 0=x, 1=y, 2=z
    int *directions;

    //! The strategy to use for incorporating multiple species in the fingerprint. Values: weighted or augmented
    string strategy;

    //! For the "weighted" strategy, what type of weight to use (e.g. atomic number, elctronegativity, etc.)
    string weight_type;

    //! In case of a local fingerprint, the cutoff to use
    double cutoff;
    double inner_cutoff;

    //! A boolean to save whether the box size is provided as part of the input file or not
    bool is_box_size_provided;

   //! A pointer to the list that holds the dimensions of the box if they are provided in input file
    double *box_size;


    //! The number of different species to consider for fingerprint generation
    int natomtypes;

    //! A pointer to the list of atomic species for fingerprint generation
    string *atomtypes;

    // Behler_Parinello Fingerprints 

    int netas;
    double *etas;

    int netas2;
    double *etas2;

    int ngammas;
    double *gammas;

    int nzetas;
    double *zetas;

    int nmax;  //Used for Zernike fingerprint calculator
    int jmax;  //Used for Bispectrum fingerprint calculator

    // AGNI Fingerprints

    int dimensionality;
    double width;
    int alpha;

    string output_file;
    string output_mode;
};

fingerprintProperties read_prop_file(string);



#endif
