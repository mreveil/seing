#ifndef AGNICALCULATOR_H
#define AGNICALCULATOR_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <time.h>
#include <vector>


#include "atomicsystem.h"
#include "inputs.h"
#include "neighborlist.h"


using namespace std;

//! An AGNI fingerprint calculator 

class AGNICalculator {

    int fpsize;
    double w;
    double *centers;
    int alpha;
    int dim;

    double cutoff;

    int natomtypes;
    
    vector<string> orderedatomtypes;
   
    AtomicSystem atomicsystem;
    fingerprintProperties fpproperties;

    double calculate_component(int, int*, double*, double*, double, string);


    public:
        //! Constructor that instantiate an AGNI fingerprint object 
        AGNICalculator(AtomicSystem&, fingerprintProperties);
        ~AGNICalculator();

        //! Returns the dimensionality of this AGNI fingerprint
        int get_size();
       // double *calculate_fingerprint(int, int, int*);

        //! Function to calculate the fingerprint of an atom given its neighbor list
        double *calculate_fingerprint(int, NeighborList&);

};


#endif
