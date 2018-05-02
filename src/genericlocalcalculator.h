#ifndef GENERICCALCULATOR_H
#define GENERICCALCULATOR_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <time.h>


#include "atomicsystem.h"
#include "inputs.h"
#include "utilities.h"
#include "neighborlist.h"

#include "zernikecalculator.h"
#include "bispectrumcalculator.h"
#include "agnicalculator.h"


#include "gaussiancalculator.h"


using namespace std;

//! Class that holds a generic calculator to switch between actual fingerprint calculators.

class GenericLocalCalculator { 

    int size;
    double cutoff;
    AtomicSystem atomicsystem;
    fingerprintProperties fpproperties;

    public:
        //! Constructor based on the atomtic system and fingerprintProperties 
        GenericLocalCalculator(AtomicSystem&, fingerprintProperties);
        ~GenericLocalCalculator();

        //! Returns the dimensionality of the fingerprint
        int get_size();

        /*! Will call the calculate_fingerprint function of the actual fingerprint calculator and returns the 
            fingerprint
        */
        double *calculate_fingerprint(int, NeighborList&);

};


#endif
