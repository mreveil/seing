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
#include "diamondcalculator.h"

#include "gaussiancalculator.h"


using namespace std;


class GenericLocalCalculator { 

    int size;
    double cutoff;
    AtomicSystem atomicsystem;
    fingerprintProperties fpproperties;

    public:

        GenericLocalCalculator(AtomicSystem&, fingerprintProperties);
        ~GenericLocalCalculator();

        int get_size();
        double *calculate_fingerprint(int, NeighborList&);

};


#endif
