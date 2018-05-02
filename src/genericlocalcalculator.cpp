// Given an atom and its neighbors, calculate an associated fingerprint
// Multiple fingerprint generation algorithm are included
//
// Author: Mardochee Reveil
// Date: 9/29/17
//
//



#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <math.h>
#include <time.h>

#include "atom.h"
#include "genericlocalcalculator.h"

#include "utilities.h"

using namespace std;

GenericLocalCalculator::GenericLocalCalculator(AtomicSystem& asys, fingerprintProperties fpprop): fpproperties(fpprop),atomicsystem(asys) {

    cutoff = fpproperties.cutoff;

    if (fpproperties.type == "gaussian"){
        GaussianCalculator gcalc = GaussianCalculator(atomicsystem,fpproperties); 
        size = gcalc.get_size();
    }

    else if (fpproperties.type == "zernike"){
        ZernikeCalculator zcalc = ZernikeCalculator(atomicsystem,fpproperties); 
        size = zcalc.get_size();
    }

    else if (fpproperties.type == "bispectrum"){
        BispectrumCalculator bcalc = BispectrumCalculator(atomicsystem,fpproperties); 
        size = bcalc.get_size();
    }
    else if (fpproperties.type == "agni"){
        AGNICalculator acalc = AGNICalculator(atomicsystem,fpproperties); 
        size = acalc.get_size();
    }
    else {
        cerr<<"ERROR: Fingerprint type "<<fpproperties.type<<"not supported.\n";
    }

}

GenericLocalCalculator::~GenericLocalCalculator(){

}

int GenericLocalCalculator::get_size(void) {

    return size;
}

double* GenericLocalCalculator::calculate_fingerprint(int atomid, NeighborList &nlist) {//int nneighbors, int* neighbors) {
 
    double* fp;
 
    if (fpproperties.type == "gaussian"){
        GaussianCalculator gcalc = GaussianCalculator(atomicsystem,fpproperties); 
        fp = gcalc.calculate_fingerprint(atomid, nlist); //nneighbors, neighbors);
    }
    else if (fpproperties.type == "zernike"){
        ZernikeCalculator zcalc = ZernikeCalculator(atomicsystem,fpproperties); 
        fp = zcalc.calculate_fingerprint(atomid, nlist);
    }
    else if (fpproperties.type == "bispectrum"){
        BispectrumCalculator bcalc = BispectrumCalculator(atomicsystem,fpproperties); 
        fp = bcalc.calculate_fingerprint(atomid,nlist);
    }
    else if (fpproperties.type == "agni"){
        AGNICalculator dcalc = AGNICalculator(atomicsystem,fpproperties); 
        fp = dcalc.calculate_fingerprint(atomid,nlist);
    }
    else {
        cerr<<"ERROR: Fingerprint type "<<fpproperties.type<<"not supported.\n";
    }

    return fp;
}

