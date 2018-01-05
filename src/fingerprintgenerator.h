#ifndef FINGERPRINTGENERATOR_H
#define FINGERPRINTGENERATOR_H

#include <string>

#include "atomicsystem.h"
#include "neighborlist.h"
#include "inputs.h"


using namespace std;

class FingerprintGenerator {

    AtomicSystem atomicsystem;
    double **fingerprints;

    int fsize, natoms, fp_natomtypes, natompairs;

    public:
        FingerprintGenerator (AtomicSystem&,fingerprintProperties);
        ~FingerprintGenerator(void);

        bool write2file(string,string);
};

#endif
