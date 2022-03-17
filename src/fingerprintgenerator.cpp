// This will generate fingerprint for each atom in a given atomic system
//
// Author: Mardochee Reveil
// Date created: 9/29/17
//

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <time.h>
#include <math.h>

#include "fingerprintgenerator.h"
#include "gaussiancalculator.h"
#include "bispectrumcalculator.h"
#include "zernikecalculator.h"
#include "genericlocalcalculator.h"

#include "inputs.h"

using namespace std;

FingerprintGenerator::FingerprintGenerator(AtomicSystem& asys, fingerprintProperties fpproperties): atomicsystem(asys), natoms(0), fp_natomtypes(1), natompairs(0), fingerprints(NULL){

    natoms = atomicsystem.get_n_atoms();

    double cutoff_skin = 0.5;
    double cutoff = fpproperties.cutoff;
    
    double xsize = asys.get_xsize();
    double ysize = asys.get_ysize();
    double zsize = asys.get_zsize();

    int x_nbins = ((int) xsize/cutoff > 0) ? (int) xsize/cutoff : 1;
    int y_nbins = ((int) ysize/cutoff > 0) ? (int) ysize/cutoff : 1;
    int z_nbins = ((int) zsize/cutoff > 0) ? (int) zsize/cutoff : 1;

    int maxneighborsperatom = 1000;

    cout<<"Generating neighborlist...";

    NeighborList neighlist = NeighborList(atomicsystem, cutoff, x_nbins, y_nbins, z_nbins, maxneighborsperatom);


   //TODO: Check if atom types from atomic system overlaps with atom types from fingerprint properties input file

    try {        
        neighlist.build();
        cout<<"..........done\n";
    } catch (const char* msg) {
        cerr<<msg<<"\n";
    }

    for (int i=0; i<fpproperties.ndirections;i++) cout<<fpproperties.directions[i]<<" ";

    GenericLocalCalculator fpcalc = GenericLocalCalculator(atomicsystem,fpproperties); 
    fsize = fpcalc.get_size();

    cout<<"Done initializing fingerprint calculator...\n";

    fingerprints = new double*[natoms]; 
    for (int i=0; i<natoms; i++)
        fingerprints[i] = new double[fsize]; //automatically includes derivatives as per user request in input file



    cout<<"\tTotal dimensionality of fingerprint is:"<<fsize<<"\n";

    cout<<"Generating fingerprints...";

    for (int atomid=0; atomid<natoms; atomid++) {

       // cout<<"Working with atom: "<<atomid<<"\n";
        //int nneighbors = neighlist.get_n_neighbors(atomid);
       // int *myneighbors = neighlist.get_sorted_neighbors(atomid);
       // fingerprints[atomid] = fpcalc.calculate_fingerprint(atomid,nneighbors,myneighbors);
       fingerprints[atomid] = fpcalc.calculate_fingerprint(atomid,neighlist);

    } 

    cout<<"All done\n";

}



FingerprintGenerator::~FingerprintGenerator (void) {

    delete [] fingerprints;
}

bool FingerprintGenerator::write2file(string filename, string mode) {
    bool success = true;
    ofstream myfile;
    
    if (mode == "append")
        myfile.open(filename,std::ios_base::app);
    else 
        myfile.open(filename);

    if (myfile.is_open()) {
        for (int i=0; i<natoms; i++){
            for (int t=0; t<fsize; t++){
                myfile << fingerprints[i][t] << " ";
            } 
 
            myfile << "\n";            
        }
        myfile.close();
    }
    else success = false;

    //return success;

    sample1 = fingerprints[3][1];
    sample2 = fingerprints[9][2];

    return success;
}

//sample1 = *fingerprints[3];
//sample2 = *fingerprints[9];
