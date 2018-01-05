#ifndef DIAMOND_H
#define DIAMOND_H

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
#include "utilities.h"
#include "neighborlist.h"

using namespace std;


class DiamondCalculator {

    int size; // number of features
    double cutoff, dx, dy, dz;
    int* ngaussians; //number of gaussians per feature
    double lattice_constant;
    double inner_cutoff;
   
    AtomicSystem atomicsystem;
    fingerprintProperties fpproperties;

    vector<vector<double>> generateGaussians(double,double,double);
    vector<vector<vector<double>>> findAtomsForEachGaussian(vector<vector<double>> gaussians, int atomid, int nneighbors, int* neighbors);

    double calculateGaussianValue(vector<vector<double>>, vector<double>); 
    double calculateDistanceToGaussian(int atomid,vector<double> gaussian);

    public:
        DiamondCalculator(AtomicSystem&, fingerprintProperties);
        ~DiamondCalculator();

        int get_size();
        double *calculate_fingerprint(int, NeighborList&);
       // double *calculate_fingerprint_prime(int, int, int*, double*, int, int);

};


#endif
