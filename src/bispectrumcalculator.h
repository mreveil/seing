#ifndef BISPECTRUM_H
#define BISPECTRUM_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <time.h>
#include <complex>


#include "atomicsystem.h"
#include "inputs.h"
#include "utilities.h"
#include "neighborlist.h"

using namespace std;


class BispectrumCalculator {

    int size;
    double cutoff;
    int jmax;
   
    AtomicSystem atomicsystem;
    fingerprintProperties fpproperties;

    double calculate_B(double j1, double j2, int j, int nneighbors, double* distances, double* psis, double* thetas, double* phis);
    complex<double> Wigner(double j,double m, double mprime,double alpha, double beta, double gamma);
    complex<double> U(double j,double m,double mprime,double omega,double theta,double phi);
    double get_CG_coefficient(double a, double alpha, double b, double beta, double c, double gamma);
    complex<double> calculate_c(double j, double mprime, double m, int nneighbors, double* distances, double* psis, double* thetas, double* phis);
     
    public:
        BispectrumCalculator(AtomicSystem&, fingerprintProperties);
        ~BispectrumCalculator();

        int get_size();
        double *calculate_fingerprint(int, NeighborList&);
       // double *calculate_fingerprint_prime(int, int, int*, double*, int, int);

};


#endif
