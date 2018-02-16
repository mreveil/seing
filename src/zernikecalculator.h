#ifndef ZERNIKE_H
#define ZERNIKE_H

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
#include "neighborlist.h"
#include "periodictable.h"

using namespace std;


class ZernikeCalculator {

    int size, subsize;
    int nderivatives;
    int ndirections;
    int *directions;
    bool include_derivatives;
    double cutoff;
    int nmax;
    int* factorial_list;

    PeriodicTable ptable;
   
    AtomicSystem atomicsystem;
    fingerprintProperties fpproperties;

    double der_position (double x0, double y0, double z0, double x1, double y1, double z1, int coef, int dir);
    complex<double> calculate_Z_prime(int n, int l, int m, double x, double y, double z, int p);
    complex<double> calculate_Z(int n, int l, int m, double  x, double y, double z);
    double calculate_q(int nu, int k, int l);
    double calculate_R(int n, int l, double rho);
    double binomial(int n, int k);

    double *get_Z_norms(int atomid, int nneighbors, int* neighbors);
    double *get_Znorms_prime (int atomid, int nneighbors, int* neighbors, int p, int direction);

    
    public:
        ZernikeCalculator (AtomicSystem&, fingerprintProperties);
        ~ZernikeCalculator ();

        int get_size();
        double *calculate_fingerprint(int, NeighborList&);



};


#endif
