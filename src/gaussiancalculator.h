#ifndef GAUSSIANCALCULATOR_H
#define GAUSSIANCALCULATOR_H

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
#include "periodictable.h"


using namespace std;


class GaussianCalculator {

    int fpsize;
    int nG4s;
    int nG2s;

    double cutoff;
    bool include_derivatives;
    int nderivatives;
    int ndirections;
    int *directions;

    int natomtypes;
    int natompairs;
    vector<vector<string>> atompairs;
    vector<string> orderedatomtypes;
   
    AtomicSystem atomicsystem;
    fingerprintProperties fpproperties;

    PeriodicTable ptable;


    double calculate_G2(int, int*, double*, double, string);
    double calculate_G4(int, int, int*, double*, double, double, double, string);
    double calculate_cos_theta(int,int,int); 

    double *dRij_dRml_vector(int, int, int, int);
    double dCos_theta_ijk_dR_ml(int i, int j, int k, double Rij, double Rik, int m, int l);
    double dRij_dRml(int, int, double, int, int);
    double calculate_G2_prime(int, int, int *, double *, double, int, int, string);
    double calculate_G4_prime(int, int, int *, double *, double, double, double, int, int, string);

    double *get_G2s(int, int, int*, double*, string);
    double *get_G4s(int, int, int*, double*, string);
    double *get_G2_primes(int, int, int*, double*, int, int, string);
    double *get_G4_primes(int, int, int*, double*, int, int, string);


    int* sort_neighbors(int id, int nneighbors, int* neighbors);
    int* get_sorted_neighbors_subset(int id, int nneighbors, int* neighbors, string atomtype);
    int* get_sorted_neighbors_subset(int id, int nneighbors, int* neighbors, vector<string> atomtypes);
    int get_n_neighbors_subset(int nneighbors, int* neighbors, string atomtype);
    int get_n_neighbors_subset(int nneighbors, int* neighbors, vector<string> atomtypes);

    public:
        GaussianCalculator(AtomicSystem&, fingerprintProperties);
        ~GaussianCalculator();

        int get_size();
       // double *calculate_fingerprint(int, int, int*);
        double *calculate_fingerprint(int, NeighborList&);

};


#endif
