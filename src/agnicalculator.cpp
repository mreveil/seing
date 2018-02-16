// Implementation of the AGNI fingerprinting scheme
// REFERENCE: https://www.nature.com/articles/s41524-017-0042-y
//
// 
// Author: Mardochee Reveil
// Creation Date: 1/15/17
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
#include <algorithm>

#include "atom.h"
#include "agnicalculator.h"
#include "utilities.h"

using namespace std;

AGNICalculator::AGNICalculator(AtomicSystem& asys, fingerprintProperties fpprop): fpproperties(fpprop),atomicsystem(asys) {

    cutoff = fpproperties.cutoff;


    natomtypes = fpproperties.natomtypes;
    string *fp_atomtypes = fpproperties.atomtypes;

    orderedatomtypes = atomicsystem.get_atom_types();
  
    w = fpproperties.width;
    alpha = fpproperties.alpha;  //direction of the force between 0 and 2 

    dim = fpproperties.dimensionality;
    fpsize   = dim*natomtypes;

    double starting_point = 0.0;
  //  double starting_point = cutoff/2;

    centers = new double[dim];
    for (int i=0; i<dim; i++){
        double v = starting_point + i*(cutoff-starting_point)/(dim-1);
        centers[i] = v;
    }
    
}

AGNICalculator::~AGNICalculator(){

}

int AGNICalculator::get_size(void) {return fpsize;}


double* AGNICalculator::calculate_fingerprint(int atomid, NeighborList &neighlist){ //returns the entire fingerprint for a given atom, including derivatives as needed


    double* fingerprint = new double[fpsize];
    int currentpos = 0;
   
    int nneighbors = neighlist.get_n_neighbors(atomid);
    int *neighbors = neighlist.get_sorted_neighbors(atomid);


    for (int atomtype=0; atomtype<natomtypes; atomtype++) {

        int nneighbors_subset = neighlist.get_n_neighbors(atomid,orderedatomtypes[atomtype]);
        int *neighborlist_subset = neighlist.get_sorted_neighbors(atomid,orderedatomtypes[atomtype]);

        double neighdist[nneighbors_subset];  
        double dist_components[nneighbors_subset];  


        for (int j=0; j<nneighbors_subset; j++) {

            double distance = sqrt(atomicsystem.get_square_distance(atomid,neighborlist_subset[j]));
            neighdist[j] = distance;    

            double r_alpha = atomicsystem.get_distance_component(atomid,neighborlist_subset[j],alpha);
            dist_components[j] = r_alpha;
        }

        
        for (int k=0; k<dim; k++){
            double a_k = centers[k];
            double component = calculate_component(nneighbors_subset,neighdist,dist_components,a_k);
            fingerprint[currentpos] = component;
            currentpos++;
        }         
    } 


    return fingerprint; 
}



double AGNICalculator::calculate_component(int nneighbors, double* distances, double* dist_components, double a_k){

    double value = 0.0;

    for (int i=0; i<nneighbors; i++){ //here we account for all neighbors no matter their types
        double R = distances[i];
        double R_alpha = dist_components[i];
        double temp = (R_alpha/R)*(1/(pow(2*M_PI,0.5)*w))*exp(-0.5*pow((R-a_k)/w,2))*cutoff_func(R,cutoff);
        value += temp;

    }
    return value;
}
