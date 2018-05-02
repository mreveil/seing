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
#include "periodictable.h"


using namespace std;

AGNICalculator::AGNICalculator(AtomicSystem& asys, fingerprintProperties fpprop): fpproperties(fpprop),atomicsystem(asys) {

    cutoff = fpproperties.cutoff;


    natomtypes = fpproperties.natomtypes;
    string *fp_atomtypes = fpproperties.atomtypes;

    orderedatomtypes = atomicsystem.get_atom_types();
  
    w = fpproperties.width;
    alpha = fpproperties.alpha;  //direction of the force between 0 and 2 

    dim = fpproperties.dimensionality;

    if (fpproperties.strategy == "augmented") {
        fpsize   = dim*natomtypes;
    } else if (fpproperties.strategy == "weighted") {
        fpsize = dim*2;
    }

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

    string weight_type = "None";

    if (fpproperties.strategy == "weighted") {

        double neighdist[nneighbors];  
        double dist_components[nneighbors];  

        for (int j=0; j<nneighbors; j++) {

            double distance = sqrt(atomicsystem.get_square_distance(atomid,neighbors[j]));
            neighdist[j] = distance;    

            double r_alpha = atomicsystem.get_distance_component(atomid,neighbors[j],alpha);
            dist_components[j] = r_alpha;

        }


        for (int w=0; w<2; w++) {

            if (w == 0) weight_type = "None";
            else if (w == 1) weight_type = fpproperties.weight_type;
        
            for (int k=0; k<dim; k++){
                double a_k = centers[k];
                double component = calculate_component(nneighbors,neighbors, neighdist,dist_components, a_k, weight_type);
                fingerprint[currentpos] = component;
                currentpos++;
            }  

        }
    } else if (fpproperties.strategy == "augmented") {

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
                double component = calculate_component(nneighbors_subset, neighborlist_subset, neighdist, dist_components,a_k,weight_type);
                fingerprint[currentpos] = component;
                currentpos++;
            }         
        } 
    }
    return fingerprint; 
}



double AGNICalculator::calculate_component(int nneighbors, int* neighbors, double* distances, double* dist_components, double a_k, string weight_type){

    double value = 0.0;
    PeriodicTable ptable;


    for (int i=0; i<nneighbors; i++){ //here we account for all neighbors no matter their types

        double weight = 1.0;
        Atom neigh_atom = atomicsystem.get_atom(neighbors[i]); 

        if (weight_type == "atomic_number") {
             weight = ptable.get_atomic_number(neigh_atom.get_atom_type());
        }
        else if (weight_type == "electronegativity") {
            weight = ptable.get_electronegativity(neigh_atom.get_atom_type());
        }
        else if (weight_type == "constant") {
            weight = 1.0;
        }


        double R = distances[i];
        double R_alpha = dist_components[i];
        double temp = (R_alpha/R)*(1/(pow(2*M_PI,0.5)*w))*exp(-0.5*pow((R-a_k)/w,2))*cutoff_func(R,cutoff);
        value += weight*temp;

    }
    return value;
}

