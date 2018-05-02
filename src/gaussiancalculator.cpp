// Given an atom and its neighbors, calculate an associated fingerprint
// Multiple fingerprint generation algorithm are included
//
// Author: Mardochee Reveil
// Creation Date: 9/29/17
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
#include "gaussiancalculator.h"
#include "utilities.h"
#include "periodictable.h"


using namespace std;

GaussianCalculator::GaussianCalculator(AtomicSystem& asys, fingerprintProperties fpprop): fpproperties(fpprop),atomicsystem(asys), ptable() {

    cutoff = fpproperties.cutoff;


    natomtypes = fpproperties.natomtypes;
    string *fp_atomtypes = fpproperties.atomtypes;

    orderedatomtypes = atomicsystem.get_atom_types();
  
    natompairs = 0;
    for (int i=0; i<natomtypes;i++){
        for (int j=i; j<natomtypes; j++){
            vector<string> pair;
            pair.push_back(fp_atomtypes[i]);
            pair.push_back(fp_atomtypes[j]);
            atompairs.push_back(pair); 
            natompairs++;
        }
    }

    nG4s = fpproperties.netas2*fpproperties.nzetas*fpproperties.ngammas;
    nG2s = fpproperties.netas;

    ndirections = fpproperties.ndirections;
    directions = fpproperties.directions;

    nderivatives = fpproperties.nderivatives;

    if (fpproperties.calculate_derivatives == "true") 
        include_derivatives = true;  
    else 
        include_derivatives = false;


    if (fpproperties.strategy == "augmented") {
        fpsize   = nG2s*natomtypes + nG4s*natompairs;    
        if (fpproperties.calculate_derivatives == "true") {
            fpsize  += fpsize*nderivatives*ndirections;
        }   
    } else if  (fpproperties.strategy == "weighted") {
        fpsize   = (nG2s+nG4s)*2;    
        if (fpproperties.calculate_derivatives == "true") {
            fpsize  += fpsize*nderivatives*ndirections;
        }  
    }

}

GaussianCalculator::~GaussianCalculator(){

}

int GaussianCalculator::get_size(void) {return fpsize;}


double* GaussianCalculator::calculate_fingerprint(int atomid, NeighborList &neighlist){ //returns the entire fingerprint for a given atom, including derivatives as needed


    double* fingerprint = new double[fpsize];
    int currentpos = 0;
    int nG2s = fpproperties.netas;

    int nneighbors = neighlist.get_n_neighbors(atomid);
    int *neighbors = neighlist.get_sorted_neighbors(atomid);

    string weight_type = "None";

    if (fpproperties.strategy == "weighted") {

        double neighdist[nneighbors];  

        for (int j=0; j<nneighbors; j++) {
            double distance = sqrt(atomicsystem.get_square_distance(atomid,neighbors[j]));
            neighdist[j] = distance;    
        }

        for (int w=0; w<2; w++) {

            if (w == 0) weight_type = "None";
            else if (w == 1) weight_type = fpproperties.weight_type;
 
            double* tempG2s  = get_G2s(atomid, nneighbors, neighbors, neighdist, weight_type);
            double* tempG4s = get_G4s(atomid, nneighbors, neighbors, neighdist, weight_type);

            for (int i=0; i<nG2s; i++) {
                fingerprint[currentpos] = tempG2s[i];
                currentpos++;
            } 
            for (int i=0; i<nG4s; i++) {
                fingerprint[currentpos] = tempG4s[i];
                currentpos++;
            }

            if (include_derivatives) {

                for (int r=0; r<ndirections; r++){

                    int dir = directions[r]; 

                    for (int d=0; d<nderivatives; d++) {
         
                        double* tempG2primes;
                        double* tempG4primes;


                        if (d==0) {
                            tempG2primes = get_G2_primes(atomid, nneighbors, neighbors, neighdist, atomid, dir, weight_type);
                            tempG4primes  = get_G4_primes(atomid, nneighbors, neighbors, neighdist, atomid, dir, weight_type);

                        }
                        else { 
                            int neigh_id = neighbors[d-1];
                            int neigh_nneighbors= neighlist.get_n_neighbors(neigh_id);
                            int *neigh_neighbors= neighlist.get_sorted_neighbors(neigh_id);
                            double neigh_neighdist[neigh_nneighbors];  
                            for (int j=0; j<neigh_nneighbors; j++) {
                                double distance = sqrt(atomicsystem.get_square_distance(neigh_id,neigh_neighbors[j]));
                                neigh_neighdist[j] = distance;    
                            }
                            tempG2primes  = get_G2_primes(neigh_id, neigh_nneighbors, neigh_neighbors, neigh_neighdist, atomid, dir, weight_type);
                            tempG4primes  = get_G4_primes(neigh_id, neigh_nneighbors, neigh_neighbors, neigh_neighdist, atomid, dir, weight_type);

                        }

                        for (int i=0; i<nG2s; i++) {
                            fingerprint[currentpos] = tempG2primes[i];
                            currentpos++;
                        } 
                        for (int i=0; i<nG4s; i++) {
                            fingerprint[currentpos] = tempG4primes[i];
                            currentpos++;
                        } 
                    }
                }
            }
        }
    }


    else if (fpproperties.strategy == "augmented") {

        for (int atomtype=0; atomtype<natomtypes; atomtype++) {

            int nneighbors_subset = neighlist.get_n_neighbors(atomid,orderedatomtypes[atomtype]);
            int *neighborlist_subset = neighlist.get_sorted_neighbors(atomid,orderedatomtypes[atomtype]);
            double neighdist[nneighbors_subset];  


           for (int j=0; j<nneighbors_subset; j++) {
                double distance = sqrt(atomicsystem.get_square_distance(atomid,neighborlist_subset[j]));
                neighdist[j] = distance;    
            }

            double* tempG2s  = get_G2s(atomid, nneighbors_subset, neighborlist_subset, neighdist, weight_type);

            for (int i=0; i<nG2s; i++) {
                fingerprint[currentpos] = tempG2s[i];
                currentpos++;
            } 

            if (include_derivatives) {

                for (int r=0; r<ndirections; r++){

                    int dir = directions[r]; 

                    for (int d=0; d<nderivatives; d++) {
         
                        double* tempG2primes;

                        if (d==0) {
                             tempG2primes = get_G2_primes(atomid, nneighbors_subset, neighborlist_subset, neighdist, atomid, dir, weight_type);
                        }
                        else { 
                            int neigh_id = neighborlist_subset[d-1];
                            int neigh_nneighbors_subset = neighlist.get_n_neighbors(neigh_id,orderedatomtypes[atomtype]);
                            int *neigh_neighborlist_subset = neighlist.get_sorted_neighbors(neigh_id,orderedatomtypes[atomtype]);
                            double neigh_neighdist[neigh_nneighbors_subset];  
                            for (int j=0; j<neigh_nneighbors_subset; j++) {
                                double distance = sqrt(atomicsystem.get_square_distance(neigh_id,neigh_neighborlist_subset[j]));
                                neigh_neighdist[j] = distance;    
                            }
                            tempG2primes  = get_G2_primes(neigh_id, neigh_nneighbors_subset, neigh_neighborlist_subset, neigh_neighdist, atomid, dir, weight_type);

                        }

                        for (int i=0; i<nG2s; i++) {
                            fingerprint[currentpos] = tempG2primes[i];
                            currentpos++;
                        } 
                    }
                }
            }
        }


        for (int atompair=0; atompair<natompairs; atompair++ ){

            int nneighbors_subset = neighlist.get_n_neighbors(atomid,atompairs[atompair]);
            int *neighborlist_subset = neighlist.get_sorted_neighbors(atomid,atompairs[atompair]);
            double neighdist[nneighbors_subset];        

            for (int j=0; j<nneighbors_subset; j++) {
                double distance = sqrt(atomicsystem.get_square_distance(atomid,neighborlist_subset[j]));
                neighdist[j] = distance;       
            }

            double* tempG4s = get_G4s(atomid, nneighbors_subset, neighborlist_subset, neighdist, weight_type);

            for (int i=0; i<nG4s; i++) {
                fingerprint[currentpos] = tempG4s[i];
                currentpos++;
            }

            if (include_derivatives) {

                for (int r=0; r<ndirections; r++){

                    int dir = directions[r]; 

                    for (int d=0; d<nderivatives; d++) {

                        double* tempG4primes;

                        if (d==0) {
                            tempG4primes  = get_G4_primes(atomid, nneighbors_subset, neighborlist_subset, neighdist, atomid, dir, weight_type);
                        }
                        else {
                            int neigh_id = neighborlist_subset[d-1];
                            int neigh_nneighbors_subset = neighlist.get_n_neighbors(neigh_id,atompairs[atompair]);
                            int *neigh_neighborlist_subset = neighlist.get_sorted_neighbors(neigh_id,atompairs[atompair]);
                            double neigh_neighdist[neigh_nneighbors_subset];  
                            for (int j=0; j<neigh_nneighbors_subset; j++) {
                                double distance = sqrt(atomicsystem.get_square_distance(neigh_id,neigh_neighborlist_subset[j]));
                                neigh_neighdist[j] = distance;    
                            }
                            tempG4primes  = get_G4_primes(neigh_id, neigh_nneighbors_subset, neigh_neighborlist_subset, neigh_neighdist, atomid, dir, weight_type);
          
                        }
    
                        for (int i=0; i<nG4s; i++) {
                            fingerprint[currentpos] = tempG4primes[i];
                            currentpos++;
                        } 
                    }
                }       
            }    
        }
    }

    return fingerprint; 
}


double* GaussianCalculator::get_G2s(int atomid, int nneighbors, int* neighbors, double* distances, string weight_type) {
 
    double* fp;
    fp = new double[fpproperties.netas];
    for (int i=0; i<fpproperties.netas; i++) fp[i] = 0.0;

    int i=0,j=0,k=0;

    int n = 0;

    for (int i=0; i<fpproperties.netas; i++){
        double eta = fpproperties.etas[i];
        double G2 = calculate_G2(nneighbors,neighbors,distances,eta, weight_type);
        fp[n] = G2;
        n++;
    }

    return fp; 
}

double* GaussianCalculator::get_G4s(int atomid,int nneighbors,int* neighbors,double* distances, string weight_type) {

    double* fp;

    int nG4s = fpproperties.netas2*fpproperties.nzetas*fpproperties.ngammas;
    fp = new double[nG4s];
    for (int i=0; i<nG4s; i++) fp[i] = 0.0;

    int n = 0;

    for (int i=0; i<fpproperties.netas2; i++){
        double eta = fpproperties.etas2[i];
        for (int j=0; j<fpproperties.nzetas; j++){
            double zeta = fpproperties.zetas[j];
            for (int k=0; k<fpproperties.ngammas; k++){
                double gamma = fpproperties.gammas[k];
                double G4 = calculate_G4(atomid,nneighbors,neighbors,distances,eta,zeta,gamma,weight_type);
                fp[n] = G4;
                n++;
            }
        }
    }

    return fp;
}

double* GaussianCalculator::get_G2_primes(int atomid2, int nneighbors2, int* neighbors2, double* distances2, int atomid, int direction, string weight_type){

    double* fpprime;
    fpprime = new double[fpsize];

    int n = 0;
    
    for (int i=0; i<fpproperties.netas; i++){
        double eta = fpproperties.etas[i];
        double G2_prime = calculate_G2_prime(atomid2,nneighbors2,neighbors2,distances2,eta,atomid,direction, weight_type);
        fpprime[n] = G2_prime;
        n++;
    }

    return fpprime;
}


double* GaussianCalculator::get_G4_primes(int atomid2, int nneighbors2, int* neighbors2, double* distances2, int atomid, int direction, string weight_type){

    double* fpprime;
    fpprime = new double[fpsize];
 
    int n = 0;

    for (int i=0; i<fpproperties.netas2; i++){
        double eta = fpproperties.etas2[i];
        for (int j=0; j<fpproperties.nzetas; j++){
            double zeta = fpproperties.zetas[j];
            for (int k=0; k<fpproperties.ngammas; k++){
                double gamma = fpproperties.gammas[k];
                double G4_prime = calculate_G4_prime(atomid2,nneighbors2,neighbors2,distances2,eta,zeta,gamma,atomid,direction,weight_type); 
                fpprime[n] = G4_prime;
                n++;
            }
        }
    }

    return fpprime;
}



double GaussianCalculator::calculate_G2(int nneighbors, int* neighbors, double* distances, double eta, string weight_type){

    double value = 0.0;

    for (int i=0; i<nneighbors; i++){ //here we account for all neighbors no matter their types
        double R = distances[i];
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

        value += weight*exp(-1.0*eta*pow(R,2)/pow(cutoff,2))*cutoff_func(R,cutoff);

    }
    return value;
}



double GaussianCalculator::calculate_cos_theta(int i, int j, int k) { //TODO: account for PBC for angle calculation too

    double result = 0.0;

    double* Rij = get_vector(atomicsystem,i,j);
    double* Rik = get_vector(atomicsystem,i,k);

    result += dot(Rij,Rik);

    return result;
}

double GaussianCalculator::calculate_G4(int atomi, int nneighbors, int* neighbors, double* distances, double eta, double zeta, double gamma, string weight_type){

    double value = 0.0, term;

    for (int j=0; j<nneighbors; j++){ //here we account for all neighbors no matter their types
        for (int k=j+1; k<nneighbors; k++){

            double weight = 1.0;
            Atom neigh_atom_j = atomicsystem.get_atom(neighbors[j]); 
            Atom neigh_atom_k = atomicsystem.get_atom(neighbors[k]); 

            if (weight_type == "atomic_number") {
                 double weight_j = ptable.get_atomic_number(neigh_atom_j.get_atom_type());
                 double weight_k = ptable.get_atomic_number(neigh_atom_k.get_atom_type());
                 weight = (weight_j + weight_k) /2;
            }
            else if (weight_type == "electronegativity") {
                double weight_j = ptable.get_electronegativity(neigh_atom_j.get_atom_type());
                double weight_k = ptable.get_electronegativity(neigh_atom_k.get_atom_type());
                weight = (weight_j + weight_k) /2;
            }
            else if (weight_type == "constant") {
                weight = 1.0;
            }

            double Rij = distances[j];
            double Rik = distances[k];
            int atomid1 = neighbors[j];
            int atomid2 = neighbors[k];
            double Rjk = sqrt(atomicsystem.get_square_distance(atomid1,atomid2));

            double cos_theta_ijk = calculate_cos_theta(atomi,neighbors[j],neighbors[k])/(Rij*Rik);

            term   = pow(1.0 + gamma*cos_theta_ijk,zeta);
            term  *= exp(-1.0*eta*(pow(Rij,2.0)+pow(Rik,2.0)+pow(Rjk,2.0))/pow(cutoff,2.0));     
            term  *= cutoff_func(Rij,cutoff)*cutoff_func(Rik,cutoff)*cutoff_func(Rjk,cutoff);     

            value += weight*term;
        }
    }
    value *= pow(2.0,1.0-zeta);
    return value;
}




double* GaussianCalculator::dRij_dRml_vector(int i, int j, int m, int l){

    double* dRij_dRml;
    dRij_dRml = new double[3]; 
    for (int r=0; r<3; r++) dRij_dRml[r] = 0.0;
  
    if ( m!= i && m != j) return dRij_dRml;
    else {

        int p = Kronecker(m,j) - Kronecker(m,i);
        dRij_dRml[0] = p*Kronecker(0,l);
        dRij_dRml[1] = p*Kronecker(1,l);
        dRij_dRml[2] = p*Kronecker(2,l);
 
        return dRij_dRml;

    }

}

double GaussianCalculator::dRij_dRml(int i, int j, double Rij, int m, int l) {

    double dRij_dRml = 0.0;

    Atom Ai = atomicsystem.get_atom(i);
    Atom Aj = atomicsystem.get_atom(j);

    double Ri[3] = {Ai.get_x(),Ai.get_y(), Ai.get_z()};
    double Rj[3] = {Aj.get_x(),Aj.get_y(), Aj.get_z()};

    if (i == j) return 0.0;

    if (m == i)
        dRij_dRml = -1.0 * (Rj[l] - Ri[l])/Rij;
    else if (m == j) 
        dRij_dRml = (Rj[l] - Ri[l])/Rij;
    else dRij_dRml = 0.0;

    return dRij_dRml;
}

double GaussianCalculator::dCos_theta_ijk_dR_ml(int i, int j, int k, double Rij_norm, double Rik_norm, int m, int l) {

    double *Rij_vector = get_vector(atomicsystem,i,j);

    double *Rik_vector = get_vector(atomicsystem,i,k);

    double dCosTheta = 0.0;

    if (Rij_norm == 0.0 || Rik_norm==0.0) return dCosTheta;

    double* dRijdRml = dRij_dRml_vector(i,j,m,l);
    dCosTheta += dot(dRijdRml,Rik_vector)/(Rij_norm*Rik_norm);

    double* dRikdRml = dRij_dRml_vector(i,k,m,l);
    dCosTheta += dot(Rij_vector,dRikdRml)/(Rij_norm*Rik_norm);

    double dRijdRml_val = dRij_dRml(i,j,Rij_norm,m,l);
    if (dRijdRml_val != 0)
        dCosTheta += -1.0 * dot(Rij_vector,Rik_vector)*dRijdRml_val/(pow(Rij_norm,2.0)*Rik_norm);

    double dRikdRml_val = dRij_dRml(i,k,Rik_norm,m,l); 
    if (dRikdRml_val != 0)
        dCosTheta += -1.0 * dot(Rij_vector,Rik_vector)*dRikdRml_val/(Rij_norm*pow(Rik_norm,2.0));

    return dCosTheta;

}




double GaussianCalculator::calculate_G2_prime(int atomi, int nneighbors, int *neighbors, double *distances, double eta, int m, int l, string weight_type) {

    double value = 0;

    for (int j=0; j<nneighbors; j++) {
        int atomj = neighbors[j];
        double Rij = distances[j];
        double dRijdRml = dRij_dRml(atomi,atomj,Rij,m,l);

        double weight = 1.0;
        Atom neigh_atom = atomicsystem.get_atom(neighbors[j]); 

        if (weight_type == "atomic_number") {
             weight = ptable.get_atomic_number(neigh_atom.get_atom_type());
        }
        else if (weight_type == "electronegativity") {
            weight = ptable.get_electronegativity(neigh_atom.get_atom_type());
        }
        else if (weight_type == "constant") {
            weight = 1.0;
        }

        if (dRijdRml != 0) {
            double y = (-2.0*eta*Rij*cutoff_func(Rij,cutoff))/pow(cutoff,2.0) + cutoff_func_prime(Rij,cutoff);
            value += weight*exp(-1.0 * eta * pow(Rij,2.0) / pow(cutoff,2.0)) * y * dRijdRml;            

        }
    }

    return value;

}


 //TODO: Maybe establish sorting by atomic number to process neighbors


double GaussianCalculator::calculate_G4_prime(int atomi, int nneighbors, int *neighbors, 
                       double *distances, double eta, double zeta, double gamma, int m, int l, string weight_type) {


    double value = 0;
    int test = 0;
    for (int j=0; j<nneighbors; j++){
        for (int k=j+1; k<nneighbors; k++){

            double Rij = distances[j];
            double Rik = distances[k];

            double weight = 1.0;
            Atom neigh_atom_j = atomicsystem.get_atom(neighbors[j]); 
            Atom neigh_atom_k = atomicsystem.get_atom(neighbors[k]); 

            if (weight_type == "atomic_number") {
                 double weight_j = ptable.get_atomic_number(neigh_atom_j.get_atom_type());
                 double weight_k = ptable.get_atomic_number(neigh_atom_k.get_atom_type());
                 weight = (weight_j + weight_k) /2;
            }
            else if (weight_type == "electronegativity") {
                double weight_j = ptable.get_electronegativity(neigh_atom_j.get_atom_type());
                double weight_k = ptable.get_electronegativity(neigh_atom_k.get_atom_type());
                weight = (weight_j + weight_k) /2;
            }
            else if (weight_type == "constant") {
                weight = 1.0;
            }


            int atomj = neighbors[j];
            int atomk = neighbors[k];
            double Rjk = sqrt(atomicsystem.get_square_distance(atomj,atomk));

            double cos_theta_ijk = calculate_cos_theta(atomi,neighbors[j],neighbors[k])/(Rij*Rik);

            double y = (1.0 + gamma * cos_theta_ijk);

            double v1 = pow(y,zeta-1.0) * exp(-1.0 * eta * (pow(Rij,2.0)+pow(Rik,2.0)+pow(Rjk,2.0))/pow(cutoff,2.0));
    
            double dCosthetadRml = dCos_theta_ijk_dR_ml(atomi,atomj,atomk,Rij,Rik,m,l);

            double v2 = gamma * zeta * dCosthetadRml; 
 
            double dRijdRml = dRij_dRml(atomi,atomj,Rij,m,l);
            v2 += -2.0 * y * eta * Rij * dRijdRml / pow(cutoff,2.0);

            double dRikdRml = dRij_dRml(atomi,atomk,Rik,m,l);
            v2 += -2.0 * y * eta * Rik * dRikdRml / pow(cutoff,2.0);

            double dRjkdRml = dRij_dRml(atomj,atomk,Rjk, m,l);
            v2 += -2.0 * y * eta * Rjk * dRjkdRml / pow(cutoff,2.0);  

            //if (m==1 && value!= value && test==0) cout<<dCosthetadRml<<" "<<dRjkdRml<<"  "<<dRikdRml<<" "<<dRijdRml<<" <-> ";


            v2 = v2 * cutoff_func(Rij,cutoff) * cutoff_func(Rik,cutoff) * cutoff_func(Rjk,cutoff);

            double v3 = cutoff_func_prime(Rij,cutoff) * dRijdRml * cutoff_func(Rik,cutoff) * cutoff_func(Rjk,cutoff);
            double v4 = cutoff_func(Rij,cutoff) * cutoff_func_prime(Rik,cutoff) * dRikdRml * cutoff_func(Rjk,cutoff);
            double v5 = cutoff_func(Rij,cutoff) * cutoff_func(Rik,cutoff) * cutoff_func_prime(Rjk,cutoff) * dRjkdRml;
            
            value += weight*v1 * (v2 + y*(v3+v4+v5)); 

         //   if (m==1 && value!= value && test==0){
         //   cout<<dCosthetadRml<<" "<<dRjkdRml<<"  "<<dRikdRml<<" "<<dRijdRml<<" <-> " <<value<<" "<<y<<" "<<v1 <<" "<<v2<<" "<<v3<<" "<<v4<<" "<<v5<<"\n"; test=1;
         //   }

        }    
    }

    value = value*pow(2.0,1.0-zeta);

    return value;

}
















