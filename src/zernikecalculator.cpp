#define _USE_MATH_DEFINES

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <math.h>
#include <time.h>
#include <complex>

//#include <boost/lambda/lambda.hpp>
//#include <boost/math/special_functions/spheric_harmonic.hpp>

#include "atom.h"
#include "zernikecalculator.h"
#include "utilities.h"
#include "periodictable.h"

#define TOLERANCE pow(10,-10)

using namespace std;

ZernikeCalculator::ZernikeCalculator(AtomicSystem& asys,fingerprintProperties fpprop): fpproperties(fpprop),atomicsystem(asys),nmax(10),ptable() {

    cutoff = fpproperties.cutoff;
    nmax   = fpproperties.nmax;
    subsize = 0;
    for (int n=0; n<nmax+1; n++) {
        for (int l=0; l<n+1; l++) {
            if ((n-l) % 2 == 0) {
                    subsize++;
            }
        }
    }

    nderivatives = fpproperties.nderivatives;
    ndirections = fpproperties.ndirections;
    directions = fpproperties.directions;

    if (fpproperties.calculate_derivatives == "true") {
        size = subsize*(nderivatives+1)*ndirections;
        include_derivatives = true;
    }
    else {
        size = subsize;
        include_derivatives = false;

    }

    ptable = PeriodicTable();

    factorial_list = new int[4*nmax+3];

    for (int i=0; i<4*nmax+3; i++) factorial_list[i] = factorial(0.5*i);

}

ZernikeCalculator::~ZernikeCalculator(){

}

int ZernikeCalculator::get_size(void) {return size;}


double* ZernikeCalculator::calculate_fingerprint(int atomid,NeighborList &neighlist) {
 
    int nneighbors = neighlist.get_n_neighbors(atomid);
    int *neighbors = neighlist.get_sorted_neighbors(atomid);

    double* fingerprint;
    fingerprint = new double[size];
    for (int i=0; i<size; i++) fingerprint[i] = 0.0;

    double* Znorms = get_Z_norms(atomid, nneighbors,neighbors);

    int currentpos = 0;

    for (int i=0; i<subsize; i++) {
        fingerprint[currentpos] = Znorms[i];
        currentpos++;
    } 
    
    
    if (include_derivatives) {

        for (int r=0; r<ndirections; r++){

            int dir = directions[r]; 

            for (int d=0; d<nderivatives; d++) {
         
                double* Znorms_prime;

                if (d==0) {
                     Znorms_prime = get_Znorms_prime(atomid, nneighbors, neighbors, atomid, dir);
                }
                else { 

                    int neigh_id = neighbors[d-1];
                    int neigh_nneighbors = neighlist.get_n_neighbors(neigh_id);
                    int *neigh_neighborlist = neighlist.get_sorted_neighbors(neigh_id);
                    double neigh_neighdist[neigh_nneighbors];  
   
                    Znorms_prime= get_Znorms_prime(neigh_id, neigh_nneighbors, neigh_neighborlist, atomid, dir);
  
                }
  
                for (int i=0; i<subsize; i++) {
                    fingerprint[currentpos] = Znorms_prime[i];
                    currentpos++;
                } 
            }
        }
    }


    return fingerprint; 
}

double *ZernikeCalculator::get_Z_norms(int atomid, int nneighbors, int* neighbors) {

    double* Z_norms;
    Z_norms = new double[size];
    for (int i=0; i<size; i++) Z_norms[i] = 0.0;

    Atom myatom = atomicsystem.get_atom(atomid); 

    double x0 = myatom.get_x(); 
    double y0 = myatom.get_y(); 
    double z0 = myatom.get_z(); 

    double* psis   = new double[nneighbors];
    double* thetas = new double[nneighbors];
    double* phis   = new double[nneighbors];


    int q = 0;

    for (int n=0; n<nmax+1; n++) {

        for (int l=0; l<n+1; l++) {

            if ((n-l) % 2 == 0) {

                complex<double> znorm(0,0);

                for (int m=0; m<l+1; m++) {

                    complex<double> c_nlm(0,0);

                    for (int neigh=0; neigh<nneighbors; neigh++){

                        Atom neigh_atom = atomicsystem.get_atom(neighbors[neigh]); 

                        double x1 = neigh_atom.get_x(); 
                        double y1 = neigh_atom.get_y(); 
                        double z1 = neigh_atom.get_z(); 

                        double dx = (x1 - x0)/cutoff;
                        double dy = (y1 - y0)/cutoff;
                        double dz = (z1 - z0)/cutoff;

                        double rho = calculate_norm(dx,dy,dz);
                        double weight = 1.0; 

                        if (fpproperties.natomtypes > 1) {

                            if (fpproperties.weight_type == "atomic_number") {
                                weight = ptable.get_atomic_number(neigh_atom.get_atom_type());
                            }
                            else if (fpproperties.weight_type == "electronegativity") {
                                weight = ptable.get_electronegativity(neigh_atom.get_atom_type()); 
                            }
                        }

                        complex<double> zernike = weight*calculate_Z(n,l,m,dx,dy,dz)*cutoff_func(rho*cutoff,cutoff);

                        double theta = 0.0;
                        if (rho > 0.0)
                            theta = acos(dz/rho);

                        double phi = 0.0;
                        if (dx < 0.0)
                            phi = M_PI + atan(dy/dx);
                        else if ((0.0 < dx) and (dy < 0.0))
                            phi = 2*M_PI + atan(dy/dx);
                        else if ((0.0 < dx) and (dy >= 0.0))
                            phi = atan(dy/dx);
                        else if ( (dx == 0.0) and (dy > 0.0) )
                            phi = 0.5*M_PI;
                        else if ( (dx == 0.0) and (dy < 0.0) )
                            phi = 1.5*M_PI;

                    //    zernike = weight*calculate_R(n,l,rho)*sph_harm(m,l,phi,theta)*cutoff_func(cutoff*rho,cutoff);

                        c_nlm += conj(zernike);
                    }

                    complex<double> t = conj(c_nlm);

                    if (m==0)
                        znorm += c_nlm*t;
                    else
                        znorm += 2.0*c_nlm*t;
                }
               // cout<<"Final value "<<real(znorm)<<"\n";
                Z_norms[q] = real(znorm);
                q++; 
            }
        }
    }
 
   return Z_norms;

}

double ZernikeCalculator::calculate_R(int n, int l, double rho) {

    if ((n-l) % 2 != 0)
        return 0.0;

    else {

        double value = 0;
        double k = (n-l)/2;

        for (int s=0; s<k+1; s++) {

            double b1 = binomial(k,s);
            double b2 = binomial(n-s-1+1.5,k);
            value += pow(-1.0,1.0*s)*b1*b2*pow(rho,n-2.0*s);
        }
        value *= pow(2*n+3.0,0.5);
        return value;
    }
}


double ZernikeCalculator::binomial(int n, int k) {

    return factorial_list[2*n]/double(factorial_list[2*k]*factorial_list[2*(n-k)]);
}

double ZernikeCalculator::calculate_q(int nu, int k, int l) {

    double result = pow(-1.0,1.0*k+nu) * pow((2.0*l+4*k+3)/3,0.5);
    result *= binomial(k,nu) * binomial(2*k,k);
    result *= binomial(2*(k+l+nu)+1,2*k)/binomial(k+l+nu,k)/pow(2,2*k);

    return result;
}

complex<double> ZernikeCalculator::calculate_Z(int n, int l, int m, double  x, double y, double z) {

    complex<double> value(0,0);
    complex<double> j(0,1);

    double temp1 = pow( (2.0*l+1)*factorial_list[2*(l+m)]*factorial_list[2*(l-m)],0.5)/factorial_list[2*l];

    temp1 = temp1*pow(2.0,-1.0*m);

    int k = int((n-l)/2.0);

    for (int nu=0; nu<k+1; nu++) {

        double q = calculate_q(nu,k,l);

        for (int alpha=0; alpha<nu+1; alpha++) {

            double b1 = binomial(nu,alpha);

            for (int beta=0; beta<nu-alpha+1; beta++) {

                double b2 = binomial(nu-alpha,beta);
                double temp2 = q*b1*b2;

                for (int u=0; u<m+1; u++) {

                    double b3 = binomial(m,u);

                    complex<double> temp3(0,0);
                    temp3 = pow(-1,m-u)*b3*pow(j,u);

                    for (int mu=0; mu<int((l-m)/2)+1; mu++) {

                        double b4 = binomial(l,mu);
                        double b5 = binomial(l-mu,m+mu);
                        double temp4 = pow(-1,mu*1.0)*pow(2.0,-2.0*mu)*b4*b5;

                        for (int eta=0; eta<mu+1; eta++) {

                            double r = 2*(eta+alpha) + u;
                            double s = 2*(mu - eta + beta) + m - u;
                            double t = 2*(nu - alpha - beta - mu) + l -m;

                            complex<double> temp(0,0);
                            temp += temp2*temp3*temp4*binomial(mu,eta);

                            temp *= pow(x,r)*pow(y,s)*pow(z,t);
                            value += temp;

                        }
                    }
                } 
            }
        }
    }

    value *= temp1*pow(j,m*1.0)/pow(4*M_PI/3,0.5);

    return value;
}

double* ZernikeCalculator::get_Znorms_prime (int atomid, int nneighbors, int* neighbors, int p, int direction) {

    double* fpprime;
    fpprime = new double[size];
    for (int i=0; i<size; i++) fpprime[i] = 0.0;


    Atom myatom = atomicsystem.get_atom(atomid); 

    double x0 = myatom.get_x(); 
    double y0 = myatom.get_y(); 
    double z0 = myatom.get_z(); 

    int q = 0;
    double weight = 14.0;

    if (fpproperties.weight_type == "atomic_number") {

        weight = ptable.get_atomic_number(myatom.get_atom_type());

    }

    for (int n=0; n<nmax+1; n++) {

        for (int l=0; l<n+1; l++) {

            if ((n-l) % 2 == 0) {

                complex<double> norm_prime(0,0);
                
                for (int m=0; m<l+1; m++) {

                    complex<double> c_nlm(0,0);
                    complex<double> c_nlm_prime(0,0);
               
                    for (int i=0; i<nneighbors; i++) {

                       Atom neigh_atom = atomicsystem.get_atom(neighbors[i]); 
                        int n_index = neighbors[i];
                        double x1 = neigh_atom.get_x(); 
                        double y1 = neigh_atom.get_y(); 
                        double z1 = neigh_atom.get_z(); 

                        double dx = (x1 - x0)/cutoff;
                        double dy = (y1 - y0)/cutoff;
                        double dz = (z1 - z0)/cutoff;

                        double rho = calculate_norm(dx,dy,dz);

                        complex<double> temp_Z_nlm = calculate_Z(n,l,m,dx,dy,dz);

                        complex<double> Z_nlm = temp_Z_nlm*cutoff_func(rho*cutoff,cutoff);

                        int coef = 1;
                        if (p == neighbors[i])
                            coef = -1;
                        else if (p == atomid)
                            coef = 1;
                        else 
                            coef = 0;

                        complex<double> Z_nlm_prime(0,0);
                        Z_nlm_prime += temp_Z_nlm*cutoff_func_prime(rho*cutoff,cutoff)*der_position(x0,y0,z0,x1,y1,z1,coef,direction);

                        complex<double> temp_Z_nlm_prime = calculate_Z_prime(n,l,m,dx,dy,dz,direction);


                        if ((Kronecker(n_index,p) - Kronecker(atomid,p)) == 1) //TODO: revisit this n_index!!!
                            Z_nlm_prime += cutoff_func(rho*cutoff,cutoff)*temp_Z_nlm_prime/cutoff;

                        else if (Kronecker(n_index,p) - Kronecker(atomid,p) == -1)
                            Z_nlm_prime -= cutoff_func(rho*cutoff,cutoff)*temp_Z_nlm_prime/cutoff;

                        c_nlm += weight*conj(Z_nlm);
                        c_nlm_prime += weight*conj(Z_nlm_prime);                      
                 
                    }                           

                    complex<double> t = conj(c_nlm_prime);

                    if (m == 0) {
                        norm_prime += 2.0*c_nlm*t;
                    }
                    else{
                        norm_prime += 4.0*c_nlm*t;
                    }
                }

                fpprime[q] = real(norm_prime);
                q++;
            }
        }
    }

    return fpprime;

}

double ZernikeCalculator::der_position (double x0, double y0, double z0, double x1, double y1, double z1, int coef, int dir) {

    double der = 0.0;
    double dist = calculate_norm(x0-x1,y0-y1,z0-z1);

    if (dir == 0)
        der = coef*(x0 - x1) / dist;
    else if (dir == 1)
        der = coef*(y0 - y1) / dist;
    else if (dir == 2)
        der = coef*(z0 - z1) / dist;

    return der;
} 



complex<double> ZernikeCalculator::calculate_Z_prime(int n, int l, int m, double x, double y, double z, int p) {

    complex<double> value(0,0);
    complex<double> j(0,1);
    double temp1 = pow(2,-1.0*m)*pow((2*l+1)*factorial_list[2*(l+m)]*factorial_list[2*(l-m)],0.5)/factorial_list[2*l];
    
    int k = (n-l)/2.0;

    for (int nu=0; nu<k+1; nu++) {

        double q = calculate_q(nu,k,l);

        for (int alpha=0; alpha<nu+1; alpha++){

            double b1 = binomial(nu,alpha);

            for (int beta=0; beta<nu - alpha + 1; beta++){

                double b2 = binomial(nu - alpha, beta);
                double temp2 = q * b1 * b2;

                for (int u=0; u<m + 1; u++) { 

                    complex<double> temp3(0,0);
                    temp3 += pow(-1,m - u) * binomial(m, u) * pow(j,u);

                    for (int mu=0; mu<int((l - m) / 2.) + 1; mu++) {

                        double temp4 = pow(-1.,mu) * pow(2.,-2. * mu);
                        temp4 *= binomial(l, mu);
                        temp4 *= binomial(l - mu, m + mu);
                        
                        for (int eta = 0; eta<mu + 1; eta++) {

                            double r = 2.0 * (eta + alpha) + u;
                            double s = 2.0 * (mu - eta + beta) + m - u;
                            double t = 2.0 * (nu - alpha - beta - mu) + l - m;

                            complex<double> c(0,0);
                         //   cout<<"mu and eta: "<<mu<<" "<<eta<<"\n";
                            double v =  binomial(mu, eta);
                            c += temp3*temp2*temp4*v;

                            if (p == 0) {
                                if (r != 0)
                                    value += c * r * pow(x, (r - 1)) * pow(y , s) * pow(z , t);
                            } else if (p == 1) {
                                if (s != 0)
                                    value += c * s * pow(x, r) * pow(y, (s - 1)) * pow(z, t);
                            } else if (p == 2){
                                if (t != 0)
                                    value += c * t * pow(x, r) * pow(y, s) * pow(z, (t - 1));
                            }
                        }
                    }
                }
            }
        }
    }
    
    value = value*temp1*pow(j,m)/pow(4*M_PI/3,0.5);
    return value;
}
