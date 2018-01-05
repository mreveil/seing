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

#include "atom.h"
#include "bispectrumcalculator.h"
#include "utilities.h"

#define TOLERANCE pow(10,-10)

using namespace std;

BispectrumCalculator::BispectrumCalculator(AtomicSystem& asys,fingerprintProperties fpprop): fpproperties(fpprop),atomicsystem(asys),jmax(10) {

    cutoff = fpproperties.cutoff;
    jmax = fpproperties.jmax;
    size   = 0;

    for (int _2j1=0; _2j1 < 2*jmax+1; _2j1++) {
        for (int j=0; j<min(_2j1,jmax) + 1; j++) {
            size++;
        }
    }

}

BispectrumCalculator::~BispectrumCalculator(){

}

int BispectrumCalculator::get_size(void) {return size;}


double* BispectrumCalculator::calculate_fingerprint(int atomid,NeighborList &nlist) {
 
    int nneighbors = nlist.get_n_neighbors(atomid);
    int* neighbors = nlist.get_sorted_neighbors(atomid);

    double* distances = new double[nneighbors];

    for (int j=0; j<nneighbors; j++) {
        double distance = sqrt(atomicsystem.get_square_distance(atomid,neighbors[j]));
        distances[j] = distance;    
    }

    double* fp;
    fp = new double[size];
    for (int i=0; i<size; i++) fp[i] = 0.0;

    Atom myatom = atomicsystem.get_atom(atomid); 

    double x0 = myatom.get_x(); 
    double y0 = myatom.get_y(); 
    double z0 = myatom.get_z(); 

    double* psis   = new double[nneighbors];
    double* thetas = new double[nneighbors];
    double* phis   = new double[nneighbors];

    for (int i=0; i<nneighbors; i++) {

        Atom neigh_atom = atomicsystem.get_atom(neighbors[i]); 

        double x1 = neigh_atom.get_x(); 
        double y1 = neigh_atom.get_y(); 
        double z1 = neigh_atom.get_z(); 

        double dx = x1 - x0;
        double dy = y1 - y0;
        double dz = z1 - z0;

        double dist = distances[i];
        double psi = asin(dist/cutoff);
        double theta = acos(dz/dist);

        if (abs(dz/dist - 1.0) < pow(10,-6)) 
            theta = 0.0;
        else if (dz/dist + 1.0 < pow(10,-6)) 
            theta = M_PI;

        double phi = 0.0;
        if (dx < 0)
            phi = M_PI + atan(dy/dx);
        else if ( (dx > 0.0) and (dy < 0.0))
            phi = 2 * M_PI + atan(dy/dx);
        else if ( (dx > 0.0) and (dy >= 0))
            phi = atan(dy/dx);
        else if ( (dx == 0) and (dy > 0.0) )
            phi = 0.5*M_PI;
        else if ( (dx == 0) and (dy < 0.0) )
            phi = 1.5*M_PI;

        psis[i] = psi;
        thetas[i] = theta;
        phis[i] = phi;    

    }

    int n = 0;
    for (int _2j1=0; _2j1 < 2*jmax+1; _2j1++) {

        double j1 = 0.5 * _2j1;
        double j2 = 0.5 * _2j1;

        for (int j=0; j<min(_2j1,jmax) + 1; j++) {

            double B = calculate_B(j1,j2,j,nneighbors,distances,psis,thetas,phis);
            fp[n] = B;
            n++;
        }
    }

    cout<<"Lenght of fp: "<<n<<"\n";
    return fp;
}


double BispectrumCalculator::calculate_B(double j1, double j2, int j, int nneighbors, double* distances, double* psis, double* thetas, double* phis){

    double* mvalues = new double[2*j+1];
    for (int i=0; i<2*j+1; i++) mvalues[i] = j-i;

    complex<double> B = (0.0, 0.0);

    for (int i=0; i<2*j+1; i++) {
        double m = mvalues[i];
        for (int k=0; k<2*j+1; k++){
            int mprime = mvalues[k];
            complex<double> c = calculate_c(j,mprime,m,nneighbors,distances,psis,thetas,phis);

            double m1lower = fmin(j1,m+j2);
            double mp1lower = fmin(j1,mprime+j2);

            double m1 = fmax(-j1,m-j2);

            while (m1 < (m1lower + 0.5)){
                double mprime1 = fmax(-j1,mprime-j2);

                while (mprime1 < mp1lower + 0.5) {

                    complex<double> c1 = calculate_c(j1,mprime1,m1,nneighbors,distances,psis,thetas,phis);
                    complex<double> c2 = calculate_c(j2,mprime-mprime1,m-m1,nneighbors,distances,psis,thetas,phis);

                    B += get_CG_coefficient(j1,m1,j2,m-m1,j,m) * get_CG_coefficient(j1,mprime1,j2,mprime-mprime1,j,mprime)*conj(c)*c1*c2; 
                    mprime1 += 1;
               }
               m1 += 1;
            }    
        }
    }

    return real(B);
}

complex<double> BispectrumCalculator::calculate_c(double j, double mprime, double m, int nneighbors, double* distances, double* psis, double* thetas, double* phis){

    complex<double> c = (0., 0.);
    double weight = 14.0;
    for (int i=0; i<nneighbors; i++) {
        double dist = distances[i];
        double psi = psis[i];
        double theta = thetas[i];
        double phi = phis[i];
        c += weight*conj(U(j,m,mprime,psi,theta,phi))*cutoff_func(dist,cutoff);

    }

    return c;
}



complex<double> BispectrumCalculator::Wigner(double j,double m, double mprime,double alpha, double beta, double gamma) {

    complex<double> value = (0,0);
    complex<double> _j(0,1);

    if (abs(beta-M_PI/2) < TOLERANCE){

        for (int i=0; i<2*j+1; i++){
        
            if (i>j+mprime or i > j-m)
                break;
            else if (i < mprime-m)
                continue;

            value += pow(-1.0,i) * 1.0 * get_binomial(j+mprime,i)*get_binomial(j-mprime,i+m-mprime);
        }

        value *= pow(-1.0,m-mprime)*pow(factorial(j+m)*factorial(j-m)*1.0/((factorial(j+mprime)*factorial(j-mprime))),0.5)/pow(2.0,j);
        value *= exp(-1.0*_j*m*alpha)*exp(-1.0*_j*mprime*gamma);
    }
    else {
        
        double* mvalues = new double[int(2*j+1)];
        for (int i=0; i<int(2*j+1); i++) 
            mvalues[i] = j-i;
        
        for (int i=0; i<int(2*j+1); i++) {

            double mpp = mvalues[i];
            double temp1 = 0.0;
            for (int k=0; k<int(2*j+1); k++){
                if ((k>j+mpp) or (k > j-m))
                    break;
                else if (k < mpp - m)
                    continue;
                temp1 += pow(-1.0,k)*1.0*get_binomial(j+mpp,k)*get_binomial(j-mpp,k+m-mpp);
            } 

            temp1 *= pow(-1.0,m-mpp)*pow(factorial(j+m)*factorial(j-m)*1.0/((factorial(j+mpp)*factorial(j-mpp))),0.5)/pow(2.0,j);

            double temp2 = 0.0;
            for (int k=0; k<int(2*j+1); k++){
                if ((k>j-mprime) || (k > j-mpp))
                    break;
                else if (k < -mprime - mpp)
                    continue;
                 
                temp2 += pow(-1.0,k)*1.0*get_binomial(j-mprime,k)*get_binomial(j+mprime,k+mprime+mpp);
            } 

            temp2 *= pow(-1.0,mprime+mpp)*pow(factorial(j+mpp)*factorial(j-mpp)*1.0/((factorial(j-mprime)*factorial(j+mprime))),0.5)/pow(2.0,j);
            value += temp1*exp(_j*(-1.0*mpp)*beta)*temp2;
        }
        value *= (pow(_j*1.0,(2 * j - m - mprime))) * (pow(-1.0,(2 * m)));
        value *= exp((-1.0*m) *_j * alpha) * exp(_j * (-1.0*mprime) * gamma);
    }

    return value;
}


complex<double> BispectrumCalculator::U(double j, double m, double mprime, double omega, double theta, double phi) {

    complex<double> value = (0.0, 0.0);
    complex<double> _j(0,1);


    double* mvalues = new double[int(2*j+1)];
    for (int i=0; i<int(2*j+1); i++) 
        mvalues[i] = j-i;

    for (int i=0; i<int(2*j+1); i++) {
        double mpp = mvalues[i];
        value += Wigner(j,m,mpp,phi,theta,-phi)*exp(_j*(-1.0*mpp)*omega)*Wigner(j,mpp,mprime,phi,-theta,-phi);
    }
 
    return value;
}


double BispectrumCalculator::get_CG_coefficient(double a, double alpha, double b, double beta, double c, double gamma) {
    double result = 0.0;
    if (alpha + beta - gamma != 0) 
        return 0.0;

    else {

        double temp = factorial(a+alpha)*factorial(a-alpha)*factorial(b+beta)*factorial(b-beta)*
                      factorial(c+gamma)*factorial(c-gamma)*(2.0*c+1)*factorial(a+b-c)*
                      factorial(a-b+c)*factorial(-a+b+c)/factorial(a+b+c+1);
        double val1 = pow(temp,0.5);

        int zmin = max(int(a+beta-c),int(b-alpha-c));
        zmin = max(zmin,0);
        int zmax = min(int(b+beta),int(a-alpha));
        zmax = min(zmax,int(a+b-c));

        double val2  = 0.0;
        for (int z=zmin; z<zmax+1; z++){
            double val3 = factorial(z)*factorial(a+b-c-z)*factorial(a-alpha-z)*
                          factorial(b+beta-z)*factorial(c-b+alpha+z)*factorial(c-a-beta+z);
            val2 += pow(-1.0,z)/val3;
        }
        result = val1*val2;

    }
    return result;
}
