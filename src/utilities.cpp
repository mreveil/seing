#define _USE_MATH_DEFINES

#include <math.h>

#include "atom.h"
#include "atomicsystem.h"
#include "utilities.h"

using namespace std;


double cutoff_func(double r,double cutoff){
    double result = 0.0;
    if (r <= cutoff) 
       result += 0.5*(1+cos(M_PI*r/cutoff));
    return result;
}

double cutoff_func_prime(double r,double cutoff){
    double result = 0.0;
    if (r <= cutoff) 
       result = (-0.5*M_PI/cutoff)*sin(M_PI*r/cutoff);
       //result += 0.5*(1+cos(M_PI*r/cutoff));
    return result;
}

double get_min_distance(double a,double b, double size, bool pbc){

    double dist = b - a;

    if (pbc == true) {

        if (dist >   size*0.5) dist = dist - size;
        if (dist <= -size*0.5) dist = dist + size;  
    }

    return dist;

}

double calculate_norm(double x, double y, double z) {

    return pow(pow(x,2)+pow(y,2)+pow(z,2),0.5);

}

double* get_vector(AtomicSystem atomicsystem,int i, int j) {

    double* Rij = new double[3];

    double xsize = atomicsystem.get_xsize();
    double ysize = atomicsystem.get_ysize();
    double zsize = atomicsystem.get_zsize();

    double xpbc = atomicsystem.get_xpbc();
    double ypbc = atomicsystem.get_ypbc();
    double zpbc = atomicsystem.get_zpbc();

    Atom Ai=atomicsystem.get_atom(i);
    Atom Aj=atomicsystem.get_atom(j);

    Rij[0] = get_min_distance(Ai.get_x(),Aj.get_x(),xsize,xpbc);
    Rij[1] = get_min_distance(Ai.get_y(),Aj.get_y(),ysize,ypbc);
    Rij[2] = get_min_distance(Ai.get_z(),Aj.get_z(),zsize,zpbc); 

    return Rij;
}


int Kronecker(int i, int j){

    if (i==j) return 1;
    else return 0;

}

double dot(double* vec1, double* vec2) {

    return (vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2]);

}

//int BispectrumCalculator::calculate_factorial(int n)
//{
//  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
//}

unsigned int factorial(unsigned int n){

    unsigned int ret = 1;
    for(unsigned int i = 1; i <= n; ++i)
        ret *= i;
    return ret;
}

double get_binomial(int n, int k) {

    return factorial(n)/(factorial(k)*factorial(n-k));
}
