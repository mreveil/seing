#ifndef ATOM_H
#define ATOM_H

#include <string>

using namespace std;

class Atom {

    string atomtype;
    double x;
    double y;
    double z;
    double charge;
    double mass;


    public:
        Atom(void);
        Atom(string,double,double,double);
        ~Atom(void);

        string get_atom_type();
        double get_x();
        double get_y();
        double get_z();

};

#endif