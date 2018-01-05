// Created by Mardochee Reveil on 9/27/18
// Class to hold info about an atom

#include <string>

#include "atom.h"

Atom::Atom(void) {

    atomtype = "O";
    x = 0.0;
    y = 0.0;
    z = 0.0; 
    charge = 0.0;

}

Atom::Atom(string attype,double cx, double cy, double cz) {

    atomtype = attype;
    x = cx;
    y = cy;
    z = cz;
    charge = 0.0;
}

Atom::~Atom(void){

}

string Atom::get_atom_type(){return atomtype;}
double Atom::get_x(){return x;}
double Atom::get_y(){return y;}
double Atom::get_z(){return z;}
