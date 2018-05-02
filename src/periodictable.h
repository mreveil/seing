#ifndef PERIDOICTABLE_H 
#define PERIDOICTABLE_H

//
// Author: Mardochee Reveil
// Date Created: 12/17/17
//


#include <string>
#include <vector>

using namespace std;


struct Element {

    string name;
    string symbol;
    int atomic_number; 
    double electronegativity;
    double mass;
    double ionization_energy; //in kJ/mol
};


class PeriodicTable {

    Element *elements;
    int nelements;

    public:

        PeriodicTable();
        ~PeriodicTable();

        int get_atomic_number(string);
        double get_electronegativity(string);
        vector<string> get_element_list();

};

#endif