
#include <string>
#include <iostream>
#include <stdio.h>



#include "periodictable.h"

using namespace std;


PeriodicTable::PeriodicTable() { //Data from https://chem.libretexts.org/Reference/Periodic_Table_of_the_Elements
    //Please use the 1st ionization energy values with caution. Highly recommended to double check values

    nelements = 36;

    elements = new Element[36];

    Element Hydrogen    = {"Hydrogen",   "H",   1, 2.2,   1.00794,  1312.0};  elements[0] = Hydrogen;
    Element Helium      = {"Helium",     "He",  2, 0.0,   4.002602, 2372.3};  elements[1] = Helium;
    Element Lithium     = {"Lithium",    "Li",  3, 0.95,  6.941,     520.2};  elements[2] = Lithium;
    Element Berylium    = {"Berylium",   "Be",  4, 1.57,  9.012182,  899.5};  elements[3] = Berylium;
    Element Boron       = {"Boron",      "B",   5, 2.04, 10.811,     800.6};  elements[4] = Boron;
    Element Carbon      = {"Carbon",     "C",   6, 2.55, 12.01,     1086.5};  elements[5] = Carbon;
    Element Nitrogen    = {"Nitrogen",   "N",   7, 3.04, 14.0067,   1402.3};  elements[6] = Nitrogen;
    Element Oxygen      = {"Oxygen",     "O",   8, 3.44, 15.9994,   1313.9};  elements[7] = Oxygen;
    Element Fluorine    = {"Fluorine",   "F",   9, 3.96, 18.998403, 1651.0};  elements[8] = Fluorine;
    Element Neon        = {"Neon",       "Ne", 10, 0.0,  20.1797,   2080.7};  elements[9] = Neon;
    Element Sodium      = {"Sodium",     "Na", 11, 0.93, 22.98976,   495.8};  elements[10] = Sodium;
    Element Magnesium   = {"Magnesium",  "Mg", 12, 1.31, 24.3050,    737.7};  elements[11] = Magnesium;
    Element Aluminum    = {"Aluminum",   "Al", 13, 1.61, 26.98,      577.5};  elements[12] = Aluminum;
    Element Silicon     = {"Silicon",    "Si", 14, 1.9,  28.0855,    785.5};  elements[13] = Silicon;
    Element Phosphorus  = {"Phosphorus", "P",  15, 2.19, 30.97696,  1011.8};  elements[14] = Phosphorus;
    Element Sulfur      = {"Sulfur",     "S",  16, 2.58, 32.065,     999.6};  elements[15] = Sulfur;
    Element Chlorine    = {"Chlorine",   "Cl", 17, 3.16, 35.453,    1251.2};  elements[16] = Chlorine;
    Element Argon       = {"Argon",      "Ar", 18, 0.0,  39.948,    1520.6};  elements[17] = Argon;
    Element Potassium   = {"Potassium",  "K",  19, 0.82, 39.0983,    418.8};  elements[18] = Potassium;
    Element Calcium     = {"Calcium",    "Ca", 20, 1.00,  40.078,    589.8};  elements[19] = Calcium;
    Element Scandium    = {"Scandium",   "Sc", 21, 1.36,  44.9591,   631.1};  elements[20] = Scandium;
    Element Titanium    = {"Titanium",   "Ti", 22, 1.54,  47.867,    658.8};  elements[21] = Titanium;
    Element Vanadium    = {"Vanadium",   "Va", 23, 1.63,  50.9415,   650.9};  elements[22] = Vanadium;
    Element Chromium    = {"Chromium",   "Cr", 24, 1.54,  47.867,    658.8};  elements[23] = Chromium;
    Element Manganese   = {"Manganese",  "Mn", 25, 1.55, 54.93804,   717.3};  elements[24] = Manganese;
    Element Iron        = {"Iron",       "Fe", 26, 1.83, 55.845,     762.5};  elements[25] = Iron;
    Element Cobalt      = {"Cobalt",     "Co", 27, 1.91, 58.93319,   760.4};  elements[26] = Cobalt;
    Element Nickel      = {"Nickel",     "Ni", 28, 1.88, 58.6934,    737.1};  elements[27] = Nickel;
    Element Copper      = {"Copper",     "Cu", 29, 1.9,  63.546,     745.5};  elements[28] = Copper;
    Element Zinc        = {"Zinc",       "Zn", 30, 1.65, 65.38,      906.4};  elements[29] = Zinc;
    Element Gallium     = {"Gallium",    "Ga", 31, 1.81, 67.723,     578.8};  elements[30] = Gallium;
    Element Germanium   = {"Germanium",  "Ge", 32, 2.01, 72.64,      762.0};  elements[31] = Germanium;
    Element Arsenic     = {"Arsenic",    "As", 33, 2.18, 74.9216,    947.0};  elements[32] = Arsenic;
    Element Selenium     = {"Selenium",  "Se", 34, 2.35, 78.96,      941.0};  elements[33] = Selenium;

    Element Indium      = {"Indium",     "In", 49, 1.78, 114.818,    558.3};  elements[34] = Indium;

    Element Tungsten    = {"Tungsten",   "W",  74, 2.36, 183.84,     770.0};  elements[35] = Tungsten;
    

}

PeriodicTable::~PeriodicTable() {

}

vector<string> PeriodicTable::get_element_list() {
    vector<string> element_list;

    for (int i=0; i<nelements; i++) {

            element_list.push_back(elements[i].symbol);
        }
    return element_list;
}


int PeriodicTable::get_atomic_number(string symbol) {

    int anumber = -1;
    bool found = false;

    for (int i=0; i<nelements; i++) {

        if (elements[i].symbol == symbol){
            anumber = elements[i].atomic_number;
            found = true;
            break;
        }
    }

    if (!found) cerr<<"ERROR: Element '"<<symbol<<"' not found in periodic table\n";

    return anumber;
}

double PeriodicTable::get_electronegativity(string symbol) {

    double enegativity = -1.0;

    bool found = false;

    for (int i=0; i<nelements; i++) {

        if (elements[i].symbol == symbol){
            enegativity = elements[i].electronegativity;
            found = true;
            break;
        }
    }

    if (!found) cerr<<"ERROR: Element '"<<symbol<<"' not found in periodic table\n";


    return enegativity;
}

