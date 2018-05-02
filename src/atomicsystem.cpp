/*!

   \file atomicsystem.cpp
   Class that holds information about the molecular system including atom types and spatial coordinates.

   Created by Mardochee Reveil on 9/27/18
   Class to hold info about any molecular sytem
 */

#include <stdlib.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include <cstdlib>
#include <ctime>
#include <algorithm>

#include "atomicsystem.h"

using namespace std;

AtomicSystem::AtomicSystem(void) : atoms(NULL) {

    natoms = 0;
    xmin = 0.0;
    ymin = 0.0;
    zmin = 0.0;
    xmax = 0.0;
    ymax = 0.0;
    zmax = 0.0;
    xpbc = true;
    ypbc = true;
    zpbc = true;
    skin = 0.0;
}

AtomicSystem::AtomicSystem(string filename, bool pbcx, bool pbcy, bool pbcz, double skin): atoms(NULL), skin(skin) {

    xpbc = pbcx;
    ypbc = pbcy;
    zpbc = pbcz;

    build_from_file(filename);

}

void AtomicSystem::set_box_size(double xmin, double ymin, double zmin, double xmax, double ymax, double zmax) {

    xmin = xmin;
    ymin = ymin;
    zmin = zmin;
    xmax = xmax;
    ymax = ymax;
    zmax = zmax;

    int bad = 0;
    for (int i=0; i<natoms; i++) {
        double badx = false,bady = false,badz=false;
        if (atoms[i].get_x() < xmin || atoms[i].get_x() > xmax) {
            badx = true;
        }
        if (atoms[i].get_y() < ymin || atoms[i].get_y() > ymax) {
            bady = true;
        }
        if (atoms[i].get_z() < zmin || atoms[i].get_z() > zmax) {
            badz = true;
        }
        if (badx || bady || badz) bad++;

    }
    if (bad > 0)
        cerr<<"ERROR: "<<bad<<"atoms found outside box\n";
    return;
}

void AtomicSystem::build_from_file(string filename){

    int n = 0;
    string atomtype;
    double x,y,z;
    double minx=1e15, miny=1e15, minz=1e15, maxx=-1e15, maxy=-1e15, maxz=-1e15;

    string line;
    ifstream myfile(filename.c_str());

    if (!myfile.is_open()) {

        throw "Cannot open xyz file";
        return;
    }

    getline(myfile,line);
    istringstream buffer(line);
    buffer >> natoms;
    atoms = new Atom[natoms];
    getline(myfile,line);

    while(getline(myfile,line)) {

        istringstream buffer(line);
        buffer >> atomtype >> x >> y >> z;
        atoms[n] = Atom(atomtype,x,y,z);
        
        if (x > maxx) maxx = x;
        if (x < minx) minx = x;

        if (y > maxy) maxy = y;
        if (y < miny) miny = y;

        if (z > maxz) maxz = z;
        if (z < minz) minz = z;

        n++;

        if (n > natoms) cerr << "ERROR: too many atoms in xyz file"<<"\n";

    }

    xmin = minx;
    ymin = miny;
    zmin = minz;

    xmax = maxx;
    ymax = maxy;
    zmax = maxz;

    myfile.close();

}



AtomicSystem::~AtomicSystem(){

}


vector<string> AtomicSystem::get_atom_types(){

    PeriodicTable ptable = PeriodicTable();
    vector<string> supportedatomlist = ptable.get_element_list(); 

    vector<string> atomtypes;

    for (int i=0; i<natoms; i++) {

        string atype = atoms[i].get_atom_type();
        if (find(atomtypes.begin(), atomtypes.end(), atype) == atomtypes.end()){
          // Element not in atom types
            if (find(supportedatomlist.begin(), supportedatomlist.end(), atype) != supportedatomlist.end()){
                //Element is found in periodic table
                atomtypes.push_back(atype);
            }
            else {

                cerr<<"ERROR: Unknown atom type found in coordinates file: '"<<atype<<"'. Please add atom to periodictable.cpp and recompile SEING.\n";
                //TODO: exit program when this error is thrown!!!
            }
        }

    }

    return atomtypes;

}

// Function used for debugging purposes only

double AtomicSystem::check_square_distance(Atom A, Atom B) {

    double dx=0.0, dy=0.0, dz=0.0, sqdistance = 0.0;

    double x1 = A.get_x();
    double y1 = A.get_y();
    double z1 = A.get_z();
 
    double x2 = B.get_x();
    double y2 = B.get_y();
    double z2 = B.get_z();

    double x_size = xmax - xmin + skin;
    double y_size = ymax - ymin + skin;
    double z_size = zmax - zmin + skin;

    dx = x2 - x1;

    if (xpbc == true) {

        if (dx >   x_size*0.5) dx = dx - x_size;
        if (dx <= -x_size*0.5) dx = dx + x_size;  
    } 

    dy = y2 - y1;

    if (ypbc == true) {

        if (dy >   y_size*0.5) dy = dy - y_size;
        if (dy <= -y_size*0.5) dy = dy + y_size;  
    } 

    dz = z2 -z1;

    if (zpbc == true) {

        if (dz >   z_size*0.5) dz = dz - z_size;
        if (dz <= -z_size*0.5) dz = dz + z_size;  
    } 

   // if (x1==0.0 && y1==0.0 && z1==0.0 && x2==1.0 && y2==0.0 && z2 ==0.0)
   // cout<<"x: "<< x_size<<" "<< skin << " " << xmin<<" "<< dx<<" " << x1<<" "<< x2<<"\n";
   // cout<<"y: "<< y_size<<" "<< ymax<< " " << ymin<<" "<< dy<<" " << y1<<" "<< y2<<"\n";
   // cout<<"z: "<< z_size<<" "<< zmax<< " " << zmin<<" "<< dz<<" " << z1<<" "<< z2<<"\n\n";

    sqdistance = pow(dx,2.0) + pow(dy,2.0) + pow(dz,2.0);

    return sqdistance;

}


double AtomicSystem::get_distance_component(Atom A, Atom B, int direction) {

    double dist_component, size, coord1, coord2;
    bool pbc = false;

    if (direction == 0) {
        coord1 = A.get_x();
        coord2 = B.get_x();
        size = xmax - xmin + skin;
        if (xpbc == true) pbc = true;
    }
    else if (direction == 1) {
        coord1 = A.get_y();
        coord2 = B.get_y();
        size = ymax - ymin + skin;
        if (ypbc == true) pbc = true;

    }
    else if (direction == 2) {
        coord1 = A.get_z();
        coord2 = B.get_z();
        size = zmax - zmin + skin;
        if (zpbc == true) pbc = true;
    }
    else {
        cerr<<"ERROR: No such component for the distance between two atoms";
        return -1;
    }

    dist_component = coord2 - coord1;

    if (pbc == true ) {

        if (dist_component >   size*0.5) dist_component = dist_component - size;
        if (dist_component <= -size*0.5) dist_component = dist_component + size;  

    }

    return dist_component;
}

double AtomicSystem::get_distance_component(int atomID1, int atomID2, int direction) {

    Atom A = get_atom(atomID1);
    Atom B = get_atom(atomID2);

    return get_distance_component(A, B, direction);

}

double AtomicSystem::get_square_distance(Atom A, Atom B) {

    double dx = get_distance_component(A, B, 0);
    double dy = get_distance_component(A, B, 1);
    double dz = get_distance_component(A, B, 2);

    return pow(dx,2.0) + pow(dy,2.0) + pow(dz,2.0);
}

double AtomicSystem::get_square_distance(int id1, int id2) {

    Atom A = get_atom(id1);
    Atom B = get_atom(id2);

   // if (id1 == 0 and id2 == 3) return check_square_distance(A,B);

    return get_square_distance(A,B);

}

Atom AtomicSystem::get_atom(int id) {

    return atoms[id];
}

int AtomicSystem::get_n_atoms() {
    return natoms;

}

double AtomicSystem::get_xmin() {return xmin;}
double AtomicSystem::get_ymin() {return ymin;}
double AtomicSystem::get_zmin() {return zmin;}

double AtomicSystem::get_xmax() {return xmax;}
double AtomicSystem::get_ymax() {return ymax;}
double AtomicSystem::get_zmax() {return zmax;}

double AtomicSystem::get_xsize() {return xmax-xmin+skin;}
double AtomicSystem::get_ysize() {return ymax-ymin+skin;}
double AtomicSystem::get_zsize() {return zmax-zmin+skin;}

double AtomicSystem::get_xpbc() {return xpbc;}
double AtomicSystem::get_ypbc() {return ypbc;}
double AtomicSystem::get_zpbc() {return zpbc;}


