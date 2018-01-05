#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <time.h>

#include "fingerprintgenerator.h"
#include "neighborlist.h"
#include "atomicsystem.h"
#include "atom.h"
#include "inputs.h"


using namespace std;

int main(int argc, char*argv[]) {

    time_t start = time(0);

    string currentdir, xyzfilename, outputfile, propfilename;
    AtomicSystem asys;
    fingerprintProperties fpproperties;
    currentdir = argv[0];

    if (argc <= 1) {
        cout<<"Not enough arguments. Usage: seing xyzfilename inputfilename\n";
        return -1;  
    }

    xyzfilename = argv[1];
    outputfile = argv[2];

  // Reading input file

    if (argc == 3) {
        propfilename = argv[2];
        try {
            cout<<"Reading fingerprint properties input file...\n\n";
            fpproperties = read_prop_file(propfilename);
        } catch (const char* msg) {
            cerr<<msg<<"\n";
            cout<<"ERROR: Could not read input file.\n";
            return -1;
        }
    }

    // Reading xyz file

    try {
        cout <<"Reading xyz file.............";
        asys = AtomicSystem(xyzfilename, true, true, true, 2.01778); //1.46 for InGaAs
        cout << "\n\tNumber of atoms: " << asys.get_n_atoms()<<"\n";

    } catch (const char* msg) {
        cerr<<msg<<"\n";
        cout<<"ERROR:\n";
        return -1;
    }

    if (fpproperties.is_box_size_provided) {
        double xmin = fpproperties.box_size[0];
        double ymin = fpproperties.box_size[1];
        double zmin = fpproperties.box_size[2];
        double xmax = fpproperties.box_size[3];
        double ymax = fpproperties.box_size[4];
        double zmax = fpproperties.box_size[5];
        asys.set_box_size(xmin, ymin, zmin, xmax, ymax, zmax);
    }
  //  else fpproperties = default_properties();

    FingerprintGenerator fpgenerator = FingerprintGenerator(asys,fpproperties);
    cout<<"Writing fingerprint to file..."<<fpproperties.output_file<<"\n";
    fpgenerator.write2file(fpproperties.output_file,fpproperties.mode);

    return 1;
}
