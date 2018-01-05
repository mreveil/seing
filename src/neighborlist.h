#ifndef NEIGHBORLIST_H
#define NEIGHBORLIST_H

#include <vector>

#include "atom.h"
#include "atomicsystem.h"

using namespace std;

class NeighborList {

    int natoms;
    int totalbins;

    int nxbins;
    int nybins;
    int nzbins;

    double xbinsize, ybinsize, zbinsize;
    double minxbox, minybox, minzbox;
    double xboxsize, yboxsize, zboxsize;

    int maxneighbors;
    
    AtomicSystem atomicsystem;
    double cutoffsq;

    int **neighboringbins;
    int *heads;
    int *binlist;
    int *atomsperbin;

    int** neighbors;
    int*  neighborsperatom;


    bool initialize_binning();
    void find_neighbors(int);
    bool is_bin_valid(double, double, double, int, int, int);
    void find_neighboring_bins();

    public:
        NeighborList();
        NeighborList(AtomicSystem&,double,int,int,int,int);
        ~NeighborList();

        void build();
        int* get_atoms_in_bin(int);
        int get_atoms_per_bin(int);
       // bool are_bins_neighbors(int,int);
        int get_bin_number(double, double, double);

        int* get_neighboring_bins(int);
        int* get_neighbors(int);

        int* get_sorted_neighbors(int);
        int* get_sorted_neighbors(int,string);
        int* get_sorted_neighbors(int,vector<string>);

        int get_n_neighbors(int);
        int get_n_neighbors(int,string);
        int get_n_neighbors(int,vector<string>);

};

#endif
