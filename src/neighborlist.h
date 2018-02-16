#ifndef NEIGHBORLIST_H
#define NEIGHBORLIST_H

#include <vector>

#include "atom.h"
#include "atomicsystem.h"

using namespace std;

//! Class to generate and hold neighbor lists 
/*! This will divide the box into bins and for a given atom, loop over 
    atoms in neighboring bins to generate its neighbor list. This speeds
    up the process because it is no longer necessary to loop over all 
    atoms in the system.
    Also, for additional speed up, the neighbor list is kept in a linked list.   
*/
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
        //! Default Constructor
        NeighborList();

        //! Creates a NeighborList object based on an AtomicSystem
        /*! \param asystem reference to the atomic system for which the neighbor lists have 
                           to be generated.
            \param cutoff maximum distance for which two atoms are considered neighbors
            \param nxb number of bins in the x direction
            \param nyb number of bins in the y direction
            \param nzb number of bins in the z direction
        */
        NeighborList(AtomicSystem&,double,int,int,int,int);
        ~NeighborList();

        //! Function that actually generates the neighbor list for each atom
        void build();

        //! Returns a pointer to the list of atoms in given bin  
        int* get_atoms_in_bin(int);

        //! Returns the number of atoms in given bins
        int get_atoms_per_bin(int);

       // bool are_bins_neighbors(int,int);

        //! Retursn the bin number in which the atom with given x, y, and z coordinates belong
        int get_bin_number(double, double, double);

        //! Returns the neighboring bins for a bin givne by its index
        int* get_neighboring_bins(int);

        //! Returns a pointer to the list of neighbors of an atom given by its index
        int* get_neighbors(int);

        //! Returns a pointer to the neighbor list of a given atom  sorted by distance
        int* get_sorted_neighbors(int);

        //! Returns only atoms of a given type from the neighbor list of a given atom 
        int* get_sorted_neighbors(int,string);

        //! Returns only atoms of specific types from the neighbor list of a given atom 
        int* get_sorted_neighbors(int,vector<string>);

        //! Returns the number of neighbors for a given atom
        int get_n_neighbors(int);

        //! Returns the number of neighbors of a given type for a given atom
        int get_n_neighbors(int,string);

        //! Returns the number of neighbors of given types for a given atom
        int get_n_neighbors(int,vector<string>);

};

#endif
