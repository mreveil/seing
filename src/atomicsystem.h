/*  Created by Mardochee Reveil on 9/27/18.
  Class to hold info about any molecular sytem.
*/

#ifndef ATOMICSYSTEM_H
#define ATOMICSYSTEM_H

#include <string>
#include <vector>

#include "atom.h"
#include "periodictable.h"

using namespace std;

//!  Class that holds information about the molecular system including atom types and spatial coordinates.

class AtomicSystem {

    private: 
        int natoms;

	Atom *atoms;

        double xmin, ymin, zmin;
        double xmax, ymax, zmax;
        double skin;
        bool xpbc, ypbc, zpbc; 

	// void build_from_file(string);
    public:

	//transferred from private:
	// Atom *atoms;
	void build_from_file(string);

        //! Default constructor 
        AtomicSystem(void);

        //! Constructor that creates the atomicsystem object from a coordinate file
        /*! \param filename name of coordinate file 
            \param pbcx boundary condition for x 
            \param pbcy boundary condition for y
            \param pbcz boundary condition for z
            \param skin size of the skin around the box for used in neighbor list generation
        */
        AtomicSystem(string,bool,bool,bool,double);
        ~AtomicSystem(void);

   


        //! Sets the boundaries of the box xmin, ymin, zmin, xmax, ymax, zmax
        void set_box_size(double, double, double, double, double, double);

        //! Calculates the x (=0), y (=1) or z (=2) component of the distance between two atoms
        /*! \param A  the first atom object.
            \param B  the second atom object.
            \param direction  the component of the distance to return 0=x, 1=y, 2=z 
        */
        double get_distance_component(Atom, Atom, int);

        double get_distance_component(int, int, int);

        //! Calculates the square distance between two atoms
        double get_square_distance(Atom,Atom);

        //! Calculates the square distance between two atoms given by their index in the AtomicSystem
        double get_square_distance(int,int);

        double check_square_distance(Atom A, Atom B); //For debugging only

        //! Finds unique atomic species and return them as a vector of atom names ordered by atomic numbers
        vector<string> get_atom_types(); 

        //! Returns an Atom given by its index in the atomic system
        Atom get_atom(int);

        //! Returns the total number of atoms in the system
        int get_n_atoms();

        double get_xsize();
        double get_ysize();
        double get_zsize();

        double get_xmin();
        double get_ymin();
        double get_zmin();

        double get_xmax();
        double get_ymax();
        double get_zmax();

        double get_xpbc();
        double get_ypbc();
        double get_zpbc();

   
};

#endif
