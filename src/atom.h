#ifndef ATOM_H
#define ATOM_H

#include <string>

using namespace std;

//!  Class that holds information about an atom and its coordinates.

class Atom {

    string atomtype;
    double x;
    double y;
    double z;
    double charge;
    double mass;


    public:
        //! Default constructor.
        Atom(void);
        //! Constructor that creates an Atom object from its name and coordinates
        /*! \param attype name of atom (e.g. H, C, Si, Au, etc.) 
            \param cx cartesian coordinates - x value
            \param cy cartesian coordinates - y value
            \param cz cartesian coordinates - z value
        */
        Atom(string,double,double,double);
        ~Atom(void);

        //! Returns the name of this atom
        string get_atom_type();
 
        //! Returns the x coordinate of this atom
        double get_x();

        //! Returns the y coordinate of this atom        
        double get_y();

        //! Returns the z coordinate of this atom
        double get_z();

};

#endif
