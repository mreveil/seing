#include "inputs.h"
#include "atom.h"
#include "periodictable.h"
#include "atomicsystem.h"
#include "fingerprintgenerator.h"
#include "neighborlist.h"
#include "gaussiancalculator.h"
#include "bispectrumcalculator.h"
#include "zernikecalculator.h"
#include "genericlocalcalculator.h"
#include "utilities.h"
#include "agnicalculator.h"
#include <pybind11/pybind11.h>

namespace py = pybind11;

PYBIND11_MODULE(_pyseing, m){
    m.doc() = "This is the Python version of the SEING library. SEING is a C/C++ package for fingerprint calculations suitable for machine learning studies of molecular systems.\n\nSEING was developed in the Clancy Group (https://clancygroup.wse.jhu.edu/) by Mardochee Reveil (mr937@cornell.edu). The Python version is implemented by Man Kit Ao (mao2@jh.edu) and Divya Sharma (dsharm23@jh.edu).\n\nFingerprints (in this context) are numerical representations of chemical environments designed to be invariant under property-perseving operations such as permutation of atoms of the same nature, geometric rotation, etc. For more information on fingerprints in general and the ones currently implemented in SEING, please see the official documentation and user-guide.";
    m.def("read_prop_file", &read_prop_file, "Reads in the option file that specifies the type of fingerprints to generate and creates an fingerprintProperties object\n\nParameters:\n(arg0: str): the name of the input file\n\nReturns:\nfingerprintProperties: a seing4python.fingerprintProperties object\n");

    py::class_<fingerprintProperties>(m, "fingerprintProperties")
      .def_readwrite("box_size", &fingerprintProperties::box_size, "Returns the box size")
      .def_readwrite("output_file", &fingerprintProperties::output_file, "Returns the name of the outputfile")
      .def_readwrite("output_mode", &fingerprintProperties::output_mode, "Returns the mode of output");
    
    //py::class_<Atom, std::shared_ptr<Atom>>(m, "Atom")
    //	    .def(py::init<const std::string &, double, double, double>())
    //	    .def(py::init<>())
    //	    .def("get_atom_type", &Atom::get_atom_type)
    //	    .def("get_x", &Atom::get_x)
    //	    .def("get_y", &Atom::get_y)
    //	    .def("get_z", &Atom::get_z);
   
    //  py::class_<Element>(m, "Element");

    //    py::class_<PeriodicTable>(m, "PeriodicTable")
    //	    .def(py::init<>())
    //	    .def("get_atomic_number", &PeriodicTable::get_atomic_number)
    //	    .def("get_electronegativity", &PeriodicTable::get_electronegativity)
    //	    .def("get_element_list", &PeriodicTable::get_element_list);

    py::class_<AtomicSystem, std::shared_ptr<AtomicSystem>>(m, "AtomicSystem")
      .def(py::init<const std::string &, bool, bool, bool, double>(), "Reads in the .xyz file and severa; parameters and creates an AtomicSystem object\n\nParameters:\n(arg0: str): the name of the .xyz file\n(arg1: bool): whether to set periodic boundary conditions for x-sides\n(arg2: bool): whether to set periodic boundary conditions for y-sides\n(arg3: bool): whether to set periodic bounardy conditions for z-sides\n(arg4: double): skin\n\nReturns:\nAtomicSystem: a seing4python.AtomicSystem object\n")
      .def(py::init<>())
	    // .def("atoms", &AtomicSystem::atoms)
	    // .def("build_from_file", &AtomicSystem::build_from_file)
      .def("set_box_size", &AtomicSystem::set_box_size, "Sets the box size of the atomic system\n\nParameters:\n(arg0: float): xmin\n(arg1: float): ymin\n(arg2: float): zmin\n(arg3: float): xmax\n(arg4: float): ymax\n(arg5: float): zmax\n")
      //.def("get_distance_component", py::overload_cast<int, int, int>(&AtomicSystem::get_distance_component))
      //.def("get_distance_component", py::overload_cast<Atom, Atom, int>(&AtomicSystem::get_distance_component))
      //.def("get_square_distance", py::overload_cast<Atom, Atom>(&AtomicSystem::get_square_distance))
      //.def("get_square_distance", py::overload_cast<int, int>(&AtomicSystem::get_square_distance))
      //.def("check_square_distance", &AtomicSystem::check_square_distance)
      .def("get_atom_types", &AtomicSystem::get_atom_types, "Returns the types of atoms in the atomic system")
      //.def("get_atom", &AtomicSystem::get_atom)
      .def("get_n_atoms", &AtomicSystem::get_n_atoms, "Returns the nnumber of atoms in the atomic system")
      .def("get_xsize", &AtomicSystem::get_xsize, "Returns the length of the x dimension")
      .def("get_ysize", &AtomicSystem::get_ysize, "Returns the length of the y dimension")
      .def("get_zsize", &AtomicSystem::get_zsize, "Returns the length of the z dimension")
      .def("get_xmin", &AtomicSystem::get_xmin, "Returns the minimum x-value where there is an atom")
      .def("get_ymin", &AtomicSystem::get_ymin, "Returns the minimum y-value where there is an atom")
      .def("get_zmin", &AtomicSystem::get_zmin, "Returns the minimum z-value where there is an atom")
      .def("get_xmax", &AtomicSystem::get_xmax, "Returns the maximum x-value where there is an atom")
      .def("get_ymax", &AtomicSystem::get_ymax, "Returns the maximum y-value where there is an atom")
      .def("get_zmax", &AtomicSystem::get_zmax, "Returns the maximum z-value where there is an atom")
      .def("get_xpbc", &AtomicSystem::get_xpbc, "Returns the current periodic boundary conditions for the x-sides of the box")
      .def("get_ypbc", &AtomicSystem::get_ypbc, "Returns the current periodic boundary conditions for the y-sides of the box")
      .def("get_zpbc", &AtomicSystem::get_zpbc, "Returns the current periodic boundary conditions for the z-sides of the box");
    
    py::class_<FingerprintGenerator>(m, "FingerprintGenerator")
      .def(py::init<AtomicSystem&, fingerprintProperties>(), "Creates a FingerprintGenerator object\n\nParameters:\n(argo: AtomicSystem): the AtomicSystem object\n(arg1: fingerprintProperties): the fingerprintProperties object\nReturns:\nFingerprintGenerator: a FingerprintGenerator object\n")
      .def("write2file", &FingerprintGenerator::write2file, "Write the fingerprints to a file\n\nParameters:\n(arg0: str): the name of the file\n(arg1: str): the mode of output\n")
	    .def_readwrite("natoms", &FingerprintGenerator::natoms)
	    .def_readwrite("fsize", &FingerprintGenerator::fsize);

    // m.def("CalculateFP", &CalculateFP, "Calculates fingerprints and writes the output to a file, the corresponding AtomicSystem object is declared with periodic boundary conditions for all sides and a skin of 2.01778\n\nParameters:\n(arg0: str): the name of the .xyz file\n(arg1: str): the name of the option file\n");

    // m.def("CalculateFP_w_fname", &CalculateFP, "Calculates fingerprints and writes the output to a file, the corresponding AtomicSystem object is declared with periodic boundary conditions for all sides and a skin of 2.01778\n\nParameters:\n(arg0: str): the name of the .xyz file\n(arg1: str): the name of the option file\n(arg2: str): the name of the output file, can be different from the one in the option file, useful in iteration to create multiple fingerprint files\n");

    
    // these two are overloaded methods of CalculateFP
    //    m.def("CalculateFP", py::overload_cast<const std::string &, const std::string &>(&CalculateFP), "Calculates fingerprints and writes the output to a file, the corresponding AtomicSystem object is declared with periodic boundary conditions for all sides and a skin of 2.01778\n\nParameters:\n(arg0: str): the name of the .xyz file\n(arg1: str): the name of the option file\n(arg2: str): (OPTIONAL) the name of the output file, can be different from the one in the option file, useful in iteration to create multiple fingerprint files\n");

    //    m.def("CalculateFP", py::overload_cast<const std::string &, const std::string &, const std::string &>(&CalculateFP),  "Calculates fingerprints and write \
s the output to a file, the corresponding AtomicSystem object is declared with periodic boundary conditions for all sides and a skin o\
f 2.01778\n\nParameters:\n(arg0: str): the name of the .xyz file\n(arg1: str): the name of the option file\n(arg2: str): \
(OPTIONAL) the name of the output file, can be different from the one in the option file, useful in iteration to create multiple finge\
rprint files\n");
}
