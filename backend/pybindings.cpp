#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
#include <pybind11/iostream.h>

#include "logging_functions.h"
#include "atom_class.h"
#include "coincidence_algorithm.h"
#include "test.h"

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(int1dvec_t);
PYBIND11_MAKE_OPAQUE(int2dvec_t);
PYBIND11_MAKE_OPAQUE(double1dvec_t);
PYBIND11_MAKE_OPAQUE(double2dvec_t);

PYBIND11_MODULE(hetbuilder_backend, m)
{
    m.doc() = "C++ implementation of the coincidence algorithm."; // optional module docstring

    // datatype casting
    py::bind_vector<int1dvec_t>(m, "int1dVector");
    py::bind_vector<int2dvec_t>(m, "int2dVector");
    py::bind_vector<double1dvec_t>(m, "double1dVector");
    py::bind_vector<double2dvec_t>(m, "double2dVector");

    // output stream
    py::add_ostream_redirect(m, "ostream_redirect");

    // combined datatype castings
    py::bind_vector<std::vector<Interface>>(m, "Interfaces");

    // class bindings
    py::class_<Atoms>(m, "CppAtomsClass")
        .def(py::init<double2dvec_t &, double2dvec_t &, int1dvec_t &>())
        .def_readwrite("lattice", &Atoms::lattice)
        .def_readwrite("positions", &Atoms::positions)
        .def_readwrite("atomic_numbers", &Atoms::atomic_numbers)
        .def("standardize", &Atoms::standardize, "Spglib standardization.")
        .def("print", &Atoms::print, py::call_guard<py::scoped_ostream_redirect, py::scoped_estream_redirect>(), "Print information.")
        .def("scale_cell", &Atoms::scale_cell, "Scales cell.");

    py::class_<Interface>(m, "CppInterfaceClass")
        .def(py::init<Atoms &, Atoms &, Atoms &, double &, int2dvec_t &, int2dvec_t &, int &>())
        .def_readonly("bottom", &Interface::bottomLayer)
        .def_readonly("top", &Interface::topLayer)
        .def_readonly("stack", &Interface::stack)
        .def_readonly("angle", &Interface::angle)
        .def_readonly("M", &Interface::M)
        .def_readonly("N", &Interface::N)
        .def_readonly("spacegroup", &Interface::spaceGroup);

    py::class_<CoincidenceAlgorithm>(m, "CppCoincidenceAlgorithmClass")
        .def(py::init<Atoms &, Atoms &, int &, int &, double1dvec_t &, double &, double &, double &, int &, double &, double &>())
        .def("run", &CoincidenceAlgorithm::run, "Runs the coincidence algorithm.");

    // function definitions
    m.def("test_coincidence_algorithm",
          &test_coincidence_algorithm,
          py::call_guard<py::scoped_ostream_redirect, py::scoped_estream_redirect>(),
          "Test function to check that hte coincidence algorithm is executed.");

    m.def("test_atom_binding", &test_atom_binding, "Test function to bind Atom class.");

    m.def("get_number_of_omp_threads", &get_number_of_threads, "Returns number of available OMP threads.");

    m.def("cpp_make_supercell", &make_supercell, "C++ implementation to make supercell");

    m.def("cpp_lattice_points_in_supercell", &lattice_points_in_supercell, py::call_guard<py::scoped_ostream_redirect, py::scoped_estream_redirect>(), "C++ implementation to find lattice points in supercell matrix.");

    m.def("cpp_rotate_atoms_around_z", &rotate_atoms_around_z, py::call_guard<py::scoped_ostream_redirect, py::scoped_estream_redirect>(), "C++ implementation to rotate cell and atomic positions.");
}