#include "Python.h"
#include "../../external/pybind11/include/pybind11/pybind11.h"

#include "coincidence_algorithm.h"
#include "test.h"

namespace py = pybind11;

PYBIND11_MODULE(pybackend, m)
{
    m.doc() = "backend c++ implementation"; // optional module docstring
    m.def("test_coincidence_algorithm", &test_coincidence_algorithm, "Test function");
    // py::bind_vector<std::vector<float>>(m, "FloatVector");
    // py::bind_vector<float_2d_vec>(m, "FloatVector2D");
    // m.def("find_coincidence", &find_coincidence, "A function that finds coincidence points.");
    // m.def("backend_routine", &backend_routine, "Full C++ routine to generate coincidence matrices.");
}