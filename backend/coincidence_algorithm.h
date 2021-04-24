#pragma once
#include "logging_functions.h"
#include "math_functions.h"
#include "atom_class.h"
#include "atom_functions.h"
#include "interface_class.h"

/**
 * Class definition of the lattice coincidence algorithm.
 * 
 * Executed by the run() method.
 */
class CoincidenceAlgorithm
{
public:
    Atoms primitive_bottom;
    Atoms primitive_top;

    CoincidenceAlgorithm(Atoms cPrimitiveBottom,
                         Atoms cPrimitiveTop)
    {
        primitive_bottom = cPrimitiveBottom;
        primitive_top = cPrimitiveTop;
    };

    int2dvec_t find_coincidences(double2dvec_t &A, double2dvec_t &B, double &theta, int &Nmin, int &Nmax, double &tolerance);

    int2dvec_t find_unique_pairs(int2dvec_t &coincidences);

    std::vector<Interface> build_all_supercells(Atoms &bottom, Atoms &top, std::map<double, int2dvec_t> &AnglesMN, double &weight, double &distance, int &no_idealize, double &symprec, double &angle_tolerance);

    std::vector<Interface> filter_supercells(std::vector<Interface> &stacks);

    std::vector<Interface> run(int cNmax = 10,
                               int cNmin = -10,
                               double1dvec_t cAngles = {0.0, 30.0, 60.0, 90.0},
                               double cTolerance = 0.01,
                               double cWeight = 0.5,
                               double cDistance = 4.0,
                               int cNoIdealize = 0,
                               double cSymPrec = 1e-5,
                               double cAngleTolerance = 5.0);
};
