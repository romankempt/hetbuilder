#pragma once
#include <vector>

typedef std::vector<std::vector<int>> int2dvec_t;
typedef std::vector<std::vector<double>> double2dvec_t;

double2dvec_t lattice_points_in_supercell(int2dvec_t &SuperCellMatrix);

/**
 * Generate a supercell by applying a SuperCellMatrix to
    the input atomic configuration prim.
*/
Atoms make_supercell(Atoms &prim, int2dvec_t &SuperCellMatrix);

Atoms rotate_atoms_around_z(Atoms &atoms, const double &theta);