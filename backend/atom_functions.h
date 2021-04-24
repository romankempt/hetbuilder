#pragma once
#include <vector>

double2dvec_t lattice_points_in_supercell(int2dvec_t &superCellMatrix);

Atoms make_supercell(Atoms &prim, int2dvec_t &superCellMatrix);

Atoms rotate_atoms_around_z(Atoms &atoms, double &theta);

Atoms stack_atoms(Atoms bottom, Atoms top, double &weight, double &distance);