#pragma once
#include <vector>

double2dvec_t lattice_points_in_supercell(int2dvec_t &superCellMatrix);

/**
 * Generate a supercell by applying an integer superCellMatrix to
    the input atomic configuration prim.
*/
Atoms make_supercell(const Atoms prim, int2dvec_t &superCellMatrix);

/**
 * Applies a rotation matrix R for the given angle theta in degree
 * to the atomic configuration and cell to rotate around the z-axis.
*/
Atoms rotate_atoms_around_z(Atoms &atoms, const double &theta);

/**
 * Stacks two Atoms bottom and top on top of eacher other
 * with the new unit cell being C = A + weight * (B-A) and an interlayer distance given by distance.
*/
Atoms stack_atoms(Atoms &bottom, Atoms &top, double &weight, double &distance);