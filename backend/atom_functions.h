#pragma once
#include <vector>

typedef std::vector<std::vector<int>> int2dvec_t;
typedef std::vector<std::vector<double>> double2dvec_t;

double2dvec_t lattice_points_in_supercell(int2dvec_t &SuperCellMatrix);