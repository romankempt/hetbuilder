#pragma once
#include <cmath>
#include <vector>
#include <iostream>
#include <array>
#include <map>
#include <algorithm>
#include <set>

#include "spglib.h"

#include "logging_functions.h"
#include "math_functions.h"
#include "atom_class.h"
#include "atom_functions.h"

using std::sin, std::cos, std::sqrt, std::pow, std::abs;

typedef std::vector<int> int1dvec_t;
typedef std::vector<double> double1dvec_t;
typedef std::vector<std::vector<int>> int2dvec_t;
typedef std::vector<std::vector<double>> double2dvec_t;

/**
 * Solves the equation |Am - R(theta)Bn| < tolerance for a given angle theta.
 * 
 * The results are stored in a 2d vector of integers containing m1, m2, n1, n2.
 * OpenMP is employed to distribute the nested loops on threads, but an ordered construct 
 * has to be used to push back the vector for thread safety. 
 * 
 * The case of m1 = m2 = n1 = n2 is already removed, including the null vector.
 */
int2dvec_t find_coincidences(const double2dvec_t &A, const double2dvec_t &B, const double &theta, const int &Nmin, const int &Nmax, const double &tolerance);

/**
 * Constructs the independent pairs (m1,m2,m3,m4) and (n1,n2,n3,n4).
 * 
 * First loop is OpenMP parallel. Second one cannot be collapsed because j > i to avoid repititions.
 * 
 * All pairs with an absolute greatest common divisor different from 1 are removed,
 * because they correspond to scalar multiples of other smaller super cells.
 */
int2dvec_t find_unique_pairs(int2dvec_t &coincidences);

/**
 * Builds all supercells, applying the supercell matrices M and N and the Rotation R(theta).
 * 
 * The unit cell of the stack (interface) is given bei C = A + weight * (B - A).
 * The interfaces are standardized via spglib for the given symprec and angle_tolerance.
 * The loop over the supecell generation and standardization is OpenMP parallel.
 * 
 * Returns a vector of interfaces.
 */
std::vector<Interface> build_all_supercells(Atoms &bottom, Atoms &top, std::map<double, int2dvec_t> &AnglesMN, double &weight, double &distance, const int &no_idealize, const double &symprec, const double &angle_tolerance);

/**
 * Filters the interfaces.
 * 
 * Interfaces are considered equal if their spacegroup, area and number of atoms matches.
 * 
 * Returns a vector of interfaces.
 */
std::vector<Interface> filter_supercells(std::vector<Interface> &stacks);
