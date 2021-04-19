#pragma once
#include <vector>
#include <iostream>
#include <map>

typedef std::vector<std::vector<int>> int2dvec_t;
typedef std::vector<std::vector<double>> double2dvec_t;

// Prints 2d vector.
template <typename T>
void print_2d_vector(const std::vector<std::vector<T>> &vec);
template <typename T>
void print_2d_vector(std::vector<std::vector<T>> &vec);

// Prints number of OpenMP threads.
void log_number_of_threads();

// Prints map of single-valued double key with a 2d vector of ints.
void print_map_key_2d_vector(std::map<double, int2dvec_t> const &m);