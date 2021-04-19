#pragma once
#include <vector>
#include <map>

using std::vector;

template <class T>
void print_2d_vector(const vector<vector<T>> &vec);

void log_number_of_threads(void);

void print_map_key_2d_vector_of_ints(std::map<double, vector<vector<int>>> const &m);